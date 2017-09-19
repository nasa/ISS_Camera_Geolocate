"""
Initial Commit 20170330
Timothy Lang (timothy.j.lang@nasa.gov)

References
----------
Schenk, T., 2005: Introduction to Photogrammetry. Department of Civil and
Environmental Engineering and Geodetic Science, The Ohio State University,
100 pp.
[Available online at http://www.mat.uc.pt/~gil/downloads/IntroPhoto.pdf]

Saplano, M. R. P., S. Bilanow, and W. Berg, 2010: SSM/I and SSMIS Stewardship
Code Geolocation Algorithm Theoretical Basis. CSU Technical Report, 33 pp.
[Available online at
http://rain.atmos.colostate.edu/FCDR/doc/CSU_FCDR_geolocation_tech_report.pdf]

"""
from __future__ import print_function  # for debug
import numpy as np
import datetime as dt
from astropy.time import Time
from sgp4.earth_gravity import wgs84
from sgp4.io import twoline2rv
import time

# WGS-84 Attributes
DEF_RE = 6378.137  # Equatorial Radius of the Earth (km)
DEF_RP = 6356.755  # Polar Radius of the Earth (km)
DEF_RM = 6371.0  # Mean Radius of the Earth (km)


def get_fov(d, f):
    """
    Helper function, if you want to figure out field of view
    for a particular camera.

    Parameters
    ----------
    d : float
        Length of camera film/sensor dimension (mm)
    f : float
        Focal length of camera lens (mm)

    returns
    -------
    fov : float
        Camera field of view (degrees)
    """
    return np.rad2deg(2.0 * np.arctan(d / (2.0 * f)))


def get_ifov(fov, pixels):
    """
    Helper function, if you want to figure out instantaneous field of view
    for a particular camera.

    Parameters
    ----------
    fov : float
        Field of view of camera, in degrees
    pixels : float
        Number of pixels in CCD

    Returns
    -------
    ifov : float
        Instantaneous field of view, in degrees
    """
    return fov / pixels


def get_camera_vector(f, hsize, vsize, hpixels, vpixels):
    """
    Given the focal length of the camera, the physical size of the focal plane,
    and the size of the image, obtain a vector of photogrammetric information
    for each pixel.

    Parameters
    -----------
    f : float
        Focal length of camera (mm)
    hsize : float
        Horizontal length of CCD imager (mm)
    vsize : float
        Vertical length of CCD imager (mm)
    hpixels : float
        Number of CCD pixels in horizontal dimension
    vpixels : float
        Number of CCD pixels in vertical dimension
    fov : float
        Field of view of camera, in degrees
    pixels : float
        Number of pixels in CCD

    Returns
    -------
    camera_vector : list
        Camera pointing vectors as a list of three-element vectors
    """
    x = np.linspace(-hsize/2.0, hsize/2.0, hpixels)
    y = np.linspace(-vsize/2.0, vsize/2.0, vpixels)
    tanx = x / (-f)
    tany = y / (-f)
    xp, yp = np.meshgrid(tanx, tany)
    xp = xp.ravel()
    yp = yp.ravel()
    return np.matrix([xp, yp, -1.0*np.ones(len(xp))])


def get_A_matrix(roll, pitch, yaw):
    """
    Compute the camera rotation matrix, given roll/pitch/yaw.
    Note: ISS attitude info is not typically publicly available, so
    the spacecraft and camera rotation matrices are combined in this
    function. If ISS attitude info is known, both rotation matrices will
    need to be multiplied together.

    Parameters
    ----------
    roll : float
        Roll angle (degrees)
    pitch : float
        Pitch angle (degrees)
    yaw : float
        Yaw angle (degrees)

    Returns
    -------
    A : numpy.matrix
        3x3 3D rotation matrix, with the multiplication ordered
        towards the largest angle adjustments
    """
    Asort = np.argsort([np.abs(roll), np.abs(pitch), np.abs(yaw)])
    roll = np.deg2rad(roll)
    pitch = np.deg2rad(pitch)
    yaw = np.deg2rad(yaw)
    Ar = np.matrix([[np.cos(roll), 0, np.sin(roll)],
                    [0, 1, 0],
                    [-np.sin(roll), 0, np.cos(roll)]])
    Ap = np.matrix([[1, 0, 0],
                    [0, np.cos(pitch), -np.sin(pitch)],
                    [0, np.sin(pitch), np.cos(pitch)]])
    Ay = np.matrix([[np.cos(yaw), -np.sin(yaw), 0],
                    [np.sin(yaw), np.cos(yaw), 0],
                    [0, 0, 1]])
    Alist = [Ar, Ap, Ay]
    return Alist[Asort[2]] * Alist[Asort[1]] * Alist[Asort[0]]


def get_vector_magnitude(vect):
    """
    Compute the magnitude of a vector.

    Parameters
    ----------
    vect : 1-dimensional array of floats
        Vector

    Returns
    -------
    vector_magnitude : float
        Magnitude of vector
    """
    return np.sqrt(np.sum(np.array(vect) * np.array(vect)))


def get_f(RE=DEF_RE, RP=DEF_RP):
    """
    Obtain the flattening factor of the Earth ellipsoid.
    Note: Default is WGS84 ellipsoid

    Parameters
    ----------
    RE : float
        Equatorial radius of the Earth (km)
    RP : float
        Polar radius of the Earth (km)

    Returns
    -------
    f : float
        Flattening factor of the Earth ellipsoid
    """
    return 1.0 - (RP / RE)


def get_fprime(f, P, Rm=DEF_RM):
    """
    Obtain the modified flattening factor of the Earth ellipsoid.

    Parameters
    ----------
    f : float
        Flattening factor of the Earth ellipsoid
    P : 1-dimensional array of floats
        Spacecraft position vector (Px, Py, Pz) in km
    Rm : float
        Mean radius of the Earth (km)

    Returns
    -------
    fprime : float
        Modified flattening factor
    """
    Pm = get_vector_magnitude(P)
    return (Rm * (1 - f)**2 + Pm - Rm) / Pm


def get_Mz(P, f, Rm=6371.0088):
    """
    Compute the z component of the nadir-to-GCI rotation matrix M

    Parameters
    ----------
    P : 1-dimensional array of floats
        Spacecraft position vector (Px, Py, Pz) in km
    f : float
        Flattening factor of the Earth ellipsoid
    Rm : float
        Mean radius of the Earth (km)

    Returns
    -------
    Mz : 1-dimensional array of floats
        Mz component of the nadir-to-GCI rotation matrix M
    """
    fp = get_fprime(f, P, Rm=Rm)
    denom = np.sqrt(P[2]**2 + fp**2 * (P[0]**2 + P[1]**2))
    Mz1 = -P[0] * fp / denom
    Mz2 = -P[1] * fp / denom
    Mz3 = -P[2] / denom
    return np.array([Mz1, Mz2, Mz3])


def get_My(Mz, V):
    """
    Compute the y component of the nadir-to-GCI rotation matrix M

    Parameters
    ----------
    Mz : 1-dimensional array of floats
        Mz component of the nadir-to-GCI rotation matrix M
    V : 1-dimensional array of floats
        Spacecraft velocity vector (Vx, Vy, Vz)

    Returns
    -------
    My : 1-dimensional array of floats
        My component of the nadir-to-GCI rotation matrix M
    """
    T = np.cross(Mz, V)
    Tm = get_vector_magnitude(T)
    return T / Tm


def get_Mx(My, Mz):
    """
    Compute the x component of the nadir-to-GCI rotation matrix M

    Parameters
    ----------
    My : 1-dimensional array of floats
        My component of the nadir-to-GCI rotation matrix M
    Mz : 1-dimensional array of floats
        Mz component of the nadir-to-GCI rotation matrix M

    Returns
    -------
    Mx : 1-dimensional array of floats
        Mx component of the nadir-to-GCI rotation matrix M
    """
    return np.cross(My, Mz)


def get_M(P, V, f, Rm=DEF_RM):
    """
    Compute the matrix M to convert from nadir to
    geocentric inertial (GCI) coordinates.

    Parameters
    ----------
    P : 1-dimensional array of floats
        Spacecraft position vector (Px, Py, Pz) in km
    V : 1-dimensional array of floats
        Spacecraft velocity vector (Vx, Vy, Vz)
    f : float
        Flattening factor of the Earth ellipsoid
    Rm : float
        Mean radius of the Earth (km)

    Returns
    -------
    M : 3x3 matrix of floats
        nadir-to-GCI rotation matrix M
    """
    Mz = get_Mz(P, f, Rm=Rm)
    My = get_My(Mz, V)
    Mx = get_Mx(My, Mz)
    # Arrangement still uncertain, needs to be [Mx, My, Mz]
    return np.matrix([Mx, My, Mz]).T


def get_Di(M, A, Ds, S=None):
    """
    Obtain the IFOV vector in GCI coordinates.

    Parameters
    ----------
    M : numpy.matrix
        3x3 nadir-to-GCI rotation matrix M
    A : numpy.matrix
        3x3 Spacecraft Attitude Matrix
    Ds : 1-dimensional array of floats
        3-element IFOV vector in sensor coordinates
    S : numpy.matrix or None
        Optional 3x3 Sensor Alignment Matrix

    Returns
    -------
    Di : 1-dimensional array of floats
        3-element IFOV vector in GCI coordinates
    """
    if S is None:
        return M * A.T * Ds
    else:
        return M * A.T * S.T * Ds


def quadratic_equation(a, b, c):
    """
    Compute the quadratic equation

    Parameters
    ---------
    a : float
        Coefficient a of a*x**2 + b*x + c
    b : float
        Coefficient a of a*x**2 + b*x + c
    c : float
        Coefficient a of a*x**2 + b*x + c

    Returns
    -------
    quadratic : tuple, float, or NoneType
        Depending on nature of square root argument, return the appropriate
        solution to the quadratic equation.
    """
    root = b**2 - 4 * a * c
    if root > 0:
        return (-b + np.sqrt(root)) / (2 * a), (-b - np.sqrt(root)) / (2 * a)
    if root == 0:
        return -b / (2 * a)
    else:
        return None


def get_G(P, Di, RE=DEF_RE, RP=DEF_RP):
    """
    Get the target position vector in GCI coordinates

    Parameters
    ----------
    P : 1-dimensional array of floats
        Spacecraft position vector (Px, Py, Pz) in km
    Di : 1-dimensional array of floats
        3-element IFOV vector in GCI coordinates
    RE : float
        Equatorial radius of the Earth (km)
    RP : float
        Polar radius of the Earth (km)

    Returns
    -------
    G : 1-dimensional array of floats
        3-element target position vector in GCI coordinates
    """
    a = (Di[0]**2 + Di[1]**2) / RE**2 + Di[2]**2 / RP**2
    b = (2.0 * P[0] * Di[0] + 2.0 * P[1] * Di[1]) / RE**2 + \
        2.0 * P[2] * Di[2] / RP**2
    c = (P[0]**2 + P[1]**2) / RE**2 + P[2]**2 / RP**2 - 1.0
    d = quadratic_equation(a, b, c)
    if d is None:
        return np.array([np.nan, np.nan, np.nan])
    if np.size(d) == 2:
        da = np.abs(d)
        d = d[np.argmin(da)]
    return np.array([P[0] + d * Di[0], P[1] + d * Di[1], P[2] + d * Di[2]])


def get_geodetic_lat_lon(G, RE=DEF_RE, RP=DEF_RP):
    """
    Get the geodetic latitude and longitude of each image pixel.
    Longitude must then be corrected for Greenwich Hour Angle.

    Parameters
    ----------
    G : 1-dimensional array of floats
        3-element target position vector in GCI coordinates
    RE : float
        Equatorial radius of the Earth (km)
    RP : float
        Polar radius of the Earth (km)

    Returns
    -------
    lat : float
        Geodetic latitude (degrees)
    lon : float
        Geodetic longitude (degrees)
    """
    f = get_f(RE=RE, RP=RP)
    lat = np.arctan2(G[2], ((1 - f)**2 * np.sqrt(G[0]**2 + G[1]**2)))
    lon = np.arctan2(G[1], G[0])
    return np.rad2deg(lat), np.rad2deg(lon)  # No correction for GHA yet


def get_tles_and_datetimes(filen):
    """
    Read a file containing a list of TLEs and convert every TLE to a
    datetime object.

    Parameters
    ----------
    filen : str
        String name for a TLE file

    Returns
    -------
    tle : ndarray
        TLE strings. Every other string is a new TLE.
    dts : ndarray
        Array of datetime objects from corresponding to every TLE in the file.
        Half the length of tle.
    """
    sec_to_micro = 1e6
    day_to_sec = 24 * 3600
    tle = []
    dts = []
    f = open(filen, 'r')
    lines = f.readlines()
    f.close()
    for i, line in enumerate(lines):
        tle.append(line.strip('\n'))
        if i % 2 == 0:
            dstr = line.split()[3]
            year = int(dstr[0:2]) + 2000  # Won't work prior to 2000
            day = np.double(dstr[2:])
            musec_of_yr = int(sec_to_micro * day_to_sec * (day - 1.0))
            dt1 = dt.datetime(year, 1, 1) + \
                dt.timedelta(microseconds=musec_of_yr)
            dts.append(dt1)
    tle = np.array(tle)
    dts = np.array(dts)
    return tle, dts


def get_closest_tle(dtobj, tle, dts, verbose=True):
    """
    Search thru a list of TLEs and find the closest in time to the input
    datetime.

    Parameters
    ----------
    dtobj : datetime.datetime object
        Datetime object to match to the closest TLE.
    tle : ndarray
        TLE strings. Every other string is a new TLE.
    dts : ndarray
        Array of datetime objects from corresponding to every TLE in the file.
        Half the length of tle.

    Returns
    -------
    tle_out : ndarray of str
        Two-Line Element (TLE) set that is closest in time to the
        user-provided datetime.
    """
    index = np.argmin(np.abs(dtobj - dts))
    if verbose:
        print(dts[index], '\n', tle[2 * index], '\n', tle[2 * index + 1])
    return dts[index], tle[2 * index], tle[2 * index + 1]


def loop_over_each_pixel_and_geolocate_it(vect, M, A, P, verbose=True):
    """
    Given a vector of camera pixels, loop over them all and geolocate
    each given understanding of camera pointing and spacecraft
    positioning. Return these positions as 1-D arrays.
    Note: This function can take a few minutes to complete,
    depending on resolution of the image.

    Parameters
    ----------
    vect : numpy.matrix
        Camera photogrammetric pixel matrix
    M : numpy.matrix
        nadir-to-GCI coordinates matrix
    A : numpy.matrix
        Camera rotation matrix
    P : numpy.ndarray
        Spacecraft position vector

    Other Parameters
    ----------------
    verbose : bool
        True - Print out helpful monitoring info

        False - Don't do this

    Returns
    -------
    lats : numpy.ndarray
        Array of geodetic latitudes for all image pixels
    lons : numpy.ndarray
        Array of uncorrected geodetic longitudes for all pixels
    """
    if verbose:
        bt = time.time()
    lats = []
    lons = []
    Di = []
    for i, Ds in enumerate(vect.T):
        Ds = np.vstack(np.asarray(Ds)[0])
        Di.append(get_Di(M, A, Ds))
        di = np.squeeze(np.asarray(Di[i]))
        G = get_G(P, di)
        lat, lon = get_geodetic_lat_lon(G)
        lats.append(lat)
        lons.append(lon)
        if i % 500000 == 0:
            if verbose:  # Plot out every 500k-th pixel
                print('i =', i, G, lat, lon)
    lats = np.array(lats)
    lons = np.array(lons)
    if verbose:
        print(time.time() - bt, 'seconds to process')
    return lats, lons


def get_GHA(dt1):
    """
    Longitude geolocation must be corrected for Greenwich Hour
    Angle. This is determined from the datetime associated
    with the image to be geolocated.

    Parameters
    ----------
    dt1 : datetime.datetime
        Input datetime object for image to be geolocated

    Returns
    -------
    GHA : float
        Greenwich Hour Angle, use this to correct image longitudes
    """
    # Convert datetime to astropy.time.Time object
    times = [dt1.strftime('%Y-%m-%dT%H:%M:%S.%f')]
    t1 = Time(times, format='isot', scale='utc')

    # Get Greenwich Mean Sidereal Time
    d_a = t1.jd[0] - 2451545.0
    t = np.double(d_a) / 36525.0
    GMST = 280.46061837 + 360.98564736629 * d_a

    # Meeus Ch21 approximate formulas
    Om = (125.04452 - 1934.136261 * t) % 360
    L = (280.4665 + 36000.7698 * t) % 360
    L1 = (218.3165 + 481267.8813 * t) % 360
    eps = (23.439 - 0.0000004 * t) % 360

    # All of these quantities must be brought into the range 0 to 360
    # using the modulo function or similar.
    dp = -17.2 * np.sin(np.deg2rad(Om)) - 1.32 * np.sin(np.deg2rad(2*L)) - \
        0.23 * np.sin(np.deg2rad(2*L1)) + 0.21 * np.sin(np.deg2rad(2*Om))
    de = 9.2 * np.cos(np.deg2rad(Om)) + 0.57 * np.cos(np.deg2rad(2*L)) + \
        0.1 * np.cos(np.deg2rad(2*L1)) - 0.09 * np.cos(np.deg2rad(2*Om))
    e = eps + de / 3600.0

    # The correction to GMST in seconds of time is given by
    correction_s = dp * np.cos(np.deg2rad(e)) / 15.0

    # and so, in degrees
    correction_d = dp * np.cos(np.deg2rad(e)) / 3600.0
    GHA = GMST % 360 + correction_d
    return GHA
