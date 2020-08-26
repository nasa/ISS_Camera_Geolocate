"""
This code enables geolocation of the ISS LIS background datasets
http://dx.doi.org/10.5067/LIS/ISSLIS/DATA206
http://dx.doi.org/10.5067/LIS/ISSLIS/DATA207

Notable required packages: xarray, cython, cartopy

Timothy Lang
timothy.j.lang@nasa.gov
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime as dt
import cartopy.crs as ccrs
import xarray
import time
import itertools
from .GeodeticFromKmCTRS import GeodeticFromKmCTRS  # cython 5x speedup

# Below are global constants based on empirically derived values
# specific to the ISS LIS optics and integration on the space station.
MAX_HT = 13.5
MAGNIFICATION_FACTOR = 1.01
ROTATE_FACTOR = 0
LEFT_RIGHT_FACTOR = -0.0022
UP_DOWN_FACTOR = 0.0205
MATRIX_VAL = np.array(
    [[-0.1322, -0.0518, 0.9911],
     [0.0268, 0.9953, 0.1055],
     [-0.9907, 0.0413, -0.1299]])


class Quaternion(object):

    """
    This class simplifies passing of quaternion information
    """

    def __init__(self, w, x, y, z):
        self.w = w
        self.x = x
        self.y = y
        self.z = z


class Geo(object):

    """
    This class simplifies passing of geolocation information
    """

    def __init__(self, lat, lon, alt):
        self.lat = lat
        self.lon = lon
        self.alt = alt


class Ephemeris(object):

    """
    This class simplifies passing of ephemeris data
    """

    def __init__(self, new_pv, new_vv, new_tm, index):
        self.state_vector = np.concatenate(
                [[], new_pv[index]])
        self.state_vector = np.concatenate(
                [self.state_vector, new_vv[index]])
        self.translation_matrix = new_tm[index]


def get_every_pointing_vector():

    """
    ; ------------------------------------------------------------------
    ;  Polynomial fit of lab measurements to get pointing vector
    ;  from x and y pixel location
    ; ------------------------------------------------------------------
    ; Look_Vector  = output = unit vector direction from CCD axis
    ;                         to geolocated position (a global variable)
    ; x_pixel      = input  = CCD pixel location in x direction (0-127)
    ; y_pixel      = input  = CCD pixel location in y direction (0-127)
    ;
    """
    Look_Vector = np.zeros((128, 128, 3), dtype='double')

    # ISS-LIS optics
    coeff = np.zeros(4, dtype='double')
    coeff[0] = 1.4754537
    coeff[1] = -0.36224695
    coeff[2] = -0.088939824
    coeff[3] = -0.28203806

    for i, j in itertools.product(range(128), range(128)):

        x = (i - 63.5) / 127.0
        y = (127 - j - 63.5) / 127.0
        xy = np.sqrt(x * x + y * y)
        convert = coeff[0] + coeff[1] * xy + coeff[2] * xy * xy + \
            coeff[3] * xy * xy * xy

        Look_Vector[i][j][0] = x * convert * MAGNIFICATION_FACTOR
        Look_Vector[i][j][1] = y * convert * MAGNIFICATION_FACTOR
        Look_Vector[i][j][2] = np.sqrt(
            1.0 - (Look_Vector[i][j][0] * Look_Vector[i][j][0] +
                   Look_Vector[i][j][1] * Look_Vector[i][j][1]))

    return Look_Vector


def QuaternionFromMatrix(m):

    tr = m[0, 0] + m[1, 1] + m[2, 2]

    if tr > 0:
        S = np.sqrt(tr + 1.0) * 2
        q = Quaternion(0.25 * S, (m[2, 1] - m[1, 2]) / S,
                       (m[0, 2] - m[2, 0]) / S, (m[1, 0] - m[0, 1]) / S)

    elif m[0, 0] > m[1, 1] and m[0, 0] > m[2, 2]:
        S = np.sqrt(1.0 + m[0, 0] - m[1, 1] - m[2, 2]) * 2
        q = Quaternion((m[2, 1] - m[1, 2]) / S, 0.25 * S,
                       (m[0, 1] + m[1, 0]) / S, (m[0, 2] + m[2, 0]) / S)

    elif m.e[1, 1] > m.e[2, 2]:
        S = np.sqrt(1.0 + m[1, 1] - m[0, 0] - m[2, 2]) * 2
        q = Quaternion((m[0, 2] - m[2, 0]) / S, (m[0, 1] + m[1, 0]) / S,
                       0.25 * S, (m[1, 2] + m[2, 1]) / S)

    else:
        S = np.sqrt(1.0 + m[2, 2] - m[0, 0] - m[1, 1]) * 2
        q = Quaternion((m[1, 0] - m[0, 1]) / S, (m[0, 2] + m[2, 0]) / S,
                       (m[1, 2] + m[2, 1]) / S, 0.25 * S)

    return q


def Vector_Add(A, B):
    """
    C-code legacy function
    """
    return A + B


def Vector_Mag(A):
    """
    C-code legacy function, could get replaced by an equivalent numpy function
    """
    return np.sqrt(
        A[0] * A[0] +
        A[1] * A[1] +
        A[2] * A[2])


def Vector_Cross(A, B):
    """
    This C-code legacy function is surprisingly faster than np.cross, go figure
    """
    C = np.zeros(3, dtype='double')
    C[0] = A[1] * B[2] - A[2] * B[1]
    C[1] = A[2] * B[0] - A[0] * B[2]
    C[2] = A[0] * B[1] - A[1] * B[0]
    return C


def MatrixFromYPR(ypr):

    cy = np.cos(ypr[0])
    sy = np.sin(ypr[0])
    cp = np.cos(ypr[1])
    sp = np.sin(ypr[1])
    cr = np.cos(ypr[2])
    sr = np.sin(ypr[2])

    MYaw = np.zeros((3, 3), dtype='double')
    MPitch = np.zeros((3, 3), dtype='double')
    MRoll = np.zeros((3, 3), dtype='double')

    # Define Yaw matrix
    MYaw[0][0] = cy
    MYaw[0][1] = -sy
    MYaw[0][2] = 0.0
    MYaw[1][0] = sy
    MYaw[1][1] = cy
    MYaw[1][2] = 0.0
    MYaw[2][0] = 0.0
    MYaw[2][1] = 0.0
    MYaw[2][2] = 1.0

    # Define Pitch matrix
    MPitch[0][0] = cp
    MPitch[0][1] = 0.0
    MPitch[0][2] = sp
    MPitch[1][0] = 0.0
    MPitch[1][1] = 1.0
    MPitch[1][2] = 0.0
    MPitch[2][0] = -sp
    MPitch[2][1] = 0.0
    MPitch[2][2] = cp

    # Define Roll matrix
    MRoll[0][0] = 1.0
    MRoll[0][1] = 0.0
    MRoll[0][2] = 0.0
    MRoll[1][0] = 0.0
    MRoll[1][1] = cr
    MRoll[1][2] = -sr
    MRoll[2][0] = 0.0
    MRoll[2][1] = sr
    MRoll[2][2] = cr

    return np.dot(np.dot(MYaw, MPitch), MRoll)


def LVLHMatrixFrom(inPos, inVel):

    Mlvlh = np.zeros((3, 3), dtype='double')
    posMag = Vector_Mag(inPos)
    speed = Vector_Mag(inVel)
    crossPr = Vector_Cross(inPos, inVel)

    for i in range(3):
        Mlvlh[0][i] = inVel[i] / speed  # LVLH X in geocentric frame
        # LVLH Y in geocentric frame
        Mlvlh[1][i] = crossPr[i] / ((-posMag) * speed)
        Mlvlh[2][i] = inPos[i] / (-posMag)  # LVLH Z in geocentric frame

    return Mlvlh


def YPRFromQuaternion(q):
    """
    ****************** DESCRIPTION ********************
    This subroutine converts quaternions to yaw, pitch, roll (Euler angles)
    """
    v = np.zeros(3, dtype='double')
    v[0] = np.arctan2(2.0 * (q.x * q.y + q.w * q.z),
                      q.w * q.w + q.x * q.x - q.y * q.y - q.z * q.z)
    v[1] = np.arcsin(-2.0 * (q.x * q.z - q.w * q.y))
    v[2] = np.arctan2(2.0 * (q.y * q.z + q.w * q.x),
                      q.w * q.w - q.x * q.x - q.y * q.y + q.z * q.z)
    return v


def get_earth_intersection(number_to_process, look_vector, Eph, ISSGD,
                           MlisOrientation, Mypr, Mlvlh, inPosCTRS):

    """
    /* ****************** DESCRIPTION ******************** */
    /* This subroutine computes the surface position of each */
    /* pixel handed to it. */
    /* It is based on CATS code. */
    /* This subroutine utilizes the subroutines above. */
    Note: Old code derived from C, could use more Python-ification.
    """

    location = np.zeros((1, 2), dtype='double')
    for i in range(number_to_process):

        PixelLOS = np.zeros(3, dtype='double')
        for j in range(3):
            PixelLOS[j] = look_vector[i][j]

        Alt = ISSGD.alt - 3.0
        RangeSpan = 5.0 * ISSGD.alt
        P2GD = Geo(-90.0, -90.0, -1.0)
        # This loop is legacy code that will be modified soon
        for z in np.arange(Alt, Alt + RangeSpan, 1.0):
            P2T = np.dot(MlisOrientation, z * PixelLOS)
            P2LVLH = np.dot(Mypr, P2T)  # Transform to LVLH
            P2g = np.dot(P2LVLH, Mlvlh)  # Transform to geocentric
            P2 = Vector_Add(inPosCTRS, P2g)
            output = GeodeticFromKmCTRS(P2)
            P2GD = Geo(output[0], output[1], output[2])
            # MAX_HT = "Cloud-Top Height"
            # P2GD.alt is nearly always negative to start, so break hits early
            if P2GD.alt < MAX_HT:
                break

        location[i][0] = P2GD.lat
        location[i][1] = P2GD.lon

    return location


def get_interpolated_matrices(lis_background, index):

    """
    ISS LIS background files contain TAI93 time stamps for each image.
    However, ephemeris information is provided at different time steps.
    This function interpolates the ephemeris data to the desired background
    image's time. The ephemeris time step is only ~1 second, so simple
    linear interpolation is used.
    """

    time1 = pd.to_datetime(
        lis_background.bg_data_summary_TAI93_time.data).to_pydatetime()
    time2 = pd.to_datetime(
        lis_background.bg_info_TAI93_time.data).to_pydatetime()
    # Work in floating point seconds since the earliest reference time.
    # np.interp does not play nicely with datetime/timedelta objects.
    min_time = np.min([time1[0], time2[0]])
    td1 = np.array([(t1 - min_time).total_seconds() for t1 in time1])
    td2 = np.array([(t2 - min_time).total_seconds() for t2 in time2])

    # Position vector
    new_pv = np.zeros(
        (time1.shape[0], lis_background.bg_info_position_vector.data.shape[1]),
        dtype='double')
    # Velocity vector
    new_vv = np.zeros_like(new_pv)
    # Transformation matrix
    new_tm = np.zeros(
        (time1.shape[0],
         lis_background.bg_info_transform_matrix.data.data.shape[1]),
        dtype='double')

    # Do the interpolation
    for i in range(9):
        if i < 3:
            new_pv[:, i] = np.interp(
                td1, td2, lis_background.bg_info_position_vector.data[:, i])
            new_vv[:, i] = np.interp(
                td1, td2, lis_background.bg_info_velocity_vector.data[:, i])
        new_tm[:, i] = np.interp(
                td1, td2, lis_background.bg_info_transform_matrix.data[:, i])

    return Ephemeris(new_pv, new_vv, new_tm, index)


def run_iss_lis_geolocation(lis_background, index, verbose=True):
    """
    Main function that drives the geolocation processing.

    Parameters
    ----------
    lis_background : xarray.Dataset object
        LIS background xarray Dataset object from xarray.open_dataset()
    index : int
        Background image index to use

    Other Parameters
    ----------------
    verbose : bool
        True - Print out helpful monitoring info

        False - Don't do this

    Returns
    -------
    locations : numpy.ndarray
        16384 x 2 array of geolocated pixels (lats = index 0, lons = index 1)
        Use numpy.reshape(locations[:, i], (128, 128)) to recover 2D structure
        Be sure to transpose the original 128x128 image data after geolocation
    """

    if verbose:
        bt = time.time()
    Eph = get_interpolated_matrices(lis_background, index)
    lookvecSC = np.zeros((1, 3), dtype='double')
    locations = []

    inPosCTRS = np.zeros(3, dtype='double')
    inVelCTRS = np.zeros(3, dtype='double')
    Mypr = np.zeros((3, 3), dtype='double')
    MlisOrientation = np.zeros((3, 3), dtype='double')
    for i in range(3):
        for j in range(3):
            Mypr[i, j] = Eph.translation_matrix[i*3 + j]
            MlisOrientation[i, j] = MATRIX_VAL[i, j]
        inPosCTRS[i] = Eph.state_vector[i] * 0.001
        inVelCTRS[i] = Eph.state_vector[i+3] * 0.001

    output = GeodeticFromKmCTRS(inPosCTRS)
    ISSGD = Geo(output[0], output[1], output[2])
    Mlvlh = LVLHMatrixFrom(inPosCTRS, inVelCTRS)

    VlisOrientation = np.zeros(3, dtype='double')
    VlisOrientation[0] = 3.182 - 0.0485 + \
        0.0485 * Eph.state_vector[5] / 6000.0 + ROTATE_FACTOR
    VlisOrientation[1] = -0.020 + 0.015 - \
        0.0042 * Eph.state_vector[5] / 6000.0 + LEFT_RIGHT_FACTOR
    VlisOrientation[2] = 0.020 - 0.051 + UP_DOWN_FACTOR
    MlisOrientation = MatrixFromYPR(VlisOrientation)

    lv = get_every_pointing_vector()
    for i, j in itertools.product(range(128), range(128)):
        lookvecSC[0][0] = lv[i][j][0]
        lookvecSC[0][1] = lv[i][j][1]
        lookvecSC[0][2] = lv[i][j][2]
        locations.append(
            get_earth_intersection(
                1, lookvecSC, Eph, ISSGD,
                MlisOrientation, Mypr, Mlvlh, inPosCTRS))

    if verbose:
        print((time.time() - bt) / 1.0, 'seconds to run')
    return np.squeeze(locations)


def plot_series_of_backgrounds(lis_background, index, save=None, cmap='bone'):
    """
    This function plots 10 panels of ISS LIS background images.
    It is a quick way to review ISS LIS backgrounds to identify images worth
    geolocating.

    Parameters
    ----------
    lis_background : xarray.Dataset object
        LIS background xarray Dataset object from xarray.open_dataset()
    locations : numpy.ndarray
        16384 x 2 array of geolocated pixels (lats = index 0, lons = index 1)
        Use numpy.reshape(locations[:, i], (128, 128)) to recover 2D structure
        Be sure to transpose the original 128x128 image data after geolocation
    index : int
        Background image index to use

    Other Parameters
    ----------------
    save : str or None
        File to save plot to (via matplotlib.pyplot.savefig)
    cmap : str
        Colormap to use

    Returns
    -------
    None

    """
    bdts = pd.to_datetime(
        lis_background.bg_data_summary_TAI93_time.data).to_pydatetime()
    fig = plt.figure(figsize=(9, 20))
    for i in range(10):
        ax = fig.add_subplot(5, 2, i+1)
        ax.imshow(lis_background.bg_data[index + i, :, :], cmap=cmap)
        ax.set_title(str(index+i) + ' - ' +
                     bdts[index+i].strftime('%Y%m%d %H:%M:%S'))
    plt.tight_layout()
    if save is not None:
        plt.savefig(save)


def plot_geolocated_quicklook(lis_background, locations, index, save=None,
                              cmap='bone', delt=5.8, alpha=0.5,
                              layer='ASTER_GDEM_Color_Shaded_Relief'):
    """
    This function uses cartopy to plot a simple quicklook of a
    non-geolocated and then geolocated background. The function also provides
    a useful template for decoding the geolocation information provided
    by run_iss_lis_geolocation().

    Parameters
    ----------
    lis_background : xarray.Dataset object
        LIS background xarray Dataset object from xarray.open_dataset()
    locations : numpy.ndarray
        16384 x 2 array of geolocated pixels (lats = index 0, lons = index 1)
        Use numpy.reshape(locations[:, i], (128, 128)) to recover 2D structure
        Be sure to transpose the original 128x128 image data after geolocation
    index : int
        Background image index to use

    Other Parameters
    ----------------
    save : str or None
        File to save plot to (via matplotlib.pyplot.savefig)
    cmap : str
        Colormap to use
    delt : float
        Half-width of geolocated quicklook panel (degrees)
    alpha : float
        Alpha value (0-1) of background image on top of WTMS imagery
    layer: str or None
        WTMS layer from https://map1c.vis.earthdata.nasa.gov/wmts-geo/wmts.cgi
        None means no layer will be plotted

    Returns
    -------
    None

    """
    url = 'https://map1c.vis.earthdata.nasa.gov/wmts-geo/wmts.cgi'
    fig = plt.figure(figsize=(16, 8))
    projection = ccrs.PlateCarree()
    lats = locations[:, 0]
    lons = locations[:, 1]

    ax = fig.add_subplot(121, projection=projection)
    ax.imshow(lis_background.bg_data.data[index, :, :], cmap='bone')
    ax.set_title('(a) Raw ISS LIS Background')

    ax = fig.add_subplot(122, projection=projection)
    ext = [np.mean(lons)-delt,
           np.mean(lons)+delt,
           np.mean(lats)-delt,
           np.mean(lats)+delt]
    ax.set_extent(ext)
    if type(layer) is str:
        ax.add_wmts(url, layer)
    ax.coastlines(resolution='10m')
    ax.pcolormesh(lons.reshape((128, 128)), lats.reshape((128, 128)),
                  lis_background.bg_data.data[index, :, :].T, cmap=cmap,
                  alpha=alpha, transform=projection)
    ax.set_title('(b) Geolocated ISS LIS Background')

    gl = ax.gridlines(linestyle='--', draw_labels=True)
    gl.ylabels_right = []
    gl.xlabels_top = []
    if save is not None:
        plt.savefig(save, bbox_inches='tight')
