from libc.math cimport sqrt, atan
cimport numpy as np


def GeodeticFromKmCTRS(double[:] CTRS_Pos):
    cdef:
        double ae = 6378.137
        double b = 6356.752
        double RAD2DEG = 57.29578
        double f, e, XX, YY, ZZ, p
        double q, r, s, t, u, v, w
        double k, D, lon, lat, h

    f = (ae - b) / ae
    e = sqrt(2.0 * f - f**2)
    XX = CTRS_Pos[0]
    YY = CTRS_Pos[1]
    ZZ = CTRS_Pos[2]
    p = ((XX**2) + (YY**2)) / (ae**2)
    q = ((1.0 - (e**2)) / (ae**2)) * (ZZ**2)
    r = (p + q - (e**4)) / 6.0
    s = ((e**4) * p * q) / (4.0 * (r**3))
    t = pow((1.0 + s + sqrt(s * (2.0 + s))), (1.0 / 3.0))
    u = r * (1.0 + t + (1.0 / t))
    v = sqrt((u**2) + ((e**4) * q))
    w = (e**2) * ((u + v - q) / (2.0 * v))
    k = sqrt(u + v + (w**2)) - w
    D = (k * sqrt((XX**2) + (YY**2))) / (k + (e**2))
    lon = RAD2DEG * (2.0 * atan(YY / (XX + sqrt((XX**2) + (YY**2)))))
    lat = RAD2DEG * (2.0 * atan(ZZ / (D + sqrt((D**2) + (ZZ**2)))))
    h = ((k + (e**2) - 1.0) / k) * sqrt((D**2) + (ZZ**2))    
    return lat, lon, h
