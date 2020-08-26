# -*- coding: utf-8 -*-
"""
Title/Version
-------------
ISS Camera Geolocate
iss_camera_geolocate v0.2
Developed & tested with Python 3.6+
Last changed 08/25/2020


Author
------
Timothy Lang
NASA MSFC
timothy.j.lang@nasa.gov
(256) 961-7861


Overview
--------
This is a Python software library that facilitates the geolocation of photographs and video frames from the International Space Station (ISS). The library provides functions that take camera and pointing information, along with publicly available ISS position information, and then will geolocate every pixel of the photograph in latitude and longitude. This enables geospatial analysis of astronaut photography from Earth, including pictures of clouds, lightning, coastlines, city lights, etc. Many images available from https://earth.jsc.nasa.gov/ can be fully geolocated using this software.

The code now also enables geolocation of the ISS Lightning Imaging Sensor (LIS) background imagery datasets. These data are available from http://dx.doi.org/10.5067/LIS/ISSLIS/DATA206 and http://dx.doi.org/10.5067/LIS/ISSLIS/DATA207.

To access this module, in the main source directory:
python setup.py install

Then in your program:
import iss_camera_geolocate as icg


Notes
-----
Dependencies: numpy, datetime, astropy, sgp4, cartopy, cython, xarray
"""

import os
import ast
import setuptools

try:
    import numpy
except ImportError:
    raise RuntimeError("Cannot find NumPy. Please install it first "
                       "or use pip >= 10, which will do so automatically.")

USE_CYTHON = True
if USE_CYTHON:
    from setuptools import setup, Extension  # analysis:ignore
    try:
        import Cython  # analysis:ignore
    except ImportError:
        raise RuntimeError("Cannot find Cython. Please install it first "
                           "or use pip >= 10, which will do so automatically.")
    from Cython.Build import cythonize  # analysis:ignore
    from Cython.Compiler import Options  # analysis:ignore
    Options.language_level = '2'

# Get current location
HERE = os.path.abspath(os.path.dirname(__file__))
# Pull the header into a variable
DOCLINES = __doc__.split('\n')
# Get packages
PACKAGES = setuptools.find_packages()

if USE_CYTHON:
    EXT = '.pyx'
EXTENSIONS = [Extension(PACKAGES[0] + '.GeodeticFromKmCTRS',
                        [PACKAGES[0] + '/GeodeticFromKmCTRS' + EXT])]
INCLUDE_DIRS = [numpy.get_include(), '.']

if USE_CYTHON:
    EXTENSIONS = cythonize(EXTENSIONS)

setup(
    name='ISS_Camera_Geolocate',
    version='0.2',
    author='Timothy Lang',
    author_email='timothy.j.lang@nasa.gov',
    packages=PACKAGES,
    license='LICENSE.md',
    description=DOCLINES[1],
    long_description=__doc__,
    install_requires=['sgp4', 'astropy', 'xarray', 'cartopy', 'cython'],
    ext_modules=EXTENSIONS,
    include_dirs=INCLUDE_DIRS,
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Environment :: Console"
        ],
)
