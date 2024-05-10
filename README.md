ISS Camera Geolocate README
----------------------------
This is a Python software library that facilitates the geolocation of photographs and video frames from the International Space Station (ISS). The library provides functions that take camera and pointing information, along with publicly available ISS position information, and then will geolocate every pixel of the photograph in latitude and longitude. This enables geospatial analysis of astronaut photography from Earth, including pictures of clouds, lightning, coastlines, city lights, etc. Many images available from https://eol.jsc.nasa.gov/  can be fully geolocated using this software.

The code also enables geolocation of the ISS Lightning Imaging Sensor (LIS) background imagery datasets. These data are available from http://dx.doi.org/10.5067/LIS/ISSLIS/DATA211. This code has been updated to handle the ISS LIS relocation after 7 July 2022, as well has the time prior to relocation. ISS LIS mission data began on 1 March 2017 and ended on 16 November 2023.

ISS Camera Geolocate Installation
---------------------------------
ISS Camera Geolocate works under Python 3.12+ on most Mac/Linux setups. Windows installation and other Python versions are currently untested.

In the main source directory:  
`python -m pip install .`

The following dependencies need to be installed first:

- A robust version of Python w/ most standard scientific packages (e.g., `numpy`, `datetime`, `astropy`, etc.) - Get one for free [here.](https://store.continuum.io/cshop/anaconda/)
- [SGP4](https://pypi.python.org/pypi/sgp4/)
-Cartopy
-Cython
-Xarray

Using ISS Camera Geolocate
--------------------------
To access everything:
```
import iss_camera_geolocate as icg
```

Demonstration notebooks are in the notebooks directory.

Latest release info:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2585824.svg)](https://doi.org/10.5281/zenodo.2585824)
