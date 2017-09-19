ISS Camera Geolocate README
----------------------------
This is a Python software library that facilitates the geolocation of photographs and video frames from the International Space Station (ISS). The library provides functions that take camera and pointing information, along with publicly available ISS position information, and then will geolocate every pixel of the photograph in latitude and longitude. This enables geospatial analysis of astronaut photography from Earth, including pictures of clouds, lightning, coastlines, city lights, etc. Many images available from https://earth.jsc.nasa.gov/ can be fully geolocated using this software.

ISS Camera Geolocate Installation
---------------------------------
ISS Camera Geolocate works under Python 2.7 and 3.6 on most Mac/Linux setups. Windows installation and other Python versions are currently untested.

In the main source directory:  
`python setup.py install`

The following dependencies need to be installed first:

- A robust version of Python w/ most standard scientific packages (e.g., `numpy`, `datetime`, `astropy`, etc.) - Get one for free [here.](https://store.continuum.io/cshop/anaconda/)
- [SGP4](https://pypi.python.org/pypi/sgp4/)


Specific import calls in the ISS Camera Geolocate source code:

```
from __future__ import print_function  # for debug
import numpy as np
import datetime as dt
from astropy.time import Time
from sgp4.earth_gravity import wgs84
from sgp4.io import twoline2rv
import time
```

Using ISS Camera Geolocate
--------------------------
To access everything:
```
import iss_camera_geolocate as icg
```
A demonstration notebook is in the notebooks directory.
