from setuptools import setup

setup(
    name='ISS_Camera_Geolocate',
    version='0.1',
    author='Timothy Lang',
    author_email='timothy.j.lang@nasa.gov',
    packages=['iss_camera_geolocate', ],
    license='LICENSE.md',
    description='ISS Camera Geolocation Library',
    long_description=open('description.txt').read(),
    install_requires=['numpy', 'sgp4', 'astropy', 'datetime'],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Environment :: Console"
        ],
)
