# GRASP 

GRASP is a highly accurate aerosol retrieval algorithm that processes properties 
of aerosol- and land-surface-reflectance in cloud free environments. It infers 
nearly 50 aerosol and surface parameters including particle size distribution, 
the spectral index of refraction, the degree of sphericity and absorption. The 
algorithm is designed for the enhanced characterization of aerosol properties 
from spectral, multiangular polarimetric remote sensing observations. GRASP works 
under different conditions, including bright surfaces like deserts, where the 
reflectance overwhelms the signal of aerosols. GRASP is highly versatile and allows
 input from a wide variety of satellite and surface measurements.

For more information about GRASP you can search official project webpage: http://www.grasp-open.com/

For user documentation please, go to: http://www.grasp-open.com/doc

For technical documentation, go to: http://www.grasp-open.com/tech-doc

The code is released under GRASP Open Source License V1.0 
Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
See LICENSE file for more information. 
For more information on commercial use or interest in a specific research project, please contact office@grasp-sas.com


# How to build

There are many different ways to get the code compiled. GRASP use cmake 
(http://www.cmake.org/) to compile the code but it is wrapper by a classical make
 script which simplified the procedure.  Please, check before continue that your 
environment fulfils all requirements. More information in "Build dependencies" section.

## Using make wrapper

```sh
make # build the project using the default build settings
sudo make install # install grasp
grasp # test the command
```


## Advanced mode: Using directly cmake

```sh
mkdir ~/grasp/build
cd ~/grasp/build
cmake .. -DCMAKE_BUILD_TYPE=Release -DADDITIONAL_DEPENDENCIES_PATH=/usr/local/grasp-deps -DCONSTANTS_SET=generic
make -j12
sudo make install
grasp # test the command
```


## Building GPGPU retrieval

GPGPU module now is an external extension. Use grasp manager to install it and then
to run the GPGPU version of GRASP the following CMake flags are required:
```sh
-DENABLE_MODULES=gpgpu -DENABLE_GPGPU=ON
```


## Build dependencies

### Ubuntu:

Required:
```sh
sudo apt-get install build-essential cmake git libsuperlu4 libsuperlu-dev \
                     gfortran libyaml-dev libglib2.0-dev libcunit1-dev 
```

User has to take into account that there are optional dependencies depending in installed extensions.


# GPGPU version

The GPGPU code is currently at the state of the v0.3.1 tag.
