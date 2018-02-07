# Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
# Licensed under the GRASP Open Source License V1.0 (see LICENSE file)

CONSTANTS_SET ?= generic
BUILD ?= Release
MPI ?= off
DEBUG_MPI ?= off
CC ?= cc
CXX ?= c++
F90 ?= gfortran
PREFIX ?= /usr/local
BUILD_DIR ?= build
CMAKE_FLAGS ?= 
BUILD_FLAGS ?= -s -j12

##### CONFIG OPTIONS #######
# CONSTANTS_SET ?= generic satellite lidar polder aeronet
# BUILD ?= Debug Release RelWithDebInfo Fast
# F90 ?= gfortran ifort
############################

.PHONY: all config mpi install dep retrieval clean

all: config
	cd $(BUILD_DIR) && make $(BUILD_FLAGS)

config: $(BUILD_DIR)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR) && cd $(BUILD_DIR) && cmake -DCMAKE_C_COMPILER=$(CC) -DCMAKE_CXX_COMPILER=$(CXX) -DCMAKE_Fortran_COMPILER=$(F90) -DADDITIONAL_DEPENDENCIES_PATH=$(PREFIX)/grasp-deps -DCMAKE_INSTALL_PREFIX=$(PREFIX) -DCONSTANTS_SET=$(CONSTANTS_SET) -DCMAKE_BUILD_TYPE=$(BUILD) -DUSE_MPI=$(MPI) -DDEBUG_MPI=$(DEBUG_MPI) -DENABLE_RETRIEVAL_MODULES=$(ENABLE_RETRIEVAL_MODULES) $(CMAKE_FLAGS) ..

mpi: 
	make $(BUILD_DIR) MPI=on && make $(BUILD_FLAGS)
	

install: $(BUILD_DIR)
	cd $(BUILD_DIR) && make PREFIX=$(PREFIX) install
	
dep:
	if test -d dependencies; then cd dependencies && make PREFIX=$(PREFIX) install; else echo "Dependencies module is not available. Install it"; fi

dep-polder:
	if test -d dependencies; then cd dependencies && make PREFIX=$(PREFIX) with-polder; else echo "Dependencies module is not available. Install it"; fi
	
retrieval:
	mkdir -p $(BUILD_DIR) && cd $(BUILD_DIR) && cmake -DCMAKE_Fortran_COMPILER=$(F90) -DADDITIONAL_DEPENDENCIES_PATH=$(PREFIX)/grasp-deps -DCMAKE_INSTALL_PREFIX=$(PREFIX) -DCONSTANTS_SET=$(CONSTANTS_SET) -DCMAKE_BUILD_TYPE=$(BUILD) -DENABLE_RETRIEVAL_MODULES=$(ENABLE_RETRIEVAL_MODULES) $(CMAKE_FLAGS) .. && make  grasp_retrieval $(BUILD_FLAGS)

clean:
	-rm -rf $(BUILD_DIR)

