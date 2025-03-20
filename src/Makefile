# ==============================================================================
# Makefile for Quantum Classical Dynamics (QCDyn).
# verison: 0.3
# Last updated: Nov. 13, 2021
# Author: Zhubin Hu @ Sun-Group @ NYU-SH
# ==============================================================================
# Before building this program from source, you will need the following:
# 1. C++ compiler (gcc >= 7.3)
# 2. CUDA library (11.2) (if OpenMM with GPU is reuired)
#
# The following libraries have been provided in QCDyn/lib/:
# 1. OpenMM 7.5.0 library (modified, dynamics shared lib, precomplied by
#	 gcc 7.3, cuda 11.2, and FFTW 3.3.8)
# 2. Eigen 3.3.9 library (template library for linear algebra)
# 3. chemfiles 0.10.2 library (static library to read/write both text (XYZ, PDB, ...)
#    and binary (NetCDF, TRR, ...) file formats, precomplied by gcc 7.3)
# 4. FFTw3.3.8 // by xiaofang
# If you build and install it in our HPC, you can load them by adding the
# following sentences to your ~/.bashrc: (The QCDyn is in home directory)
# module load gcc/7.3 cuda/11.2
# export QCDYN=~/QCDyn
# export PATH=$PATH:$QCDYN/bin
# export OPENMM_INSTALL_DIR=$QCDYN/lib/openmm-7.5.0
# export OPENMM_PLUGIN_DIR=$OPENMM_INSTALL_DIR/lib/plugins
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OPENMM_INSTALL_DIR/lib
# export OPENMM_CUDA_COMPILER=/gpfsnyu/packages/cuda/11.2/bin/nvcc
# ------------------------------------------------------------------------------

# Check the install directory of OpenMM library
OpenMM_INSTALL_DIR=../lib/openmm-7.5.0
LIB_DIR=$(OpenMM_INSTALL_DIR)/lib
INCLUDE_DIR=$(OpenMM_INSTALL_DIR)/include
LIBS= -lOpenMM

# Check the install directory of Eigen library
Eigen_INC_DIR=../lib/eigen-3.3.9

# Check the install directory of chemfiles library
chemfiles_INC_DIR=../lib/chemfiles-0.10.2/include
chemfiles_LIBS=../lib/chemfiles-0.10.2/lib/libchemfiles.a

# Check the install directory of fftw3 library
fftw_INC_DIR=../lib/fftw_3.3.8_double_intel_avx2/include
fftw__LIBS=../lib/fftw_3.3.8_double_intel_avx2/lib/libfftw3.a
# g++ flags
#CPPFLAGS = -std=c++11 -O3

# Define build type: Debug or Release
BUILD_TYPE ?= Release
 # Compiler flags for different build types
 # use this commond:  make BUILD_TYPE=Debug or make BUILD_TYPE=Release to choose the build type
ifeq ($(BUILD_TYPE), Debug)
     CPPFLAGS = -std=c++11 -g -O0
else ifeq ($(BUILD_TYPE), Release)
     CPPFLAGS = -std=c++11 -O3
else
     $(error Unknown BUILD_TYPE: $(BUILD_TYPE))
endif

# The name of program
PROG=QCDyn

# The location of program
PROG_DIR=../../$(PROG)

# The location of source files (.cpp)
SRC_DIR=$(PROG_DIR)/src

# The location of header files (.h)
INC_DIR=$(PROG_DIR)/src

# The location of binaray file
BIN_DIR=$(PROG_DIR)/bin

# All source files (.cpp)
ALL_CPP = $(SRC_DIR)/$(wildcard *.cpp)

# All header files (.h)
ALL_H = $(INC_DIR)/$(wildcard *.h)

# All object files (.o)
ALL_O = $(ALL_H:.h=.o)

# Build main program and link
# This will build object files firstly if they don't exist
all : $(PROG).cpp $(ALL_O)
	g++ $(CPPFLAGS) -I$(INCLUDE_DIR) -L$(LIB_DIR) $(LIBS) -I$(Eigen_INC_DIR) -I$(chemfiles_INC_DIR) -I$(fftw_INC_DIR) -o $(BIN_DIR)/$(PROG) $(PROG).cpp $(ALL_O) $(chemfiles_LIBS) $(fftw__LIBS)

# Build object files (Treat all .cpp source files (except main.cpp) the same way)
%.o : %.cpp %.h
	g++ $(CPPFLAGS) -I$(INCLUDE_DIR) -L$(LIB_DIR) $(LIBS) -I$(Eigen_INC_DIR) -I$(chemfiles_INC_DIR) -I$(fftw_INC_DIR) -c $<

# Build main program from all source and header files directly
QCDyn : $(ALL_CPP) $(ALL_H)
	g++ $(CPPFLAGS) -I$(INCLUDE_DIR) -L$(LIB_DIR) $(LIBS) -I$(Eigen_INC_DIR) -I$(chemfiles_INC_DIR) -I$(fftw_INC_DIR) -o $(BIN_DIR)/$(PROG) $(ALL_CPP) $(chemfiles_LIBS) $(fftw__LIBS)

# Delete binary file and object files
clean :
	rm $(ALL_O) $(BIN_DIR)/$(PROG)

