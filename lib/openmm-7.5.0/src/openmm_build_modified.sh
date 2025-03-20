#!/bin/bash
## build OpenMM lib from source code with gcc 7.3 and cuda 11.2

# load modules
module load cmake/3.15.0 gcc/7.3 cuda/11.2

# build FFTW lib firstly with gcc/7.3 (done)
#wget ftp://ftp.fftw.org/pub/fftw/fftw-3.3.8.tar.gz
#tar xvf fftw-3.3.8.tar.gz
#cd fftw-3.3.8
#./configure --prefix=/xspace/sungroup/software/fftw_3.3.8 --enable-sse2 --enable-avx --enable-avx2 --enable-float --enable-shared --enable-threads
#make install -j 8

# specify the path of source code and installation
SRC=/xspace/sungroup/software/source/openmm-7.5.0_mod
INSTALL=/xspace/sungroup/software/openmm_7.5.0_mod
FFTW=/xspace/sungroup/software/fftw_3.3.8_openmm
# Swig and Doxygen for Python API (optional, -DOPENMM_BUILD_PYTHON_WRAPPERS=ON)
# But for the modified source code for QCDyn, doesn't support Python API 
SWIG=/xspace/sungroup/software/swig_4.0.2
DOXY=/xspace/sungroup/software/doxygen_1.9.0
GCC=/gpfsnyu/packages/gcc/7.3/bin/gcc
GXX=/gpfsnyu/packages/gcc/7.3/bin/g++
BUILD=openmm_build_mod
NPROCES=24

# configure, build, and install OpenMM
rm -rf $INSTALL
rm -rf $BUILD
mkdir $BUILD
cd $BUILD
export CC=$GCC
export CXX=$GXX
export CMAKE_PREFIX_PATH=$FFTW:${SWIG}:${DOXY}
cmake $SRC -DCMAKE_INSTALL_PREFIX=$INSTALL -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF \
        -DOPENMM_BUILD_AMOEBA_CUDA_LIB=ON \
        -DOPENMM_BUILD_AMOEBA_PLUGIN=ON \
        -DOPENMM_BUILD_CPU_LIB=ON \
        -DOPENMM_BUILD_CUDA_COMPILER_PLUGIN=ON \
        -DOPENMM_BUILD_CUDA_DOUBLE_PRECISION_TESTS=OFF \
        -DOPENMM_BUILD_CUDA_LIB=ON \
        -DOPENMM_BUILD_CUDA_TESTS=OFF \
        -DOPENMM_BUILD_C_AND_FORTRAN_WRAPPERS=OFF \
        -DOPENMM_BUILD_DRUDE_CUDA_LIB=ON \
        -DOPENMM_BUILD_DRUDE_OPENCL_LIB=ON \
        -DOPENMM_BUILD_DRUDE_PLUGIN=ON \
        -DOPENMM_BUILD_OPENCL_DOUBLE_PRECISION_TESTS=OFF \
        -DOPENMM_BUILD_OPENCL_LIB=ON \
        -DOPENMM_BUILD_OPENCL_TESTS=OFF \
        -DOPENMM_BUILD_PME_PLUGIN=ON \
        -DOPENMM_BUILD_PYTHON_WRAPPERS=OFF \
        -DOPENMM_BUILD_RPMD_CUDA_LIB=ON \
        -DOPENMM_BUILD_RPMD_OPENCL_LIB=ON \
        -DOPENMM_BUILD_RPMD_PLUGIN=ON \
        -DOPENMM_BUILD_SHARED_LIB=ON \
        -DOPENMM_BUILD_STATIC_LIB=OFF
make -j $NPROCES
make install
# install python API (optional, -DOPENMM_BUILD_PYTHON_WRAPPERS=ON)
#make PythonInstall
