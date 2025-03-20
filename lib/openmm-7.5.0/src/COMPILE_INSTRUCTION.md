A compile instruction for the QCDyn modified openmm-7.5.0

Modified By Zhubin Hu @ Nov 2020. 
Original source directory : 
/xspace/sungroup/software/source/openmm-7.5.0_mod

To compile:

Set the correct directories for the source and everything
 
# specify the path of source code and installation
SRC=$path_to_source/src
INSTALL=$path_to_source/src
FFTW=$path_to_source/fftw_3.3.8_openmm

# Swig and Doxygen for Python API (optional, -DOPENMM_BUILD_PYTHON_WRAPPERS=ON)
# But for the modified source code for QCDyn, doesn't support Python API 
SWIG=$path_to_source/swig_4.0.2
DOXY=$path_to_source/doxygen_1.9.0
GCC=$compiler_directory/packages/gcc/7.3/bin/gcc
GXX=$compiler_directory/packages/gcc/7.3/bin/g++

Run the following script 
./openmm_build_modified.sh

Example script is given in this directory too. 


