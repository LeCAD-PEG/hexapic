#!/bin/bash

module purge
module load OpenMPI/4.1.6-GCC-13.2.0 
module load OpenBLAS/0.3.24-GCC-13.2.0
module load FFmpeg/6.0-GCCcore-13.2.0

export OMP_NUM_THREADS=1

# OpenPMD
export LD_LIBRARY_PATH=$HOME/openpmd/lib:$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=$HOME/openpmd/lib/pkgconfig:$PKG_CONFIG_PATH
# change path to your python MAJOR.MINOR version
export PYTHONPATH=$HOME/openpmd/lib/python3.13/site-packages:$PYTHONPATH


cp makefile_HPCFS makefile
make clean
make
