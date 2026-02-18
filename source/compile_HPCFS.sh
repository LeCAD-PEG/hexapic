#!/bin/bash

module purge
module load OpenMPI/5.0.8-GCC-14.3.0
module load OpenBLAS/0.3.30-GCC-14.3.0
module load FFmpeg/7.1.2-GCCcore-14.3.0

export OMP_NUM_THREADS=1

# OpenPMD
export LD_LIBRARY_PATH=$HOME/openpmd/lib64:$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=$HOME/openpmd/lib64/pkgconfig:$PKG_CONFIG_PATH
# change path to your python MAJOR.MINOR version
export PYTHONPATH=$HOME/openpmd/lib64/python3.13/site-packages:$PYTHONPATH


cp makefile_HPCFS makefile
make clean
make
