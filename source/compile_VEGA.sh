#!/bin/bash

module --force purge
module load OpenMPI/4.1.6-GCC-13.2.0
module load Python/3.11.5-GCCcore-13.2.0

export OMP_NUM_THREADS=1

# OpenPMD
export LD_LIBRARY_PATH=$HOME/openpmd/lib64:$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=$HOME/openpmd/lib64/pkgconfig:$PKG_CONFIG_PATH
# change path to your python MAJOR.MINOR version
export PYTHONPATH=$HOME/openpmd/lib64/python3.11/site-packages:$PYTHONPATH

# HYPRE install directory
export LD_LIBRARY_PATH="$HOME/hypre/install/lib:${LD_LIBRARY_PATH}"


cp makefile_VEGA makefile
make clean
make
