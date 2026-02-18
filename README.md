# [HEXAPIC](https://github.com/lecad-peg/hexapic) #

[**H**eterogenous **EXA**scale **P**article-**I**n-**C**ell](https://github.com/lecad-peg/hexapic)

OpenMPI version with parallel I/O through openPMD with ADIOS2 backend.

Feature: 3D visualisation with pyqtgraph.

### Dependency installation ###

See source/compile_HPCFS.sh and source/compile_VEGA.sh for the required modules before the installation of libraries on different HPCs.

#### ADIOS2 ####

```bash
cd ~/
git clone https://github.com/ornladios/ADIOS2.git ADIOS2
cd ADIOS2
mkdir adios2-build
cd adios2-build
cmake ../ -DADIOS2_USE_MPI=ON -DADIOS2_USE_HDF5=OFF -DADIOS2_USE_Python=ON -DADIOS2_USE_SST=ON -DADIOS2_USE_BZip2=ON -DCMAKE_INSTALL_PREFIX=~/adios2/ # optional: -DBUILD_TESTING=ON
make -j8
# optional: ctest
make install
```

#### openPMD ####
```bash
cd ~/
git clone https://github.com/openPMD/openPMD-api.git
cd openPMD-api
mkdir openPMD-api-build
cd openPMD-api-build
cmake ../ -DopenPMD_USE_MPI=ON -DopenPMD_USE_HDF5=OFF -DopenPMD_USE_ADIOS2=ON -DopenPMD_USE_PYTHON=ON -DCMAKE_INSTALL_PREFIX=~/openpmd/ -DADIOS2_DIR=~/adios2/lib/cmake/adios2/ # optional: -DPYTHON_EXECUTABLE=/software/Python/3.11.3-GCCcore-12.3.0/bin/python -DPYTHON_INCLUDE_DIRS=/software/Python/3.11.3-GCCcore-12.3.0/include/python3.11
cmake --build .
# optional: ctest
cmake --build . --target install
```

#### HYPRE ####
```bash
cd ~/
git clone https://github.com/hypre-space/hypre.git
cd hypre/src
./configure --prefix=$HOME/hypre/install --with-MPI
make -j$(nproc)
# optional: ctest
make install
```


#### Python modules ####

```bash
pip install numpy
pip install scipy
pip install pyqtgraph
pip install picmistandard
```


### Compilation ###

#### HPCFS ####

```bash
cd source
source compile_HPCFS.sh
cd ../
```
#### VEGA ####

```bash
cd source
source compile_VEGA.sh
cd ../
```


### Pre-commit ###

We use [pre-commit](https://pre-commit.com/) to enforce trailing-whitespace and end-of-file fixes, YAML checks, and clang-format for C/C++ in the repo. Install pre-commit (e.g. `pip install pre-commit`), then from the repo root run `pre-commit install`. Hooks run automatically on `git commit`; you can run them on all files anytime with `pre-commit run --all-files`.

### Usage ###

Run the executable:

```bash
# optional: cp picmi/input_file.inp test/
mpiexec -n 4 source/HEXAPIC test/input_file.inp -steps Nsteps -dsteps Nskip
```
* If `-steps Nsteps` is not provided, the simulation will run for 100 steps.
* Create and save  diagnostics output every `Nskip` steps. Default is `1`.

### Visualisation ###

#### Custom Python script

```bash
cd test/
python ../vis/vis3D-class.py input_file.inp.bp4
```

#### ParaView (version >= 5.13)

1. Open ParaView
2. Open > select `.bp4` output.
3. On dialog *Open Data With...* select `ADIOS2CoreImageReader` and press `OK`.
4. On the lfet side, under *Proeperties* where arrays are located, select the arrays to plot. Do not mix cell-centered meshes (n, T, v) with node-centered meshes (V).
5. Under *Properties*, on *Image Dimension* select one of your selected arrays.
6. On *Time step array* select `/data/meshes/time`.
7. Press `Apply` and wait for the data to load.
8. Change `Solid Color` to any of the desired loaded arrays.
9. Change `Outline` to `Surface`.
10. Optional: press the *play* button to see the time evolution of the simulation.

### Documentation ###

Documentation is available at: https://lecad-peg.github.io/hexapic/
