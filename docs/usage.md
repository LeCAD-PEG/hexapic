# Usage

## Build executable

1. Navigate to `source/` directory.

2. Source (execute in the current shell):

	```bash
	# OpenPMD install directory
	export LD_LIBRARY_PATH=$HOME/openpmd/lib:$LD_LIBRARY_PATH
	export PKG_CONFIG_PATH=$HOME/openpmd/lib/pkgconfig:$PKG_CONFIG_PATH
	# change path to your python MAJOR.MINOR version
	export PYTHONPATH=$HOME/openpmd/lib/python3.13/site-packages:$PYTHONPATH

	# HYPRE install directory
	export LD_LIBRARY_PATH="$HOME/hypre/install/lib:${LD_LIBRARY_PATH}"

	make clean # optional
	make
	```

## Create HEXAPIC input file

1. Navigate to `picmi` directory.

2. Modify `create_input_file.py` according to your simulation needs.

3. Run:
	```bash
	python create_input_file.py`
	``` 
	to create `input_file.inp`.


## Run the simulation

1. Copy the previously created input file into the directory where you want to run the simulation, e.g. `test/`.

2. Run the HEXAPIC executable:
	```bash
	mpiexec --bind-to core -n 4 source/HEXAPIC test/input_file.inp -steps Nsteps
	```
	**_NOTE_**: If `-steps Nsteps` is not provided, the simulation will run for **100** steps.

## Visualise code output

1. Install required Python modules: *numpy*, *scipy*, *pyqtgraph*, *pyside6*, *opengl*
 
2. Run:

	```bash
	python vis/vis3D-class.py test/input_file.inp.bp4
	```
	This will bring up the GUI:

	<p align="center">
	  <img src="../images/gui.png" width="90%">
	</p>

3. Click on the fields to the left to switch between the physical quantities to visualise.

	GUI field description:
	
	- n&#8203;*x* : Particle density of specie *x* in *m<sup>−3</sup>*
	- T&#8203;*x*: Temperature of specie *x* in *eV*
	- v&#8203;*x*: Velocity of specie *x* in *m/s*
	- V   : Space potential in V