# Overview

**H**eterogenous **EXA**scale **P**article-**I**n-**C**ell

- **Language**: C++

- **Dependencies**: 
	- MPI
	- openPMD with ADIOS2 backend and Python support
	- PICMI (Python module)
	- pyqtgraph (Python module)

- **Target**: Fusion plasma, Scrape-Off-Layer of tokamaks.

- **Features**:
	- Full-orbit, sheath-resolving electrostatic PIC
	- Multi-node
	- Domain decomposition
	- XOOPIC-like input file
	- Particle and heat sources: planar, volumetric
	- Multi-species
	- Monte-Carlo Collisions
	- Plasma-surface interactions:
		- ion-recycling (into neutrals)
		- secondary-electron emmission
		- particle-impact erosion and impurity injection
	- 3D visualisation
