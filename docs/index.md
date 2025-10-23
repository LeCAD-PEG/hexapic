# Overview

**H**&#8203;eterogenous **EXA**&#8203;scale **P**&#8203;article-&#8203;**I**&#8203;n-&#8203;**C**&#8203;ell

- **Language**: C++

- **Dependencies**: 
	- MPI
	- openPMD with ADIOS2 backend and Python support
	- PICMI (Python module)
	- pyqtgraph (Python module)

- **Target**: Fusion plasma, Scrape-Off-Layer of tokamaks.

- **Features**:
	- Full-orbit, Debye-sheath-resolving electrostatic PIC
	- Multi-node MPI application
	- Domain decomposition
	- XOOPIC-like input file
	- Particle and heat sources: planar, volumetric
	- Multi-species
	- Monte-Carlo Collisions: elastic, excitation, ionization, charge-exchange
	- Plasma-surface interactions:
		- ion-recycling (into neutrals)
		- secondary-electron emmission
		- particle-impact erosion and impurity injection
	- 3D visualisation

<p align="center">
  <img src="images/V.png" alt="space potential" width="28%">
  <img src="images/n.png" alt="electron density" width="27%">
  <img src="images/T.png" alt="electron temperature" width="30%">
</p>

<p align="center">
  &nbsp; Space potential [V] &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Electron density [m<sup>−3</sup>] &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp; Electron temperature [eV]
</p>


- **Workflow**:

```plantuml
@startuml
start
:init; <<input>>
:initial load; <<load>>
repeat :inject particles; <<task>>
:charge to grid; <<task>>
:compute field; <<task>>
:move particles; <<task>>
:check boundaries; <<task>>
:collisions; <<task>>
split
	if (output?) then (yes)
		:diagnostics; <<output>>
	else (no)
	endif
split again
	if (checkpoint?) then (yes)
		:save state; <<save>>
	else (no)
	endif 
end split
repeat while (particles and/or steps?) is (Yes) not (No)
:free; <<task>>
stop
@enduml
```
