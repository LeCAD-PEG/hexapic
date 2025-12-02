# CHANGELOG


## 0.8.0 - 14.11.2025

### Added
* Persistent MPI

### Changed
* Formulas for secondary electron emission

### Fixed
* Reserve memory to prevent reallocation
* Increase MPI communication buffer


## 0.7.0 - 16.10.2025

### Changed

* PETSc solver changed to HYPRE PFMG solver 


## 0.6.0 - 27.02.2025

### Added

* Monte-Carlo Collisions: elastic, excitation, ionisation, charge-exchange.
* Plasma-Surface Interactions: electron-induced and ion-induced secondary electron emission, ion recycling into neutrals.


## 0.5.0 - 17.12.2024

### Added

* Domain decomposition using recursive bisection method
* Diagnostics: density, temperature and fluid velocity on cells, space potential on grid nondes.
* In-situ visualisation with pyqtgraph


## 0.4.0 - 14.10.2024

### Added

* Parallel I/O with openPMD and ADIOS2 backend
* Particle pusher subcycling
* Periodic and Neumann boundary conditions


## 0.3.0 - 16.09.2024

### Added

* **Parallel** version using MPI

### Fixed

* Parallel update of the field source term
* Reduce calls to PETSc vector assembly

### Changed

* Field solver changed from Algebraic Multigrid to Biconjugate Gradient with Geometric Multigrid preconditioner
* Type for indexing changed from `int` to `size_t` for large grids.


## 0.2.0 - 20.08.2024

### Added

* PICMI standard
* XOOPIC-compatible input file
* Dirichlet boundaries with external electric potential
* Uniform per-specie particle load in a custom region
* Volumetric particle source based on source temperature
* Plane particle injection source at boundaries

### Changed

* Field solver from Jacobi to Algebraic Multigrig


## 0.1.0 - 29.07.2024

### Added 

* **Serial** 2D PIC routine
* Uniform initial particle load
* Boris particle pusher
* PETSc field solver (Jacobi)
* Absorbing boundaries
