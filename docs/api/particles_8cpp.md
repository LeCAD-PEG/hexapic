

# File particles.cpp



[**FileList**](files.md) **>** [**hexapic**](dir_b8175b39202167732592fd9c59b14660.md) **>** [**source**](dir_4d88c249dac475c11ca7037572002e86.md) **>** [**particles.cpp**](particles_8cpp.md)



[_**Particle**_](structParticle.md) _pusher, weighting/interpolation, and particle management (injection, movement, deposition)._[More...](#detailed-description)

* `#include "hexapic.hpp"`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**add\_received\_particles**](#function-add_received_particles) () <br> |
|  void | [**cell2part**](#function-cell2part) (PetscScalar px, PetscScalar py, PetscScalar \* V\_cell, PetscScalar \* Ef) <br> |
|  void | [**cells\_init**](#function-cells_init) () <br> |
|  void | [**grid2cell**](#function-grid2cell) (PetscInt J, PetscScalar \* V, PetscScalar \* V\_cell) <br> |
|  void | [**grid2part**](#function-grid2part) (PetscScalar Px, PetscScalar Py, PetscScalar \* V, PetscScalar \* Ef) <br> |
|  void | [**grid\_init**](#function-grid_init) (int argc, char \*\* args) <br> |
|  void | [**initial\_particle\_load**](#function-initial_particle_load) () <br> |
|  void | [**inject\_particle**](#function-inject_particle) (int sp, PetscScalar x, PetscScalar y) <br> |
|  void | [**inject\_particles**](#function-inject_particles) () <br> |
|  PetscScalar | [**normvel\_MaxBol**](#function-normvel_maxbol) () <br> |
|  PetscScalar | [**normvel\_vFv**](#function-normvel_vfv) () <br> |
|  void | [**num\_param\_init**](#function-num_param_init) () <br> |
|  void | [**part2grid**](#function-part2grid) () <br> |
|  void | [**particle\_boundaries**](#function-particle_boundaries) () <br> |
|  void | [**particle\_boundaries\_dd**](#function-particle_boundaries_dd) () <br> |
|  void | [**particle\_mover\_boris**](#function-particle_mover_boris) () <br> |
|  void | [**particles\_init**](#function-particles_init) () <br> |
|  void | [**remove\_particle**](#function-remove_particle) (int cell\_i, int sp, int p) <br> |
|  void | [**send\_receive\_particles**](#function-send_receive_particles) () <br> |




























## Detailed Description




**Copyright:**

Copyright (c) 2025 




    
## Public Functions Documentation




### function add\_received\_particles 

```C++
void add_received_particles () 
```




<hr>



### function cell2part 

```C++
void cell2part (
    PetscScalar px,
    PetscScalar py,
    PetscScalar * V_cell,
    PetscScalar * Ef
) 
```




<hr>



### function cells\_init 

```C++
void cells_init () 
```




<hr>



### function grid2cell 

```C++
void grid2cell (
    PetscInt J,
    PetscScalar * V,
    PetscScalar * V_cell
) 
```




<hr>



### function grid2part 

```C++
void grid2part (
    PetscScalar Px,
    PetscScalar Py,
    PetscScalar * V,
    PetscScalar * Ef
) 
```




<hr>



### function grid\_init 

```C++
void grid_init (
    int argc,
    char ** args
) 
```




<hr>



### function initial\_particle\_load 

```C++
void initial_particle_load () 
```




<hr>



### function inject\_particle 

```C++
void inject_particle (
    int sp,
    PetscScalar x,
    PetscScalar y
) 
```




<hr>



### function inject\_particles 

```C++
void inject_particles () 
```




<hr>



### function normvel\_MaxBol 

```C++
PetscScalar normvel_MaxBol () 
```




<hr>



### function normvel\_vFv 

```C++
PetscScalar normvel_vFv () 
```




<hr>



### function num\_param\_init 

```C++
void num_param_init () 
```




<hr>



### function part2grid 

```C++
void part2grid () 
```




<hr>



### function particle\_boundaries 

```C++
void particle_boundaries () 
```




<hr>



### function particle\_boundaries\_dd 

```C++
void particle_boundaries_dd () 
```




<hr>



### function particle\_mover\_boris 

```C++
void particle_mover_boris () 
```




<hr>



### function particles\_init 

```C++
void particles_init () 
```




<hr>



### function remove\_particle 

```C++
void remove_particle (
    int cell_i,
    int sp,
    int p
) 
```




<hr>



### function send\_receive\_particles 

```C++
void send_receive_particles () 
```




<hr>

------------------------------


