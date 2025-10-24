

# File plasma\_wall.cpp



[**FileList**](files.md) **>** [**hexapic**](dir_b8175b39202167732592fd9c59b14660.md) **>** [**source**](dir_4d88c249dac475c11ca7037572002e86.md) **>** [**plasma\_wall.cpp**](plasma__wall_8cpp.md)



[_**Boundary**_](structBoundary.md) _conditions and plasma–wall interaction (sheath models, secondary emission, reflection)._[More...](#detailed-description)

* `#include "hexapic.hpp"`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**PSI\_init**](#function-psi_init) () <br> |
|  bool | [**boundary\_in\_rank**](#function-boundary_in_rank) (int bound\_i) <br> |
|  PetscScalar | [**gamma**](#function-gamma) (PetscScalar E, PetscScalar sigmd, PetscScalar gamma0, PetscScalar E0, PetscScalar Ep, int sec\_type) <br> |
|  void | [**inject\_secondary**](#function-inject_secondary) (int sp1, int sp2, int inj\_type, PetscScalar & E, PetscScalar T0, int rank, int dir, PetscScalar xy\_along, PetscScalar z) <br> |
|  void | [**particle\_wall\_interaction**](#function-particle_wall_interaction) (int cell\_i, int sp, int p\_i, int lrbt) <br> |
|  void | [**random\_order**](#function-random_order) (std::vector&lt; int &gt; & order, int n) <br> |
|  void | [**reflect\_particle**](#function-reflect_particle) (int cell\_i, int sp, int p\_i, int rank, int dir) <br> |
|  void | [**secondary\_emission**](#function-secondary_emission) (PetscScalar & E, PetscScalar sigmd, int sp, int sec\_i, int rank, int dir, PetscScalar xy\_along, PetscScalar z) <br> |
|  void | [**sigmd\_E**](#function-sigmd_e) (PetscScalar & sigmd, PetscScalar & E, int cell\_i, int sp, int p\_i, int dir) <br> |
|  void | [**wall\_interaction**](#function-wall_interaction) (int cell\_i, int sp, int p\_i, int rank, int dir) <br> |




























## Detailed Description




**Copyright:**

Copyright (c) 2025 




    
## Public Functions Documentation




### function PSI\_init 

```C++
void PSI_init () 
```




<hr>



### function boundary\_in\_rank 

```C++
bool boundary_in_rank (
    int bound_i
) 
```




<hr>



### function gamma 

```C++
PetscScalar gamma (
    PetscScalar E,
    PetscScalar sigmd,
    PetscScalar gamma0,
    PetscScalar E0,
    PetscScalar Ep,
    int sec_type
) 
```




<hr>



### function inject\_secondary 

```C++
void inject_secondary (
    int sp1,
    int sp2,
    int inj_type,
    PetscScalar & E,
    PetscScalar T0,
    int rank,
    int dir,
    PetscScalar xy_along,
    PetscScalar z
) 
```




<hr>



### function particle\_wall\_interaction 

```C++
void particle_wall_interaction (
    int cell_i,
    int sp,
    int p_i,
    int lrbt
) 
```




<hr>



### function random\_order 

```C++
void random_order (
    std::vector< int > & order,
    int n
) 
```




<hr>



### function reflect\_particle 

```C++
void reflect_particle (
    int cell_i,
    int sp,
    int p_i,
    int rank,
    int dir
) 
```




<hr>



### function secondary\_emission 

```C++
void secondary_emission (
    PetscScalar & E,
    PetscScalar sigmd,
    int sp,
    int sec_i,
    int rank,
    int dir,
    PetscScalar xy_along,
    PetscScalar z
) 
```




<hr>



### function sigmd\_E 

```C++
void sigmd_E (
    PetscScalar & sigmd,
    PetscScalar & E,
    int cell_i,
    int sp,
    int p_i,
    int dir
) 
```




<hr>



### function wall\_interaction 

```C++
void wall_interaction (
    int cell_i,
    int sp,
    int p_i,
    int rank,
    int dir
) 
```




<hr>

------------------------------


