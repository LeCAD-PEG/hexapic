

# File collisions.cpp



[**FileList**](files.md) **>** [**hexapic**](dir_b8175b39202167732592fd9c59b14660.md) **>** [**source**](dir_4d88c249dac475c11ca7037572002e86.md) **>** [**collisions.cpp**](collisions_8cpp.md)



_Monte-Carlo Collisions (MCC): elastic, ionization, excitation, and charge-exchange processes._ [More...](#detailed-description)

* `#include "hexapic.hpp"`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**COM\_velocity**](#function-com_velocity) (PetscScalar v\_com, int cell\_i, int sp\_1, int sp\_2, int p\_1, int p\_2) <br> |
|  void | [**MCC**](#function-mcc) () <br> |
|  void | [**add\_new\_empty\_collision**](#function-add_new_empty_collision) (int sp1, int sp2) <br> |
|  void | [**add\_sp\_out**](#function-add_sp_out) (int coll\_i, int reaction\_index) <br> |
|  void | [**charge\_exchange\_collision**](#function-charge_exchange_collision) (int coll\_i, int react\_i, int cell\_i, int sp\_1, int sp\_2, int p\_1, int p\_2, PetscScalar E) <br> |
|  void | [**collision\_in\_cell**](#function-collision_in_cell) (int coll\_i, int cell\_i) <br> |
|  void | [**collisions\_init**](#function-collisions_init) () <br> |
|  PetscScalar | [**density\_computed\_particles**](#function-density_computed_particles) (int cell\_i, int sp) <br> |
|  void | [**elastic\_scattering\_collision**](#function-elastic_scattering_collision) (int coll\_i, int react\_i, int cell\_i, int sp\_1, int sp\_2, int p\_1, int p\_2, PetscScalar E) <br> |
|  void | [**excitation\_collision**](#function-excitation_collision) (int coll\_i, int react\_i, int cell\_i, int sp\_1, int sp\_2, int p\_1, int p\_2, PetscScalar E) <br> |
|  reactionFunction | [**get\_reaction**](#function-get_reaction) (int reaction\_index) <br> |
|  void | [**ionisation\_collision**](#function-ionisation_collision) (int coll\_i, int react\_i, int cell\_i, int sp\_1, int sp\_2, int p\_1, int p\_2, PetscScalar E) <br> |
|  void | [**isotropic\_scatter**](#function-isotropic_scatter) (PetscScalar v\_scatt) <br> |
|  void | [**isotropic\_scatter**](#function-isotropic_scatter) (PetscScalar v\_scatt, PetscScalar v) <br> |
|  void | [**null\_collision**](#function-null_collision) (int coll\_i, int react\_i, int cell\_i, int sp\_1, int sp\_2, int p\_1, int p\_2, PetscScalar E) <br> |
|  void | [**read\_add\_collision**](#function-read_add_collision) (int i) <br> |
|  void | [**read\_add\_cross\_sections**](#function-read_add_cross_sections) (int coll\_i, std::string filename) <br> |
|  PetscScalar | [**read\_cross\_section**](#function-read_cross_section) (int coll\_i, int j, PetscScalar E) <br> |
|  void | [**recombination\_collision**](#function-recombination_collision) (int coll\_i, int react\_i, int cell\_i, int sp\_1, int sp\_2, int p\_1, int p\_2, PetscScalar E) <br> |
|  PetscScalar | [**reduced\_mass**](#function-reduced_mass) (PetscScalar m1, PetscScalar m2) <br> |
|  PetscScalar | [**relative\_velocity**](#function-relative_velocity) (int cell\_i, int sp1, int sp2, int p1, int p2) <br> |
|  void | [**set\_maximal\_sigma\_v**](#function-set_maximal_sigma_v) (int coll\_i) <br> |
|  int | [**specie\_index**](#function-specie_index) (std::string specie\_name) <br> |
|  void | [**test\_sp\_reaction\_compatibility**](#function-test_sp_reaction_compatibility) (int sp1, int sp2, int reaction\_index) <br> |




























## Detailed Description




**Copyright:**

Copyright (c) 2025 




    
## Public Functions Documentation




### function COM\_velocity 

```C++
void COM_velocity (
    PetscScalar v_com,
    int cell_i,
    int sp_1,
    int sp_2,
    int p_1,
    int p_2
) 
```




<hr>



### function MCC 

```C++
void MCC () 
```




<hr>



### function add\_new\_empty\_collision 

```C++
void add_new_empty_collision (
    int sp1,
    int sp2
) 
```




<hr>



### function add\_sp\_out 

```C++
void add_sp_out (
    int coll_i,
    int reaction_index
) 
```




<hr>



### function charge\_exchange\_collision 

```C++
void charge_exchange_collision (
    int coll_i,
    int react_i,
    int cell_i,
    int sp_1,
    int sp_2,
    int p_1,
    int p_2,
    PetscScalar E
) 
```




<hr>



### function collision\_in\_cell 

```C++
void collision_in_cell (
    int coll_i,
    int cell_i
) 
```




<hr>



### function collisions\_init 

```C++
void collisions_init () 
```




<hr>



### function density\_computed\_particles 

```C++
PetscScalar density_computed_particles (
    int cell_i,
    int sp
) 
```




<hr>



### function elastic\_scattering\_collision 

```C++
void elastic_scattering_collision (
    int coll_i,
    int react_i,
    int cell_i,
    int sp_1,
    int sp_2,
    int p_1,
    int p_2,
    PetscScalar E
) 
```




<hr>



### function excitation\_collision 

```C++
void excitation_collision (
    int coll_i,
    int react_i,
    int cell_i,
    int sp_1,
    int sp_2,
    int p_1,
    int p_2,
    PetscScalar E
) 
```




<hr>



### function get\_reaction 

```C++
reactionFunction get_reaction (
    int reaction_index
) 
```




<hr>



### function ionisation\_collision 

```C++
void ionisation_collision (
    int coll_i,
    int react_i,
    int cell_i,
    int sp_1,
    int sp_2,
    int p_1,
    int p_2,
    PetscScalar E
) 
```




<hr>



### function isotropic\_scatter 

```C++
void isotropic_scatter (
    PetscScalar v_scatt
) 
```




<hr>



### function isotropic\_scatter 

```C++
void isotropic_scatter (
    PetscScalar v_scatt,
    PetscScalar v
) 
```




<hr>



### function null\_collision 

```C++
void null_collision (
    int coll_i,
    int react_i,
    int cell_i,
    int sp_1,
    int sp_2,
    int p_1,
    int p_2,
    PetscScalar E
) 
```




<hr>



### function read\_add\_collision 

```C++
void read_add_collision (
    int i
) 
```




<hr>



### function read\_add\_cross\_sections 

```C++
void read_add_cross_sections (
    int coll_i,
    std::string filename
) 
```




<hr>



### function read\_cross\_section 

```C++
PetscScalar read_cross_section (
    int coll_i,
    int j,
    PetscScalar E
) 
```




<hr>



### function recombination\_collision 

```C++
void recombination_collision (
    int coll_i,
    int react_i,
    int cell_i,
    int sp_1,
    int sp_2,
    int p_1,
    int p_2,
    PetscScalar E
) 
```




<hr>



### function reduced\_mass 

```C++
PetscScalar reduced_mass (
    PetscScalar m1,
    PetscScalar m2
) 
```




<hr>



### function relative\_velocity 

```C++
PetscScalar relative_velocity (
    int cell_i,
    int sp1,
    int sp2,
    int p1,
    int p2
) 
```




<hr>



### function set\_maximal\_sigma\_v 

```C++
void set_maximal_sigma_v (
    int coll_i
) 
```




<hr>



### function specie\_index 

```C++
int specie_index (
    std::string specie_name
) 
```




<hr>



### function test\_sp\_reaction\_compatibility 

```C++
void test_sp_reaction_compatibility (
    int sp1,
    int sp2,
    int reaction_index
) 
```




<hr>

------------------------------


