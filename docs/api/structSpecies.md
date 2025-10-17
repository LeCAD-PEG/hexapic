

# Struct Species



[**ClassList**](annotated.md) **>** [**Species**](structSpecies.md)


























## Public Attributes

| Type | Name |
| ---: | :--- |
|  std::vector&lt; PetscScalar &gt; | [**T\_store**](#variable-t_store)  <br>_temperature diagnostic storage_  |
|  int | [**charge**](#variable-charge)  <br>_specie charge/elementary charge_  |
|  PetscScalar | [**chom**](#variable-chom)  <br>_specie mass and charge/mass ratio_  |
|  std::vector&lt; [**Species\_inject**](structSpecies__inject.md) &gt; | [**inject**](#variable-inject)  <br>_particle injection at boundary_  |
|  std::vector&lt; [**Species\_load**](structSpecies__load.md) &gt; | [**load**](#variable-load)  <br>_allows multiple loads with same species_  |
|  PetscScalar | [**mass**](#variable-mass)  <br> |
|  unsigned int | [**n**](#variable-n)  <br>_number of particles of a specie_  |
|  std::vector&lt; PetscScalar &gt; | [**n\_store**](#variable-n_store)  <br>_particle density diagnostic storage_  |
|  std::string | [**name**](#variable-name)  <br>_specie name, e.g. electron_  |
|  PetscScalar | [**sb**](#variable-sb)  <br>_specie coefficients for Boris mover_  |
|  std::vector&lt; [**SourceSpecie**](structSourceSpecie.md) &gt; | [**source**](#variable-source)  <br>_volumetric specie sources_  |
|  unsigned int | [**subcycle**](#variable-subcycle)  <br> |
|  PetscScalar | [**tb**](#variable-tb)  <br> |
|  std::vector&lt; PetscScalar &gt; | [**v\_store**](#variable-v_store)  <br>_fluid velocity diagnostic storage_  |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**correct\_inject**](#function-correct_inject) (int Nx, int Ny, int Nx\_old, int Ny\_old) <br> |




























## Public Attributes Documentation




### variable T\_store 

_temperature diagnostic storage_ 
```C++
std::vector< PetscScalar > Species::T_store;
```




<hr>



### variable charge 

_specie charge/elementary charge_ 
```C++
int Species::charge;
```




<hr>



### variable chom 

_specie mass and charge/mass ratio_ 
```C++
PetscScalar Species::chom;
```




<hr>



### variable inject 

_particle injection at boundary_ 
```C++
std::vector<Species_inject> Species::inject;
```




<hr>



### variable load 

_allows multiple loads with same species_ 
```C++
std::vector<Species_load> Species::load;
```




<hr>



### variable mass 

```C++
PetscScalar Species::mass;
```




<hr>



### variable n 

_number of particles of a specie_ 
```C++
unsigned int Species::n;
```




<hr>



### variable n\_store 

_particle density diagnostic storage_ 
```C++
std::vector< PetscScalar > Species::n_store;
```




<hr>



### variable name 

_specie name, e.g. electron_ 
```C++
std::string Species::name;
```




<hr>



### variable sb 

_specie coefficients for Boris mover_ 
```C++
PetscScalar Species::sb[3];
```




<hr>



### variable source 

_volumetric specie sources_ 
```C++
std::vector<SourceSpecie> Species::source;
```




<hr>



### variable subcycle 

```C++
unsigned int Species::subcycle;
```




<hr>



### variable tb 

```C++
PetscScalar Species::tb[3];
```




<hr>



### variable v\_store 

_fluid velocity diagnostic storage_ 
```C++
std::vector< PetscScalar > Species::v_store;
```




<hr>
## Public Functions Documentation




### function correct\_inject 

```C++
inline void Species::correct_inject (
    int Nx,
    int Ny,
    int Nx_old,
    int Ny_old
) 
```




<hr>

------------------------------


