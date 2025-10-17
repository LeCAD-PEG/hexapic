

# Struct Boundary



[**ClassList**](annotated.md) **>** [**Boundary**](structBoundary.md)


























## Public Attributes

| Type | Name |
| ---: | :--- |
|  PetscScalar | [**C**](#variable-c)  <br>_fixed potential on boundary, in V_  |
|  PetscScalar | [**QuseFlag**](#variable-quseflag)  <br>_charge accumulation (for neumann=0)_  |
|  PetscInt | [**j1**](#variable-j1)  <br> |
|  PetscInt | [**j2**](#variable-j2)  <br>_start-end indexes of the boundary along X_  |
|  PetscInt | [**k1**](#variable-k1)  <br> |
|  PetscInt | [**k2**](#variable-k2)  <br>_start-end indexes of the boundary along Y_  |
|  std::string | [**name**](#variable-name)  <br>_boundary name, e.g. Equipotential_  |
|  PetscScalar | [**reflection**](#variable-reflection)  <br>_part of reflected particles (for neumann=1)_  |
|  std::vector&lt; std::vector&lt; [**SecondaryEmission**](structSecondaryEmission.md) &gt; &gt; | [**sp\_sec**](#variable-sp_sec)  <br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**add\_parameters**](#function-add_parameters) (std::string NAME, PetscInt J1, PetscInt J2, PetscInt K1, PetscInt K2) <br> |




























## Public Attributes Documentation




### variable C 

_fixed potential on boundary, in V_ 
```C++
PetscScalar Boundary::C;
```




<hr>



### variable QuseFlag 

_charge accumulation (for neumann=0)_ 
```C++
PetscScalar Boundary::QuseFlag;
```




<hr>



### variable j1 

```C++
PetscInt Boundary::j1;
```




<hr>



### variable j2 

_start-end indexes of the boundary along X_ 
```C++
PetscInt Boundary::j2;
```




<hr>



### variable k1 

```C++
PetscInt Boundary::k1;
```




<hr>



### variable k2 

_start-end indexes of the boundary along Y_ 
```C++
PetscInt Boundary::k2;
```




<hr>



### variable name 

_boundary name, e.g. Equipotential_ 
```C++
std::string Boundary::name;
```




<hr>



### variable reflection 

_part of reflected particles (for neumann=1)_ 
```C++
PetscScalar Boundary::reflection;
```




<hr>



### variable sp\_sec 

```C++
std::vector<std::vector<SecondaryEmission> > Boundary::sp_sec;
```




<hr>
## Public Functions Documentation




### function add\_parameters 

```C++
inline void Boundary::add_parameters (
    std::string NAME,
    PetscInt J1,
    PetscInt J2,
    PetscInt K1,
    PetscInt K2
) 
```




<hr>

------------------------------


