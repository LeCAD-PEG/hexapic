

# Struct Species\_load



[**ClassList**](annotated.md) **>** [**Species\_load**](structSpecies__load.md)


























## Public Attributes

| Type | Name |
| ---: | :--- |
|  PetscScalar | [**density**](#variable-density)  <br>_initial load density, in particles/m^-3_  |
|  int | [**method**](#variable-method)  <br>_load method, 0=uniform_  |
|  PetscInt | [**np2c**](#variable-np2c)  <br>_ratio of real to computed particles_  |
|  PetscScalar | [**v1drift**](#variable-v1drift)  <br> |
|  PetscScalar | [**v1thermal**](#variable-v1thermal)  <br> |
|  PetscScalar | [**v2drift**](#variable-v2drift)  <br> |
|  PetscScalar | [**v2thermal**](#variable-v2thermal)  <br> |
|  PetscScalar | [**v3drift**](#variable-v3drift)  <br>_components of drift velocity_  |
|  PetscScalar | [**v3thermal**](#variable-v3thermal)  <br>_components of thermal velocity_  |
|  PetscScalar | [**x1f**](#variable-x1f)  <br>_X-coordinate interval to load in, in m._  |
|  PetscScalar | [**x1s**](#variable-x1s)  <br> |
|  PetscScalar | [**x2f**](#variable-x2f)  <br>_Y-coordinate interval to load in, in m._  |
|  PetscScalar | [**x2s**](#variable-x2s)  <br> |












































## Public Attributes Documentation




### variable density 

_initial load density, in particles/m^-3_ 
```C++
PetscScalar Species_load::density;
```




<hr>



### variable method 

_load method, 0=uniform_ 
```C++
int Species_load::method;
```




<hr>



### variable np2c 

_ratio of real to computed particles_ 
```C++
PetscInt Species_load::np2c;
```




<hr>



### variable v1drift 

```C++
PetscScalar Species_load::v1drift;
```




<hr>



### variable v1thermal 

```C++
PetscScalar Species_load::v1thermal;
```




<hr>



### variable v2drift 

```C++
PetscScalar Species_load::v2drift;
```




<hr>



### variable v2thermal 

```C++
PetscScalar Species_load::v2thermal;
```




<hr>



### variable v3drift 

_components of drift velocity_ 
```C++
PetscScalar Species_load::v3drift;
```




<hr>



### variable v3thermal 

_components of thermal velocity_ 
```C++
PetscScalar Species_load::v3thermal;
```




<hr>



### variable x1f 

_X-coordinate interval to load in, in m._ 
```C++
PetscScalar Species_load::x1f;
```




<hr>



### variable x1s 

```C++
PetscScalar Species_load::x1s;
```




<hr>



### variable x2f 

_Y-coordinate interval to load in, in m._ 
```C++
PetscScalar Species_load::x2f;
```




<hr>



### variable x2s 

```C++
PetscScalar Species_load::x2s;
```




<hr>

------------------------------


