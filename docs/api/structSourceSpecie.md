

# Struct SourceSpecie



[**ClassList**](annotated.md) **>** [**SourceSpecie**](structSourceSpecie.md)


























## Public Attributes

| Type | Name |
| ---: | :--- |
|  PetscScalar | [**heat\_frac**](#variable-heat_frac)  <br>_freq\*heat\_frac to get specie specific heating_  |
|  int | [**index**](#variable-index)  <br>_source index (ordering done after reading input file)_  |
|  PetscScalar | [**source\_frac**](#variable-source_frac)  <br>_n\_inj \* source\_frac to get specie specific injecting_  |
|  PetscScalar | [**vtpar**](#variable-vtpar)  <br> |
|  PetscScalar | [**vtper**](#variable-vtper)  <br>_parallel and perpendicular thermal velocity_  |












































## Public Attributes Documentation




### variable heat\_frac 

_freq\*heat\_frac to get specie specific heating_ 
```C++
PetscScalar SourceSpecie::heat_frac;
```




<hr>



### variable index 

_source index (ordering done after reading input file)_ 
```C++
int SourceSpecie::index;
```




<hr>



### variable source\_frac 

_n\_inj \* source\_frac to get specie specific injecting_ 
```C++
PetscScalar SourceSpecie::source_frac;
```




<hr>



### variable vtpar 

```C++
PetscScalar SourceSpecie::vtpar;
```




<hr>



### variable vtper 

_parallel and perpendicular thermal velocity_ 
```C++
PetscScalar SourceSpecie::vtper;
```




<hr>

------------------------------


