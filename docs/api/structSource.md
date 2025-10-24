

# Struct Source



[**ClassList**](annotated.md) **>** [**Source**](structSource.md)


























## Public Attributes

| Type | Name |
| ---: | :--- |
|  PetscScalar | [**Fx0**](#variable-fx0)  <br> |
|  PetscScalar | [**Fx1**](#variable-fx1)  <br>_from cumulative distribution part on this rank_  |
|  std::set&lt; std::pair&lt; PetscScalar, int &gt; &gt; | [**Fxcell**](#variable-fxcell)  <br> |
|  PetscScalar | [**Fy0**](#variable-fy0)  <br> |
|  PetscScalar | [**Fy1**](#variable-fy1)  <br>_same for y_  |
|  std::set&lt; std::pair&lt; PetscScalar, int &gt; &gt; | [**Fycell**](#variable-fycell)  <br>_for normal distribution_  |
|  PetscInt | [**Lx**](#variable-lx)  <br> |
|  PetscInt | [**Ly**](#variable-ly)  <br>_x, y full length (global)_  |
|  PetscScalar | [**Tpar**](#variable-tpar)  <br> |
|  PetscScalar | [**Tper**](#variable-tper)  <br>_temperature parallel and perpendicular to B [eV]_  |
|  int | [**active**](#variable-active)  <br>_if 0 not active_  |
|  PetscScalar | [**cx**](#variable-cx)  <br> |
|  PetscScalar | [**cy**](#variable-cy)  <br>_center position of global source (in local coordinates)_  |
|  PetscScalar | [**heat\_freq**](#variable-heat_freq)  <br>_frequency of collisions of a particle inside region_  |
|  int | [**index**](#variable-index)  <br>_source index_  |
|  PetscInt | [**j1**](#variable-j1)  <br> |
|  PetscInt | [**j2**](#variable-j2)  <br> |
|  PetscInt | [**k1**](#variable-k1)  <br>_bottom left corner_  |
|  PetscInt | [**k2**](#variable-k2)  <br>_top right corner_  |
|  PetscInt | [**lx**](#variable-lx)  <br> |
|  PetscInt | [**ly**](#variable-ly)  <br>_x, y local length_  |
|  PetscScalar | [**n\_heat**](#variable-n_heat)  <br>_number of heat collisions of a particle in timestep_  |
|  PetscScalar | [**n\_inj**](#variable-n_inj)  <br>_[s^-1 m^-3] /nc2p when injecting_  |
|  PetscScalar | [**n\_inj\_step**](#variable-n_inj_step)  <br>_injected computed particles per timestep_  |
|  int | [**shape\_x**](#variable-shape_x)  <br>_0-uniform, 1-cos, 2-normal, details in implementation_  |
|  int | [**shape\_y**](#variable-shape_y)  <br>_0-uniform, 1-cos, 2-normal_  |
|  PetscScalar | [**t0**](#variable-t0)  <br> |
|  PetscScalar | [**t1**](#variable-t1)  <br> |
|  PetscScalar | [**t2**](#variable-t2)  <br> |
|  PetscScalar | [**t3**](#variable-t3)  <br>_[t0, t1] source 0-&gt;1, [t2, t3] source 1-&gt;0 (linear)_  |
|  int | [**t\_dep**](#variable-t_dep)  <br>_time-dependend (source and heating)_  |












































## Public Attributes Documentation




### variable Fx0 

```C++
PetscScalar Source::Fx0;
```




<hr>



### variable Fx1 

_from cumulative distribution part on this rank_ 
```C++
PetscScalar Source::Fx1;
```




<hr>



### variable Fxcell 

```C++
std::set<std::pair< PetscScalar , int> > Source::Fxcell;
```




<hr>



### variable Fy0 

```C++
PetscScalar Source::Fy0;
```




<hr>



### variable Fy1 

_same for y_ 
```C++
PetscScalar Source::Fy1;
```




<hr>



### variable Fycell 

_for normal distribution_ 
```C++
std::set<std::pair< PetscScalar , int> > Source::Fycell;
```




<hr>



### variable Lx 

```C++
PetscInt Source::Lx;
```




<hr>



### variable Ly 

_x, y full length (global)_ 
```C++
PetscInt Source::Ly;
```




<hr>



### variable Tpar 

```C++
PetscScalar Source::Tpar;
```




<hr>



### variable Tper 

_temperature parallel and perpendicular to B [eV]_ 
```C++
PetscScalar Source::Tper;
```




<hr>



### variable active 

_if 0 not active_ 
```C++
int Source::active;
```




<hr>



### variable cx 

```C++
PetscScalar Source::cx;
```




<hr>



### variable cy 

_center position of global source (in local coordinates)_ 
```C++
PetscScalar Source::cy;
```




<hr>



### variable heat\_freq 

_frequency of collisions of a particle inside region_ 
```C++
PetscScalar Source::heat_freq;
```




<hr>



### variable index 

_source index_ 
```C++
int Source::index;
```




<hr>



### variable j1 

```C++
PetscInt Source::j1;
```




<hr>



### variable j2 

```C++
PetscInt Source::j2;
```




<hr>



### variable k1 

_bottom left corner_ 
```C++
PetscInt Source::k1;
```




<hr>



### variable k2 

_top right corner_ 
```C++
PetscInt Source::k2;
```




<hr>



### variable lx 

```C++
PetscInt Source::lx;
```




<hr>



### variable ly 

_x, y local length_ 
```C++
PetscInt Source::ly;
```




<hr>



### variable n\_heat 

_number of heat collisions of a particle in timestep_ 
```C++
PetscScalar Source::n_heat;
```




<hr>



### variable n\_inj 

_[s^-1 m^-3] /nc2p when injecting_ 
```C++
PetscScalar Source::n_inj;
```




<hr>



### variable n\_inj\_step 

_injected computed particles per timestep_ 
```C++
PetscScalar Source::n_inj_step;
```




<hr>



### variable shape\_x 

_0-uniform, 1-cos, 2-normal, details in implementation_ 
```C++
int Source::shape_x;
```




<hr>



### variable shape\_y 

_0-uniform, 1-cos, 2-normal_ 
```C++
int Source::shape_y;
```




<hr>



### variable t0 

```C++
PetscScalar Source::t0;
```




<hr>



### variable t1 

```C++
PetscScalar Source::t1;
```




<hr>



### variable t2 

```C++
PetscScalar Source::t2;
```




<hr>



### variable t3 

_[t0, t1] source 0-&gt;1, [t2, t3] source 1-&gt;0 (linear)_ 
```C++
PetscScalar Source::t3;
```




<hr>



### variable t\_dep 

_time-dependend (source and heating)_ 
```C++
int Source::t_dep;
```




<hr>

------------------------------


