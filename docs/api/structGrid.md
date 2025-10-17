

# Struct Grid



[**ClassList**](annotated.md) **>** [**Grid**](structGrid.md)


























## Public Attributes

| Type | Name |
| ---: | :--- |
|  PetscInt | [**N0x**](#variable-n0x)  <br> |
|  PetscInt | [**N0y**](#variable-n0y)  <br>_shift of grid cells of local grid in global grid_  |
|  PetscInt | [**Ncx**](#variable-ncx)  <br> |
|  PetscInt | [**Ncy**](#variable-ncy)  <br>_number of local grid cells in X and Y directions_  |
|  PetscInt | [**Nx**](#variable-nx)  <br> |
|  PetscInt | [**Nxg**](#variable-nxg)  <br> |
|  PetscInt | [**Ny**](#variable-ny)  <br>_number of local grid edges in X and Y directions_  |
|  PetscInt | [**Nyg**](#variable-nyg)  <br>_number of global grid edges in X and Y directions_  |
|  std::vector&lt; PetscScalar &gt; | [**V\_store**](#variable-v_store)  <br>_electric potential diagnostic storage_  |
|  std::vector&lt; [**Boundary**](structBoundary.md) &gt; | [**boundaries**](#variable-boundaries)  <br>_also includes periodic boundaries_  |
|  PetscScalar | [**dx**](#variable-dx)  <br> |
|  PetscScalar | [**dy**](#variable-dy)  <br>_cell dimensions along X and Y coordinates, in m_  |
|  [**Neighbours**](structNeighbours.md) | [**neighbours**](#variable-neighbours)  <br> |
|  PetscScalar | [**p2n**](#variable-p2n)  <br>_particle-to-density constant_  |
|  PetscScalar | [**x1f**](#variable-x1f)  <br>_start and end extension along X, in m_  |
|  PetscScalar | [**x1s**](#variable-x1s)  <br> |
|  PetscScalar | [**x2f**](#variable-x2f)  <br>_start and end extension along Y, in m_  |
|  PetscScalar | [**x2s**](#variable-x2s)  <br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**correct\_boundaries**](#function-correct_boundaries) () <br> |
|  void | [**sort\_boundaries**](#function-sort_boundaries) () <br> |




























## Public Attributes Documentation




### variable N0x 

```C++
PetscInt Grid::N0x;
```




<hr>



### variable N0y 

_shift of grid cells of local grid in global grid_ 
```C++
PetscInt Grid::N0y;
```




<hr>



### variable Ncx 

```C++
PetscInt Grid::Ncx;
```




<hr>



### variable Ncy 

_number of local grid cells in X and Y directions_ 
```C++
PetscInt Grid::Ncy;
```




<hr>



### variable Nx 

```C++
PetscInt Grid::Nx;
```




<hr>



### variable Nxg 

```C++
PetscInt Grid::Nxg;
```




<hr>



### variable Ny 

_number of local grid edges in X and Y directions_ 
```C++
PetscInt Grid::Ny;
```




<hr>



### variable Nyg 

_number of global grid edges in X and Y directions_ 
```C++
PetscInt Grid::Nyg;
```




<hr>



### variable V\_store 

_electric potential diagnostic storage_ 
```C++
std::vector< PetscScalar > Grid::V_store;
```




<hr>



### variable boundaries 

_also includes periodic boundaries_ 
```C++
std::vector<Boundary> Grid::boundaries;
```




<hr>



### variable dx 

```C++
PetscScalar Grid::dx;
```




<hr>



### variable dy 

_cell dimensions along X and Y coordinates, in m_ 
```C++
PetscScalar Grid::dy;
```




<hr>



### variable neighbours 

```C++
Neighbours Grid::neighbours;
```




<hr>



### variable p2n 

_particle-to-density constant_ 
```C++
PetscScalar Grid::p2n;
```




<hr>



### variable x1f 

_start and end extension along X, in m_ 
```C++
PetscScalar Grid::x1f;
```




<hr>



### variable x1s 

```C++
PetscScalar Grid::x1s;
```




<hr>



### variable x2f 

_start and end extension along Y, in m_ 
```C++
PetscScalar Grid::x2f;
```




<hr>



### variable x2s 

```C++
PetscScalar Grid::x2s;
```




<hr>
## Public Functions Documentation




### function correct\_boundaries 

```C++
inline void Grid::correct_boundaries () 
```




<hr>



### function sort\_boundaries 

```C++
inline void Grid::sort_boundaries () 
```




<hr>

------------------------------


