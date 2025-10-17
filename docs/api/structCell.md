

# Struct Cell



[**ClassList**](annotated.md) **>** [**Cell**](structCell.md)


























## Public Attributes

| Type | Name |
| ---: | :--- |
|  PetscInt | [**n0x**](#variable-n0x)  <br> |
|  PetscInt | [**n0y**](#variable-n0y)  <br>_position of bottom left corrner in local grid_  |
|  std::vector&lt; int &gt; | [**neighbours**](#variable-neighbours)  <br>_all neighbouring cells_  |
|  std::vector&lt; std::vector&lt; [**Particle**](structParticle.md) &gt; &gt; | [**part**](#variable-part)  <br>_particle of each specie in cell_  |












































## Public Attributes Documentation




### variable n0x 

```C++
PetscInt Cell::n0x;
```




<hr>



### variable n0y 

_position of bottom left corrner in local grid_ 
```C++
PetscInt Cell::n0y;
```




<hr>



### variable neighbours 

_all neighbouring cells_ 
```C++
std::vector<int> Cell::neighbours;
```




<hr>



### variable part 

_particle of each specie in cell_ 
```C++
std::vector<std::vector<Particle> > Cell::part;
```




<hr>

------------------------------


