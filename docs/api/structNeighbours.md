

# Struct Neighbours



[**ClassList**](annotated.md) **>** [**Neighbours**](structNeighbours.md)


























## Public Attributes

| Type | Name |
| ---: | :--- |
|  std::set&lt; int &gt; | [**all**](#variable-all)  <br>_set of all neighbours (to expect data)_  |
|  std::vector&lt; int &gt; | [**corners**](#variable-corners)  <br>_ranks of corner neighbours_  |
|  int | [**max\_recv\_size**](#variable-max_recv_size)  <br>_maximal size of received data_  |
|  std::vector&lt; std::vector&lt; PetscScalar &gt; &gt; | [**recv\_data**](#variable-recv_data)  <br>_receive buffer_  |
|  std::vector&lt; int &gt; | [**rxmax**](#variable-rxmax)  <br> |
|  std::vector&lt; int &gt; | [**rxmin**](#variable-rxmin)  <br> |
|  std::vector&lt; int &gt; | [**rymax**](#variable-rymax)  <br>_ranks of neighbours_  |
|  std::vector&lt; int &gt; | [**rymin**](#variable-rymin)  <br> |
|  std::vector&lt; std::vector&lt; PetscScalar &gt; &gt; | [**send\_data**](#variable-send_data)  <br>_send buffer_  |
|  std::vector&lt; std::vector&lt; std::vector&lt; PetscScalar &gt; &gt; &gt; | [**send\_data\_3d**](#variable-send_data_3d)  <br> |
|  std::vector&lt; std::vector&lt; PetscScalar &gt; &gt; | [**shift\_xy**](#variable-shift_xy)  <br>_shift between local grids_  |












































## Public Attributes Documentation




### variable all 

_set of all neighbours (to expect data)_ 
```C++
std::set<int> Neighbours::all;
```




<hr>



### variable corners 

_ranks of corner neighbours_ 
```C++
std::vector<int> Neighbours::corners;
```




<hr>



### variable max\_recv\_size 

_maximal size of received data_ 
```C++
int Neighbours::max_recv_size;
```




<hr>



### variable recv\_data 

_receive buffer_ 
```C++
std::vector<std::vector< PetscScalar > > Neighbours::recv_data;
```




<hr>



### variable rxmax 

```C++
std::vector<int> Neighbours::rxmax;
```




<hr>



### variable rxmin 

```C++
std::vector<int> Neighbours::rxmin;
```




<hr>



### variable rymax 

_ranks of neighbours_ 
```C++
std::vector<int> Neighbours::rymax;
```




<hr>



### variable rymin 

```C++
std::vector<int> Neighbours::rymin;
```




<hr>



### variable send\_data 

_send buffer_ 
```C++
std::vector<std::vector< PetscScalar > > Neighbours::send_data;
```




<hr>



### variable send\_data\_3d 

```C++
std::vector<std::vector<std::vector< PetscScalar > > > Neighbours::send_data_3d;
```




<hr>



### variable shift\_xy 

_shift between local grids_ 
```C++
std::vector<std::vector< PetscScalar > > Neighbours::shift_xy;
```




<hr>

------------------------------


