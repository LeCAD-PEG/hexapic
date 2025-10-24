

# Struct Collision



[**ClassList**](annotated.md) **>** [**Collision**](structCollision.md)


























## Public Attributes

| Type | Name |
| ---: | :--- |
|  std::vector&lt; PetscScalar &gt; | [**E\_min**](#variable-e_min)  <br>_minimal energy required for reaction_  |
|  int | [**n\_reactions**](#variable-n_reactions)  <br>_number of reactions between this two species_  |
|  std::vector&lt; reactionFunction &gt; | [**processCollision**](#variable-processcollision)  <br> |
|  std::vector&lt; std::string &gt; | [**reaction\_names**](#variable-reaction_names)  <br>_collision reaction names_  |
|  std::vector&lt; std::set&lt; std::pair&lt; PetscScalar, PetscScalar &gt; &gt; &gt; | [**sigma**](#variable-sigma)  <br> |
|  PetscScalar | [**sigma\_v\_max**](#variable-sigma_v_max)  <br>_maximal over all energies (sigma\_total\*v)_  |
|  int | [**sp**](#variable-sp)  <br>_colliding species indexes_  |
|  std::vector&lt; std::vector&lt; int &gt; &gt; | [**sp\_out**](#variable-sp_out)  <br>_new species created in collision_  |












































## Public Attributes Documentation




### variable E\_min 

_minimal energy required for reaction_ 
```C++
std::vector< PetscScalar > Collision::E_min;
```




<hr>



### variable n\_reactions 

_number of reactions between this two species_ 
```C++
int Collision::n_reactions;
```




<hr>



### variable processCollision 

```C++
std::vector<reactionFunction> Collision::processCollision;
```




<hr>



### variable reaction\_names 

_collision reaction names_ 
```C++
std::vector<std::string> Collision::reaction_names;
```




<hr>



### variable sigma 

```C++
std::vector<std::set<std::pair< PetscScalar , PetscScalar > > > Collision::sigma;
```




<hr>



### variable sigma\_v\_max 

_maximal over all energies (sigma\_total\*v)_ 
```C++
PetscScalar Collision::sigma_v_max;
```




<hr>



### variable sp 

_colliding species indexes_ 
```C++
int Collision::sp[2];
```




<hr>



### variable sp\_out 

_new species created in collision_ 
```C++
std::vector<std::vector<int> > Collision::sp_out;
```




<hr>

------------------------------


