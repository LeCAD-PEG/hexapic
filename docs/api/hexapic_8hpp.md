

# File hexapic.hpp



[**FileList**](files.md) **>** [**hexapic**](dir_b8175b39202167732592fd9c59b14660.md) **>** [**source**](dir_4d88c249dac475c11ca7037572002e86.md) **>** [**hexapic.hpp**](hexapic_8hpp.md)



_Core data structures, types, and global parameters for the HEXAPIC Particle-in-Cell (PIC) code._ [More...](#detailed-description)

* `#include <mpi.h>`
* `#include <math.h>`
* `#include <stdio.h>`
* `#include <unistd.h>`
* `#include <stdlib.h>`
* `#include <string.h>`
* `#include <limits.h>`
* `#include <time.h>`
* `#include <iostream>`
* `#include <vector>`
* `#include <set>`
* `#include <fstream>`
* `#include <sstream>`
* `#include <petsc.h>`
* `#include <iomanip>`
* `#include <openPMD/openPMD.hpp>`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**openPMD**](namespaceopenPMD.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| struct | [**Boundary**](structBoundary.md) <br> |
| struct | [**Cell**](structCell.md) <br> |
| struct | [**Collision**](structCollision.md) <br> |
| struct | [**Grid**](structGrid.md) <br> |
| struct | [**Neighbours**](structNeighbours.md) <br> |
| struct | [**Particle**](structParticle.md) <br> |
| struct | [**SecondaryEmission**](structSecondaryEmission.md) <br> |
| struct | [**Source**](structSource.md) <br> |
| struct | [**SourceSpecie**](structSourceSpecie.md) <br> |
| struct | [**Species**](structSpecies.md) <br> |
| struct | [**Species\_inject**](structSpecies__inject.md) <br> |
| struct | [**Species\_load**](structSpecies__load.md) <br> |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef void(\* | [**reactionFunction**](#typedef-reactionfunction)  <br> |




## Public Attributes

| Type | Name |
| ---: | :--- |
|  Mat | [**A**](#variable-a)  <br>_Matrix._  |
|  PetscScalar | [**Bb**](#variable-bb)  <br>_normalised basis for B_  |
|  PetscScalar | [**Bf**](#variable-bf)  <br>_external magnetic field_  |
|  std::vector&lt; std::vector&lt; std::string &gt; &gt; | [**MCC\_init\_data**](#variable-mcc_init_data)  <br> |
|  std::vector&lt; std::vector&lt; std::string &gt; &gt; | [**PSI\_init\_data**](#variable-psi_init_data)  <br> |
|  PetscScalar | [**Sconst**](#variable-sconst)  <br> |
|  PetscScalar \* | [**V\_global**](#variable-v_global)  <br> |
|  PetscScalar \* | [**V\_local**](#variable-v_local)  <br> |
|  Vec | [**b**](#variable-b)  <br>_Solution and source term._  |
|  PetscInt | [**boris**](#variable-boris)  <br> |
|  std::vector&lt; [**Cell**](structCell.md) &gt; | [**cells**](#variable-cells)  <br> |
|  std::vector&lt; std::vector&lt; PetscScalar &gt; &gt; | [**chunks\_float**](#variable-chunks_float)  <br> |
|  std::vector&lt; [**Collision**](structCollision.md) &gt; | [**collisions**](#variable-collisions)  <br> |
|  VecScatter | [**ctx**](#variable-ctx)  <br> |
|  PetscScalar | [**dV**](#variable-dv)  <br> |
|  PetscInt | [**davg**](#variable-davg)  <br> |
|  bool | [**digital\_smoothing**](#variable-digital_smoothing)  <br> |
|  PetscInt | [**dsteps**](#variable-dsteps)  <br> |
|  PetscScalar | [**dt**](#variable-dt)  <br> |
|  PetscScalar | [**dx**](#variable-dx)  <br> |
|  PetscScalar | [**dx2\_inv**](#variable-dx2_inv)  <br> |
|  PetscScalar | [**dy**](#variable-dy)  <br> |
|  PetscScalar | [**dy2\_inv**](#variable-dy2_inv)  <br> |
|  double | [**et**](#variable-et)  <br>_elapsed time for profiling_  |
|  FILE \* | [**f**](#variable-f)  <br> |
|  [**Grid**](structGrid.md) | [**grid**](#variable-grid)  <br> |
|  PetscScalar | [**hdt**](#variable-hdt)  <br> |
|  PetscErrorCode | [**ierr**](#variable-ierr)  <br>_Error catcher._  |
|  KSP | [**ksp**](#variable-ksp)  <br>_Krylov subspace solver context._  |
|  std::vector&lt; MeshRecordComponent &gt; | [**meshes**](#variable-meshes)  <br> |
|  int | [**nmpi**](#variable-nmpi)  <br> |
|  PC | [**pc**](#variable-pc)  <br>_Preconditioner context._  |
|  std::vector&lt; PetscScalar &gt; | [**rho**](#variable-rho)  <br> |
|  std::vector&lt; PetscInt &gt; | [**rho\_index**](#variable-rho_index)  <br> |
|  int | [**rmpi**](#variable-rmpi)  <br> |
|  Series | [**series**](#variable-series)  <br> |
|  std::vector&lt; [**Source**](structSource.md) &gt; | [**sources**](#variable-sources)  <br> |
|  std::vector&lt; [**Species**](structSpecies.md) &gt; | [**species**](#variable-species)  <br> |
|  std::vector&lt; PetscScalar &gt; | [**t\_store**](#variable-t_store)  <br> |
|  std::vector&lt; PetscScalar &gt; | [**trajx**](#variable-trajx)  <br> |
|  std::vector&lt; PetscScalar &gt; | [**trajy**](#variable-trajy)  <br> |
|  std::vector&lt; PetscScalar &gt; | [**trajz**](#variable-trajz)  <br> |
|  PetscInt | [**tstep**](#variable-tstep)  <br> |
|  Vec | [**x**](#variable-x)  <br> |
|  Vec | [**x\_array**](#variable-x_array)  <br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**MCC**](#function-mcc) () <br> |
|  void | [**PSI\_init**](#function-psi_init) () <br> |
|  void | [**add\_received\_particles**](#function-add_received_particles) () <br> |
|  void | [**cells\_init**](#function-cells_init) () <br> |
|  void | [**collisions\_init**](#function-collisions_init) () <br> |
|  void | [**create\_output**](#function-create_output) () <br> |
|  void | [**decompose\_domain**](#function-decompose_domain) () <br> |
|  void | [**field\_solver\_petsc**](#function-field_solver_petsc) () <br> |
|  void | [**field\_source\_update**](#function-field_source_update) () <br> |
|  void | [**grid2part**](#function-grid2part) (PetscScalar, PetscScalar, PetscScalar \*, PetscScalar \*) <br> |
|  void | [**grid\_init**](#function-grid_init) (int, char \*\*) <br> |
|  void | [**initial\_particle\_load**](#function-initial_particle_load) () <br> |
|  void | [**inject\_particle**](#function-inject_particle) (int, PetscScalar, PetscScalar) <br> |
|  void | [**inject\_particles**](#function-inject_particles) () <br> |
|  PetscScalar | [**normvel\_MaxBol**](#function-normvel_maxbol) () <br> |
|  PetscScalar | [**normvel\_vFv**](#function-normvel_vfv) () <br> |
|  void | [**num\_param\_init**](#function-num_param_init) () <br> |
|  void | [**part2grid**](#function-part2grid) () <br> |
|  void | [**particle\_boundaries**](#function-particle_boundaries) () <br> |
|  void | [**particle\_boundaries\_dd**](#function-particle_boundaries_dd) () <br> |
|  void | [**particle\_mover\_boris**](#function-particle_mover_boris) () <br> |
|  void | [**particle\_wall\_interaction**](#function-particle_wall_interaction) (int, int, int, int) <br> |
|  void | [**particles\_init**](#function-particles_init) () <br> |
|  void | [**petsc\_cleanup**](#function-petsc_cleanup) () <br> |
|  void | [**petsc\_init**](#function-petsc_init) (int, char \*\*) <br> |
|  void | [**petsc\_plot\_finalize**](#function-petsc_plot_finalize) () <br> |
|  void | [**petsc\_plot\_init**](#function-petsc_plot_init) () <br> |
|  void | [**petsc\_plot\_update**](#function-petsc_plot_update) (int) <br> |
|  void | [**plot\_profiling**](#function-plot_profiling) () <br> |
|  void | [**read\_input\_file**](#function-read_input_file) (char \*\*) <br> |
|  void | [**remove\_particle**](#function-remove_particle) (int, int, int) <br> |
|  void | [**save\_output**](#function-save_output) () <br> |
|  void | [**send\_receive\_particles**](#function-send_receive_particles) () <br> |
|  void | [**source**](#function-source) () <br> |
|  void | [**source\_heating**](#function-source_heating) () <br> |
|  void | [**source\_init**](#function-source_init) () <br> |
|  int | [**specie\_index**](#function-specie_index) (std::string) <br> |



























## Macros

| Type | Name |
| ---: | :--- |
| define  | [**EPS0**](hexapic_8hpp.md#define-eps0)  `8.8542e-12`<br> |
| define  | [**EXTERNAL**](hexapic_8hpp.md#define-external)  `extern`<br> |
| define  | [**FALSE**](hexapic_8hpp.md#define-false)  `0`<br> |
| define  | [**PI**](hexapic_8hpp.md#define-pi)  `3.141592653589793`<br> |
| define  | [**REAL**](hexapic_8hpp.md#define-real)  `PetscScalar`<br> |
| define  | [**TRUE**](hexapic_8hpp.md#define-true)  `1`<br> |
| define  | [**TWOPI**](hexapic_8hpp.md#define-twopi)  `6.28318530717959`<br> |
| define  | [**comm**](hexapic_8hpp.md#define-comm)  `MPI\_COMM\_WORLD`<br> |
| define  | [**frand**](hexapic_8hpp.md#define-frand) () `((REAL) rand() / (RAND\_MAX+1.0))`<br> |
| define  | [**np2c\_global**](hexapic_8hpp.md#define-np2c_global)  `2e6`<br> |
| define  | [**plot\_3D**](hexapic_8hpp.md#define-plot_3d)  `FALSE`<br> |
| define  | [**q\_e**](hexapic_8hpp.md#define-q_e)  `1.602e-19`<br> |

## Detailed Description




**Copyright:**

Copyright (c) 2025 




    
## Public Types Documentation




### typedef reactionFunction 

```C++
typedef void(* reactionFunction) (int, int, int, int, int, int, int, PetscScalar);
```




<hr>
## Public Attributes Documentation




### variable A 

_Matrix._ 
```C++
Mat A;
```




<hr>



### variable Bb 

_normalised basis for B_ 
```C++
PetscScalar Bb[3][3];
```




<hr>



### variable Bf 

_external magnetic field_ 
```C++
PetscScalar Bf[3];
```




<hr>



### variable MCC\_init\_data 

```C++
std::vector<std::vector<std::string> > MCC_init_data;
```




<hr>



### variable PSI\_init\_data 

```C++
std::vector<std::vector<std::string> > PSI_init_data;
```




<hr>



### variable Sconst 

```C++
PetscScalar Sconst;
```




<hr>



### variable V\_global 

```C++
PetscScalar * V_global;
```




<hr>



### variable V\_local 

```C++
PetscScalar* V_local;
```




<hr>



### variable b 

_Solution and source term._ 
```C++
Vec b;
```




<hr>



### variable boris 

```C++
PetscInt boris;
```




<hr>



### variable cells 

```C++
std::vector<Cell> cells;
```




<hr>



### variable chunks\_float 

```C++
std::vector<std::vector< PetscScalar > > chunks_float;
```




<hr>



### variable collisions 

```C++
std::vector<Collision> collisions;
```




<hr>



### variable ctx 

```C++
VecScatter ctx;
```




<hr>



### variable dV 

```C++
PetscScalar dV;
```




<hr>



### variable davg 

```C++
PetscInt davg;
```




<hr>



### variable digital\_smoothing 

```C++
bool digital_smoothing;
```




<hr>



### variable dsteps 

```C++
PetscInt dsteps;
```




<hr>



### variable dt 

```C++
PetscScalar dt;
```




<hr>



### variable dx 

```C++
PetscScalar dx;
```




<hr>



### variable dx2\_inv 

```C++
PetscScalar dx2_inv;
```




<hr>



### variable dy 

```C++
PetscScalar dy;
```




<hr>



### variable dy2\_inv 

```C++
PetscScalar dy2_inv;
```




<hr>



### variable et 

_elapsed time for profiling_ 
```C++
double et[10];
```




<hr>



### variable f 

```C++
FILE* f;
```




<hr>



### variable grid 

```C++
Grid grid;
```




<hr>



### variable hdt 

```C++
PetscScalar hdt;
```




<hr>



### variable ierr 

_Error catcher._ 
```C++
PetscErrorCode ierr;
```




<hr>



### variable ksp 

_Krylov subspace solver context._ 
```C++
KSP ksp;
```




<hr>



### variable meshes 

```C++
std::vector<MeshRecordComponent> meshes;
```




<hr>



### variable nmpi 

```C++
int nmpi;
```




<hr>



### variable pc 

_Preconditioner context._ 
```C++
PC pc;
```




<hr>



### variable rho 

```C++
std::vector<PetscScalar> rho;
```




<hr>



### variable rho\_index 

```C++
std::vector<PetscInt> rho_index;
```




<hr>



### variable rmpi 

```C++
int rmpi;
```




<hr>



### variable series 

```C++
Series series;
```




<hr>



### variable sources 

```C++
std::vector<Source> sources;
```




<hr>



### variable species 

```C++
std::vector<Species> species;
```




<hr>



### variable t\_store 

```C++
std::vector< PetscScalar > t_store;
```




<hr>



### variable trajx 

```C++
std::vector< PetscScalar > trajx;
```




<hr>



### variable trajy 

```C++
std::vector< PetscScalar > trajy;
```




<hr>



### variable trajz 

```C++
std::vector< PetscScalar > trajz;
```




<hr>



### variable tstep 

```C++
PetscInt tstep;
```




<hr>



### variable x 

```C++
Vec x;
```




<hr>



### variable x\_array 

```C++
Vec x_array;
```




<hr>
## Public Functions Documentation




### function MCC 

```C++
void MCC () 
```




<hr>



### function PSI\_init 

```C++
void PSI_init () 
```




<hr>



### function add\_received\_particles 

```C++
void add_received_particles () 
```




<hr>



### function cells\_init 

```C++
void cells_init () 
```




<hr>



### function collisions\_init 

```C++
void collisions_init () 
```




<hr>



### function create\_output 

```C++
void create_output () 
```




<hr>



### function decompose\_domain 

```C++
void decompose_domain () 
```




<hr>



### function field\_solver\_petsc 

```C++
void field_solver_petsc () 
```




<hr>



### function field\_source\_update 

```C++
void field_source_update () 
```




<hr>



### function grid2part 

```C++
void grid2part (
    PetscScalar,
    PetscScalar,
    PetscScalar *,
    PetscScalar *
) 
```




<hr>



### function grid\_init 

```C++
void grid_init (
    int,
    char **
) 
```




<hr>



### function initial\_particle\_load 

```C++
void initial_particle_load () 
```




<hr>



### function inject\_particle 

```C++
void inject_particle (
    int,
    PetscScalar,
    PetscScalar
) 
```




<hr>



### function inject\_particles 

```C++
void inject_particles () 
```




<hr>



### function normvel\_MaxBol 

```C++
PetscScalar normvel_MaxBol () 
```




<hr>



### function normvel\_vFv 

```C++
PetscScalar normvel_vFv () 
```




<hr>



### function num\_param\_init 

```C++
void num_param_init () 
```




<hr>



### function part2grid 

```C++
void part2grid () 
```




<hr>



### function particle\_boundaries 

```C++
void particle_boundaries () 
```




<hr>



### function particle\_boundaries\_dd 

```C++
void particle_boundaries_dd () 
```




<hr>



### function particle\_mover\_boris 

```C++
void particle_mover_boris () 
```




<hr>



### function particle\_wall\_interaction 

```C++
void particle_wall_interaction (
    int,
    int,
    int,
    int
) 
```




<hr>



### function particles\_init 

```C++
void particles_init () 
```




<hr>



### function petsc\_cleanup 

```C++
void petsc_cleanup () 
```




<hr>



### function petsc\_init 

```C++
void petsc_init (
    int,
    char **
) 
```




<hr>



### function petsc\_plot\_finalize 

```C++
void petsc_plot_finalize () 
```




<hr>



### function petsc\_plot\_init 

```C++
void petsc_plot_init () 
```




<hr>



### function petsc\_plot\_update 

```C++
void petsc_plot_update (
    int
) 
```




<hr>



### function plot\_profiling 

```C++
void plot_profiling () 
```




<hr>



### function read\_input\_file 

```C++
void read_input_file (
    char **
) 
```




<hr>



### function remove\_particle 

```C++
void remove_particle (
    int,
    int,
    int
) 
```




<hr>



### function save\_output 

```C++
void save_output () 
```




<hr>



### function send\_receive\_particles 

```C++
void send_receive_particles () 
```




<hr>



### function source 

```C++
void source () 
```




<hr>



### function source\_heating 

```C++
void source_heating () 
```




<hr>



### function source\_init 

```C++
void source_init () 
```




<hr>



### function specie\_index 

```C++
int specie_index (
    std::string
) 
```




<hr>
## Macro Definition Documentation





### define EPS0 

```C++
#define EPS0 `8.8542e-12`
```




<hr>



### define EXTERNAL 

```C++
#define EXTERNAL `extern`
```




<hr>



### define FALSE 

```C++
#define FALSE `0`
```




<hr>



### define PI 

```C++
#define PI `3.141592653589793`
```




<hr>



### define REAL 

```C++
#define REAL `PetscScalar`
```




<hr>



### define TRUE 

```C++
#define TRUE `1`
```




<hr>



### define TWOPI 

```C++
#define TWOPI `6.28318530717959`
```




<hr>



### define comm 

```C++
#define comm `MPI_COMM_WORLD`
```




<hr>



### define frand 

```C++
#define frand (
    
) `((REAL) rand() / (RAND_MAX+1.0))`
```




<hr>



### define np2c\_global 

```C++
#define np2c_global `2e6`
```




<hr>



### define plot\_3D 

```C++
#define plot_3D `FALSE`
```




<hr>



### define q\_e 

```C++
#define q_e `1.602e-19`
```




<hr>

------------------------------


