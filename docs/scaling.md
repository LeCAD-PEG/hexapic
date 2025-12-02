# Scaling v0.8.0

## VEGA HPC

128 cores/node.

See the machine website for more specifications: https://doc.vega.izum.si/general-spec/.

### Strong scaling

* Mesh size 2048x2048, single node RAM usage ~180 GB.

<p align="center">
  <img src="../images/strong_scaling.png" alt="strong_scaling" width="100%">
</p>

<p align="center">
  <img src="../images/speedup.png" alt="speedup" width="100%">
</p>


### Weak scaling

* Single node mesh size 2048x2048, single node RAM usage ~180 GB. 
* Mesh size increase 2048\*nodes, e.g.:
	- 1 node  : 2048x2048
	- 2 nodes : 2048x4096
	- 4 nodes : 4096x4096
	- etc...

<p align="center">
  <img src="../images/weak_scaling.png" alt="weak_scaling" width="100%">
</p>


<p align="center">
  <img src="../images/parallel_efficiency.png" alt="parallel_efficiency" width="100%">
</p>