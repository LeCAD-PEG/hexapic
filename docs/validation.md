# vs. XOOPIC

## Two-stream instability

* **Description**: Formation of two-stream instability when counter-propagating electron beams have velocities much larger than their thermal velocity, in a background of cold ions.

* **Input file:** original XOOPIC `two_stream_ee_es.inp` with modified grid size 256$\times$8 (from 256$\times$2)

<video controls loop muted width="100%">
  <source src="../videos/hexapic_vs_xoopic_twostream.mp4" type="video/mp4">
</video>

<p></p>

## Electron Beam Spreading

* **Description**: Spreading of an electron beam when propagating along metallic walls.

* **Input file:** original XOOPIC `xybeam.inp` with modified grid size 48$\times$24 (from 20$\times$20)

<video controls loop muted width="100%">
  <source src="../videos/hexapic_vs_xoopic_xybeam.mp4" type="video/mp4">
</video>

<p></p>

## Deposition

* **Description**: Ion deposition on a substrate surface, including a trench.

* **Input file:** original XOOPIC `deposition.inp` with the following modifications:
	- grid size from 30$\times$100 to 32$\times$64 but keeping the same spatial resolution `dx=dy=5e-8m`
	- time step from `1e-13s` to `0.5e-13s` (better time resolution)
	- `np2c` from `1e4` to `1e3` (10$\times$ simulated particles for better statistics) 

<video controls loop muted width="100%">
  <source src="../videos/hexapic_vs_xoopic_deposition.mp4" type="video/mp4">
</video>