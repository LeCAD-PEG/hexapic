# Physics Validation Tests

## Running the Tests

Tests are run with `bash run_physics_validation.slm` (or `sbatch` on HPC); use `--particle-pusher`, `--field-solver`, `--particle-mesh-coupling`, or `--all` to choose groups. On first run the script creates a local `.venv` in this directory and installs dependencies from `requirements.txt` (numpy, openpmd-api); later runs reuse that venv. The HEXAPIC executable needs `libopenPMD.so` at runtime — the script sets `LD_LIBRARY_PATH` from `OPENPMD_ROOT` (default `../openpmd` relative to the repo); override with `export OPENPMD_ROOT=/path/to/openpmd` if your openPMD install is elsewhere.

---

## I. Particle Pusher Verification

*These tests verify the Boris integrator and 3V vector logic. Maxwell solver and Deposition are disabled (`FieldSolverFlag = 0`).*

*__Threshold framework:__ All thresholds are adaptive, computed as*

$$
\text{threshold} = K \times (\text{physics\_err} + \text{float\_err})
$$

*where $K=10$ is a safety factor. Physics error captures Boris integrator $O(\Delta t^2)$ phase drift, and float error captures $O(\epsilon_m \sqrt{N})$ accumulation over $N$ output steps ($\epsilon_m \approx 2.22 \times 10^{-16}$).*

### 1) Multi-Species Larmor Orbit (In-Plane)

**Description:** Verifies basic cyclotron rotation in the simulation plane (x, y) and charge sign handling using an electron–positron pair.

* **Field E:** Disabled
* **Field B:** $B=(0,\,0,\,0.1)$ T
* **Particles:**
    - **Species 1 (e−):** $q=-e,\; m=m_e,\; v_0=(10^6,\; -8794,\; 0)$ m/s
    - **Species 2 (e+):** $q=+e,\; m=m_e,\; v_0=(10^6,\; +8794,\; 0)$ m/s

    *(The non-zero $v_y$ is the Boris half-step back correction: $v_y \approx \mp v_\perp \omega_c \Delta t / 2$.)*

* **Simulation Time:** $\Delta t = 10^{-12}$ s, $T_{\text{total}}= 3 \times T_{\text{Larmor}}$. Only the last 1/3 of output steps is analyzed (initial transient skipped).

* **Metrics & Pass Criteria:**

    1. **Energy Conservation (Absolute):**

        $$
        \epsilon_E = \max|E_k(t)-E_k(0)|, \quad E_k = \tfrac{1}{2}m(v_x^2+v_y^2+v_z^2)
        $$

        **Pass if:** $\epsilon_E < K \cdot \epsilon_m \cdot \sqrt{N}$

    2. **4-Radius Separation:**

        $$
        \Delta_{\text{sep}} = \max_t \left|\mathbf{x}_{e^-}(t) - \mathbf{x}_{e^+}(t)\right|
        $$

        Expected $4 R_L$ where $R_L = \dfrac{m\, v_\perp}{|q|\, B_0}$.

        **Pass if:**

        $$
        |\Delta_{\text{sep}} - 4 R_L| < K \cdot \!\left(4 R_L \cdot \frac{20(\omega_c \Delta t)^2}{12} + 4 R_L \cdot \epsilon_m \sqrt{N}\right)
        $$

    3. **Y-Projection Trajectory:**

        $$
        \text{Extent} = \bigl(\max(y) - \min(y)\bigr) \cdot \Delta y
        $$

        Expected $2 R_L$.

        **Pass if:**

        $$
        |\text{extent} - 2 R_L| < K \cdot \!\left(2 R_L \cdot \frac{10(\omega_c \Delta t)^2}{12} + 2 R_L \cdot \epsilon_m \sqrt{N}\right)
        $$

    4. **Z-Axis Stability:**

        $$
        \Delta z = \max|z(t)-z(0)|
        $$

        **Pass if:** $\Delta z < K \cdot R_L \cdot \epsilon_m \cdot \sqrt{N}$


### 2) Multi-Species Larmor Orbit (Out-of-Plane / 3V)

**Description:** Verifies the coupling between the spatial grid (x) and the virtual velocity component ($v_z$) using a $B_y$ field. Only the electron trajectory is analyzed.

* **Field E:** Disabled
* **Field B:** $B=(0,\,0.1,\,0)$ T
* **Particles:**
    - **Species 1 (e−) & Species 2 (e+):** $v_0=(0,\;0,\;10^6)$ m/s
* **Simulation Time:** $\Delta t = 10^{-12}$ s, $T_{\text{total}}= 3 \times T_{\text{Larmor}}$. Only the last 1/3 of output steps is analyzed.

* **Metrics & Pass Criteria:**

    1. **Energy Conservation (Absolute):**

        $$
        \epsilon_E = \max|E_k(t)-E_k(0)|
        $$

        **Pass if:** $\epsilon_E < K \cdot \epsilon_m \cdot \sqrt{N}$

    2. **X-Projection Trajectory:**

        $$
        \text{Extent} = \bigl(\max(x) - \min(x)\bigr) \cdot \Delta x
        $$

        Expected $2 R_L$.

        **Pass if:**

        $$
        |\text{extent} - 2 R_L| < K \cdot \!\left(2 R_L \cdot \frac{10(\omega_c \Delta t)^2}{12} + 2 R_L \cdot \epsilon_m \sqrt{N}\right)
        $$

    3. **Y-Axis Stability:**

        $$
        \Delta y = \max|y(t)-y(0)|
        $$

        **Pass if:** $\Delta y < K \cdot R_L \cdot \epsilon_m \cdot \sqrt{N}$


### 3) E×B Drift

**Description:** Verifies the Lorentz force combination by creating a drift along the Y-axis via crossed E and B fields.

* **Field E:** From boundary potentials: $V=0$ at $x=0$, $V=100$ V at $x=L_x=0.01$ m. Gives $E_x = -10^4$ V/m.
* **Field B:** $B=(0,\,0,\,0.1)$ T
* **Particles:**
    - **Species 1 (e−):** $v_0=(0,\; 10^5,\; 0)$ m/s. Placed at domain center.
* **Simulation Time:** $\Delta t = 10^{-11}$ s, 400 steps, output every 20 steps.

* **Metrics & Pass Criteria:**

    1. **Hamiltonian Conservation:**

        $$
        H = \tfrac{1}{2}m|\mathbf{v}|^2 + q\phi
        $$

        where $\phi = -E_x \cdot x$.

        **Pass if:**

        $$
        \max|H(t) - H(0)| < K \cdot \!\left(|H_0| \cdot \frac{(\omega_c \Delta t)^2}{8} + |H_0| \cdot \epsilon_m \sqrt{N}\right)
        $$

    2. **Drift Velocity (Y):** Mean $v_y$ over the second half of the simulation, expected $v_{\text{drift}} = -E_x / B_z$.

        **Pass if:**

        $$
        |\bar{v}_y - v_{\text{drift}}| < K \cdot \!\left(|v_{\text{drift}}| \cdot \frac{(\omega_c \Delta t)^2}{8} + |v_{\text{drift}}| \cdot \epsilon_m \sqrt{N}\right)
        $$


---

## II. Field Solver Verification

*These tests verify the Poisson solver logic*

$$
\nabla^2\phi = -\frac{\rho}{\epsilon_0}
$$

*and boundary condition implementation. The solver uses an iterative method with target relative residual $\text{rtol} \sim 10^{-6}$.*

*__Threshold framework:__ Thresholds combine solver residual error ($O(\text{rtol})$), floating point accumulation ($O(\epsilon_m \sqrt{N_{\text{cells}}})$), and grid truncation error ($O(1/N)$), scaled by safety factor $K=10$.*

### 4) Mixed Boundaries (Capacitor / Dirichlet-Neumann)

**Description:** Verifies the transition between fixed potential (electrodes) and open space (dielectric). Simulates a finite parallel-plate capacitor.

* **Domain:** $32 \times 32$ grid, $L = 0.001$ m
* **Boundary Conditions:**
    - **Y-Axis (Electrodes):** $\phi(x, 0) = -10$ V, $\phi(x, L_y) = +10$ V (Dirichlet)
    - **X-Axis (Open):** $\dfrac{\partial \phi}{\partial x} = 0$ at $x=0,\, x=L_x$ (Neumann / Dielectric)
* **Source:** No particles (vacuum solve)

* **Metrics & Pass Criteria:**

    1. **Boundary Value Integrity:**

        $$
        \epsilon_{\text{BC}} = \max\!\bigl(|\phi_{y=0} - (-10)|,\; |\phi_{y=L_y} - 10|\bigr)
        $$

        **Pass if:** $\epsilon_{\text{BC}} < K \cdot \epsilon_m \cdot V_0$

    2. **Field Uniformity (Interior):** Standard deviation of $E_y$ in the interior (excluding 2-cell border near Neumann/Dirichlet corners).

        **Pass if:**

        $$
        \sigma(E_y) < K \cdot \!\left(|\bar{E}_y| \cdot \text{rtol} + |\bar{E}_y| \cdot \epsilon_m \sqrt{N_{\text{cells}}}\right)
        $$


### 5) Fully Open Boundaries (Isolated System)

**Description:** Verifies solver stability for a floating potential system using all-Neumann boundaries. Uses a neutral dipole to ensure $\int \rho \, dV = 0$.

* **Domain:** $32 \times 32$ grid, $L = 0.001$ m
* **Boundary Conditions:**
    - **All boundaries:** $\dfrac{\partial \phi}{\partial n} = 0$ (Neumann everywhere)
* **Source:** Neutral Dipole — electron cloud at $x \approx 0.35\, L$, positron cloud at $x \approx 0.65\, L$.

* **Metrics & Pass Criteria:**

    1. **Boundary Flux Leakage:** Maximum of mean $|E_n|$ across all four boundaries.

        **Pass if:**

        $$
        \max\!\bigl(\overline{|E_x|}_{\text{walls}},\; \overline{|E_y|}_{\text{walls}}\bigr) < K \cdot \!\left(E_{\text{scale}} \cdot \text{rtol} + E_{\text{scale}} \cdot \epsilon_m \sqrt{N} + \frac{E_{\text{scale}}}{N_{\min}}\right)
        $$

        where $E_{\text{scale}} = \max|\mathbf{E}|$ over the domain, $N_{\min} = \min(N_x, N_y)$.


### 6) Closed Metal Box (All-Dirichlet)

**Description:** Verifies the basic "grounded box" scenario. Tests solver matrix construction and solution symmetry.

* **Domain:** $32 \times 32$ grid, $L = 0.001$ m
* **Boundary Conditions:**
    - **All boundaries:** $\phi = 0$ (Ground)
* **Source:** Electron charge cloud centered at $(L_x/2,\, L_y/2)$.

* **Metrics & Pass Criteria:**

    1. **Boundary Zero Check:**

        $$
        \epsilon_{\text{wall}} = \max_{\text{all boundaries}} |\phi|
        $$

        **Pass if:** $\epsilon_{\text{wall}} < K \cdot \epsilon_m \cdot \max(|\phi_{\text{boundary}}|, 1)$

    2. **Symmetry:**

        $$
        \epsilon_{\text{sym}} = P_{99.9}\!\bigl(|\phi(x,y) - \phi(L_x-x,\, L_y-y)|\bigr)
        $$

        over interior points (excluding 2-cell border and 3-cell radius around center singularity).

        **Pass if:**

        $$
        \epsilon_{\text{sym}} < K \cdot \!\left(V_{\text{scale}} \cdot \text{rtol} + V_{\text{scale}} \cdot \epsilon_m \sqrt{N} + \frac{V_{\text{scale}}}{N^{1.5}}\right)
        $$


### 7) Periodic Boundaries

**Description:** Verifies matrix topology connectivity. Signals leaving the right boundary must wrap to the left boundary. Periodic in X direction with Neumann in Y.

* **Domain:** $32 \times 32$ grid, $L = 0.001$ m
* **Boundary Conditions:**
    - **X-Axis:** $\phi(0, y) = \phi(L_x, y)$ (Periodic, `PeriodicFlagX1 = 1`)
    - **Y-Axis:** $\dfrac{\partial \phi}{\partial y} = 0$ (Neumann)
* **Source:** Neutral dipole — electron and positron clouds at symmetric x-positions.

* **Metrics & Pass Criteria:**

    1. **Periodicity (Derivative Continuity):** Checks that the first derivative of $\phi$ is continuous across the periodic boundary (left–right wrap-around). Automatically detects whether the solver duplicates the endpoint ($\phi[0] = \phi[N_x]$) or not.

        **Pass if:**

        $$
        \text{derivative jump} < K \cdot \!\left(\text{max internal curvature} + V_{\text{range}} \cdot \text{rtol} + V_{\text{range}} \cdot \epsilon_m \sqrt{N}\right)
        $$


### 8) Analytical Accuracy (Green's Function)

**Description:** Quantifies numerical error by comparing the grid solution to the analytical $\ln(r)$ profile for a 2D point charge.

* **Domain:** $64 \times 64$ grid, $L = 0.001$ m
* **Boundary Conditions:** All Dirichlet $\phi = 0$.
* **Source:** Electron charge cloud at center $(L_x/2,\, L_y/2)$.
* **Comparison:** Fit $\phi_{\text{num}}(r)$ to $A \ln(r) + B$ via linear regression in the annular region $3\Delta x < r < 0.35 \min(L_x, L_y)$.

* **Metrics & Pass Criteria:**

    1. **Analytical Fit ($1 - R^2$):** Coefficient of determination for $\phi$ vs $\ln(r)$ regression. Measures how well the numerical potential follows the expected logarithmic profile.

        **Pass if:**

        $$
        1 - R^2 < K \cdot \!\left((\Delta x / r_{\min})^2 + 100 \cdot \epsilon_m\right)
        $$


---

## III. Particle-Mesh Coupling Verification

*These tests verify the projection of particle quantities to the grid (Deposition) and the sampling of grid fields to the particle (Interpolation/Gathering).*

*__Threshold framework:__ All thresholds are adaptive, computed as*

$$
\text{threshold} = K \times (\text{physics\_err} + \text{float\_err})
$$

*where $K=10$ is a safety factor. Physics error captures solver residual ($O(\text{rtol})$), discretization ($O(\Delta x^2)$), and CIC shape-function effects. Float error captures $O(\epsilon_m \sqrt{N})$ accumulation ($\epsilon_m \approx 2.22 \times 10^{-16}$). Grid charge is inferred via Gauss's law:*

$$
Q_{\text{grid}} = -\epsilon_0 \sum_{\text{interior}} (\nabla^2 \phi)_{i,j}\, \Delta x \, \Delta y
$$


### 9) Charge Deposition & Potential Topology

**Description:** Verifies that deposited charge creates a quantitatively correct and symmetric potential well. Ensures coupling between CIC particle weighting and the Poisson solver, and that physical constants ($\epsilon_0$) are correctly applied.

* **Domain:** $32 \times 32$ grid, $L = 0.001$ m, All-Dirichlet ($\phi = 0$).
* **Source:** Charge cloud at center $(L_x/2,\, L_y/2)$. Total charge $Q_{\text{total}} = -q_e \cdot n \cdot \Delta x \cdot \Delta y$ where $n = 2.1 \times 10^{15}$ m$^{-3}$ (from `static_deposit.inp`).
* **Simulation:** 5 steps, $\Delta t = 10^{-12}$ s. Field read from final iteration.

* **Metrics & Pass Criteria:**

    1. **Integrated Charge Conservation:**

        $$
        \epsilon_Q = |Q_{\text{grid}} - Q_{\text{total}}|
        $$

        **Pass if:**

        $$
        \epsilon_Q < K \cdot |Q_{\text{total}}| \cdot \!\left(\text{rtol} + \frac{1}{N} + \epsilon_m \sqrt{N_{\text{cells}}}\right)
        $$

        *(Combines solver residual, $O(1/N)$ Laplacian discretization error, and float accumulation.)*

    2. **Symmetry:**

        $$
        \epsilon_{\text{sym}} = P_{99.9}\!\bigl(|\phi(x,y) - \phi(L_x-x,\, L_y-y)|\bigr)
        $$

        over interior points (excluding 2-cell border and 3-cell radius around center singularity).

        **Pass if:**

        $$
        \epsilon_{\text{sym}} < K \cdot V_{\text{scale}} \cdot \!\left(\text{rtol} + \epsilon_m \sqrt{N_{\text{pts}}} + \frac{1}{N_{\text{grid}}}\right)
        $$

        *(Combines solver residual, float accumulation, and $O(1/N)$ grid discretization asymmetry from CIC and Laplacian.)*

    3. **L2 Profile Error:** Relative RMS between simulated $\phi$ and analytic Green's function $\phi_{\text{ana}}$ (2D Fourier series for point charge in grounded box, 50 modes).

        **Pass if:**

        $$
        \sqrt{\frac{\sum (\phi - \phi_{\text{ana}})^2}{\sum \phi_{\text{ana}}^2}} < K \cdot \!\left(\text{rtol} + \frac{1}{N_{\text{grid}}} + \epsilon_m \sqrt{N_{\text{pts}}}\right)
        $$

        *(Combines solver residual, $O(1/N)$ Fourier truncation and grid discretization error, and float accumulation.)*


### 10) Field Interpolation (Acceleration)

**Description:** Verifies that a particle in a uniform E-field experiences the correct Lorentz acceleration $a = qE/m$ by measuring velocity change over multiple timesteps.

* **Domain:** $32 \times 32$ grid, $L_y$ from input file (default $L = 0.001$ m).
* **Boundary Conditions:**
    - **Y-Axis (Electrodes):** $\phi(x, 0) = -10$ V, $\phi(x, L_y) = +10$ V (Dirichlet)
    - **X-Axis (Open):** $\dfrac{\partial \phi}{\partial x} = 0$ (Neumann / Dielectric)
* **Particle:** Single electron ($q = -e$, $m = m_e$) at center, $v_0 = 0$.
* **Simulation:** 20 steps, $\Delta t = 10^{-12}$ s, output every step.

* **Metrics & Pass Criteria:**

    1. **Acceleration Accuracy:** Linear fit $v_y(t) = a_{\text{meas}} \cdot t + v_0$ gives measured acceleration.

        Expected: $a_{\text{expected}} = \dfrac{q\, E_y}{m}$ where $E_y = -(V_{\text{top}} - V_{\text{bot}}) / L_y$.

        **Pass if:**

        $$
        |a_{\text{meas}} - a_{\text{expected}}| < K \cdot |a_{\text{expected}}| \cdot \!\left(\text{rtol} + \frac{1}{N_{\text{grid}}} + \epsilon_m \sqrt{N_{\text{pts}}}\right)
        $$

        *(Combines solver residual, $O(1/N)$ grid error from CIC and finite-difference E, and float accumulation in the linear fit.)*


### 11) Self-Force Stability (Stationary Particle)

**Description:** Integration test for the full PIC loop: Deposit $\rho$ → Solve Poisson → Interpolate $\mathbf{E}$ → Push. A stationary particle at the center of a symmetric grounded box should feel near-zero net force; any displacement is parasitic (self-force).

* **Domain:** $32 \times 32$ grid, $L = 0.001$ m, All-Dirichlet ($\phi = 0$).
* **Particle:** Single electron at rest ($v_0 = 0$) centered at $(L_x/2,\, L_y/2)$.
* **Simulation:** 20 steps, $\Delta t = 10^{-12}$ s, output every step.

* **Metrics & Pass Criteria:**

    1. **Parasitic Displacement:**

        $$
        \Delta r = \sqrt{(x_{\text{final}} - x_0)^2 + (y_{\text{final}} - y_0)^2}
        $$

        **Pass if:**

        $$
        \Delta r < K \cdot \!\left(N_{\text{steps}} \cdot \Delta x + \epsilon_m \cdot \Delta x \cdot \sqrt{N_{\text{steps}}}\right)
        $$

        *(PIC codes exhibit parasitic motion from CIC shape asymmetry, field gradient noise, and pusher-solver coupling. Threshold is grid-based: displacement bounded by $O(N \cdot \Delta x)$ over the simulation.)*


### 12) Current Deposition (Charge Conservation)

**Description:** Verifies that moving a particle and re-depositing charge conserves total charge on the grid. This ensures Gauss's law ($\nabla \cdot \mathbf{E} = \rho / \epsilon_0$) holds at every timestep, critical for Esirkepov-type current deposition schemes.

* **Domain:** $32 \times 32$ grid, $L = 0.001$ m, All-Dirichlet ($\phi = 0$).
* **Particle:** Single electron with $v_x = 10^4$ m/s (traverses ${\sim}\,0.2$ cells over 10 steps).
* **Simulation:** 10 steps, $\Delta t = 10^{-12}$ s, output every step. Poisson solver active (`FieldSolverFlag = 1`).

* **Metrics & Pass Criteria:**

    1. **Charge Conservation (Gauss):**

        $$
        \Delta Q = |Q_{\text{grid}}(t_{\text{final}}) - Q_{\text{grid}}(t_0)|
        $$

        where $Q_{\text{grid}}$ is computed via discrete Laplacian of $\phi$.

        **Pass if:**

        $$
        \Delta Q < K \cdot |Q| \cdot \!\left(\text{rtol} + \epsilon_m \sqrt{N_{\text{cells}}} + \frac{1}{N}\right)
        $$

        *(Combines solver residual at two iterations, float accumulation, and $O(1/N)$ CIC spreading discretization error.)*


---

## Output

```bash
ROOT_DIR: /xxxxx/hexapic
HEXAPIC:  /xxxxx/hexapic/source/HEXAPIC
NPROCS:   1

=======================================================================================================
                                    HEXAPIC PHYSICS VALIDATION SUITE                                   
                 Date: 2026-02-18 14:25 | Host: xxxxxx | Procs: 1 | Commit: 8fc15a8                 
=======================================================================================================

-------------------------------------------------------------------------------------------------------
                                         PARTICLE PUSHER TESTS                                         
-------------------------------------------------------------------------------------------------------

[TEST 01] In-Plane Larmor Orbit
   Context: B=0.1T, v=1e+06 m/s, dt=1ps, Steps=1071
   ----------------------------------------------------------------------------------------------------
   METRIC                              MEASURED        EXPECTED         ERROR       MAX ERR   STATUS
   ----------------------------------------------------------------------------------------------------
   Energy Conservation                 1.16e-33        0.00e+00      1.16e-33      4.20e-14   [PASS]
   4R Separation                     2.2802e-04      2.2743e-04    5.9871e-07    1.1725e-06   [PASS]
   Y-projection trajectory           1.1373e-04      1.1371e-04    1.9311e-08    2.9314e-07   [PASS]
   Z-Axis Drift                        0.00e+00        0.00e+00      0.00e+00      2.39e-18   [PASS]
   ----------------------------------------------------------------------------------------------------
   Result: OK  Time: 7.62s

[TEST 02] 3V Larmor Orbit
   Context: B_y=0.1T, v_z=1e+06 m/s, 3D Rotation Check
   ----------------------------------------------------------------------------------------------------
   METRIC                              MEASURED        EXPECTED         ERROR       MAX ERR   STATUS
   ----------------------------------------------------------------------------------------------------
   Energy Conservation                 6.74e-34        0.00e+00      6.74e-34      1.88e-14   [PASS]
   X-Projection Trajectory           1.1370e-04      1.1371e-04    1.0182e-08    2.9314e-07   [PASS]
   Y-Axis Drift                        0.00e+00        0.00e+00      0.00e+00      1.07e-18   [PASS]
   ----------------------------------------------------------------------------------------------------
   Result: OK  Time: 2.64s

[TEST 03] ExB Drift
   Context: E=-1e+04 V/m, B=0.1T, Expected v_drift=1e+05 m/s
   ----------------------------------------------------------------------------------------------------
   METRIC                              MEASURED        EXPECTED         ERROR       MAX ERR   STATUS
   ----------------------------------------------------------------------------------------------------
   Hamiltonian Conservation            2.10e-22        0.00e+00      2.10e-22      1.87e-19   [PASS]
   Drift Velocity                    1.0023e+05      1.0000e+05    2.3437e+02    3.8668e+03   [PASS]
   ----------------------------------------------------------------------------------------------------
   Result: OK  Time: 2.29s

-------------------------------------------------------------------------------------------------------
                                           FIELD SOLVER TESTS                                          
-------------------------------------------------------------------------------------------------------

[TEST 04] Mixed Boundaries
   Context: Dirichlet Y, Neumann X
   ----------------------------------------------------------------------------------------------------
   METRIC                              MEASURED        EXPECTED         ERROR       MAX ERR   STATUS
   ----------------------------------------------------------------------------------------------------
   Boundary Integrity                  0.00e+00        0.00e+00      0.00e+00      2.22e-14   [PASS]
   Field Uniformity                  3.9426e-02        0.00e+00    3.9426e-02    2.0000e-01   [PASS]
   ----------------------------------------------------------------------------------------------------
   Result: OK  Time: 0.01s

[TEST 05] Fully Open Boundaries
   Context: All Neumann, Neutral Dipole
   ----------------------------------------------------------------------------------------------------
   METRIC                              MEASURED        EXPECTED         ERROR       MAX ERR   STATUS
   ----------------------------------------------------------------------------------------------------
   Boundary Flux Leakage             4.4524e+02        0.00e+00    4.4524e+02    1.5967e+04   [PASS]
   ----------------------------------------------------------------------------------------------------
   Result: OK  Time: 0.01s

[TEST 06] Closed Metal Box
   Context: All phi=0, Centered Charge
   ----------------------------------------------------------------------------------------------------
   METRIC                              MEASURED        EXPECTED         ERROR       MAX ERR   STATUS
   ----------------------------------------------------------------------------------------------------
   Boundary Zero Check                 0.00e+00        0.00e+00      0.00e+00      2.22e-15   [PASS]
   Symmetry Error                    1.8462e-02        0.00e+00    1.8462e-02    6.0076e-01   [PASS]
   ----------------------------------------------------------------------------------------------------
   Result: OK  Time: 0.02s

[TEST 07] Periodic X Boundaries
   Context: Periodic X, Neumann Y
   ----------------------------------------------------------------------------------------------------
   METRIC                              MEASURED        EXPECTED         ERROR       MAX ERR   STATUS
   ----------------------------------------------------------------------------------------------------
   Periodicity                       1.2711e+02        0.00e+00    1.2711e+02    6.3732e+02   [PASS]
   ----------------------------------------------------------------------------------------------------
   Result: OK  Time: 0.01s

[TEST 08] Analytical Accuracy
   Context: Line Charge vs ln
   ----------------------------------------------------------------------------------------------------
   METRIC                              MEASURED        EXPECTED         ERROR       MAX ERR   STATUS
   ----------------------------------------------------------------------------------------------------
   Analytical Fit                    2.0674e-04        0.00e+00    2.0674e-04    1.1111e+00   [PASS]
   ----------------------------------------------------------------------------------------------------
   Result: OK  Time: 0.01s

-------------------------------------------------------------------------------------------------------
                                      PARTICLE MESH COUPLING TESTS                                     
-------------------------------------------------------------------------------------------------------

[TEST 09] Charge Deposition & Potential Topology
   Context: 32x32, L=0.001 m, grounded box, Q_total=-3.29e-13 C
   ----------------------------------------------------------------------------------------------------
   METRIC                              MEASURED        EXPECTED         ERROR       MAX ERR   STATUS
   ----------------------------------------------------------------------------------------------------
   Integrated Charge Conserv...        8.18e-15        0.00e+00      8.18e-15      9.96e-14   [PASS]
   Symmetry                          1.2351e-03        0.00e+00    1.2351e-03    3.2676e-03   [PASS]
   L2 Profile Error                  4.3836e-02        0.00e+00    4.3836e-02    3.0304e-01   [PASS]
   ----------------------------------------------------------------------------------------------------
   Result: OK  Time: 2.11s

[TEST 10] Field Interpolation
   Context: Capacitor dV=20V, L=0.001m
   ----------------------------------------------------------------------------------------------------
   METRIC                              MEASURED        EXPECTED         ERROR       MAX ERR   STATUS
   ----------------------------------------------------------------------------------------------------
   Acceleration Accuracy             4.3808e+13        0.00e+00    4.3808e+13    1.0993e+15   [PASS]
   ----------------------------------------------------------------------------------------------------
   Result: OK  Time: 2.06s

[TEST 11] Self-Force Stability
   Context: Stationary particle, Centered
   ----------------------------------------------------------------------------------------------------
   METRIC                              MEASURED        EXPECTED         ERROR       MAX ERR   STATUS
   ----------------------------------------------------------------------------------------------------
   Parasitic Displacement            1.6419e-04        0.00e+00    1.6419e-04    6.2500e-03   [PASS]
   ----------------------------------------------------------------------------------------------------
   Result: OK  Time: 2.06s

[TEST 12] Current Deposition
   Context: Moving particle, Gauss law charge check
   ----------------------------------------------------------------------------------------------------
   METRIC                              MEASURED        EXPECTED         ERROR       MAX ERR   STATUS
   ----------------------------------------------------------------------------------------------------
   Charge Conservation                 3.02e-17        0.00e+00      3.02e-17      9.71e-14   [PASS]
   ----------------------------------------------------------------------------------------------------
   Result: OK  Time: 2.04s

=======================================================================================================
                                                SUMMARY                                                
=======================================================================================================
PARTICLE PUSHER
  1. In-Plane Larmor Orbit               [PASS]
  2. 3V Larmor Orbit                     [PASS]
  3. ExB Drift                           [PASS]

FIELD SOLVER
  4. Mixed Boundaries                    [PASS]
  5. Fully Open Boundaries               [PASS]
  6. Closed Metal Box                    [PASS]
  7. Periodic X Boundaries               [PASS]
  8. Analytical Accuracy                 [PASS]

PARTICLE MESH COUPLING
  9. Charge Deposition & Potential ...   [PASS]
  10. Field Interpolation                 [PASS]
  11. Self-Force Stability                [PASS]
  12. Current Deposition                  [PASS]

-------------------------------------------------------------------------------------------------------
OVERALL: PASSED (12/12 tests passed, 22/22 metrics)
Total Runtime: 20.89s
=======================================================================================================


```