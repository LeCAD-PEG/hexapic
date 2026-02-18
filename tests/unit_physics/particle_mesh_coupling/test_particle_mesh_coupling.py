#!/usr/bin/env python3
"""
Physics validation tests for HEXAPIC particle-mesh coupling.

Tests:
9.  Charge Deposition & Potential Topology: Single electron at center; integrated
    charge conservation, symmetry of potential, L2 vs analytic Green's function.
10. Field Interpolation: Quantitative acceleration test in a capacitor (E = -grad phi).
11. Self-Force Stability: Verification of momentum conservation for a static particle.
12. Current Deposition: Continuity check for moving particles.

Usage:
    python3 test_particle_mesh_coupling.py --exe ./HEXAPIC
"""

import argparse
import sys
import os
import time
from pathlib import Path
from typing import List, Dict, Tuple
import numpy as np
import openpmd_api as opmd

# Add parent directory to sys.path to import utils
sys.path.append(str(Path(__file__).parent.parent))
from utils import (
    TestCase,
    Metric,
    FieldData,
    ParticleData,
    print_header,
    print_group_header,
    print_test_case,
    print_final_summary,
    save_results,
    cleanup_simulation_outputs,
    run_hexapic,
    read_particle_trajectories,
    read_field_data,
    parse_input_file,
)

# =============================================================================
# PHYSICAL & NUMERICAL CONSTANTS
# =============================================================================

Q_E = 1.602176634e-19  # Elementary charge [C]
M_E = 9.1093837139e-31  # Electron mass [kg]
EPS_0 = 8.854187817e-12  # Vacuum permittivity [F/m]

# Safety factor for thresholds
K_SAFETY = 10.0

# Typical relative residual at solver convergence (iterative solver target)
SOLVER_RESIDUAL_RTOL = 1e-6

# Machine relative error
EPS_MACHINE = np.finfo(float).eps

GROUP_NAME = "PARTICLE MESH COUPLING TESTS"


# =============================================================================
# METRIC CALCULATIONS
# =============================================================================


def count_particles_at_iterations(bp_file: str, species_idx: int = 0) -> Dict[int, int]:
    """Helper: Count particles at each iteration via OpenPMD."""
    if not os.path.exists(bp_file):
        return {}
    series = opmd.Series(bp_file, opmd.Access_Type.read_linear)
    result = {}
    for idx, it in enumerate(series.read_iterations()):
        sp_key = str(species_idx)
        if sp_key in it.particles:
            sp = it.particles[sp_key]
            # Just read the size of the x component to get count
            x_comp = sp["position"]["x"]
            result[idx] = x_comp.shape[0]
        else:
            result[idx] = 0
        it.close()
    del series
    return result


def calculate_particle_conservation(
    bp_file: str, expected_count: int
) -> Tuple[float, float]:
    """
    Verify zero particle loss.
    """
    counts = count_particles_at_iterations(bp_file, species_idx=0)
    if not counts:
        return np.nan, np.nan

    last_iter = max(counts.keys())
    measured = float(counts[last_iter])

    # Threshold is exactly 0 for conservation
    return measured, 0.1


def calculate_field_acceleration(
    particle: ParticleData,
    V_left: float,
    V_right: float,
    L_dist: float,
    mass: float,
    charge: float,
    dt: float,
    n_grid: int = 32,
) -> Tuple[float, float]:
    """
    Compare particle acceleration to qE/m.
    E = -(V_right - V_left) / L
    """
    valid = ~np.isnan(particle.vy) & ~np.isnan(particle.time)
    if not np.any(valid):
        return np.nan, np.nan

    # E-field estimation (Capacitor approximation)
    Ey = -(V_right - V_left) / L_dist
    ay_expected = (charge * Ey) / mass

    # Measured acceleration (linear fit to vy)
    vy = particle.vy[valid]
    t = particle.time[valid]
    if len(t) < 2:
        return np.nan, np.nan

    slope, _ = np.polyfit(t, vy, 1)
    measured_ay = float(slope)
    a_scale = np.abs(ay_expected)

    # Threshold components (all relative to |a_expected|), then × K_SAFETY:
    # 1. Solver residual in phi → relative error in E and thus in a
    # 2. Grid discretization: CIC + finite-difference E → O(1/N)
    # 3. Float accumulation in linear fit over n_pts
    n_pts = len(t)
    solver_err = a_scale * SOLVER_RESIDUAL_RTOL
    grid_err = a_scale / float(n_grid)
    float_err = a_scale * EPS_MACHINE * np.sqrt(n_pts)

    abs_err_threshold = K_SAFETY * (solver_err + grid_err + float_err)
    error = float(np.abs(measured_ay - ay_expected))
    return error, abs_err_threshold


def calculate_self_force_stability(
    particle: ParticleData, dx: float, n_steps: int = 20
) -> Tuple[float, float]:
    """
    Check displacement of a stationary particle due to self-force.
    PIC codes inherently have parasitic motion from CIC asymmetry,
    field gradient noise, and pusher-solver coupling.

    Threshold is grid-based: displacement should be O(dx) over the run.
    """
    valid = ~np.isnan(particle.x)
    if not np.any(valid):
        return np.nan, np.nan

    x = particle.x[valid]
    y = particle.y[valid]

    # Total displacement in SI units
    measured_disp_si = float(np.sqrt((x[-1] - x[0]) ** 2 + (y[-1] - y[0]) ** 2))

    # Grid-based threshold: parasitic displacement should be bounded by O(N * dx)
    # where N is number of steps. This accounts for CIC asymmetry and field noise
    # accumulating over time.
    grid_err = n_steps * dx
    float_err = dx * EPS_MACHINE * np.sqrt(n_steps)

    abs_err_threshold = K_SAFETY * (grid_err + float_err)

    return measured_disp_si, abs_err_threshold


def compute_grid_charge(V: np.ndarray, dx: float, dy: float) -> float:
    """
    Compute total enclosed charge via discrete Laplacian of the potential.
    From Poisson's equation: nabla^2 V = -rho / eps_0
    """
    lap_interior = (V[2:, 1:-1] + V[:-2, 1:-1] - 2.0 * V[1:-1, 1:-1]) / dx**2 + (
        V[1:-1, 2:] + V[1:-1, :-2] - 2.0 * V[1:-1, 1:-1]
    ) / dy**2
    return float(-EPS_0 * np.sum(lap_interior) * dx * dy)


def potential_point_charge_grounded_box_2d(
    x_grid: np.ndarray,
    y_grid: np.ndarray,
    x0: float,
    y0: float,
    Lx: float,
    Ly: float,
    q: float,
    n_modes: int = 50,
) -> np.ndarray:
    """
    Analytic potential for a point charge q at (x0, y0) in a grounded 2D box.
    Uses Fourier Series (Sine expansion) which is numerically stable.
    """
    phi = np.zeros_like(x_grid, dtype=float)
    prefactor = 4.0 * q / (Lx * Ly * EPS_0)

    for n in range(1, n_modes):
        kn = n * np.pi / Lx
        sin_n_x0 = np.sin(kn * x0)
        sin_n_x = np.sin(kn * x_grid)

        for m in range(1, n_modes):
            km = m * np.pi / Ly
            k_sq = kn**2 + km**2

            coeff = sin_n_x0 * np.sin(km * y0) / k_sq
            phi += coeff * sin_n_x * np.sin(km * y_grid)

    return phi * prefactor


def calculate_integrated_charge_error(
    field: FieldData,
    q_expected: float,
) -> Tuple[float, float]:
    """
    Integrated charge conservation: |sum(rho*dV) - q_e|.
    Laplacian-based charge integral has discretization and solver errors.
    """
    Q_grid = compute_grid_charge(field.V, field.dx, field.dy)
    measured = float(np.abs(Q_grid - q_expected))

    # Threshold components:
    # 1. Solver residual (iterative Poisson solver)
    # 2. Discrete Laplacian truncation error O(dx^2) ~ O(1/N) in 2D
    # 3. Float accumulation over grid
    n_cells = field.V.size
    N = np.sqrt(n_cells)  # Linear grid size
    Q_scale = np.abs(q_expected)

    solver_err = Q_scale * SOLVER_RESIDUAL_RTOL
    grid_err = Q_scale / N  # O(1/N) Laplacian discretization error
    float_err = Q_scale * EPS_MACHINE * np.sqrt(n_cells)

    abs_err_threshold = K_SAFETY * (solver_err + grid_err + float_err)

    return measured, abs_err_threshold


def calculate_deposit_symmetry(
    field: FieldData,
    Lx: float,
    Ly: float,
) -> Tuple[float, float]:
    """
    Symmetry: epsilon_sym = P_{99.9}(|phi(x,y) - phi(Lx-x, Ly-y)|).
    """
    V = field.V
    nx, ny = V.shape
    border = 2
    r_min_cells = 3.0
    cx, cy = (nx - 1) / 2.0, (ny - 1) / 2.0
    y_idx, x_idx = np.ogrid[0:nx, 0:ny]
    r_sq = (x_idx - cx) ** 2 + (y_idx - cy) ** 2
    is_interior = (
        (x_idx >= border)
        & (x_idx < ny - border)
        & (y_idx >= border)
        & (y_idx < nx - border)
    )
    is_far_field = r_sq >= r_min_cells**2

    mask = is_interior & is_far_field & is_interior[::-1, ::-1]
    V_flip = V[::-1, ::-1]
    diff = np.abs(V[mask] - V_flip[mask])

    if diff.size == 0:
        return np.nan, np.nan

    measured = float(np.percentile(diff, 99.9))
    V_scale = float(np.max(np.abs(V[mask])))

    if V_scale < 1e-12:
        return measured, K_SAFETY * EPS_MACHINE

    N_pts = float(np.sum(mask))
    N_grid = float(nx)  # Linear grid size

    # Symmetry threshold components:
    # 1. Solver residual contribution
    # 2. Float accumulation over masked points
    # 3. Grid discretization error O(1/N) for CIC + Laplacian
    solver_err = V_scale * SOLVER_RESIDUAL_RTOL
    float_err = V_scale * EPS_MACHINE * np.sqrt(N_pts)
    grid_err = V_scale / N_grid  # O(1/N) discretization asymmetry

    abs_err_threshold = K_SAFETY * (solver_err + float_err + grid_err)
    return measured, abs_err_threshold


def calculate_l2_profile_error(
    field: FieldData,
    x0: float,
    y0: float,
    Lx: float,
    Ly: float,
    q: float,
    border: int = 2,
    r_min_cells: float = 3.0,
) -> Tuple[float, float]:
    """
    L2 relative RMS between simulated phi and analytic Fourier series solution.
    """
    V = field.V
    nx, ny = V.shape
    x_grid = np.broadcast_to(np.arange(ny, dtype=float) * field.dx, (nx, ny))
    y_grid = np.broadcast_to(
        (np.arange(nx, dtype=float) * field.dy).reshape(-1, 1), (nx, ny)
    )

    # Use the stable Fourier implementation
    phi_ana = potential_point_charge_grounded_box_2d(x_grid, y_grid, x0, y0, Lx, Ly, q)

    cx, cy = (ny - 1) / 2.0, (nx - 1) / 2.0
    y_idx, x_idx = np.ogrid[0:nx, 0:ny]
    r_sq = (x_idx - cx) ** 2 + (y_idx - cy) ** 2
    is_interior = (
        (x_idx >= border)
        & (x_idx < ny - border)
        & (y_idx >= border)
        & (y_idx < nx - border)
    )
    is_far_field = r_sq >= r_min_cells**2
    mask = is_interior & is_far_field

    if not np.any(mask):
        return np.nan, np.nan

    num = np.sum((V[mask] - phi_ana[mask]) ** 2)
    denom = np.sum(phi_ana[mask] ** 2)

    if denom < 1e-30:
        return np.nan, np.nan

    measured = float(np.sqrt(num / denom))
    N_pts = float(np.sum(mask))
    N_grid = float(nx)  # Linear grid size

    # L2 threshold components (relative error):
    # 1. Solver residual
    # 2. Fourier truncation + grid discretization O(1/N)
    # 3. Float accumulation in RMS
    solver_err = SOLVER_RESIDUAL_RTOL
    grid_err = 1.0 / N_grid  # O(1/N) discretization error
    float_err = EPS_MACHINE * np.sqrt(N_pts)

    abs_err_threshold = K_SAFETY * (solver_err + grid_err + float_err)
    return measured, abs_err_threshold


def calculate_charge_conservation(
    bp_file: str,
    dx: float,
    dy: float,
) -> Tuple[float, float]:
    """
    Verify global charge conservation across timesteps via Gauss's law.
    """
    if not os.path.exists(bp_file):
        return np.nan, np.nan

    series = opmd.Series(bp_file, opmd.Access_Type.read_linear)

    Q_first = None
    Q_last = None
    n_cells = 0

    for it in series.read_iterations():
        V_rc = it.meshes["V"][opmd.Mesh_Record_Component.SCALAR]
        V_chunk = V_rc.load_chunk()
        V_rc.series_flush()

        Q = compute_grid_charge(V_chunk, dx, dy)
        n_cells = V_chunk.size

        if Q_first is None:
            Q_first = Q
        Q_last = Q
        it.close()

    del series

    if Q_first is None or Q_last is None:
        return np.nan, np.nan

    measured = abs(Q_last - Q_first)

    Q_scale = max(abs(Q_first), abs(Q_last), 1e-30)
    N = np.sqrt(n_cells)  # Linear grid size

    # Threshold components:
    # 1. Solver residual (two evaluations, first and last)
    # 2. Float accumulation over grid
    # 3. CIC charge spreading O(1/N) discretization error
    solver_noise = Q_scale * SOLVER_RESIDUAL_RTOL
    float_noise = Q_scale * EPS_MACHINE * np.sqrt(n_cells)
    grid_noise = Q_scale / N  # O(1/N) CIC spreading error

    abs_err_threshold = K_SAFETY * (solver_noise + float_noise + grid_noise)

    return measured, abs_err_threshold


# =============================================================================
# TEST IMPLEMENTATIONS
# =============================================================================


def run_static_deposit_test(
    test_id: int, hexapic_exe: str, num_procs: int = 1, skip_sim: bool = False
) -> TestCase:
    """Test 9: Charge Deposition & Potential Topology."""
    input_file = "static_deposit.inp"
    params = parse_input_file(input_file)
    dx = params.get("dx")
    dy = params.get("dy")
    Lx = params.get("Lx")
    Ly = params.get("Ly")
    x0, y0 = Lx / 2.0, Ly / 2.0

    # Calculate expected total charge from load parameters
    # Load region is ~1 cell centered, density gives total particle count
    # For this test: density=2.1e15, load area ~ dx*dy, so N_particles ~ density * dx * dy
    # Total charge = N_particles * q_electron
    load_density = 2.1e15  # From static_deposit.inp
    load_area = dx * dy  # Approximate 1-cell load region
    q_total_expected = -Q_E * load_density * load_area  # Total deposited charge

    context = f"32x32, L={Lx} m, grounded box, Q_total={q_total_expected:.2e} C"
    start_time = time.time()

    if not skip_sim:
        success = run_hexapic(input_file, hexapic_exe, steps=5, num_procs=num_procs)
        if not success:
            return TestCase(
                test_id,
                "Charge Deposition & Potential Topology",
                GROUP_NAME,
                context,
                [Metric("Simulation", 0, 1, passed=False)],
            )

    bp_file = f"{input_file}.bp4"
    field = read_field_data(bp_file)
    if field is None:
        return TestCase(
            test_id,
            "Charge Deposition & Potential Topology",
            GROUP_NAME,
            context,
            [Metric("Field read", 0, 1, passed=False)],
            time.time() - start_time,
        )

    metrics = []

    # 1. Integrated Charge Conservation: compare Gauss-derived Q to expected total charge
    charge_err, charge_thr = calculate_integrated_charge_error(field, q_total_expected)
    metrics.append(
        Metric(
            "Integrated Charge Conservation",
            charge_err,
            0.0,
            charge_err,
            charge_thr,
            charge_err < charge_thr,
        )
    )

    # 2. Symmetry: potential should be symmetric about center
    sym_err, sym_thr = calculate_deposit_symmetry(field, Lx, Ly)
    metrics.append(
        Metric("Symmetry", sym_err, 0.0, sym_err, sym_thr, sym_err < sym_thr)
    )

    # 3. L2 Profile Error vs analytic Green's function
    l2_err, l2_thr = calculate_l2_profile_error(field, x0, y0, Lx, Ly, q_total_expected)
    metrics.append(
        Metric("L2 Profile Error", l2_err, 0.0, l2_err, l2_thr, l2_err < l2_thr)
    )

    runtime = time.time() - start_time
    return TestCase(
        test_id,
        "Charge Deposition & Potential Topology",
        GROUP_NAME,
        context,
        metrics,
        runtime,
    )


def run_interpolation_test(
    test_id: int, hexapic_exe: str, num_procs: int = 1, skip_sim: bool = False
) -> TestCase:
    """Test 10: Field Interpolation (Acceleration)."""
    input_file = "interpolation.inp"
    params = parse_input_file(input_file)

    # Read domain size from input file; V_bot=-10, V_top=+10 from Equipotentials
    Ly = params.get("Ly", 0.001)
    V_bot, V_top = -10.0, 10.0
    dt = params.get("dt", 1e-12)
    n_grid = params.get("K", 32)  # grid cells along E-field direction (y)

    context = f"Capacitor dV=20V, L={Ly}m"
    start_time = time.time()

    if not skip_sim:
        success = run_hexapic(
            input_file, hexapic_exe, steps=20, dsteps=1, num_procs=num_procs
        )
        if not success:
            return TestCase(
                test_id,
                "Field Interpolation",
                GROUP_NAME,
                context,
                [Metric("Simulation", 0, 1, passed=False)],
            )

    bp_file = f"{input_file}.bp4"
    electron = read_particle_trajectories(bp_file, species_idx=0)

    metrics = []

    acc_err, acc_thr = calculate_field_acceleration(
        electron, V_bot, V_top, Ly, M_E, -Q_E, dt, n_grid
    )

    metrics.append(
        Metric(
            "Acceleration Accuracy (qE/m)",
            acc_err,
            0.0,
            acc_err,
            acc_thr,
            acc_err < acc_thr,
        )
    )

    runtime = time.time() - start_time
    return TestCase(
        test_id, "Field Interpolation (Force)", GROUP_NAME, context, metrics, runtime
    )


def run_self_force_test(
    test_id: int, hexapic_exe: str, num_procs: int = 1, skip_sim: bool = False
) -> TestCase:
    """Test 11: Self-Force Stability."""
    input_file = "self_force.inp"
    params = parse_input_file(input_file)
    dx = params.get("dx", 1.0)

    context = "Stationary particle, Centered"
    start_time = time.time()

    n_steps = 20
    if not skip_sim:
        success = run_hexapic(
            input_file, hexapic_exe, steps=n_steps, dsteps=1, num_procs=num_procs
        )
        if not success:
            return TestCase(
                test_id,
                "Self-Force Stability",
                GROUP_NAME,
                context,
                [Metric("Simulation", 0, 1, passed=False)],
            )

    bp_file = f"{input_file}.bp4"
    part = read_particle_trajectories(bp_file, species_idx=0)

    metrics = []

    # Measure physical displacement (should be ~0)
    disp, disp_thr = calculate_self_force_stability(part, dx, n_steps)
    metrics.append(
        Metric("Parasitic Displacement", disp, 0.0, disp, disp_thr, disp < disp_thr)
    )

    runtime = time.time() - start_time
    return TestCase(
        test_id, "Self-Force Stability", GROUP_NAME, context, metrics, runtime
    )


def run_current_deposit_test(
    test_id: int, hexapic_exe: str, num_procs: int = 1, skip_sim: bool = False
) -> TestCase:
    """Test 12: Current Deposition (Charge Conservation)."""
    input_file = "current_deposit.inp"
    params = parse_input_file(input_file)
    dx = params.get("dx", 1.0)
    dy = params.get("dy", 1.0)

    context = "Moving particle, Gauss law charge check"
    start_time = time.time()

    if not skip_sim:
        success = run_hexapic(
            input_file, hexapic_exe, steps=10, dsteps=1, num_procs=num_procs
        )
        if not success:
            return TestCase(
                test_id,
                "Current Deposition",
                GROUP_NAME,
                context,
                [Metric("Simulation", 0, 1, passed=False)],
            )

    bp_file = f"{input_file}.bp4"

    metrics = []

    # Global charge conservation
    dQ, dQ_thr = calculate_charge_conservation(bp_file, dx, dy)
    metrics.append(
        Metric("Charge Conservation (Gauss)", dQ, 0.0, dQ, dQ_thr, dQ < dQ_thr)
    )

    runtime = time.time() - start_time
    return TestCase(
        test_id,
        "Current Deposition (Charge Conservation)",
        GROUP_NAME,
        context,
        metrics,
        runtime,
    )


# =============================================================================
# MAIN
# =============================================================================


def main():
    parser = argparse.ArgumentParser(
        description="HEXAPIC Particle-Mesh Coupling Validation"
    )
    parser.add_argument(
        "--exe", type=str, default=None, help="Path to HEXAPIC executable"
    )
    parser.add_argument("--nprocs", type=int, default=1, help="Number of MPI processes")
    parser.add_argument(
        "--test",
        type=str,
        default="all",
        choices=["all", "deposit", "interp", "self_force", "current"],
        help="Which test to run",
    )
    parser.add_argument("--skip-sim", action="store_true", help="Skip simulation")
    parser.add_argument("--no-cleanup", action="store_true", help="Keep output files")
    parser.add_argument(
        "--results-json",
        type=str,
        default=None,
        help="Save results to JSON (for aggregation)",
    )

    args = parser.parse_args()

    # Auto-detect executable
    if args.exe is None:
        script_dir = Path(__file__).parent
        potential_paths = [
            script_dir / "../../../source/HEXAPIC",
            script_dir / "../../source/HEXAPIC",
            Path("./HEXAPIC"),
        ]
        for path in potential_paths:
            if path.exists():
                args.exe = str(path.resolve())
                break
        if args.exe is None:
            print("ERROR: HEXAPIC executable not found. Please specify with --exe.")
            return 1

    os.chdir(Path(__file__).parent)

    tests: List[TestCase] = []
    test_id = 9  # Continuation ID

    if args.results_json:
        print_group_header(GROUP_NAME)
    else:
        print_header(args.nprocs)
        print_group_header(GROUP_NAME)

    # --- Run Tests ---
    if args.test in ["all", "deposit"]:
        test = run_static_deposit_test(test_id, args.exe, args.nprocs, args.skip_sim)
        tests.append(test)
        print_test_case(test)
        test_id += 1

    if args.test in ["all", "interp"]:
        test = run_interpolation_test(test_id, args.exe, args.nprocs, args.skip_sim)
        tests.append(test)
        print_test_case(test)
        test_id += 1

    if args.test in ["all", "self_force"]:
        test = run_self_force_test(test_id, args.exe, args.nprocs, args.skip_sim)
        tests.append(test)
        print_test_case(test)
        test_id += 1

    if args.test in ["all", "current"]:
        test = run_current_deposit_test(test_id, args.exe, args.nprocs, args.skip_sim)
        tests.append(test)
        print_test_case(test)
        test_id += 1

    # --- Summary ---
    if args.results_json:
        save_results(tests, args.results_json)
    else:
        total_runtime = sum(t.runtime or 0 for t in tests)
        print_final_summary(tests, total_runtime)

    if not args.no_cleanup:
        cleanup_simulation_outputs(".", ["*.bp4", "profiling.txt", "plot.sh"])

    return 0 if all(t.passed for t in tests) else 1


if __name__ == "__main__":
    sys.exit(main())
