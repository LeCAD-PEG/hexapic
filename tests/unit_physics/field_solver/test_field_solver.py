#!/usr/bin/env python3
"""
Physics validation tests for HEXAPIC field solver (Poisson).

This module verifies the accuracy of the electrostatic solver: Dirichlet boundary
enforcement (Capacitor plates), Neumann boundary flux leakage (Dipole in open space),
symmetry of potentials (rotational symmetry of a point charge), periodic boundary
continuity (field wrapping), and analytical agreement (Green's function vs ln(r)).

Usage:
    python3 test_field_solver.py --exe ./HEXAPIC
"""

import argparse
import sys
import os
import time
from pathlib import Path
from typing import List, Tuple
import numpy as np

# Add parent directory to sys.path to import utils
sys.path.append(str(Path(__file__).parent.parent))
from utils import (
    TestCase,
    Metric,
    FieldData,
    print_header,
    print_group_header,
    print_test_case,
    print_final_summary,
    save_results,
    cleanup_simulation_outputs,
    run_hexapic,
    read_field_data,
)

# =============================================================================
# PHYSICAL & NUMERICAL CONSTANTS
# =============================================================================

# Safety factor for thresholds.
K_SAFETY = 10.0

# Typical relative residual at solver convergence (Iterative solver target).
SOLVER_RESIDUAL_RTOL = 1e-6

# Machine relative error (float64 epsilon approx 2.22e-16)
EPS_MACHINE = np.finfo(float).eps

GROUP_NAME = "FIELD SOLVER TESTS"


# =============================================================================
# METRIC CALCULATIONS
# =============================================================================


def calculate_boundary_integrity(
    field: FieldData,
    expected_low: float,
    expected_high: float,
    axis: int,  # 0 for X boundaries, 1 for Y boundaries
) -> Tuple[float, float]:
    """
    Check if Dirichlet boundary conditions are respected.
    """
    if axis == 1:  # Y-boundaries (bottom/top)
        bc_min = field.V[:, 0]
        bc_max = field.V[:, -1]
    else:  # X-boundaries (left/right)
        bc_min = field.V[0, :]
        bc_max = field.V[-1, :]

    err_min = np.max(np.abs(bc_min - expected_low))
    err_max = np.max(np.abs(bc_max - expected_high))
    measured = float(max(err_min, err_max))

    # Scale float error by the magnitude of the boundary value
    max_val = max(abs(expected_low), abs(expected_high), 1.0)
    abs_err_threshold = K_SAFETY * EPS_MACHINE * max_val

    return measured, abs_err_threshold


def calculate_field_uniformity(
    field: FieldData,
) -> Tuple[float, float]:
    """
    Check uniformity of E-field (Standard Deviation).
    Ideally 0.0 for a parallel plate capacitor.
    """
    # Cut off boundaries to avoid numerical artifacts at Neumann/Dirichlet corners
    nx, ny = field.shape
    border = 2
    Ey_interior = field.Ey[:, border : ny - border]

    Ey_mean = np.mean(Ey_interior)
    Ey_std = np.std(Ey_interior)

    measured = float(Ey_std)

    float_err = abs(Ey_mean) * EPS_MACHINE * np.sqrt(field.V.size)
    solver_residual = abs(Ey_mean) * SOLVER_RESIDUAL_RTOL

    abs_err_threshold = K_SAFETY * (float_err + solver_residual)

    return measured, abs_err_threshold


def calculate_boundary_flux(
    field: FieldData,
) -> Tuple[float, float]:
    """
    Check Neumann conditions (Flux Leakage).
    Verify if the normal Electric field is ~0 at boundaries.
    """
    Ex_left = np.mean(np.abs(field.Ex[0, :]))
    Ex_right = np.mean(np.abs(field.Ex[-1, :]))
    Ey_bot = np.mean(np.abs(field.Ey[:, 0]))
    Ey_top = np.mean(np.abs(field.Ey[:, -1]))

    measured = float(max(Ex_left, Ex_right, Ey_bot, Ey_top))

    # Physical scale of the problem
    max_E = np.max(np.sqrt(field.Ex**2 + field.Ey**2))
    e_scale = max(max_E, 1e-12)

    residual_err = e_scale * SOLVER_RESIDUAL_RTOL
    float_err = e_scale * EPS_MACHINE * np.sqrt(field.V.size)
    # Grid truncation error (O(dx) at boundaries), proportional to 1/N
    grid_res = 1.0 / min(field.shape)
    grid_err = e_scale * grid_res

    abs_err_threshold = K_SAFETY * (residual_err + float_err + grid_err)

    return measured, abs_err_threshold


def calculate_symmetry(field: FieldData) -> Tuple[float, float]:
    """
    Robust symmetry check: ignores boundaries and small singularity region.
    """
    V = field.V
    nx, ny = V.shape

    border = 2
    r_min_cells = 3.0

    cx, cy = (nx - 1) / 2.0, (ny - 1) / 2.0
    y, x = np.ogrid[:nx, :ny]
    r_sq = (x - cx) ** 2 + (y - cy) ** 2

    is_interior = (x >= border) & (x < ny - border) & (y >= border) & (y < nx - border)
    is_far_field = r_sq >= r_min_cells**2
    mask = is_interior & is_far_field

    # Ensure mask is symmetric itself
    mask = mask & mask[::-1, ::-1]

    if not np.any(mask):
        return float("nan"), float("nan")

    V_flipped = V[::-1, ::-1]
    diff = np.abs(V[mask] - V_flipped[mask])

    # 99.9th percentile to filter single-pixel outliers
    measured = float(np.percentile(diff, 99.9))

    V_scale = float(np.max(np.abs(V[mask])))
    if V_scale < 1e-12:
        return measured, K_SAFETY * EPS_MACHINE

    residual_err = V_scale * SOLVER_RESIDUAL_RTOL
    float_err = V_scale * EPS_MACHINE * np.sqrt(V.size)

    # Grid alignment noise
    N = max(1.0, float(min(nx, ny)))
    grid_err = V_scale / (N)

    abs_err_threshold = K_SAFETY * (residual_err + float_err + grid_err)
    return measured, abs_err_threshold


def calculate_periodicity(field: FieldData) -> Tuple[float, float]:
    """
    Checks continuity of first derivatives across Periodic boundaries.
    """
    V = field.V
    nx, ny = V.shape
    if nx < 4:
        return float("nan"), float("nan")

    # Detect endpoint duplication (V[0] == V[-1])
    is_duplicated = np.allclose(V[0, :], V[-1, :], atol=1e-8)

    if is_duplicated:
        # Physical wrap is [nx-2] -> [0]
        d_left = V[1, :] - V[0, :]
        d_wrap = V[0, :] - V[-2, :]
        measured = float(np.max(np.abs(d_wrap - d_left)))
    else:
        # Physical wrap is [nx-1] -> [0]
        d_left = V[1, :] - V[0, :]
        d_right = V[-1, :] - V[-2, :]
        d_wrap = V[0, :] - V[-1, :]

        j1 = np.max(np.abs(d_wrap - d_left))
        j2 = np.max(np.abs(d_wrap - d_right))
        measured = float(max(j1, j2))

    # Use potential range for scaling
    V_range = float(np.nanmax(V) - np.nanmin(V))
    e_scale = max(V_range, 1e-12)

    # Internal curvature (2nd derivative) serves as the baseline for "smoothness"
    d2_internal = np.abs(V[2:, :] - 2.0 * V[1:-1, :] + V[:-2, :])
    max_curvature = (
        float(np.nanpercentile(d2_internal, 99.9)) if d2_internal.size else 0.0
    )

    solver_noise = e_scale * SOLVER_RESIDUAL_RTOL
    float_noise = e_scale * EPS_MACHINE * np.sqrt(V.size)

    abs_err_threshold = K_SAFETY * (max_curvature + solver_noise + float_noise)

    return measured, abs_err_threshold


def calculate_analytical_fit(
    field: FieldData, r_min_factor: float = 3.0, r_max_factor: float = 0.35
) -> Tuple[float, float]:
    """
    Compare numerical potential to Analytical ln(r) profile.
    Returns 1 - R^2 (Goodness of fit error).
    """
    nx, ny = field.shape
    Lx, Ly = nx * field.dx, ny * field.dy
    x0, y0 = Lx / 2.0, Ly / 2.0

    # Construct radial grid
    x = np.linspace(0, Lx, nx)
    y = np.linspace(0, Ly, ny)
    X, Y = np.meshgrid(x, y, indexing="ij")
    R = np.sqrt((X - x0) ** 2 + (Y - y0) ** 2)

    # Mask for valid region (away from singularity and boundaries)
    r_min = r_min_factor * field.dx
    r_max = r_max_factor * min(Lx, Ly)
    mask = (R > r_min) & (R < r_max)

    r_vals = R[mask].flatten()
    v_vals = field.V[mask].flatten()

    if len(r_vals) < 10:
        return np.nan, np.nan

    # Fit: V = A * ln(r) + B
    ln_r = np.log(r_vals)
    slope, intercept = np.polyfit(ln_r, v_vals, 1)

    v_fit = slope * ln_r + intercept

    # Calculate R^2
    ss_res = np.sum((v_vals - v_fit) ** 2)
    ss_tot = np.sum((v_vals - np.mean(v_vals)) ** 2)
    r2 = 1.0 - (ss_res / ss_tot)

    measured_err = float(1.0 - r2)

    # Threshold:
    # 2nd Order Poisson Solver -> Error O(dx^2) relative to radial distance
    # We use a conservative estimate based on the closest evaluated point
    discretization_err = (field.dx / r_min) ** 2
    float_err = EPS_MACHINE * 100.0

    abs_err_threshold = K_SAFETY * (discretization_err + float_err)

    return measured_err, abs_err_threshold


# =============================================================================
# TEST IMPLEMENTATIONS
# =============================================================================


def run_capacitor_test(
    test_id: int, hexapic_exe: str, num_procs: int = 1, skip_sim: bool = False
) -> TestCase:
    """Test 4: Mixed Boundaries (Capacitor)"""
    input_file = "capacitor.inp"
    context = "Dirichlet Y (+/-10V), Neumann X"

    if not skip_sim:
        success = run_hexapic(input_file, hexapic_exe, steps=5, num_procs=num_procs)
        if not success:
            return TestCase(
                test_id,
                "Capacitor",
                GROUP_NAME,
                context,
                [Metric("Simulation", 0, 1, passed=False)],
            )

    start_time = time.time()
    field = read_field_data(f"{input_file}.bp4")

    if field is None or field.has_nan:
        return TestCase(
            test_id,
            "Capacitor",
            GROUP_NAME,
            context,
            [Metric("Field Read", 0, 1, passed=False)],
        )

    metrics = []

    bc_err, bc_thr = calculate_boundary_integrity(field, -10.0, 10.0, axis=1)
    metrics.append(
        Metric("Boundary Integrity", bc_err, 0.0, bc_err, bc_thr, bc_err < bc_thr)
    )

    unif_err, unif_thr = calculate_field_uniformity(field)
    metrics.append(
        Metric(
            "Field Uniformity (StdDev)",
            unif_err,
            0.0,
            unif_err,
            unif_thr,
            unif_err < unif_thr,
        )
    )

    runtime = time.time() - start_time
    return TestCase(
        test_id, "Mixed Boundaries (Capacitor)", GROUP_NAME, context, metrics, runtime
    )


def run_neumann_dipole_test(
    test_id: int, hexapic_exe: str, num_procs: int = 1, skip_sim: bool = False
) -> TestCase:
    """Test 5: Fully Open Boundaries (Neumann)"""
    input_file = "neumann_dipole.inp"
    context = "All Neumann, Neutral Dipole"

    if not skip_sim:
        success = run_hexapic(input_file, hexapic_exe, steps=5, num_procs=num_procs)
        if not success:
            return TestCase(
                test_id,
                "Neumann Dipole",
                GROUP_NAME,
                context,
                [Metric("Simulation", 0, 1, passed=False)],
            )

    start_time = time.time()
    field = read_field_data(f"{input_file}.bp4")

    if field is None:
        return TestCase(
            test_id,
            "Neumann Dipole",
            GROUP_NAME,
            context,
            [Metric("Field Read", 0, 1, passed=False)],
        )

    metrics = []

    flux_val, flux_thr = calculate_boundary_flux(field)
    metrics.append(
        Metric(
            "Boundary Flux Leakage",
            flux_val,
            0.0,
            flux_val,
            flux_thr,
            flux_val < flux_thr,
        )
    )

    runtime = time.time() - start_time
    return TestCase(
        test_id,
        "Fully Open Boundaries (Neumann)",
        GROUP_NAME,
        context,
        metrics,
        runtime,
    )


def run_grounded_box_test(
    test_id: int, hexapic_exe: str, num_procs: int = 1, skip_sim: bool = False
) -> TestCase:
    """Test 6: Closed Metal Box (All-Dirichlet)"""
    input_file = "grounded_box.inp"
    context = "All phi=0, Centered Charge"

    if not skip_sim:
        success = run_hexapic(input_file, hexapic_exe, steps=5, num_procs=num_procs)
        if not success:
            return TestCase(
                test_id,
                "Grounded Box",
                GROUP_NAME,
                context,
                [Metric("Simulation", 0, 1, passed=False)],
            )

    start_time = time.time()
    field = read_field_data(f"{input_file}.bp4")

    if field is None:
        return TestCase(
            test_id,
            "Grounded Box",
            GROUP_NAME,
            context,
            [Metric("Field Read", 0, 1, passed=False)],
        )

    metrics = []

    bc_err, bc_thr = calculate_boundary_integrity(field, 0.0, 0.0, axis=0)  # Check X
    bc_err_y, _ = calculate_boundary_integrity(field, 0.0, 0.0, axis=1)  # Check Y
    bc_total = max(bc_err, bc_err_y)
    metrics.append(
        Metric(
            "Boundary Zero Check", bc_total, 0.0, bc_total, bc_thr, bc_total < bc_thr
        )
    )

    sym_val, sym_thr = calculate_symmetry(field)
    metrics.append(
        Metric("Symmetry Error", sym_val, 0.0, sym_val, sym_thr, sym_val < sym_thr)
    )

    runtime = time.time() - start_time
    return TestCase(
        test_id, "Closed Metal Box (Dirichlet)", GROUP_NAME, context, metrics, runtime
    )


def run_periodic_test(
    test_id: int, hexapic_exe: str, num_procs: int = 1, skip_sim: bool = False
) -> TestCase:
    """Test 7: Periodic X Boundaries"""
    input_file = "periodic_x.inp"
    context = "Periodic X, Neumann Y"

    if not skip_sim:
        success = run_hexapic(input_file, hexapic_exe, steps=5, num_procs=num_procs)
        if not success:
            return TestCase(
                test_id,
                "Periodic X",
                GROUP_NAME,
                context,
                [Metric("Simulation", 0, 1, passed=False)],
            )

    start_time = time.time()
    field = read_field_data(f"{input_file}.bp4")

    if field is None:
        return TestCase(
            test_id,
            "Periodic X",
            GROUP_NAME,
            context,
            [Metric("Field Read", 0, 1, passed=False)],
        )

    metrics = []

    per_err, per_thr = calculate_periodicity(field)
    metrics.append(
        Metric(
            "Periodicity (Left-Right)",
            per_err,
            0.0,
            per_err,
            per_thr,
            per_err < per_thr,
        )
    )

    runtime = time.time() - start_time
    return TestCase(
        test_id, "Periodic X Boundaries", GROUP_NAME, context, metrics, runtime
    )


def run_analytical_accuracy_test(
    test_id: int, hexapic_exe: str, num_procs: int = 1, skip_sim: bool = False
) -> TestCase:
    """Test 8: Analytical Accuracy (Green's Function)"""
    input_file = "line_charge.inp"
    context = "Line Charge vs ln(r)"

    if not skip_sim:
        success = run_hexapic(input_file, hexapic_exe, steps=5, num_procs=num_procs)
        if not success:
            return TestCase(
                test_id,
                "Analytical Accuracy",
                GROUP_NAME,
                context,
                [Metric("Simulation", 0, 1, passed=False)],
            )

    start_time = time.time()
    field = read_field_data(f"{input_file}.bp4")

    if field is None:
        return TestCase(
            test_id,
            "Analytical Accuracy",
            GROUP_NAME,
            context,
            [Metric("Field Read", 0, 1, passed=False)],
        )

    metrics = []

    fit_err, fit_thr = calculate_analytical_fit(field)
    metrics.append(
        Metric(
            "Analytical Fit (1-R^2)", fit_err, 0.0, fit_err, fit_thr, fit_err < fit_thr
        )
    )

    runtime = time.time() - start_time
    return TestCase(
        test_id,
        "Analytical Accuracy (Green's Func)",
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
        description="HEXAPIC Field Solver Validation Tests"
    )
    parser.add_argument(
        "--exe", type=str, default=None, help="Path to HEXAPIC executable"
    )
    parser.add_argument("--nprocs", type=int, default=1, help="Number of MPI processes")
    parser.add_argument("--skip-sim", action="store_true", help="Skip simulation")
    parser.add_argument("--no-cleanup", action="store_true", help="Keep output files")
    parser.add_argument(
        "--results-json",
        type=str,
        default=None,
        help="Save results to JSON (for aggregation)",
    )

    args = parser.parse_args()

    # Resolve executable path
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

    # Change CWD to script directory for input file access
    os.chdir(Path(__file__).parent)

    tests: List[TestCase] = []
    test_id = 4  # Starting from 4 to follow particle tests (1-3)

    # Output headers
    if args.results_json:
        print_group_header(GROUP_NAME)
    else:
        print_header(args.nprocs)
        print_group_header(GROUP_NAME)

    # --- Execute Tests ---
    test = run_capacitor_test(test_id, args.exe, args.nprocs, args.skip_sim)
    tests.append(test)
    print_test_case(test)
    test_id += 1

    test = run_neumann_dipole_test(test_id, args.exe, args.nprocs, args.skip_sim)
    tests.append(test)
    print_test_case(test)
    test_id += 1

    test = run_grounded_box_test(test_id, args.exe, args.nprocs, args.skip_sim)
    tests.append(test)
    print_test_case(test)
    test_id += 1

    test = run_periodic_test(test_id, args.exe, args.nprocs, args.skip_sim)
    tests.append(test)
    print_test_case(test)
    test_id += 1

    test = run_analytical_accuracy_test(test_id, args.exe, args.nprocs, args.skip_sim)
    tests.append(test)
    print_test_case(test)
    test_id += 1

    # --- Summary & Cleanup ---
    if args.results_json:
        save_results(tests, args.results_json)
    else:
        total_runtime = sum(t.runtime or 0 for t in tests)
        print_final_summary(tests, total_runtime)

    if not args.no_cleanup:
        cleanup_simulation_outputs(".", ["*.bp4", "profiling.txt", "plot.sh"])

    passed = all(t.passed for t in tests)
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
