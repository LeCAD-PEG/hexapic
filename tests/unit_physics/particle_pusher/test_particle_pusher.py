#!/usr/bin/env python3
"""
Physics validation tests for HEXAPIC particle pusher.

Tests: Multi-Species Larmor Orbit (In-Plane; gyroradius, energy conservation, drift stability),
Larmor Orbit (Out-of-Plane / 3V; Boris rotator, electron only), ExB Drift (velocity in crossed E and B).
Thresholds account for Boris O(dt^2) phase errors, floating point accumulation over N steps,
and sampling error from discrete output steps.
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
    ParticleData,
    print_header,
    print_group_header,
    print_test_case,
    print_final_summary,
    save_results,
    cleanup_simulation_outputs,
    run_hexapic,
    read_particle_trajectories,
    parse_input_file,
)

# =============================================================================
# PHYSICAL & NUMERICAL CONSTANTS
# =============================================================================

Q_E = 1.602176634e-19  # Elementary charge [C]
M_E = 9.1093837139e-31  # Electron mass [kg]

# Safety factor for thresholds.
K_SAFETY = 10.0

# Machine relative error (float64 epsilon approx 2.22e-16)
EPS_MACHINE = np.finfo(float).eps

GROUP_NAME = "PARTICLE PUSHER TESTS"


# =============================================================================
# METRIC CALCULATIONS
# =============================================================================


def calculate_energy_conservation(
    particle: ParticleData, mass: float
) -> Tuple[float, float]:
    """
    Calculate max kinetic energy deviation.
    Boris pusher is not perfectly symplectic but conserves energy very well.
    Error is dominated by floating point accumulation.
    """
    valid = ~np.isnan(particle.vx) & ~np.isnan(particle.vy) & ~np.isnan(particle.vz)
    if not np.any(valid):
        return np.nan, np.nan

    vx, vy, vz = particle.vx[valid], particle.vy[valid], particle.vz[valid]
    E_k = 0.5 * mass * (vx**2 + vy**2 + vz**2)
    E_k0 = E_k[0]

    measured = float(np.max(np.abs(E_k - E_k0)))

    float_err = EPS_MACHINE * np.sqrt(len(E_k))

    abs_err_threshold = K_SAFETY * float_err

    return measured, abs_err_threshold


def calculate_separation(
    e: ParticleData,
    p: ParticleData,
    dx: float,
    dy: float,
    R_larmor: float,
    omega_c: float,
    dt: float,
) -> Tuple[float, float]:
    """
    Calculate max separation between electron and positron.
    Expected to be 4*R_larmor for opposite initial velocities.
    """
    valid = ~np.isnan(e.x) & ~np.isnan(e.y) & ~np.isnan(p.x) & ~np.isnan(p.y)
    if not np.any(valid):
        return np.nan, np.nan

    x_e, y_e = e.x[valid] * dx, e.y[valid] * dy
    x_p, y_p = p.x[valid] * dx, p.y[valid] * dy

    measured = float(np.max(np.sqrt((x_e - x_p) ** 2 + (y_e - y_p) ** 2)))
    N_steps = len(x_e)

    phase_err_rel = (omega_c * dt) ** 2 / 12.0
    phase_err_rel *= 20.0  # safety factor
    boris_err = 4.0 * R_larmor * phase_err_rel
    float_err = 4.0 * R_larmor * EPS_MACHINE * np.sqrt(N_steps)

    abs_err_threshold = K_SAFETY * (boris_err + float_err)

    return measured, abs_err_threshold


def calculate_extent(
    arr: np.ndarray,
    dx: float,
    R_larmor: float,
    omega_c: float,
    dt: float,
) -> Tuple[float, float]:
    """
    Calculate trajectory extent (max - min).
    Expected to be 2*R_larmor.
    """
    valid = arr[~np.isnan(arr)]
    if len(valid) == 0:
        return np.nan, np.nan

    measured = float((np.max(valid) - np.min(valid)) * dx)
    N_steps = len(valid)

    phase_err_rel = (omega_c * dt) ** 2 / 12.0
    phase_err_rel *= 10.0  # safety factor
    boris_err = 2.0 * R_larmor * phase_err_rel
    float_err = 2.0 * R_larmor * EPS_MACHINE * np.sqrt(N_steps)

    abs_err_threshold = K_SAFETY * (boris_err + float_err)

    return measured, abs_err_threshold


def calculate_stability(
    arr: np.ndarray,
    scale: float,
    ref_scale: float,
    dx: float,
) -> Tuple[float, float]:
    """
    Calculate max deviation from initial value (Drift).
    Should be 0 for axes perpendicular to B-field force.
    """
    valid = arr[~np.isnan(arr)]
    if len(valid) == 0:
        return np.nan, np.nan

    measured = float(np.max(np.abs(valid - valid[0])) * scale)
    N_steps = len(valid)

    float_err = ref_scale * EPS_MACHINE * np.sqrt(N_steps)

    abs_err_threshold = K_SAFETY * float_err

    return measured, abs_err_threshold


def calculate_hamiltonian_conservation(
    particle: ParticleData,
    mass: float,
    charge: float,
    Ex: float,
    dx: float,
    dt: float,
    B: float,
) -> Tuple[float, float]:
    """
    Calculate max Hamiltonian deviation.
    H = (1/2)mv^2 + q*phi, where phi = -Ex * x.
    """
    valid = (
        ~np.isnan(particle.vx)
        & ~np.isnan(particle.vy)
        & ~np.isnan(particle.vz)
        & ~np.isnan(particle.x)
    )
    if not np.any(valid):
        return np.nan, np.nan

    vx = particle.vx[valid]
    vy = particle.vy[valid]
    vz = particle.vz[valid]
    x = particle.x[valid] * dx

    E_k = 0.5 * mass * (vx**2 + vy**2 + vz**2)
    phi = -Ex * x
    H = E_k + charge * phi

    H0 = H[0]
    ref_H = abs(H0) if abs(H0) > 1e-30 else 1.0
    measured = float(np.max(np.abs(H - H0)))

    omega_c = np.abs(charge) * B / mass
    boris_err = ref_H * (omega_c * dt) ** 2 / 12.0
    float_err = ref_H * EPS_MACHINE * np.sqrt(len(H))

    abs_err_threshold = K_SAFETY * (boris_err + float_err)

    return measured, abs_err_threshold


def calculate_drift_velocity(
    particle: ParticleData,
    v_expected: float,
    mass: float,
    charge: float,
    B: float,
    dt: float,
) -> Tuple[float, float]:
    """
    Calculate drift velocity (mean of vy over second half of sim).
    """
    valid = ~np.isnan(particle.vy) & ~np.isnan(particle.time)
    if not np.any(valid):
        return np.nan, np.nan

    vy = particle.vy[valid]
    half_idx = len(vy) // 2
    if half_idx < 1:
        return np.nan, np.nan

    measured = float(np.mean(vy[half_idx:]))

    omega_c = np.abs(charge) * B / mass
    boris_err = np.abs(v_expected) * ((omega_c * dt) ** 2 / 8.0)
    float_err = EPS_MACHINE * np.abs(v_expected) * np.sqrt(len(vy))

    abs_err_threshold = K_SAFETY * (boris_err + float_err)

    return measured, abs_err_threshold


# =============================================================================
# TEST IMPLEMENTATIONS
# =============================================================================


def run_in_plane_larmor_test(
    test_id: int, hexapic_exe: str, num_procs: int = 1, skip_sim: bool = False
) -> TestCase:
    """Test 1: In-Plane Larmor Orbit"""
    input_file = "in_plane_larmor.inp"
    params = parse_input_file(input_file)

    dt = params.get("dt", 1e-12)
    B0 = params["B"][2]
    dx, dy = params.get("dx", 1.0), params.get("dy", 1.0)

    # Safety check if parsing failed
    v_perp = 1e6
    if params["loads"] and "v1drift" in params["loads"][0]:
        v_perp = params["loads"][0]["v1drift"]

    # Theoretical Larmor parameters
    omega = np.abs(Q_E) * B0 / M_E
    R_larmor = M_E * v_perp / (np.abs(Q_E) * B0)
    T_larmor = 2 * np.pi / omega
    T_total = 3 * T_larmor
    steps = int(T_total / dt)

    context = f"B={B0}T, v={v_perp:.0e} m/s, dt={dt*1e12:.0f}ps, Steps={steps}"

    start_time = time.time()

    if not skip_sim:
        success = run_hexapic(
            input_file, hexapic_exe, steps, dsteps=1, num_procs=num_procs
        )
        if not success:
            return TestCase(
                test_id,
                "Larmor In-Plane",
                GROUP_NAME,
                context,
                [Metric("Simulation", 0, 1, passed=False)],
            )

    bp_file = f"{input_file}.bp4"
    electron = read_particle_trajectories(bp_file, species_idx=0)
    positron = read_particle_trajectories(bp_file, species_idx=1)

    # Skip transient (first 2/3)
    n = len(electron.time)
    start = (2 * n) // 3
    if start < n:
        electron = ParticleData(
            x=electron.x[start:],
            y=electron.y[start:],
            z=electron.z[start:],
            vx=electron.vx[start:],
            vy=electron.vy[start:],
            vz=electron.vz[start:],
            time=electron.time[start:],
        )
        positron = ParticleData(
            x=positron.x[start:],
            y=positron.y[start:],
            z=positron.z[start:],
            vx=positron.vx[start:],
            vy=positron.vy[start:],
            vz=positron.vz[start:],
            time=positron.time[start:],
        )

    runtime = time.time() - start_time
    metrics = []

    energy_abs_err, err_threshold = calculate_energy_conservation(electron, M_E)
    metrics.append(
        Metric(
            "Energy Conservation",
            energy_abs_err,
            0.0,
            energy_abs_err,
            err_threshold,
            energy_abs_err < err_threshold,
        )
    )

    expected_sep = 4 * R_larmor
    measured_sep, sep_threshold = calculate_separation(
        electron, positron, dx, dy, R_larmor, omega, dt
    )
    sep_error = np.abs(measured_sep - expected_sep)
    metrics.append(
        Metric(
            "4R Separation",
            measured_sep,
            expected_sep,
            sep_error,
            sep_threshold,
            sep_error < sep_threshold,
        )
    )

    expected_y = 2 * R_larmor
    measured_y, y_threshold = calculate_extent(electron.y, dy, R_larmor, omega, dt)
    y_error = np.abs(measured_y - expected_y)
    metrics.append(
        Metric(
            "Y-projection trajectory",
            measured_y,
            expected_y,
            y_error,
            y_threshold,
            y_error < y_threshold,
        )
    )

    z_drift, z_threshold = calculate_stability(electron.z, 1.0, R_larmor, dx)
    metrics.append(
        Metric(
            "Z-Axis Drift", z_drift, 0.0, z_drift, z_threshold, z_drift <= z_threshold
        )
    )

    return TestCase(
        test_id,
        "In-Plane Larmor Orbit (In-Plane)",
        GROUP_NAME,
        context,
        metrics,
        runtime,
    )


def run_3v_larmor_test(
    test_id: int, hexapic_exe: str, num_procs: int = 1, skip_sim: bool = False
) -> TestCase:
    """Test 2: 3V / Out-of-Plane Larmor Orbit"""
    input_file = "3v_larmor.inp"
    params = parse_input_file(input_file)

    dt = params.get("dt", 1e-12)
    B0 = params["B"][1]  # By != 0
    dx = params.get("dx", 1.0)

    v_perp = 1e6
    if params["loads"] and "v3drift" in params["loads"][0]:
        v_perp = params["loads"][0]["v3drift"]

    omega = np.abs(Q_E) * B0 / M_E
    R_larmor = M_E * v_perp / (np.abs(Q_E) * B0)
    T_larmor = 2 * np.pi / omega
    T_total = 3 * T_larmor
    steps = int(T_total / dt)
    dsteps = max(1, steps // 200)

    context = f"B_y={B0}T, v_z={v_perp:.0e} m/s, 3D Rotation Check"

    start_time = time.time()

    if not skip_sim:
        success = run_hexapic(
            input_file, hexapic_exe, steps, dsteps=dsteps, num_procs=num_procs
        )
        if not success:
            return TestCase(
                test_id,
                "Larmor 3V",
                GROUP_NAME,
                context,
                [Metric("Simulation", 0, 1, passed=False)],
            )

    bp_file = f"{input_file}.bp4"
    electron = read_particle_trajectories(bp_file, species_idx=0)

    n = len(electron.time)
    start = (2 * n) // 3
    if start < n:
        electron = ParticleData(
            x=electron.x[start:],
            y=electron.y[start:],
            z=electron.z[start:],
            vx=electron.vx[start:],
            vy=electron.vy[start:],
            vz=electron.vz[start:],
            time=electron.time[start:],
        )

    runtime = time.time() - start_time
    metrics = []

    energy_abs_err, err_threshold = calculate_energy_conservation(electron, M_E)
    metrics.append(
        Metric(
            "Energy Conservation",
            energy_abs_err,
            0.0,
            energy_abs_err,
            err_threshold,
            energy_abs_err < err_threshold,
        )
    )

    expected_x = 2 * R_larmor
    measured_x, x_threshold = calculate_extent(electron.x, dx, R_larmor, omega, dt)
    x_error = np.abs(measured_x - expected_x)
    metrics.append(
        Metric(
            "X-Projection Trajectory",
            measured_x,
            expected_x,
            x_error,
            x_threshold,
            x_error < x_threshold,
        )
    )

    y_drift, y_threshold = calculate_stability(electron.y, 1.0, R_larmor, dx)
    metrics.append(
        Metric(
            "Y-Axis Drift", y_drift, 0.0, y_drift, y_threshold, y_drift <= y_threshold
        )
    )

    return TestCase(test_id, "3V Larmor Orbit", GROUP_NAME, context, metrics, runtime)


def run_exb_drift_test(
    test_id: int, hexapic_exe: str, num_procs: int = 1, skip_sim: bool = False
) -> TestCase:
    """Test 3: ExB Drift"""
    input_file = "exb_drift.inp"
    params = parse_input_file(input_file)

    dt = params.get("dt", 1e-12)
    Bz = params["B"][2]
    dx = params.get("dx", 1.0)

    # E-field from boundaries: V=0 at x=0, V=100 at x=0.01
    # Assuming standard test case setup
    Ex = -100.0 / 0.01
    v_exb_expected = -Ex / Bz

    steps = 400
    dsteps = 20
    context = f"E={Ex:.0e} V/m, B={Bz}T, Expected v_drift={v_exb_expected:.0e} m/s"

    start_time = time.time()

    if not skip_sim:
        success = run_hexapic(
            input_file, hexapic_exe, steps, dsteps=dsteps, num_procs=num_procs
        )
        if not success:
            return TestCase(
                test_id,
                "ExB Drift",
                GROUP_NAME,
                context,
                [Metric("Simulation", 0, 1, passed=False)],
            )

    bp_file = f"{input_file}.bp4"
    electron = read_particle_trajectories(bp_file, species_idx=0)

    runtime = time.time() - start_time
    metrics = []

    H_dev, err_threshold = calculate_hamiltonian_conservation(
        electron, M_E, -Q_E, Ex, dx, dt, Bz
    )
    metrics.append(
        Metric(
            "Hamiltonian Conservation",
            H_dev,
            0.0,
            H_dev,
            err_threshold,
            H_dev < err_threshold,
        )
    )

    measured_vy, vy_threshold = calculate_drift_velocity(
        electron, v_exb_expected, M_E, -Q_E, Bz, dt
    )
    vy_error = np.abs(measured_vy - v_exb_expected)
    metrics.append(
        Metric(
            "Drift Velocity (Vy)",
            measured_vy,
            v_exb_expected,
            vy_error,
            vy_threshold,
            vy_error < vy_threshold,
        )
    )

    return TestCase(test_id, "ExB Drift", GROUP_NAME, context, metrics, runtime)


# =============================================================================
# MAIN
# =============================================================================


def main():
    parser = argparse.ArgumentParser(
        description="HEXAPIC Particle Pusher Validation Tests"
    )
    parser.add_argument(
        "--exe", type=str, default=None, help="Path to HEXAPIC executable"
    )
    parser.add_argument("--nprocs", type=int, default=1, help="Number of MPI processes")
    parser.add_argument(
        "--test",
        type=str,
        default="all",
        choices=["all", "in_plane", "3v", "exb"],
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
    test_id = 1

    # Output headers
    if args.results_json:
        print_group_header(GROUP_NAME)
    else:
        print_header(args.nprocs)
        print_group_header(GROUP_NAME)

    # --- Execute Tests ---
    if args.test in ["all", "in_plane"]:
        test = run_in_plane_larmor_test(test_id, args.exe, args.nprocs, args.skip_sim)
        tests.append(test)
        print_test_case(test)
        test_id += 1

    if args.test in ["all", "3v"]:
        test = run_3v_larmor_test(test_id, args.exe, args.nprocs, args.skip_sim)
        tests.append(test)
        print_test_case(test)
        test_id += 1

    if args.test in ["all", "exb"]:
        test = run_exb_drift_test(test_id, args.exe, args.nprocs, args.skip_sim)
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
