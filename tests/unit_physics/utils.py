#!/usr/bin/env python3
"""
Shared utilities for HEXAPIC physics validation tests.

This module provides:
- Data structures for test results and metrics
- Formatted output printing
- Temporary directory management
- Result serialization for aggregation
"""

import json
import os
import re
import shutil
import socket
import subprocess
import sys
import tempfile
import numpy as np
import openpmd_api as opmd
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Dict, Any, Tuple
from contextlib import contextmanager


# =============================================================================
# DATA STRUCTURES
# =============================================================================


@dataclass
class Metric:
    """A single metric within a test case."""

    name: str
    measured: float
    expected: Optional[float] = None
    error: Optional[float] = None
    threshold: float = 0.0
    passed: bool = True
    unit: str = ""
    is_percentage: bool = False  # Display threshold/error as percentage
    warn_only: bool = False  # WARN instead of FAIL

    @property
    def status_str(self) -> str:
        if self.passed:
            return "[PASS]"
        elif self.warn_only:
            return "[WARN]"
        else:
            return "[FAIL]"


@dataclass
class TestCase:
    """A complete test case with multiple metrics."""

    id: int
    name: str
    group: str  # e.g., "PARTICLE PUSHER TESTS"
    context: str  # Brief description of test setup
    metrics: List[Metric] = field(default_factory=list)
    runtime: Optional[float] = None

    @property
    def passed(self) -> bool:
        # Warnings don't count as failures
        return all(m.passed or m.warn_only for m in self.metrics)

    @property
    def total_metrics(self) -> int:
        return len(self.metrics)

    @property
    def passed_metrics(self) -> int:
        return sum(1 for m in self.metrics if m.passed)


@dataclass
class FieldData:
    """Container for field data (Potential and E-field)."""

    V: np.ndarray  # Scalar potential [V]
    Ex: np.ndarray  # Electric field X [V/m]
    Ey: np.ndarray  # Electric field Y [V/m]
    dx: float
    dy: float
    shape: Tuple[int, int]
    iterations: int
    has_nan: bool


@dataclass
class ParticleData:
    """Container for particle trajectory data."""

    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    vx: np.ndarray
    vy: np.ndarray
    vz: np.ndarray
    time: np.ndarray


# =============================================================================
# I/O & PARSING HELPER FUNCTIONS
# =============================================================================


def parse_input_file(input_file: str) -> Dict[str, Any]:
    """
    Robustly parses HEXAPIC input file.
    Handles comments, whitespace, and extracts grid/load parameters.
    """
    params = {
        "dt": 0.0,
        "B": [0.0, 0.0, 0.0],
        "dx": 1.0,
        "dy": 1.0,
        "Lx": 0.0,
        "Ly": 0.0,
        "J": 0,
        "K": 0,
        "loads": [],
    }

    if not os.path.exists(input_file):
        return params

    with open(input_file, "r") as f:
        content = f.read()

    current_load = {}

    for line in content.split("\n"):
        # Strip comments (# or //) and whitespace
        line = line.split("#")[0].split("//")[0].strip()
        if not line or "=" not in line:
            continue

        parts = line.split("=")
        key = parts[0].strip()
        try:
            val_str = parts[1].strip()
            val = float(val_str)

            if key == "dt":
                params["dt"] = val
            elif key == "B01":
                params["B"][0] = val
            elif key == "B02":
                params["B"][1] = val
            elif key == "B03":
                params["B"][2] = val
            elif key == "x1s":
                params["x1s"] = val
            elif key == "x1f":
                params["x1f"] = val
            elif key == "x2s":
                params["x2s"] = val
            elif key == "x2f":
                params["x2f"] = val
            elif key == "J":
                params["J"] = int(val)
            elif key == "K":
                params["K"] = int(val)

            # Load parameters
            elif key == "v1drift":
                current_load["v1drift"] = val
            elif key == "v2drift":
                current_load["v2drift"] = val
            elif key == "v3drift":
                current_load["v3drift"] = val
                params["loads"].append(current_load.copy())
                current_load = {}  # Reset for next load

        except (ValueError, IndexError):
            continue

    # Derived parameters
    if "x1f" in params and "x1s" in params:
        params["Lx"] = params["x1f"] - params["x1s"]
    if "x2f" in params and "x2s" in params:
        params["Ly"] = params["x2f"] - params["x2s"]

    if params["J"] > 0 and params["Lx"] > 0:
        params["dx"] = params["Lx"] / params["J"]
    if params["K"] > 0 and params["Ly"] > 0:
        params["dy"] = params["Ly"] / params["K"]

    return params


def read_field_data(bp_file: str) -> Optional[FieldData]:
    """Read the final field state from OpenPMD output."""
    if not os.path.exists(bp_file):
        return None

    series = opmd.Series(bp_file, opmd.Access_Type.read_linear)

    V_final = None
    iterations = 0
    has_nan = False

    # Iterate to find the last step
    for it in series.read_iterations():
        iterations += 1
        V_rc = it.meshes["V"][opmd.Mesh_Record_Component.SCALAR]
        V_chunk = V_rc.load_chunk()
        V_rc.series_flush()

        if np.any(np.isnan(V_chunk)):
            has_nan = True

        V_final = V_chunk.copy()
        it.close()

    if V_final is None:
        return None
    # Get grid info from params (fallback if OpenPMD metadata is missing)
    input_file = bp_file.replace(".bp4", "")
    params = parse_input_file(input_file)
    dx, dy = params.get("dx", 1.0), params.get("dy", 1.0)

    # Compute E-field from potential: E = -grad(V)
    # np.gradient uses 2nd order central differences in interior
    dV_dx, dV_dy = np.gradient(V_final, dx, dy)
    Ex, Ey = -dV_dx, -dV_dy
    return FieldData(
        V=V_final,
        Ex=Ex,
        Ey=Ey,
        dx=dx,
        dy=dy,
        shape=V_final.shape,
        iterations=iterations,
        has_nan=has_nan,
    )


def read_particle_trajectories(bp_file: str, species_idx: int = 0) -> ParticleData:
    """Read particle trajectory data from OpenPMD bp4 file."""
    if not os.path.exists(bp_file):
        # Return empty data structure if file doesn't exist to allow graceful failure later
        return ParticleData(
            np.array([]),
            np.array([]),
            np.array([]),
            np.array([]),
            np.array([]),
            np.array([]),
            np.array([]),
        )

    series = opmd.Series(bp_file, opmd.Access_Type.read_linear)

    x_list, y_list, z_list = [], [], []
    vx_list, vy_list, vz_list = [], [], []
    time_list = []

    for it in series.read_iterations():
        time_rc = it.meshes["time"][opmd.Mesh_Record_Component.SCALAR]
        time_data = time_rc.load_chunk()
        time_rc.series_flush()
        time_list.append(float(time_data[0]))

        sp_key = str(species_idx)
        if sp_key in it.particles:
            sp = it.particles[sp_key]
            # Load chunks
            try:
                pos_x = sp["position"]["x"].load_chunk()
                pos_y = sp["position"]["y"].load_chunk()
                pos_z = sp["position"]["z"].load_chunk()
                vel_x = sp["velocity"]["x"].load_chunk()
                vel_y = sp["velocity"]["y"].load_chunk()
                vel_z = sp["velocity"]["z"].load_chunk()
                it.close()

                if len(pos_x) > 0:
                    x_list.append(float(pos_x[0]))
                    y_list.append(float(pos_y[0]))
                    z_list.append(float(pos_z[0]))
                    vx_list.append(float(vel_x[0]))
                    vy_list.append(float(vel_y[0]))
                    vz_list.append(float(vel_z[0]))
                else:
                    # Particle might have left the domain or wasn't saved
                    for lst in [x_list, y_list, z_list, vx_list, vy_list, vz_list]:
                        lst.append(np.nan)
            except Exception:
                it.close()
                for lst in [x_list, y_list, z_list, vx_list, vy_list, vz_list]:
                    lst.append(np.nan)
        else:
            it.close()
            for lst in [x_list, y_list, z_list, vx_list, vy_list, vz_list]:
                lst.append(np.nan)

    return ParticleData(
        x=np.array(x_list),
        y=np.array(y_list),
        z=np.array(z_list),
        vx=np.array(vx_list),
        vy=np.array(vy_list),
        vz=np.array(vz_list),
        time=np.array(time_list),
    )


# =============================================================================
# OUTPUT FORMATTING
# =============================================================================

# Column widths for metrics table (header and rows must match)
_W_METRIC = 28
_W_MEASURED = 14
_W_EXPECTED = 14
_W_ERROR = 12
_W_THRESH = 12
_W_STATUS = 7
_TABLE_WIDTH = (
    3
    + _W_METRIC
    + 2
    + _W_MEASURED
    + 2
    + _W_EXPECTED
    + 2
    + _W_ERROR
    + 2
    + _W_THRESH
    + 2
    + _W_STATUS
)
# All separators (===, ---) and table borders use this width so lines align
WIDTH = 3 + _TABLE_WIDTH

METRIC_COL_WIDTHS = {
    "name": 26,
    "measured": 14,
    "expected": 14,
    "error": 10,
    "threshold": 10,
    "status": 8,
}


def _fmt_value(val: Optional[float], unit: str = "", is_pct: bool = False) -> str:
    """Format a numeric value for display."""
    if val is None:
        return "---"
    if is_pct:
        return f"{val*100:.2f}%"
    if abs(val) < 1e-10 and val != 0:
        return f"{val:.2e}"
    if abs(val) >= 1e4 or (abs(val) < 0.01 and val != 0):
        return f"{val:.4e}"
    if val == 0:
        return "0.00e+00"
    return f"{val:.4e}"


def _center(text: str, width: int = WIDTH) -> str:
    return text.center(width)


def _line(char: str = "=", width: int = WIDTH) -> str:
    return char * width


def _strip_parens(s: str) -> str:
    """Remove all parenthetical content for clean report formatting."""
    return re.sub(r"\s*\([^)]*\)", "", s).strip()


def print_header(nprocs: int = 1):
    """Print the main test suite header."""
    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    host = socket.gethostname()
    commit = _get_git_commit()

    print(_line("="))
    print(_center("HEXAPIC PHYSICS VALIDATION SUITE"))
    info = f"Date: {now} | Host: {host} | Procs: {nprocs} | Commit: {commit}"
    print(info.center(WIDTH) if len(info) <= WIDTH else info)
    print(_line("="))
    print()


def print_group_header(group_name: str):
    """Print a test group header (e.g., PARTICLE PUSHER TESTS)."""
    print(_line("-"))
    print(_center(group_name.upper()))
    print(_line("-"))
    print()


def print_test_case(test: TestCase):
    """Print a complete test case with all its metrics."""
    # Test header (no parentheses for clean formatting)
    print(f"[TEST {test.id:02d}] {_strip_parens(test.name)}")
    print(f"   Context: {_strip_parens(test.context)}")
    print(f"   {'-' * _TABLE_WIDTH}")

    # Header: metric name | measured value | expected value | error (measured−expected) | max allowed error | status
    hdr = (
        f"   {'METRIC':<{_W_METRIC}}  "
        f"{'MEASURED':>{_W_MEASURED}}  "
        f"{'EXPECTED':>{_W_EXPECTED}}  "
        f"{'ERROR':>{_W_ERROR}}  "
        f"{'MAX ERR':>{_W_THRESH}}  "
        f"{'STATUS':>{_W_STATUS}}"
    )
    print(hdr)
    print(f"   {'-' * _TABLE_WIDTH}")

    # Rows: same column widths
    for m in test.metrics:
        measured_str = _fmt_value(m.measured, m.unit)
        expected_str = (
            _fmt_value(m.expected, m.unit) if m.expected is not None else "---"
        )
        error_str = (
            _fmt_value(m.error, is_pct=m.is_percentage)
            if m.error is not None
            else "---"
        )
        thresh_str = _fmt_value(m.threshold, is_pct=m.is_percentage)
        name_display = _strip_parens(m.name)
        if len(name_display) > _W_METRIC:
            name_display = name_display[: _W_METRIC - 3] + "..."

        row = (
            f"   {name_display:<{_W_METRIC}}  "
            f"{measured_str:>{_W_MEASURED}}  "
            f"{expected_str:>{_W_EXPECTED}}  "
            f"{error_str:>{_W_ERROR}}  "
            f"{thresh_str:>{_W_THRESH}}  "
            f"{m.status_str:>{_W_STATUS}}"
        )
        print(row)

    # Result line
    print(f"   {'-' * _TABLE_WIDTH}")
    result_str = "OK" if test.passed else "FAILED"
    time_str = f"  Time: {test.runtime:.2f}s" if test.runtime else ""
    print(f"   Result: {result_str}{time_str}")
    print()


def print_final_summary(tests: List[TestCase], total_runtime: Optional[float] = None):
    """Print the final summary of all tests."""
    print(_line("="))
    print(_center("SUMMARY"))
    print(_line("="))

    # Group tests by their group
    groups: Dict[str, List[TestCase]] = {}
    for t in tests:
        if t.group not in groups:
            groups[t.group] = []
        groups[t.group].append(t)

    # Print each group
    for group_name, group_tests in groups.items():
        # Clean up group name for display
        display_name = group_name.replace(" TESTS", "").replace("TESTS", "").strip()
        print(display_name)
        for t in group_tests:
            status = "[PASS]" if t.passed else "[FAIL]"
            # Create a short name for summary (no parentheses)
            name_clean = _strip_parens(t.name)
            short_name = name_clean[:30] + "..." if len(name_clean) > 33 else name_clean
            print(f"  {t.id}. {short_name:<35} {status}")
        print()

    # Overall stats
    total = len(tests)
    passed = sum(1 for t in tests if t.passed)
    total_metrics = sum(t.total_metrics for t in tests)
    passed_metrics = sum(t.passed_metrics for t in tests)

    print(_line("-"))
    status_str = "PASSED" if passed == total else "FAILED"
    print(
        f"OVERALL: {status_str} ({passed}/{total} tests passed, {passed_metrics}/{total_metrics} metrics)"
    )
    if total_runtime:
        print(f"Total Runtime: {total_runtime:.2f}s")
    print(_line("="))
    print()


# =============================================================================
# SERIALIZATION (for aggregating results across scripts)
# =============================================================================


def _convert_for_json(obj: Any) -> Any:
    """Convert numpy types and other non-serializable objects."""
    if hasattr(obj, "item"):
        return obj.item()
    if isinstance(obj, (list, tuple)):
        return [_convert_for_json(x) for x in obj]
    if isinstance(obj, dict):
        return {k: _convert_for_json(v) for k, v in obj.items()}
    return obj


def save_results(tests: List[TestCase], filepath: str):
    """Save test results to JSON for later aggregation."""
    data = []
    for t in tests:
        t_dict = {
            "id": t.id,
            "name": t.name,
            "group": t.group,
            "context": t.context,
            "runtime": t.runtime,
            "metrics": [_convert_for_json(asdict(m)) for m in t.metrics],
        }
        data.append(_convert_for_json(t_dict))

    with open(filepath, "w") as f:
        json.dump(data, f, indent=2)


def load_results(filepath: str) -> List[TestCase]:
    """Load test results from JSON."""
    if not os.path.exists(filepath):
        return []
    with open(filepath, "r") as f:
        data = json.load(f)

    tests = []
    for t_dict in data:
        metrics = [Metric(**m) for m in t_dict["metrics"]]
        tests.append(
            TestCase(
                id=t_dict["id"],
                name=t_dict["name"],
                group=t_dict["group"],
                context=t_dict["context"],
                metrics=metrics,
                runtime=t_dict.get("runtime"),
            )
        )
    return tests


def aggregate_results_main(argv: List[str]) -> int:
    """
    Aggregate test results from multiple JSON files and print the final summary.
    argv: list of arguments (argv[0] is program name, argv[1:] are JSON file paths).
    Returns: 0 if all tests passed, 1 otherwise or on usage/empty error.
    """
    if len(argv) < 2:
        print("Usage: utils.py result1.json result2.json ...")
        return 1

    all_tests: List[TestCase] = []
    total_runtime = 0.0

    for filepath in argv[1:]:
        tests = load_results(filepath)
        all_tests.extend(tests)
        total_runtime += sum(t.runtime or 0 for t in tests)

    if not all_tests:
        print("No results found to aggregate.")
        return 1

    for i, test in enumerate(all_tests, 1):
        test.id = i

    print_final_summary(all_tests, total_runtime)

    return 0 if all(t.passed for t in all_tests) else 1


# =============================================================================
# TEMPORARY DIRECTORY MANAGEMENT
# =============================================================================


@contextmanager
def temp_test_directory(base_dir: str, prefix: str = "run_", cleanup: bool = True):
    """
    Context manager for creating and optionally cleaning up a temporary test directory.

    Usage:
        with temp_test_directory("/path/to/base", "test_") as tmp_dir:
            # Run tests in tmp_dir
            pass
        # Directory is automatically cleaned up
    """
    tmp_dir = tempfile.mkdtemp(dir=base_dir, prefix=prefix)
    try:
        yield tmp_dir
    finally:
        if cleanup and os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir, ignore_errors=True)


def cleanup_simulation_outputs(directory: str, patterns: List[str] = None):
    """
    Clean up simulation output files.

    Args:
        directory: Directory to clean
        patterns: List of patterns to remove (default: common simulation outputs)
    """
    if patterns is None:
        patterns = ["*.bp4", "profiling.txt", "plot.sh"]

    for pattern in patterns:
        for path in Path(directory).glob(pattern):
            if path.is_dir():
                shutil.rmtree(path, ignore_errors=True)
            else:
                path.unlink(missing_ok=True)


# =============================================================================
# CONSTANTS
# =============================================================================

# Grid noise floor - numerical noise from position-to-grid interpolation
# Typical value: 0.001% of cell size (1e-5)
GRID_NOISE_FLOOR = 1e-5

# Machine precision for float64
EPS_MACHINE = np.finfo(float).eps


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================


def _get_git_commit() -> str:
    """Get current git commit hash (short form)."""
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode == 0:
            return result.stdout.strip()
    except Exception:
        pass
    return "unknown"


def run_hexapic(
    input_file: str,
    hexapic_exe: str,
    steps: int,
    dsteps: int = 1,
    num_procs: int = 1,
    quiet: bool = True,
    extra_args: List[str] = None,
    force_mpi: bool = False,
) -> bool:
    """
    Run HEXAPIC simulation.
    For num_procs==1, tries direct execution first (avoids mpiexec/PMIx in containers or as root).
    If force_mpi=True and num_procs==1, use mpiexec -n 1 so MPI/OpenPMD output is written correctly.

    Returns:
        True if simulation completed successfully
    """
    base = [hexapic_exe, input_file, "-steps", str(steps), "-dsteps", str(dsteps)]
    if quiet:
        base.append("-quiet")
    if extra_args:
        base.extend(extra_args)

    def _run(cmd: List[str]) -> bool:
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            return True
        except subprocess.CalledProcessError as e:
            if e.stderr:
                print(f"Simulation failed: {e.stderr}")
            return False

    if num_procs == 1:
        if force_mpi:
            cmd = ["mpiexec", "-n", "1"] + base
            if _run(cmd):
                return True
            try:
                if os.geteuid() == 0:
                    out = subprocess.run(
                        ["mpiexec", "--help"], capture_output=True, text=True, timeout=5
                    )
                    if "allow-run-as-root" in (out.stdout or "") + (out.stderr or ""):
                        if _run(["mpiexec", "--allow-run-as-root", "-n", "1"] + base):
                            return True
            except (AttributeError, FileNotFoundError, subprocess.TimeoutExpired):
                pass
            return False
        if _run(base):
            return True
        # Fallback to mpiexec (e.g. MPI_Init requires launcher)
        cmd = ["mpiexec", "-n", "1"] + base
        if _run(cmd):
            return True
        # Only add --allow-run-as-root when launcher supports it (Open MPI), not e.g. mpiexec.slurm
        try:
            if os.geteuid() == 0:
                out = subprocess.run(
                    ["mpiexec", "--help"], capture_output=True, text=True, timeout=5
                )
                if "allow-run-as-root" in (out.stdout or "") + (out.stderr or ""):
                    if _run(["mpiexec", "--allow-run-as-root", "-n", "1"] + base):
                        return True
        except (AttributeError, FileNotFoundError, subprocess.TimeoutExpired):
            pass
        return False
    else:
        # Under SLURM use srun so allocated tasks are used; otherwise mpiexec
        if os.environ.get("SLURM_JOB_ID"):
            cmd = ["srun", "-n", str(num_procs)] + base
            if _run(cmd):
                return True
        cmd = ["mpiexec", "-n", str(num_procs)] + base
        if _run(cmd):
            return True
        try:
            if os.geteuid() == 0:
                out = subprocess.run(
                    ["mpiexec", "--help"], capture_output=True, text=True, timeout=5
                )
                if "allow-run-as-root" in (out.stdout or "") + (out.stderr or ""):
                    if _run(
                        ["mpiexec", "--allow-run-as-root", "-n", str(num_procs)] + base
                    ):
                        return True
        except (AttributeError, FileNotFoundError, subprocess.TimeoutExpired):
            pass
        return False


# =============================================================================
# LEGACY COMPATIBILITY (for gradual migration)
# =============================================================================


@dataclass
class TestResult:
    """Legacy TestResult for backwards compatibility."""

    group: str
    name: str
    passed: bool
    metric_value: float
    threshold: float
    message: str


def print_test_result(result: TestResult):
    """Legacy print function for backwards compatibility."""
    status = "PASS" if result.passed else "FAIL"
    print(f"[{status}] {result.name}: {result.message}")


if __name__ == "__main__":
    sys.exit(aggregate_results_main(sys.argv))
