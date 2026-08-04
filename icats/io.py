#!/usr/bin/env python3
"""Input/output helpers for icats.

Sections:
- Molecule-level I/O validation
- Run-level I/O validation and template generation
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple

from .functions import File2InputList


# ----------------------------- Molecule-Level I/O -----------------------------
def validate_molecule_file_refs(base_dir: Path, mol_file: str) -> List[str]:
    """Validate referenced files inside a molecule input file.

    Returns warning strings for non-fatal issues.
    """
    warnings: List[str] = []
    mol_path = (base_dir / mol_file).resolve()
    if not mol_path.is_file():
        raise ValueError(f"Referenced molecule file does not exist: {mol_file}")

    entries = File2InputList(str(mol_path))
    for key, vals in entries:
        if key in ("xyz", "hess", "w"):
            if not vals:
                raise ValueError(f"Missing value for '{key}' in {mol_file}")
            fpath = (mol_path.parent / vals[0]).resolve()
            if key == "hess" and not fpath.exists():
                warnings.append(f"Optional hessian file missing for {mol_file}: {vals[0]}")
            elif key in ("xyz", "w") and not fpath.exists():
                raise ValueError(f"Referenced file missing in {mol_file}: {vals[0]}")
    return warnings


# ------------------------------- Run-Level I/O --------------------------------
KNOWN_KEYS = {
    "mol",
    "tvel",
    "relative-velocity",
    "relative-velocity-fwhm",
    "collision-energy",
    "incoming-p0",
    "incoming-k",
    "fileout",
    "dirout",
    "seed",
    "trot",
    "tvib",
    "maxb",
    "fixed-b",
    "impact-phi",
    "output-frame",
    "orbital-sampling",
    "vib-mode",
    "maxj",
    "maxl",
    "chi",
    "maxv",
    "nsamp",
    "workers",
    "rz",
    "beam-angle",
    "ordist",
    "orientation-mode",
    "orientation-frame",
    "orientation-thin",
    "printout",
    "rot-param",
    "phisample",
    "plothist",
    "hist_initial",
    "hist_sampled",
    "wang",
    "keepinfo",
    "plotinit",
    "continue",
    "wlmode",
    "wl-target",
    "wl-ff",
    "wl-nstep",
    "wl-flatness",
    "wl-wn-factor",
    "wl-wn",
    "wl-j-bins",
    "wl-j-range",
    "wl-j-min",
    "wl-low-j-scale",
    "wl-l-cap",
    "wl-angular-sampler",
    "wl-audit-angular-sampler",
    "audit-initial-sample",
    "audit-initial-energy-tol",
    "audit-initial-angular-tol",
    "audit-initial-vib-tol",
    "audit-initial-velocity-tol",
    "seed-mode",
    "run-mode",
    "run-tag",
    "wl-tol",
    "wl-max-iter",
    "wl-log-every",
    "progress",
    "dry-run",
    "check-input",
    "save-frequency",
    "output-format",
    "units-out",
}


def _parse_bool(raw: str, key: str) -> None:
    if raw not in ("True", "False"):
        raise ValueError(f"Invalid boolean for '{key}': {raw} (use True/False)")


def _ensure_numeric(raw: str, key: str, as_int: bool = False) -> None:
    try:
        _ = int(raw) if as_int else float(raw)
    except ValueError as exc:
        kind = "integer" if as_int else "number"
        raise ValueError(f"Invalid {kind} for '{key}': {raw}") from exc


def load_and_validate(input_file: str) -> Dict[str, object]:
    """Parse and validate top-level input, returning a resolved summary dict."""
    input_path = Path(input_file).resolve()
    if not input_path.is_file():
        raise ValueError(f"Input file does not exist: {input_file}")

    entries: List[Tuple[str, List[str]]] = File2InputList(str(input_path))
    warnings: List[str] = []
    mol_entries: List[Tuple[str, List[str]]] = []

    for key, vals in entries:
        if key not in KNOWN_KEYS:
            warnings.append(f"Unknown key (ignored by current parser): {key}")
            continue
        if key == "mol":
            if len(vals) < 2:
                raise ValueError("Each 'mol' entry must be: mol = <0|1> <molecule_file>")
            if vals[0] not in ("0", "1"):
                raise ValueError(f"Invalid molecule index in 'mol': {vals[0]} (expected 0 or 1)")
            mol_entries.append((key, vals))
            continue
        if key in ("wang", "keepinfo", "continue", "plothist", "phisample", "dry-run", "check-input", "hist_initial", "hist_sampled", "wl-audit-angular-sampler", "audit-initial-sample"):
            _parse_bool(vals[0], key)
        if key in ("workers", "seed", "maxj", "maxl", "maxv", "wl-nstep", "wl-wn", "wl-j-bins", "wl-max-iter", "wl-log-every", "save-frequency"):
            _ensure_numeric(vals[0], key, as_int=True)
        if key in ("tvel", "relative-velocity", "relative-velocity-fwhm", "collision-energy", "incoming-p0", "incoming-k", "trot", "tvib", "maxb", "fixed-b", "impact-phi", "rz", "beam-angle", "wl-ff", "wl-flatness", "wl-wn-factor", "wl-j-range", "wl-l-cap", "wl-tol", "audit-initial-energy-tol", "audit-initial-angular-tol", "audit-initial-vib-tol", "audit-initial-velocity-tol"):
            _ensure_numeric(vals[0], key, as_int=False)

    if len(mol_entries) < 2:
        raise ValueError("Expected two 'mol' entries (for molecule 0 and 1).")

    by_idx = {vals[0] for _, vals in mol_entries}
    if by_idx != {"0", "1"}:
        raise ValueError("Input must include both 'mol = 0 ...' and 'mol = 1 ...' entries.")

    for _, vals in mol_entries:
        warnings.extend(validate_molecule_file_refs(input_path.parent, vals[1]))

    return {
        "input_path": str(input_path),
        "entries": entries,
        "warnings": warnings,
    }


def write_run_log(path: str, lines: List[str]) -> None:
    Path(path).write_text("".join(lines))


def write_templates(base_dir: str = ".", overwrite: bool = False) -> List[str]:
    """Write starter template files for a standard NH3 + H2O workflow.

    Returns the list of written file paths.
    """
    out = Path(base_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)

    input_txt = """# icats input template (defaults-oriented)
# Molecule inputs (required)
mol = 0 ammonia_dat.txt
mol = 1 h2o_dat.txt

# Core sampling controls
Nsamp = 10000
# For Wang-Landau runs use maxj; for non-WL runs use maxl.
maxl = 20
# maxj = 20
maxb = 10

# Thermodynamic conditions
Tvib = 500.0
Trot = 500.0
# Tvel = 500.0
# Tvel = -500.0 50.0
# relative-velocity = 1000.0
# collision-energy = 0.050
# incoming-k = 12.0

# Geometry / beam
Rz = 15
# beam-angle = 90.0

# Output
fileout = out
dirout = outputs
printout = 0 0 0 0
plothist = False
hist_initial = False
hist_sampled = False

# Runtime
workers = 1
seed = 400
continue = False

# Wang-Landau (disabled by default)
wang = False
wlmode = default
# wl-target is automatic by default:
#   geometric orbital-sampling -> linear-j
#   flat-l orbital-sampling    -> flat-j
# wl-target = linear-j
# wl-target = flat-j
# wlmode = fast
# wlmode = accurate
# wl-ff = 1.05
# wl-nstep = 500
# wl-flatness = 0.90
# wl-wn-factor = 4.0
# wl-wn = 80
# wl-j-bins = 80
# wl-j-range = 60
# wl-j-min = 0.0
# wl-low-j-scale = 0.25
# wl-l-cap = 60
# wl-tol = 1.000001
# wl-max-iter = 0
# wl-log-every = 1
# wl-angular-sampler = fast
# wl-audit-angular-sampler = False
# audit-initial-sample = False
# audit-initial-energy-tol = 0.02
# audit-initial-angular-tol = 0.0
# audit-initial-vib-tol = 0.0
# audit-initial-velocity-tol = 0.0

# User-friendly run behavior
run-mode = fresh
# run-tag = my_run_name
seed-mode = fixed
progress = normal
dry-run = False
check-input = False
save-frequency = 0
output-format = xyzvel
units-out = ang-fs
output-frame = internal
"""

    h2o_txt = """# H2O molecule input template
# Required file links (update paths as needed)
xyz  = h2o_geom.xyz
# hess is the full mass-weighted Cartesian Hessian in atomic units.
# Its eigenvalues should be omega^2; ICATS diagonalizes this matrix directly.
hess = h2o_hessian.txt
w    = h2o_freq.txt

# Optional beam velocity model for this molecule
vel  = 1000 100 3

# Notes:
# - If hess is unavailable, harmonic vibrational sampling is disabled.
# - Keep file names relative to this file for portability.
"""

    nh3_txt = """# NH3 molecule input template
# Required file links (update paths as needed)
xyz  = ammonia_geom.xyz
# hess is the full mass-weighted Cartesian Hessian in atomic units.
# Its eigenvalues should be omega^2; ICATS diagonalizes this matrix directly.
hess = ammonia_hessian.txt
w    = ammonia_freq.txt

# Optional beam velocity model for this molecule
vel  = 1000 100 3

# Notes:
# - If hess is unavailable, harmonic vibrational sampling is disabled.
# - Keep file names relative to this file for portability.
"""

    files = {
        "input.template": input_txt,
        "h2o_dat.template.txt": h2o_txt,
        "ammonia_dat.template.txt": nh3_txt,
    }

    written: List[str] = []
    for name, text in files.items():
        path = out / name
        if path.exists() and not overwrite:
            continue
        path.write_text(text)
        written.append(str(path))
    return written


__all__ = [
    "KNOWN_KEYS",
    "validate_molecule_file_refs",
    "load_and_validate",
    "write_run_log",
    "write_templates",
]
