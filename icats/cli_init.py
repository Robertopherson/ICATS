#!/usr/bin/env python3
import argparse
import os
import sys
import tempfile
from pathlib import Path


def _set_default_cache_dirs() -> None:
    user = os.environ.get("USER", "user")
    tmp = tempfile.gettempdir()
    numba_dir = os.path.join(tmp, "numba_cache_" + user)
    mpl_dir = os.path.join(tmp, "mpl_cache_" + user)
    os.makedirs(numba_dir, exist_ok=True)
    os.makedirs(mpl_dir, exist_ok=True)
    os.environ.setdefault("NUMBA_CACHE_DIR", numba_dir)
    os.environ.setdefault("MPLCONFIGDIR", mpl_dir)


def _ensure_repo_parent_on_path() -> None:
    # Keep script executable from any working directory.
    this_file = Path(__file__).resolve()
    repo_parent = str(this_file.parent.parent)
    if repo_parent not in sys.path:
        sys.path.insert(0, repo_parent)

REQUIRED_INPUT_KEYS = [
    ("mol", "Molecule definition lines. Expected twice (mi file), e.g. 'mol = 0 ammonia_dat.txt'."),
]

OPTIONAL_INPUT_KEYS = [
    ("Nsamp", "Number of generated samples."),
    ("maxj", "Use for Wang-Landau runs: total angular-momentum cap used in WL/J setup."),
    ("maxl", "Use for non-WL runs: orbital angular-momentum cap L."),
    ("maxb", "Maximum impact parameter (Angstrom in input, converted internally)."),
    ("Trot", "Rotational temperature."),
    ("Tvib", "Vibrational temperature."),
    ("Tvel", "Intermolecular velocity temperature or beam-centered velocity form."),
    ("workers", "Number of parallel workers."),
    ("seed", "Random seed."),
    ("seed-mode", "Seed strategy: fixed, time, per-worker."),
    ("run-mode", "Run behavior: fresh or continue. Existing wang.pkl is reused only if compatible."),
    ("wang", "Enable Wang-Landau weighting (True/False)."),
    ("wlmode", "Wang-Landau profile: fast, default, accurate."),
    ("wl-ff", "Override Wang-Landau initial f (>1)."),
    ("wl-nstep", "Override Wang-Landau nstep multiplier per wn bin."),
    ("wl-flatness", "Override Wang-Landau flatness criterion."),
    ("wl-wn-factor", "Override Wang-Landau wn density factor (wn ~ PeakJab * factor)."),
    ("wl-wn", "Override Wang-Landau wn bins directly (absolute)."),
    ("wl-tol", "Wang-Landau stopping tolerance on f."),
    ("wl-max-iter", "Maximum WL iterations (0 means no cap)."),
    ("wl-log-every", "WL status print period in iterations."),
    ("wl-angular-sampler", "WL angular sampler implementation: fast or legacy."),
    ("wl-audit-angular-sampler", "Compare fast and legacy WL angular samplers on initial draws (True/False)."),
    ("audit-initial-sample", "Verify generated sample bookkeeping against immediate coordinate analysis (True/False)."),
    ("audit-initial-energy-tol", "Initial-sample audit energy tolerance in eV."),
    ("audit-initial-angular-tol", "Optional angular-momentum magnitude tolerance; 0 disables angular checks."),
    ("progress", "Console verbosity: quiet, normal, verbose."),
    ("dry-run", "Validate/load setup but skip generation (True/False)."),
    ("check-input", "Check input/setup only and exit before sampling (True/False)."),
    ("save-frequency", "Worker checkpoint frequency in samples (0 disables)."),
    ("output-format", "Output format: xyzvel, npz, both."),
    ("units-out", "Output units: ang-fs, au."),
    ("write-templates", "Generate starter template files and exit."),
    ("continue", "Continue from cached pickles if available (True/False)."),
    ("rz", "Initial intermolecular Z distance (Angstrom)."),
    ("beam-angle", "Cross-beam angle in degrees."),
    ("fileout", "Output file prefix."),
    ("dirout", "Output directory."),
    ("printout", "Print flags (4 ints)."),
    ("plothist", "Generate histograms (True/False)."),
    ("keepinfo", "Store extra per-sample info (True/False)."),
    ("plotinit", "Number of pre-samples for distribution plotting."),
    ("ordist", "Orientation distribution function and parameters."),
    ("rot-param", "Rotation parameterization override."),
    ("phisample", "Sample orbital azimuthal coordinate phi (True/False)."),
    ("maxv", "Maximum vibrational state."),
    ("chi", "Azimuthal scattering angle."),
]


def _print_input_spec() -> None:
    print("Input specification for icats.init")
    print("")
    print("Required keys:")
    for key, desc in REQUIRED_INPUT_KEYS:
        print(f"  - {key}: {desc}")
    print("")
    print("Optional keys:")
    for key, desc in OPTIONAL_INPUT_KEYS:
        print(f"  - {key}: {desc}")


def main() -> int:
    _set_default_cache_dirs()
    _ensure_repo_parent_on_path()
    parser = argparse.ArgumentParser(
        prog="icats.init",
        description="Generate initial conditions from an icats input file.",
    )
    parser.add_argument("input_file", nargs="?", help="Input configuration file (e.g. input)")
    parser.add_argument(
        "--show-input-options",
        action="store_true",
        help="Print required and optional input-file keys, then exit.",
    )
    parser.add_argument(
        "--write-templates",
        action="store_true",
        help="Write starter input/molecule template files in current directory and exit.",
    )
    parser.add_argument(
        "--overwrite-templates",
        action="store_true",
        help="When used with --write-templates, overwrite existing template files.",
    )
    args = parser.parse_args()
    if args.show_input_options:
        _print_input_spec()
        return 0
    if args.write_templates:
        from icats.io import write_templates

        written = write_templates(".", overwrite=args.overwrite_templates)
        if written:
            print("Wrote templates:")
            for p in written:
                print(" -", p)
        else:
            print("No template files written (already exist). Use --overwrite-templates to replace.")
        return 0
    if not args.input_file:
        parser.error("input_file is required unless --show-input-options is used")

    from icats.iscattering import icats
    from icats.io import load_and_validate, write_run_log

    cfg = load_and_validate(args.input_file)
    for w in cfg["warnings"]:
        print("WARNING:", w)
    sc = icats()
    sc.ReadInput(args.input_file)
    # Flush parsed/setup information before sampling starts.
    write_run_log(args.input_file + ".logfile", sc.log)
    sc.GenSamples()
    write_run_log(args.input_file + ".logfile", sc.log)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
