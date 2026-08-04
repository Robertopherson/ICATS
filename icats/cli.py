#!/usr/bin/env python3
import argparse
import shutil
import subprocess
import sys
from pathlib import Path


TUTORIAL_ORDER = [
    "quickstart",
    "diatomic_n2_n2_fast",
    "mixed_h2o_n2",
    "methane_methane",
    "single_atom_he_he",
    "single_atom_diatom_he_n2",
    "polarized_orientation_he_no",
    "fixed_plane_atom_diatom_ar_no",
    "flat_l_atom_diatom_ar_no",
    "single_atom_diatom_he_n2_wl",
    "flat_l_atom_diatom_he_n2_wl",
    "flat_l_diatom_diatom_n2_n2_wl",
    "wang_landau_nh3_h2o",
    "paper_nh3_h2o_100k",
    "npz_output_co2_co2",
]

TUTORIALS = {
    "quickstart": {
        "desc": "NH3 + H2O minimal end-to-end (recommended first run)",
        "mol0": "ammonia_dat.txt",
        "mol1": "h2o_dat.txt",
        "nsamp": 10,
        "ntraj": 3,
        "steps": 200,
        "maxl": 182,
        "maxb": 16,
        "extra": [
            "wang = False",
            "phisample = True",
        ],
        "notes": [
            "Baseline tutorial with the smallest practical settings.",
            "Shows the standard flow: init -> cheap dynamics -> analysis.",
        ],
    },
    "diatomic_n2_n2_fast": {
        "desc": "N2 + N2 fast sanity check for environment/runtime",
        "mol0": "n2_dat.txt",
        "mol1": "n2_dat.txt",
        "nsamp": 20,
        "ntraj": 3,
        "steps": 150,
        "maxl": 174,
        "maxb": 16,
        "extra": [
            "wang = False",
            "phisample = True",
        ],
        "notes": [
            "Fastest physical template for quick checks.",
            "Useful to confirm libraries and command flow are healthy.",
        ],
    },
    "mixed_h2o_n2": {
        "desc": "H2O + N2 mixed rotor types, moderate sampling",
        "mol0": "h2o_dat.txt",
        "mol1": "n2_dat.txt",
        "nsamp": 40,
        "ntraj": 3,
        "steps": 200,
        "maxl": 187,
        "maxb": 16,
        "extra": [
            "wang = False",
            "phisample = True",
        ],
        "notes": [
            "Good intermediate template with hetero collision partners.",
            "Uses common flags you will likely keep for production runs.",
        ],
    },
    "methane_methane": {
        "desc": "CH4 + CH4 heavier symmetric-top test case",
        "mol0": "meth_dat.txt",
        "mol1": "meth_dat.txt",
        "nsamp": 30,
        "ntraj": 3,
        "steps": 200,
        "maxl": 99,
        "maxb": 16,
        "extra": [
            "wang = False",
            "phisample = True",
        ],
        "notes": [
            "Useful to test a heavier polyatomic pair.",
            "Same simple workflow, slightly richer internal structure.",
        ],
    },
    "single_atom_he_he": {
        "desc": "He + He single-atom edge case (no internal vib/rot modes)",
        "mol0": "he_dat.txt",
        "mol1": "he_dat.txt",
        "nsamp": 20,
        "ntraj": 3,
        "steps": 120,
        "maxl": 41,
        "maxb": 16,
        "extra": [
            "wang = False",
            "phisample = True",
        ],
        "notes": [
            "Exercises single-atom path with no molecular internal modes.",
            "Good edge-case sanity check for setup and analysis plumbing.",
        ],
    },
    "single_atom_diatom_he_n2": {
        "desc": "He + N2 atom-diatom template without Wang-Landau",
        "mol0": "he_dat.txt",
        "mol1": "n2_dat.txt",
        "nsamp": 30,
        "ntraj": 3,
        "steps": 150,
        "maxl": 59,
        "maxb": 16,
        "extra": [
            "wang = False",
            "phisample = True",
        ],
        "notes": [
            "Atom-diatom baseline; good before trying WL.",
            "Simple and fast while including one molecular rotor.",
        ],
    },
    "polarized_orientation_he_no": {
        "desc": "He + NO toy polarized-orientation PDF tutorial",
        "mol0": "he_dat.txt",
        "mol1": "no_polarized_dat.txt",
        "extra_files": ["orientation_pdfs.py"],
        "nsamp": 2000,
        "ntraj": 1,
        "steps": 120,
        "maxl": 80,
        "maxb": 8,
        "tvib": 0.0,
        "trot": 0.0,
        "printout": "0 1 0 0",
        "hist_by_default": True,
        "extra": [
            "vib-mode = rigid",
            "maxv = 0",
            "incoming-k = 10.0",
            "wang = False",
            "phisample = True",
        ],
        "notes": [
            "Demonstrates a user-supplied molecular orientation PDF.",
            "NO is kept rigid with Trot=0 so the orientation histogram is the main signal.",
            "The example PDF treats the NO body-z axis as a dipole in a tilted field in the ICATS scattering frame.",
            "The tilted field makes the sampled distribution depend on both alpha and beta.",
            "Run rd_tutorial_input/histograms/plot_polarization_check.py after icats.init to compare sampled angles with the expected trends.",
        ],
    },
    "fixed_plane_atom_diatom_ar_no": {
        "desc": "Ar + NO constrained atom-diatom setup with fixed b and impact plane",
        "mol0": "ar_dat.txt",
        "mol1": "no_dat.txt",
        "nsamp": 100,
        "ntraj": 1,
        "steps": 120,
        "maxl": 120,
        "maxb": 4.5,
        "tvib": 0.0,
        "trot": 0.0,
        "rz": 58.20949319933,
        "printout": "1 1 0 0",
        "extra": [
            "vib-mode = rigid",
            "maxv = 0",
            "incoming-k = 13.615392",
            "fixed-b = 4.5",
            "impact-phi = 0.0",
            "wang = False",
            "phisample = True",
        ],
        "notes": [
            "Constrained diagnostic setup: rigid NO, zero initial rotor angular momentum.",
            "Fixes impact parameter and lab impact-plane azimuth while sampling NO orientation.",
            "Writes combined out_full.xyz, out_full.vel, and out_full.info for inspection.",
            "PySCF MINDO/3 does not support Ar in this setup, so do not use run_cheap_dynamics.sh for this tutorial.",
            "Use these files as initial conditions for inspection or for an external dynamics/QM code.",
        ],
    },
    "flat_l_atom_diatom_ar_no": {
        "desc": "Ar + NO diagnostic setup with uniform L proposals",
        "mol0": "ar_dat.txt",
        "mol1": "no_dat.txt",
        "nsamp": 200,
        "ntraj": 1,
        "steps": 120,
        "maxl": 205,
        "maxb": 8.0,
        "tvib": 0.0,
        "trot": 0.0,
        "rz": 58.20949319933,
        "printout": "0 1 0 0",
        "extra": [
            "vib-mode = rigid",
            "maxv = 0",
            "incoming-k = 13.615392",
            "impact-phi = 0.0",
            "orbital-sampling = flat-l",
            "wang = False",
            "phisample = True",
        ],
        "notes": [
            "Diagnostic tutorial for non-geometric orbital proposals.",
            "Samples L uniformly instead of using the geometric P(L) proportional to L measure.",
            "With rigid NO and Trot=0, Jab=0 and J=L, so the histogram check is transparent.",
            "Reweight by L or b before using this ensemble for geometric cross-section averages.",
            "PySCF MINDO/3 does not support Ar in this setup, so do not use run_cheap_dynamics.sh for this tutorial.",
        ],
    },
    "single_atom_diatom_he_n2_wl": {
        "desc": "He + N2 atom-diatom template with Wang-Landau enabled",
        "mol0": "he_dat.txt",
        "mol1": "n2_dat.txt",
        "nsamp": 40,
        "ntraj": 3,
        "steps": 150,
        "maxl": 59,
        "maxb": 16,
        "extra": [
            "wang = True",
            "wlmode = fast",
            "run-mode = fresh",
            "phisample = True",
        ],
        "notes": [
            "Advanced atom-diatom variant to inspect WL effects.",
            "Use run_continue.sh with run-mode=continue for accumulation.",
        ],
    },
    "flat_l_atom_diatom_he_n2_wl": {
        "desc": "He + N2 Wang-Landau diagnostic with uniform L proposals",
        "mol0": "he_dat.txt",
        "mol1": "n2_dat.txt",
        "nsamp": 40,
        "ntraj": 3,
        "steps": 150,
        "maxl": 59,
        "maxb": 16,
        "extra": [
            "wang = True",
            "wlmode = default",
            "wl-target = flat-j",
            "wl-j-range = 60",
            "wl-j-bins = 60",
            "wl-l-cap = 59",
            "wl-flatness = 0.95",
            "wl-nstep = 1000",
            "run-mode = fresh",
            "orbital-sampling = flat-l",
            "phisample = True",
        ],
        "notes": [
            "Diagnostic companion to single_atom_diatom_he_n2_wl.",
            "Uses the same one-dimensional WL-on-J machinery, but proposes L uniformly.",
            "Uses wl-target = flat-j so WL does not add the linear 2J+1 target factor.",
            "Useful for checking whether flat-L proposals improve low-L coverage before reweighting.",
            "The final L and J marginals can be only approximately flat because the one-dimensional WL acceptance still couples L, J, and Jab.",
            "Uses wl-j-range, wl-j-bins, and wl-l-cap to make the WL diagnostic range explicit.",
            "Inspect the WL umbrella plus sampled L, J, and Jab histograms before trusting production settings.",
        ],
    },
    "flat_l_diatom_diatom_n2_n2_wl": {
        "desc": "N2 + N2 flat-L Wang-Landau diagnostic with two molecular rotors",
        "mol0": "n2_dat.txt",
        "mol1": "n2_dat.txt",
        "nsamp": 40,
        "ntraj": 3,
        "steps": 150,
        "maxl": 100,
        "maxb": 10,
        "extra": [
            "wang = True",
            "wlmode = default",
            "wl-target = flat-j",
            "wl-j-range = 100",
            "wl-j-bins = 50",
            "wl-j-min = 5.0",
            "wl-low-j-scale = 0.15",
            "wl-l-cap = 100",
            "wl-flatness = 0.95",
            "wl-nstep = 1500",
            "run-mode = fresh",
            "orbital-sampling = flat-l",
            "phisample = True",
        ],
        "notes": [
            "Diatom-diatom companion to the He + N2 flat-L WL diagnostic.",
            "Uses two molecular rotors, so Jab has more angular freedom than in atom-diatom tests.",
            "This is intended as the clearer flat-J WL tutorial; He + N2 is kept as a stress case.",
            "Uses wl-j-min = 5.0 and wl-low-j-scale = 0.15 so rare near-zero-J cancellation bins are not forced too strongly.",
            "Inspect L, J, Jab, and cos(theta_L,Jab) diagnostics after running.",
            "If you change WL range, bins, target, low-J boundary, or L cap, move rd_tutorial_input/wang.pkl before rerunning.",
        ],
    },
    "wang_landau_nh3_h2o": {
        "desc": "NH3 + H2O with Wang-Landau weighting enabled",
        "mol0": "ammonia_dat.txt",
        "mol1": "h2o_dat.txt",
        "nsamp": 60,
        "ntraj": 3,
        "steps": 200,
        "maxl": 182,
        "maxb": 16,
        "extra": [
            "wang = True",
            "wlmode = normal",
            "wl-flatness = 0.90",
            "wl-tol = 1.00001",
            "run-mode = fresh",
            "phisample = True",
        ],
        "notes": [
            "Template to demonstrate WL workflow and restart behavior.",
            "For continuation, switch run-mode to continue in tutorial_input.txt.",
        ],
    },
    "paper_nh3_h2o_100k": {
        "desc": "NH3 + H2O 100k validation ensemble for manuscript Figures 6, 7, and 9",
        "mol0": "paper_ammonia_dat.txt",
        "mol1": "paper_h2o_dat.txt",
        "nsamp": 100010,
        "workers": 4,
        "ntraj": 1,
        "steps": 200,
        "maxl": 70,
        "maxb": 16,
        "tvib": 500.0,
        "trot": 500.0,
        "rz": 20,
        "printout": "0 0 0 0",
        "extra": [
            "beam-angle = 90.0",
            "wang = True",
            "wlmode = normal",
            "wl-ff = 1.105170185988091",
            "wl-nstep = 500",
            "wl-flatness = 0.90",
            "wl-wn-factor = 4.0",
            "wl-tol = 1.00001",
            "wl-angular-sampler = fast",
            "wl-audit-angular-sampler = True",
            "run-mode = fresh",
            "phisample = True",
        ],
        "notes": [
            "Paper-validation tutorial rather than a short introductory calculation.",
            "Uses 100,010 collision samples, 500 K vibrational and rotational temperatures, and a 90 degree crossed-beam geometry.",
            "The NH3 and H2O beam-speed distributions are centred at 600 and 800 m/s, respectively, with 100 m/s FWHM.",
            "Uses an orbital cap of 70, giving the J/L support and upper-edge roll-off shown in manuscript Figure 9.",
            "Freezes the normal WL profile at ff=exp(0.1), 500 trial steps per active bin, flatness 0.90, and tolerance 1.00001.",
            "Run export_paper_histogram_data.py after icats.init to write compact CSV data for manuscript Figures 6, 7, and 9.",
            "Run plot_paper_histograms.py to turn those tables into PNG and PDF validation composites.",
            "Figure 8 is a separate fixed-J asymmetric-top benchmark and is not generated from this collision ensemble.",
            "This calculation can take substantial time; run_paper_slurm.sh is supplied for a configured cluster environment.",
        ],
    },
    "npz_output_co2_co2": {
        "desc": "CO2 + CO2 with dual output format (xyzvel + npz)",
        "mol0": "co2_dat.txt",
        "mol1": "co2_dat.txt",
        "nsamp": 30,
        "ntraj": 3,
        "steps": 180,
        "maxl": 273,
        "maxb": 16,
        "extra": [
            "wang = False",
            "phisample = True",
            "output-format = both",
            "units-out = ang-fs",
        ],
        "notes": [
            "Shows output-format flag usage while keeping run size small.",
            "Generates both plain text outputs and npz trajectory files.",
        ],
    },
}


def _root_dir() -> Path:
    return Path(__file__).resolve().parent.parent


def _examples_dir() -> Path:
    return _root_dir() / "icats" / "examples"


def _precomputed_dir() -> Path:
    return _root_dir() / "icats" / "tutorials" / "precomputed"

def _dyn_runner_path() -> Path:
    return _root_dir() / "icats" / "run_tutorial_dyn.py"


def _write_file(path: Path, text: str, executable: bool = False) -> None:
    path.write_text(text)
    if executable:
        path.chmod(0o755)

def _write_initial_histogram_helper(out_dir: Path, run_tag: str = "tutorial_input") -> None:
    hist_dir = out_dir / f"rd_{run_tag}" / "histograms"
    hist_dir.mkdir(parents=True, exist_ok=True)

    _write_file(
        hist_dir / "plot_initial.sh",
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "SCRIPT_DIR=\"$(cd \"$(dirname \"${BASH_SOURCE[0]}\")\" && pwd)\"\n"
        "HROOT=\"${HIST_ROOT:-$SCRIPT_DIR}\"\n"
        "export MPLCONFIGDIR=\"${MPLCONFIGDIR:-/tmp/mpl_cache_${USER:-user}}\"\n"
        "mkdir -p \"$MPLCONFIGDIR\"\n"
        "mkdir -p \"$HROOT/plots/initial\"\n"
        "find \"$HROOT/initial\" -type f -name 'hist_*.py' 2>/dev/null | sort | while read -r f; do\n"
        "  stem=\"$(basename \"${f%.py}\")\"\n"
        "  python \"$f\" --no-show --outfile \"$HROOT/plots/initial/${stem}\"\n"
        "done\n"
        "echo \"Initial histogram plots written to $HROOT/plots/initial/\"\n",
        executable=True,
    )


def _write_histogram_helpers(
    out_dir: Path,
    run_tag: str = "tutorial_input",
    include_paper: bool = False,
) -> None:
    _write_initial_histogram_helper(out_dir, run_tag=run_tag)
    hist_dir = out_dir / f"rd_{run_tag}" / "histograms"
    exporter = r'''#!/usr/bin/env python3
"""Export the compact validation data used by manuscript Figures 6, 7, and 9."""

from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path

import numpy as np


def write_columns(path: Path, header: str, *columns) -> None:
    arrays = [np.asarray(column) for column in columns]
    if not arrays or any(len(column) != len(arrays[0]) for column in arrays):
        raise ValueError("CSV columns must have equal lengths: " + str(path))
    data = np.column_stack(arrays)
    np.savetxt(path, data, delimiter=",", header=header, comments="")


def write_histogram(path: Path, values, bins) -> None:
    counts, edges = np.histogram(np.asarray(values, dtype=float), bins=bins)
    widths = np.diff(edges)
    density = counts / max(1, counts.sum()) / widths
    write_columns(path, "bin_left,bin_right,count,density", edges[:-1], edges[1:], counts, density)


def state_histogram(path: Path, values) -> None:
    states, counts = np.unique(np.asarray(values, dtype=int), return_counts=True)
    write_columns(path, "state,count,probability", states, counts, counts / counts.sum())


def load_wang(path: Path):
    with path.open("rb") as handle:
        payload = pickle.load(handle)
    if isinstance(payload, dict):
        return payload["iwld"], payload["td"], payload.get("metadata", {})
    uu, iwld, td = payload
    return iwld, td, {"format": "legacy tuple"}


def projection_arrays(rotor, j_value: int):
    entry = rotor.get(j_value, rotor.get(str(j_value), {}))
    return (
        np.asarray(entry.get("jz", []), dtype=float),
        np.asarray(entry.get("sjz", []), dtype=float),
        np.asarray(entry.get("qjz", []), dtype=float),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", default="rd_tutorial_input")
    parser.add_argument("--output-dir", default="paper_histogram_data")
    parser.add_argument("--vibrational-seed", type=int, default=20260619)
    parser.add_argument("--vibrational-samples", type=int, default=100000)
    args = parser.parse_args()

    run_dir = Path(args.run_dir)
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    data_path = run_dir / "dat_''' + run_tag + r'''.txt.pkl"
    wang_path = run_dir / "wang.pkl"
    if not data_path.exists():
        raise SystemExit("Missing sampled data: " + str(data_path) + ". Run icats.init first.")
    if not wang_path.exists():
        raise SystemExit("Missing Wang-Landau umbrella: " + str(wang_path) + ". Run icats.init first.")

    with data_path.open("rb") as handle:
        sampled = pickle.load(handle)

    # Figure 6 is a one-mode component validation. All harmonic modes share this
    # reduced distribution after mass-frequency scaling, so it is generated
    # separately at the manuscript's NH3 umbrella-mode frequency and temperature.
    nsamp = len(sampled["orb"]["sJ"])
    vib_nsamp = args.vibrational_samples
    if vib_nsamp <= 0:
        raise SystemExit("--vibrational-samples must be positive")
    rng = np.random.default_rng(args.vibrational_seed)
    temp_k = 500.0
    freq_cm = 652.27
    kbt_cm = 0.69503476 * temp_k
    states = np.arange(25)
    probabilities = np.exp(-freq_cm * states / kbt_cm)
    probabilities /= probabilities.sum()
    vib_n = rng.choice(states, size=vib_nsamp, p=probabilities)
    r2 = rng.gamma(shape=2 * vib_n + 1, scale=1.0)
    phase = rng.random(vib_nsamp) * 2.0 * np.pi
    vib_q = np.sqrt(r2) * np.cos(phase)
    vib_p = np.sqrt(r2) * np.sin(phase)
    vib_energy = 0.5 * (vib_q**2 + vib_p**2)
    write_columns(
        output / "figure6_mode1_vibrations.csv",
        "state,Q,P,dimensionless_energy",
        vib_n,
        vib_q,
        vib_p,
        vib_energy,
    )
    state_histogram(output / "figure6_state_populations.csv", vib_n)
    for state in (1, 2, 3):
        write_histogram(output / f"figure6_n{state}_Q_histogram.csv", vib_q[vib_n == state], 100)

    # Figure 7 uses the sampled molecular J distributions and the J=5
    # body-fixed projection records. For asymmetric tops, `jz` is the selected
    # eigenstate projection centre, whereas `sjz` is the vector-model projection
    # reconstructed from the sampled Cartesian coordinates.
    nh3_rot = sampled["rot"]["m0"]
    h2o_rot = sampled["rot"]["m1"]
    nh3_j = np.asarray(nh3_rot["J"], dtype=int)
    h2o_j = np.asarray(h2o_rot["J"], dtype=int)
    write_columns(output / "figure7_rotor_J_samples.csv", "sample,NH3_J,H2O_J", np.arange(nsamp), nh3_j, h2o_j)
    state_histogram(output / "figure7_NH3_J_histogram.csv", nh3_j)
    state_histogram(output / "figure7_H2O_J_histogram.csv", h2o_j)
    for label, rotor in (("NH3", nh3_rot), ("H2O", h2o_rot)):
        jz, analysed_jz, rms_jz = projection_arrays(rotor, 5)
        nproj = min(len(jz), len(analysed_jz), len(rms_jz))
        write_columns(
            output / f"figure7_{label}_J5_projection.csv",
            "projection_centre_Jz,reconstructed_vector_Jz,quantum_rms_Jz",
            jz[:nproj],
            analysed_jz[:nproj],
            rms_jz[:nproj],
        )
        if nproj:
            write_histogram(output / f"figure7_{label}_J5_sampled_projection_histogram.csv", analysed_jz[:nproj], 30)

    # Figure 9 is taken directly from the accepted bimolecular ensemble.
    orb = sampled["orb"]
    relative_velocity = np.asarray(sampled["vel"]["ivel"], dtype=float)
    write_columns(
        output / "figure9_intermolecular_samples.csv",
        "sample,J,L,Jab,relative_velocity_m_per_s",
        np.arange(nsamp),
        orb["sJ"],
        orb["sL"],
        orb["sJab"],
        relative_velocity,
    )
    angular_max = max(float(np.max(orb["sJ"])), float(np.max(orb["sL"])))
    angular_edges = np.arange(0.0, 4.0 * np.ceil(angular_max / 4.0) + 8.0, 4.0)
    write_histogram(output / "figure9_J_histogram.csv", orb["sJ"], angular_edges)
    write_histogram(output / "figure9_L_histogram.csv", orb["sL"], angular_edges)
    write_histogram(output / "figure9_relative_velocity_histogram.csv", relative_velocity, 100)
    omega, acceptance, wang_metadata = load_wang(wang_path)
    write_columns(
        output / "figure9_wang_landau.csv",
        "J,estimated_Omega_J,acceptance_weight",
        np.arange(len(acceptance)),
        np.pad(np.asarray(omega), (0, max(0, len(acceptance) - len(omega))), constant_values=np.nan)[:len(acceptance)],
        acceptance,
    )

    metadata = {
        "source_sample_pickle": str(data_path),
        "source_wang_pickle": str(wang_path),
        "accepted_samples": nsamp,
        "figure6_samples": vib_nsamp,
        "molecule_order": ["NH3", "H2O"],
        "temperature_K": 500.0,
        "beam_centres_m_per_s": [600.0, 800.0],
        "beam_fwhm_m_per_s": [100.0, 100.0],
        "beam_speed_power": [3, 3],
        "beam_angle_degrees": 90.0,
        "rotational_constants_hartree": {
            "NH3_C": 5.0407924631164556e-05,
            "NH3_Cz": -2.3489116084713895e-05,
            "H2O_Cbar": 8.533146923414483e-05,
        },
        "figure6_frequency_cm_inverse": freq_cm,
        "figure6_seed": args.vibrational_seed,
        "wang_landau": wang_metadata,
        "scope": "Histogram data for manuscript Figures 6, 7, and 9; Figure 8 is a separate fixed-J benchmark.",
    }
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print("Wrote paper validation data to", output.resolve())


if __name__ == "__main__":
    main()
'''
    if include_paper:
        _write_file(out_dir / "export_paper_histogram_data.py", exporter, executable=True)

    paper_plotter = r'''#!/usr/bin/env python3
"""Plot the compact validation data for manuscript Figures 6, 7, and 9."""

from __future__ import annotations

import argparse
import json
from math import comb, gamma
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def trapezoid(y, x):
    implementation = getattr(np, "trapezoid", None)
    if implementation is None:
        implementation = np.trapz
    return implementation(y, x)


def read_csv(path: Path):
    if not path.exists():
        raise SystemExit("Missing " + str(path) + ". Run export_paper_histogram_data.py first.")
    return np.genfromtxt(path, delimiter=",", names=True)


def normalize(y, x):
    area = trapezoid(y, x)
    return y / area if area > 0.0 else y


def hermite_values(n: int, x):
    if n == 0:
        return np.ones_like(x)
    if n == 1:
        return 2.0 * x
    hm2 = np.ones_like(x)
    hm1 = 2.0 * x
    for k in range(1, n):
        hn = 2.0 * x * hm1 - 2.0 * k * hm2
        hm2, hm1 = hm1, hn
    return hn


def hermite_density(n: int, x):
    return normalize(hermite_values(n, x) ** 2 * np.exp(-(x**2)), x)


def leading_wigner_projection(n: int, x):
    order = 2 * n
    poly = np.zeros_like(x)
    for j in range(order + 1):
        poly += comb(order, j) * gamma(j + 0.5) * x ** (2 * (order - j))
    return normalize(np.exp(-(x**2)) * poly / (np.pi * gamma(order + 1)), x)


def plot_figure6(data_dir: Path, output: Path) -> None:
    samples = read_csv(data_dir / "figure6_mode1_vibrations.csv")
    grid = np.linspace(-5.5, 5.5, 1400)
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.2))
    for ax, state in zip(axes, (1, 2, 3)):
        q = samples["Q"][samples["state"].astype(int) == state]
        if len(q):
            ax.hist(q, bins=100, density=True, color="tab:blue", edgecolor="black", linewidth=0.3, alpha=0.55)
        ax.plot(grid, leading_wigner_projection(state, grid), color="tab:red", lw=1.8, label="leading-Wigner")
        ax.plot(grid, hermite_density(state, grid), color="black", lw=1.5, label="HO density")
        label = rf"$n={state}$" + "\n" + rf"$N_{{MC}}={len(q)}$"
        ax.text(0.04, 0.95, label, transform=ax.transAxes, va="top")
        ax.set_xlim(-5.0, 5.0)
        ax.set_xlabel(r"$Q$")
    axes[0].set_ylabel(r"$\rho(Q)$")
    axes[-1].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(output / "figure6_validation.png", dpi=180)
    fig.savefig(output / "figure6_validation.pdf")
    plt.close(fig)


def rotational_target(states, rotational_constant, temperature):
    kb_hartree_per_k = 3.166811563e-6
    weights = (2.0 * states + 1.0) * np.exp(
        -rotational_constant * states * (states + 1.0) / (kb_hartree_per_k * temperature)
    )
    return weights / weights.sum()


def state_bars(ax, table, title, rotational_constant, temperature):
    ax.bar(table["state"], table["probability"], width=0.9, color="tab:blue", edgecolor="black", linewidth=0.35)
    states = np.asarray(table["state"], dtype=float)
    ax.plot(states, rotational_target(states, rotational_constant, temperature), color="tab:red", lw=1.5)
    ax.set_title(title)
    ax.set_xlabel(r"$J$")
    ax.set_ylabel("probability")


def projection_panel(ax, data_dir: Path, molecule: str, temperature: float, nh3_cz: float):
    values = read_csv(data_dir / f"figure7_{molecule}_J5_projection.csv")
    if molecule == "NH3":
        bins = np.arange(-5.5, 6.5, 1.0)
        projection = np.asarray(values["reconstructed_vector_Jz"])
        ax.hist(projection, bins=bins, density=True, alpha=0.65, color="tab:blue", edgecolor="black", linewidth=0.35)
        k = np.arange(-5.0, 6.0)
        kb_hartree_per_k = 3.166811563e-6
        target = np.exp(-nh3_cz * k**2 / (kb_hartree_per_k * temperature))
        target /= target.sum()
        ax.plot(k, target, color="tab:red", lw=1.5)
        ax.set_xlabel(r"body-fixed $K$")
    else:
        bins = np.linspace(-np.sqrt(30.0), np.sqrt(30.0), 12)
        projection = np.asarray(values["reconstructed_vector_Jz"])
        ax.hist(projection, bins=bins, density=True, alpha=0.65, color="tab:blue", edgecolor="black", linewidth=0.35)
        ax.set_xlabel("selected body-fixed projection")
    ax.set_title(molecule + r" projection, $J=5$")
    ax.set_ylabel("density")


def plot_figure7(data_dir: Path, output: Path, metadata) -> None:
    temperature = float(metadata["temperature_K"])
    constants = metadata["rotational_constants_hartree"]
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.8))
    state_bars(axes[0, 0], read_csv(data_dir / "figure7_NH3_J_histogram.csv"), r"NH$_3$ rotational states", constants["NH3_C"], temperature)
    projection_panel(axes[0, 1], data_dir, "NH3", temperature, constants["NH3_Cz"])
    state_bars(axes[1, 0], read_csv(data_dir / "figure7_H2O_J_histogram.csv"), r"H$_2$O rotational states", constants["H2O_Cbar"], temperature)
    projection_panel(axes[1, 1], data_dir, "H2O", temperature, constants["NH3_Cz"])
    fig.tight_layout()
    fig.savefig(output / "figure7_validation.png", dpi=180)
    fig.savefig(output / "figure7_validation.pdf")
    plt.close(fig)


def histogram_bars(ax, table, title, xlabel):
    centres = 0.5 * (table["bin_left"] + table["bin_right"])
    widths = table["bin_right"] - table["bin_left"]
    ax.bar(centres, table["density"], width=widths, color="tab:blue", edgecolor="black", linewidth=0.35, alpha=0.65)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("density")
    return centres


def weighted_beam_speed(v, centre, fwhm, power):
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    return np.where(v >= 0.0, v**power * np.exp(-0.5 * ((v - centre) / sigma) ** 2), 0.0)


def relative_speed_target(grid, metadata):
    centres = metadata["beam_centres_m_per_s"]
    fwhm = metadata["beam_fwhm_m_per_s"]
    powers = metadata["beam_speed_power"]
    theta = np.linspace(0.0, 0.5 * np.pi, 2401)
    normalization_grid = np.linspace(0.0, max(centres) + 10.0 * max(fwhm), 30001)
    norm_a = trapezoid(weighted_beam_speed(normalization_grid, centres[0], fwhm[0], powers[0]), normalization_grid)
    norm_b = trapezoid(weighted_beam_speed(normalization_grid, centres[1], fwhm[1], powers[1]), normalization_grid)
    va = grid[:, None] * np.cos(theta)[None, :]
    vb = grid[:, None] * np.sin(theta)[None, :]
    density = grid * trapezoid(
        weighted_beam_speed(va, centres[0], fwhm[0], powers[0])
        * weighted_beam_speed(vb, centres[1], fwhm[1], powers[1]),
        theta,
    ) / (norm_a * norm_b)
    return normalize(density, grid)


def plot_figure9(data_dir: Path, output: Path, metadata) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 6.8))
    l_table = read_csv(data_dir / "figure9_L_histogram.csv")
    j_table = read_csv(data_dir / "figure9_J_histogram.csv")
    j_grid = histogram_bars(axes[0, 0], j_table, "Total angular momentum", r"$J$")
    l_grid = histogram_bars(axes[0, 1], l_table, "Orbital angular momentum", r"$L$")
    for ax, grid, label in ((axes[0, 0], j_grid, r"$1+2J$"), (axes[0, 1], l_grid, r"$1+2L$")):
        target = 1.0 + 2.0 * grid
        target /= trapezoid(target, grid)
        ax.plot(grid, target, color="tab:red", lw=1.5, label=label)
        ax.legend(frameon=False, fontsize=8)
    v_table = read_csv(data_dir / "figure9_relative_velocity_histogram.csv")
    v_grid = histogram_bars(axes[1, 0], v_table, "Relative beam speed", r"$|v_{AB}|$ / m s$^{-1}$")
    axes[1, 0].plot(v_grid, relative_speed_target(v_grid, metadata), color="tab:red", lw=1.5, label="analytic target")
    axes[1, 0].legend(frameon=False, fontsize=8)
    wl = read_csv(data_dir / "figure9_wang_landau.csv")
    learned = np.isfinite(wl["estimated_Omega_J"])
    axes[1, 1].bar(
        wl["J"][learned],
        wl["acceptance_weight"][learned],
        width=0.9,
        color="tab:blue",
        edgecolor="black",
        linewidth=0.35,
    )
    axes[1, 1].set_title("Wang-Landau acceptance weight")
    axes[1, 1].set_xlabel(r"$J$")
    axes[1, 1].set_ylabel(r"$W(J)$")
    fig.tight_layout()
    fig.savefig(output / "figure9_validation.png", dpi=180)
    fig.savefig(output / "figure9_validation.pdf")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", default="paper_histogram_data")
    parser.add_argument("--output-dir", default="paper_histogram_plots")
    args = parser.parse_args()
    data_dir = Path(args.data_dir)
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    metadata = json.loads((data_dir / "metadata.json").read_text())
    plot_figure6(data_dir, output)
    plot_figure7(data_dir, output, metadata)
    plot_figure9(data_dir, output, metadata)
    print("Wrote validation plots for", metadata["accepted_samples"], "accepted samples to", output.resolve())


if __name__ == "__main__":
    main()
'''
    if include_paper:
        _write_file(out_dir / "plot_paper_histograms.py", paper_plotter, executable=True)

    if include_paper:
        _write_file(
            out_dir / "run_paper_slurm.sh",
        "#!/usr/bin/env bash\n"
        "#SBATCH --job-name=icats-paper-100k\n"
        "#SBATCH --export=ALL\n"
        "#SBATCH --partition=nodes\n"
        "#SBATCH --nodes=1\n"
        "#SBATCH --ntasks=1\n"
        "#SBATCH --cpus-per-task=16\n"
        "#SBATCH --mem=32G\n"
        "#SBATCH --time=12:00:00\n"
        "#SBATCH --output=icats-paper-%j.out\n"
        "#SBATCH --error=icats-paper-%j.err\n\n"
        "set -euo pipefail\n"
        "source \"${flight_ROOT:-/opt/flight}/etc/setup.sh\" >/dev/null 2>&1 || true\n"
        "if [ -n \"${SLURM_SUBMIT_DIR:-}\" ]; then\n"
        "  cd \"$SLURM_SUBMIT_DIR\"\n"
        "else\n"
        "  cd \"$(dirname \"${BASH_SOURCE[0]}\")\"\n"
        "fi\n"
        "export OMP_NUM_THREADS=1\n"
        "export OPENBLAS_NUM_THREADS=1\n"
        "export MKL_NUM_THREADS=1\n"
        "export NUMEXPR_NUM_THREADS=1\n"
        "export MPLBACKEND=Agg\n"
        "export NUMBA_CACHE_DIR=\"${NUMBA_CACHE_DIR:-/tmp/numba_cache_${USER:-user}}\"\n"
        "export MPLCONFIGDIR=\"${MPLCONFIGDIR:-/tmp/mpl_cache_${USER:-user}}\"\n"
        "mkdir -p \"$NUMBA_CACHE_DIR\" \"$MPLCONFIGDIR\"\n"
        "ICATS_ENV=\"${ICATS_ENV:-icats_clean}\"\n"
        "if command -v conda >/dev/null 2>&1; then\n"
        "  set +u\n"
        "  eval \"$(conda shell.bash hook)\"\n"
        "  conda activate \"$ICATS_ENV\"\n"
        "  set -u\n"
        "fi\n"
        "if ! command -v icats.init >/dev/null 2>&1; then\n"
        "  echo \"icats.init is not available. Activate the environment containing the frozen ICATS release.\" >&2\n"
        "  exit 1\n"
        "fi\n"
        "WORKERS=\"${SLURM_CPUS_PER_TASK:-1}\"\n"
        "sed -i -E \"s/^[[:space:]]*workers[[:space:]]*=.*/workers = ${WORKERS}/\" tutorial_input.txt\n"
        "printf 'ICATS commit: '; python -c 'import pathlib, subprocess; p=pathlib.Path(__import__(\"icats\").__file__).resolve().parent.parent; print(subprocess.run([\"git\", \"-C\", str(p), \"rev-parse\", \"HEAD\"], text=True, capture_output=True).stdout.strip() or p)'\n"
        "echo \"ICATS workers: $WORKERS\"\n"
        "icats.init tutorial_input.txt\n"
        "python export_paper_histogram_data.py --vibrational-samples 100000\n"
        "python plot_paper_histograms.py\n",
            executable=True,
        )

    _write_file(
        hist_dir / "plot_sampled.sh",
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "SCRIPT_DIR=\"$(cd \"$(dirname \"${BASH_SOURCE[0]}\")\" && pwd)\"\n"
        "HROOT=\"${HIST_ROOT:-$SCRIPT_DIR}\"\n"
        "export MPLCONFIGDIR=\"${MPLCONFIGDIR:-/tmp/mpl_cache_${USER:-user}}\"\n"
        "mkdir -p \"$MPLCONFIGDIR\"\n"
        "mkdir -p \"$HROOT/plots/sampled\"\n"
        "find \"$HROOT/sampled\" -type f -name 'hist_*.py' 2>/dev/null | sort | while read -r f; do\n"
        "  stem=\"$(basename \"${f%.py}\")\"\n"
        "  python \"$f\" --no-show --outfile \"$HROOT/plots/sampled/${stem}\"\n"
        "done\n"
        "echo \"Sampled histogram plots written to $HROOT/plots/sampled/\"\n",
        executable=True,
    )

    _write_file(
        hist_dir / "plot_orbital_jljab.sh",
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "SCRIPT_DIR=\"$(cd \"$(dirname \"${BASH_SOURCE[0]}\")\" && pwd)\"\n"
        "HROOT=\"${HIST_ROOT:-$SCRIPT_DIR}\"\n"
        "export MPLCONFIGDIR=\"${MPLCONFIGDIR:-/tmp/mpl_cache_${USER:-user}}\"\n"
        "mkdir -p \"$MPLCONFIGDIR\"\n"
        "mkdir -p \"$HROOT/plots/sampled\"\n"
        "for metric in sl sj sjab cosljab_thet; do\n"
        "  f=\"$HROOT/sampled/system/hist_sam_sys_orb_${metric}.py\"\n"
        "  if [ ! -f \"$f\" ]; then\n"
        "    if [ \"$metric\" = cosljab_thet ]; then\n"
        "      echo \"Skipping missing optional $f (no nonzero L/Jab angles recorded).\" >&2\n"
        "      continue\n"
        "    fi\n"
        "    echo \"Missing $f. Run icats.init with hist_sampled = True first.\" >&2\n"
        "    exit 1\n"
        "  fi\n"
        "  python \"$f\" --no-show --outfile \"$HROOT/plots/sampled/hist_sam_sys_orb_${metric}\"\n"
        "done\n"
        "echo \"Orbital J/L/Jab/cos(theta) plots written to $HROOT/plots/sampled/\"\n",
        executable=True,
    )

    _write_file(
        hist_dir / "plot_orbital_correlation.sh",
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "SCRIPT_DIR=\"$(cd \"$(dirname \"${BASH_SOURCE[0]}\")\" && pwd)\"\n"
        "HROOT=\"${HIST_ROOT:-$SCRIPT_DIR}\"\n"
        "RROOT=\"$(dirname \"$HROOT\")\"\n"
        "CONV=\"$RROOT/convergence\"\n"
        "export MPLCONFIGDIR=\"${MPLCONFIGDIR:-/tmp/mpl_cache_${USER:-user}}\"\n"
        "mkdir -p \"$MPLCONFIGDIR\"\n"
        "mkdir -p \"$HROOT/plots/sampled\"\n"
        "f=\"$HROOT/sampled/system/hist_sam_sys_orb_cosljab_thet.py\"\n"
        "if [ -f \"$f\" ]; then\n"
        "  python \"$f\" --no-show --outfile \"$HROOT/plots/sampled/hist_sam_sys_orb_cosljab_thet\"\n"
        "else\n"
        "  echo \"Missing optional $f (no nonzero L/Jab angles recorded).\" >&2\n"
        "fi\n"
        "if [ -f \"$CONV/costheta_j_vs_costheta_plot.py\" ]; then\n"
        "  (cd \"$CONV\" && python costheta_j_vs_costheta_plot.py --no-show --outfile \"$HROOT/plots/sampled/costheta_j_vs_costheta\")\n"
        "else\n"
        "  echo \"Missing $CONV/costheta_j_vs_costheta_plot.py. Run icats.init first.\" >&2\n"
        "fi\n"
        "if [ -f \"$CONV/wl_costheta_last_iter_j_vs_costheta_plot.py\" ]; then\n"
        "  (cd \"$CONV\" && python wl_costheta_last_iter_j_vs_costheta_plot.py --no-show --outfile \"$HROOT/plots/sampled/wl_costheta_last_iter_j_vs_costheta\")\n"
        "fi\n"
        "echo \"Orbital correlation plots written to $HROOT/plots/sampled/\"\n",
        executable=True,
    )

    _write_file(
        hist_dir / "plot_all.sh",
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "cd \"$(dirname \"${BASH_SOURCE[0]}\")\"\n"
        "./plot_initial.sh\n"
        "./plot_sampled.sh\n",
        executable=True,
    )

    _write_file(
        hist_dir / "plot_compare_pairs.py",
        """#!/usr/bin/env python3
\"\"\"Generate init-vs-sampled comparison plots and stats.

Default behavior:
- Compare only intermolecular/system quantities (J, L, vel, and when available b, phi).
- Skip rotational/vibrational pair rendering unless explicitly requested.
- Skip unmatched rendering unless requested.
\"\"\"
from __future__ import annotations

import argparse
import ast
import os
import re
import subprocess
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("HIST_ROOT", str(SCRIPT_DIR))).resolve()
INI_DIR = ROOT / "initial"
SAM_DIR = ROOT / "sampled"
OUT = ROOT / "plots" / "compare"
OUT_INI = OUT / "init_single"
OUT_SAM = OUT / "sam_single"
OUT_PAIR = OUT / "pairs"
OUT_UNM_INI = OUT / "unmatched_init"
OUT_UNM_SAM = OUT / "unmatched_sam"
for d in (OUT, OUT_INI, OUT_SAM, OUT_PAIR, OUT_UNM_INI, OUT_UNM_SAM):
    d.mkdir(parents=True, exist_ok=True)

ap = argparse.ArgumentParser()
ap.add_argument("--include-unmatched", action="store_true",
                help="Also render unmatched initial/sampled histograms after paired comparisons.")
ap.add_argument("--include-rot", action="store_true",
                help="Include molecule rotational pair comparisons (off by default).")
ap.add_argument("--include-vib", action="store_true",
                help="Include molecule vibrational pair comparisons (off by default).")
args = ap.parse_args()

rx = re.compile(r"^hist_(ini|sam)_([^_]+)_(.+)\\.py$")

def parse_scripts(base: Path):
    out = {}
    for p in sorted(base.rglob("hist_*.py")):
        m = rx.match(p.name)
        if not m:
            continue
        _stage, scope, metric = m.groups()
        out[(scope, metric)] = p
    return out

def render_hist(script_path: Path, out_base: Path):
    subprocess.run(
        ["python", str(script_path), "--no-show", "--outfile", str(out_base)],
        check=True,
        cwd=str(ROOT),
    )
    return out_base.with_suffix(".png")

def side_by_side(left_png: Path, right_png: Path, title: str, out_png: Path):
    li = mpimg.imread(left_png)
    ri = mpimg.imread(right_png)
    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    axs[0].imshow(li)
    axs[0].axis("off")
    axs[0].set_title("Initial")
    axs[1].imshow(ri)
    axs[1].axis("off")
    axs[1].set_title("Sampled")
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)

def extract_data(script_path: Path):
    txt = script_path.read_text()
    m = re.search(r"data\\s*=\\s*np\\.array\\((\\[[\\s\\S]*?\\])\\)", txt)
    if not m:
        return None
    return np.asarray(ast.literal_eval(m.group(1)), dtype=float).ravel()

def dist_stats(a, b):
    a = np.asarray(a, dtype=float).ravel()
    b = np.asarray(b, dtype=float).ravel()
    if len(a) < 2 or len(b) < 2:
        return {"bins": 0, "jsd": float("nan"), "emd": float("nan")}
    mix = np.concatenate([a, b])
    edges = np.histogram_bin_edges(mix, bins="fd")
    if len(edges) < 3:
        edges = np.histogram_bin_edges(mix, bins=20)
    ha, _ = np.histogram(a, bins=edges, density=False)
    hb, _ = np.histogram(b, bins=edges, density=False)
    eps = 1e-12
    pa = (ha + eps) / np.sum(ha + eps)
    pb = (hb + eps) / np.sum(hb + eps)
    mm = 0.5 * (pa + pb)
    jsd = np.sqrt(0.5 * (np.sum(pa * np.log(pa / mm)) + np.sum(pb * np.log(pb / mm))))
    cdfa = np.cumsum(pa)
    cdfb = np.cumsum(pb)
    emd = np.sum(np.abs(cdfa - cdfb) * np.diff(edges))
    return {"bins": int(len(edges) - 1), "jsd": float(jsd), "emd": float(emd)}

ini = parse_scripts(INI_DIR)
sam = parse_scripts(SAM_DIR)

if len(ini) == 0 and len(sam) == 0:
    print(f"No histogram scripts found under {ROOT}.")
    print("Run icats.init first (with histogram mode enabled), then re-run this script.")
    raise SystemExit(0)

sys_map = {
    "j": "orb_ij",
    "l": "orb_il",
    "b": "orb_sb",
    "phi": "orb_sphi",
    "vel": "vel_ivel",
}
sys_fallback = {
    "j": "orb_sj",
    "l": "orb_sl",
}
sys_default = {"j", "l", "vel", "b", "phi"}

pairs = []
used_ini = set()
used_sam = set()
for (scope, metric), ipath in ini.items():
    if scope == "sys":
        if metric not in sys_default:
            continue
        skey = (scope, sys_map.get(metric, metric))
        if skey not in sam and metric in sys_fallback:
            skey = (scope, sys_fallback[metric])
    elif scope in ("m0", "m1"):
        if metric.startswith("vi_mode"):
            if not args.include_vib:
                continue
        elif metric == "j":
            if not args.include_rot:
                continue
        else:
            continue
        skey = (scope, metric)
    else:
        continue
    if skey in sam:
        pairs.append(((scope, metric), skey, ipath, sam[skey]))
        used_ini.add((scope, metric))
        used_sam.add(skey)

manifest = []
for (iscope, imetric), (_sscope, smetric), ipath, spath in pairs:
    key = f"{iscope}_{imetric}"
    ipng = render_hist(ipath, OUT_INI / f"{key}_ini")
    spng = render_hist(spath, OUT_SAM / f"{key}_sam")
    cpng = OUT_PAIR / f"{key}__ini_vs_sam.png"
    desc = f"Compare {iscope}:{imetric} (init) vs {iscope}:{smetric} (sampled)"
    side_by_side(ipng, spng, desc, cpng)
    ida = extract_data(ipath)
    sda = extract_data(spath)
    st = dist_stats(ida, sda) if ida is not None and sda is not None else {"bins": 0, "jsd": float("nan"), "emd": float("nan")}
    manifest.append((key, desc, str(st["bins"]), f"{st['jsd']:.6g}", f"{st['emd']:.6g}", str(cpng)))

if args.include_unmatched:
    for (scope, metric), p in ini.items():
        if (scope, metric) in used_ini:
            continue
        render_hist(p, OUT_UNM_INI / f"{scope}_{metric}_ini")
    for (scope, metric), p in sam.items():
        if (scope, metric) in used_sam:
            continue
        render_hist(p, OUT_UNM_SAM / f"{scope}_{metric}_sam")

mf = OUT / "pair_manifest.tsv"
hdr = "key\\tdescription\\tbins\\tjsd\\temd\\tcomparison_png\\n"
mf.write_text(hdr + "\\n".join("\\t".join(x) for x in manifest) + "\\n")
print(f"Paired comparisons: {len(manifest)}")
print(f"Manifest: {mf}")
print(f"Pair plots: {OUT_PAIR}")
if args.include_unmatched:
    print(f"Unmatched init: {OUT_UNM_INI}")
    print(f"Unmatched sampled: {OUT_UNM_SAM}")
else:
    print("Unmatched histograms skipped (use --include-unmatched to render them).")
""",
        executable=True,
    )

    _write_file(
        hist_dir / "plot_compare_pairs.sh",
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "SCRIPT_DIR=\"$(cd \"$(dirname \"${BASH_SOURCE[0]}\")\" && pwd)\"\n"
        "export HIST_ROOT=\"${HIST_ROOT:-$SCRIPT_DIR}\"\n"
        "export MPLCONFIGDIR=\"${MPLCONFIGDIR:-/tmp/mpl_cache_${USER:-user}}\"\n"
        "mkdir -p \"$MPLCONFIGDIR\"\n"
        "python \"$SCRIPT_DIR/plot_compare_pairs.py\" \"$@\"\n",
        executable=True,
    )

    _write_file(
        hist_dir / "HISTOGRAMS_README.txt",
        "Histogram helpers\n"
        "=================\n\n"
        "Run from tutorial root:\n"
        "- ./rd_tutorial_input/histograms/plot_initial.sh\n"
        "- ./rd_tutorial_input/histograms/plot_sampled.sh\n"
        "- ./rd_tutorial_input/histograms/plot_all.sh\n\n"
        "- ./rd_tutorial_input/histograms/plot_compare_pairs.sh\n"
        "- ./rd_tutorial_input/histograms/plot_compare_pairs.sh --include-rot --include-vib\n"
        "- ./rd_tutorial_input/histograms/plot_compare_pairs.sh --include-unmatched\n\n"
        "What to expect:\n"
        "- Initial histograms (hist_ini_*) show the distributions used to draw initial conditions.\n"
        "- Sampled histograms (hist_sam_*) summarize the generated sample set.\n"
        "- Default paired comparisons are intermolecular/system only: j, l, vel, and if available b, phi.\n"
        "- Add --include-rot and/or --include-vib for molecular rotational/vibrational comparisons.\n"
        "- For shared quantities, initial and sampled shapes\n"
        "  should be qualitatively similar when sampling is healthy.\n"
        "- `plot_compare_pairs.sh` is fast by default (paired comparisons only).\n"
        "- Add `--include-unmatched` when you also want every remaining histogram rendered.\n"
        "- pair_manifest.tsv reports bins, JSD, and EMD for each paired comparison.\n"
        "- Histogram data and plots are stored under rd_tutorial_input/histograms/ by default.\n"
        "- Differences are expected from finite-sample noise and downstream conditioning.\n",
    )


def _write_polarization_histogram_helper(out_dir: Path, run_tag: str = "tutorial_input") -> None:
    hist_dir = out_dir / f"rd_{run_tag}" / "histograms"
    hist_dir.mkdir(parents=True, exist_ok=True)
    _write_file(
        hist_dir / "plot_polarization_check.py",
        """#!/usr/bin/env python3
from __future__ import annotations

import ast
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def trapezoid(y, x):
    implementation = getattr(np, "trapezoid", None)
    if implementation is None:
        implementation = np.trapz
    return implementation(y, x)

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "plots" / "polarization"
OUT.mkdir(parents=True, exist_ok=True)

# These values match icats/examples/no_polarized_dat.txt.
A = 0.75
field_theta = np.pi / 4.0
field_phi = 0.0


def _read_hist_data(candidates):
    for path in candidates:
        if path.exists():
            txt = path.read_text()
            m = re.search(r"data\\s*=\\s*np\\.array\\((\\[[\\s\\S]*?\\])\\)", txt)
            if not m:
                raise SystemExit(f"Could not read embedded data from {path}")
            return np.asarray(ast.literal_eval(m.group(1)), dtype=float), path
    raise SystemExit("Missing orientation histogram scripts. Run icats.init tutorial_input.txt first.")


alpha, alpha_path = _read_hist_data(
    [
        ROOT / "sampled" / "molecule_m1" / "hist_sam_m1_salpha.py",
        ROOT / "sampled" / "molecule_m1" / "hist_sam_m1_alpha.py",
        ROOT / "sampled" / "molecule_m1" / "hist_sam_m1_ori_alpha.py",
        ROOT / "initial" / "molecule_m1" / "hist_ini_m1_alpha.py",
        ROOT / "initial" / "molecule_m1" / "hist_ini_m1_ori_alpha.py",
    ]
)
beta, beta_path = _read_hist_data(
    [
        ROOT / "sampled" / "molecule_m1" / "hist_sam_m1_sbeta.py",
        ROOT / "sampled" / "molecule_m1" / "hist_sam_m1_beta.py",
        ROOT / "sampled" / "molecule_m1" / "hist_sam_m1_ori_beta.py",
        ROOT / "initial" / "molecule_m1" / "hist_ini_m1_beta.py",
        ROOT / "initial" / "molecule_m1" / "hist_ini_m1_ori_beta.py",
    ]
)

u = np.cos(beta)

fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))

alpha_grid = np.linspace(-np.pi, np.pi, 500)
alpha_pdf = 1.0 + A * (np.pi / 4.0) * np.sin(field_theta) * np.cos(alpha_grid - field_phi)
axs[0].hist(alpha, bins=40, density=True, alpha=0.72, edgecolor="k", label="sampled")
axs[0].plot(alpha_grid, alpha_pdf / trapezoid(alpha_pdf, alpha_grid), "r-", lw=2, label="expected")
axs[0].set_xlabel("alpha / rad")
axs[0].set_ylabel("density")
axs[0].set_title("Azimuthal bias from tilted field")
axs[0].legend()

u_grid = np.linspace(-1.0, 1.0, 500)
u_pdf = 1.0 + A * np.cos(field_theta) * u_grid
axs[1].hist(u, bins=40, density=True, alpha=0.72, edgecolor="k", label="sampled")
axs[1].plot(u_grid, 0.5 * u_pdf, "r-", lw=2, label="expected")
axs[1].set_xlabel("cos(beta)")
axs[1].set_ylabel("density")
axs[1].set_title("Polar bias of body-z axis")
axs[1].legend()

fig.suptitle("NO orientation PDF in the ICATS scattering frame")
fig.tight_layout()
out = OUT / "polarized_orientation_check.png"
fig.savefig(out, dpi=160)
print(f"Read alpha from {alpha_path}")
print(f"Read beta from {beta_path}")
print(f"Wrote {out}")
""",
        executable=True,
    )


def _dat_references(dat_path: Path) -> list[str]:
    refs = []
    for raw in dat_path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        parts = line.split()
        if len(parts) >= 3 and parts[0].lower() in ("xyz", "hess", "w"):
            refs.append(parts[2])
    return refs


def _tutorial_input_text(name: str, cfg: dict, hist_samples: int | None = None) -> str:
    hist_enabled = (hist_samples is not None and hist_samples > 0) or bool(cfg.get("hist_by_default", False))
    nsamp_value = int(hist_samples) if hist_samples is not None and hist_samples > 0 else int(cfg["nsamp"])
    # Histogram diagnostics should not heavily clip orbital-L via b cutoff.
    maxb_value = float(cfg["maxb"])
    is_wl = any(line.strip().lower() == "wang = true" for line in cfg.get("extra", []))
    ang_key = "maxj" if is_wl else "maxl"
    # Histogram-heavy tutorials are typically distribution diagnostics; avoid writing
    # per-sample xyz/vel files by default in this mode.
    printout_line = "printout = 0 0 0 0" if hist_enabled else f"printout = {cfg.get('printout', '0 0 1 0')}"
    lines = [
        f"# Tutorial: {name}",
        f"# {cfg['desc']}",
        "",
        "# Molecules (required)",
        f"mol = 0 {cfg['mol0']}",
        f"mol = 1 {cfg['mol1']}",
        "",
        "# Core sampling controls",
        f"Nsamp = {nsamp_value}",
        f"workers = {int(cfg.get('workers', 1))}",
        "seed = 400",
        "continue = False",
        "run-mode = fresh",
        "run-tag = tutorial_input",
        "seed-mode = fixed",
        "",
        "# Thermal/beam settings",
        f"Tvib = {float(cfg.get('tvib', 500.0))}",
        f"Trot = {float(cfg.get('trot', 500.0))}",
        f"Rz = {cfg.get('rz', 15)}",
        "",
        f"# Orbital controls ({ang_key} and maxb)",
        f"{ang_key} = {cfg['maxl']}",
        f"maxb = {maxb_value:g}",
        "",
        "# Output",
        "fileout = out",
        "dirout = outputs",
        printout_line,
        f"plothist = {'True' if hist_enabled else 'False'}",
        f"hist_initial = {'True' if hist_enabled else 'False'}",
        f"hist_sampled = {'True' if hist_enabled else 'False'}",
        "progress = normal",
        "output-format = xyzvel",
        "units-out = ang-fs",
        "output-frame = internal",
        "",
        "# Tutorial-specific options",
    ]
    if hist_enabled:
        lines += [
            f"plotinit = {nsamp_value}",
            "# plotinit draws MC samples for distribution histograms",
            "# Nsamp is matched to plotinit in histogram mode for fair init-vs-sampled comparisons",
            "# maxb is kept fixed in histogram mode to match the chosen angular-momentum cap",
            "# hist_initial/hist_sampled control initial and sampled histogram generation",
            "# printout is set to 0 0 0 0 in histogram mode to suppress xyz/vel outputs",
        ]
    lines.extend(cfg["extra"])
    lines.extend([
        "",
        "# Notes:",
    ])
    lines.extend([f"# - {n}" for n in cfg["notes"]])
    return "\n".join(lines) + "\n"


def _generate_tutorial(
    name: str,
    out_dir: Path,
    nsamp_override: int | None = None,
    ntraj_override: int | None = None,
    hist_samples: int | None = None,
) -> None:
    cfg = TUTORIALS[name]
    nsamp = int(nsamp_override) if nsamp_override is not None else int(cfg["nsamp"])
    ntraj = int(ntraj_override) if ntraj_override is not None else int(cfg["ntraj"])
    hist_enabled = (hist_samples is not None and hist_samples > 0) or bool(cfg.get("hist_by_default", False))
    is_wl = any(line.strip().lower() == "wang = true" for line in cfg.get("extra", []))
    is_paper = name == "paper_nh3_h2o_100k"
    ex_dir = _examples_dir()
    out_dir.mkdir(parents=True, exist_ok=True)

    # Copy molecule inputs and their referenced files.
    dats = sorted({cfg["mol0"], cfg["mol1"]})
    for dat_name in dats:
        dat_src = ex_dir / dat_name
        if not dat_src.exists():
            raise SystemExit(f"Missing example data file: {dat_src}")
        shutil.copy2(dat_src, out_dir / dat_name)
        for ref in _dat_references(dat_src):
            ref_src = ex_dir / ref
            if ref_src.exists():
                shutil.copy2(ref_src, out_dir / ref)
            else:
                raise SystemExit(f"Missing referenced file '{ref}' for {dat_name}")
    for extra_name in cfg.get("extra_files", []):
        extra_src = ex_dir / extra_name
        if not extra_src.exists():
            raise SystemExit(f"Missing example companion file: {extra_src}")
        shutil.copy2(extra_src, out_dir / extra_name)

    cfg_for_file = dict(cfg)
    cfg_for_file["nsamp"] = nsamp
    cfg_for_file["ntraj"] = ntraj
    _write_file(out_dir / "tutorial_input.txt", _tutorial_input_text(name, cfg_for_file, hist_samples=hist_samples))
    dyn_runner_src = _dyn_runner_path()
    if not dyn_runner_src.exists():
        raise SystemExit(f"Missing tutorial dynamics runner: {dyn_runner_src}")
    shutil.copy2(dyn_runner_src, out_dir / "run_tutorial_dyn.py")
    (out_dir / "run_tutorial_dyn.py").chmod(0o755)

    wl_umbrella_note = ""
    pre_src = _precomputed_dir() / name / "rd_tutorial_input"
    if not pre_src.exists():
        pre_src = _precomputed_dir() / name / "run_data"
    if pre_src.exists():
        if pre_src.name == "rd_tutorial_input":
            shutil.copytree(pre_src, out_dir / "rd_tutorial_input", dirs_exist_ok=True)
        else:
            legacy_input = pre_src / "tutorial_input"
            if legacy_input.exists():
                shutil.copytree(legacy_input, out_dir / "rd_tutorial_input", dirs_exist_ok=True)
            else:
                shutil.copytree(pre_src, out_dir / "rd_tutorial_input", dirs_exist_ok=True)
        run_input_dir = out_dir / "rd_tutorial_input"
        wl_dir = run_input_dir / "histograms" / "wl"
        wl_dir.mkdir(parents=True, exist_ok=True)
        for fn in ("wl_td_plot.py", "wl_wl_plot.py"):
            old = run_input_dir / fn
            if old.exists():
                shutil.move(str(old), str(wl_dir / fn))
        wl_umbrella_note = (
            "Precomputed WL umbrella copied:\n"
            "- rd_tutorial_input/wang.pkl\n"
            "- rd_tutorial_input/histograms/wl/wl_td_plot.py\n"
            "- rd_tutorial_input/histograms/wl/wl_wl_plot.py\n"
            "Visualize with:\n"
            "- python rd_tutorial_input/histograms/wl/wl_td_plot.py\n"
            "- python rd_tutorial_input/histograms/wl/wl_wl_plot.py\n\n"
            "If input settings change, move or rename rd_tutorial_input/wang.pkl before rerunning.\n"
            "icats.init will refuse incompatible WL umbrellas instead of overwriting them.\n\n"
        )

    _write_file(
        out_dir / "run_cheap_dynamics.sh",
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "cd \"$(dirname \"${BASH_SOURCE[0]}\")\"\n"
        "# Cluster/sandbox compatibility: keep caches on writable paths.\n"
        "export NUMBA_CACHE_DIR=\"${NUMBA_CACHE_DIR:-/tmp/numba_cache_${USER:-user}}\"\n"
        "export MPLCONFIGDIR=\"${MPLCONFIGDIR:-/tmp/mpl_cache_${USER:-user}}\"\n"
        "mkdir -p \"$NUMBA_CACHE_DIR\" \"$MPLCONFIGDIR\"\n"
        "CHARGE=\"${CHARGE:-0}\"\n"
        "SPIN=\"${SPIN:-0}\"\n"
        f"./run_tutorial_dyn.py --outdir rd_tutorial_input/outputs --prefix out --ntraj {ntraj} --steps {cfg['steps']} --charge \"$CHARGE\" --spin \"$SPIN\"\n",
        executable=True,
    )

    _write_file(
        out_dir / "run_analysis.sh",
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "cd \"$(dirname \"${BASH_SOURCE[0]}\")\"\n"
        "# Cluster/sandbox compatibility: keep caches on writable paths.\n"
        "export NUMBA_CACHE_DIR=\"${NUMBA_CACHE_DIR:-/tmp/numba_cache_${USER:-user}}\"\n"
        "export MPLCONFIGDIR=\"${MPLCONFIGDIR:-/tmp/mpl_cache_${USER:-user}}\"\n"
        "mkdir -p \"$NUMBA_CACHE_DIR\" \"$MPLCONFIGDIR\"\n"
        f"icats.analyse tutorial_input.txt --dir rd_tutorial_input/outputs --prefix out --ntraj {ntraj}\n",
        executable=True,
    )

    _write_file(
        out_dir / "run_continue.sh",
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "cd \"$(dirname \"${BASH_SOURCE[0]}\")\"\n"
        "# Cluster/sandbox compatibility: keep caches on writable paths.\n"
        "export NUMBA_CACHE_DIR=\"${NUMBA_CACHE_DIR:-/tmp/numba_cache_${USER:-user}}\"\n"
        "export MPLCONFIGDIR=\"${MPLCONFIGDIR:-/tmp/mpl_cache_${USER:-user}}\"\n"
        "mkdir -p \"$NUMBA_CACHE_DIR\" \"$MPLCONFIGDIR\"\n"
        "ADD_NSAMP=\"${1:-10}\"\n"
        "export ADD_NSAMP\n"
        "cp tutorial_input.txt tutorial_input.txt.bak\n"
        "trap 'mv -f tutorial_input.txt.bak tutorial_input.txt' EXIT\n"
        "python3 - <<'PY'\n"
        "from pathlib import Path\n"
        "import os\n"
        "p = Path('tutorial_input.txt')\n"
        "ns = os.environ.get('ADD_NSAMP', '10')\n"
        "lines = p.read_text().splitlines()\n"
        "out = []\n"
        "seen_nsamp = False\n"
        "seen_mode = False\n"
        "for ln in lines:\n"
        "    s = ln.strip().lower()\n"
        "    if s.startswith('nsamp') and '=' in ln:\n"
        "        out.append(f'Nsamp = {ns}')\n"
        "        seen_nsamp = True\n"
        "    elif s.startswith('run-mode') and '=' in ln:\n"
        "        out.append('run-mode = continue')\n"
        "        seen_mode = True\n"
        "    else:\n"
        "        out.append(ln)\n"
        "if not seen_nsamp:\n"
        "    out.append(f'Nsamp = {ns}')\n"
        "if not seen_mode:\n"
        "    out.append('run-mode = continue')\n"
        "p.write_text('\\n'.join(out) + '\\n')\n"
        "PY\n"
        "ADD_NSAMP=\"$ADD_NSAMP\" icats.init tutorial_input.txt\n"
        "mv -f tutorial_input.txt.bak tutorial_input.txt\n"
        "trap - EXIT\n",
        executable=True,
    )

    _write_file(
        out_dir / "setup_conda_env.sh",
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n\n"
        "ENV_NAME=\"${1:-icats_clean}\"\n\n"
        "echo \"Creating conda env: $ENV_NAME\"\n"
        "conda create -y -n \"$ENV_NAME\" -c conda-forge \\\n"
        "  python=3.11 numpy=1.26 scipy=1.11 h5py=3.10 pyscf=2.4 \\\n"
        "  joblib tqdm matplotlib numba ffmpeg\n\n"
        "echo \"Activating env and installing pip extras...\"\n"
        "eval \"$(conda shell.bash hook)\"\n"
        "conda activate \"$ENV_NAME\"\n"
        "pip install --upgrade pip\n"
        "pip install pyscf-semiempirical spherical quaternionic\n\n"
        "echo \"Verifying core imports...\"\n"
        "python - <<'PYCHECK'\n"
        "import numpy, h5py, pyscf\n"
        "print(\"numpy\", numpy.__version__)\n"
        "print(\"h5py\", h5py.__version__)\n"
        "print(\"pyscf\", pyscf.__version__)\n"
        "PYCHECK\n\n"
        "echo \"\"\n"
        "echo \"Done. Activate with:\"\n"
        "echo \"  conda activate $ENV_NAME\"\n",
        executable=True,
    )

    notes_text = "".join(f"- {note}\n" for note in cfg["notes"])
    inputs_text = (
        "Input files copied into this tutorial:\n"
        f"- tutorial_input.txt\n"
        f"- {cfg['mol0']}\n"
        f"- {cfg['mol1']}\n\n"
        "Frame convention:\n"
        "- tutorial_input.txt includes output-frame = internal by default.\n"
        "- Keep this for the historical ICATS/tutorial convention.\n"
        "- Change it to output-frame = incoming-k-plus-z before icats.init if an external scattering/QM code expects incoming k along space-fixed +Z.\n"
        "- Changing the frame changes vector components, SF/BF Euler angles, and exported xyz/vel coordinates, so regenerate samples after changing it.\n\n"
    )

    if is_paper:
        run_order = (
            "Run sequence:\n"
            "1) Read this file and inspect tutorial_input.txt.\n"
            "2) Submit to a configured SLURM compute node:\n"
            "   sbatch run_paper_slurm.sh\n"
            "   Alternatively, on a suitable interactive compute node run:\n"
            "   icats.init tutorial_input.txt\n"
            "3) For an interactive run, export the paper data after sampling:\n"
            "   python export_paper_histogram_data.py\n"
            "4) Create the paper-validation plots:\n"
            "   python plot_paper_histograms.py\n"
            "   The SLURM wrapper performs all three commands automatically.\n\n"
        )
    else:
        run_order = (
            "Run sequence:\n"
            "1) Read this file and inspect tutorial_input.txt.\n"
            "2) Run initial-condition generation:\n"
            "   icats.init tutorial_input.txt\n"
        )
    if is_wl and not is_paper:
        run_order += (
            "3) For Wang-Landau tutorials, the first run may mainly build rd_tutorial_input/wang.pkl.\n"
            "   If sampled histogram scripts are missing after the first run, run the same command again:\n"
            "   icats.init tutorial_input.txt\n"
            "   The second run should reuse the compatible wang.pkl and generate sampled histograms.\n"
        )
        dynamics_step = 4
        analysis_step = 5
    elif not is_paper:
        dynamics_step = 3
        analysis_step = 4
    if not is_paper:
        run_order += (
            f"{dynamics_step}) Optional cheap dynamics demonstration:\n"
            "   ./run_cheap_dynamics.sh\n"
            f"{analysis_step}) Optional trajectory analysis after cheap dynamics:\n"
            "   ./run_analysis.sh\n\n"
        )

    expected_text = (
        "Main outputs to expect:\n"
        "- tutorial_input.txt.logfile\n"
        "- rd_tutorial_input/\n"
    )
    if is_wl:
        expected_text += (
            "- rd_tutorial_input/wang.pkl\n"
            "- rd_tutorial_input/histograms/wl/wl_td_plot.py\n"
            "- rd_tutorial_input/histograms/wl/wl_wl_plot.py\n"
        )
    if hist_enabled:
        expected_text += (
            "- rd_tutorial_input/histograms/initial/\n"
            "- rd_tutorial_input/histograms/sampled/\n"
            "- rd_tutorial_input/histograms/plots/\n"
        )
    if is_paper:
        expected_text += (
            "- rd_tutorial_input/dat_tutorial_input.txt.pkl\n"
            "- paper_histogram_data/figure6_*.csv\n"
            "- paper_histogram_data/figure7_*.csv\n"
            "- paper_histogram_data/figure9_*.csv\n"
            "- paper_histogram_data/metadata.json\n\n"
            "- paper_histogram_plots/figure6_validation.png and .pdf\n"
            "- paper_histogram_plots/figure7_validation.png and .pdf\n"
            "- paper_histogram_plots/figure9_validation.png and .pdf\n\n"
        )
    else:
        expected_text += (
            "- rd_tutorial_input/outputs/out_*.xyz and out_*.vel when printout requests trajectory files\n"
            "- rd_tutorial_input/outputs/dynamics*.analinfo after ./run_analysis.sh\n\n"
        )

    after_init_note = ""
    if hist_enabled:
        after_init_note += (
            "After icats.init: histogram checks\n"
            "- Histogram scripts live under rd_tutorial_input/histograms/.\n"
            "- The broad helpers are available, but may plot many files:\n"
            "  ./rd_tutorial_input/histograms/plot_initial.sh\n"
            "  ./rd_tutorial_input/histograms/plot_sampled.sh\n"
            "  ./rd_tutorial_input/histograms/plot_compare_pairs.sh\n\n"
            "- For a quick intermolecular check, plot only system L, J, and Jab:\n"
            "  ./rd_tutorial_input/histograms/plot_orbital_jljab.sh\n  ./rd_tutorial_input/histograms/plot_orbital_correlation.sh\n\n"
            "- Open the resulting PNGs in:\n"
            "  rd_tutorial_input/histograms/plots/sampled/\n\n"
        )
    if is_wl:
        after_init_note += (
            "After icats.init: Wang-Landau checks\n"
            "- Confirm the stored umbrella exists:\n"
            "  rd_tutorial_input/wang.pkl\n"
            "- List WL plotting helpers:\n"
            "  find rd_tutorial_input/histograms/wl -maxdepth 1 -type f | sort\n"
            "- Common WL plots, when present:\n"
            "  cd rd_tutorial_input/histograms/wl\n"
            "  python wl_td_plot.py\n"
            "  python wl_wl_plot.py\n"
            "  cd ../../..\n\n"
        )

    if is_paper:
        general_advice = (
            "General Advice\n"
            "==============\n\n"
            "Environment setup:\n"
            "1) Install and activate the frozen ICATS release used for the paper.\n"
            "2) Confirm that `icats.init --help` resolves inside that environment.\n\n"
            "Cluster note:\n"
            "- Do not run the 100,000-sample calculation on a login node.\n"
            "- Inspect the CPU, memory, wall-time, account, and partition directives in run_paper_slurm.sh before submission.\n"
            "- The wrapper places numba and matplotlib caches under writable /tmp paths.\n\n"
            "Rebuilding the umbrella:\n"
            "- A compatible rd_tutorial_input/wang.pkl is reused on rerun.\n"
            "- Move or rename it before rerunning only when the WL settings have changed.\n"
            "- ICATS refuses metadata-incompatible umbrellas rather than silently overwriting them.\n\n"
            "Interpreting the validation plots:\n"
            "- Figure 6: sampled Q histograms should follow the smooth leading-Wigner curves; the black curves are the oscillatory harmonic-oscillator densities.\n"
            "- Figure 7: the rotational J distributions should follow their thermal reference curves. The NH3 K distribution has an analytic symmetric-top reference; the H2O panel shows the sampled projection along the Wang-state-selected body axis.\n"
            "- Figure 9: L and J should follow the triangular reference over most of the requested range, the relative-speed histogram should follow the transformed pair of weighted molecular-beam distributions, and the WL correction should become comparatively flat once L dominates Jab.\n"
            "- Retain wang.pkl and metadata.json. A visually poor or strongly structured result should be diagnosed rather than hidden by replotting.\n\n"
        )
    else:
        general_advice = (
            "General Advice\n"
            "==============\n\n"
            "Environment setup:\n"
            "1) ./setup_conda_env.sh icats_clean\n"
            "2) conda activate icats_clean\n\n"
            "Cluster/sandbox note:\n"
            "- ICATS entry points and tutorial scripts set writable /tmp cache paths by default.\n"
            "- This avoids numba/matplotlib cache errors on restricted/shared filesystems.\n\n"
            "Spin/charge note for dynamics:\n"
            "- Defaults are CHARGE=0 and SPIN=0.\n"
            "- For open-shell systems set env vars, e.g.:\n"
            "- CHARGE=0 SPIN=1 ./run_cheap_dynamics.sh\n\n"
            "Continue sampling (adds more initial conditions):\n"
            "- ./run_continue.sh 10\n\n"
        )

    readme_text = (
        f"Tutorial: {name}\n"
        f"{cfg['desc']}\n\n"
        "This Tutorial\n"
        "=============\n\n"
        "Purpose:\n"
        f"{notes_text}\n"
        f"{inputs_text}"
        f"{run_order}"
        f"{wl_umbrella_note}"
        f"{after_init_note}"
        f"{expected_text}"
        f"{general_advice}"
    )
    if not hist_enabled and not is_paper:
        readme_text += (
            "Optional histogram-heavy setup:\n"
            "- Regenerate this tutorial with --histograms [--hist-samples 10000].\n"
            "- Then run icats.init tutorial_input.txt to populate rd_tutorial_input/histograms/.\n\n"
        )
    _write_file(out_dir / "tutorial_README.txt", readme_text)
    _write_histogram_helpers(
        out_dir,
        run_tag="tutorial_input",
        include_paper=name == "paper_nh3_h2o_100k",
    )
    if name == "polarized_orientation_he_no":
        _write_polarization_histogram_helper(out_dir, run_tag="tutorial_input")
    if is_paper:
        for helper in (
            "run_analysis.sh",
            "run_cheap_dynamics.sh",
            "run_continue.sh",
            "run_tutorial_dyn.py",
            "setup_conda_env.sh",
        ):
            (out_dir / helper).unlink(missing_ok=True)


def _default_tutorial_dir(name: str) -> str:
    if name == "quickstart":
        return "tutorial_quickstart"
    return f"tutorial_{name}"


def main() -> int:
    ap = argparse.ArgumentParser(prog="icats", description="ICATS helper CLI")
    ap.add_argument("--list-tutorials", action="store_true", help="List available tutorials")
    ap.add_argument(
        "--tutorial",
        nargs="?",
        const="quickstart",
        help="Generate tutorial directory (name or index; default: quickstart)",
    )
    ap.add_argument("--tutorial-dir", default=None, help="Output directory for generated tutorial")
    ap.add_argument("--nsamp", type=int, default=None, help="Override Nsamp in generated tutorial input")
    ap.add_argument("--ntraj", type=int, default=None, help="Override dynamics/analysis trajectory count in helper scripts")
    ap.add_argument("--histograms", action="store_true", help="Enable histogram generation in tutorial_input.txt")
    ap.add_argument("--hist-samples", type=int, default=10000, help="plotinit count used with --histograms (default: 10000)")
    ap.add_argument("--setup-only", action="store_true", help="Generate files only, skip icats.init run")
    ap.add_argument("--run-dynamics", action="store_true", help="Also run cheap dynamics after generation")
    args = ap.parse_args()

    if args.list_tutorials:
        print("Available tutorials:")
        for i, name in enumerate(TUTORIAL_ORDER):
            print(f"{i}: {name} - {TUTORIALS[name]['desc']}")
        return 0

    if not args.tutorial:
        ap.print_help()
        return 0

    if args.setup_only and args.run_dynamics:
        raise SystemExit("Invalid options: --setup-only cannot be combined with --run-dynamics")
    if args.nsamp is not None and args.nsamp <= 0:
        raise SystemExit("--nsamp must be a positive integer")
    if args.ntraj is not None and args.ntraj <= 0:
        raise SystemExit("--ntraj must be a positive integer")
    if args.hist_samples <= 0:
        raise SystemExit("--hist-samples must be a positive integer")

    tut = args.tutorial
    if tut.isdigit():
        idx = int(tut)
        if idx < 0 or idx >= len(TUTORIAL_ORDER):
            raise SystemExit(f"Unknown tutorial index '{tut}'. Use --list-tutorials")
        tut = TUTORIAL_ORDER[idx]

    if tut not in TUTORIALS:
        raise SystemExit(f"Unknown tutorial '{args.tutorial}'. Use --list-tutorials")

    out_name = args.tutorial_dir if args.tutorial_dir else _default_tutorial_dir(tut)
    out_dir = Path(out_name).resolve()

    hist_samples = args.hist_samples if args.histograms else None
    _generate_tutorial(
        tut,
        out_dir,
        nsamp_override=args.nsamp,
        ntraj_override=args.ntraj,
        hist_samples=hist_samples,
    )
    print(f"[1/2] Tutorial '{tut}' generated in: {out_dir}")

    if not args.setup_only:
        subprocess.run([sys.executable, "-m", "icats.cli_init", "tutorial_input.txt"], cwd=str(out_dir), check=True)
        print("[2/2] Generated initial conditions via icats.init")
    else:
        print("[2/2] Setup-only mode: skipped icats.init run")

    if args.run_dynamics:
        subprocess.run([str(out_dir / "run_cheap_dynamics.sh")], cwd=str(out_dir), check=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
