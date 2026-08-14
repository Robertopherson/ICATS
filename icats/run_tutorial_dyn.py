#!/usr/bin/env python3
import argparse
import glob
import os
import warnings
from pathlib import Path

import numpy as np

gto = None
md = None
se_mindo3 = None


# ---- unit conversions ----
FMT2AU = 41.3413745
AU2FMT = 1.0 / FMT2AU
ANG2AU = 1.88973
AU2ANG = 1.0 / ANG2AU
AU2EV = 1.0 / 0.0367494


def read_velocities_xyz_like(path):
    """Read per-atom velocities from XYZ-like text: Elem vx vy vz (Ang/fs)."""
    lines = Path(path).read_text().splitlines()
    data = []
    for ln in lines[2:]:
        parts = ln.split()
        if len(parts) >= 4:
            data.append([float(parts[-3]), float(parts[-2]), float(parts[-1])])
    if not data:
        return None
    v = np.array(data, dtype=float) * (ANG2AU / FMT2AU)  # Ang/fs -> bohr/a.u.
    return v


def build_mindo3_scanner(mol):
    if se_mindo3 is None:
        raise RuntimeError("Missing pyscf-semiempirical. Install: pip install pyscf-semiempirical")
    if hasattr(se_mindo3, "RMINDO3"):
        mf = se_mindo3.RMINDO3(mol) if mol.spin == 0 else se_mindo3.UMINDO3(mol)
    elif hasattr(se_mindo3, "MINDO3"):
        mf = se_mindo3.MINDO3(mol)
    else:
        raise RuntimeError("Could not find RMINDO3/UMINDO3/MINDO3 in pyscf.semiempirical.")
    mf.conv_tol = 1e-6
    mf.max_cycle = 50
    mf.verbose = 0
    mf.kernel()
    return mf.nuc_grad_method().as_scanner()


def discover_pairs(outdir, prefix):
    xyzs = glob.glob(os.path.join(outdir, f"{prefix}_*.xyz"))

    def _idx(path):
        base = os.path.basename(path)
        stem = base[:-4]  # strip .xyz
        try:
            return int(stem.rsplit("_", 1)[1])
        except Exception:
            return 10**9

    xyzs = sorted(xyzs, key=_idx)
    pairs = []
    for xf in xyzs:
        vf = xf[:-4] + ".vel"
        if os.path.exists(vf):
            pairs.append((xf, vf))
    return pairs


def run_one(geom_file, vel_file, dt_au, steps, charge, spin, ref_en):
    mol = gto.Mole()
    mol.atom = geom_file
    mol.basis = "sto-3g"
    mol.charge = charge
    mol.spin = spin
    mol.symmetry = False
    mol.verbose = 0
    mol.build()

    scanner = build_mindo3_scanner(mol)
    init_veloc = read_velocities_xyz_like(vel_file)

    velocs = {}

    def _cb(local):
        t = local["self"].time
        velocs[t] = local["self"].veloc.copy()

    stem = os.path.splitext(geom_file)[0]
    md_data = stem + ".md.data"
    md_xyz = stem + ".md.xyz"
    integ = md.NVE(
        scanner,
        dt=dt_au,
        steps=steps,
        veloc=init_veloc,
        callback=_cb,
        data_output=md_data,
        trajectory_output=md_xyz,
    )
    with open(os.devnull, "w") as quiet:
        integ.stdout = quiet
        integ.run(verbose=0)
    integ.data_output.close()
    integ.trajectory_output.close()

    # Convert md.data to .dat in eV, optionally subtracting reference energy
    lines = Path(md_data).read_text().splitlines()
    out = [lines[0] + "\n"] if lines else []
    total_energies = []
    for ln in lines[1:]:
        toks = ln.split()
        if len(toks) < 2:
            continue
        step = toks[0]
        vals = np.array([float(x) for x in toks[1:]], dtype=float)
        if len(vals) >= 3:
            vals[0] -= ref_en
            vals[2] -= ref_en
        vals *= AU2EV
        if len(vals) >= 3:
            total_energies.append(vals[2])
        out.append(step + "".join(f"{v:20.9f}" for v in vals) + "\n")
    Path(stem + ".dat").write_text("".join(out))

    # Dump per-step velocities in Ang/fs
    na = len(mol._atom)
    el = [atm[0] for atm in mol._atom]
    vout = []
    for t in sorted(velocs.keys()):
        vout.append(f"{na}\n")
        vout.append(f"Time = {t}\n")
        vv = velocs[t] * (AU2ANG / AU2FMT)
        for i in range(na):
            vout.append(f"{el[i]}  {vv[i,0]:13.9f}{vv[i,1]:13.9f}{vv[i,2]:13.9f}\n")
    Path(stem + ".md.vel").write_text("".join(vout))

    if not total_energies:
        return None
    total_energies = np.asarray(total_energies)
    energy_drift = total_energies - total_energies[0]
    return float(energy_drift[-1]), float(np.max(np.abs(energy_drift)))


def main():
    global gto, md, se_mindo3
    ap = argparse.ArgumentParser(
        description="Tutorial cheap MD runner for ISCatter outputs (MINDO/3 only)."
    )
    ap.add_argument("--outdir", default="outputs", help="Directory with ISCatter output .xyz/.vel files")
    ap.add_argument("--prefix", default="out", help="File prefix, e.g. out -> out_0.xyz")
    ap.add_argument("--ntraj", type=int, default=10, help="How many trajectories to run")
    ap.add_argument("--dt", type=float, default=10.0, help="MD timestep in a.u.")
    ap.add_argument("--steps", type=int, default=500, help="Number of MD steps per trajectory")
    ap.add_argument("--charge", type=int, default=0, help="Total charge")
    ap.add_argument("--spin", type=int, default=0, help="2S spin multiplicity offset (0 singlet)")
    ap.add_argument("--ref-energy-file", default="ref.en", help="Optional reference energy file (Hartree)")
    args = ap.parse_args()

    warnings.filterwarnings(
        "ignore",
        message=r"Since PySCF-2\.3, B3LYP.*",
        category=UserWarning,
    )
    try:
        from pyscf import gto as _gto, md as _md
        gto, md = _gto, _md
        try:
            from pyscf.semiempirical import mindo3 as _se_mindo3
            se_mindo3 = _se_mindo3
        except Exception:
            se_mindo3 = None
    except Exception as exc:
        raise SystemExit(
            "PySCF import failed. Ensure compatible pyscf/h5py/numpy are installed in this environment.\n"
            f"Import error: {exc}"
        )

    pairs = discover_pairs(args.outdir, args.prefix)
    if not pairs:
        raise SystemExit(f"No matching trajectory inputs found in '{args.outdir}' for prefix '{args.prefix}'.")

    ref_en = 0.0
    if os.path.exists(args.ref_energy_file):
        try:
            ref_en = float(Path(args.ref_energy_file).read_text().splitlines()[0].strip())
        except Exception:
            ref_en = 0.0

    nrun = min(args.ntraj, len(pairs))
    print(f"Running {nrun} cheap trajectories from {args.outdir} (prefix={args.prefix})")
    for i, (xf, vf) in enumerate(pairs[:nrun]):
        print(f"[{i+1}/{nrun}] {xf}")
        energy_check = run_one(xf, vf, args.dt, args.steps, args.charge, args.spin, ref_en)
        if energy_check is not None:
            final_drift, max_drift = energy_check
            print(
                f"    total-energy drift: final={final_drift:+.6f} eV; "
                f"max |drift|={max_drift:.6f} eV"
            )


if __name__ == "__main__":
    main()
