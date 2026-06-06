#!/usr/bin/env python3
import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np


def _set_default_cache_dirs() -> None:
    user = os.environ.get("USER", "user")
    tmp = tempfile.gettempdir()
    os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(tmp, "numba_cache_" + user))
    os.environ.setdefault("MPLCONFIGDIR", os.path.join(tmp, "mpl_cache_" + user))
    os.makedirs(os.environ["NUMBA_CACHE_DIR"], exist_ok=True)
    os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)


def _replace_or_insert(text: str, key: str, value: str) -> str:
    rx = re.compile(rf"^{re.escape(key)}\s*=.*$", re.MULTILINE)
    line = f"{key} = {value}"
    if rx.search(text):
        return rx.sub(line, text)
    marker = "\n# Tutorial-specific options"
    if marker in text:
        return text.replace(marker, "\n" + line + "\n" + marker, 1)
    return text.rstrip() + "\n" + line + "\n"


def _patch_input(path: Path, frame: str) -> None:
    text = path.read_text()
    for key, value in (
        ("Nsamp", "1"),
        ("workers", "1"),
        ("printout", "1 1 0 0"),
        ("progress", "quiet"),
        ("plothist", "False"),
        ("hist_initial", "False"),
        ("hist_sampled", "False"),
        ("output-frame", frame),
        ("audit-initial-sample", "True"),
        ("audit-initial-energy-tol", "1.0e-5"),
        ("audit-initial-angular-tol", "1.0e-5"),
        ("audit-initial-vib-tol", "0.0"),
        ("audit-initial-velocity-tol", "1.0e-5"),
        ("wang", "False"),
    ):
        text = _replace_or_insert(text, key, value)
    path.write_text(text)


def _read_xyz_like(path: Path):
    lines = path.read_text().splitlines()
    n = int(lines[0])
    rows = []
    for line in lines[2:2 + n]:
        parts = line.split()
        rows.append((parts[0], np.array([float(x) for x in parts[1:4]], dtype=float)))
    return rows


def _max_rxpi_mismatch(internal: Path, plusz: Path) -> float:
    a = _read_xyz_like(internal)
    b = _read_xyz_like(plusz)
    if len(a) != len(b):
        raise ValueError(f"Different atom counts: {internal} vs {plusz}")
    worst = 0.0
    for (ela, va), (elb, vb) in zip(a, b):
        if ela != elb:
            raise ValueError(f"Element mismatch: {ela} vs {elb}")
        expected = np.array([va[0], -va[1], -va[2]])
        worst = max(worst, float(np.max(np.abs(vb - expected))))
    return worst


def _run(cmd, cwd: Path, log: Path) -> int:
    with log.open("w") as fh:
        proc = subprocess.run(cmd, cwd=cwd, stdout=fh, stderr=subprocess.STDOUT)
    return proc.returncode


def main() -> int:
    _set_default_cache_dirs()
    parser = argparse.ArgumentParser(
        prog="icats.frame-smoke",
        description="Check internal vs incoming-k-plus-z output-frame consistency on a fixed-b tutorial.",
    )
    parser.add_argument("--out-dir", default="", help="Output directory. Defaults to a temporary directory.")
    parser.add_argument("--tolerance", type=float, default=1.0e-10, help="Maximum allowed Cartesian Rx(pi) mismatch.")
    parser.add_argument("--keep", action="store_true", help="Keep temporary directory when --out-dir is omitted.")
    args = parser.parse_args()

    if args.out_dir:
        base = Path(args.out_dir).resolve()
        base.mkdir(parents=True, exist_ok=True)
        cleanup = False
    else:
        base = Path(tempfile.mkdtemp(prefix="icats_frame_smoke_"))
        cleanup = not args.keep

    try:
        template = base / "template"
        gen_rc = _run(
            [
                sys.executable,
                "-m",
                "icats.cli",
                "--tutorial",
                "fixed_plane_atom_diatom_ar_no",
                "--tutorial-dir",
                str(template),
                "--setup-only",
            ],
            base,
            base / "generate.log",
        )
        if gen_rc != 0:
            print(f"generation failed; see {base / 'generate.log'}")
            return 1

        cases = {"internal": base / "internal", "incoming-k-plus-z": base / "incoming-k-plus-z"}
        for frame, tdir in cases.items():
            if tdir.exists():
                shutil.rmtree(tdir)
            shutil.copytree(template, tdir)
            _patch_input(tdir / "tutorial_input.txt", frame)
            rc = _run([sys.executable, "-m", "icats.cli_init", "tutorial_input.txt"], tdir, tdir / "init.log")
            if rc != 0:
                print(f"{frame} run failed; see {tdir / 'init.log'}")
                return 1
            info = (tdir / "out_full.info").read_text(errors="ignore")
            if "Initial sample audit: OK" not in info:
                print(f"{frame} audit marker missing; see {tdir / 'out_full.info'}")
                return 1

        xyz_mismatch = _max_rxpi_mismatch(cases["internal"] / "out_full.xyz", cases["incoming-k-plus-z"] / "out_full.xyz")
        vel_mismatch = _max_rxpi_mismatch(cases["internal"] / "out_full.vel", cases["incoming-k-plus-z"] / "out_full.vel")
        print(f"frame smoke directory: {base}")
        print(f"out_full.xyz Rx(pi) max mismatch: {xyz_mismatch:.3e}")
        print(f"out_full.vel Rx(pi) max mismatch: {vel_mismatch:.3e}")
        if xyz_mismatch > args.tolerance or vel_mismatch > args.tolerance:
            print(f"mismatch exceeds tolerance {args.tolerance:.3e}")
            return 1
        print("frame smoke: OK")
        return 0
    finally:
        if cleanup:
            shutil.rmtree(base, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
