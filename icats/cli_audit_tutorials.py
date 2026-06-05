#!/usr/bin/env python3
import argparse
import os
import re
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path

from .cli import TUTORIAL_ORDER


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


def _patch_input(path: Path, args: argparse.Namespace) -> None:
    text = path.read_text()
    for key, value in (
        ("Nsamp", str(args.nsamp)),
        ("workers", str(args.workers)),
        ("progress", "quiet"),
        ("printout", "0 1 0 0"),
        ("plothist", "False"),
        ("hist_initial", "False"),
        ("hist_sampled", "False"),
        ("audit-initial-sample", "True"),
        ("audit-initial-energy-tol", str(args.energy_tol)),
        ("audit-initial-angular-tol", str(args.angular_tol)),
        ("audit-initial-vib-tol", str(args.vib_tol)),
        ("audit-initial-velocity-tol", str(args.velocity_tol)),
    ):
        text = _replace_or_insert(text, key, value)
    if not args.include_wl:
        text = _replace_or_insert(text, "wang", "False")
    path.write_text(text)


def _parse_audit(path: Path):
    text = path.read_text(errors="ignore")
    ok_count = text.count("Initial sample audit: OK")
    energy_rx = re.compile(
        r"Audit energy\s+(\w+): generation\s+([-+0-9.eE]+) eV, "
        r"analysis\s+([-+0-9.eE]+) eV, diff\s+([-+0-9.eE]+)"
    )
    generic_rx = re.compile(
        r"Audit\s+(angular|vector|mol scalar|mol vector|scalar|vib|angle)\s+(.+?):.*?"
        r"(?:diff|vector-norm diff|component-rms diff|circular diff)\s+([-+0-9.eE]+)"
    )
    worst_energy = (0.0, "")
    worst_state = (0.0, "")
    for match in energy_rx.finditer(text):
        diff = float(match.group(4))
        if diff > worst_energy[0]:
            worst_energy = (diff, match.group(1))
    for match in generic_rx.finditer(text):
        diff = float(match.group(3))
        if diff > worst_state[0]:
            worst_state = (diff, match.group(1).strip() + ":" + match.group(2).strip())
    return ok_count, worst_energy, worst_state


def main() -> int:
    _set_default_cache_dirs()
    parser = argparse.ArgumentParser(
        prog="icats.audit-tutorials",
        description="Generate all tutorials and run the initial-sample audit without dynamics.",
    )
    parser.add_argument("--out-dir", default="", help="Output directory for the audit run.")
    parser.add_argument("--nsamp", type=int, default=2, help="Samples per tutorial.")
    parser.add_argument("--workers", type=int, default=1, help="Workers per tutorial.")
    parser.add_argument("--energy-tol", type=float, default=0.02, help="Energy audit tolerance in eV.")
    parser.add_argument("--angular-tol", type=float, default=2.0, help="Angular audit tolerance in au; 0 disables angular checks.")
    parser.add_argument("--vib-tol", type=float, default=2.0, help="Vibrational Q/P RMS-per-component audit tolerance; 0 records diagnostics without failing.")
    parser.add_argument("--velocity-tol", type=float, default=5.0, help="Relative-velocity audit tolerance in m/s.")
    parser.add_argument("--include-wl", action="store_true", help="Keep Wang-Landau enabled in tutorials that request it.")
    parser.add_argument("--keep-going", action="store_true", help="Continue after a failed tutorial.")
    args = parser.parse_args()

    if args.out_dir:
        base = Path(args.out_dir).resolve()
    else:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        base = Path.cwd() / "smoke_results" / f"initial_audit_{stamp}"
    base.mkdir(parents=True, exist_ok=True)

    rows = ["tutorial\tgenerate_rc\tinit_rc\taudit_ok\tmax_energy_diff_eV\tworst_energy\tmax_state_diff\tworst_state\tnote\n"]
    failures = []
    for name in TUTORIAL_ORDER:
        tdir = base / name
        gen_log = base / f"{name}.generate.log"
        with gen_log.open("w") as fh:
            gen = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "icats.cli",
                    "--tutorial",
                    name,
                    "--tutorial-dir",
                    str(tdir),
                    "--nsamp",
                    str(args.nsamp),
                    "--ntraj",
                    "1",
                    "--setup-only",
                ],
                stdout=fh,
                stderr=subprocess.STDOUT,
            )
        if gen.returncode != 0:
            rows.append(f"{name}\t{gen.returncode}\tNA\tno\tNA\tNA\tNA\tNA\tgeneration failed\n")
            failures.append(name)
            if not args.keep_going:
                break
            continue

        input_file = tdir / "tutorial_input.txt"
        _patch_input(input_file, args)
        init_log = tdir / "smoke_init.log"
        with init_log.open("w") as fh:
            init = subprocess.run(
                [sys.executable, "-m", "icats.cli_init", "tutorial_input.txt"],
                cwd=tdir,
                stdout=fh,
                stderr=subprocess.STDOUT,
            )

        out_info = tdir / "out_full.info"
        note = ""
        max_energy = "NA"
        worst_energy = "NA"
        max_state = "NA"
        worst_state = "NA"
        audit_ok = "no"
        if out_info.exists():
            ok_count, worst_e, worst_s = _parse_audit(out_info)
            max_energy = f"{worst_e[0]:.6g}"
            worst_energy = worst_e[1] or "none"
            max_state = f"{worst_s[0]:.6g}"
            worst_state = worst_s[1] or "none"
            audit_ok = "yes" if init.returncode == 0 and ok_count == args.nsamp else "no"
            note = f"{ok_count}/{args.nsamp} audit OK"
        else:
            note = "missing out_full.info"
        rows.append(
            f"{name}\t{gen.returncode}\t{init.returncode}\t{audit_ok}\t"
            f"{max_energy}\t{worst_energy}\t{max_state}\t{worst_state}\t{note}\n"
        )
        if audit_ok != "yes":
            failures.append(name)
            if not args.keep_going:
                break

    summary = base / "summary.tsv"
    summary.write_text("".join(rows))
    readme = base / "README.md"
    readme.write_text(
        "# ICATS initial-sample audit\n\n"
        f"- Samples per tutorial: `{args.nsamp}`\n"
        f"- Workers: `{args.workers}`\n"
        f"- Energy tolerance: `{args.energy_tol}` eV\n"
        f"- Angular tolerance: `{args.angular_tol}`\n"
        f"- Vibrational Q/P RMS-per-component tolerance: `{args.vib_tol}`\n"
        f"- Relative-velocity tolerance: `{args.velocity_tol}` m/s\n"
        f"- Wang-Landau enabled: `{args.include_wl}`\n"
        "- Dynamics were not run.\n\n"
        "See `summary.tsv` for per-tutorial results.\n"
    )
    print(f"Audit directory: {base}")
    print(summary.read_text(), end="")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
