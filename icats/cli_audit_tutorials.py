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
    ):
        text = _replace_or_insert(text, key, value)
    path.write_text(text)


def _parse_audit(path: Path):
    text = path.read_text(errors="ignore")
    ok_count = text.count("Initial sample audit: OK")
    rx = re.compile(
        r"Audit energy\s+(\w+): generation\s+([-+0-9.eE]+) eV, "
        r"analysis\s+([-+0-9.eE]+) eV, diff\s+([-+0-9.eE]+)"
    )
    worst = (0.0, "")
    for match in rx.finditer(text):
        diff = float(match.group(4))
        if diff > worst[0]:
            worst = (diff, match.group(1))
    return ok_count, worst


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
    parser.add_argument("--angular-tol", type=float, default=0.0, help="Angular audit tolerance; 0 disables angular checks.")
    parser.add_argument("--keep-going", action="store_true", help="Continue after a failed tutorial.")
    args = parser.parse_args()

    if args.out_dir:
        base = Path(args.out_dir).resolve()
    else:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        base = Path.cwd() / "smoke_results" / f"initial_audit_{stamp}"
    base.mkdir(parents=True, exist_ok=True)

    rows = ["tutorial\tgenerate_rc\tinit_rc\taudit_ok\tmax_diff_eV\tworst_component\tnote\n"]
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
            rows.append(f"{name}\t{gen.returncode}\tNA\tno\tNA\tNA\tgeneration failed\n")
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
        max_diff = "NA"
        worst_comp = "NA"
        audit_ok = "no"
        if out_info.exists():
            ok_count, worst = _parse_audit(out_info)
            max_diff = f"{worst[0]:.6g}"
            worst_comp = worst[1] or "none"
            audit_ok = "yes" if init.returncode == 0 and ok_count == args.nsamp else "no"
            note = f"{ok_count}/{args.nsamp} audit OK"
        else:
            note = "missing out_full.info"
        rows.append(f"{name}\t{gen.returncode}\t{init.returncode}\t{audit_ok}\t{max_diff}\t{worst_comp}\t{note}\n")
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
        "- Dynamics were not run.\n\n"
        "See `summary.tsv` for per-tutorial results.\n"
    )
    print(f"Audit directory: {base}")
    print(summary.read_text(), end="")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
