#!/usr/bin/env python3
import argparse
import glob
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
    this_file = Path(__file__).resolve()
    repo_parent = str(this_file.parent.parent)
    if repo_parent not in sys.path:
        sys.path.insert(0, repo_parent)


def _parse_index(prefix: str, base: str) -> int:
    if not base.startswith(prefix + "_"):
        return 10**9
    tail = base[len(prefix) + 1 :]
    token = tail.split(".", 1)[0]
    try:
        return int(token)
    except Exception:
        return 10**9


def discover_pairs(directory: str, prefix: str):
    by_idx = {}

    # Dynamics files: out_*.md.xyz + out_*.md.vel
    for xf in glob.glob(os.path.join(directory, f"{prefix}_*.md.xyz")):
        vf = xf[:-7] + ".md.vel"
        if os.path.exists(vf):
            base = os.path.basename(xf)
            idx = _parse_index(prefix, base)
            by_idx[idx] = (0, xf, vf)

    # Plain files: out_*.xyz + out_*.vel
    for xf in glob.glob(os.path.join(directory, f"{prefix}_*.xyz")):
        if xf.endswith(".md.xyz"):
            continue
        vf = xf[:-4] + ".vel"
        if os.path.exists(vf):
            base = os.path.basename(xf)
            idx = _parse_index(prefix, base)
            if idx not in by_idx:
                by_idx[idx] = (1, xf, vf)

    pairs = [(idx, rank, xf, vf) for idx, (rank, xf, vf) in by_idx.items()]
    pairs.sort(key=lambda t: (t[0], t[1], t[2]))
    return [(x, v) for _, _, x, v in pairs]


def main() -> int:
    _set_default_cache_dirs()
    _ensure_repo_parent_on_path()
    ap = argparse.ArgumentParser(
        prog="icats.analyse",
        description="Analyse concordant *xyz/*vel pairs for a prefix.",
    )
    ap.add_argument("input_file", help="icats input file used to define molecule metadata")
    ap.add_argument("--dir", default=".", help="Directory to scan for pairs (default: current dir)")
    ap.add_argument("--prefix", default="out", help="File prefix (default: out)")
    ap.add_argument("--ntraj", type=int, default=0, help="Max number of pairs to analyse (0 = all)")
    ap.add_argument("--out-prefix", default="dynamics", help="Output prefix for analysis logs")
    args = ap.parse_args()

    from icats.iscattering import icats

    sc = icats()
    sc.ReadInput(args.input_file)

    scan_dir = Path(args.dir).resolve()
    if not scan_dir.is_dir():
        raise SystemExit(f"Directory not found: {scan_dir}")

    pairs = discover_pairs(str(scan_dir), args.prefix)
    if not pairs:
        raise SystemExit(f"No concordant pairs found in {scan_dir} for prefix '{args.prefix}'")

    if args.ntraj > 0:
        pairs = pairs[: args.ntraj]

    analysed = 0
    for i, (xf, vf) in enumerate(pairs):
        sa = sc.InitializeWorker(i)
        sc.ReadSamples(sa, xf, vf)
        (scan_dir / f"{args.out_prefix}{i}.analinfo").write_text("".join(sa.slog))
        analysed += 1

    print(f"Analysed {analysed} trajectories in {scan_dir}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
