#!/usr/bin/env python3
import os


def write_histogram_helpers_runtime(hist_root):
    """Emit histogram helper scripts in rd_<tag>/histograms for runtime runs."""
    hroot = str(hist_root)
    os.makedirs(hroot, exist_ok=True)

    def _w(path, text, exe=False):
        with open(path, "w") as f:
            f.write(text)
        if exe:
            os.chmod(path, 0o755)

    _w(
        os.path.join(hroot, "plot_initial.sh"),
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "SCRIPT_DIR=\"$(cd \"$(dirname \"${BASH_SOURCE[0]}\")\" && pwd)\"\n"
        "HROOT=\"${HIST_ROOT:-$SCRIPT_DIR}\"\n"
        "mkdir -p \"$HROOT/plots/initial\"\n"
        "find \"$HROOT/initial\" -type f -name 'hist_*.py' 2>/dev/null | sort | while read -r f; do\n"
        "  stem=\"$(basename \"${f%.py}\")\"\n"
        "  python \"$f\" --no-show --outfile \"$HROOT/plots/initial/${stem}\"\n"
        "done\n"
        "echo \"Initial histogram plots written to $HROOT/plots/initial/\"\n",
        exe=True,
    )

    _w(
        os.path.join(hroot, "plot_sampled.sh"),
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "SCRIPT_DIR=\"$(cd \"$(dirname \"${BASH_SOURCE[0]}\")\" && pwd)\"\n"
        "HROOT=\"${HIST_ROOT:-$SCRIPT_DIR}\"\n"
        "mkdir -p \"$HROOT/plots/sampled\"\n"
        "find \"$HROOT/sampled\" -type f -name 'hist_*.py' 2>/dev/null | sort | while read -r f; do\n"
        "  stem=\"$(basename \"${f%.py}\")\"\n"
        "  python \"$f\" --no-show --outfile \"$HROOT/plots/sampled/${stem}\"\n"
        "done\n"
        "echo \"Sampled histogram plots written to $HROOT/plots/sampled/\"\n",
        exe=True,
    )

    _w(
        os.path.join(hroot, "plot_all.sh"),
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "SCRIPT_DIR=\"$(cd \"$(dirname \"${BASH_SOURCE[0]}\")\" && pwd)\"\n"
        "\"$SCRIPT_DIR/plot_initial.sh\"\n"
        "\"$SCRIPT_DIR/plot_sampled.sh\"\n",
        exe=True,
    )

    _w(
        os.path.join(hroot, "plot_compare_pairs.py"),
        "#!/usr/bin/env python3\n"
        "from __future__ import annotations\n"
        "import argparse, ast, os, re, subprocess\n"
        "from pathlib import Path\n"
        "import matplotlib.image as mpimg\n"
        "import matplotlib.pyplot as plt\n"
        "import numpy as np\n"
        "SCRIPT_DIR = Path(__file__).resolve().parent\n"
        "ROOT = Path(os.environ.get('HIST_ROOT', str(SCRIPT_DIR))).resolve()\n"
        "INI_DIR = ROOT / 'initial'\n"
        "SAM_DIR = ROOT / 'sampled'\n"
        "OUT = ROOT / 'plots' / 'compare'\n"
        "OUT_INI = OUT / 'init_single'\n"
        "OUT_SAM = OUT / 'sam_single'\n"
        "OUT_PAIR = OUT / 'pairs'\n"
        "OUT_UNM_INI = OUT / 'unmatched_init'\n"
        "OUT_UNM_SAM = OUT / 'unmatched_sam'\n"
        "for d in (OUT, OUT_INI, OUT_SAM, OUT_PAIR, OUT_UNM_INI, OUT_UNM_SAM):\n"
        "    d.mkdir(parents=True, exist_ok=True)\n"
        "ap = argparse.ArgumentParser()\n"
        "ap.add_argument('--include-unmatched', action='store_true')\n"
        "ap.add_argument('--include-rot', action='store_true')\n"
        "ap.add_argument('--include-vib', action='store_true')\n"
        "args = ap.parse_args()\n"
        "rx = re.compile(r'^hist_(ini|sam)_([^_]+)_(.+)\\.py$')\n"
        "def parse_scripts(base: Path):\n"
        "    out = {}\n"
        "    for p in sorted(base.rglob('hist_*.py')):\n"
        "        m = rx.match(p.name)\n"
        "        if not m:\n"
        "            continue\n"
        "        _st, scp, met = m.groups()\n"
        "        out[(scp, met)] = p\n"
        "    return out\n"
        "def render_hist(script_path: Path, out_base: Path):\n"
        "    subprocess.run(['python', str(script_path), '--no-show', '--outfile', str(out_base)], check=True, cwd=str(ROOT))\n"
        "    return out_base.with_suffix('.png')\n"
        "def side_by_side(lp: Path, rp: Path, title: str, out_png: Path):\n"
        "    li = mpimg.imread(lp)\n"
        "    ri = mpimg.imread(rp)\n"
        "    fig, axs = plt.subplots(1, 2, figsize=(12, 5))\n"
        "    axs[0].imshow(li); axs[0].axis('off'); axs[0].set_title('Initial')\n"
        "    axs[1].imshow(ri); axs[1].axis('off'); axs[1].set_title('Sampled')\n"
        "    fig.suptitle(title); fig.tight_layout(); fig.savefig(out_png, dpi=150); plt.close(fig)\n"
        "def extract_data(script_path: Path):\n"
        "    txt = script_path.read_text()\n"
        "    m = re.search(r'data\\s*=\\s*np\\.array\\((\\[[\\s\\S]*?\\])\\)', txt)\n"
        "    if not m: return None\n"
        "    return np.asarray(ast.literal_eval(m.group(1)), dtype=float).ravel()\n"
        "def dist_stats(a, b):\n"
        "    a = np.asarray(a, dtype=float).ravel(); b = np.asarray(b, dtype=float).ravel()\n"
        "    if len(a) < 2 or len(b) < 2: return {'bins':0, 'jsd':float('nan'), 'emd':float('nan')}\n"
        "    mix = np.concatenate([a,b])\n"
        "    edges = np.histogram_bin_edges(mix, bins='fd')\n"
        "    if len(edges) < 3: edges = np.histogram_bin_edges(mix, bins=20)\n"
        "    ha,_ = np.histogram(a, bins=edges, density=False); hb,_ = np.histogram(b, bins=edges, density=False)\n"
        "    eps = 1e-12\n"
        "    pa = (ha+eps)/np.sum(ha+eps); pb = (hb+eps)/np.sum(hb+eps); mm = 0.5*(pa+pb)\n"
        "    jsd = np.sqrt(0.5*(np.sum(pa*np.log(pa/mm))+np.sum(pb*np.log(pb/mm))))\n"
        "    cdfa = np.cumsum(pa); cdfb = np.cumsum(pb)\n"
        "    emd = np.sum(np.abs(cdfa-cdfb)*np.diff(edges))\n"
        "    return {'bins':int(len(edges)-1),'jsd':float(jsd),'emd':float(emd)}\n"
        "ini = parse_scripts(INI_DIR); sam = parse_scripts(SAM_DIR)\n"
        "if len(ini) == 0 and len(sam) == 0:\n"
        "    print(f'No histogram scripts found under {ROOT}.')\n"
        "    print('Run icats.init first (with histogram mode enabled), then re-run this script.')\n"
        "    raise SystemExit(0)\n"
        "sys_map = {'j':'orb_ij','l':'orb_il','b':'orb_sb','phi':'orb_sphi','vel':'vel_ivel'}\n"
        "sys_fallback = {'j':'orb_sj','l':'orb_sl'}\n"
        "sys_default = {'j','l','vel','b','phi'}\n"
        "pairs=[]; used_ini=set(); used_sam=set()\n"
        "for (scp,met), ip in ini.items():\n"
        "    if scp=='sys':\n"
        "        if met not in sys_default: continue\n"
        "        skey=(scp, sys_map.get(met, met))\n"
        "        if skey not in sam and met in sys_fallback:\n"
        "            skey=(scp, sys_fallback[met])\n"
        "    elif scp in ('m0','m1'):\n"
        "        if met.startswith('vi_mode') and not args.include_vib: continue\n"
        "        if met=='j' and not args.include_rot: continue\n"
        "        if not (met.startswith('vi_mode') or met=='j'): continue\n"
        "        skey=(scp,met)\n"
        "    else:\n"
        "        continue\n"
        "    if skey in sam:\n"
        "        pairs.append(((scp,met), skey, ip, sam[skey])); used_ini.add((scp,met)); used_sam.add(skey)\n"
        "manifest=[]\n"
        "for (iscp,imet),(_ssc,smet), ipath, spath in pairs:\n"
        "    key=f'{iscp}_{imet}'\n"
        "    ipng=render_hist(ipath, OUT_INI / f'{key}_ini')\n"
        "    spng=render_hist(spath, OUT_SAM / f'{key}_sam')\n"
        "    cpng=OUT_PAIR / f'{key}__ini_vs_sam.png'\n"
        "    desc=f'Compare {iscp}:{imet} (init) vs {iscp}:{smet} (sampled)'\n"
        "    side_by_side(ipng, spng, desc, cpng)\n"
        "    ida=extract_data(ipath); sda=extract_data(spath)\n"
        "    st=dist_stats(ida,sda) if ida is not None and sda is not None else {'bins':0,'jsd':float('nan'),'emd':float('nan')}\n"
        "    manifest.append((key, desc, str(st['bins']), f\"{st['jsd']:.6g}\", f\"{st['emd']:.6g}\", str(cpng)))\n"
        "if args.include_unmatched:\n"
        "    for (scp,met), p in ini.items():\n"
        "        if (scp,met) in used_ini: continue\n"
        "        render_hist(p, OUT_UNM_INI / f'{scp}_{met}_ini')\n"
        "    for (scp,met), p in sam.items():\n"
        "        if (scp,met) in used_sam: continue\n"
        "        render_hist(p, OUT_UNM_SAM / f'{scp}_{met}_sam')\n"
        "mf = OUT / 'pair_manifest.tsv'\n"
        "hdr = 'key\\tdescription\\tbins\\tjsd\\temd\\tcomparison_png\\n'\n"
        "mf.write_text(hdr + '\\n'.join('\\t'.join(x) for x in manifest) + '\\n')\n"
        "print(f'Paired comparisons: {len(manifest)}')\n"
        "print(f'Manifest: {mf}')\n"
        "print(f'Pair plots: {OUT_PAIR}')\n"
        "print('Unmatched rendered.' if args.include_unmatched else 'Unmatched skipped (use --include-unmatched).')\n",
        exe=True,
    )

    _w(
        os.path.join(hroot, "plot_compare_pairs.sh"),
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "SCRIPT_DIR=\"$(cd \"$(dirname \"${BASH_SOURCE[0]}\")\" && pwd)\"\n"
        "export HIST_ROOT=\"${HIST_ROOT:-$SCRIPT_DIR}\"\n"
        "python \"$SCRIPT_DIR/plot_compare_pairs.py\" \"$@\"\n",
        exe=True,
    )
