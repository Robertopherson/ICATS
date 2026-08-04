#!/usr/bin/env python3
import os
import numpy as np


def write_costheta_diagnostics(conv_dir, vals, js=None, prefix="costheta"):
    """Write L-Jab angular-correlation diagnostics and optional heatmap helper.

    Returns a small dict with summary values and a warning flag.
    """
    if vals is None or len(vals) == 0:
        return {"written": False, "warning": False}
    raw_arr = np.asarray(vals, dtype=float)
    mask = np.isfinite(raw_arr) & (np.abs(raw_arr) <= 1.0)
    arr = raw_arr[mask]
    if arr.size == 0:
        return {"written": False, "warning": False}

    os.makedirs(conv_dir, exist_ok=True)
    n = int(arr.size)
    mean = float(np.mean(arr))
    std = float(np.std(arr, ddof=1)) if n > 1 else 0.0
    sem = float(std / np.sqrt(n)) if n > 0 else 0.0
    mean_abs = float(np.mean(np.abs(arr)))
    warn_limit = max(0.05, 5.0 * sem)
    warning = abs(mean) > warn_limit

    summ = [
        "# cos(theta_{L,Jab}) angular-correlation summary\n",
        "# For an isotropic uncorrelated relative angle, mean -> 0 and mean_abs -> 0.5.\n",
        f"n = {n}\n",
        f"mean = {mean:.8e}\n",
        f"std = {std:.8e}\n",
        f"sem = {sem:.8e}\n",
        f"mean_abs = {mean_abs:.8e}\n",
        f"warning_limit_abs_mean = {warn_limit:.8e}\n",
        f"warning = {warning}\n",
    ]
    if warning:
        summ.append(
            "WARNING: |mean cos(theta_{L,Jab})| is larger than the heuristic "
            "limit; inspect J-vs-cos(theta) plots for WL-induced angular correlation.\n"
        )
    with open(os.path.join(conv_dir, prefix + "_summary.txt"), "w") as f:
        f.writelines(summ)

    csum = np.cumsum(arr)
    csum2 = np.cumsum(arr * arr)
    abs_csum = np.cumsum(np.abs(arr))
    lines = ["sample_count\tmean\tstd\tsem\tmean_abs\n"]
    for i in range(1, n + 1):
        mu = csum[i - 1] / i
        if i > 1:
            var = max((csum2[i - 1] - i * mu * mu) / (i - 1), 0.0)
            si = np.sqrt(var)
        else:
            si = 0.0
        sei = si / np.sqrt(i)
        mai = abs_csum[i - 1] / i
        lines.append(f"{i}\t{mu:.8e}\t{si:.8e}\t{sei:.8e}\t{mai:.8e}\n")
    with open(os.path.join(conv_dir, prefix + "_cumulative.tsv"), "w") as f:
        f.writelines(lines)

    if js is not None and len(js) == len(vals):
        jarr = np.asarray(js, dtype=float)
        pair_mask = mask & np.isfinite(jarr)
        if np.count_nonzero(pair_mask) > 1:
            table = np.column_stack([jarr[pair_mask], raw_arr[pair_mask]])
            np.savetxt(
                os.path.join(conv_dir, prefix + "_j_vs_costheta.tsv"),
                table,
                delimiter="\t",
                header="J\tcos_theta_L_Jab",
                comments="",
            )
            write_costheta_heatmap_script(conv_dir, prefix)

    return {
        "written": True,
        "warning": warning,
        "n": n,
        "mean": mean,
        "std": std,
        "sem": sem,
        "mean_abs": mean_abs,
        "warning_limit": warn_limit,
        "summary_path": os.path.join(conv_dir, prefix + "_summary.txt"),
    }


def write_costheta_heatmap_script(conv_dir, prefix="costheta"):
    script = os.path.join(conv_dir, prefix + "_j_vs_costheta_plot.py")
    lines = [
        "#!/usr/bin/env python3\n",
        "import argparse\n",
        "import pathlib\n",
        "import numpy as np\n",
        "import matplotlib.pyplot as plt\n",
        "\n",
        "p = argparse.ArgumentParser()\n",
        "p.add_argument(\"--no-show\", action=\"store_true\")\n",
        f"p.add_argument(\"--outfile\", default=\"{prefix}_j_vs_costheta\")\n",
        "p.add_argument(\"--bins-j\", type=int, default=50)\n",
        "p.add_argument(\"--bins-cos\", type=int, default=50)\n",
        "args = p.parse_args()\n",
        "root = pathlib.Path(__file__).resolve().parent\n",
        f"data = np.loadtxt(root / \"{prefix}_j_vs_costheta.tsv\", skiprows=1)\n",
        "if data.ndim == 1:\n",
        "    data = data.reshape(1, -1)\n",
        "J = data[:, 0]\n",
        "C = data[:, 1]\n",
        "fig, ax = plt.subplots(figsize=(7, 5))\n",
        "h = ax.hist2d(J, C, bins=[args.bins_j, args.bins_cos], cmap=\"viridis\")\n",
        "fig.colorbar(h[3], ax=ax, label=\"counts\")\n",
        "ax.axhline(0.0, color=\"w\", lw=1, alpha=0.8)\n",
        "ax.set_xlabel(\"J\")\n",
        "ax.set_ylabel(r\"$\\cos\\theta_{L,Jab}$\")\n",
        "ax.set_title(r\"$J$ vs. $\\cos\\theta_{L,Jab}$\")\n",
        "fig.tight_layout()\n",
        "out = pathlib.Path(args.outfile)\n",
        "fig.savefig(out.with_suffix(\".png\"), dpi=160)\n",
        "print(\"Wrote\", out.with_suffix(\".png\"))\n",
        "if not args.no_show:\n",
        "    plt.show()\n",
    ]
    with open(script, "w") as f:
        f.writelines(lines)
    os.chmod(script, 0o755)
