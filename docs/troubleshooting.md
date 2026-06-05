---
layout: default
title: Troubleshooting
---

# Troubleshooting

This page lists common problems a new user may see while learning ICATS.

## `icats` command not found

The package is probably not installed in the active Python environment. From the
repository root, run:

```bash
python -m pip install -e .
```

Then check:

```bash
icats --list-tutorials
```

## PySCF or MINDO errors

The initial-condition generator does not require the cheap dynamics backend. If
`run_cheap_dynamics.sh` fails, first check that the initial-condition stage still
works:

```bash
icats.init tutorial_input.txt
```

If only dynamics fails, install the optional dynamics dependencies:

```bash
python -m pip install -e ".[dynamics]"
```

If a binary compatibility error mentions NumPy, create a fresh environment and
reinstall ICATS there.

## Matplotlib or numba cache warnings

ICATS command-line tools set temporary default cache directories when needed.
If running Python modules directly, set:

```bash
export NUMBA_CACHE_DIR=/tmp/numba_cache_$USER
export MPLCONFIGDIR=/tmp/mpl_cache_$USER
```

This is most useful on clusters or filesystems where home-directory cache writes
are slow or restricted.

## Wang-Landau umbrella is incompatible

If ICATS refuses to use `wang.pkl`, the stored umbrella does not match the
requested calculation. This is intentional. Move the old file aside and rerun:

```bash
mv rd_tutorial_input/wang.pkl rd_tutorial_input/wang_old.pkl
icats.init tutorial_input.txt
```

Do not overwrite a useful `wang.pkl` unless you are sure you no longer need it.

## No histogram plots appear

Histogram scripts are generated only when histogram output is enabled. Check the
input file for:

```text
plothist = True
hist_initial = True
hist_sampled = True
```

Then run:

```bash
./rd_tutorial_input/histograms/plot_initial.sh
./rd_tutorial_input/histograms/plot_sampled.sh
```

For only the sampled system `L`, `J`, and `Jab` plots, use:

```bash
./rd_tutorial_input/histograms/plot_orbital_jljab.sh
```

## Audit fails

The initial-sample audit checks whether generated model energies are recovered
from the immediate Cartesian analysis. If it fails:

1. Confirm you are using the current code.
2. Rerun with `workers = 1` to simplify the log.
3. Inspect `out_full.info` for the failed sample.
4. Compare generated energies with the `Audit energy ...` lines.

For vibrating polyatomics, compare sampled rotor energy with
`vector rot. energy`, not `full rot. energy`.
