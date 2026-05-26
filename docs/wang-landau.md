---
layout: default
title: Wang-Landau Sampling
---

# Wang-Landau Sampling

Wang-Landau sampling is used when direct sampling does not populate the desired
total angular momentum range well enough. This is common when molecular angular
momenta and orbital angular momentum combine to produce a nontrivial density of
states for `J`.

In a molecular collision, `J = L + Jab` as a vector relation. Sampling `L` from
the impact-parameter distribution does not automatically produce the desired
triangular total-`J` distribution, especially at low `J` where the molecular
rotor sum `Jab` is comparable to the orbital angular momentum. Wang-Landau is
used to learn the density-of-states correction needed for this acceptance step.

Enable Wang-Landau sampling with:

```text
wang = True
```

The generated umbrella is stored as:

```text
rd_<run-tag>/wang.pkl
```

The file stores the umbrella arrays plus metadata describing the settings used
to generate them. If the requested input range is incompatible with the stored
file, ICATS refuses to continue rather than silently reusing the wrong
umbrella.

Useful controls:

```text
wlmode = fast
wl-angular-sampler = fast
wl-audit-angular-sampler = False
wl-tol = 1.000001
wl-flatness = 0.85
```

`wlmode` is the Wang-Landau parameter preset. `wl-angular-sampler` controls the
numerical implementation used to sample angular momenta during umbrella
generation.

`wlmode = fast` changes the Wang-Landau convergence settings. It is not the same
thing as `wl-angular-sampler = fast`, which controls the numerical shortcut used
to draw angular momenta during the repeated Wang-Landau trial steps.

## Typical Workflow

1. Start from a tutorial that already has Wang-Landau enabled.
2. Run `icats.init tutorial_input.txt`.
3. If `wang.pkl` is absent, ICATS generates it.
4. If `wang.pkl` is present and compatible, ICATS reuses it.
5. Inspect the Wang-Landau plots and sampled `J`, `L` histograms.

For a new molecular system, first run a small calculation and check the plots
before committing to a large sample count.

## Stored Umbrellas

`wang.pkl` is not just a raw array. It includes metadata such as the requested
angular momentum range and Wang-Landau settings. ICATS refuses to reuse an
incompatible umbrella.

If you need to regenerate the umbrella, move the old file aside:

```bash
mv rd_tutorial_input/wang.pkl rd_tutorial_input/wang_previous.pkl
icats.init tutorial_input.txt
```

This avoids accidentally overwriting a useful precomputed umbrella.

## Parameter Guide

| Key | Meaning | Beginner advice |
| --- | --- | --- |
| `wang` | Enables Wang-Landau weighting | Use `True` only when needed |
| `wlmode` | Parameter preset | Start with tutorial defaults |
| `wl-angular-sampler` | Numerical implementation | Keep `fast` unless debugging |
| `wl-audit-angular-sampler` | Compare fast and legacy WL angular samplers | Use only for debugging |
| `wl-flatness` | Histogram flatness target | Higher is stricter |
| `wl-tol` | Stopping tolerance for modification factor | Lower tolerance takes longer |
| `wl-wn` | Number of bins, if set directly | Leave unset at first |

## Plotting Checks

After a Wang-Landau run, inspect the learned distribution and weights:

```bash
cd rd_tutorial_input/histograms/wl
python wl_td_plot.py
python wl_wl_plot.py
```

Initial and sampled angular-momentum histograms can be generated with:

```bash
./rd_tutorial_input/histograms/plot_initial.sh
./rd_tutorial_input/histograms/plot_sampled.sh
```

These plots are a practical sanity check that the requested `J` and `L` ranges
are populated sensibly.

## What To Look For

A useful Wang-Landau run should show:

- sampled `J` values spanning the requested useful range,
- no sharp unexplained cutoff inside the intended range,
- `L` values broad enough for the intended impact parameters,
- an umbrella that changes smoothly enough to be physically plausible.

If the run refuses to start because `wang.pkl` is too short or incompatible,
that is a successful safety check rather than a crash. Regenerate the umbrella
for the new requested range.

For a tutorial-quality stored umbrella, the aim is not perfect production
convergence. The aim is that the resulting sampled `J` and `L` histograms are
smooth, cover the requested region, and do not show obvious boundary artefacts
that would hide mistakes in the rest of the workflow.
