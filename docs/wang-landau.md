---
layout: default
title: Wang-Landau Sampling
---

# Wang-Landau Sampling

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
