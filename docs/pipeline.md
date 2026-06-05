---
layout: default
title: Dynamics and Analysis Pipeline
---

# Dynamics and Analysis Pipeline

![ICATS workflow](assets/pipeline.svg)

The ICATS workflow has three stages:

1. Generate initial conditions with `icats.init`.
2. Optionally run cheap demonstration dynamics.
3. Analyse trajectory coordinate and velocity files with `icats.analyse`.

The tutorial helper scripts show the full flow:

```bash
icats.init tutorial_input.txt
./run_cheap_dynamics.sh
./run_analysis.sh
```

The cheap dynamics are intended as a toy method for learning the interface and
file flow. Energy drift in those trajectories reflects the dynamics/integrator
and semiempirical potential choices, not necessarily an initial-condition
problem.

The initial-condition audit described in [Diagnostics](diagnostics.md) checks
the sampling and analysis round trip before any dynamics are run.

This separation is important. If the audit passes, the initial Cartesian sample
is internally consistent with the model used to generate it. A later energy
drift in MINDO dynamics is then a dynamics, integration, or potential-energy
surface issue, not direct evidence that the initial-condition sampler is broken.

For a quick regression check of the tutorial initial-condition pipeline:

```bash
icats.audit-tutorials --nsamp 3 --keep-going
```

This command does not call the cheap dynamics backend. It only generates each
tutorial, runs `icats.init`, and checks that sampled model energies and
comparable state variables are recovered from the immediate Cartesian analysis.
Wang-Landau is kept off by default so the command remains a quick regression
check after code changes.

## Files Produced By Each Stage

After `icats --tutorial ... --setup-only`:

```text
tutorial_input.txt
run_cheap_dynamics.sh
run_analysis.sh
molecule data files
```

After `icats.init tutorial_input.txt`:

```text
rd_tutorial_input/dat_tutorial_input.txt.pkl
rd_tutorial_input/work_sys_tutorial_input.txt.pkl
rd_tutorial_input/outputs/out_*.xyz or out_*.md.xyz
rd_tutorial_input/outputs/out_*.vel or out_*.md.vel
```

For Wang-Landau tutorials, the run directory can also contain:

```text
rd_tutorial_input/wang.pkl
rd_tutorial_input/histograms/wl/
```

After `./run_analysis.sh`:

```text
rd_tutorial_input/outputs/dynamics*.analinfo
```

These `.analinfo` files are plain text and are the first place to look when
checking whether a trajectory or initial condition is sensible.

## What The Dynamics Step Means

The bundled dynamics scripts are deliberately cheap. They demonstrate how ICATS
coordinate/velocity files can be passed to a dynamics engine and then read back
for analysis.

Do not judge the initial-condition generator by long-time energy drift in the
toy dynamics alone. First use the initial-sample audit to check that the
generated model state is internally consistent before dynamics.
