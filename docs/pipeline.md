# Dynamics and Analysis Pipeline

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

For a quick regression check of the tutorial initial-condition pipeline:

```bash
icats.audit-tutorials
```

This command does not call the cheap dynamics backend. It only generates each
tutorial, runs `icats.init`, and checks that sampled model energies are
recovered from the immediate Cartesian analysis.
