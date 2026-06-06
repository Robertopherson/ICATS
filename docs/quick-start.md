---
layout: default
title: Quick Start
---

# Quick Start

This page walks through the smallest useful workflow. It assumes ICATS has
already been installed; see [Installation](installation.md) first if the
commands are not available.

## 1. See Available Tutorials

```bash
icats --list-tutorials
```

Expected output is a numbered list of examples, including `quickstart`,
`single_atom_he_he`, and Wang-Landau tutorials.

## 2. Generate The Default Tutorial

```bash
icats --tutorial quickstart --setup-only
cd tutorial_quickstart
```

This creates a self-contained tutorial directory with:

```text
tutorial_input.txt
run_cheap_dynamics.sh
run_analysis.sh
*_dat.txt
*_geom.xyz
*_hessian.txt
*_freq.txt
```

The file to inspect first is `tutorial_input.txt`. It contains the two molecule
definitions and the sampling options.

## 3. Generate Initial Conditions

```bash
icats.init tutorial_input.txt
```

The most important generated directory is:

```text
rd_tutorial_input/
```

For the default tutorial, the generated Cartesian samples appear under:

```text
rd_tutorial_input/outputs/
```

The tutorial input normally uses:

```text
printout = 0 0 1 0
```

which writes one coordinate/velocity pair per sample:

```text
rd_tutorial_input/outputs/out_0.xyz
rd_tutorial_input/outputs/out_0.vel
```

If a downstream code prefers one combined coordinate file and one combined
velocity file, use:

```text
printout = 1 0 0 0
```

which writes `out_full.xyz` and `out_full.vel`.

The initial-condition stage is independent of the cheap dynamics backend.

## 4. Run The Initial-Condition Audit

For a quick check of all bundled tutorials:

```bash
icats.audit-tutorials --nsamp 3 --keep-going
```

To repeat the audit in the `incoming-k-plus-z` scattering convention:

```bash
icats.audit-tutorials --nsamp 2 --keep-going --output-frame incoming-k-plus-z
```

For a focused frame-convention regression check:

```bash
icats.frame-smoke
```

For a single tutorial, add the following lines to `tutorial_input.txt` before
running `icats.init`:

```text
audit-initial-sample = True
audit-initial-energy-tol = 0.02
audit-initial-angular-tol = 2.0
audit-initial-vib-tol = 2.0
audit-initial-velocity-tol = 5.0
```

The audit verifies that sampled model energies, angular momenta, impact
parameter, relative velocity, and defined normal coordinates are recovered from
the immediate Cartesian-coordinate analysis at `t = 0`. The tutorial sweep keeps
Wang-Landau off by default, so it is meant as a quick regression check after
code changes rather than as a full umbrella-convergence test.

To understand the first generated log, read
[Annotated First Output](annotated-output.md) with `out_full.info` or
`dynamics*.analinfo` open.

## 5. Run The Cheap Dynamics Demonstration

The tutorial can now pass the generated samples to the simple demonstration
dynamics script:

```bash
./run_cheap_dynamics.sh
```

This step is optional. It exists to teach the file pipeline; serious production
work may use a different dynamics code.

## 6. Analyse Generated Trajectories

```bash
./run_analysis.sh
```

Important outputs include:

```text
rd_tutorial_input/outputs/out_*.md.xyz
rd_tutorial_input/outputs/out_*.md.vel
rd_tutorial_input/outputs/dynamics*.analinfo
```

The `dynamics*.analinfo` files contain orientation, rotational, vibrational,
intermolecular, and energy-decomposition diagnostics for each saved frame.

## What Success Looks Like

A successful first run should produce:

- no Python traceback,
- `rd_tutorial_input/dat_tutorial_input.txt.pkl`,
- `rd_tutorial_input/work_sys_tutorial_input.txt.pkl`,
- one or more `out_*.md.xyz` and `out_*.md.vel` files after dynamics,
- one or more `dynamics*.analinfo` files after analysis.

If the dynamics step fails but `icats.init` succeeds, the initial-condition
generator is probably working. Check [Troubleshooting](troubleshooting.md) for
PySCF and environment notes.
