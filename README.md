# ICATS

ICATS generates initial conditions and analysis diagnostics for quasi-classical
molecular scattering calculations.

This repository contains the code package, small example molecule files, and
tutorial-generation helpers. Generated tutorial outputs and run directories are
not intended to be committed.

## Install for Development

From this directory:

```bash
python -m pip install -e .
```

This installs the command-line tools:

```bash
icats
icats.init
icats.analyse
```

## Basic Use

Generate initial conditions from an input file:

```bash
icats.init tutorial_input.txt
```

List and generate tutorials:

```bash
icats --list-tutorials
icats --tutorial quickstart --setup-only
```

Analyse trajectory output pairs:

```bash
icats.analyse tutorial_input.txt --dir rd_tutorial_input/outputs --prefix out
```

Run the tutorial initial-condition audit without dynamics:

```bash
icats.audit-tutorials
```

## Manual

The user manual is being developed in [`docs/`](docs/index.md). It covers
installation, the tutorial workflow, input files, Wang-Landau umbrellas,
diagnostics, and the theory behind the initial-condition sampling.

## Wang-Landau Umbrellas

For Wang-Landau runs, `wang.pkl` stores the generated Wang-Landau umbrella and
metadata describing the input settings used to create it. If the input settings
change, move or rename the existing `wang.pkl` before rerunning. ICATS refuses
incompatible stored umbrellas rather than silently reusing them.

## Notes

The CLI entry points set default temporary locations for numba and matplotlib
cache directories when the user has not already configured them. This helps on
cluster filesystems, conda environments, and restricted workspaces.
