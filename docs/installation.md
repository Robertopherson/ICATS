---
layout: default
title: Installation
---

# Installation

## Requirements

ICATS is a Python package. A new user should start from a clean Python
environment so that numerical dependencies do not conflict with packages already
installed for other projects.

One suitable workflow is:

```bash
conda create -n icats python=3.11
conda activate icats
```

Other environment managers are fine, provided they supply a recent Python and
can install NumPy/SciPy-based packages.

## Install from GitHub

Clone the repository and install ICATS in editable mode:

```bash
git clone https://github.com/Robertopherson/ICATS.git
cd ICATS
python -m pip install -e .
```

Editable mode is useful while the package is still under active development,
because updates to the source tree are immediately visible to the command-line
tools.

## Optional Dynamics Dependencies

The initial-condition generator and analysis tools are the core of ICATS. The
cheap demonstration dynamics use PySCF semiempirical functionality, which may
need additional packages:

```bash
python -m pip install -e ".[dynamics]"
```

This extra currently installs:

```text
pyscf
pyscf-semiempirical
```

If the dynamics dependencies are not available, users can still generate initial
conditions, inspect histograms, and analyse compatible trajectory files produced
elsewhere.

## Command-Line Tools

The install provides:

```bash
icats
icats.init
icats.analyse
icats.audit-tutorials
```

## Checking the Install

```bash
icats --list-tutorials
icats.init --show-input-options
```

For a more complete check of the initial-condition pipeline, run:

```bash
icats.audit-tutorials --nsamp 3 --keep-going
```

This generates each bundled tutorial and verifies that the generated sample
energies and comparable state variables are recovered by the immediate
Cartesian-coordinate analysis. It does not run the cheap dynamics backend, and
it keeps Wang-Landau disabled unless `--include-wl` is requested.

## Cache Directories

The command-line entry points set temporary default cache directories for numba
and matplotlib when the user has not already configured them. This helps avoid
cache-write problems on cluster filesystems and in restricted environments.
