---
layout: default
title: Installation
---

# Installation

## Development Install

From the root of the ICATS repository:

```bash
python -m pip install -e .
```

This installs the command-line tools:

```bash
icats
icats.init
icats.analyse
```

For the current local development copy, the source directory is:

```bash
/home/chris/Programs/ICATS
```

## Environment Notes

The cheap dynamics examples use PySCF semiempirical functionality. On the local
test machine, the clean environment was:

```bash
conda activate icats_clean
```

If ICATS is not installed into that environment, it can still be used during
development with:

```bash
export PATH=/home/chris/Programs/ICATS/scripts:$PATH
export PYTHONPATH=/home/chris/Programs/ICATS
```

The command-line entry points set temporary default cache directories for numba
and matplotlib when the user has not already configured them. This helps avoid
cache-write problems on cluster filesystems and in restricted environments.

## Checking the Install

```bash
icats --list-tutorials
icats.init --show-input-options
```
