---
layout: default
title: Input Files
---

# Input Files

ICATS input files are plain text key-value files. Show the recognized options:

```bash
icats.init --show-input-options
```

Core entries include:

```text
mol = 0 molecule_a_dat.txt
mol = 1 molecule_b_dat.txt
Nsamp = 100
workers = 4
seed = 400
run-mode = fresh
run-tag = tutorial_input
```

Common sampling controls:

```text
Tvib = 500.0
Trot = 500.0
Rz = 15
maxl = 182
maxb = 16
phisample = True
```

Common output controls:

```text
fileout = out
dirout = outputs
printout = 0 0 1 0
output-format = xyzvel
units-out = ang-fs
```

`run-tag` controls the generated run directory name:

```text
rd_<run-tag>/
```

When `run-tag` is absent, ICATS derives the run directory from the input file
stem.
