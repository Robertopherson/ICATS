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

## Minimal Shape Of An Input File

Every scattering input needs two molecule entries:

```text
mol = 0 ammonia_dat.txt
mol = 1 h2o_dat.txt
```

The `0` and `1` labels identify the two collision partners. The molecule files
point to molecular data files, which in turn reference geometries, Hessians, and
frequencies when needed.

## Run Directory Behaviour

ICATS writes generated files under:

```text
rd_<run-tag>/
```

For example:

```text
run-tag = tutorial_input
```

creates:

```text
rd_tutorial_input/
```

This directory contains restart files, Wang-Landau umbrellas, histograms, and
outputs. If two calculations use the same `run-tag`, they use the same run
directory. Choose a new `run-tag` when comparing different setups.

## Common Keys

| Key | Meaning | Typical use |
| --- | --- | --- |
| `Nsamp` | Number of initial-condition samples | Increase for production statistics |
| `workers` | Number of parallel workers | Use a small number first |
| `seed` | Random seed | Use fixed seed for reproducible tests |
| `Trot` | Rotational temperature | Controls rotor state distribution |
| `Tvib` | Vibrational temperature | Controls vibrational state distribution |
| `Rz` | Initial separation along the collision axis | Should place molecules far apart |
| `maxl` | Maximum orbital angular momentum | Often more important than `maxb` |
| `maxb` | Maximum impact parameter estimate | Useful for physical intuition |
| `wang` | Enable Wang-Landau weighting | Use for difficult `J` distributions |
| `plothist` | Generate histogram scripts | Use for diagnostics |

## Output Controls

The common tutorial setting:

```text
printout = 0 0 1 0
```

requests Cartesian coordinate/velocity output suitable for dynamics. To write a
full generation log for debugging, use:

```text
printout = 0 1 0 0
```

This creates `out_full.info`, which is useful for reading the generation and
immediate analysis blocks for each sample.

## Safe First Values

For a quick test:

```text
Nsamp = 2
workers = 1
progress = quiet
```

For a small tutorial run:

```text
Nsamp = 10
workers = 2
```

Increase `Nsamp` only after the audit and histograms look sensible.
