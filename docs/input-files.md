---
layout: default
title: Input Files
---

# Input Files

ICATS input files are plain text key-value files. Show the recognized options:

```bash
icats.init --show-input-options
```

Keys are parsed case-insensitively. The examples use `Nsamp`, `Trot`, `Tvib`,
and `Rz` because that is the older human-readable style, but internally these
are read as lowercase keys.

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

## Option Reference

The tables below group the current input keys by what they control. They are
intended as a readable reference, not as a replacement for the tutorials.

### Molecule Selection

These options appear in the top-level scattering input file.

| Key | Form | Meaning |
| --- | --- | --- |
| `mol` | `mol = 0 file.txt` and `mol = 1 file.txt` | Required. Points to the two molecule definition files. |

The molecule index must be `0` or `1`. File paths are resolved relative to the
top-level input file.

### Molecule Definition Files

These options appear inside each molecule file referenced by `mol`.

| Key | Form | Meaning |
| --- | --- | --- |
| `name` | `name = ammonia` | Optional label used in logs. |
| `xyz` | `xyz = geom.xyz` | Reference geometry in Angstrom. Required for normal molecule setup. |
| `hess` | `hess = hessian.txt` | Mass-weighted or code-compatible Hessian used to construct normal modes. |
| `w` | `w = freq.txt` | Frequency-file reference used by templates/validation; current normal-mode sampling is driven by `hess`. |
| `trot` | `trot = 50.0` | Molecule-level rotational temperature override. |
| `tvib` | `tvib = 50.0` | Molecule-level vibrational temperature override. |
| `vel` | `vel = 1000 100 3` | Molecule beam-speed distribution: centre, FWHM, and velocity power. |
| `nfreeze` | `nfreeze = 1 2 3` | Atom indices to freeze/exclude from normal-mode treatment. |
| `ordist` | `ordist = read file function ...` | Read a molecule orientation distribution from a user-supplied function. |
| `ordist` | `ordist = fixed` | Fixed-orientation mode. |
| `rot-param` | `rot-param = xyz` or `euler` | Rotation-angle parameterisation for this molecule. |

For atoms, ICATS forces molecular `Trot` and `Tvib` to zero because there are no
internal rotor or vibrational degrees of freedom.

### Ensemble Size and Runtime

| Key | Form | Meaning |
| --- | --- | --- |
| `Nsamp` | `Nsamp = 1000` | Number of generated initial-condition samples. |
| `workers` | `workers = 4` | Number of parallel workers. Start with `1` or `2` while debugging. |
| `seed` | `seed = 400` | Base random seed. |
| `seed-mode` | `fixed`, `time`, `per-worker` | How the seed is chosen or distributed. |
| `progress` | `quiet`, `normal`, `verbose` | Console verbosity. |
| `save-frequency` | `save-frequency = 100` | Worker checkpoint frequency; `0` disables periodic checkpointing. |
| `continue` | `True` or `False` | Legacy restart flag using existing worker pickle files. |
| `run-mode` | `fresh`, `continue`, `rebuild-wang` | Higher-level run behaviour. `rebuild-wang` refuses to overwrite an existing `wang.pkl`. |
| `run-tag` | `run-tag = my_run` | Names the run directory `rd_<run-tag>/`. |
| `dry-run` | `True` or `False` | Validate/load setup and skip sample generation. |
| `check-input` | `True` or `False` | Check input/setup only and exit before sampling. |

For normal use, prefer `run-mode` over the older `continue` flag.

When `wang = True`, the first run may spend most of its time building
`rd_<run-tag>/wang.pkl`. This Wang-Landau umbrella build can take many minutes,
and for larger polyatomic/high-`maxj` cases it can take an hour or more.
Increasing `workers` can help, but use `workers = 1` while debugging input
files and increase to a few cores only after the setup is known to be correct.

### Temperatures and Velocity

| Key | Form | Meaning |
| --- | --- | --- |
| `Trot` | `Trot = 50.0` | System rotational temperature. |
| `Tvib` | `Tvib = 50.0` | System vibrational temperature. |
| `Tvel` | `Tvel = 300.0` | Intermolecular velocity temperature. |
| `Tvel` | `Tvel = -500.0 50.0` | Explicit beam-centred velocity form: centre and optional FWHM in m/s. |

The top-level `Trot` and `Tvib` are copied into molecules unless molecule-level
values override them. Atom-only systems force both to zero.

### Intermolecular Geometry and Angular Momentum

| Key | Form | Meaning |
| --- | --- | --- |
| `Rz` | `Rz = 15.0` | Initial separation along the collision axis, in Angstrom. |
| `beam-angle` | `beam-angle = 90.0` | Crossed-beam angle in degrees. |
| `maxb` | `maxb = 10.0` | Maximum impact parameter estimate, in Angstrom. |
| `maxl` | `maxl = 80` | Maximum orbital angular momentum quantum number. |
| `maxj` | `maxj = 80` | Maximum total angular momentum for Wang-Landau/J setup. If `maxl` is omitted, ICATS can use `maxj` as the orbital cap. |
| `chi` | `chi = 0.0` | Azimuthal scattering angle/control used by the intermolecular setup. |
| `phisample` | `True` or `False` | Whether to sample the orbital azimuthal coordinate `phi`. |

For non-Wang-Landau runs, `maxl` is usually the direct cap. For Wang-Landau
runs, `maxj` is the important requested total-`J` range, but `L` still enters
the actual Cartesian collision geometry.

### Rotations and Orientations

| Key | Form | Meaning |
| --- | --- | --- |
| `ordist` | `ordist = ...` | System-level orientation distribution and parameters. |
| `rot-param` | `rot-param = xyz` or `euler` | System-level rotation parameterisation; overrides molecule-level values. |

Use isotropic/default orientations first. When using custom orientation
distributions, enable histograms and inspect the resulting angular distributions.

### Vibrations

| Key | Form | Meaning |
| --- | --- | --- |
| `maxv` | `maxv = 10` | Maximum vibrational quantum number considered in Boltzmann sampling. |
| `Tvib` | `Tvib = 500.0` | Vibrational temperature, also listed above because it controls vibrational populations. |

The normal modes come from the molecule Hessian. If a molecule has no Hessian,
the current model can still handle rigid/atomic pieces but has no harmonic
vibrational mode sampling for that molecule.

### Wang-Landau

| Key | Form | Meaning |
| --- | --- | --- |
| `wang` | `True` or `False` | Enables Wang-Landau rejection weighting. |
| `wlmode` | `fast`, `default`, `accurate` | Preset for Wang-Landau convergence parameters. |
| `wl-ff` | `wl-ff = 1.05` | Initial modification factor override. Must be greater than 1. |
| `wl-nstep` | `wl-nstep = 500` | Number of steps per WL bin multiplier. |
| `wl-flatness` | `wl-flatness = 0.90` | Histogram flatness criterion. Larger is stricter. |
| `wl-wn-factor` | `wl-wn-factor = 4.0` | Sets WL bins from `PeakJab * factor`. |
| `wl-wn` | `wl-wn = 80` | Directly sets the number of WL bins. |
| `wl-tol` | `wl-tol = 1.000001` | Stopping tolerance for the modification factor. |
| `wl-max-iter` | `wl-max-iter = 0` | Maximum WL iterations; `0` means no explicit cap. |
| `wl-log-every` | `wl-log-every = 1` | WL status print period. |
| `wl-angular-sampler` | `fast` or `legacy` | Numerical implementation used inside WL trial sampling. |
| `wl-audit-angular-sampler` | `True` or `False` | Compares fast and legacy angular samplers for debugging. |

The stored umbrella lives at `rd_<run-tag>/wang.pkl`. ICATS checks metadata
before reuse and refuses incompatible files.

### Output Files

| Key | Form | Meaning |
| --- | --- | --- |
| `fileout` | `fileout = out` | Output filename prefix. |
| `dirout` | `dirout = outputs` | Output directory. Relative paths are placed inside the run directory. |
| `printout` | `printout = 0 0 1 0` | Four output flags: single xyz, full info file, xyz/vel directory, info directory. |
| `output-format` | `xyzvel`, `npz`, `both` | Coordinate/velocity output format. |
| `units-out` | `ang-fs`, `au` | Output unit system. |
| `keepinfo` | `True` or `False` | Store extra per-sample information. |

The common tutorial setting `printout = 0 0 1 0` writes a directory of Cartesian
coordinate and velocity files suitable for the demonstration dynamics.

### Histograms and Diagnostics

| Key | Form | Meaning |
| --- | --- | --- |
| `plothist` | `True` or `False` | Convenience switch for sampled histograms. |
| `hist_initial` | `True` or `False` | Generate initial/pre-sampling histogram helpers. |
| `hist_sampled` | `True` or `False` | Generate sampled-output histogram helpers. |
| `plotinit` | `plotinit = 1000` | Number of pre-samples used for initial distribution plots. |
| `audit-initial-sample` | `True` or `False` | Run generation-vs-analysis audit at `t = 0`. |
| `audit-initial-energy-tol` | `audit-initial-energy-tol = 0.02` | Energy tolerance in eV for the initial audit. |
| `audit-initial-angular-tol` | `audit-initial-angular-tol = 0.0` | Angular-momentum tolerance; `0` disables angular pass/fail checks. |

For development and tutorials, `audit-initial-sample = True` is often the most
important safety check. For production ensembles, histograms are the quickest
way to catch accidental boundary effects or inappropriate `J`/`L` ranges.

## Choosing `maxb`, `maxl`, and `J`

The impact parameter `b` is the geometrical collision variable. A larger
`maxb` includes more glancing trajectories, but those trajectories may not
interact strongly if the molecules are already outside the useful range of the
potential. A very small `maxb` can make a tutorial run quickly while hiding the
large-`L` part of the collision ensemble.

The orbital angular momentum `L` is connected to `b` through the relative
momentum. For atom-atom scattering this is almost the whole angular-momentum
story. For molecules, the total angular momentum `J` also includes the rotor
sum:

```text
J = L + J_A + J_B
```

This means a sensible `maxb` or `maxl` does not by itself guarantee a sensible
total-`J` distribution. When `wang = True`, always inspect the sampled `J` and
`L` histograms after the run.

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
