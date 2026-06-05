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
orbital-sampling = geometric
phisample = True
```

Useful constrained-sampling controls:

```text
vib-mode = rigid
fixed-b = 3.5
impact-phi = 0.0
wang = False
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
| `vel` | `vel = 1000 100 3` | Molecule beam-speed distribution: centre speed in m/s, FWHM in m/s, and velocity-power weight. |
| `nfreeze` | `nfreeze = 1 2 3` | Normal-mode indices to leave at zero during vibrational sampling. |
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

To write starter files from the command line, use `icats.init
--write-templates`. Add `--overwrite-templates` only when replacing existing
template files intentionally.

When `wang = True`, the first run may spend most of its time building
`rd_<run-tag>/wang.pkl`. This Wang-Landau umbrella build can take many minutes,
and for larger polyatomic/high-`maxj` cases it can take an hour or more.
Increasing `workers` can help, but use `workers = 1` while debugging input
files and increase to a few cores only after the setup is known to be correct.
During production sampling, Wang-Landau acceptance/rejection counters are
appended to `<input-file>.logfile` so that long runs can be monitored with
`tail -f`.

### Temperatures and Velocity

| Key | Form | Meaning |
| --- | --- | --- |
| `Trot` | `Trot = 50.0` | System rotational temperature. |
| `Tvib` | `Tvib = 50.0` | System vibrational temperature. |
| `Tvel` | `Tvel = 300.0` | Intermolecular velocity temperature. |
| `Tvel` | `Tvel = -500.0 50.0` | Legacy direct relative-speed form: centre and optional FWHM in m/s. |
| `relative-velocity` | `relative-velocity = 1000.0` | Fixed direct relative speed in m/s. |
| `relative-velocity` | `relative-velocity = 1000.0 50.0` | Direct relative-speed distribution with centre and FWHM in m/s. |
| `relative-velocity-fwhm` | `relative-velocity-fwhm = 50.0` | Optional FWHM for a direct relative-speed input, in m/s. |
| `collision-energy` | `collision-energy = 0.050` | Fixed direct collision energy in eV. |
| `incoming-p0` | `incoming-p0 = 12.0` | Fixed direct relative momentum in atomic units. |
| `incoming-k` | `incoming-k = 12.0` | Fixed direct incoming wave number in bohr^-1. In atomic units this is numerically the same as `p0`. |

The top-level `Trot` and `Tvib` are copied into molecules unless molecule-level
values override them. Atom-only systems force both to zero.

There are two physically different ways to set the incoming translational
motion.

**Crossed-beam mode** uses the molecule-file `vel` entries. Each beam speed is
sampled independently, and `beam-angle` determines the relative speed:

```text
# molecule A file
vel = 600.0 100.0 3

# molecule B file
vel = 800.0 100.0 3

# top-level scattering input
beam-angle = 90.0
```

This is the closest model of a crossed-beam experiment. The third number in
`vel` is the velocity-power weight; `3` corresponds to the flux-weighted form
often used for effusive or supersonic molecular beams.

**Direct relative-channel mode** bypasses separate beam speeds and sets the
Jacobi relative speed directly. Use this for idealised fixed-energy scattering
or when thinking in partial-wave variables:

```text
relative-velocity = 1000.0
```

or equivalently:

```text
collision-energy = 0.050
```

or:

```text
incoming-k = 12.0
```

For direct-channel inputs ICATS converts internally between

```text
E_coll = 1/2 * mu * v_rel^2
p0 = mu * v_rel
k = p0 / hbar
```

where `hbar = 1` in atomic units. Do not combine molecule-file `vel`
distributions with direct-channel inputs in the same production run unless you
are intentionally testing precedence; the top-level direct-channel input takes
over the intermolecular relative speed.

### Intermolecular Geometry and Angular Momentum

| Key | Form | Meaning |
| --- | --- | --- |
| `Rz` | `Rz = 15.0` | Initial separation along the collision axis, in Angstrom. |
| `beam-angle` | `beam-angle = 90.0` | Crossed-beam angle in degrees. |
| `maxb` | `maxb = 10.0` | Maximum impact parameter estimate, in Angstrom. |
| `fixed-b` | `fixed-b = 3.5` | Fixed impact parameter in Angstrom. This bypasses impact-parameter sampling and derives the corresponding orbital angular momentum from the sampled relative velocity. Currently use with `wang = False`. |
| `impact-phi` | `impact-phi = 0.0` | Fixed impact-parameter azimuth in the lab frame, in radians. Omit this key to sample the azimuth. This controls the collision plane, not the full Jacobi/Euler `phi, beta, chi` transformation reconstructed later from the Cartesian geometry. |
| `maxl` | `maxl = 80` | Maximum orbital angular momentum quantum number. |
| `maxj` | `maxj = 80` | Maximum total angular momentum for Wang-Landau/J setup. If `maxl` is omitted, ICATS can use `maxj` as the orbital cap. |
| `orbital-sampling` | `geometric` or `flat-l` | Proposal distribution for non-fixed orbital angular momentum. `geometric` is the default/current behavior and samples the classical impact-parameter measure, approximately `P(L) ~ L`. `flat-l` samples `L` uniformly over the requested cap and lets `J = L + Jab` emerge. |
| `chi` | `chi = 0.0` | Legacy azimuthal scattering-angle/control variable. Do not use this to fix the lab impact plane; use `impact-phi` for that. |
| `phisample` | `True` or `False` | Legacy standard-azimuth switch. In isotropic runs, `True` rotates the sampled system to a standard azimuth unless `impact-phi` is explicitly set. |

For non-Wang-Landau runs, `maxl` is usually the direct cap. For Wang-Landau
runs, `maxj` is the important requested total-`J` range, but `L` still enters
the actual Cartesian collision geometry.

Use `orbital-sampling = flat-l` for diagnostic or trajectory-budgeting runs
where the small-impact-parameter region should be sampled as heavily as the
large-impact-parameter region. The resulting ensemble is not geometrically
weighted by default; reweight by the desired measure, for example by `L` or
`b`, when reconstructing a geometric cross-section-like average.

#### Fixed Impact-Parameter Runs

For diagnostic calculations, it is often useful to hold the collision geometry
fixed while still sampling internal rotations and orientations. For example, an
atom-diatom setup such as Ar + NO can be run with a fixed impact parameter,
fixed lab-frame impact plane, and rigid NO bond:

```text
mol = 0 ar_dat.txt
mol = 1 no_dat.txt

Trot = 300.0
Tvib = 0.0
vib-mode = rigid

fixed-b = 3.5
impact-phi = 0.0
wang = False
```

In this mode ICATS samples the relative velocity first, then sets the orbital
angular momentum from

```text
|L| = mu * v_rel * b
```

and uses `impact-phi` as the lab-frame azimuth of the impact parameter. If
`impact-phi` is absent, the azimuth is sampled. The later Jacobi/Euler
`phi, beta, chi` values are still reconstructed from the generated Cartesian
geometry; they are not direct input controls for the full body-frame
transformation. Older sample-data keys may call the polar angle `theta`; in the
manual and current logs this same angle is labelled `beta`.

The packaged tutorial `fixed_plane_atom_diatom_ar_no` is the runnable version
of this idea. It uses `Trot = 0.0`, `vib-mode = rigid`, `incoming-k`,
`fixed-b`, and `impact-phi` to create a fixed-plane Ar + NO ensemble while
still sampling the NO orientation.

### Rotations and Orientations

| Key | Form | Meaning |
| --- | --- | --- |
| `ordist` | `ordist = ...` | System-level orientation distribution and parameters. |
| `rot-param` | `rot-param = xyz` or `euler` | System-level rotation parameterisation; overrides molecule-level values. |

Use isotropic/default orientations first. When using custom orientation
distributions, enable histograms and inspect the resulting angular distributions.

Set `Trot = 0.0` when the molecule should have no initial rigid-rotor angular
momentum but should still have a sampled orientation. This is useful for
fixed-plane atom-diatom tests where the diatom bond direction is sampled but
the diatom starts non-rotating.

For molecular orientations, ICATS uses the manuscript Euler labels
`alpha, beta, gamma`. Linear molecules have an arbitrary final spin about the
bond axis, so ICATS reports `gamma = 0` for that gauge coordinate; the physical
axis direction is controlled by `alpha` and `beta`.

### Vibrations

| Key | Form | Meaning |
| --- | --- | --- |
| `maxv` | `maxv = 10` | Maximum vibrational quantum number considered in Boltzmann sampling. |
| `Tvib` | `Tvib = 500.0` | Vibrational temperature, also listed above because it controls vibrational populations. |
| `vib-mode` | `vib-mode = sample` or `rigid` | `sample` uses the harmonic normal-mode sampler. `rigid` skips vibrational sampling and leaves molecules at the reference geometry while still sampling rotations and orientations. |

The normal modes come from the molecule Hessian. If a molecule has no Hessian,
the current model can still handle rigid/atomic pieces but has no harmonic
vibrational mode sampling for that molecule.

`vib-mode = rigid` is the simplest way to freeze all intramolecular vibrational
motion for a run without editing the molecule file. Molecule-level `nfreeze`
can be used for finer control when only selected normal modes should be kept at
zero. `maxv = 0` is different: it still samples the ground-state Husimi/Wigner
width, so it is not a rigid-bond setting.

### Wang-Landau

| Key | Form | Meaning |
| --- | --- | --- |
| `wang` | `True` or `False` | Enables Wang-Landau rejection weighting. |
| `wlmode` | `fast`, `default`, `accurate` | Preset for Wang-Landau convergence parameters. |
| `wl-target` | `linear-j` or `flat-j` | Target distribution in total `J` after dividing by the WL density estimate. If omitted, ICATS uses `linear-j` with `orbital-sampling = geometric` and `flat-j` with `orbital-sampling = flat-l`. Mixed combinations are rejected. |
| `wl-ff` | `wl-ff = 1.05` | Initial modification factor override. Must be greater than 1. |
| `wl-nstep` | `wl-nstep = 500` | Number of steps per WL bin multiplier. |
| `wl-flatness` | `wl-flatness = 0.90` | Histogram flatness criterion. Larger is stricter. |
| `wl-wn-factor` | `wl-wn-factor = 4.0` | Sets WL bins from `PeakJab * factor`. |
| `wl-j-range` | `wl-j-range = 60` | Upper `J` range for the explicit WL density estimate. |
| `wl-j-bins` | `wl-j-bins = 80` | Directly sets the number of WL bins. |
| `wl-l-cap` | `wl-l-cap = 60` | Orbital `L` proposal cap during WL construction. |
| `wl-wn` | `wl-wn = 80` | Old alias for `wl-j-bins`. |
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
| `printout` | `printout = 0 0 1 0` | Four output flags: combined xyz/vel files, full info file, xyz/vel directory, info directory. |
| `output-format` | `xyzvel`, `npz`, `both` | Coordinate/velocity output format. |
| `units-out` | `ang-fs`, `au` | Output unit system. |
| `keepinfo` | `True` or `False` | Store extra per-sample information. |

The `printout` line has four integer switches:

| Position | Example value | Files written |
| --- | --- | --- |
| 1 | `1 0 0 0` | Combined coordinate/velocity files: `out_full.xyz`, `out_full.vel`. |
| 2 | `0 1 0 0` | Combined generation and analysis log: `out_full.info`. |
| 3 | `0 0 1 0` | Per-sample files in `dirout`: `outputs/out_0.xyz`, `outputs/out_0.vel`, ... |
| 4 | `0 0 0 1` | Per-sample info logs in `dirout`: `outputs/out_0.info`, ... |

The common tutorial setting `printout = 0 0 1 0` writes a directory of Cartesian
coordinate and velocity files suitable for the demonstration dynamics. The
combined-file setting `printout = 1 0 0 0` writes the same samples into
`out_full.xyz` and `out_full.vel`.

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
| `audit-initial-vib-tol` | `audit-initial-vib-tol = 0.0` | Normal-coordinate `Q/P` RMS-per-component tolerance; `0` records diagnostics without failing. |
| `audit-initial-velocity-tol` | `audit-initial-velocity-tol = 0.0` | Relative-velocity tolerance in m/s. |

For development and tutorials, `audit-initial-sample = True` is often the most
important safety check. It compares the generated bookkeeping with the immediate
Cartesian-coordinate analysis before any dynamics are run. For production
ensembles, histograms are the quickest way to catch accidental boundary effects
or inappropriate `J`/`L` ranges.

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

The `printout` setting is the main control over how generated samples are
written. It is a four-number switch:

```text
printout = combined_xyzvel combined_info per_sample_xyzvel per_sample_info
```

Each number is either `0` or `1`. Several switches can be enabled at once, but
large runs should avoid writing unnecessary logs.

### Dynamics-Friendly Output

The common tutorial setting is:

```text
printout = 0 0 1 0
```

This writes one coordinate file and one velocity file per sampled initial
condition:

```text
rd_<run-tag>/outputs/out_0.xyz
rd_<run-tag>/outputs/out_0.vel
rd_<run-tag>/outputs/out_1.xyz
rd_<run-tag>/outputs/out_1.vel
```

This is usually the easiest output for trajectory scripts, because each sample
is already a separate coordinate/velocity pair.

### Single Combined Coordinate/Velocity Pair

To write all initial conditions into one coordinate file and one matching
velocity file, use:

```text
printout = 1 0 0 0
```

This writes:

```text
out_full.xyz
out_full.vel
```

The two files contain matching sample blocks in the same order.

### Debugging A Sample

To write the detailed generation and immediate analysis log, use:

```text
printout = 0 1 0 0
```

This creates `out_full.info`, which is useful for reading the generation and
immediate analysis blocks for each sample.

For a very small debugging run, it is often useful to write both the Cartesian
files and the full log:

```text
Nsamp = 2
printout = 1 1 1 0
```

For large runs, avoid `printout = 0 1 0 0` or `printout = 0 0 0 1` unless you
really need the text logs; they can become bulky.

### Histogram-Only Or Setup Checks

If you only want histograms or input checking, it can be useful to suppress
coordinate/velocity output:

```text
printout = 0 0 0 0
plothist = True
hist_initial = True
hist_sampled = True
```

This is the pattern used by some histogram-focused tutorials.

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
