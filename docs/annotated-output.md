---
layout: default
title: Annotated First Output
---

# Annotated First Output

This page explains how to read `out_full.info` or a per-sample `.info` file.
The `.info` file is deliberately split into two matched halves:

```text
===================================================================
Sample 0 | generation
===================================================================

...

===================================================================
Sample 0 | analysis
===================================================================
```

`generation` is the model-side draw: the vibrational, rotational,
orientational, velocity, and impact-parameter values ICATS intended to sample.
`analysis` is the round-trip check: ICATS reconstructs the same kinds of
quantities from the Cartesian coordinates and velocities it just wrote. The
sections are ordered and named so that corresponding information can be
compared directly.

## Output Conventions

Scalar lines use:

```text
label                  =          value  unit
```

Vector lines use:

```text
label                  = [          x,          y,          z ]  unit, |name| =      value
```

This spacing is intentional. The commas prevent negative signs from merging
with the previous number, and the explicit unit keeps text parsing simple.

Angles in Euler triples are printed in units of `pi rad`:

```text
system Euler           = [    -0.0000,     0.0246,     0.0000 ]  pi rad   # phi, beta, chi
alpha,beta,gamma       = [     0.6287,     0.3324,     0.0000 ]  pi rad   # molecular Euler
```

## 1. Generation Header

A generated sample begins with:

```text
===================================================================
Sample 0 | generation
===================================================================
```

Everything under this header is sampled-model information. It is not a later
dynamics result and it is not read from a potential-energy surface.

## 2. Intermolecular Setup

The intermolecular block gives the relative incoming speed, collision energy,
and later the angular-momentum geometry:

```text
[intermolecular]
relative velocity      =    953.5976463  m/s
collision energy       =      0.0807475  eV
ar_dat z velocity      =    409.0352371  m/s
no_dat z velocity      =   -544.5624092  m/s
ar_dat kinetic         =      0.0346357  eV
no_dat kinetic         =      0.0461117  eV

fixed b                =        4.50000  Ang
implied L              =         115.28
total J                =         115.28
```

For direct-channel runs, `relative velocity`, `collision energy`, or
`incoming-k` determines the incoming two-body channel. For crossed-beam runs,
the molecule-level beam velocities are first converted into this relative
channel.

## 3. Orientation

For a diatom with sampled orientation and zero rotational angular momentum:

```text
[orientation]
ar_dat                 = no orientational state
no_dat                 = sampled orientation
alpha,beta,gamma       = [  0.6287488,  0.3323750,  0.0000000 ]  pi rad   # molecular Euler
wx,wy,wz               = [ -0.3205039, -0.1105013,  0.0000000 ]  pi rad   # XYZ rotation angles
```

For a linear molecule, `gamma` is an arbitrary spin about the molecular axis.
ICATS reports it as zero. The physical bond direction is carried by
`alpha,beta`.

## 4. Rotation

Generation shows the sampled rotor model. Analysis reconstructs several
rotational diagnostics from the Cartesian coordinates:

```text
[rotation]
no_dat                 = linear, symmetry constant =  -1.00
full J, space          = [      0.00000,      0.00000,      0.00000 ]  au, |J| =      0.00000  (QM:   0.00)
full J, Eckart         = [      0.00000,      0.00000,      0.00000 ]  au, |J| =      0.00000  (QM:   0.00)
vector J, space        = [      0.00000,      0.00000,      0.00000 ]  au, |J| =      0.00000  (QM:   0.00)
vector J, Eckart       = [      0.00000,      0.00000,      0.00000 ]  au, |J| =      0.00000  (QM:   0.00)
vibrational J          = [      0.00000,      0.00000,      0.00000 ]  au, |J| =      0.00000  (QM:   0.00)
vector rot. energy     = [      0.00000,      0.00000,      0.00000 ]  eV, |E| =      0.00000
full rot. energy       = [      0.00000,      0.00000,      0.00000 ]  eV, |E| =      0.00000
```

The most important comparison for sampler validation is the generated
rigid-rotor energy against `vector rot. energy`. `full rot. energy` uses the
instantaneous realised geometry, so for vibrating polyatomics it is a useful
diagnostic but not always the exact sampled model energy.

## 5. Vibration

If `vib-mode = rigid`, generation says so directly:

```text
[vibration]
vib-mode              = rigid; harmonic oscillator sampling skipped
```

For a vibrationally active molecule, ICATS prints mode-by-mode harmonic
oscillator information. The key columns are:

- `vstat`: sampled harmonic-oscillator state.
- `Q`, `P`: sampled dimensionless normal coordinate and momentum.
- `QE`, `PE`: coordinate and momentum energy contributions.
- `EE`: total harmonic energy for that mode, in eV.

In analysis, ICATS projects the Cartesian coordinates and velocities back onto
the same normal-mode coordinates. That reconstructed table should agree with
the generated one within the accuracy expected from the model and numerical
projection.

## 6. Intermolecular Analysis

The analysis-side intermolecular block reconstructs the two-body Jacobi
quantities from the Cartesian sample:

```text
[intermolecular]
molecules              = ar_dat x no_dat
angular energy         =    4.79710e-04  eV
radial energy          =    8.02678e-02  eV
total energy           =        0.08075  eV
Jab, vector model      = [      0.00000,      0.00000,      0.00000 ]  au, |Jab| =      0.00000  (QM:   0.00)
L                      = [      0.00000,    115.78237,      0.00000 ]  au, |L| =    115.78237  (QM: 115.28)
P_R                    = [      1.04631,     -0.00000,     13.53450 ]  au, |P_R| =     13.57489  (QM:  13.08)
J = L + Jab            = [      0.00000,    115.78237,      0.00000 ]  au, |J| =    115.78237  (QM: 115.28)
b                      =        4.50000  Ang
phi                    =        1.00000  pi rad
COM 1                  = [      1.93023,      0.00000,     24.96832 ]  Ang
COM 2                  = [     -2.56977,      0.00000,    -33.24117 ]  Ang
```

The sign of a vector component may differ between generation and analysis if
the coordinate convention or equivalent body-frame representation has changed.
For initial-condition sanity checks, compare the magnitudes and the energy
summary first, then inspect signs and frames if the system is anisotropic.

## 7. System Angles

The system-angle block is the two-vector embedding reconstruction:

```text
[system angles]
system Euler           = [    -0.0000,     0.0246,     0.0000 ]  pi rad   # phi, beta, chi
mol 1 BF Euler         = [     0.0000,     0.0000,     0.0000 ]  pi rad   # alpha, beta, gamma
mol 2 BF Euler         = [     0.6413,     0.3425,    -0.5382 ]  pi rad   # alpha, beta, gamma
v1-v2 dihedral         =         0.0000  pi rad
Jacobi R               =       58.38317  Ang
```

Do not confuse `impact-phi` in the input with the reconstructed system Euler
`phi`. The input fixes the lab impact-parameter azimuth. The system Euler
triple is reconstructed later from the full Cartesian geometry.

## 8. Energy Summary

Generation and analysis both print the same table shape:

```text
[energy summary]
component             ammonia_dat/eV     h2o_dat/eV       total/eV
vibrational                 1.0268         0.3820         1.4088
rotational                  0.0082         0.0498         0.0580
velocity                    0.0492         0.0465         0.0957
total                       1.0842         0.4783         1.5625
```

Read the generation table as:

```text
What did ICATS intend to sample?
```

Read the analysis table as:

```text
What do the final Cartesian coordinates and velocities contain when analysed back?
```

For a first check, the generation and analysis totals should be close. For a
strict initial-condition round trip, enable `audit-initial-sample = True`,
which compares the correct vector-model terms directly.

## Worked Example: NH3 + H2O

For a polyatomic pair such as ammonia + water, the `.info` file is most useful
when read section by section rather than top to bottom. A representative
generation energy summary might look like:

```text
[energy summary]
component             ammonia_dat/eV     h2o_dat/eV       total/eV
vibrational                 1.0268         0.3820         1.4088
rotational                  0.0082         0.0498         0.0580
velocity                    0.0492         0.0465         0.0957
total                       1.0842         0.4783         1.5625
```

Start with the `velocity` row. It tells you the incoming translational energy
assigned to the two molecular centres of mass. This should match the beam or
direct-channel setup you intended.

Next read the `rotational` row. For water and ammonia, the sampled rotor energy
comes from the rigid-rotor/vector-model state selected for each molecule. In
the analysis block, compare this to `vector rot. energy`, not blindly to
`full rot. energy`, because the full instantaneous line uses the realised
vibrationally distorted geometry.

Then read the `vibrational` row and the mode tables. If the vibrational
temperature is high, a single sample can carry substantial oscillator energy.
The ensemble should reproduce the intended Husimi/harmonic-oscillator
population, while one sample is just one phase-space point.

Finally check the intermolecular block:

```text
Jab, vector model      = [     -1.43210,      0.28744,      2.91013 ]  au, |Jab| =      3.25660  (QM:   2.81)
L                      = [      0.00000,    114.07766,      0.00000 ]  au, |L| =    114.07766  (QM: 113.58)
J = L + Jab            = [     -1.43210,    114.36510,      2.91013 ]  au, |J| =    114.41125  (QM: 113.91)
b                      =        5.69860  Ang
```

This is the part to inspect for Wang-Landau or angular-momentum-window
questions. In a molecular collision, `J` is not simply `L`; the rotor vector
contribution `Jab = JA + JB` shifts the total angular momentum. Wang-Landau
sampling is useful when you want the accepted ensemble to balance the intended
total-`J` window with a physically broad impact-parameter/orbital-`L`
distribution.

## Audit Block

With `audit-initial-sample = True`, ICATS also prints a compact check:

```text
Audit energy   vib: generation 1.408804 eV, analysis 1.410228 eV, diff 1.424e-03 eV [OK]
Audit energy   rot: generation 0.058004 eV, analysis 0.058009 eV, diff 4.764e-06 eV [OK]
Audit energy   vel: generation 0.095675 eV, analysis 0.095720 eV, diff 4.514e-05 eV [OK]
Audit energy total: generation 1.562482 eV, analysis 1.563956 eV, diff 1.474e-03 eV [OK]
Initial sample audit: OK
```

If this passes, ICATS has generated Cartesian initial conditions consistent
with its own harmonic-oscillator, rigid-rotor, and Jacobi model. It does not
prove that a later dynamics method conserves energy or that the potential
energy surface is appropriate.

## Ensemble Histograms

A single `.info` block checks one sample. For any production calculation,
inspect histograms as well:

```bash
./rd_tutorial_input/histograms/plot_initial.sh
./rd_tutorial_input/histograms/plot_sampled.sh
./rd_tutorial_input/histograms/compare_initial_sampled.py --render-unmatched
```

For Wang-Landau runs, inspect both the Wang-Landau umbrella and the sampled
`L`/`J` distributions. A good-looking single sample does not prove that the
ensemble populated the intended angular-momentum range.
