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
system Euler           = [    -0.0000,     0.0246,     0.0000 ]  pi rad   # SFF -> Jacobi BF; phi, beta, chi
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

If a molecule has no Hessian, ICATS cannot assign internal motion to harmonic
normal modes. The analysis then reports a residual outside the translation and
rotation space:

```text
ar_dat                 = vibrational analysis
vibrational space      = none (atom)
no_dat                 = vibrational analysis
vibrational modes      = unavailable; no Hessian/normal modes
internal residual      = projected outside translation/rotation
residual |dx|          =      0.0000000  Ang
residual |p|           =      0.0000000  au
residual kinetic       =      0.0000000  eV
```

For rigid tests this is a useful check: the residual should be zero or very
small. It is not a harmonic vibrational energy table, because no Hessian or
frequencies were supplied.

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
frame note             = SFF is the lab Cartesian frame; Jacobi BF is the collision frame
system Euler           = [    -0.0000,     0.0246,     0.0000 ]  pi rad   # SFF -> Jacobi BF; phi, beta, chi
mol 1 BF Euler         = [     0.0000,     0.0000,     0.0000 ]  pi rad   # molecule 1 in Jacobi BF; alpha, beta, gamma
mol 2 BF Euler         = [     0.6413,     0.3425,     0.0000 ]  pi rad   # molecule 2 in Jacobi BF; alpha, beta, gamma
v1-v2 dihedral         =         0.0000  pi rad
Jacobi R               =       58.38317  Ang
SFF->Jacobi BF         = [    1.0000,   -0.0000,    0.0000 ]; [    0.0000,    0.9969,    0.0772 ]; [   -0.0000,   -0.0772,    0.9969 ]
Jacobi BF->mol1 BF     = [    1.0000,    0.0000,    0.0000 ]; [    0.0000,    1.0000,    0.0000 ]; [    0.0000,    0.0000,    1.0000 ]
```

Do not confuse `impact-phi` in the input with the reconstructed system Euler
`phi`. The input fixes the lab impact-parameter azimuth. The system Euler
triple is reconstructed later from the full Cartesian geometry.
For a linear molecule, `gamma` in `mol i BF Euler` is an arbitrary spin about
the molecular axis, so ICATS reports it as zero.
The printed matrices are row-vector rotation matrices used by the analysis:
`SFF->Jacobi BF` takes lab Cartesian components into the collision frame, and
`Jacobi BF->mol i BF` takes Jacobi-frame components into the molecule's
Eckart/body-fixed frame.

## 8. Energy Summary

Generation and analysis both print the same table shape:

```text
[energy summary]
component              =   ammonia_dat/eV     h2o_dat/eV       total/eV
vibrational            =           1.0268         0.3820         1.4088
rotational             =           0.0082         0.0498         0.0580
velocity               =           0.0492         0.0465         0.0957
total                  =           1.0842         0.4783         1.5625
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

A small NH3 + H2O run using the quickstart molecule files produced the
following representative `.info` sample:

```text
Nsamp = 12
Tvib = 500.0
Trot = 500.0
Rz = 15
maxl = 182
maxb = 16
wang = False
printout = 0 1 0 0
```

Sample 3 is a useful teaching example because both molecules rotate, both have
normal-mode vibrations, and the analysis reconstructs similar values from the
Cartesian coordinates and velocities.

### Generation: Rotor Draws

```text
Sample 3 | generation
 - Orbital Angular Q.N. L = : 169.52

[orientation]
ammonia_dat            = angular velocity from sampled rotor
quantum J              =           11.0
quantum J_z            =          9.000  au
classical J            = [   -7.0537000,   -1.1159377,    9.0000000 ]  au, |J| =     11.48913  (QM:  11.00)
classical energy       =        0.12929  eV

h2o_dat                = angular velocity from sampled rotor
quantum J              =            2.0
classical J            = [   -2.2262829,    0.0424868,    1.0207152 ]  au, |J| =      2.44949  (QM:   2.00)
classical energy       =        0.02153  eV
```

`L` is the sampled orbital angular-momentum quantum number. The molecular
`quantum J` lines are the rotor states selected for each monomer. The
`classical J` vectors are the vector-model realisation of those quantum labels:
the magnitude is approximately `sqrt(J(J+1))`, which is why `J = 11` appears
as `|J| = 11.48913`.

```text
ammonia_dat            = sampled orientation
alpha,beta,gamma       = [  0.7452936,  0.7338185,  0.7312064 ]  pi rad
h2o_dat                = sampled orientation
alpha,beta,gamma       = [  0.6461051,  0.6799557, -0.2812087 ]  pi rad
total J                =         172.24
```

These `alpha,beta,gamma` triples are molecular Euler angles. They describe
each molecule's body-frame orientation before the Cartesian sample is written.
`total J` is the vector sum of sampled orbital `L` and molecular
`Jab = Ja + Jb`.

### Generation: Vibrations

```text
[vibration]
  :ammonia_dat :
    mode   freq  vstat    Q         P         QE        PE        EE (eV)
    0     652.3  1     1.075534  0.926277  0.046775  0.034693  0.081468
    1    1816.3  0     0.220889  0.101438  0.005494  0.001159  0.006652
...
  :h2o_dat :
    mode   freq  vstat    Q         P         QE        PE        EE (eV)
    0    1756.9  0     0.486248  0.110618  0.025751  0.001333  0.027084
```

`freq` is the harmonic normal-mode frequency. `vstat` is the sampled oscillator
state used by the Husimi-like normal-mode sampler. `Q` and `P` are the sampled
mass/frequency-scaled coordinate and momentum. `QE`, `PE`, and `EE` are the
coordinate, momentum, and total harmonic energy for that mode.

### Generation: Intermolecular Channel

```text
[intermolecular]
relative velocity      =   1435.6131231  m/s
collision energy       =      0.0934994  eV
ammonia_dat kinetic    =      0.0480626  eV
h2o_dat kinetic        =      0.0454368  eV
Ja, vector model       = [      5.42409,      0.00708,    -10.12814 ]  au, |Ja| =     11.48913
Jb, vector model       = [      0.14867,      2.36552,      0.61824 ]  au, |Jb| =      2.44949
L                      = [      0.00000,    170.01775,      0.00000 ]  au, |L| =    170.01775
J = L + Jab            = [      5.57276,    172.39035,     -9.50990 ]  au, |J| =    172.74237
b                      =        8.59126  Ang
impact phi             =        1.00000  pi rad
```

The relative velocity and collision energy are the incoming two-body channel.
`Ja` and `Jb` are molecular angular momenta in the lab frame, `Jab` is their
sum, `L` is orbital angular momentum, and `J = L + Jab` is the total. The
impact parameter follows from `L = mu v_rel b`.

```text
[energy summary]
component              =   ammonia_dat/eV     h2o_dat/eV       total/eV
vibrational            =           0.2478         0.3656         0.6135
rotational             =           0.1293         0.0215         0.1508
velocity               =           0.0481         0.0454         0.0935
total                  =           0.4252         0.4326         0.8578
```

Read this as the sampled model budget: harmonic vibration, vector-model
rotation, and incoming COM motion.

### Analysis: Reconstructing the Cartesian Sample

```text
Sample 3 | analysis

[rotation]
ammonia_dat            = oblate, symmetry constant =   1.00
full J, space          = [      6.78943,     -0.28659,     -9.65061 ]  au, |J| =     11.80309
vector J, space        = [      5.42514,      0.00672,    -10.12783 ]  au, |J| =     11.48935
vibrational J          = [     -0.51418,      1.38125,      0.05324 ]  au, |J| =      1.47481
vector rot. energy     = [      0.06826,      0.00171,      0.05933 ]  eV, |E| =      0.12929
full rot. energy       = [      0.00171,      0.07625,      0.05964 ]  eV, |E| =      0.13759
```

The analysis block is reconstructed from the actual Cartesian coordinates and
velocities. `vector J` is the part corresponding to the intended rigid-rotor
model. `full J` uses the instantaneous molecular geometry and therefore also
sees vibrational angular momentum. Their difference is reported as
`vibrational J`. For sampler checks, compare the generated rotor energy with
`vector rot. energy`; inspect `full rot. energy` to see how the instantaneous
distorted geometry changes the decomposition.

```text
[vibration]
ammonia_dat            = vibrational analysis
 mode    freq     ~vstat       Q         P         QE        PE       EE (eV)
   0    652.3     0.5073  1.075533  0.926191  0.046775  0.034687  0.081461
...
h2o_dat                = vibrational analysis
   0   1756.9    -0.3757  0.486248  0.110606  0.025751  0.001332  0.027083
```

The reconstructed `Q`, `P`, `QE`, and `PE` values should match the generated
ones closely. The `~vstat` value is not resampling a quantum number; it is the
classical reconstructed oscillator energy expressed as `E/omega - 0.5`.

```text
[intermolecular]
angular energy         =    2.30266e-02  eV
radial energy          =    7.02990e-02  eV
total energy           =        0.09333  eV
L                      = [     -0.47993,   -169.76306,     -0.27488 ]  au, |L| =    169.76396
J = L + Jab            = [      5.09247,   -167.39197,     -9.78332 ]  au, |J| =    167.75493
b                      =        8.58642  Ang
```

The analysed collision energy, impact parameter, and angular-momentum
magnitudes are close to the generated values. Signs may differ because the
analysis reports the reconstructed geometry after standard orientation and
frame choices; start by comparing magnitudes and energies.

```text
[energy summary]
component              =   ammonia_dat/eV     h2o_dat/eV       total/eV
vibrational            =           0.2500         0.3643         0.6143
rotational             =           0.1376         0.0205         0.1581
velocity               =           0.0480         0.0454         0.0933
total                  =           0.4356         0.4302         0.8658
```

This is the final round-trip check. The generation total was `0.8578 eV`; the
analysis total is `0.8658 eV`. The small difference comes from reconstructing
the Cartesian geometry into instantaneous rotational and vibrational
components. Large unexpected differences are a sign to inspect the vibrational
table, the vector/full rotational split, and the intermolecular `L`, `Jab`,
and `J` lines.

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
