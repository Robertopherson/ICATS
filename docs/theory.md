---
layout: default
title: Initial-Condition Theory
---

# Initial-Condition Theory

![ICATS analysis map](assets/angular-decomposition.svg)

This page explains the model ICATS uses and, just as importantly, how the code
turns that model into the quantities printed in the analysis files. The aim is
not to reproduce the full manuscript. The aim is to make the program output
readable.

The central approximation is the entrance-channel separation used in the
manuscript. Before collision, the two molecules are assumed to be far enough
apart that the intramolecular potentials are independent of the intermolecular
interaction. ICATS therefore builds each initial condition from four pieces:

- molecular vibrations, treated as harmonic oscillator normal modes,
- molecular rotations, treated with a rigid-rotor/vector-model description,
- molecular orientation in the space-fixed frame,
- intermolecular Jacobi motion, described by relative velocity, separation,
  impact parameter, orbital angular momentum, and total angular momentum.

The code first samples these model variables, then converts them into Cartesian
atomic positions and velocities. The analysis code then performs the reverse
projection:

```text
model variables -> Cartesian x/v -> reconstructed model variables
```

This round trip is why the analysis output is useful even before dynamics. It
checks whether the sampled internal state survived the conversion to atom-wise
Cartesian coordinates and velocities.

## Code-Level Workflow

For each sample, `iscattering.GenerateSample` follows this order:

1. `InitializeSample`
   Creates a fresh sample object containing two molecule sample objects.

2. `SampleOrbitalL`
   Samples the orbital angular momentum quantum number `L` and a classical
   vector whose length is `sqrt(L(L+1))`.

3. `SampleRigidRotorState0`
   Samples the molecular rotor states for molecule 0 and molecule 1.

4. `SampleOrientat0`
   Samples molecular orientations and sets molecular angular velocities.

5. `SampleJ`
   Forms the total angular momentum vector from the orbital and molecular
   rotor vectors.

6. Optional Wang-Landau rejection
   If `wang = True`, the trial is accepted or rejected using the stored or newly
   calculated `wang.pkl` umbrella.

7. `SampleHOVibrState`
   Samples harmonic oscillator vibrational states and applies the normal-mode
   displacements and velocities.

8. `SetInterZDist`, `SampleInterMolZVeloc`, and `SetImpactParam`
   Place the molecules at the requested initial separation, assign their
   relative velocity, and shift them by the sampled impact parameter.

9. `SummarizeLogEnergy(..., False)`
   Prints the generation-side energy bookkeeping.

10. `AnalyseSample`
    Reconstructs orientation, rotational, vibrational, intermolecular, and
    energy-decomposition quantities from the Cartesian sample.

The generation and analysis blocks are intentionally both printed. When they
disagree, the disagreement tells you where the model-to-Cartesian conversion or
the analysis assumptions need attention.

## Frames Used By The Analysis

Several frames appear in the output:

- `space`: the external Cartesian frame used for the trajectory files.
- molecular COM frame: each molecule after removing its centre-of-mass position
  and velocity.
- Eckart frame: the instantaneous molecule rotated onto the reference geometry
  by the mass-weighted best-fit/Eckart alignment.
- instantaneous principal-axis frame: the frame that diagonalises the inertia
  tensor of the realised, possibly vibrationally distorted geometry.
- Jacobi/intermolecular frame: the two molecular centres of mass treated as a
  two-body problem.

The distinction between the Eckart frame and the instantaneous principal-axis
frame matters. The vector-model rotor is defined with the reference geometry.
The full instantaneous rotational analysis is defined from the actual distorted
geometry at that snapshot.

## Vibrational Sampling

The molecule setup constructs normal modes from the reference geometry and
Hessian. In code, the normal-mode transformation matrices are stored on the
molecule object and used by `molecules.SetHOVibrState` and
`molecules.CalcInterEner`.

For each normal mode, ICATS samples:

- a vibrational quantum number `vstat` from a Boltzmann distribution,
- dimensionless normal coordinate `Q`,
- dimensionless normal momentum `P`.

The phase-space `Q, P` values are drawn from the Husimi-style distributions used
by the sampler. The harmonic oscillator energy for one mode is then

```text
E_mode = 0.5 * omega * (Q^2 + P^2)
```

where `omega` is the mode frequency in atomic units. This is the quantity
printed as `EE (eV)` in the vibrational table. The coordinate and momentum
parts are printed separately as:

```text
QE = 0.5 * omega * Q^2
PE = 0.5 * omega * P^2
EE = QE + PE
```

The approximate state label `~vstat` printed by the analysis is

```text
~vstat = E_mode / omega - 0.5
```

This is not a new quantum measurement. It is a convenient way to express the
reconstructed classical oscillator energy in units of one vibrational quantum.

## Vibrational Analysis

`molecules.CalcInterEner` reconstructs the vibrational table from Cartesian
coordinates and velocities:

1. remove the molecular COM position and velocity,
2. rotate the instantaneous molecule into the Eckart frame,
3. subtract the reference geometry,
4. project out translation and rotation using the reference
   translation/rotation vectors,
5. transform the remaining displacement and momentum into normal-mode
   coordinates,
6. compute `Q`, `P`, `QE`, `PE`, and `EE`.

This is why the initial-sample audit is meaningful. The generated `Q, P`
variables are converted to Cartesian coordinates and velocities, then the
analysis tries to recover the same `Q, P` from the final Cartesian sample.

If the molecule is an atom, or if no vibrational normal-mode space exists, the
analysis prints `No vibrational space`.

## Rotational Sampling

Molecular rotations are sampled with the quasi-classical vector model. For a
sampled rotor state, the code builds a classical angular momentum vector with
length approximately

```text
|J| = sqrt(J(J+1))
```

and with a body-fixed projection determined by the sampled rotor state. For
symmetric-top-like cases this resembles the usual `J, K` vector-model picture.
For asymmetric tops, ICATS uses the asymmetric-top/Wang-state machinery in the
molecule setup to decide how the projection information is represented.

After the angular momentum vector is selected, `molecules.SetAngularVeloc`
constructs the rotational velocity field as

```text
omega = I_ref^-1 J_vector
v_rot = omega x x_ref
```

where `I_ref` and `x_ref` belong to the reference/equilibrium Eckart geometry.
This is deliberate. The vector model describes the persistent rigid-rotor
motion of the reference structure. It should not be reconstructed from the
instantaneous vibrationally distorted geometry.

The corresponding reference-geometry rotor energy is

```text
E_vector = 0.5 * sum_i J_vector_i^2 / I_ref_i
```

This is the energy that should be compared to the sampled rigid-rotor energy
when validating the initial condition.

## Rotational Analysis

`molecules.CalcRotEner` starts from the Cartesian positions and velocities of
one molecule. It removes the molecular COM and calculates several related but
not identical angular momenta.

### `Full Ang. Mom. (space)`

This is the total instantaneous molecular angular momentum in the space-fixed
frame:

```text
J_full_space = sum_i m_i r_i x v_i
```

where `r_i` and `v_i` are measured after removing the molecular COM.

### `Full Ang. Mom. (Eckart)`

The same instantaneous molecule is rotated into the Eckart frame and the
angular momentum is recomputed:

```text
J_full_eckart = sum_i m_i r_i^E x v_i^E
```

The magnitude should match the space-frame value apart from numerical noise,
but the components are expressed in the Eckart molecular frame.

### `Vector Model. Ang. (Eckart)`

This is the reference-geometry rotor component extracted from the same Eckart
velocities:

```text
J_vector_eckart = sum_i m_i x_ref_i x v_i^E
```

The key difference is that the cross product uses the reference geometry
`x_ref`, not the realised distorted geometry `r_i^E`. This is the analysis-side
version of the rotor angular momentum assumed by the sampler.

### `Vector Model. Ang. (space)`

This is the same vector-model angular momentum rotated back to the space-fixed
frame. It is useful when comparing molecular rotor vectors to the
intermolecular `L` and total `J` vectors.

### `Vibr. Ang. Mom. (Eckart)`

The code defines this residual as:

```text
J_vib_eckart = J_full_eckart - J_vector_eckart
```

It measures how much instantaneous molecular angular momentum is carried by
the vibrational distortion and vibrational velocity in the realised snapshot.
In the small-amplitude Eckart approximation, this term is expected to be small
or to average away over vibrational phase, but it does not have to be exactly
zero in a single Cartesian sample.

This label should not be read as a separately sampled conserved quantum number.
It is a diagnostic of the difference between the full instantaneous angular
momentum and the reference-geometry vector-model angular momentum.

Do not read this residual as implying a simple scalar energy identity such as
`E_full = E_vector + E_vib`. The angular momenta are useful diagnostic vectors,
but the energy expressions use different inertia tensors and can contain
cross/higher-order geometry effects.

### `Vector Model Rotational (ev)`

This is computed from `J_vector_eckart` and the reference principal moments of
inertia:

```text
E_vector_i = 0.5 * J_vector_i^2 / I_ref_i
```

The printed line gives the three principal-axis contributions and their sum.
This is the correct comparison to the sampled rotor energy at `t = 0`.

### `Full Rotational Energy (ev)`

For this quantity, the code diagonalises the inertia tensor of the realised
instantaneous geometry. It then computes a rigid-body rotational energy from
the full instantaneous angular momentum in that instantaneous principal-axis
frame:

```text
E_full_i = 0.5 * J_full_i^2 / I_inst_i
```

This is a real diagnostic of the Cartesian snapshot, but it is not the same
quantity as the sampled vector-model rotor energy for a vibrating molecule. It
can include instantaneous vibrational angular momentum and small higher-order
geometry effects.

## Intermolecular Sampling

The intermolecular variables describe the motion of the two molecular centres
of mass. ICATS samples an orbital angular momentum `L`, a relative speed, an
initial separation `Rz`, and an impact parameter `b`.

The sampled orbital vector has length

```text
|L| = sqrt(L(L+1))
```

The two molecular vector-model angular momenta are rotated into the space frame
and added:

```text
Jab = Ja + Jb
J_total = L + Jab
```

For atom-atom scattering, `J` and `L` are essentially the same. For molecular
scattering, `Jab` can be comparable to `L`, especially at low collision angular
momentum. This is why a good `L` or impact-parameter distribution does not
automatically guarantee a good total-`J` distribution.

The semiclassical map used throughout the log is:

```text
|A|^2 = q(q+1)
q = 0.5 * (-1 + sqrt(1 + 4 |A|^2))
```

The printed `(QM: ...)` values are this inverse map from a classical vector
length back to an effective quantum number.

## Intermolecular Analysis

`iscattering.CalcInterMolMomentum` reconstructs the two-body Jacobi quantities
from the Cartesian sample:

1. compute the molecular COM positions and velocities,
2. remove the total system COM velocity,
3. form the relative vector `R = R_0 - R_1`,
4. form the relative velocity,
5. project the velocity into radial and angular components,
6. compute orbital angular momentum, radial momentum, impact parameter, and
   total angular momentum.

The printed `Init. Angular Energy` is

```text
E_L = |L|^2 / (2 mu |R|^2)
```

where `mu` is the reduced mass. The printed `Init. Radial Energy` is

```text
E_R = |P_R|^2 / (2 mu)
```

The line `Cylindrical Coords` reports the reconstructed impact parameter and
azimuthal angle:

```text
b = |L| / (mu |v_rel|)
phi = atan2(L_x, -L_y)
```

The lines labelled `Tot. Mol Ja/Jb/Jab` use the full molecular angular momenta
from the Cartesian analysis. The lines labelled `Vec Model Ja/Jb/Jab` use the
reference-geometry vector-model molecular angular momenta. For checking the
sampled total angular momentum distribution, the vector-model `Jab` is usually
the more relevant quantity, because it corresponds to how the sample was
generated.

## Energy Decomposition Blocks

The analysis file contains two energy decomposition blocks during generation
and one block when analysing existing trajectory files. Their names are easy to
misread.

### `Energy Decomposition (From Generation)`

This block is the bookkeeping from the sampled model variables. It uses:

- harmonic oscillator energies sampled for each molecule,
- rotor energies assigned by the vector-model/reference-geometry construction,
- intermolecular COM kinetic energies from the sampled relative velocity.

This is the intended model energy before a dynamics code has acted on the
sample. It does not include a later trajectory-code potential energy.

### `Energy Decomposition (From Sample)`

This block is reconstructed from the Cartesian coordinates and velocities. It
uses:

- vibrational energies from `molecules.CalcInterEner`,
- rotational energies from `molecules.CalcRotEner`,
- intermolecular kinetic energies from `iscattering.CalcInterMolMomentum`.

At present, the printed rotational total in this summary follows the full
instantaneous rotational-energy analysis. For vibrating polyatomics, this means
it may not be identical to the sampled vector-model rotor energy. To check the
sampler itself, compare the sampled rotor energy to the explicit `Vector Model
Rotational (ev)` line, not to the full instantaneous rotational energy.

The `Velocity` row in these summaries is the translational/intermolecular COM
kinetic contribution split between the two molecular partners. It is not a
separate intramolecular vibrational velocity term.

### Initial-Sample Audit

The optional audit exists because of this distinction. When
`audit-initial-sample = True`, ICATS compares:

- generated vibrational energy against reconstructed vibrational energy,
- generated rotor energy against reconstructed vector-model rotational energy,
- generated intermolecular velocity energy against reconstructed velocity
  energy,
- the corresponding model total.

This is a test of the initial-condition round trip, not a test of a later
dynamics method. If the audit passes and a trajectory later drifts in energy,
the likely source is the dynamics/integrator/potential-energy-surface side, not
the initial-condition conversion itself.

## What To Use For Which Question

Use `Vector Model Rotational (ev)` when asking:

```text
Did ICATS recover the rigid-rotor energy it sampled?
```

Use `Full Rotational Energy (ev)` when asking:

```text
What is the instantaneous rigid-body rotational energy of this Cartesian
snapshot if I use the realised distorted geometry?
```

Use `Vibr. Ang. Mom. (Eckart)` when asking:

```text
How much molecular angular momentum is not described by the reference-geometry
vector-model component in this snapshot?
```

Use `Init. Angular Energy`, `Init. Radial Energy`, `Init. Ang. Momentum`, and
`Tot Ang. Momentum` when asking:

```text
Is the collision geometry and total angular momentum range sensible?
```

Use the histogram scripts when asking:

```text
Does the ensemble, not just one sample, populate the intended J, L, b,
vibrational, rotational, and orientation ranges?
```

## Limits Of The Model

The harmonic oscillator and rigid rotor approximations assume an incoming
molecule near one reference structure. They are most appropriate for small
amplitude vibrations and well separated collision partners. Large-amplitude
motion, floppy modes, strong Coriolis coupling, centrifugal distortion, or a
poor potential-energy surface can all make the decomposition less clean.

ICATS tries to expose those limitations rather than hiding them. The separate
vector-model, full-rotational, vibrational-angular-momentum, and
intermolecular-analysis lines are there so that a user can see when the simple
model remains self-consistent and when the Cartesian sample contains extra
structure that the idealised model does not fully absorb.
