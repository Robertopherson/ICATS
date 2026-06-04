---
layout: default
title: Initial-Condition Theory
---

# Initial-Condition Theory

![ICATS analysis map](assets/angular-decomposition.svg)

For the main coordinate and angular-momentum diagrams, see the
[Visual Guide](visual-guide.md).

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

## Why These Sampling Tools

ICATS uses three specialised ideas because the three parts of the entrance
channel have different natural variables.

For vibrations, the natural model is a harmonic oscillator in normal-mode
phase space. ICATS therefore samples a vibrational quantum number and then
samples dimensionless `Q, P` coordinates from a Husimi-style phase-space
distribution. This gives Cartesian displacements and velocities that retain the
harmonic oscillator energy and uncertainty-like spread expected from the chosen
state.

$$
E_{\mathrm{HO},k}
  = \frac{1}{2}\omega_k\left(Q_k^2 + P_k^2\right),
\qquad
E_{v,k}
  = \omega_k\left(v_k+\frac{1}{2}\right).
$$

The first expression is the classical phase-point energy reconstructed from the
sampled `Q, P` values. The second is the harmonic-oscillator level energy used
when assigning Boltzmann populations. A single Husimi phase-space draw does not
have to make these two numbers identical mode by mode.

For molecular rotations, the natural variables are angular momenta rather than
Cartesian velocities. ICATS therefore uses a quasi-classical vector model:
sample a rotational state, construct a classical angular momentum vector with
length close to `sqrt(J(J+1))`, then convert that vector into an angular
velocity field on the reference geometry.

$$
|\vec{J}| = \sqrt{J(J+1)}, \qquad
\vec{\omega}_{\mathrm{ref}} = I_{\mathrm{ref}}^{-1}\vec{J},
\qquad
\vec{v}_{\mathrm{rot},i} = \vec{\omega}_{\mathrm{ref}}\times \vec{x}_{\mathrm{ref},i}.
$$

For intermolecular scattering, the natural classical variable is the impact
parameter, while the partial-wave view is organised by total angular momentum.
For atoms these almost coincide because `J = L`. For molecules, however,
`J = L + Jab`, where `Jab` comes from the molecular rotors. Wang-Landau sampling
is used to balance these two views without manually forcing the angle between
`L` and `Jab`.

$$
d\sigma \propto b db,
\qquad
P(L) dL \propto (2L+1) dL,
\qquad
\vec{J} = \vec{L}+\vec{J}_{AB}.
$$

$$
\vec{J}_{AB}=\vec{J}_A+\vec{J}_B,
\qquad
P_{\mathrm{acc}}(J)\propto \frac{2J+1}{\Omega_{\mathrm{WL}}(J)}.
$$

Here $\Omega_{\mathrm{WL}}(J)$ denotes the Wang-Landau estimate of the density
of trial samples produced at total angular momentum `J`.

## Distributions Used In ICATS

Most ICATS input options eventually define one of the following distributions.
The names in the middle column are the code-level ideas a user will see in
logs, histogram file names, or input options.

| Type | Sampled quantity | Distribution or weight |
| --- | --- | --- |
| Vibrational state | normal-mode quantum numbers `vstat` | Boltzmann harmonic-oscillator populations, $P(v)\propto \exp[-E_v/(k_B T_{\mathrm{vib}})]$ |
| Vibrational phase space | normal-mode `Q, P` | Husimi-style harmonic-oscillator phase-space density for the sampled `vstat` |
| Rigid-rotor total angular momentum | molecular `J` | Boltzmann rotor populations, with $(2J+1)$ degeneracy for isotropic tops and state-specific weights for anisotropic/asymmetric tops |
| Symmetric-top projection | body-fixed `K` or projection-like coordinate | Boltzmann projection distribution at fixed `J` |
| Asymmetric-top state | Wang-basis eigenstate at fixed `J` | Boltzmann distribution over asymmetric-rotor eigenenergies and symmetry labels |
| Asymmetric vector model | projection spread and unresolved azimuth | rejection-sampled Gaussian-sine, azimuthal, or Bingham-like auxiliary distributions chosen from the Wang-state expectation values |
| Molecular orientation | Euler angles | isotropic orientations by default, or user-supplied/fixed/read orientation distributions when requested |
| Impact parameter | `b` | geometric incoming-flux measure, $P(b) db \propto b db$ on the requested interval |
| Orbital angular momentum | `L` | semiclassical partial-wave measure, $P(L) dL \propto (2L+1) dL$ |
| Total angular momentum | `J` | target partial-wave measure, $P(J) dJ \propto (2J+1) dJ$ |
| Orbital azimuth | `phi` | uniform angle when `phisample = True` |
| Relative speed | intermolecular velocity | Crossed-beam molecular speed distributions, direct relative-speed/channel input, or Maxwell-Boltzmann relative-speed sampling, depending on input mode |
| Wang-Landau correction | accepted total `J` | rejection weight proportional to $(2J+1)/\Omega_{\mathrm{WL}}(J)$ |

The important practical point is that these distributions are not all sampled
at the same level. Vibrational and rotor state distributions define internal
molecular states. Orientational and phase-space distributions turn those states
into Cartesian coordinates and velocities. Intermolecular distributions define
the incoming collision geometry. Wang-Landau then acts on the combined trial
sample, after `L`, `J_A`, and `J_B` have already produced a candidate total
angular momentum.

## Intermolecular Velocity Ensembles

The incoming translational motion can be interpreted in two different ways.

In a crossed-beam interpretation, each molecule has its own beam speed
distribution. ICATS samples the two molecular centre-of-mass speeds separately
from the molecule-file `vel` entries, then combines them using the requested
`beam-angle`. This produces the Jacobi relative speed and therefore the
collision energy:

$$
v_{\mathrm{rel}}^2
= v_A^2 + v_B^2 - 2 v_A v_B \cos\varphi ,
\qquad
E_{\mathrm{coll}} = \frac{1}{2}\mu v_{\mathrm{rel}}^2 .
$$

This is the natural model for a crossed-beam experiment where the beam
velocities and crossing angle are known.

In a direct-channel interpretation, the user specifies the relative incoming
channel itself. The equivalent quantities are:

$$
p_0 = \mu v_{\mathrm{rel}},
\qquad
E_{\mathrm{coll}} = \frac{p_0^2}{2\mu},
\qquad
k = \frac{p_0}{\hbar}.
$$

Since ICATS uses atomic units internally, \(\hbar=1\), so the numerical values
of `incoming-p0` and `incoming-k` are the same when both are expressed in
atomic units. This direct-channel setup is closest to the usual partial-wave
language, where an incoming channel is labelled by \(k\) and semiclassically
\(L \simeq k b\).

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

## Euler Angles And Output Names

ICATS uses the same rotation convention as the manuscript, but older internal
keys used a less careful naming scheme. The current logs and manual use these
names:

| Quantity | Manuscript-style output | Meaning |
| --- | --- | --- |
| molecular orientation | `alpha, beta, gamma` | Euler rotation from the space-fixed frame into the molecule/Eckart frame. |
| system Jacobi embedding | `phi, beta, chi` | Euler rotation from the space-fixed frame into the two-vector body-fixed frame of the collision complex. |
| impact-parameter azimuth | `impact-phi` | Input/sample coordinate fixing the lab collision plane; not the full system Euler `phi`. |

The molecular Euler triple follows the usual `Z-Y-Z`/line-of-nodes
construction:

$$
R(\alpha,\beta,\gamma)
  = R_Z(\alpha)R_Y(\beta)R_Z(\gamma).
$$

In the molecular output, `alpha` is the first azimuthal rotation about the
space-fixed axis, `beta` is the polar line-of-nodes rotation, and `gamma` is the
final body-fixed spin about the molecule's own axis. For atoms there is no
molecular orientation space. For linear molecules, the final spin about the
bond is arbitrary; ICATS reports `gamma = 0` for that gauge coordinate.

For the two-body collision analysis, the first embedding vector is the Jacobi
COM-to-COM vector:

$$
\vec{R}_{AB} = \vec{R}_A-\vec{R}_B.
$$

ICATS computes the system angles from the generated Cartesian coordinates as

$$
\phi = \mathrm{atan2}(R_y,R_x),
\qquad
\beta = \cos^{-1}\left(\frac{R_z}{|\vec{R}_{AB}|}\right),
$$

then uses a second molecular embedding vector to define the final azimuth
`chi` about the system body-fixed `z` axis. The sampled data still includes a
legacy key `theta`; it is an alias for this same system polar angle `beta`.

The practical output mapping is:

| Output block/key | Read as |
| --- | --- |
| `Orientations (alpha,beta,gamma)` | generated molecular orientation. |
| `Orientation, Euler Angles (alpha,beta,gamma)` | reconstructed molecular/Eckart orientation from Cartesian coordinates. |
| `System Euler (phi,beta,chi)` | reconstructed two-vector Jacobi/system embedding. |
| `Molecule (1) BF Euler`, `Molecule (2) BF Euler` | molecule frame angles after rotating into the system body-fixed frame. |
| `2bJac.theta` | legacy alias for `2bJac.beta`. |

This distinction matters in constrained runs. For example, `impact-phi = 0`
fixes the incoming impact-parameter plane, but it does not force a diatomic
molecule to be coplanar. The diatom orientation is still sampled through its
molecular `alpha,beta` angles, while its molecular `gamma` spin is arbitrary
and therefore reported as zero.

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

or, in the same notation as the manuscript,

$$
E_{\mathrm{mode}}
  = \frac{1}{2}\omega\left(Q^2+P^2\right)
  = Q_E + P_E.
$$

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

equivalently,

$$
\tilde{v} = \frac{E_{\mathrm{mode}}}{\omega}-\frac{1}{2}.
$$

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

or

$$
|\vec{J}| = \sqrt{J(J+1)}.
$$

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

or

$$
\vec{\omega}_{\mathrm{ref}}
  = I_{\mathrm{ref}}^{-1}\vec{J}_{\mathrm{vector}},
\qquad
\vec{v}_{\mathrm{rot},i}
  = \vec{\omega}_{\mathrm{ref}}\times\vec{x}_{\mathrm{ref},i}.
$$

where `I_ref` and `x_ref` belong to the reference/equilibrium Eckart geometry.
This is deliberate. The vector model describes the persistent rigid-rotor
motion of the reference structure. It should not be reconstructed from the
instantaneous vibrationally distorted geometry.

The corresponding reference-geometry rotor energy is

```text
E_vector = 0.5 * sum_i J_vector_i^2 / I_ref_i
```

or

$$
E_{\mathrm{vector}}
  = \frac{1}{2}\sum_i
    \frac{J_{\mathrm{vector},i}^2}{I_{\mathrm{ref},i}}.
$$

This is the energy that should be compared to the sampled rigid-rotor energy
when validating the initial condition.

## Rotor Types

ICATS classifies each molecule from the principal moments of inertia of the
reference geometry. The code records a rotor type such as `atom`, `linear`,
`prolate`, `oblate`, `spherical`, `asym-prolate`, `asym-oblate`, or
`asym-spherical`.

Atoms have no molecular rotational or vibrational space in this model, so their
molecular `Trot` and `Tvib` are forced to zero.

Linear rotors have one vanishing principal moment. ICATS samples the total
rotational angular momentum and treats the body-fixed projection consistently
with the linear-rotor limit.

Near-symmetric tops, labelled `prolate` or `oblate`, are handled with the
standard vector-model picture: sample total `J`, sample the body-fixed
projection, choose the unresolved azimuthal angle, and build the classical
angular momentum vector on the corresponding cone.

Spherical tops are close to isotropic in their rotational constants. Their
rotational energy is essentially controlled by `J(J+1)`, while the direction of
the angular momentum is sampled isotropically.

Asymmetric tops require the most machinery. ICATS diagonalises the asymmetric
rotor Hamiltonian in a Wang-state basis for each `J`. It stores the eigenstate
energies, symmetry labels, and expectation values such as `<J_x^2>`, `<J_y^2>`,
and `<J_z^2>`. The extended vector model then samples a classical vector whose
projection spread is consistent with those eigenstate expectation values. For
near-prolate or near-oblate asymmetric tops the sampling remains close to a
single projection axis. For more spheroidal asymmetric tops, the code uses the
Wang-state symmetry information to choose the relevant principal axis and then
rotates the sampled vector back into the molecular frame.

The point of this machinery is practical: the user can give a general
polyatomic geometry and Hessian, and ICATS will choose the closest rotor model
instead of requiring the user to hand-code a separate sampler for each top.

## Rotational Analysis

`molecules.CalcRotEner` starts from the Cartesian positions and velocities of
one molecule. It removes the molecular COM and calculates several related but
not identical angular momenta.

### `full J, space`

This is the total instantaneous molecular angular momentum in the space-fixed
frame:

```text
J_full_space = sum_i m_i r_i x v_i
```

or

$$
\vec{J}_{\mathrm{full}}^{\mathrm{space}}
  = \sum_i m_i \vec{r}_i\times\vec{v}_i.
$$

where `r_i` and `v_i` are measured after removing the molecular COM.

### `full J, Eckart`

The same instantaneous molecule is rotated into the Eckart frame and the
angular momentum is recomputed:

```text
J_full_eckart = sum_i m_i r_i^E x v_i^E
```

or

$$
\vec{J}_{\mathrm{full}}^{E}
  = \sum_i m_i \vec{r}_i^{E}\times\vec{v}_i^{E}.
$$

The magnitude should match the space-frame value apart from numerical noise,
but the components are expressed in the Eckart molecular frame.

### `vector J, Eckart`

This is the reference-geometry rotor component extracted from the same Eckart
velocities:

```text
J_vector_eckart = sum_i m_i x_ref_i x v_i^E
```

or

$$
\vec{J}_{\mathrm{vector}}^{E}
  = \sum_i m_i \vec{x}_{\mathrm{ref},i}\times\vec{v}_i^{E}.
$$

The key difference is that the cross product uses the reference geometry
`x_ref`, not the realised distorted geometry `r_i^E`. This is the analysis-side
version of the rotor angular momentum assumed by the sampler.

### `vector J, space`

This is the same vector-model angular momentum rotated back to the space-fixed
frame. It is useful when comparing molecular rotor vectors to the
intermolecular `L` and total `J` vectors.

### `vibrational J`

The code defines this residual as:

```text
J_vib_eckart = J_full_eckart - J_vector_eckart
```

or

$$
\vec{J}_{\mathrm{vib}}^{E}
  = \vec{J}_{\mathrm{full}}^{E}
  - \vec{J}_{\mathrm{vector}}^{E}.
$$

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

### `vector rot. energy`

This is computed from `J_vector_eckart` and the reference principal moments of
inertia:

```text
E_vector_i = 0.5 * J_vector_i^2 / I_ref_i
```

or

$$
E_{\mathrm{vector},i}
  = \frac{1}{2}\frac{J_{\mathrm{vector},i}^2}{I_{\mathrm{ref},i}}.
$$

The printed line gives the three principal-axis contributions and their sum.
This is the correct comparison to the sampled rotor energy at `t = 0`.

### `full rot. energy`

For this quantity, the code diagonalises the inertia tensor of the realised
instantaneous geometry. It then computes a rigid-body rotational energy from
the full instantaneous angular momentum in that instantaneous principal-axis
frame:

```text
E_full_i = 0.5 * J_full_i^2 / I_inst_i
```

or

$$
E_{\mathrm{full},i}
  = \frac{1}{2}\frac{J_{\mathrm{full},i}^2}{I_{\mathrm{inst},i}}.
$$

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

or

$$
|\vec{L}| = \sqrt{L(L+1)}.
$$

The two molecular vector-model angular momenta are rotated into the space frame
and added:

```text
Jab = Ja + Jb
J_total = L + Jab
```

or

$$
\vec{J}_{AB} = \vec{J}_A+\vec{J}_B,
\qquad
\vec{J} = \vec{L}+\vec{J}_{AB}.
$$

For atom-atom scattering, `J` and `L` are essentially the same. For molecular
scattering, `Jab` can be comparable to `L`, especially at low collision angular
momentum. This is why a good `L` or impact-parameter distribution does not
automatically guarantee a good total-`J` distribution.

This is also why Wang-Landau sampling appears in the code. The geometric
incoming-flux measure naturally weights impact parameter, and therefore `L`.
The partial-wave view naturally weights total angular momentum, and therefore
`J`. For molecules these two one-dimensional distributions are coupled through
`J = L + Jab`. ICATS uses Wang-Landau to estimate the `J` density of states
generated by independent `L` and `Jab` sampling, then reweights accepted samples
toward the desired `2J + 1` distribution without manually forcing the angle
between `L` and `Jab`.

More practically, the larger-`L` part of a crossed-beam ensemble is often
dominated by orbital angular momentum and already resembles the geometric
impact-parameter measure. Wang-Landau is most useful in the intermediate region
where `L` and `Jab` mix appreciably. The intended outcome is a balanced accepted
ensemble: the `L` and `J` histograms both remain close to their intended
degeneracy-weighted forms, without creating a strong artificial correlation
between `L` and `Jab`.

The semiclassical map used throughout the log is:

```text
|A|^2 = q(q+1)
q = 0.5 * (-1 + sqrt(1 + 4 |A|^2))
```

or

$$
|\vec{A}|^2 = q(q+1),
\qquad
q = \frac{-1+\sqrt{1+4|\vec{A}|^2}}{2}.
$$

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

or

$$
E_L = \frac{|\vec{L}|^2}{2\mu |\vec{R}|^2}.
$$

where `mu` is the reduced mass. The printed `Init. Radial Energy` is

```text
E_R = |P_R|^2 / (2 mu)
```

or

$$
E_R = \frac{|\vec{P}_R|^2}{2\mu}.
$$

The lines `b` and `phi` in the `[intermolecular]` block report the reconstructed
impact parameter and azimuthal angle:

```text
b = |L| / (mu |v_rel|)
phi = atan2(L_x, -L_y)
```

or

$$
b = \frac{|\vec{L}|}{\mu |\vec{v}_{\mathrm{rel}}|},
\qquad
\phi = \mathrm{atan2}(L_x,-L_y).
$$

The lines labelled `Tot. Mol Ja/Jb/Jab` use the full molecular angular momenta
from the Cartesian analysis. The lines labelled `Vec Model Ja/Jb/Jab` use the
reference-geometry vector-model molecular angular momenta. For checking the
sampled total angular momentum distribution, the vector-model `Jab` is usually
the more relevant quantity, because it corresponds to how the sample was
generated.

## Energy Summary Blocks

The analysis file contains matched `[energy summary]` blocks during generation
and analysis. Their position in the file tells you whether the values came from
the sampled model variables or from the reconstructed Cartesian sample.

### `Sample N | generation`

This block is the bookkeeping from the sampled model variables. It uses:

- harmonic oscillator energies sampled for each molecule,
- rotor energies assigned by the vector-model/reference-geometry construction,
- intermolecular COM kinetic energies from the sampled relative velocity.

This is the intended model energy before a dynamics code has acted on the
sample. It does not include a later trajectory-code potential energy.

### `Sample N | analysis`

This block is reconstructed from the Cartesian coordinates and velocities. It
uses:

- vibrational energies from `molecules.CalcInterEner`,
- rotational energies from `molecules.CalcRotEner`,
- intermolecular kinetic energies from `iscattering.CalcInterMolMomentum`.

At present, the printed rotational total in this summary follows the full
instantaneous rotational-energy analysis. For vibrating polyatomics, this means
it may not be identical to the sampled vector-model rotor energy. To check the
sampler itself, compare the sampled rotor energy to the explicit
`vector rot. energy` line, not to the full instantaneous rotational energy.

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

Use `vector rot. energy` when asking:

```text
Did ICATS recover the rigid-rotor energy it sampled?
```

Use `full rot. energy` when asking:

```text
What is the instantaneous rigid-body rotational energy of this Cartesian
snapshot if I use the realised distorted geometry?
```

Use `vibrational J` when asking:

```text
How much molecular angular momentum is not described by the reference-geometry
vector-model component in this snapshot?
```

Use `angular energy`, `radial energy`, `L`, `P_R`, and `J = L + Jab` in the
`[intermolecular]` block when asking:

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
