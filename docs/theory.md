---
layout: default
title: Initial-Condition Theory
---

# Initial-Condition Theory

![Angular momentum decomposition](assets/angular-decomposition.svg)

ICATS samples molecular-scattering initial conditions using separable model
degrees of freedom:

- harmonic oscillator vibrational states,
- rigid-rotor molecular angular momenta,
- molecular orientations,
- intermolecular separation and relative velocity,
- orbital angular momentum or impact parameter.

The generated Cartesian positions and velocities are then analysed back into
the same components. This round trip is useful because it checks whether the
sampled bookkeeping survives conversion to atom-resolved coordinates.

The approximation is the usual entrance-channel one: before collision, the two
molecules are far enough apart that their intramolecular potentials can be
treated independently from the intermolecular interaction. ICATS therefore
builds a sample in internal coordinates first, then embeds it into Cartesian
coordinates for dynamics.

## Vibrations

Vibrational coordinates are sampled in the harmonic-oscillator normal-mode
basis. For polyatomic molecules, the vibrational energy includes zero-point
energy, so total vibrational energies can be several eV even at moderate
temperatures when high-frequency modes are present.

The practical check is direct. At `t = 0`, the analysis should reconstruct the
normal-mode coordinates and momenta and recover the sampled vibrational energy
within the audit tolerance. If that fails, inspect the reference geometry,
normal modes, frequencies, and mass conventions before trusting trajectory
propagation.

## Rotations

The sampled rigid-rotor state defines a vector-model angular momentum and
rotational energy. The analysis also reports a full instantaneous rotational
energy from the realized geometry and velocities.

For vibrating polyatomics, these are not always identical because the realised
geometry can carry vibrational angular momentum. The analysis therefore prints:

```text
Vector Model. Ang.
Vibr. Ang. Mom. (Eckart)
Vector Model Rotational (ev)
Full Rotational Energy (ev)
```

The vector-model rotational energy is the apples-to-apples comparison with the
sampled rigid-rotor energy.

Linear, symmetric, and asymmetric tops differ in how their body-fixed projection
information is represented. For asymmetric tops, ICATS uses the Wang-state
description internally when sampling rotational states. The operational point is
that the sampled rotor state should be checked through the reconstructed
`Vector Model Rotational (ev)` and angular-momentum diagnostics.

## Orientations

The molecular angular momentum vector and the molecular frame are separate
pieces of information. ICATS samples the rotor state and then samples the
orientation of the molecule in the space-fixed frame. Under isotropic
conditions this corresponds to uniform sampling over the relevant Euler-angle
measure. More specialised orientation or polarisation distributions should be
checked with histograms when used.

## Intermolecular Motion

The relative motion of the two molecules is sampled using a separation,
relative velocity, and either an impact parameter or orbital angular momentum.
For many scattering observables, the orbital angular momentum range controls how
large an impact parameter is represented in the ensemble.

The total angular momentum combines:

```text
J = L + J_A + J_B
```

where `L` is the orbital angular momentum and `J_A`, `J_B` are the molecular
rotor angular momenta.

For atom-atom scattering, `J` and `L` are effectively the same angular
momentum. For molecular scattering, total angular momentum is the vector sum of
orbital and molecular rotor contributions. This is why a calculation can have a
sensible impact-parameter distribution but still require extra care to obtain
the desired total-`J` distribution.

The impact parameter is the geometrical variable connected to cross-sectional
area: larger annuli contribute more trajectories through the familiar `b db`
weighting. The angular-momentum view highlights the statistical `2J+1`
degeneracy. ICATS provides tools to sample these consistently enough for
practical quasi-classical trajectory ensembles.

## Why The Round Trip Matters

The code samples model quantities, but dynamics programs usually need atom-wise
Cartesian coordinates and velocities. The round trip is:

```text
model variables -> Cartesian x/v -> reconstructed model variables
```

If this round trip fails, the sampled ensemble is not being represented
faithfully in the files sent to dynamics. The initial-sample audit exists to
catch this before a user spends time on trajectory propagation.

## What Is Approximate

The harmonic oscillator and rigid rotor approximations assume an incoming
molecule near one equilibrium structure. They are usually most reasonable for
small-amplitude vibrations and well separated collision partners. Large
amplitude motion, floppy modes, strong Coriolis coupling, or a poor potential
energy surface can all make the toy model less reliable.

ICATS exposes these approximations in the analysis rather than hiding them. For
example, vibrational angular momentum is reported separately from the
reference-geometry vector-model rotor energy.
