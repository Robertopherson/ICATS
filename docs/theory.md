# Initial-Condition Theory

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

## Vibrations

Vibrational coordinates are sampled in the harmonic-oscillator normal-mode
basis. For polyatomic molecules, the vibrational energy includes zero-point
energy, so total vibrational energies can be several eV even at moderate
temperatures when high-frequency modes are present.

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
