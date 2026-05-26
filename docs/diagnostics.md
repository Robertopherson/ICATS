---
layout: default
title: Diagnostics and Audits
---

# Diagnostics and Audits

Diagnostics are not just for debugging code. They are how a user decides
whether an initial-condition ensemble is physically and numerically sensible
before running many trajectories.

## Initial-Sample Audit

Enable the generation-time audit with:

```text
audit-initial-sample = True
audit-initial-energy-tol = 0.02
audit-initial-angular-tol = 0.0
```

This compares generated sample bookkeeping against the immediate coordinate and
velocity analysis at `t = 0`, before any dynamics are run.

The audit compares:

- vibrational energy,
- vector-model rotational energy,
- intermolecular velocity energy,
- total sampled model energy.

Run the audit over every bundled tutorial with:

```bash
icats.audit-tutorials
```

By default this creates:

```text
smoke_results/initial_audit_<timestamp>/
```

and writes a `summary.tsv` table with the pass/fail status and largest
generation-vs-analysis energy difference for each tutorial.

Angular quantities are currently captured but not enforced by default. The
analysis already reports vibrational angular momentum, but a pass/fail angular
criterion should be defined carefully because frame conventions matter for atoms
and linear/symmetric molecules.

## Analysis Output

The `dynamics*.analinfo` files contain:

- orientation analysis,
- rotational analysis,
- vibrational analysis,
- intermolecular analysis,
- energy decomposition.

For vibrating polyatomics, compare sampled rigid-rotor energy with the analysis
`Vector Model Rotational (ev)` line. The `Full Rotational Energy (ev)` line is a
different decomposition that includes the realised instantaneous geometry.

## Reading Rotational Analysis

A typical polyatomic rotational block contains several related angular momenta.
They are all useful, but they are not interchangeable:

```text
Full Ang. Mom. (Eckart)
Vector Model. Ang.  (Eckart)
Vibr. Ang. Mom. (Eckart)
Vector Model Rotational (ev)
Full Rotational Energy (ev)
```

`Full Ang. Mom. (space)` is the total instantaneous angular momentum of the
molecule after removing its centre-of-mass motion, expressed in the external
Cartesian frame.

`Full Ang. Mom. (Eckart)` is the same instantaneous angular momentum after the
molecule has been rotated into the Eckart frame. The magnitude should be the
same as the space-frame value, but the components are now molecular-frame
components.

`Vector Model. Ang. (Eckart)` is reconstructed with the reference geometry
rather than the distorted instantaneous geometry. This is the component that
corresponds most directly to the rigid-rotor vector sampled by ICATS.

`Vibr. Ang. Mom. (Eckart)` is the residual:

```text
Full Ang. Mom. (Eckart) - Vector Model. Ang. (Eckart)
```

It is an instantaneous diagnostic of vibrational angular momentum and
higher-order geometry effects. It is not a separately sampled rotor state.

`Vector Model Rotational (ev)` uses the reference moments of inertia and the
vector-model angular momentum. Use this line when checking the sampled rigid
rotor.

`Full Rotational Energy (ev)` uses the realised instantaneous geometry and its
instantaneous inertia tensor. This can differ from the vector-model energy in
vibrating polyatomics.

## Reading Vibrational Analysis

The vibrational block reports each normal mode:

```text
mode  freq  ~vstat  Q  P  QE  PE  EE
```

`Q` and `P` are the reconstructed normal-mode coordinate and momentum. `EE` is
the corresponding harmonic oscillator energy. In a good initial-condition round
trip, these values should agree with the generated sample within the audit
tolerance.

The columns mean:

- `freq`: harmonic normal-mode frequency.
- `~vstat`: reconstructed classical oscillator energy expressed as
  `E / omega - 0.5`.
- `Q`: mass/frequency-scaled normal coordinate.
- `P`: mass/frequency-scaled normal momentum.
- `QE`: coordinate part of the harmonic energy.
- `PE`: momentum part of the harmonic energy.
- `EE`: `QE + PE`.

## Reading Intermolecular Analysis

The intermolecular block reports:

```text
Init. Angular Energy
Init. Radial Energy
Init. Ang. Momentum
Init. Rad. Momentum
Tot Ang. Momentum
Jacobi R
```

These quantities describe the collision geometry and relative motion. They are
especially useful when checking whether `maxl`, `maxb`, and the velocity
settings are physically plausible.

`Init. Angular Energy` is the centrifugal/orbital part of the two-body Jacobi
motion, computed as `L^2 / (2 mu R^2)`.

`Init. Radial Energy` is the radial relative-motion part, computed as
`P_R^2 / (2 mu)`.

`Init. Ang. Momentum` is the reconstructed orbital angular momentum `L`.

`Tot Ang. Momentum` is formed from the orbital angular momentum plus the
molecular vector-model angular momenta:

```text
J = L + J_a + J_b
```

The `Tot. Mol Ja/Jb/Jab` lines use the full molecular angular momenta from the
Cartesian analysis. The `Vec Model Ja/Jb/Jab` lines use the reference-geometry
vector-model molecular angular momenta and are usually the better comparison
for the sampled total-`J` distribution.

## Reading Energy Decomposition

During initial-condition generation, `out_full.info` can contain both:

```text
Energy Decomposition (From Generation)
Energy Decomposition (From Sample)
```

`From Generation` is the bookkeeping of the model variables ICATS sampled:
harmonic vibrational energy, vector-model rotor energy, and sampled
intermolecular velocity energy.

`From Sample` is reconstructed from the Cartesian coordinates and velocities.
Its rotational row follows the full instantaneous rotational-energy analysis.
For vibrating polyatomics this row can differ from the sampled rotor energy
because the realised geometry can carry vibrational angular momentum.

The `Velocity` row is the intermolecular centre-of-mass kinetic contribution.
It is not the momentum part of the normal-mode vibrational energy; that is the
`PE` column in the vibrational table.

For the cleanest t=0 sampler check, use the initial-sample audit. It compares
generated rotor energy to `Vector Model Rotational (ev)`, not to the full
instantaneous rotational-energy row.

## Histogram Checks

When histogram generation is enabled, plotting helpers are written under:

```text
rd_<run-tag>/histograms/
```

Use these to inspect sampled `J`, `L`, vibrational coordinates, rotor
projections, and Wang-Landau distributions.

Useful first plots are:

```text
J distribution
L distribution
intermolecular velocity distribution
Wang-Landau weights, if wang = True
```

If a histogram shows almost all samples at a boundary, the requested range is
probably too narrow or the stored Wang-Landau umbrella is not appropriate for
the calculation.
