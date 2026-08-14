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
audit-initial-angular-tol = 2.0
audit-initial-vib-tol = 2.0
audit-initial-velocity-tol = 5.0
```

This compares generated sample bookkeeping against the immediate coordinate and
velocity analysis at `t = 0`, before any dynamics are run.

For a guided explanation of the first generation and analysis blocks, see
[Annotated First Output](annotated-output.md).

To run the bundled tutorial audit in the alternative scattering-axis
convention:

```bash
icats.audit-tutorials --nsamp 2 --keep-going --output-frame incoming-k-plus-z
```

This keeps Wang-Landau disabled by default, so it is a quick initial-condition
regression check rather than an umbrella-convergence test.

For a sharper frame-convention smoke test:

```bash
icats.frame-smoke
```

This generates the fixed-plane Ar + NO tutorial once in `output-frame =
internal` and once in `output-frame = incoming-k-plus-z`. It checks that both
initial-sample audits pass and that the plus-Z `.xyz` and `.vel` files are
exactly the `Rx(pi)` transform of the internal-frame files.

The audit compares:

- vibrational energy,
- vector-model rotational energy,
- intermolecular velocity energy,
- total sampled model energy,
- orbital and total angular momenta `L`, `Jab`, and `J`,
- per-molecule vector-model angular momenta,
- impact parameter and relative velocity,
- vibrational normal coordinates where they are defined,
- angular variables using circular angle differences.

Not every printed value is a fair generation-vs-analysis pair. For example,
linear-molecule spin about the molecular axis is gauge-like, atom Euler angles
are arbitrary, and the full instantaneous molecular angular momentum includes
vibrational angular momentum that was not directly sampled as a rigid-rotor
vector. Those quantities remain useful diagnostics, but they are not all strict
pass/fail checks.

Run the audit over every bundled tutorial with:

```bash
icats.audit-tutorials --nsamp 3 --keep-going
```

By default this creates:

```text
smoke_results/initial_audit_<timestamp>/
```

and writes a `summary.tsv` table with the pass/fail status, largest energy
difference, and largest audited state residual for each tutorial. Wang-Landau is
disabled in this sweep by default so the command checks the initial-condition
round trip quickly. Add `--include-wl` only when you deliberately want the
tutorials that request Wang-Landau to build or use their umbrella.

## Analysis Output

The `dynamics*.analinfo` files contain:

- orientation analysis,
- rotational analysis,
- vibrational analysis,
- intermolecular analysis,
- energy summary.

For an initial-condition round trip, compare the sampled rigid-rotor energy
with `rigid rot. energy`. For a propagated snapshot, `Eckart rot. energy`
uses the angular velocity and realised instantaneous geometry.

## Reading Rotational Analysis

A typical polyatomic rotational block contains several related angular momenta.
They are all useful, but they are not interchangeable:

```text
full J, Eckart
rigid J0, Eckart
geometry J(1), Eckart
geometry J(2), Eckart
intrinsic pi, Eckart
closure J, Eckart
rigid rot. energy
Eckart rot. energy
```

`full J, space` is the total instantaneous angular momentum of the
molecule after removing its centre-of-mass motion, expressed in the external
Cartesian frame.

`full J, Eckart` is the same instantaneous angular momentum after the
molecule has been rotated into the Eckart frame. The magnitude should be the
same as the space-frame value, but the components are now molecular-frame
components.

`rigid J0, Eckart` is the reference-geometry contribution
`I0 Omega_Eck`. For generated initial conditions it independently
reconstructs the vector-model angular momentum sampled by ICATS.

`geometry J(1), Eckart` and `geometry J(2), Eckart` are the first- and
second-order changes caused by applying the rotational velocity to a displaced
geometry. `intrinsic pi, Eckart` is

```text
pi = sum_i m_i u_i cross du_i/dt
```

These are distinct physical contributions, not separately sampled rotor
states. `closure J, Eckart` checks that their sum with `rigid J0, Eckart`
recovers `full J, Eckart`.

`rigid rot. energy` uses the reference moments of inertia and `J0`. Use this
line when checking the sampled rigid rotor.

For propagated snapshots, `Eckart rot. energy` is
`0.5 Omega_Eck^T I(u) Omega_Eck`. It excludes the intrinsic angular momentum
`pi`, which must not be treated as rigid rotation.

## Reading Vibrational Analysis

The vibrational block reports each normal mode:

```text
mode  freq  ~vstat  Q  P  QE  PE  EE
```

For propagated snapshots ICATS first evaluates `Omega_Eck` from the angular
decomposition and forms `du_i/dt = v_i^b - Omega_Eck cross x_i`. Residual
linear rotational and translational leakage is projected out before the
normal-mode transformation. `Q` and `P` are the resulting reconstructed
normal-mode coordinate and momentum. `EE` is
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

If no harmonic normal-mode basis is available, the block is shorter. For an
atom it reports:

```text
vibrational space      = none (atom)
```

For a non-atomic molecule without a Hessian it reports an internal residual:

```text
vibrational modes      = unavailable; no Hessian/normal modes
internal residual      = projected outside translation/rotation
residual |dx|          =      0.0000000  Ang
residual |p|           =      0.0000000  au
residual kinetic       =      0.0000000  eV
```

These residuals are useful in rigid or constrained tests. They say whether the
Cartesian sample contains motion outside the reference translation/rotation
space. They are not harmonic vibrational energies, because no Hessian or
frequencies were supplied.

## Reading Intermolecular Analysis

The intermolecular block reports:

```text
angular energy
radial energy
L
P_R
J = L + Jab
Jacobi R
```

These quantities describe the collision geometry and relative motion. They are
especially useful when checking whether `maxl`, `maxb`, and the velocity
settings are physically plausible.

`angular energy` is the centrifugal/orbital part of the two-body Jacobi
motion, computed as `L^2 / (2 mu R^2)`.

`radial energy` is the radial relative-motion part, computed as
`P_R^2 / (2 mu)`.

`L` is the reconstructed orbital angular momentum.

`J = L + Jab` is formed from the orbital angular momentum plus the
molecular vector-model angular momenta:

```text
J = L + J_a + J_b
```

The `Ja/Jb/Jab, full` lines use the full molecular angular momenta from the
Cartesian analysis. The `Ja/Jb/Jab, vector model` lines use the
reference-geometry vector-model molecular angular momenta and are usually the
better comparison
for the sampled total-`J` distribution.

## Reading Energy Summary

During initial-condition generation, `out_full.info` can contain both:

```text
Sample N | generation
Sample N | analysis
```

The `[energy summary]` block under `generation` is the bookkeeping of the model
variables ICATS sampled:
harmonic vibrational energy, vector-model rotor energy, and sampled
intermolecular velocity energy.

The `[energy summary]` block under `analysis` is reconstructed from the
Cartesian coordinates and velocities. For generated samples its rotational
row is the recovered rigid-rotor energy. For propagated snapshots it uses the
Eckart rotational velocity and instantaneous geometry.

The `Velocity` row is the intermolecular centre-of-mass kinetic contribution.
It is not the momentum part of the normal-mode vibrational energy; that is the
`PE` column in the vibrational table.

For the cleanest `t = 0` sampler check, use the initial-sample audit. It also
enforces the angular-decomposition closure and velocity Eckart condition.

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
Jab distribution
intermolecular velocity distribution
Wang-Landau weights, if wang = True
```

Generated tutorial directories include a focused helper for the common system
angular-momentum check:

```bash
./rd_tutorial_input/histograms/plot_orbital_jljab.sh
./rd_tutorial_input/histograms/plot_orbital_correlation.sh
```

This renders the sampled `L`, `J`, and `Jab` histograms without plotting every
molecular internal-coordinate histogram.

If a histogram shows almost all samples at a boundary, the requested range is
probably too narrow or the stored Wang-Landau umbrella is not appropriate for
the calculation.

During long Wang-Landau production sampling, the run logfile also receives live
acceptance diagnostics such as:

```text
Sampling worker 0: accepted=500 trials=8231 trial_acceptance=0.06074 J_cap_rejects=120 WL_rejects=7611 WL_range_rejects=0
```

Large `WL_rejects` indicate a low umbrella acceptance probability, large
`J_cap_rejects` indicate hard total-`J` cutoff pressure, and nonzero
`WL_range_rejects` indicate that generated `J` values are falling outside the
stored WL table.
