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

A typical polyatomic rotational block contains:

```text
Full Ang. Mom. (Eckart)
Vector Model. Ang.  (Eckart)
Vibr. Ang. Mom. (Eckart)
Vector Model Rotational (ev)
Full Rotational Energy (ev)
```

For checking the sampled rigid rotor, use `Vector Model Rotational (ev)`.
`Vibr. Ang. Mom. (Eckart)` tells you how much instantaneous angular momentum is
carried by the vibrational motion in that snapshot.

## Reading Vibrational Analysis

The vibrational block reports each normal mode:

```text
mode  freq  ~vstat  Q  P  QE  PE  EE
```

`Q` and `P` are the reconstructed normal-mode coordinate and momentum. `EE` is
the corresponding harmonic oscillator energy. In a good initial-condition round
trip, these values should agree with the generated sample within the audit
tolerance.

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
