---
layout: default
title: Diagnostics and Audits
---

# Diagnostics and Audits

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

## Histogram Checks

When histogram generation is enabled, plotting helpers are written under:

```text
rd_<run-tag>/histograms/
```

Use these to inspect sampled `J`, `L`, vibrational coordinates, rotor
projections, and Wang-Landau distributions.
