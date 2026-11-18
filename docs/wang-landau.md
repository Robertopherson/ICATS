---
layout: default
title: Retired Wang-Landau Workflow
---

# Retired Wang-Landau Workflow

The revised tutorials do not use Wang-Landau reweighting. Set `wang = False`
and generate a fresh ensemble rather than reusing old checkpoints or profiles.

Geometric impact-parameter sampling supplies the incident area measure.
The total-J reference also includes the allowed molecular and orbital channels;
it is not generally proportional to 2J+1 alone. For a flat orbital proposal,
apply the corresponding area weights when estimating cross sections.

See [Tutorials](tutorials.md), [Theory](theory.md) and
[Diagnostics](diagnostics.md) for the revised workflow. Historical WL options
remain in the code for compatibility, but are not recommended for reproducing
the revised manuscript.
