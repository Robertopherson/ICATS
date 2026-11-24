---
layout: default
title: Historical Sampling Methods
---

# Historical Sampling Methods

## Wang-Landau Reweighting

Earlier versions of ICATS used Wang-Landau (WL) flat-histogram reweighting to
broaden the sampled total-J distribution beyond the geometric orbital proposal.
The idea was to populate a wide range of total-J values approximately uniformly
in log density, allowing cross-section estimates at a range of initial states from
a single ensemble. The WL convergence was controlled by the modification factor
`f`, the flatness criterion, and the number of workers sharing a common umbrella.

This approach was retired for two reasons. First, the target distribution was
identified incorrectly: it was constructed as a histogram proportional to `2J+1`
rather than the correct distribution obtained by coupling the molecular populations
with the accessible orbital and molecular channels. Second, the WL reweighting
inflates or suppresses entire regions of configuration space, which can disguise
sampling errors in the rotor and vibrational distributions. The revised workflow
avoids these issues by using the correct thermal populations and geometric orbital
sampling directly.

Historical WL input keys (`wang`, `wl-target`, `wl-maxf`, `wl-flatfrac`,
`wl-workers`) remain in the code for compatibility with old inputs, but are not
recommended for new calculations. Set `wang = False` (or omit it; False is the
default) and generate a fresh ensemble.

Old Wang-Landau checkpoints, stored umbrellas, or reweighting profiles must not
be reused with the revised sampler. Start in a new directory.

## Old Gaussian Rotor Broadening

Earlier versions used a Gaussian broadening in the projection quantum number
`K_c` (and a heuristic azimuthal distribution) when constructing the vector-model
angular momentum for asymmetric-top molecules. This scheme did not in general
preserve the correct discretised `K` distribution or the second moment
`<K_c^2> = K_0^2`, so the mean rigid-rotor energy was not conserved in the
ensemble.

This was replaced with the interval/von Mises scheme described in the manuscript
(Section III C 1):

- the continuous projection `K_c` is drawn from a unit-width interval centred on
  the selected integer `K_0`;
- a relative-entropy correction enforces `<K_c^2> = K_0^2`;
- the azimuthal angle follows a von Mises distribution that similarly preserves
  the transverse second moments.

The old Gaussian scheme is not in the current release. If old inputs set
`rotor-broadening = gaussian`, that key will be ignored or raise an error; use
the defaults (`rotor-polar = interval`, `rotor-azimuth = vonmises`) instead.

## Current Approach

The current recommended workflow for a thermal ensemble is:

```text
wang = False
rotor-polar = interval
rotor-azimuth = vonmises
```

These are the defaults. For the corrected canonical rotor populations and
geometric orbital sampling see the `revision_nh3_h2o` tutorial and
[Theory](theory.md). For cross-section estimation see [Diagnostics](diagnostics.md).
