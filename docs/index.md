---
layout: default
title: ICATS Manual
---

# ICATS Manual

![ICATS logo](assets/icats-logo.svg)

![ICATS workflow](assets/pipeline.svg)

ICATS generates molecular-scattering initial conditions and analysis diagnostics
for quasi-classical trajectory workflows. This manual is a working guide to the
whole pipeline: install the code, generate initial conditions, run the provided
tutorials, inspect Wang-Landau sampling, run cheap MINDO dynamics, and interpret
the analysis output.

The intended reader is a new research student who wants to run the examples,
understand what files are being produced, and learn enough of the theory to make
sensible checks. ICATS does not replace judgement about the potential energy
surface or dynamics method. Its job is to generate and analyse the initial
conditions cleanly.

## Contents

- [Installation](installation.md)
- [Quick Start](quick-start.md)
- [Tutorials](tutorials.md)
- [Input Files](input-files.md)
- [Wang-Landau Sampling](wang-landau.md)
- [Initial-Condition Theory](theory.md)
- [Dynamics and Analysis Pipeline](pipeline.md)
- [Diagnostics and Audits](diagnostics.md)
- [Troubleshooting](troubleshooting.md)
- [School and Outreach Videos](videos.md)

## What ICATS Does

ICATS starts from two molecule definitions and a plain-text scattering input
file. It samples vibrational, rotational, orientational, and intermolecular
degrees of freedom, converts those samples into Cartesian coordinates and
velocities, and writes files that can be passed to a trajectory code.

The same code can analyse Cartesian trajectory output and reconstruct the
internal components: vibrational mode energies, molecular angular momenta,
intermolecular angular momentum, and energy decomposition.

## What To Read First

1. Install ICATS using [Installation](installation.md).
2. Run the [Quick Start](quick-start.md).
3. Use [Diagnostics and Audits](diagnostics.md) to check that the generated
   initial conditions make sense.
4. Read [Initial-Condition Theory](theory.md) when the output terms need
   interpretation.
