---
layout: default
title: Tutorials
---

# Tutorials

Tutorials are small, reproducible examples that exercise different parts of the
program. They are not meant to be polished production scattering calculations.
Their purpose is to teach the workflow and catch common setup mistakes.

List available tutorials:

```bash
icats --list-tutorials
```

Generate a named tutorial:

```bash
icats --tutorial single_atom_diatom_he_n2 --setup-only
```

Current tutorials cover:

- `quickstart`: NH3 + H2O baseline workflow.
- `diatomic_n2_n2_fast`: fast N2 + N2 environment sanity check.
- `mixed_h2o_n2`: mixed rotor types.
- `methane_methane`: heavier symmetric-top example.
- `single_atom_he_he`: atom-atom edge case.
- `single_atom_diatom_he_n2`: atom-diatom template.
- `single_atom_diatom_he_n2_wl`: atom-diatom with Wang-Landau weighting.
- `wang_landau_nh3_h2o`: NH3 + H2O with Wang-Landau weighting.
- `npz_output_co2_co2`: dual xyz/vel and NPZ output.

Tutorials are feature demonstrations. Some use inexpensive mock or
semiempirical settings so that users can learn the workflow before replacing the
toy dynamics with more serious calculations.

## Recommended Learning Order

1. `single_atom_he_he`
   This is the simplest possible case. There are no molecular rotations or
   vibrations, so the output focuses on intermolecular motion.

2. `single_atom_diatom_he_n2`
   Adds one linear rotor. This is the easiest place to learn molecular angular
   momentum output.

3. `quickstart`
   Introduces two polyatomic molecules and the full file pipeline.

4. `methane_methane`
   Shows a heavier symmetric-top case. This is useful for seeing vibrational
   angular momentum in the analysis.

5. `single_atom_diatom_he_n2_wl`
   Adds Wang-Landau weighting in a relatively cheap system.

6. `wang_landau_nh3_h2o`
   Demonstrates Wang-Landau sampling in a more demanding polyatomic system.

## Running A Tutorial

The standard pattern is:

```bash
icats --tutorial quickstart --setup-only
cd tutorial_quickstart
icats.init tutorial_input.txt
./run_cheap_dynamics.sh
./run_analysis.sh
```

To generate fewer samples for a fast test:

```bash
icats --tutorial quickstart --nsamp 2 --ntraj 1 --setup-only
```

## What To Inspect

After `icats.init`, inspect:

```text
tutorial_input.txt.logfile
rd_tutorial_input/
rd_tutorial_input/outputs/
```

After analysis, inspect:

```text
rd_tutorial_input/outputs/dynamics0.analinfo
```

Search for:

```text
[energy summary]
[rotation]
[vibration]
[intermolecular]
```

These blocks are the fastest way to understand what ICATS sampled and what the
analysis reconstructed.
