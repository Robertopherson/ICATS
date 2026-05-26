---
layout: default
title: Tutorials
---

# Tutorials

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
