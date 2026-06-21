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
- `polarized_orientation_he_no`: atom-diatom toy polarization example using a
  user-supplied molecular orientation PDF.
- `fixed_plane_atom_diatom_ar_no`: constrained Ar + NO diagnostic setup. It
  uses `incoming-k = 13.615392`, `fixed-b = 4.5`, `impact-phi = 0.0`,
  `vib-mode = rigid`, `Trot = 0.0`, and `Tvib = 0.0`; it samples only the NO
  orientation and writes combined `out_full.info`, `out_full.xyz`, and
  `out_full.vel` files.
- `flat_l_atom_diatom_ar_no`: Ar + NO diagnostic setup for
  `orbital-sampling = flat-l`. It samples `L` uniformly instead of using the
  geometric `P(L) ~ L` proposal. With rigid NO and `Trot = 0.0`, `Jab = 0` and
  `J = L`, so the L/J histogram check is intentionally easy to interpret.
- `single_atom_diatom_he_n2_wl`: atom-diatom with Wang-Landau weighting.
- `flat_l_atom_diatom_he_n2_wl`: atom-diatom Wang-Landau diagnostic using
  `orbital-sampling = flat-l`. This is the direct companion to
  `single_atom_diatom_he_n2_wl`, but the base orbital proposal is uniform in
  `L` rather than geometric.
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

3. `fixed_plane_atom_diatom_ar_no`
   Shows how to hold selected intermolecular variables fixed: impact parameter
   `b`, lab impact-plane azimuth `impact-phi`, rigid vibration, and zero
   initial diatom angular momentum. This is mainly an initial-condition
   diagnostic tutorial.

4. `polarized_orientation_he_no`
   Shows how to provide a molecule-level orientation PDF. The example keeps NO
   rigid and non-rotating, then biases the NO body axis with a simple tilted
   field model.

5. `flat_l_atom_diatom_ar_no`
   Shows how to sample orbital angular momentum uniformly rather than
   geometrically. Use it to learn the difference between a proposal ensemble and
   a weighted/geometric cross-section ensemble.

6. `quickstart`
   Introduces two polyatomic molecules and the full file pipeline.

7. `methane_methane`
   Shows a heavier symmetric-top case. This is useful for seeing vibrational
   angular momentum in the analysis.

8. `single_atom_diatom_he_n2_wl`
   Adds Wang-Landau weighting in a relatively cheap system.

9. `flat_l_atom_diatom_he_n2_wl`
   Uses the same one-dimensional Wang-Landau correction on total `J`, but with
   uniform `L` proposals. This is a diagnostic for low-`L` coverage and
   reweighting strategy, not a two-dimensional Wang-Landau calculation.

10. `wang_landau_nh3_h2o`
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

Before running `icats.init`, inspect the frame convention in
`tutorial_input.txt`:

```text
output-frame = internal
```

This is the historical ICATS/tutorial convention. If the generated xyz/vel
files must match a scattering or QM input convention where the incoming
relative wave vector is along space-fixed `+Z`, change the line before
generation:

```text
output-frame = incoming-k-plus-z
```

Changing `output-frame` changes Cartesian components and SF/BF Euler angles, so
regenerate the samples after changing it. Magnitudes and energies are unchanged.

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

## Constrained Ar + NO Tutorial

Use this tutorial when you want a small setup that fixes selected variables
instead of sampling the full intermolecular geometry:

```bash
icats --tutorial fixed_plane_atom_diatom_ar_no --setup-only
cd tutorial_fixed_plane_atom_diatom_ar_no
icats.init tutorial_input.txt
```

The generated input uses:

```text
Trot = 0.0
Tvib = 0.0
vib-mode = rigid
incoming-k = 13.615392
fixed-b = 4.5
impact-phi = 0.0
output-frame = internal
printout = 1 1 0 0
```

This freezes the NO bond at the reference geometry, starts the diatom with no
rigid-rotor angular momentum, fixes the impact parameter and lab impact-plane
azimuth, and still samples the NO bond orientation. Inspect `out_full.info`,
`out_full.xyz`, and `out_full.vel` after `icats.init`.

This tutorial is an initial-condition/export example, not a cheap-MINDO
dynamics example. The PySCF MINDO/3 backend used by `run_cheap_dynamics.sh`
does not support Ar for this setup, so use the generated files for inspection
or pass them to an external dynamics/QM code.

## Polarized NO Orientation Tutorial

Use this tutorial to test a user-supplied molecular orientation PDF. The
example is intentionally simple: He is an atom, NO is rigid and non-rotating,
and the main diagnostic is the biased NO orientation.

```bash
icats --tutorial polarized_orientation_he_no --setup-only
cd tutorial_polarized_orientation_he_no
icats.init tutorial_input.txt
python rd_tutorial_input/histograms/plot_polarization_check.py
```

The generated NO molecule file contains:

```text
orientation-mode = pdf orientation_pdfs.py dipole_field_tilted 0.75 0.7853981633974483 0.0
orientation-frame = scattering
orientation-thin = 100
rot-param = euler
```

The referenced function has the form:

```python
def dipole_field_tilted(alpha, beta, gamma, strength, field_theta, field_phi):
    ...
```

It treats the NO body-fixed `z` axis as a toy dipole and the field direction as
already expressed in the ICATS scattering frame. The example field is tilted by
`pi/4`, so the sampled distribution depends on both `alpha` and `beta`. ICATS
multiplies the user PDF by the Euler measure `sin(beta)`, so the user function
should return only the physical angular weight. In particular, the example
function computes the body-axis direction from the supplied Euler angles, takes
its dot product with the field axis, and returns a non-negative value
proportional to `1 + A cos(theta_muE)`.

The plotting helper writes:

```text
rd_tutorial_input/histograms/plots/polarization/polarized_orientation_check.png
```

![Polarized NO tutorial orientation check](assets/figures/polarized-orientation-check.png)

The left panel compares sampled `alpha` with the expected azimuthal trend. The
right panel compares sampled `cos(beta)` with the expected linear trend. With
the default tutorial sample count the curves are not production-smooth, but the
tilt and sign should be visible. Increase `Nsamp` if you want a cleaner
diagnostic plot.

## Flat-L Ar + NO Tutorial

Use this tutorial when you want to oversample the small-impact-parameter region
without running Wang-Landau. The default geometric orbital proposal samples
impact parameter area, so `P(b) db ~ b db` and equivalently `P(L) dL ~ L dL`
at fixed incoming momentum. That is the natural choice for direct geometric
averages, but it can leave relatively few samples at small `b`.

The flat-L tutorial deliberately changes only that proposal:

```bash
icats --tutorial flat_l_atom_diatom_ar_no --setup-only
cd tutorial_flat_l_atom_diatom_ar_no
icats.init tutorial_input.txt
```

The generated input contains:

```text
Trot = 0.0
Tvib = 0.0
vib-mode = rigid
incoming-k = 13.615392
maxb = 8.0
impact-phi = 0.0
output-frame = internal
orbital-sampling = flat-l
wang = False
printout = 0 1 0 0
```

As above, this Ar + NO tutorial is not compatible with the bundled
`run_cheap_dynamics.sh` MINDO helper, because PySCF MINDO/3 does not support Ar
for this setup. It is intended for initial-condition diagnostics and histogram
inspection.

Here `orbital-sampling = flat-l` means ICATS proposes `L` uniformly between
zero and the requested upper limit. It does not impose a flat physical
cross-section distribution. If you use these samples to estimate a geometric
quantity, reweight the samples by the missing area factor, proportional to `L`
or `b` at fixed incoming momentum.

This tutorial is intentionally simple. The NO bond is rigid and `Trot = 0.0`,
so the internal rotor contribution is zero: `Jab = 0`. Therefore the total
angular momentum is just `J = L`. That makes it a clean user check: sampled `L`
and sampled `J` should both look roughly flat.

To make histogram helper scripts at setup time:

```bash
icats --tutorial flat_l_atom_diatom_ar_no --histograms --hist-samples 5000 --setup-only
cd tutorial_flat_l_atom_diatom_ar_no
icats.init tutorial_input.txt
mkdir -p rd_tutorial_input/histograms/plots/sampled
python rd_tutorial_input/histograms/sampled/system/hist_sam_sys_orb_sl.py \
  --no-show --outfile rd_tutorial_input/histograms/plots/sampled/hist_sam_sys_orb_sl
python rd_tutorial_input/histograms/sampled/system/hist_sam_sys_orb_sj.py \
  --no-show --outfile rd_tutorial_input/histograms/plots/sampled/hist_sam_sys_orb_sj
```

Inspect the resulting system `L` and `J` plots under:

```text
rd_tutorial_input/histograms/plots/sampled/
```

The broader `plot_sampled.sh` helper is still available, but it plots every
sampled histogram. For this diagnostic tutorial the two system orbital plots
above are usually the useful check.

For this tutorial the important lesson is the shape of the proposal
distribution. It is not a Wang-Landau run, and it is not a replacement for
checking whether the requested `maxb`, `maxl`, and incoming momentum are
appropriate for the physical calculation.

## Flat-L He + N2 Wang-Landau Diagnostic

Use this tutorial to test the combination of uniform-`L` proposals and the
existing one-dimensional Wang-Landau correction on total angular momentum `J`.
It is intentionally a companion to `single_atom_diatom_he_n2_wl`: the molecule
pair and cheap dynamics settings are the same, but the generated input adds

```text
orbital-sampling = flat-l
wang = True
wlmode = default
wl-target = flat-j
wl-j-range = 60
wl-j-bins = 80
wl-l-cap = 59
wl-flatness = 0.90
wl-nstep = 500
run-mode = fresh
```

This is not a two-dimensional Wang-Landau calculation in `(L,J)`. ICATS still
builds a one-dimensional WL profile for the sampled total `J`. In this tutorial
`wl-target = flat-j`, so the WL acceptance uses `1/Omega(J)` rather than the
usual `(2J+1)/Omega(J)` target. The practical question is whether proposing `L`
uniformly gives better coverage of the low-`L` region while the WL acceptance
keeps the total-`J` ensemble from becoming badly concentrated.

The final sampled `L` and `J` histograms should be broadly even in the useful
range, but they do not have to be mathematically flat bin-by-bin. This is still
a one-dimensional WL correction on `J`, and the accepted ensemble can retain
small structure where `L`, `J`, and the atom-diatom `Jab` distribution overlap
strongly, especially near `J,L = 0`.

The WL controls are deliberately explicit in this tutorial. `wl-j-range = 60`
sets the upper `J` range for the explicit WL density estimate, `wl-j-bins = 80`
sets the resolution within that range, and `wl-l-cap = 59` makes the WL-building
orbital proposal cover the same scale as the production `maxl`. The tutorial
uses the `default` WL profile, with explicit flatness and step settings for a
moderately stricter diagnostic run.

Generate the test directory with histogram helpers:

```bash
icats --tutorial flat_l_atom_diatom_he_n2_wl --histograms --hist-samples 5000 --setup-only
cd tutorial_flat_l_atom_diatom_he_n2_wl
```

Then run the initial-condition generator. This may take noticeably longer than
the non-WL tutorials because it first builds `wang.pkl`:

```bash
icats.init tutorial_input.txt
```

`icats.init` sets writable default cache directories for numba and matplotlib
when the user has not already configured them. The generated tutorial keeps
`workers = 1` for maximum portability. On your laptop or a small workstation,
edit `tutorial_input.txt` to `workers = 4` before running `icats.init` if you
want the WL build to use four parallel workers.

After it finishes, inspect the generated WL file and the log:

```text
rd_tutorial_input/wang.pkl
tutorial_input.txt.logfile
rd_tutorial_input/histograms/wl/
```

The WL helper scripts are written under `rd_tutorial_input/histograms/wl/`.
The exact filenames are generated by the code; list them with:

```bash
find rd_tutorial_input/histograms/wl -maxdepth 1 -type f | sort
```

For the sampled orbital histograms, the most useful quick check is:

```bash
./rd_tutorial_input/histograms/plot_orbital_jljab.sh
```

Inspect:

```text
rd_tutorial_input/histograms/plots/sampled/hist_sam_sys_orb_sl.png
rd_tutorial_input/histograms/plots/sampled/hist_sam_sys_orb_sj.png
rd_tutorial_input/histograms/plots/sampled/hist_sam_sys_orb_sjab.png
```

For a sensible diagnostic run:

- `L` should reflect the flat-L proposal modified by WL acceptance.
- `J` should not show an obvious missing-window or hard cutoff inside the
  requested range.
- `Jab` should remain a molecular-rotor distribution, not something manually
  forced by the WL machinery.
- If `J` only looks acceptable because of rare cancellations between `L` and
  `Jab`, inspect `cos(theta_{L,Jab})` with
  `hist_sam_sys_orb_cosljab_thet.py` before treating the run as useful.

If you want to compare against the default geometric proposal, generate and run
the sibling tutorial:

```bash
icats --tutorial single_atom_diatom_he_n2_wl --histograms --hist-samples 5000 --setup-only
```

Compare the sampled `L`, `J`, and `Jab` plots between the two directories. The
flat-L tutorial is mainly a coverage and weighting diagnostic; it should not be
used as a production scattering ensemble without an explicit weighting strategy.
