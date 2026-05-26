---
layout: default
title: Annotated First Output
---

# Annotated First Output

This page shows how to read the first useful ICATS output. It is meant to be
read with a small tutorial run open in another terminal.

For the most explicit first run, ask ICATS to write the generation log and the
initial-sample audit:

```text
printout = 0 1 0 0
audit-initial-sample = True
audit-initial-energy-tol = 0.02
```

Then run:

```bash
icats.init tutorial_input.txt
```

The generation log is written as:

```text
out_full.info
```

If the tutorial suppresses `out_full.info`, either enable the second `printout`
flag as above or read the later trajectory-analysis files:

```text
rd_<run-tag>/outputs/dynamics*.analinfo
```

The numbers below are representative excerpts. The exact values change with
the molecule, random seed, temperature, and angular-momentum range.

## 1. Sample Header

Each generated initial condition starts with a sample number:

```text
###################################################
############## Sample Number 0
###################################################
 - Orbital Angular Q.N. L = : 113.58
```

This is the model-side draw before ICATS converts the sampled variables into
Cartesian atom positions and velocities. `L` is the orbital angular-momentum
quantum number associated with the incoming collision geometry.

If Wang-Landau is enabled, the program may try many trial angular momenta before
one sample is accepted. That rejection step happens before the final Cartesian
sample is written.

## 2. Molecular Rotor State

A molecular rotational block looks like:

```text
 :h2o_dat :
    Quantum    J   :       4, J_z :  3.847, Equilib. Energy : 0.02879 eV
    Classical |J|  :    4.47, J_z :  1.259, Energy : 0.04977 eV
    Cartesian Vec. : -2.6489804  3.3759998  1.2591777
```

The first line names the sampled rotor state. The second line is the
quasi-classical vector-model angular momentum used to build the rotational
velocity field. The Cartesian vector is the same angular momentum expressed in
the molecular frame used at that point in the setup.

For checking the initial-condition sampler, the important later comparison is
not the full instantaneous rotational energy. It is:

```text
Vector Model Rotational (ev)
```

That line is reconstructed from the final Cartesian coordinates and velocities
using the same reference-geometry vector model.

## 3. Vibrational State

A vibrational generation block looks like:

```text
mode   freq  vstat    Q         P         QE        PE        EE (eV)
0     652.3  1     0.308314 -1.653056  0.003844  0.110494  0.114338
```

The columns mean:

- `vstat`: harmonic-oscillator state chosen from the Boltzmann population.
- `Q`, `P`: sampled normal-mode phase-space coordinate and momentum.
- `QE`, `PE`: coordinate and momentum contributions to the harmonic energy.
- `EE`: `QE + PE`, in eV.

The sampled Husimi phase point does not have to have exactly the energy
`omega(v + 1/2)` in every single draw. The ensemble, not one row, is what should
represent the chosen vibrational temperature.

## 4. Intermolecular Geometry

The intermolecular setup contains the collision separation, velocity, angular
momenta, and impact parameter:

```text
-Z coordinate distance (Ang)      :   15.00000
-Z Inter-mol Velocity Sample (m/s):    746.5005542   -705.7184817 = 1452.2190359
-Orbital Angular mometnum         :    0.00000 114.07766 0.00000, |L| = 114.07766
-Total   Angular mometnum         :   -2.40033 111.44411 3.24126, |J| = 111.51707
-InterM. Cylindrical Coor         :    5.69860 Ang, 1.00000 pi rad
```

These are the quantities to inspect when asking whether `Rz`, `maxb`, `maxl`,
`maxj`, and the velocity settings are physically sensible.

In molecular scattering, `J` is not just `L`. ICATS forms:

$$
\vec{J}=\vec{L}+\vec{J}_A+\vec{J}_B.
$$

This is why a clean impact-parameter distribution does not automatically give a
clean total-`J` distribution for polyatomic systems.

## 5. Generation Energy Summary

The generation-side energy block is bookkeeping from the sampled model
variables:

```text
########## Energy Decomposition (From Generation) ##############
                       ammonia_dat               h2o_dat               Total (eV)
Vibrational  :              1.0268               0.3820               1.4088
Rotational   :              0.0082               0.0498               0.0580
Velocity     :              0.0492               0.0465               0.0957
Total Energy :              1.0842               0.4783               1.5625
```

This block answers:

```text
What did ICATS intend to sample?
```

It is the cleanest statement of the model energy before any dynamics code or
potential-energy surface has acted on the sample.

## 6. Reconstructed Rotational Analysis

After ICATS has built Cartesian coordinates and velocities, it analyses those
same coordinates:

```text
Full Ang. Mom. (space)      :    1.79275 -0.43980  1.61681, |J| = 2.45387
Full Ang. Mom. (Eckart)     :    2.12232 -1.18324 -0.34230, |J| = 2.45387
Vector Model. Ang. (Eckart) :    1.89920 -1.54788  0.00001, |J| = 2.45008
Vibr. Ang. Mom. (Eckart)    :    0.22312  0.36464 -0.34231, |J| = 0.54765
Vector Model Rotational (ev):    0.00495  0.00329  0.00000, Sum = 0.00823 eV
Full Rotational Energy (ev) :    0.00187  0.00629  0.00008, Sum = 0.00823 eV
```

Read these lines as different projections of the same Cartesian snapshot:

- `Full Ang. Mom.` uses the realised instantaneous geometry.
- `Vector Model. Ang.` uses the reference geometry used by the sampler.
- `Vibr. Ang. Mom.` is the residual difference between those two angular
  momenta in the Eckart frame.
- `Vector Model Rotational (ev)` is the right comparison to the sampled
  rigid-rotor energy.
- `Full Rotational Energy (ev)` is useful, but it is not the quantity sampled
  by the rigid-rotor vector model for a vibrating molecule.

For atoms, there is no molecular rotational space. For linear molecules, some
principal-axis components are zero or nearly zero by construction.

## 7. Reconstructed Vibrational Analysis

The reconstructed vibrational table should resemble the generated `Q`, `P`,
and energy values:

```text
mode    freq     ~vstat       Q         P         QE        PE       EE (eV)
0      652.3     0.9138  0.308314 -1.653030  0.003844  0.110491  0.114334
```

This block answers:

```text
Can the Cartesian coordinates and velocities be projected back onto the
normal-mode variables that ICATS intended to sample?
```

Small differences are normal because the code is doing a finite-precision
coordinate transformation, Eckart alignment, and projection. Large differences
mean the molecule setup, normal modes, frame reconstruction, or initial
condition may need attention.

## 8. Reconstructed Intermolecular Analysis

The intermolecular analysis reconstructs the Jacobi two-body quantities:

```text
Init. Angular Energy   : 1.20012e-02
Init. Radial  Energy   : 8.37188e-02
Tot Ang.   Momentum au :   -2.15919 -116.39845 3.33320, |J| = 116.46618
Cylindrical Coords     :    5.68170 Ang    0.00067 pi rad
Jacobi R (Ang)         :   16.04600
```

This block is the place to check:

- whether the impact parameter is in the intended range,
- whether the initial separation is still large enough,
- whether `L` and `J` are in sensible windows,
- whether the relative velocity is plausible for the intended collision.

## 9. Sample Energy Summary

The sample-side energy block is reconstructed from the Cartesian sample:

```text
########## Energy Decomposition (From Sample) ################
                       ammonia_dat               h2o_dat               Total (eV)
Vibrational  :              1.0280               0.3822               1.4102
Rotational   :              0.0082               0.0496               0.0578
Velocity     :              0.0492               0.0465               0.0957
Total Energy :              1.0854               0.4784               1.5638
```

This block answers:

```text
What does the final Cartesian sample contain when ICATS analyses it back?
```

For a first check, compare this block to `From Generation`. For the strict
sampler check, rely on the audit block because it compares rotational energy to
the vector-model reconstruction rather than to the full instantaneous rotor
diagnostic.

## 10. Initial-Sample Audit

The audit is the fastest way to decide whether the initial-condition round trip
is behaving:

```text
########## Initial Sample Audit ###############################
Audit energy   vib: generation 1.408804 eV, analysis 1.410228 eV, diff 1.424e-03 eV [OK]
Audit energy   rot: generation 0.058004 eV, analysis 0.058009 eV, diff 4.764e-06 eV [OK]
Audit energy   vel: generation 0.095675 eV, analysis 0.095720 eV, diff 4.514e-05 eV [OK]
Audit energy total: generation 1.562482 eV, analysis 1.563956 eV, diff 1.474e-03 eV [OK]
Initial sample audit: OK
```

If this passes, ICATS has generated Cartesian initial conditions that are
consistent with its own harmonic-oscillator, rigid-rotor, and Jacobi model. It
does not prove that a later dynamics method conserves energy or that the
potential-energy surface is appropriate.

If this fails, reduce the problem:

1. set `Nsamp = 1`,
2. set `workers = 1`,
3. keep `audit-initial-sample = True`,
4. inspect the first failing component: `vib`, `rot`, `vel`, or `total`.

## 11. Ensemble Histograms

A single sample can look reasonable while the ensemble is wrong. For any real
calculation, inspect histograms:

```bash
./rd_tutorial_input/histograms/plot_initial.sh
./rd_tutorial_input/histograms/plot_sampled.sh
```

For Wang-Landau runs also inspect:

```bash
cd rd_tutorial_input/histograms/wl
python wl_td_plot.py
python wl_wl_plot.py
```

The most important first plots are:

- sampled `J`,
- sampled `L`,
- impact parameter `b`,
- intermolecular velocity,
- Wang-Landau weights, if `wang = True`.

## Wang-Landau Runtime

Wang-Landau is often the expensive part of an ICATS setup. When `wang = True`
and `rd_<run-tag>/wang.pkl` is absent, ICATS must first estimate the
total-`J` density of states before it can generate the final accepted samples.

For small atom-diatom examples this can be quick. For polyatomic systems, large
`maxj`, strict `wlmode` settings, or many WL bins, the first umbrella build can
take many minutes and can plausibly take an hour or more. Using a few workers,
for example:

```text
workers = 4
```

can help, but the scaling is not perfect. Use one or two workers while
debugging input files, then use four or more only once the setup is known to be
correct.

Once a compatible `wang.pkl` exists, later runs reuse it and skip the expensive
umbrella build. If you need to rebuild it, move the old file aside rather than
overwriting it:

```bash
mv rd_tutorial_input/wang.pkl rd_tutorial_input/wang_previous.pkl
icats.init tutorial_input.txt
```
