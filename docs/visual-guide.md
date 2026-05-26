---
layout: default
title: Visual Guide
---

# Visual Guide

These figures are visual anchors for the main ICATS ideas. They are not a
substitute for the equations, but they help connect the output files to the
physical coordinates being sampled.

## Lab And Collision Frames

![Lab-frame collision geometry](assets/figures/labdiag.png)

The lab-frame diagram shows how the two incoming molecular beams define the
experimental geometry. ICATS ultimately converts this into a relative
intermolecular velocity, an initial separation, and a collision angle.

## Intermolecular Angular Momenta

![Intermolecular angular momentum vectors](assets/figures/interdiag.png)

This is the diagram to keep in mind when reading Wang-Landau output. The
orbital angular momentum `L` comes from the incoming collision geometry, while
the molecular rotors contribute `J_A` and `J_B`. ICATS samples and analyses the
vector relation:

$$
\vec{J}=\vec{L}+\vec{J}_A+\vec{J}_B.
$$

## Body-Fixed And Euler Angles

![Body-fixed Euler-angle construction](assets/figures/bffdiag.png)

Molecular rotations and orientations are easiest to understand as rotations
between the space-fixed and body-fixed frames. The same distinction appears in
the logs as space-frame, Eckart-frame, and body-frame angular quantities.

## Vibrational Phase Space

![Husimi vibrational distributions](assets/figures/husimi.png)

ICATS samples harmonic vibrational states and then draws phase-space `Q, P`
values from Husimi-style distributions. The histogram does not mean every
sample has exactly the quantum level energy; it means the phase-space ensemble
represents the chosen oscillator state.

## Rotor State Sampling

![Rotor-state distribution examples](assets/figures/rotato.png)

The rotational state distributions show the two-step character of rotor
sampling: choose a total rotational state and then choose projection or
asymmetric-top vector-model information. The same logic is reconstructed later
in the rotational-analysis block.
