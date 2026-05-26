# Quick Start

Generate the default tutorial:

```bash
icats --tutorial quickstart --setup-only
cd tutorial_quickstart
```

Generate initial conditions:

```bash
icats.init tutorial_input.txt
```

Run the cheap dynamics demonstration:

```bash
./run_cheap_dynamics.sh
```

Analyse generated trajectories:

```bash
./run_analysis.sh
```

The main generated directory is:

```text
rd_tutorial_input/
```

Important outputs include:

```text
rd_tutorial_input/outputs/out_*.md.xyz
rd_tutorial_input/outputs/out_*.md.vel
rd_tutorial_input/outputs/dynamics*.analinfo
```

The `dynamics*.analinfo` files contain orientation, rotational, vibrational,
intermolecular, and energy-decomposition diagnostics for each saved frame.
