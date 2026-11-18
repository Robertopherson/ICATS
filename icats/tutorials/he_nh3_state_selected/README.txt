He + NH3: state-selected initial conditions and scattering opacities
=================================================================

Quick walkthrough (from this directory):
  python run.py --samples 128 --workers 4 --output quick
  python plot.py quick --output quick_plots
128 samples only illustrate the workflow; they are not converged opacities.
Use --prepare-only to generate initial conditions without dynamics.
The driver writes ordinary input.txt, nh3.txt and he.txt in each batch.
You can inspect/edit these and invoke icats.init input.txt yourself in a
new directory. run.py verifies its own inputs and refuses incompatible reuse.

Compare sharp K and width-free interval smoothing, for example:
  python run.py --j 3 --k 2 --width 0 --samples 8192 --workers 4 --output j3_s0
  python run.py --j 3 --k 2 --samples 8192 --workers 4 --output j3_interval
  python plot.py j3_s0 j3_interval --output comparison

The default --polar interval preserves each integer-K population and second
moment using a corrected uniform prior on its unit-width interval. Positive
--width values do not change this law; --width 0 selects sharp K.
Historical Gaussian runs require --polar broadened --width 0.8 explicitly.
Never combine historical Gaussian and interval batches into one ensemble.

Optional larger production calculation:
  python run.py --paper-task 0 --workers 8 --output paper
Run task numbers 0 through 255 (paper.sbatch groups these into 16 Slurm jobs).
Four cases: (j,|K|)=(3,2),(5,4), each sharp and interval-smoothed.
The existing CASES width values are switches (0 sharp, positive interval),
not adjustable widths for the interval method. With --polar broadened the
same tasks reproduce the earlier s=0/0.8 Gaussian comparison instead.
16 one-bohr shells, four batches of 4096 trajectories per shell and case:
262144 trajectories/case, 1048576 in total. Seed base 2026111017;
batch seed = base + state_index*10000 + shell*100 + batch.
Sampling is uniform b within each shell, with identical seeds across methods.
Initial-condition generation is single-worker per batch; changing the number
of dynamics workers does not change the generated ensemble.
  python plot.py paper --output paper_plots
Plots use whatever completed batches are present; inspect their coverage report
before treating a partial run as a paper result. Retain full b=0..16 coverage.

This larger configuration is not the exact trajectory count used in the revised
article. Its interval-smoothed comparisons use 32768 trajectories per incoming
state: 2048 in each of sixteen one-bohr shells. The separate paper-data bundle
is being updated to provide the frozen data and figure commands. The plot.py
here demonstrates ordinary binning with the filters described below, not the
paper's additional energy-transfer cuts and locally renormalized blue curves.

Model and limitations
---------------------
The model.py potential is the analytic He-NH3 MOLSCAT example, Eq. (5.7) of
its manual. B=9.9441, C=6.2294 cm^-1; reduced mass 3.24 amu.
Collision energy 93.8265 cm^-1. The input geometry and explicit masses are
matched to these model parameters, not a new ab initio ammonia structure.
ICATS writes full-precision atomic Cartesian positions/velocities. The driver
recovers COM, rotation, and angular momentum from those exports without
replacing them by a second random orientation or rotational sample.
Only initial/final compact data and audits are retained, not trajectory paths.
Dynamics requires NumPy/SciPy and ICATS's ordinary dependencies; no electronic
structure program, PySCF, or MOLSCAT installation is needed.

Both incoming levels are allowed para-NH3 levels. The input explicitly requests
nh3-para and equal +/-K populations. Broadening preserves the second moment
in the ensemble, not the energy of each individual trajectory. Initial
rotational energies therefore have a spread when smoothing is enabled.
The rigid dynamics has no umbrella inversion or quantum spin symmetry.
The quantum reference averages the two initial inversion partners equally
and sums outgoing K/inversion channels for each j'. See reference.json.

Default plots: nearest-bin j,K assignment, para-allowed final K only, with
assigned rotational energy <= the individual trajectory's initial total
energy. All incoming trajectories remain in the denominator. This rejection
is an approximate analysis prescription, NOT a restoration of quantum dynamics.
Use --no-energy-filter for the symmetry-only comparison, or --no-symmetry-filter
to retain all nearest-bin assignments. The latter still uses an energy filter
unless --no-energy-filter is also supplied. Never silently renormalize.
The substantial j'=6 discrepancy is a limitation, not a hidden missing channel.
Cross sections shown are j-changing: K-only/inversion-only transitions are
not included in that sum. Individual j' panels sum their accepted final K.

For each interval, sigma = mean[2*pi*b*interval_width*outcome] (area units).
Sum over intervals; the opacity is that interval contribution/annular area.
The statistical bars are 1 standard error. Zero observed events do not prove
zero probability. The bundled quantum CSVs are already averaged and are not
rescaled to fit the classical data. Here, 0.5-bohr plotting bins
subdivide each generated one-bohr shell; no data is truncated in integration.

Input options are documented in docs/state_selected_rotors.txt in the source.
The scripts record seeds, inputs, code hashes, geometry checks, and per-trajectory
energy/angular-momentum conservation. Keep the same release and environment
for exact reproduction; cross-platform floating-point integration can differ.
