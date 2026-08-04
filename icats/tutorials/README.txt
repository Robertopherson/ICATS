ICATS Tutorials

This directory is reserved for named tutorial definitions.

Current tutorial index:
0 -> quickstart
1 -> diatomic_n2_n2_fast
2 -> mixed_h2o_n2
3 -> methane_methane
4 -> single_atom_he_he
5 -> single_atom_diatom_he_n2
6 -> polarized_orientation_he_no
7 -> fixed_plane_atom_diatom_ar_no
8 -> flat_l_atom_diatom_ar_no
9 -> single_atom_diatom_he_n2_wl
10 -> flat_l_atom_diatom_he_n2_wl
11 -> flat_l_diatom_diatom_n2_n2_wl
12 -> wang_landau_nh3_h2o
13 -> paper_nh3_h2o_100k
14 -> npz_output_co2_co2

Generate with:
  icats --tutorial quickstart --setup-only
or
  icats --tutorial 0 --setup-only

Override sample/trajectory counts:
  icats --tutorial quickstart --nsamp 200 --ntraj 10 --setup-only

List options:
  icats --list-tutorials

Prefer tutorial names in saved commands and job scripts because numeric indices
can change when new tutorials are added.
