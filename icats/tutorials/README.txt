ICATS Tutorials

This directory is reserved for named tutorial definitions.

Use `icats --list-tutorials` for the current tutorial names and indices.
The revised NH3+H2O example is revision_nh3_h2o. Wang-Landau tutorials and
their precomputed profiles are retired and are not installed.
The separate he_nh3_state_selected tutorial supplies analytic rigid-top
dynamics and saved quantum reference tables.

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
