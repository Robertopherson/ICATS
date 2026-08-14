import json
import pickle
import subprocess
import sys

import numpy as np

from icats.cli import _generate_tutorial


def test_paper_tutorial_generation_and_export(tmp_path):
    tutorial = tmp_path / "paper"
    _generate_tutorial("paper_nh3_h2o_100k", tutorial)

    input_text = (tutorial / "tutorial_input.txt").read_text()
    assert "Nsamp = 100000" in input_text
    assert "maxj = 70" in input_text
    assert "beam-angle = 90.0" in input_text
    assert "workers = 4" in input_text
    assert "wl-nstep = 4000" in input_text
    assert "wl-flatness = 0.95" in input_text
    assert "wl-tol = 1.0000001" in input_text
    readme_text = (tutorial / "tutorial_README.txt").read_text()
    assert "incoming relative momentum points along space-fixed -Z" in readme_text
    assert "incoming relative momentum points along space-fixed +Z" in readme_text
    assert "historical ICATS/tutorial convention" not in readme_text
    assert "vel  = 600 100 3" in (tutorial / "paper_ammonia_dat.txt").read_text()
    assert "vel  = 800 100 3" in (tutorial / "paper_h2o_dat.txt").read_text()
    assert (tutorial / "export_paper_histogram_data.py").exists()
    assert (tutorial / "plot_paper_histograms.py").exists()
    assert (tutorial / "run_paper_tutorial.sh").exists()
    assert (tutorial / "tutorial_environment.sh").exists()
    assert (tutorial / "find_python.sh").exists()
    runner_text = (tutorial / "run_paper_tutorial.sh").read_text()
    assert "SBATCH" not in runner_text
    assert "SLURM" not in runner_text
    assert "icats.init tutorial_input.txt" in runner_text
    assert "source ./tutorial_environment.sh" in runner_text
    assert "source ./find_python.sh" in runner_text
    assert "NUMBA_CACHE_DIR" not in runner_text
    assert '"$PYTHON" export_paper_histogram_data.py' in runner_text
    assert '"$PYTHON" plot_paper_histograms.py' in runner_text
    assert not (tutorial / "run_paper_slurm.sh").exists()
    assert not (tutorial / "run_cheap_dynamics.sh").exists()
    assert "NUMBA_CACHE_DIR" in (tutorial / "tutorial_environment.sh").read_text()
    assert "command -v python3" in (tutorial / "find_python.sh").read_text()

    run_dir = tutorial / "rd_tutorial_input"
    with (run_dir / "wang.pkl").open("rb") as handle:
        supplied_wang = pickle.load(handle)
    assert supplied_wang["metadata"]["wl_nstep_mult"] == 4000
    assert supplied_wang["metadata"]["wl_flatness"] == 0.95
    assert supplied_wang["metadata"]["wl_tol"] == 1.0000001
    assert len(supplied_wang["td"]) == 77
    sampled = {
        "vel": {"ivel": [990.0, 1000.0, 1010.0, 1020.0]},
        "orb": {
            "sJ": [1.0, 2.0, 3.0, 4.0],
            "sL": [1.1, 2.1, 3.1, 4.1],
            "sJab": [0.2, 0.3, 0.4, 0.5],
        },
        "rot": {
            "m0": {"J": [5, 5, 4, 3], 5: {"jz": [1, 2], "sjz": [1, 2], "qjz": [1, 2]}},
            "m1": {"J": [5, 4, 5, 3], 5: {"jz": [0, 1], "sjz": [0, 1], "qjz": [1, 1]}},
        },
    }
    with (run_dir / "dat_tutorial_input.txt.pkl").open("wb") as handle:
        pickle.dump(sampled, handle)
    with (run_dir / "wang.pkl").open("wb") as handle:
        pickle.dump(
            {
                "metadata": {"format": "icats-wang-umbrella", "version": 2},
                "uu": np.ones(5),
                "iwld": np.arange(1.0, 6.0),
                "td": np.linspace(0.2, 1.0, 5),
            },
            handle,
        )

    subprocess.run(
        [sys.executable, "export_paper_histogram_data.py", "--vibrational-samples", "20"],
        cwd=tutorial,
        check=True,
    )
    metadata = json.loads((tutorial / "paper_histogram_data" / "metadata.json").read_text())
    assert metadata["accepted_samples"] == 4
    assert metadata["figure6_samples"] == 20
    projection_path = tutorial / "paper_histogram_data" / "figure7_NH3_J5_projection.csv"
    assert projection_path.exists()
    assert projection_path.read_text().splitlines()[0] == (
        "projection_centre_Jz,reconstructed_vector_Jz,quantum_rms_Jz"
    )
    assert (tutorial / "paper_histogram_data" / "figure9_wang_landau.csv").exists()

    subprocess.run([sys.executable, "plot_paper_histograms.py"], cwd=tutorial, check=True)
    for figure in (6, 7, 9):
        assert (tutorial / "paper_histogram_plots" / f"figure{figure}_validation.png").exists()
        assert (tutorial / "paper_histogram_plots" / f"figure{figure}_validation.pdf").exists()


def test_standard_tutorial_keeps_general_histogram_helpers(tmp_path):
    tutorial = tmp_path / "quickstart"
    _generate_tutorial("quickstart", tutorial)
    histograms = tutorial / "rd_tutorial_input" / "histograms"
    assert (histograms / "plot_initial.sh").exists()
    assert (histograms / "plot_sampled.sh").exists()
    assert (histograms / "plot_orbital_jljab.sh").exists()
    assert not (tutorial / "export_paper_histogram_data.py").exists()
    assert not (tutorial / "plot_paper_histograms.py").exists()
    dynamics = (tutorial / "run_cheap_dynamics.sh").read_text()
    assert "--ntraj 1 --steps 20" in dynamics
    setup = (tutorial / "setup_conda_env.sh").read_text()
    assert "ICATS.git@v0.1.3" in setup
    assert "import icats, numpy, h5py, pyscf" in setup
