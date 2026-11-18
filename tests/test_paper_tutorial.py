import pytest

from icats.cli import RETIRED_TUTORIALS, TUTORIAL_ORDER, _generate_tutorial


@pytest.mark.parametrize('name', sorted(RETIRED_TUTORIALS))
def test_retired_tutorials_fail_without_writing(tmp_path, name):
    destination = tmp_path / name
    assert name not in TUTORIAL_ORDER
    with pytest.raises(SystemExit, match='retired Wang-Landau'):
        _generate_tutorial(name, destination)
    assert not destination.exists()


def test_revised_tutorial_uses_current_sampling(tmp_path):
    destination = tmp_path / 'revised'
    _generate_tutorial('revision_nh3_h2o', destination, nsamp_override=4)
    text = (destination / 'tutorial_input.txt').read_text()
    assert 'Nsamp = 4' in text
    assert 'wang = False' in text
    assert 'orbital-sampling = geometric' in text
    assert 'beam-angle = 90.0' in text
    notes = (destination / 'tutorial_README.txt').read_text()
    assert 'within-K-interval' in notes
    assert 'K width 0.8' not in notes
    assert not list(destination.rglob('wang.pkl'))


def test_standard_tutorial_keeps_general_histogram_helpers(tmp_path):
    tutorial = tmp_path / 'quickstart'
    _generate_tutorial('quickstart', tutorial)
    histograms = tutorial / 'rd_tutorial_input' / 'histograms'
    assert (histograms / 'plot_initial.sh').exists()
    assert (histograms / 'plot_sampled.sh').exists()
    assert (histograms / 'plot_orbital_jljab.sh').exists()
    assert not (tutorial / 'export_paper_histogram_data.py').exists()
    dynamics = (tutorial / 'run_cheap_dynamics.sh').read_text()
    assert '--ntraj 1 --steps 20' in dynamics
    setup = (tutorial / 'setup_conda_env.sh').read_text()
    assert 'import icats, numpy, h5py, pyscf' in setup
    assert 'ICATS_SOURCE' in setup
    assert 'pip install "$ICATS_SOURCE"' in setup
    assert 'revision/no-ar-handoff' not in setup
