import numpy as np
import pytest
from scipy.special import gammainc

from icats.ensemble_statistics import (
    distribution_statistics, discrete_statistics, uniform_cdf_statistics,
    quantum_j_probabilities, orbital_channel_probabilities, mean_statistics,
)


def test_structureless_channel_count():
    pl = np.arange(1, 22, 2, dtype=float)
    pl /= pl.sum()
    np.testing.assert_allclose(quantum_j_probabilities([1], [1], pl), pl)


def test_channel_count_matches_explicit_couplings():
    pa, pb, pl = [.2, .8], [.3, .2, .5], [.1, .3, .6]
    expected = np.zeros(6)
    for a, a_p in enumerate(pa):
        for b, b_p in enumerate(pb):
            for s in range(abs(a-b), a+b+1):
                for l, l_p in enumerate(pl):
                    for j in range(abs(l-s), l+s+1):
                        expected[j] += a_p*b_p*l_p*(2*j+1)/((2*a+1)*(2*b+1)*(2*l+1))
    np.testing.assert_allclose(quantum_j_probabilities(pa, pb, pl), expected, atol=1e-14)


def test_total_j_tail_survives():
    p = quantum_j_probabilities([0, 0, 1], [1], [0, 0, 1])
    assert p[3:].sum() > 0
    assert p.sum() == pytest.approx(1)


def test_matched_orbital_cells():
    cap = (np.sqrt(1+4*71**2)-1)/2
    np.testing.assert_allclose(orbital_channel_probabilities(cap), (2*np.arange(71)+1)/71**2)
    assert orbital_channel_probabilities(70).sum() == pytest.approx(1)


def test_statistics_keep_out_of_support_samples():
    result = discrete_statistics([0, 1, 2], [.5, .5])
    assert result["outside_reference_count"] == 1
    assert result["total_variation"] == pytest.approx(1/3)
    assert uniform_cdf_statistics([-.1, .2, 1, 1.1])["outside_reference_count"] == 2
    assert mean_statistics([1])["standard_error"] is None
    assert distribution_statistics([], [1]) == {"n": 0}


def test_leading_wigner_reference():
    rng = np.random.default_rng(123)
    r2 = rng.gamma(5, 1, 30000)
    result = uniform_cdf_statistics(gammainc(5, r2))
    assert result["total_variation"] < .025
    assert mean_statistics(r2/2, 2.5)["mean_error"] == pytest.approx(0, abs=.04)


def test_unobserved_tail_has_nonzero_reference_uncertainty():
    result = discrete_statistics(np.zeros(100, dtype=int), [.999, .001])
    assert result["moments"]["standard_error"] == 0
    assert result["moments"]["reference_standard_error"] > 0


def test_vibrational_worker_merge_preserves_modes():
    from types import SimpleNamespace
    from copy import deepcopy
    from icats.iscattering import icats
    first = {"vib": {"m0": {"vi": [[0], [1]], **{k: [[1], [2]] for k in ("Q", "sQ", "P", "sP")}}}}
    second = deepcopy(first)
    sc = SimpleNamespace(mol=[SimpleNamespace(ip=SimpleNamespace(MaxV=1))])
    merged = icats.MergeSdats(sc, [SimpleNamespace(sdat=first), SimpleNamespace(sdat=second)])
    assert merged["vib"]["m0"]["vi"] == [[0, 0], [1, 1]]
    assert first["vib"]["m0"]["vi"] == [[0], [1]]


def test_report_skips_biased_reference(tmp_path):
    from types import SimpleNamespace
    from icats.ensemble_statistics import write_ensemble_statistics
    import json
    sc = SimpleNamespace(ip=SimpleNamespace(usewang=True),
                         sdat={"orb": {"cosLJab_thet": []}})
    write_ensemble_statistics(sc, [], tmp_path)
    report = json.loads((tmp_path/"ensemble_statistics.json").read_text())
    assert "unconditional references" in report["skipped"]
