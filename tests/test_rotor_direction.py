import numpy as np
import pytest

from icats.dist import reject_bingham
from icats.functions import iq2R


@pytest.mark.parametrize("variance", [[2., 8., 20.], [1., 1., 1.], [0., 1., 0.]])
def test_normalized_gaussian_direction_roundtrip(variance):
    reference = np.random.default_rng(712)
    sampled = np.random.default_rng(712)
    for _ in range(1000):
        direction = reference.multivariate_normal(np.zeros(3), np.diag(variance))
        direction /= np.linalg.norm(direction)
        angles, weight = reject_bingham([sampled], variance)
        rebuilt = iq2R(angles).T @ np.array([0., 0., 1.])
        np.testing.assert_allclose(rebuilt, direction, atol=2e-13, rtol=0)
        assert weight == 1
