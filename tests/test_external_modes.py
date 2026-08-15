import numpy as np

from icats.functions import GetRotTransVec


def test_external_mode_count_and_orthonormality():
    cases = [
        (np.array([[0.0, 0.0, 0.0]]), np.array([4.0]), 3),
        (
            np.array([[0.0, 0.0, -0.55], [0.0, 0.0, 0.55]]),
            np.array([14.0, 14.0]),
            5,
        ),
        (
            np.array(
                [[0.0, 0.0, 0.0], [0.76, 0.59, 0.0], [-0.76, 0.59, 0.0]]
            ),
            np.array([16.0, 1.0, 1.0]),
            6,
        ),
    ]

    for coordinates, masses, expected_modes in cases:
        modes = GetRotTransVec(coordinates, masses, [])
        assert modes.shape == (expected_modes, coordinates.size)
        np.testing.assert_allclose(
            modes @ modes.T, np.eye(expected_modes), atol=1.0e-12
        )


def test_linear_molecule_retains_one_vibrational_coordinate():
    coordinates = np.array([[0.0, 0.0, -0.55], [0.0, 0.0, 0.55]])
    modes = GetRotTransVec(coordinates, np.array([14.0, 14.0]), [])

    assert coordinates.size - modes.shape[0] == 1
