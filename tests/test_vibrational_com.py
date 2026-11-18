from types import SimpleNamespace

import numpy as np
import pytest

from icats.molecules import imolecule


@pytest.mark.parametrize("angle", [0.0, 0.7])
def test_vibration_preserves_molecular_com_velocity(angle):
    mass = np.array([16., 1., 1.])
    x = np.array([[0., 0., 0.1], [0.8, 0., -0.8], [-0.3, 0.7, -0.8]])
    x -= np.average(x, axis=0, weights=mass)
    r = np.array([[np.cos(angle), -np.sin(angle), 0.],
                  [np.sin(angle), np.cos(angle), 0.], [0., 0., 1.]])
    # Small translational contamination must not alter the prepared beam velocity.
    du = np.array([[.01, .02, .03], [-.02, .03, .01], [.04, -.01, .02]])
    p = du * mass[:, None]
    mol = SimpleNamespace(sp=SimpleNamespace(mass=mass, na=3, xxe=x,
                          n2x=np.zeros((1, 9)), n2p=p.reshape(1, 9)), ip=None)
    v = np.array([.2, -.3, .4])
    sa = SimpleNamespace(sxx=x @ r.T, svv=np.tile(v, (3, 1)))
    imolecule.SetHOVibrState(mol, sa, np.array([0.]), np.array([1.]))
    np.testing.assert_allclose(np.average(sa.svv, axis=0, weights=mass), v, atol=1e-14)
    # Removing a uniform translation must leave relative vibrational motion intact.
    np.testing.assert_allclose(sa.svv[1] - sa.svv[0], (du[1] - du[0]) @ r.T, atol=1e-14)
