import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import i0e
from icats.dist import RotorAzimuthVonMises, _azimuth_concentration


@pytest.mark.parametrize('a,b', [(1.,1.), (1.,4.), (20.,1.), (1.,1000.)])
def test_azimuth_moments(a, b):
    r = (a-b)/(a+b)
    k = _azimuth_concentration(abs(r))
    centre = 0 if r >= 0 else np.pi
    density = lambda g: np.exp(k*(np.cos(2*g-centre)-1))/(2*np.pi*i0e(k))
    points = np.linspace(-np.pi, np.pi, 33)
    assert quad(density, -np.pi, np.pi, points=points, limit=200, epsabs=1e-10)[0] == pytest.approx(1)
    assert quad(lambda g: np.cos(g)**2*density(g), -np.pi, np.pi,
                points=points, limit=200, epsabs=1e-10)[0] == pytest.approx(a/(a+b), abs=1e-10)


@pytest.mark.parametrize('a,b', [(1.,1.), (1.,4.), (20.,1.), (0.,1.), (1.,0.), (0.,0.)])
def test_azimuth_draws(a, b):
    rng = [np.random.default_rng(154)]
    angles = np.array([RotorAzimuthVonMises(rng, a, b)[0][0] for _ in range(30000)])
    target = a/(a+b) if a+b else .5
    assert np.mean(np.cos(angles)**2) == pytest.approx(target, abs=.012)
    assert abs(np.mean(np.cos(angles))) < .02
    assert abs(np.mean(np.sin(angles))) < .02


def test_invalid_moments():
    with pytest.raises(ValueError):
        RotorAzimuthVonMises([np.random.default_rng(1)], -1, 2)
