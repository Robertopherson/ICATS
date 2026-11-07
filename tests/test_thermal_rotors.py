from types import SimpleNamespace

import numpy as np
import pytest

from icats.constants import kboltz
from icats.dist import RotorBoltzmannProbabilities, RotorTotBoltzICDF, SymRotorProjBoltzICDF, JcrossICDFc, IsotropicDistICDF
from icats.mc import InitICDF, ICDFscalar
from icats.molecules import imolecule
from icats import wang


def molecule(kind, a=1e-4, c=0.0, maximum=25, rr=None):
    mol = imolecule()
    mol.sp = SimpleNamespace(rsym=kind, J2c=a, Jzc=c, RR=rr,
                             rA=a, rB=a, rC=a+c)
    mol.ip.MaxR = maximum
    return mol


@pytest.mark.parametrize("kind,power", [("linear", 1), ("spherical", 2)])
def test_multiplicity(kind, power):
    mol = molecule(kind)
    j = np.arange(mol.ip.MaxR+1)
    expected = (2*j+1)**power*np.exp(-mol.sp.J2c*j*(j+1)/(kboltz*500))
    np.testing.assert_allclose(mol.ThermalRotorProbabilities(500), expected/expected.sum())


@pytest.mark.parametrize("c", [-5e-5, 0.0, 2e-4])
def test_symmetric_joint_distribution(c):
    mol = molecule("oblate" if c < 0 else "prolate", c=c, maximum=10)
    p = mol.ThermalRotorProbabilities(500)
    joint = []
    expected = []
    for j in range(len(p)):
        k = np.arange(-j, j+1)
        conditional = np.exp(-c*k*k/(kboltz*500))
        joint.extend(p[j]*conditional/conditional.sum())
        expected.extend((2*j+1)*np.exp(-(mol.sp.J2c*j*(j+1)+c*k*k)/(kboltz*500)))
    expected = np.asarray(expected)
    np.testing.assert_allclose(joint, expected/expected.sum(), rtol=1e-13)


def test_asymmetric_sublevel_sum():
    a, b, c = 3e-4, 1.4e-4, 0.9e-4
    base = (a+b+c)/3
    rr, full = [], []
    for j in range(13):
        k = np.arange(-j, j+1)
        jp = np.diag(np.sqrt(j*(j+1)-k[:-1]*(k[:-1]+1)), -1)
        jx, jy = (jp+jp.T)/2, (jp-jp.T)/(2j)
        h = a*(jx@jx)+b*(jy@jy)+c*np.diag(k*k)
        energies = np.linalg.eigvalsh(h)
        full.append(energies)
        rr.append([np.arange(2*j+1), energies-base*j*(j+1)])
    mol = molecule("asym-spherical", a=base, maximum=12, rr=rr)
    weights = np.array([(2*j+1)*np.exp(-ee/(kboltz*500)).sum() for j, ee in enumerate(full)])
    np.testing.assert_allclose(mol.ThermalRotorProbabilities(500), weights/weights.sum())


@pytest.mark.parametrize("c", [-5e-5, 0.0, 2e-4])
def test_discrete_k_samples(c):
    j, n = 4, 40000
    rng = [np.random.default_rng(164)]
    samples = np.array([SymRotorProjBoltzICDF(rng, 500, c, j)[0] for _ in range(n)])
    k = np.arange(-j, j+1)
    prob = np.exp(-c*k*k/(kboltz*500))
    prob /= prob.sum()
    measured = np.bincount(samples+j, minlength=len(k))/n
    assert np.all(np.abs(measured-prob) < 6*np.sqrt(prob*(1-prob)/n)+1/n)


def test_j_sampler_and_restart_state():
    p = molecule("spherical").ThermalRotorProbabilities(500)
    first = InitICDF(1, RotorTotBoltzICDF, [p], seed=12)
    second = InitICDF(1, RotorTotBoltzICDF, [p], seed=12)
    assert [ICDFscalar(first) for _ in range(100)] == [ICDFscalar(second) for _ in range(100)]
    rng = [np.random.default_rng(87)]
    n = 40000
    samples = np.array([RotorTotBoltzICDF(rng, p)[0] for _ in range(n)])
    assert samples.min() >= 0 and samples.max() <= 25
    measured = np.bincount(samples, minlength=len(p))/n
    assert np.all(np.abs(measured-p) < 6*np.sqrt(p*(1-p)/n)+2/n)


def test_cold_rotor_is_stable():
    p = molecule("oblate", c=-5e-5).ThermalRotorProbabilities(1e-5)
    assert p[0] == 1 and np.isfinite(p).all()
    for c in (-1e-4, 1e-4):
        k, weight = SymRotorProjBoltzICDF([np.random.default_rng(5)], 1e-5, c, 10)
        assert abs(k) == (10 if c < 0 else 0)
        assert np.isfinite(weight)


def test_conservative_cutoff_independent_of_polarization():
    mol = molecule("oblate", c=-5e-5)
    mol.ip.isotropic = True
    cutoff = mol.EstimateMaxR(500)
    mol.ip.isotropic = False
    assert cutoff == mol.EstimateMaxR(500)
    mol.ip.MaxR = cutoff+30
    assert mol.ThermalRotorProbabilities(500)[cutoff+1:].sum() < 1e-12


def test_old_umbrella_cannot_be_reused():
    with pytest.raises(ValueError, match="predates canonical"):
        wang.validate("old.pkl", [], [], np.ones(5), None,
                      {"maxj": 5, "rotor_sampling": "canonical-sublevels-v1"})


@pytest.mark.parametrize("temperature", [0, -1, np.inf, np.nan])
def test_invalid_thermal_temperature(temperature):
    with pytest.raises(ValueError):
        RotorBoltzmannProbabilities(temperature, [[0]])


def test_orbital_inverse_cdf():
    first, reference = np.random.default_rng(123), np.random.default_rng(123)
    for cap in (1, 5, 70):
        for _ in range(100):
            ell, _ = JcrossICDFc([first], True, cap)
            np.testing.assert_allclose(ell*(ell+1)/(cap*(cap+1)), reference.random(), atol=1e-14)


def test_spherical_vector_generation_does_not_require_k_sampler():
    mol = molecule("spherical")
    mol.sp.Irep = 3
    cont = InitICDF(1, IsotropicDistICDF,
                    ["eul", np.array([2*np.pi, 1, 2*np.pi]), np.array([np.pi, 0, np.pi])], seed=16)
    sa = SimpleNamespace(sJ=5, SampInfo={}, dist={"rotJ": {}, "arotJz": {"cont": cont}})
    mol.SampleRigidRotorState(sa)
    np.testing.assert_allclose(np.dot(sa.srpar[-1], sa.srpar[-1]), 30)
    np.testing.assert_allclose(sa.srpar[1], mol.sp.J2c*30)
