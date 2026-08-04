import numpy as np

from icats.dist import HusimiFuncICDF as LeadingWignerFuncICDF, uniform
from icats.mc import ICDFsample, ICDFscalar, InitICDF


def test_scalar_icdf_returns_python_scalar():
    cont = InitICDF(1, uniform, [0.0, 1.0], seed=8123)
    sampled = ICDFscalar(cont)
    assert isinstance(sampled, float)
    assert 0.0 <= sampled <= 1.0


def test_leading_wigner_icdf_preserves_oscillator_energy_moment():
    """The sampled radial moment should reproduce E/w = n + 1/2."""
    nsamp = 20000
    for n in range(5):
        cont = InitICDF(1, LeadingWignerFuncICDF, [n], seed=12007 + n)
        energies = np.empty(nsamp)
        for i in range(nsamp):
            q, p = ICDFsample(cont)[0]
            energies[i] = 0.5 * (q * q + p * p)
        assert abs(float(energies.mean()) - (n + 0.5)) < 0.06
