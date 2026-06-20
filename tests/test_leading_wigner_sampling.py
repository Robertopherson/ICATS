import numpy as np

from icats.dist import HusimiFuncICDF as LeadingWignerFuncICDF
from icats.mc import ICDFsample, InitICDF


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
