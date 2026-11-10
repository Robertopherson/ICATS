"""Discrete state selection before continuous rotor-direction sampling."""
import numpy as np
from scipy.special import logsumexp

from .constants import kboltz


def configured(ip):
    return (getattr(ip, 'rotor_state', None) is not None or
            getattr(ip, 'rotor_k_rule', 'all') != 'all' or
            getattr(ip, 'rotor_allowed_k', None) is not None)


def state_tables(mol, temperature):
    """Return P(J), P(signed K|J); NH3 inversion splitting is neglected."""
    sp, ip = mol.sp, mol.ip
    if sp.rsym not in ('linear', 'prolate', 'oblate'):
        raise ValueError('rotor-state and K restrictions require a linear or symmetric top')
    fixed = getattr(ip, 'rotor_state', None)
    rule = getattr(ip, 'rotor_k_rule', 'all')
    allowed = getattr(ip, 'rotor_allowed_k', None)
    if rule not in ('all', 'nh3-para', 'nh3-ortho'):
        raise ValueError('rotor-k-rule must be all, nh3-para or nh3-ortho')
    if rule != 'all' and sorted(str(e).upper() for e in sp.el) != ['H','H','H','N']:
        raise ValueError('NH3 K rules require NH3; they are not generic isotope/spin rules')
    if not np.isfinite(temperature) or temperature < 0:
        raise ValueError('Rotational temperature must be finite and nonnegative')
    maximum = int(ip.MaxR)
    if fixed is not None:
        j, k = fixed
        if j < 0 or k < 0 or k > j or int(j) != j or int(k) != k:
            raise ValueError('rotor-state requires integer J >= |K| >= 0')
        maximum = max(maximum, j)
    energies, weights = [], []
    for j in range(maximum+1):
        ks = np.arange(-j, j+1)
        keep = np.ones(len(ks), dtype=bool)
        if sp.rsym == 'linear':keep &= ks == 0
        if rule == 'nh3-para':keep &= ks % 3 != 0
        if rule == 'nh3-ortho':keep &= ks % 3 == 0
        if allowed is not None:keep &= np.isin(abs(ks), allowed)
        if fixed is not None:keep &= (j == fixed[0]) & (abs(ks) == fixed[1])
        w = keep.astype(float)
        # K=0 has one allowed inversion component; |K|>0 has two.
        # The common spin multiplicity cancels within a specified species.
        if rule == 'nh3-ortho':w[ks == 0] *= .5
        weights.append(w)
        energies.append(sp.J2c*j*(j+1)+sp.Jzc*ks**2)
    valid = [ee[ww > 0] for ee, ww in zip(energies, weights) if np.any(ww > 0)]
    if not valid:
        raise ValueError('No rotational states satisfy the requested state/K restrictions')
    ground = min(ee.min() for ee in valid)
    logs=[]
    for ee, ww in zip(energies, weights):
        logw=np.full(len(ww), -np.inf)
        keep=ww>0
        if temperature == 0:keep &= np.isclose(ee,ground,rtol=1e-12,atol=1e-15)
        logw[keep]=np.log(ww[keep])
        if temperature > 0:logw[keep] -= (ee[keep]-ground)/(kboltz*temperature)
        logs.append(logw)
    # Fixed-level preparation overrides the thermal distribution.
    if fixed is not None:
        logs=[np.where(w>0,0.,-np.inf) for w in weights]
    marg=np.array([logsumexp(v)+np.log(2*j+1) for j,v in enumerate(logs)])
    pj=np.exp(marg-logsumexp(marg))
    pk=[np.exp(v-logsumexp(v)) if np.isfinite(v).any() else np.zeros(len(v)) for v in logs]
    return pj, pk


def projection_icdf(rng, tables, j):
    p=tables[int(j)]
    if not np.isclose(np.sum(p),1):
        raise ValueError('Cannot sample a forbidden J manifold')
    return np.int64(rng[0].choice(len(p), p=p))-int(j), 1.
