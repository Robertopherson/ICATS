"""Discrete state selection before continuous rotor-direction sampling."""
import numpy as np
from scipy.special import logsumexp

from .constants import kboltz


def configured(ip):
    return (getattr(ip, 'rotor_state', None) is not None or
            getattr(ip, 'rotor_k_rule', 'all') != 'all' or
            getattr(ip, 'rotor_allowed_k', None) is not None or
            getattr(ip, 'rotor_k_weights', None) is not None)


def residue_weights(ip):
    """One relative weight per signed-K residue; zero excludes that class."""
    values = getattr(ip, 'rotor_k_weights', None)
    if values is None:
        return None
    w = np.asarray(values, dtype=float)
    if w.ndim != 1 or not len(w) or not np.isfinite(w).all() or np.any(w < 0) or not np.any(w > 0):
        raise ValueError('rotor-k-weights requires finite nonnegative weights, at least one positive')
    if not np.allclose(w, w[(-np.arange(len(w))) % len(w)], rtol=1e-12, atol=0):
        raise ValueError('rotor-k-weights must give equal weights to +K and -K')
    if getattr(ip, 'rotor_k_rule', 'all') != 'all' or getattr(ip, 'rotor_allowed_k', None) is not None:
        raise ValueError('Do not combine rotor-k-weights with NH3 presets or rotor-allowed-k')
    return w


def state_tables(mol, temperature):
    """Return P(J) and conditional K/eigenstate weights; NH3 inversion is unresolved."""
    sp, ip = mol.sp, mol.ip
    residues = residue_weights(ip)
    asymmetric = 'asym' in sp.rsym
    if sp.rsym not in ('linear', 'prolate', 'oblate') and not (asymmetric and residues is not None):
        raise ValueError('rotor-state and K restrictions require a linear or symmetric top')
    fixed = getattr(ip, 'rotor_state', None)
    rule = getattr(ip, 'rotor_k_rule', 'all')
    allowed = getattr(ip, 'rotor_allowed_k', None)
    if asymmetric and fixed is not None:
        raise ValueError('rotor-state J,K requires a linear or symmetric top, not an asymmetric eigenstate')
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
        if residues is not None:
            wk = residues[ks % len(residues)]
            if asymmetric:
                probabilities = np.asarray(sp.rotor_projection_probabilities[j])
                # Do not project/truncate an energy eigenvector to impose a rule.
                w = np.empty(2*j+1)
                for state in range(2*j+1):
                    support = wk[probabilities[:, state] > 1e-10]
                    if not len(support) or not np.allclose(support, support[0], rtol=1e-12, atol=0):
                        raise ValueError(f'Incompatible rotor-k-weights: J={j}, state={state} mixes differently weighted K residues')
                    w[state] = support[0]
            else:
                w *= wk
        # Signed K already counts nonzero |K| twice versus K=0 once.
        # For unresolved ortho NH3 this supplies the inversion count;
        # the common spin multiplicity cancels within the species.
        weights.append(w)
        projection = np.asarray(sp.RR[j][1]) if asymmetric else sp.Jzc*ks**2
        energies.append(sp.J2c*j*(j+1)+projection)
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


def eigenstate_icdf(rng, tables, rr, j):
    p = tables[int(j)]
    if not np.isclose(np.sum(p), 1):
        raise ValueError('Cannot sample a forbidden J manifold')
    state = rng[0].choice(len(p), p=p)
    return [column[state] for column in rr[int(j)]], 1.
