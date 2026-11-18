"""Projection smoothing with a minimum-relative-entropy moment correction."""

from functools import lru_cache

import numpy as np
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.optimize import brentq


@lru_cache(maxsize=128)
def projection_cdf(j, probabilities, width=0.8):
    """Tabulate the corrected Gaussian mixture; cache by eigenstate and width.

    The constraint is on the whole mixture, not separately on each K peak.
    Zero width and zero second moment are handled directly by the sampler.
    """
    p = np.asarray(probabilities, dtype=float)
    if j < 1 or int(j) != j or p.shape != (2*int(j)+1,):
        raise ValueError("K broadening requires integer J and 2J+1 probabilities")
    if not np.isfinite(p).all() or np.any(p < 0) or p.sum() <= 0:
        raise ValueError("Invalid projection probabilities")
    if not np.isfinite(width) or width <= 0:
        raise ValueError("CDF width must be positive and finite")
    p = p/p.sum()
    k = np.arange(-j, j+1)
    radius = np.sqrt(j*(j+1))
    target = np.dot(p, k*k)/radius**2
    if target <= 0:
        raise ValueError("Zero second moment is a point mass, not a continuous CDF")
    # Resolve narrow peaks as well as the default width; fail rather than
    # silently replacing an unresolved continuous distribution by sharp K.
    count = max(32769, int(np.ceil(64*radius/width))+1)
    if count > 262145:
        raise ValueError("K width too small to resolve; use width 0 for discrete K")
    u = np.linspace(-1, 1, count)
    logbase = np.full(count, -np.inf)
    for peak, weight in zip(k, p):
        if weight > 0:
            logbase = np.logaddexp(logbase, np.log(weight)-0.5*((radius*u-peak)/width)**2)

    def density(lam):
        exponent = logbase+lam*u*u
        return np.exp(exponent-exponent.max())

    def residual(lam):
        d = density(lam)
        return trapezoid(d*u*u, u)/trapezoid(d, u)-target

    lo, hi = -1., 1.
    for _ in range(60):
        if residual(lo) <= 0 <= residual(hi):
            break
        lo *= 2
        hi *= 2
    else:
        raise RuntimeError("Cannot bracket the K second-moment correction")
    lam = brentq(residual, lo, hi, xtol=1e-11)
    d = density(lam)
    cdf = cumulative_trapezoid(d, u, initial=0)
    cdf /= cdf[-1]
    grid = radius*u
    grid.setflags(write=False)
    cdf.setflags(write=False)
    return grid, cdf, lam


@lru_cache(maxsize=512)
def interval_cdf(j, k):
    """Uniform prior in a positive K interval, tilted to preserve <Kc^2>=K^2.

    Negative projections use the mirrored kernel; K=0 stays a point mass.
    The scaled, centred statistic avoids cancellation at large J.
    """
    if int(j) != j or int(k) != k or not 1 <= k <= j:
        raise ValueError("Interval CDF requires integer 1 <= K <= J")
    radius = np.sqrt(j*(j+1))
    grid = np.linspace(k-.5, min(k+.5, radius), 8193)
    statistic = (grid-k)*(grid+k)/k

    def density(lam):
        exponent = lam*statistic
        return np.exp(exponent-exponent.max())

    def residual(lam):
        d = density(lam)
        return trapezoid(d*statistic, grid)/trapezoid(d, grid)

    lo, hi = -1., 1.
    for _ in range(60):
        if residual(lo) <= 0 <= residual(hi):
            break
        lo *= 2
        hi *= 2
    else:
        raise RuntimeError("Cannot bracket the interval moment correction")
    lam = brentq(residual, lo, hi, xtol=1e-12)
    cdf = cumulative_trapezoid(density(lam), grid, initial=0)
    cdf /= cdf[-1]
    grid.setflags(write=False)
    cdf.setflags(write=False)
    return grid, cdf, lam


def RotorProjectionICDF(rng, width, j, probabilities, method='interval'):
    """Sample a parent K, then its corrected interval; width=0 keeps sharp K.

    Positive width is only used by the retained Gaussian comparison method.
    """
    p = np.asarray(probabilities, dtype=float)
    if int(j) != j or j < 0 or p.shape != (2*int(j)+1,):
        raise ValueError("Expected integer J and 2J+1 projection probabilities")
    if not np.isfinite(width) or width < 0:
        raise ValueError("rotor-k-width must be finite and nonnegative")
    if not np.isfinite(p).all() or np.any(p < 0) or p.sum() <= 0:
        raise ValueError("Invalid projection probabilities")
    p = p/p.sum()
    k = np.arange(-j, j+1)
    if method not in ('interval', 'gaussian'):
        raise ValueError("Unknown projection smoothing method")
    if np.dot(p, k*k) == 0:
        value = 0.
    elif width == 0 or method == 'interval':
        parent = int(rng[0].choice(k, p=p))
        value = float(parent)
        if parent != 0 and width != 0:
            grid, cdf, _ = interval_cdf(int(j), abs(parent))
            value = np.sign(parent)*float(np.interp(rng[0].random(), cdf, grid))
    else:
        grid, cdf, _ = projection_cdf(int(j), tuple(p), float(width))
        value = float(np.interp(rng[0].random(), cdf, grid))
    return np.float64(value), 1.


def RotorGaussianProjectionICDF(rng, width, j, probabilities):
    """Historical Gaussian-mixture law for explicitly requested comparisons."""
    return RotorProjectionICDF(rng, width, j, probabilities, method='gaussian')
