"""Orbital-channel quadrature in atomic units; no quantum propagation implied."""
import math


def channel_impact(channel, momentum, rule, uniform_draw=0.5):
    """Return b and its shell-integration weight (bohr, bohr squared)."""
    if not math.isfinite(channel) or channel < 0 or channel != int(channel):
        raise ValueError('orbital-channel must be a nonnegative integer')
    if not math.isfinite(momentum) or momentum <= 0:
        raise ValueError('orbital-channel requires positive incoming momentum')
    ell = int(channel)
    area = math.pi*(2*ell+1)/momentum**2
    if rule == 'sqrt':
        return math.sqrt(ell*(ell+1))/momentum, area
    if rule == 'langer':
        return (ell+0.5)/momentum, area
    if rule == 'uniform-b':
        if not 0 <= uniform_draw <= 1:
            raise ValueError('uniform draw must lie in [0,1]')
        b = (ell+uniform_draw)/momentum
        return b, 2*math.pi*b/momentum
    raise ValueError('orbital-channel-rule must be sqrt, langer or uniform-b')
