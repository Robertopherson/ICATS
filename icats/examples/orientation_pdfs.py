#!/usr/bin/env python3
"""Example molecular orientation PDFs for ICATS tutorials.

Each function receives molecular Euler angles (alpha, beta, gamma) in radians
in the ICATS scattering frame and returns a non-negative physical density.
ICATS applies the sin(beta) Euler measure during sampling, so the function
should return the density *without* that factor.

To supply your own PDF, copy this file to your working directory, add a
function with the same signature, and reference it in the ICATS input::

    orientation-mode = pdf orientation_pdfs.py my_function [param1 param2 ...]

Extra positional arguments after the function name are passed as strings and
should be cast to float inside the function.
"""
import numpy as np


def _unit(v):
    v = np.asarray(v, dtype=float)
    n = np.linalg.norm(v)
    if n <= 0.0:
        raise ValueError("zero-length axis in orientation PDF")
    return v / n


def _rz(alpha):
    c, s = np.cos(alpha), np.sin(alpha)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def _ry(beta):
    c, s = np.cos(beta), np.sin(beta)
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])


def _rotation(alpha, beta, gamma):
    return _rz(alpha) @ _ry(beta) @ _rz(gamma)


def dipole_field_z(alpha, beta, gamma, strength=0.6):
    """Dipole orientation PDF for a body-z dipole and field along ICATS Z.

    P = 1 + A cos(theta_muE).  Use |A| <= 1 to keep the PDF non-negative.
    """
    strength = float(strength)
    return max(0.0, 1.0 + strength * np.cos(beta))


def dipole_field_tilted(alpha, beta, gamma, strength=0.6, field_theta=0.5, field_phi=0.0):
    """Dipole orientation PDF for a tilted field in the ICATS scattering frame.

    The molecular dipole is assumed to lie along body z.  The field direction
    is specified by spherical angles (field_theta, field_phi) in the ICATS
    scattering frame.  This generally depends on both alpha and beta.
    """
    strength = float(strength)
    field_theta = float(field_theta)
    field_phi = float(field_phi)
    field = _unit(
        [
            np.sin(field_theta) * np.cos(field_phi),
            np.sin(field_theta) * np.sin(field_phi),
            np.cos(field_theta),
        ]
    )
    mu_lab = _rotation(alpha, beta, gamma) @ np.array([0.0, 0.0, 1.0])
    cos_mu_e = float(np.dot(mu_lab, field))
    return max(0.0, 1.0 + strength * cos_mu_e)


def alignment_field_z(alpha, beta, gamma, strength=0.6):
    """Laser-like alignment example, P = 1 + A P2(cos beta)."""
    strength = float(strength)
    p2 = 0.5 * (3.0 * np.cos(beta) ** 2 - 1.0)
    return max(0.0, 1.0 + strength * p2)


def _wigner_d_small(j, m, k, beta):
    """Wigner small-d matrix element d^j_{m,k}(beta) for integer or half-integer j."""
    from scipy.special import gammaln
    j, m, k = float(j), float(m), float(k)
    cb2 = np.cos(beta / 2.0)
    sb2 = np.sin(beta / 2.0)
    half_log_norm = 0.5 * (
        gammaln(j + m + 1) + gammaln(j - m + 1)
        - gammaln(j + k + 1) - gammaln(j - k + 1)
    )
    s_min = int(round(max(0.0, k - m)))
    s_max = int(round(min(j - m, j + k)))
    result = 0.0
    for s in range(s_min, s_max + 1):
        p1 = 2 * j - 2 * s - m + k
        p2 = 2 * s + m - k
        if p1 < 0 or p2 < 0:
            continue
        log_binom = (
            gammaln(j + k - s + 1) - gammaln(s + 1)
            - gammaln(j - m - s + 1) - gammaln(m - k + s + 1)
        )
        sign = (-1) ** int(round(m - k + s))
        result += sign * np.exp(log_binom + half_log_norm) * (cb2 ** p1) * (sb2 ** p2)
    return result


def rotor_state(alpha, beta, gamma, j=1, m=0, k=0):
    """State-selected symmetric-top (or linear-rotor) orientation PDF.

    Returns |d^j_{m,k}(beta)|^2, the body-axis angular distribution for a
    molecule prepared in a definite (j, m, k) rotational eigenstate with the
    quantization axis along the ICATS Z-axis (relative velocity direction).

    Parameters
    ----------
    j : int or half-integer float
        Total angular momentum quantum number.
    m : int or half-integer float
        Space-fixed projection along the quantization axis.
    k : int or half-integer float
        Body-fixed projection (0 for a linear molecule).

    Examples
    --------
    Linear molecule in |j=2, m=0>::

        orientation-mode = pdf orientation_pdfs.py rotor_state 2 0 0

    Symmetric top in |j=2, m=2, k=2>::

        orientation-mode = pdf orientation_pdfs.py rotor_state 2 2 2
    """
    j, m, k = float(j), float(m), float(k)
    d = _wigner_d_small(j, m, k, beta)
    return float(d * d)


def hexapole_no(alpha, beta, gamma, m=0.5, field_theta=0.0, field_phi=0.0):
    """Hexapole-selected NO orientation PDF for the X2Pi ground state.

    For NO X(2)Pi (j=1/2, Omega=1/2), an electrostatic hexapole state-selector
    transmits only |M|=1/2, quantized along the orienting electric field.
    The body-axis angular distribution is the Wigner-d probability:

        P = |d^{1/2}_{M, 1/2}(theta_E)|^2 = (1 +/- cos theta_E) / 2

    where theta_E is the angle between the molecular body axis and the
    orienting field.

    Parameters
    ----------
    m : float
        Space-fixed projection along the field.  Use m=+0.5 for N-end-first
        orientation (N-atom toward the field) or m=-0.5 for O-end-first.
    field_theta, field_phi : float
        Direction of the orienting electric field in the ICATS scattering
        frame F, in radians.  field_theta=0 (default) puts the field along
        the ICATS Z-axis (relative velocity).  In a real experiment the
        field direction in F must be computed from the beam geometry (see
        Eqn. for phi in the manuscript); passing field_theta != 0 makes the
        distribution depend on both alpha and beta.

    With field_theta=0 the PDF reduces to p(cos beta) = (1+-cos beta)/2 and
    depends only on beta (1-D); with field_theta != 0 it is genuinely 2-D.

    Usage in an ICATS data file::

        # field along beam axis (1-D):
        orientation-mode = pdf orientation_pdfs.py hexapole_no  0.5
        orientation-mode = pdf orientation_pdfs.py hexapole_no -0.5
        # field tilted 45 deg from beam axis in F (2-D):
        orientation-mode = pdf orientation_pdfs.py hexapole_no  0.5 0.7854
    """
    m = float(m)
    field_theta = float(field_theta)
    field_phi = float(field_phi)
    if abs(field_theta) < 1e-9:
        cos_e = np.cos(beta)
    else:
        field = _unit([
            np.sin(field_theta) * np.cos(field_phi),
            np.sin(field_theta) * np.sin(field_phi),
            np.cos(field_theta),
        ])
        mol_axis = np.array([
            np.sin(beta) * np.cos(alpha),
            np.sin(beta) * np.sin(alpha),
            np.cos(beta),
        ])
        cos_e = float(np.dot(field, mol_axis))
    sign = 1.0 if m > 0.0 else -1.0
    return max(0.0, (1.0 + sign * cos_e) / 2.0)
