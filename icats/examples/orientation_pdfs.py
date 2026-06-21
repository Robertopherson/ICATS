#!/usr/bin/env python3
"""Example molecular orientation PDFs for ICATS tutorials.

Functions receive molecular Euler angles in radians in the ICATS scattering
frame and return a non-negative physical density. ICATS applies the Euler
measure sin(beta) during sampling.
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
