#!/usr/bin/env python3
import sys
import os
import re
import numpy as np
from numpy.linalg import norm, eigh, det, svd, pinv, inv
from numpy import reshape, matmul, array, zeros, identity
from numpy import random, cross, prod, diag, copy
from numpy import cos, sin, arccos, arcsin, arctan2, sqrt, exp
from numpy.linalg import qr, eigh, norm, eig, svd, solve
from math  import acos
import pathlib, textwrap, json, inspect
from .constants import *
import matplotlib.pyplot as plt
from copy import deepcopy
from collections.abc import Mapping, MutableSequence

np.set_printoptions(suppress=True)
_HIST_BASE_DIR = pathlib.Path("histograms")

def get_size(x):
    if isinstance(x, np.ndarray):
        if x.shape == ():  # zero-dimensional scalar array
            return 1
        else:
            return x.size
    else:
        # Assume plain Python scalar
        return 1

def block_diag(blocks):
    n = len(blocks)
    grid = []
    for i in range(n):
        row = []
        for j in range(n):
            if i == j:
                row.append(blocks[i])
            else:
                row.append(
                    np.zeros((blocks[i].shape[0], blocks[j].shape[1]),
                             dtype=blocks[i].dtype)
                )
        grid.append(row)
    return np.block(grid)

def svdsolve(A,b):
  U, s, VT = svd(A, full_matrices=False)
  tol = (np.finfo(float).eps * max(A.shape)) * s[0]
  r = np.sum(s > tol)  
  x = (VT[:r].T @ ( (U[:, :r].T @ b) / s[:r] ))
  return x 

def gamma_mt(n):
    """Generate one sample from Gamma(n+1, 1) using Marsaglia & Tsang's method."""
    k = n+1
    d = k - 1.0 / 3.0
    c = 1.0 / np.sqrt(9.0 * d)
    while True:
        x = np.random.normal()
        v = (1 + c * x)**3
        if v <= 0:
            continue
        u = np.random.rand()
        if u < 1 - 0.0331 * (x**4):
            return d * v
        if np.log(u) < 0.5 * x**2 + d * (1 - v + np.log(v)):
            return d * v


def rot_match_vec(r1, r2):
    """
    Rotate the vector r1 to match the direction of vector r2.

    Parameters:
    r1 (numpy.ndarray): The input vector to be rotated.
    r2 (numpy.ndarray): The target vector to match.

    Returns:
    numpy.ndarray: The rotated vector matching the direction of r2.
    numpy.ndarray: The rotation matrix used for the transformation.
    """
    U, s, V = svd(matmul(r1.T, r2))
    if det(matmul(U, V.T)) < 0:
        ON = identity(3)
        ON[2, 2] = -1
        RO = matmul(U, matmul(ON, V)).T
    else:
        RO = matmul(U, V).T
    return matmul(RO, r1.T).T, RO


# from the body frame we have : Ix w_x^2 + Iy w_y^2 + Iz w_z^2 = 2*E
# defines a ellipsoid:
#   1) w_z = sqrt(2E/Iz) cos(theta)
#   2) w_x = sqrt(2E/Ix) sin(theta) * sin(phi)
#   3) w_y = sqrt(2E/Iy) sin(theta) * cos(phi)
#   we know E and w_z = Lz/Iz, so theta = cos-1(w_z sqrt(Iz/2E)), uniform sample of phi and use 2) and 3) to calculate w_x and w_y.
# returns the classical angular momentum vector.
# note that even in the ekart frame, we are still using the energy obtained from the equlibrium geometry.
def ClassicalAngularMomentum(Ib, Lz, E, **dic):
    """
    Calculate the classical angular momentum vector for a rotating body.

    Parameters:
    Ib (numpy.ndarray): A numpy array containing the moments of inertia along the principal axes.
    Lz (float): The z-component of angular momentum.
    E (float): The total energy of the system.
    dic (dict): Additional keyword arguments, including 'phi' for specifying the azimuthal angle.

    Returns:
    numpy.ndarray: The classical angular momentum vector.
    """
    if "phi" in dic.keys():
        phi = dic["phi"]
    else:
        phi = tpi * random.random()
    Ix, Iy, Iz = Ib
    if abs(Iz) < 1.0e-8:
        wz = 0.0
        theta = hpi
    else:
        wz = Lz / Iz
        if E > 0.0:
            dm = wz * sqrt(Iz / (2 * E))
            if dm < -1.0 or dm > 1.0:
                print('WARNING, inverse cos not numerically right: is this |dm| a lot greater than 1? :',dm)
                dm = 1.0*(dm/abs(dm))
            theta = acos(dm)
        else:
            theta = 0.0
    wx = sqrt(2 * E / Ix) * sin(theta) * cos(phi)
    wy = sqrt(2 * E / Iy) * sin(theta) * sin(phi)
    return array([wx, wy, wz]) * Ib


def polar2xyz(R, psi, theta):
    """
    Convert spherical coordinates to Cartesian coordinates.

    Parameters:
    R (float): The radial distance.
    psi (float): The azimuthal angle.
    theta (float): The polar angle.

    Returns:
    numpy.ndarray: The Cartesian coordinates (x, y, z) corresponding to the given spherical coordinates.
    """
    zi = R * cos(theta)
    xi = R * sin(theta) * cos(psi)
    yi = R * sin(theta) * sin(psi)
    return array([xi, yi, zi])

def xyz2polar(xyz): 
    R = norm(xyz) 
    if R > 0.0:
      theta = np.arccos(xyz[2]/R)
      phi = arctan2(xyz[0],-xyz[1])
      return np.array([R, theta, phi])
    else:
      return np.array([0.0, 0.0, 0.0])

def GetRotTransVec(xx0, ms, el):
    """
    Calculate the rotation and translation vectors for a set of coordinates.

    Parameters:
    xx0 (numpy.ndarray): The input coordinates.
    ms (numpy.ndarray): The masses of the atoms.
    el (list): The list of element symbols.

    Returns:
    numpy.ndarray: The rotation and translation vectors.
    """
    x0 = xx0 - COM(xx0, ms).T
    # need to work in eckart frame
    RR, Ibm, Is = iI(iX(x0), iX(x0), ms)
    x0 = matmul(RR.T, x0.T).T
    na = x0.shape[0]
    ve = []
    for d in [x, y, z]:
        rr = cross(d, x0)
        ve.append(rr)
    for d in [x, y, z]:
        tt = zeros((na, 3)) + d
        ve.append(tt)
    vee = []
    for v in ve:
        # move back to origical frame
        vi = matmul(RR, v.T).T
        # mass scale and normalize
        vi = reshape(mscale(vi, ms, 1), (na * 3,))
        # orthonormalize the translation and rotations separately 
        #if norm(vi) > 0.0001:
        #    vi = vi / norm(vi)
        vee.append(vi.tolist())
    ve = vee
    debug = False
    # debug = True
    if debug:
        out = []
        for v in ve:
            v1 = mscale(reshape(array(v), (na, 3)), ms, -1)
            v1 = 4 * v1 / norm(v1)
            out += XYZlist(el, np.concatenate((xx0.T, v1.T)).T)
        open("6mod.xyz", "w").writelines(out)
    ve = array(ve)
    # debug=True
    if debug:
        print("OVER = ")
        print(matmul(ve, ve.T))
    # smooth out blocks
    Qr, _ = qr(ve[0:3,:].T)   
    Qt, _ = qr(ve[3:6,:].T)   
    ve = np.vstack([Qr.T, Qt.T]) 
    return ve

def ScaleTransform2(VV, **dic):
    """
    Scale transformation matrices for mass-scaled normal modes.

    Parameters:
    VV  (numpy.ndarray): Eigenvectors of mass-scaled hessian 
    freq (list): List of frequencies.
    mass (list): List of atomic masses.

    Returns:
    numpy.ndarray: Scaled transformation matrix for Cartesian to normal mode coordinates.
    """
    nmodes = VV.shape[1]
    if 'omega' in dic.keys():
      ws = abs(array(dic['omega']))
      for i in range(len(ws)):
          if abs(ws[i]) < 1.0e-8:
              ws[i] = 1.0
      if len(ws) != nmodes:
          raise ValueError(f"ScaleTransform2: len(omega)={len(ws)} does not match nmodes={nmodes}")
      W = np.diag(np.sqrt(ws))
      iW = np.diag(1/np.sqrt(ws))
    else:
      W = np.eye(nmodes)
      iW = W
    if 'mass' in dic.keys(): 
      M = np.sqrt(np.diag(dic['mass']))
      iM = np.sqrt(np.diag(1/dic['mass'])) 
    else:
      M = np.eye(VV.shape[0])
      iM = M
    x2Q = matmul(W,matmul(VV.T,M)) 
    Q2x = matmul(iM,matmul(VV,iW))
    p2P = matmul(iW,matmul(VV.T,iM))
    P2p = matmul(M,matmul(VV,W))
    return x2Q.T, Q2x.T, p2P.T, P2p.T


# to go from UN-mass-scaled to frequency mass scaled normal modes...
def ScaleTransform(x2n, n2x, freq, mass, dirr, fscale):
    """
    Scale transformation matrices for mass-scaled normal modes.

    Parameters:
    x2n (numpy.ndarray): Transformation matrix from Cartesian to normal mode coordinates.
    n2x (numpy.ndarray): Transformation matrix from normal mode to Cartesian coordinates.
    freq (list): List of frequencies.
    mass (list): List of atomic masses.
    dirr (bool): Direction of scaling (True for scaling up, False for scaling down).
    fscale (float): Scaling factor.

    Returns:
    numpy.ndarray: Scaled transformation matrix for Cartesian to normal mode coordinates.
    numpy.ndarray: Scaled transformation matrix for normal mode to Cartesian coordinates.
    """
    natoms = len(mass)
    nm     = n2x.shape[0]
    x2n2 = x2n.copy()
    n2x2 = n2x.copy()
    ws = abs(array(freq))
    for i in range(len(ws)):
        if abs(ws[i]) < 1.0e-8:
            ws[i] = 1.0
    if dirr:
        for i in range(nm):
            for j in range(natoms * 3):
                x2n2[j, i] = x2n2[j, i] * (sqrt(ws[i] ** fscale * mass[int(j / 3)]))
                n2x2[i, j] = n2x2[i, j] / (sqrt(ws[i] ** fscale * mass[int(j / 3)]))
    else:
        for i in range(nm):
            for j in range(natoms * 3):
                x2n2[j, i] = x2n2[j, i] / (sqrt(ws[i] ** fscale * mass[int(j / 3)]))
                n2x2[i, j] = n2x2[i, j] * (sqrt(ws[i] ** fscale * mass[int(j / 3)]))
    return x2n2, n2x2


def el2Mass(el):
    """
    Map element symbols to atomic masses.

    Parameters:
    el (list): List of element symbols.

    Returns:
    umpy.ndarray: List of atomic masses corresponding to the element symbols.
    """
    return array([el2m[e] for e in el])


def COM(r, mass):
    """
    Calculate the center of mass of a system.

    Parameters:
    r (numpy.ndarray): The array of coordinates.
    mass (numpy.ndarray): The masses of the atoms.

    Returns:
    numpy.ndarray: The center of mass coordinates.
    """
    com = np.sum(r.T * mass, axis=1, keepdims=True) / sum(mass)
    if False: # debug 
     m = np.diag(mass)
     on = np.ones(m.shape[1])
     mm = np.dot(mass,on)
     cc = np.matmul(np.matmul(r.T,m),on)/mm
     print('COM close? = ', np.allclose(cc,com))
    return com


def ReadXYZ(filn):
    """
    Read atomic coordinates from an XYZ file.

    Parameters:
    filn (str): The path to the XYZ file.

    Returns:
    list: List of element symbols.
    numpy.ndarray: Array of atomic coordinates.
    """
    dat = open(filn, "r").readlines()
    na = int(dat[0].strip())
    xcoo = []
    el = []
    for i in range(2, na + 2):
        lns = dat[i].strip().split()
        el.append(lns[0])
        xcoo.append([float(f) for f in lns[1:]])
    return el, array(xcoo)


def ReadXYZs(filn):
    """
    Read XYZ coordinate data from a file and return element symbols and atomic positions.

    Parameters:
    filn (str): The path to the XYZ file.

    Returns:
    tuple: A tuple containing two elements:
        - el (list): A list of element symbols in the same order as the atomic positions.
        - xcoos (list): A list of NumPy arrays, each containing the atomic coordinates for a molecule.
        - mess (list): A list of messages in the description part of the data

    Example:
    el, xcoos, mess = read_xyzs('molecule.xyz')
    """
    dat = open(filn, 'r').readlines()
    na = int(dat[0].strip())
    xcoos, xx, ell = [], [], []
    mess = []
    for ln in [la.strip() for la in dat]:
        lns = ln.split()
        if len(lns) == 1:
            if len(xx) > 0:
                if len(xcoos) == 0:
                    el = ell.copy()
                xcoos.append(np.array(xx))
            xx = []
        elif len(lns) == 4:
            xx.append([float(ff) for ff in lns[1:]])
            ell.append(lns[0])
        else:
            mess.append(ln)
    if len(xx) > 0:
        xcoos.append(np.array(xx))

    return el, xcoos, mess

def XYZlist(el, xcoo, **dic):
    """
    Generate an XYZ format file from element symbols and atomic coordinates.

    Parameters:
    el (list): List of element symbols.
    xcoo (numpy.ndarray): Array of atomic coordinates.
    dic (dict): Additional keyword arguments, including 'mess' for a message in the file.

    Returns:
    list: List of strings containing the XYZ file content.
    """
    if "mess" in dic.keys():
        message = dic["mess"]
    else:
        message = " "
    na, ni = xcoo.shape
    xout = [str(na) + "\n"]
    xout.append(message + "\n")
    xout += [
        el[i]
        + " "
        + "".join(["{0:14.9f}".format(xcoo[i, j]) + " " for j in range(ni)])
        + "\n"
        for i in range(na)
    ]
    return xout


def mscale(mat, mass, dr):
    """
    Mass-scale the given matrix.

    Parameters:
    mat (numpy.ndarray): The input matrix.
    mass (numpy.ndarray): The masses of the atoms.
    dr (int): Scaling direction (1 for scaling up, -1 for scaling down).

    Returns:
    numpy.ndarray: The mass-scaled matrix.
    """
    na = mat.shape[0]
    mo = mat.copy()
    for i in range(na):
        if dr == 1:
            mo[i, :] = mo[i, :] * sqrt(mass[i])
        else:
            mo[i, :] = mo[i, :] / sqrt(mass[i])
    return mo


def mscale2(arr, masses, dr=1):
    """
    Vectorized mass scaling. arr: shape (3N,) or (k, 3N).
    dr = +1 multiply by sqrt(m), dr = -1 divide.
    """
    arr = np.asarray(arr)
    N = masses.shape[0]
    scale = np.repeat(np.sqrt(masses), 3)
    if dr == -1:
        scale = 1.0 / scale
    if arr.ndim == 1:
        return arr * scale
    elif arr.ndim == 2:
        return arr * scale  # broadcast over rows
    else:
        raise ValueError("Unsupported shape")

# mass scale the coordinates as 1d tensor.
def mscale3(vi, mass, dr):
    """
    Mass-scale a 1D tensor of coordinates.

    Parameters:
    vi (numpy.ndarray): The 1D tensor of coordinates.
    mass (numpy.ndarray): The masses of the corresponding atoms.
    dr (int): Scaling direction (1 for scaling up, -1 for scaling down).

    Returns:
    numpy.ndarray: The mass-scaled 1D tensor of coordinates.
    """
    mo = vi.copy()
    if len(list(vi.shape)) == 1:
        nd = vi.shape[0]
        for k, i in enumerate(range(0, nd, 3)):
            if dr == 1:
                mo[i : i + 3] = mo[i : i + 3] * sqrt(mass[k])
            else:
                mo[i : i + 3] = mo[i : i + 3] / sqrt(mass[k])
    else:
        nd = vi.shape[1]
        for j in range(vi.shape[0]):
            for k, i in enumerate(range(0, nd, 3)):
                if dr == 1:
                    mo[j, i : i + 3] = mo[j, i : i + 3] * sqrt(mass[k])
                else:
                    mo[j, i : i + 3] = mo[j, i : i + 3] / sqrt(mass[k])
    return mo


def ProjectRTSpace(vvi,xxi,x0,mass,sp): 
    #vvi can be position or momenta. xxi is instantaneous position, and x0 is equilibrium position
    #turn to eckart 
    na,_ = vvi.shape 
    cx, cv =  COM(xxi, mass).T, COM(vvi, mass).T
    xx = xxi - cx
    vv = vvi - cv 
    U = EckartFrameTrans(x0,xx,mass)
    xx = matmul(U,xx.T).T    
    vv = matmul(U,vv.T).T    
    #get rotations and translations 
    rt = GetRotTransVec(xx, mass, [])  
    # build projectors 
    Pr, Pt = matmul(rt[:3,:].T,rt[:3,:]), matmul(rt[3:,:].T,rt[3:,:])
    Pi = np.eye(na*3) - Pr - Pt 
    #make final projector 
    PP = 0
    if 't' in sp:
      PP += Pt
    if 'r' in sp: 
      PP += Pr
    if 'i' in sp: 
      PP += Pi
    # scale and project and unscale   
    mv = mscale2(reshape(vv,(na*3,)), mass, +1)
    mv = matmul(PP,mv)
    vvo = reshape(mscale2(mv, mass, -1),(na,3)) 
    # return to original frame
    vvo = matmul(U.T,vvo.T).T  
    if 'r' or 'i' in sp: 
     return vvo  
    else:
     return vvo+cv




# XX matrix, cross product like matrix that represents the coordinate in the matrix=>vector mapping
def iX(r):
    """
    Generate the XX matrix for a given set of coordinates.

    Parameters:
    r (numpy.ndarray): The input coordinates.

    Returns:
    numpy.ndarray: The XX matrix.
    """
    na = r.shape[0]
    xo = zeros((na, 3, 3))
    for i in range(na):
        xo[i, 0, 1], xo[i, 0, 2], xo[i, 1, 2] = -r[i, 2], +r[i, 1], -r[i, 0]
        xo[i, :, :] = xo[i, :, :] - xo[i, :, :].T
    return xo


# uses the quaterion approach to eckart frame as discussed in
# stepanov et al, jcp, 140, 2014
# note that this is, at least most of the time, the same as doing:
# mx, me = mscale(self.sxx,self.mass,1), mscale(self.xxe,self.mass,1)
# U4 = rot_match_vec(mx,me)[1]
def EckartFrameTrans(xxe, sxx, mass):
    """
    Transform coordinates into the Eckart frame.

    Parameters:
    xxe (numpy.ndarray): Mass-scaled Cartesian coordinates.
    sxx (numpy.ndarray): Vectorial coordinates.
    mass (numpy.ndarray): The masses of the atoms.

    Returns:
    numpy.ndarray: Transformed coordinates in the Eckart frame.
    """
    C = zeros((4, 4))
    xp = xxe + sxx
    xm = xxe - sxx
    m = mass
    C[0, 1] = sum(m * (xp[:, 1] * xm[:, 2] - xm[:, 1] * xp[:, 2]))
    C[0, 2] = sum(m * (xm[:, 0] * xp[:, 2] - xp[:, 0] * xm[:, 2]))
    C[0, 3] = sum(m * (xp[:, 0] * xm[:, 1] - xm[:, 0] * xp[:, 1]))
    C[1, 2] = sum(m * (xm[:, 0] * xm[:, 1] - xp[:, 0] * xp[:, 1]))
    C[1, 3] = sum(m * (xm[:, 0] * xm[:, 2] - xp[:, 0] * xp[:, 2]))
    C[2, 3] = sum(m * (xm[:, 1] * xm[:, 2] - xp[:, 1] * xp[:, 2]))
    C += C.T
    C[0, 0] = sum(m * (xm[:, 0] ** 2 + xm[:, 1] ** 2 + xm[:, 2] ** 2))
    C[1, 1] = sum(m * (xm[:, 0] ** 2 + xp[:, 1] ** 2 + xp[:, 2] ** 2))
    C[2, 2] = sum(m * (xp[:, 0] ** 2 + xm[:, 1] ** 2 + xp[:, 2] ** 2))
    C[3, 3] = sum(m * (xp[:, 0] ** 2 + xp[:, 1] ** 2 + xm[:, 2] ** 2))
    ee, VV = eigh(C)
    # Find Orientat closest to equlibrium
    q = VV[:, 0]
    U = zeros((3, 3))
    U[0, 0] = q[0] ** 2 + q[1] ** 2 - q[2] ** 2 - q[3] ** 2
    U[1, 1] = q[0] ** 2 - q[1] ** 2 + q[2] ** 2 - q[3] ** 2
    U[2, 2] = q[0] ** 2 - q[1] ** 2 - q[2] ** 2 + q[3] ** 2
    U[0, 1] = 2 * (q[1] * q[2] + q[0] * q[3])
    U[1, 0] = 2 * (q[1] * q[2] - q[0] * q[3])
    U[0, 2] = 2 * (q[1] * q[3] - q[0] * q[2])
    U[2, 0] = 2 * (q[1] * q[3] + q[0] * q[2])
    U[1, 2] = 2 * (q[2] * q[3] + q[0] * q[1])
    U[2, 1] = 2 * (q[2] * q[3] - q[0] * q[1])
    if False:
        sxx = matmul(U, sxx.T).T
        print("U  =, angles = ", iR2q(U))
        print(U)
    return U


# get inertia matrix body-fixed frame 1) transformation, 2) I in body fixed fram (I diagonal), 3) in space fixed frame
def iI(XX1, XX2, mass):
    """
    Calculate the inertia matrix and eigenvectors for a set of coordinates.

    Parameters:
    XX1 (numpy.ndarray): XX matrix for coordinates.
    XX2 (numpy.ndarray): XX matrix for transformed coordinates.
    mass (numpy.ndarray): The masses of the atoms.

    Returns:
    numpy.ndarray: Eigenvectors of the inertia matrix.
    numpy.ndarray: Eigenvalues of the inertia matrix.
    numpy.ndarray: The inertia matrix.
    """
    debug = False
    # debug = True
    na = len(mass)
    II = zeros((3, 3))
    for i in range(na):
        II += mass[i] * matmul(XX1[i, :, :].T, XX2[i, :, :])
    evl, evc = eigh(II)
    if det(evc) < 0.0:
        evc = matmul(evc, diag([1.0, 1.0, -1.0]))
    if debug:
        print("iEVL = ", evl)
        print("EVC = ")
        print(evc)
    return evc, diag(evl), II


# rotate about some cartesian axis
def Rabout(ang, d):
    """
    Generate a rotation matrix for a rotation about a Cartesian axis.

    Parameters:
    ang (float): The rotation angle (in radians).
    d (int): Cartesian axis of rotation (0 for x, 1 for y, 2 for z).

    Returns:
    numpy.ndarray: The rotation matrix.
    """
    if d == 1:
        ro = array([[cos(ang), sin(ang)], [-sin(ang), cos(ang)]])
    else:
        ro = array([[cos(ang), -sin(ang)], [sin(ang), cos(ang)]])
    oo = np.eye(3)
    ii = array([i for i in range(3) if i != d]).reshape(2, 1)
    oo[ii, ii.T] = ro
    return oo

def iang2R(q,typ):
    if typ == 'eul': 
      return iq2R(q)
    else:
      return ixyz2R(q)
def iR2ang(q,typ):
    if typ == 'eul': 
      return iR2q(q)
    else:
      return iR2xyz(q)

# rotation matrix in euler parametrization
def iq2R(q):
    """
    Convert a vector to a rotation matrix in Euler parametrization.

    Parameters:
    q (numpy.ndarray): alpha, beta, gamma angles. range: -pi < alpha < pi, 0 < beta < pi , -pi < gam < pi

    Returns:
    numpy.ndarray: The rotation matrix.
    """
    alp, bet, gam = q
    return matmul(matmul(Rabout(alp, 2), Rabout(bet, 1)), Rabout(gam, 2))

def iR2q(R):
    """
    Convert a rotation matrix to angles in Euler parametrization. range: -pi < alpha < pi, 0 < beta < pi , -pi < gam < pi

    Parameters:
    R (numpy.ndarray): The rotation matrix.

    Returns:
    numpy.ndarray: Angles represented as (alpha, beta, gamma) angles.
    """
    if R[2, 2] < +1:
        if R[2, 2] > -1:
            bet = arccos(R[2, 2])
            alp = arctan2(R[1, 2], R[0, 2])
            gam = arctan2(R[2, 1], -R[2, 0])
        else:  # r22=−1
            bet = pi
            alp = -arctan2(R[1, 0], R[1, 1])
            gam = 0.0
    else:  # r22=+1
        bet = 0.0
        alp = arctan2(R[1, 0], R[1, 1])
        gam = 0.0
    return array([alp, bet, gam])

# rotation matrix in cartesian parametrization RxRyRz
def ixyz2R(w):
    """
    Convert a set of Euler angles to a rotation matrix in Cartesian parametrization.

    Parameters:
    w (numpy.ndarray): Euler angles represented as (alpha, beta, gamma) angles.
        
    Returns:
    numpy.ndarray: The rotation matrix.
    """
    return matmul(matmul(Rabout(w[0], 0), Rabout(w[1], 1)), Rabout(w[2], 2))
          
def iR2xyz(R): # RxRyRz
    """
    Convert a rotation matrix to a xyz rotation parametrization. 
    Note this assumes the range: -pi < w0 < pi, -pi/2 < w1 < pi/2 , -pi < w2 < pi 
        
    Parameters:
    R (numpy.ndarray): The rotation matrix.
        
    Returns:
    numpy.ndarray: angles represented as (wx, wy, wz) angles.
    """
    w = zeros(3)
    if R[0,2] < 1.000: 
      if R[0,2] >= -1:
        w[1] = arcsin( R[0,2])
        w[0] = arctan2(-R[1,2],R[2,2])
        w[2] = arctan2(-R[0,1],R[0,0])
      else:
        w[1] = -pi*0.5 
        w[2] = -arctan2(-R[1,0],R[1,1]) 
        w[0] = 0.0 
    else: 
      w[1] = pi*0.5 
      w[2] = arctan2(R[1,0],R[1,1]) 
      w[0] = 0.0 
    return w


# rotation matrix in cartesian parametrization  RzRyRx
def ixyz2R2(w):
    """
    Convert a set of Euler angles to a rotation matrix in Cartesian parametrization.

    Parameters:
    w (numpy.ndarray): Euler angles represented as (alpha, beta, gamma) angles.

    Returns:
    numpy.ndarray: The rotation matrix.
    """
    return matmul(matmul(Rabout(w[2], 2), Rabout(w[1], 1)), Rabout(w[0], 0))

def iR2xyz2(R): # RzRyRx
    """
    Convert a rotation matrix to a xyz rotation parametrization. 
    Note this assumes the range: -pi < w0 < pi, -pi/2 < w1 < pi/2 , -pi < w2 < pi 

    Parameters:
    R (numpy.ndarray): The rotation matrix.

    Returns:
    numpy.ndarray: angles represented as (wx, wy, wz) angles.
    """
    w = np.zeros(3)
    if R[2,0] < 1.0: 
      if R[2,0] > -1:
        w[1] = np.arcsin(-R[2,0])
        w[2] = np.arctan2(R[1,0],R[0,0])
        w[0] = np.arctan2(R[2,1],R[2,2])
      else:
        w[1] = np.pi*0.5 
        w[2] = -np.arctan2(-R[1,2],R[1,1]) 
        w[0] = 0.0 
    else: 
      w[1] = -np.pi*0.5 
      w[2] = -np.arctan2(-R[1,2],R[1,1]) 
      w[0] = 0.0 
    return w

def load_function_from_file(filepath, func_name):  
    import importlib.util  
    # Generate a module name (e.g., based on the filename)  
    module_name = os.path.splitext(os.path.basename(filepath))[0]  
    # Load the module                               
    spec = importlib.util.spec_from_file_location(module_name, filepath)  
    module = importlib.util.module_from_spec(spec)  
    sys.modules[module_name] = module  
    spec.loader.exec_module(module)  
    # Get the function by name                      
    if not hasattr(module, func_name):  
        raise AttributeError(f"Function '{func_name}' not found in '{filepath}'")  
      
    return getattr(module, func_name)  



def stdq(q):
    """
    Standardize a quaternion to the range [0, 2π].

    Parameters:
    q (numpy.ndarray): Quaternion represented as (alpha, beta, gamma) angles.

    Returns:
    numpy.ndarray: Standardized quaternion.
    """
    for j in [0, 2]:
        # temp shift to 0:2pi
        q[j] += pi
        if q[j] < 0:
            q[j] = tpi + (q[j] - ceil(q[j] / tpi) * tpi)
        elif q[j] >= 2 * pi:
            q[j] = q[j] - floor(q[j] / tpi) * tpi
        q[j] -= pi
    if q[1] < 0:
        q[1] = pi + (q[1] - ceil(q[1] / pi) * pi)
    elif q[1] > pi:
        q[1] = q[1] - floor(q[1] / pi) * pi
    return q


def File2InputList(fnam):
    """
    Read input parameters from a text file.

    Parameters:
    fnam (str): Path to the input file.

    Returns:
    list: List of input parameter pairs.
    """
    lines = open(fnam, "r").readlines()
    inpd = []
    for line in lines:
        # Support inline comments while preserving valid key/value content.
        line = line.split("#", 1)[0].strip()
        if len(line) == 0:
            continue
        if "=" in line:
            key, rhs = line.split("=", 1)
            inpd.append([key.strip().lower(), rhs.strip().split()])
    return inpd
def _suggest_bins(x):
    """Robust binning heuristic."""
    n = len(x)
    if n < 2 or np.all(x == x[0]):
      return "auto", 1
    q25, q75 = np.percentile(x, [25, 75])
    iqr = q75 - q25
    sigma = np.std(x)

    if np.issubdtype(x.dtype, np.integer) and np.ptp(x) <= 50:
        return "integer", int(np.ptp(x)) + 1

    # Freedman-Diaconis (robust for most data, good default)
    h_fd = 2 * iqr / np.cbrt(n) if iqr > 0 else 0
    # Scott (optimal for Gaussian)
    h_scott = 3.5 * sigma / np.cbrt(n) if sigma > 0 else 0

    # Choose the smaller bin-width (higher resolution), but not too small
    if h_fd > 0 and h_scott > 0:
        h = min(h_fd, h_scott)
    elif h_fd > 0:
        h = h_fd
    elif h_scott > 0:
        h = h_scott
    else:
        # fallback to sqrt rule if everything else fails
        return "sqrt", None

    # Limit bins to reasonable numbers (max ~100 bins)
    n_bins = int(np.clip(np.ceil((x.max() - x.min()) / h), 10, 100))
    return "auto", n_bins

def _hist_token(name):
    txt = str(name).strip().lower()
    txt = re.sub(r"[^a-z0-9]+", "_", txt).strip("_")
    return txt or "x"

def _hist_stage_alias(stage):
    stg = _hist_token(stage)
    return {"initial": "ini", "sampled": "sam"}.get(stg, stg)

def _hist_scope_alias(scope):
    scp = _hist_token(scope)
    return {
        "system": "sys",
        "molecule_m0": "m0",
        "molecule_m1": "m1",
    }.get(scp, scp)

def hist_filename(stage, scope, metric):
    stg = _hist_stage_alias(stage)
    scp = _hist_scope_alias(scope)
    met = _hist_token(metric)
    return f"hist_{stg}_{scp}_{met}.py"

def hist_outdir(base_dir="histograms", stage="initial", scope="system"):
    return pathlib.Path(base_dir) / _hist_token(stage) / _hist_token(scope)

def set_hist_base_dir(base_dir):
    global _HIST_BASE_DIR
    _HIST_BASE_DIR = pathlib.Path(base_dir)

def hist_script(data, name, out_dir="hist_scripts"):
    x = np.asarray(data).ravel()
    rule, k = _suggest_bins(x)

    out_dir = pathlib.Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    #script_path = out_dir / f"make_{name}_hist.py"
    script_path = out_dir / f"hist_{name}.py"

    # Format data as multiline JSON array for readability
    data_json = json.dumps(x.tolist(), indent=4)

    default_bins_literal = repr(k if k is not None else rule)

    tmpl = rf'''#!/usr/bin/env python3
"""Histogram of '{name}' auto-generated by export_histogram_script.

Default bin choice: {rule.upper()} {f"(bins={k})" if k else ""}.
"""

import argparse, matplotlib.pyplot as plt, numpy as np, pathlib, json

# Embedded data (modify or replace with file-loading if needed)
data = np.array({data_json})

# CLI arguments
p = argparse.ArgumentParser()
p.add_argument("--bins", default=None,
               help="Bin rule name (scott, fd, sqrt, sturges, auto) "
                    "or explicit integer.")
p.add_argument("--edges", default=None,
               help="Comma-separated bin edges (overrides --bins).")
p.add_argument("--no-show", action="store_true")
p.add_argument("--outfile", default="{name}.png",
               help="Output filename base (extension ignored)")
p.add_argument("--pdf", action="store_true",
               help="Also write PDF in addition to PNG")
args = p.parse_args()

# Bin determination
if args.edges:
    bins = np.fromstring(args.edges, sep=",")
elif args.bins:
    bins = int(args.bins) if args.bins.isdigit() else args.bins
else:
    bins = {default_bins_literal}

# Plot histogram
fig, ax = plt.subplots(figsize=(7,5))
counts, edges, _ = ax.hist(data, bins=bins, color="tab:blue",
                           alpha=0.7, edgecolor="k")

# Mean and standard deviation annotation
if len(data) == 0:
    mu, sigma = float("nan"), float("nan")
elif len(data) == 1:
    mu, sigma = float(data[0]), 0.0
else:
    mu, sigma = np.mean(data), np.std(data, ddof=1)
ax.annotate(rf"$\mu={{mu:.3g}},\;\sigma={{sigma:.3g}}$",
            xy=(0.98, 0.95), xycoords="axes fraction",
            ha="right", va="top", fontsize=10,
            bbox=dict(boxstyle="round", fc="w", ec="0.7", alpha=0.8))

ax.set(title="Histogram of {name}",
       xlabel="{name}", ylabel="Counts")

fig.tight_layout()
out = pathlib.Path(args.outfile)
fig.savefig(out.with_suffix(".png"), dpi=150)
if args.pdf:
    fig.savefig(out.with_suffix(".pdf"))
    print("Wrote", out.with_suffix(".png"), "and", out.with_suffix(".pdf"))
else:
    print("Wrote", out.with_suffix(".png"))

if not args.no_show:
    plt.show()
'''

    script_path.write_text(textwrap.dedent(tmpl))
    script_path.chmod(0o755)
    return script_path

def hist_emit(data, metric, *, stage="initial", scope="system", base_dir=None):
    if base_dir is None:
        base_dir = _HIST_BASE_DIR
    out_dir = hist_outdir(base_dir=base_dir, stage=stage, scope=scope)
    name = f"{_hist_stage_alias(stage)}_{_hist_scope_alias(scope)}_{_hist_token(metric)}"
    return hist_script(data, name, out_dir=out_dir)




'''
Symmetric orthogonalization leakage correction
Paper:
Colclough, G. L., Brookes, M., Smith, S. M. and Woolrich, M. W., "A symmetric multivariate leakage correction for MEG connectomes," NeuroImage 117, pp. 439-448 (2015)
Translated from MATLAB:
https://github.com/OHBA-analysis/MEG-ROI-nets/blob/master/%2BROInets/symmetric_orthogonalise.m
Main function:
closest_ortho_matrix(dat)
dat: np.array with k x n shape
k: number of regions or sensors or sources of interest
n: number of samples
'''
def symmetric_ortho(dat):
    U,S,V = svd(dat,full_matrices=0)
    #rank checking
    S = S #this is different from matlab, as the diagonal component is automatically obtained in np svd, for matlab, we need diag(S)
    tol = max(dat.shape)*S[0]*(np.finfo(dat.dtype).eps) #tolerance level
    r = np.sum(S>tol) #number of S larger than tolerance
    isFullRank = (r >= dat.shape[0]) #dat.shape[0] here is number of ROIs
    #in matlab -> [U,S,V] = svd(a)
    #in python U, S, Vh = linalg.svd(a) and V = Vh.T
    '''if isFullRank == False:
        print('Warning: The input ts matrix is not full rank.')
        print(r)
        print(dat.shape[0])'''
    L = U.dot(np.conj(V))
    #W = V.T.dot(np.diag(1/S)).dot(V) #working weights, but not using
    return(L,isFullRank)

#fast_svd.py assumes dat is already transposed
def fast_svd(dat,N):
    #N = 1
    if N < dat.shape[1]:
        eigs2 = eig(dat.dot(dat.T))
        #eigs2 = eig(dat.dot(dat))
        d = max(eigs2[0])

        U = eigs2[1][:,0]

        S = np.sqrt(np.abs(d))
        V = dat.T.dot(U.dot(1/S))

        #U = dat.dot(V.dot(1/S))
    return(S) # for the purpose of tolerance finding, only need S, need a constant

def scale_cols(dat,s):
    newdat = dat * s
    return(newdat)

def reldiff(a,b):
    if a == 0 or b == 0:
        outcome = 0
    else:
        outcome = (2*np.abs(a-b) / (np.abs(a)+np.abs(b)))
    return(outcome)

def closest_ortho_matrix(dat):
    print('Starting symmetric orthogonalization leakage correction')
    #dat = dat.T # data has to be transposed before processing as per their matlab code... weird practice
    itere = 0
    #dat = dat.astype(np.float64) #use double precision
    MAX_ITER = 2e2
    #slightly different because of single precision float32
    tol = np.finfo(dat.dtype).eps
    print(tol)
    A_b = np.conj(dat)
    d = np.sqrt(np.sum(dat.conj()*A_b,axis=0))
    rho = []
    Ls = []
    '''
    dot(A,B) of same size is simply in matlab:
    sum(conj(A).*B)

    in python it is:
    np.sum(A.conj()*B, axis=0)
    '''
    isFullRank = True
    while itere < MAX_ITER:
        V, isFullRank = symmetric_ortho(scale_cols(dat,d))
        d = np.sum(A_b.conj()*V,axis=0)
        L = scale_cols(V,d)
        Ls.append(L)
        if isFullRank == False:
            print('  No longer full rank. Optimal matrix reached at iteration %s' % (str(itere)))
            break
        E = dat - L
        rho.append(np.sqrt(np.sum(np.sum(E.conj()*np.conj(E),axis=0))))
        if itere > 0:
            val = reldiff(rho[itere],rho[itere-1])
            print('  Iteration: %s\n   Tolerance: %s\n   Relative difference: %s\n   Rhos: %s' % (str(itere+1),str(tol),str(val),str(rho[itere])))
            if val <= tol:
                print('  Optimal matrix reached at iteration %s\n  Tolerance: %s\n  Relative difference: %s\n  Rhos: %s' % (str(itere+1),str(tol),str(val),str(rho[itere])))
                break
        itere+=1
    if isFullRank == False:
        return(Ls[-1])
    else:
        return(L)


def plot_wl_weights(weights,
                    bin_edges=None,
                    J_range=(0.0, 40),
                    #J_range=(0.0, 100),
                    is_log=False,
                    #is_log=True,
                    normalize=True,
                    title="",
                    ylabel="target weight (2J+1)/Omega(J)",
                    #title="Wang–Landau weights vs. J",
                    outfile="wl_weights.png"):

    """
    Plot Wang–Landau weights as bars across φ and save to file.

    Parameters
    ----------
    weights : array-like
        Vector of weights per φ-bin. If `is_log=True`, these are log-weights (ln w).
    bin_edges : array-like, optional
        Bin edges for φ (length = len(weights)+1). If None, assume uniform bins over J_range.
    J_range : (float, float)
        Range for φ when bin_edges is None. Default (0, 2π).
    is_log : bool
        True if `weights` are log-weights (ln w). False if already linear weights.
    normalize : bool
        If True, rescale so the tallest bar is 1.
    title : str
        Plot title.
    ylabel : str
        Label for the plotted WL quantity.
    outfile : str
        Filename (with path) where the figure is saved.
    """
    w = np.asarray(weights, dtype=float)
    n = w.size

    # Bin edges
    edges = (np.linspace(J_range[0], J_range[1], n + 1)
             if bin_edges is None else np.asarray(bin_edges, dtype=float))
    if edges.size != n + 1:
        raise ValueError("bin_edges must have length len(weights)+1")

    # Convert to linear weights if needed
    y = np.exp(w - (w.max() if normalize else 0.0)) if is_log else w.copy()
    if normalize and not is_log and y.max() > 0:
        y /= y.max()

    centers = 0.5 * (edges[:-1] + edges[1:])
    widths  = np.diff(edges)

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.bar(centers, y, width=widths, align="center",
    color="tab:blue",   # fill color
    alpha=0.6,          # transparency
    edgecolor="k"       # bar border color
           )
    ax.set_xlim(edges[0], edges[-1])

    # Pretty φ ticks
    if np.isclose(edges[0], 0.0) and np.isclose(edges[-1], 2*np.pi):
        xt = [0, 0.5*np.pi, np.pi, 1.5*np.pi, 2*np.pi]
        ax.set_xticks(xt)
        ax.set_xticklabels([r"$0$", r"$\tfrac{\pi}{2}$", r"$\pi$",
                            r"$\tfrac{3\pi}{2}$", r"$2\pi$"])

    font = {'size'   : 12}
    plt.rc('font', **font)
    ax.tick_params(axis='both', which='major', labelsize=11)
    ax.set_xlabel("Total angular momentum J",fontsize=12)
    ax.set_ylabel(ylabel,fontsize=11)
    #ax.set_ylabel("weight" + (" (normalized)" if normalize else ""))
    ax.set_title("")
    #ax.set_title(title)
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout(pad=1.2)
    #plt.show()
    fig.savefig(outfile, dpi=300)
    plt.close(fig)
    return outfile


def write_wl_plot_script(weights,
                         bin_edges=None,
                         *,
                         J_range=(0.0, 100),
                         is_log=False,
                         normalize=True,
                         title="Wang–Landau weights vs. J",
                         ylabel="target weight (2J+1)/Omega(J)",
                         outfile="wl_weights.png",
                         script_path="make_wl_plot.py"):
  
    """
    Emit a helper script that reproduces the WL-weight histogram by
    importing *this* module and re-calling `plot_wl_weights`.
    --------------------------------------------------------------------
    Parameters are identical to `plot_wl_weights`, plus
    --------------------------------------------------------------------
    script_path : str or Path
        Where to write the burner script (default: make_wl_plot.py)
    """
    # ––––– assemble the call-string exactly as Python would see it –––––
    modname = inspect.getmodule(write_wl_plot_script).__name__

    arg_lines = [
        f"weights=np.asarray({json.dumps(np.asarray(weights).tolist())}, dtype=float)",
        ("bin_edges=np.asarray(" + json.dumps(np.asarray(bin_edges).tolist()) + ", dtype=float)"
         if bin_edges is not None else "bin_edges=None"),
        f"J_range=({J_range[0]}, {J_range[1]})",
        f"is_log={is_log}",
        f"normalize={normalize}",
        f"title={title!r}",
        f"ylabel={ylabel!r}",
        f"outfile={outfile!r}"
    ]
    call_string = ", ".join(arg_lines)

   
    tpl = f"""\
    #!/usr/bin/env python
    \"\"\"Helper script – auto-generated by `{modname}.write_wl_plot_script`.

    Edit the arrays at the top or tweak any keyword arguments, then run

        python {pathlib.Path(script_path).name}

    to regenerate the figure (`{outfile}`).
    \"\"\"

    import numpy as np
    import {modname} as wlmod

    # ----- raw data -----  (feel free to edit) -------------------------
    {arg_lines[0]}
    {arg_lines[1]}
    # ------------------------------------------------------------------

    wlmod.plot_wl_weights({call_string})
    print("Plot saved to", {outfile!r})
    """

    pathlib.Path(script_path).write_text(textwrap.dedent(tpl))
    return script_path
