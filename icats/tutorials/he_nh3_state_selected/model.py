"""Rigid He-NH3 analytic MOLSCAT example, atomic units (manual Eq. 5.7).

This omits ammonia inversion dynamics. It is a tutorial potential, not an
ab initio replacement for the full molecular system.
"""
import numpy as np

CM=219474.6313705
MU=3.24*1822.888486209
CONSTANTS=np.array([9.9441,9.9441,6.2294])/CM
INERTIA=1/(2*CONSTANTS)


def cross_matrix(w):
    x,y,z=w
    return np.array([[0.,-z,y],[z,0.,-x],[-y,x,0.]])


def potential(r):
    radius=np.sqrt(r@r)
    x,y,z=n=r/radius
    a,b,c,d,e=(np.sqrt(1/(4*np.pi)),np.sqrt(3/(4*np.pi)),np.sqrt(5/(16*np.pi)),
               np.sqrt(7/(16*np.pi)),-np.sqrt(35/(16*np.pi)))
    harmonics=np.array([a,b*z,c*(3*z*z-1),d*(5*z**3-3*z),e*(x**3-3*x*y*y)])
    dn=np.array([[0,0,0],[0,0,b],[0,0,6*c*z],[0,0,d*(15*z*z-3)],
                 [3*e*(x*x-y*y),-6*e*x*y,0]])
    decay=np.exp(-2.4*radius)
    radial=np.array([-1.2e7/radius**6+6.2e8*decay,-5.4e7*decay,2.2e6/radius**6,
                     3.6e7*decay,1.e8*decay])
    derivative=np.array([7.2e7/radius**7-2.4*6.2e8*decay,2.4*5.4e7*decay,
                         -6*2.2e6/radius**7,-2.4*3.6e7*decay,-2.4*1.e8*decay])
    angular_gradient=radial@dn
    gradient=((derivative@harmonics)*n+(angular_gradient-n*(n@angular_gradient))/radius)
    return radial@harmonics/CM,gradient/CM


def rhs(t,state):
    r,p,ell=state[:3],state[3:6],state[6:9]
    q=state[9:].reshape(3,3)
    rb=q.T@r
    _,gradient=potential(rb)
    omega=ell/INERTIA
    torque=np.cross(rb,gradient)
    return np.concatenate((p/MU,-q@gradient,torque-np.cross(omega,ell),
                           (q@cross_matrix(omega)).ravel()))


def invariants(state):
    r,p,ell=state[:3],state[3:6],state[6:9]
    q=state[9:].reshape(3,3)
    energy=p@p/(2*MU)+np.sum(CONSTANTS*ell**2)+potential(q.T@r)[0]
    return energy,np.cross(r,p)+q@ell,np.max(np.abs(q.T@q-np.eye(3)))
