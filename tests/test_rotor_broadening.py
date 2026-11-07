from types import SimpleNamespace

import numpy as np
import pytest
from scipy.integrate import cumulative_trapezoid

from icats.dist import RotorAzimuthVonMises, uniform
from icats.mc import InitICDF
from icats.molecules import imolecule
from icats.rotor_sampling import RotorProjectionICDF, projection_cdf


@pytest.mark.parametrize('j,k,width', [(1,1,.8), (5,2,.8), (5,5,.8),
                                      (40,1,.8), (40,40,.8), (5,3,.2), (5,3,2.)])
def test_bounded_cdf_and_independent_moment(j,k,width):
    p = np.zeros(2*j+1)
    p[j-k] = p[j+k] = .5
    grid,cdf,lam = projection_cdf(j,tuple(p),width)
    assert cdf[0] == 0 and cdf[-1] == 1
    assert np.all(np.diff(cdf) >= 0)
    # Independent finer quadrature of the actual analytical density.
    r = np.sqrt(j*(j+1))
    x = np.linspace(-r,r,131073)
    logd = np.logaddexp(-.5*((x-k)/width)**2,-.5*((x+k)/width)**2)+lam*(x/r)**2
    d = np.exp(logd-logd.max())
    assert abs(np.trapz(d*x*x,x)/np.trapz(d,x)-k*k)/(j*(j+1)) < 1e-7
    finer = cumulative_trapezoid(d,x,initial=0)
    finer /= finer[-1]
    assert np.max(np.abs(np.interp(grid,x,finer)-cdf)) < 1e-5


def test_zero_moment_and_discrete_limit():
    rng = [np.random.default_rng(24)]
    p = [0,0,1,0,0]
    assert RotorProjectionICDF(rng,.8,2,p)[0] == 0
    assert RotorProjectionICDF(rng,.8,0,[1])[0] == 0
    p = [.2,0,0,0,.8]
    values = [RotorProjectionICDF(rng,0,2,p)[0] for _ in range(2000)]
    assert set(values) == {-2.,2.}
    assert abs(np.mean(np.array(values)==2)-.8)<.04


def selected_state(rng, record, j=None):
    return np.asarray(record,dtype=float if np.isscalar(record) else object), 1.


def make_asymmetric(tmp_path, constants, rep, j, state):
    mol = imolecule()
    a,b,c = constants
    mol.sp = SimpleNamespace(rA=a,rB=b,rC=c,rsym='asym-prolate' if rep==1 else 'asym-oblate')
    mol.ip.diagdir = str(tmp_path)
    mol.ip.name = 'test'
    mol.SetOrientatConvention(np.eye(3),.5/np.array(constants),rep)
    mol.AsymRigidRotorProjEnergies(j)
    record = [col[state] for col in mol.sp.RR[j]]
    sample = SimpleNamespace(sJ=j,SampInfo={},dist={
        'rotJ':{},
        'arotJz':{'cont':InitICDF(1,selected_state,[record],seed=23)},
        'rotK':{'cont':InitICDF(1,RotorProjectionICDF,[.8],par0=[0,[1]],seed=41)},
        'arotGam':{'cont':InitICDF(1,RotorAzimuthVonMises,[],par0=[1.,1.],seed=52)}})
    return mol,sample


@pytest.mark.parametrize('constants,rep,state', [([3.,1.4,1.],1,3),
    ([3.,2.9,1.],3,5),([1.,.9999,.9998],1,2)])
def test_production_moments_and_energy(tmp_path,constants,rep,state):
    constants = np.array(constants)*1e-4
    mol,sa = make_asymmetric(tmp_path,constants,rep,5,state)
    sp = mol.sp
    target = np.array([col[state] for col in sp.RR[5][3:6]])
    body = np.array([sp.J2c+2*sp.Jcc,sp.J2c-2*sp.Jcc,sp.J2c+sp.Jzc])
    energy = sp.J2c*30+sp.RR[5][1][state]
    np.testing.assert_allclose(body@target,energy,rtol=1e-9)
    p = sp.rotor_projection_probabilities[5][:,state]
    np.testing.assert_allclose(p.sum(),1.,atol=1e-14)
    np.testing.assert_allclose(p@np.arange(-5,6)**2,target[2],atol=1e-8)
    vectors,energies = [],[]
    for _ in range(12000):
        mol.SampleRigidRotorState(sa)
        vectors.append(sa.SampInfo['rot']['vecJ'].copy())
        energies.append(sa.SampInfo['rot']['cenergy'])
    squared = np.array(vectors)**2
    np.testing.assert_allclose(squared.sum(axis=1),30,atol=1e-12)
    assert np.all(abs(squared.mean(axis=0)-target) < 6*squared.std(axis=0)/np.sqrt(len(squared))+1e-6)
    np.testing.assert_allclose(energies,squared@body,rtol=1e-13)
    assert abs(np.mean(energies)-energy) < 6*np.std(energies)/np.sqrt(len(energies))+1e-10


def test_linear_rotor_needs_no_broadening():
    mol = imolecule()
    mol.sp = SimpleNamespace(rsym='linear',J2c=1e-4,Jzc=0.,Jcc=0.,Irep=0)
    sa = SimpleNamespace(sJ=5,SampInfo={},dist={'rotJ':{},
        'srotGam':{'cont':InitICDF(1,uniform,[-np.pi,np.pi],seed=52)}})
    mol.SampleRigidRotorState(sa)
    assert abs(sa.srpar[-1][2]) < 1e-14
    np.testing.assert_allclose(sa.SampInfo['rot']['cenergy'],30e-4)


@pytest.mark.parametrize('k', [0,2,5])
def test_symmetric_top_uses_broadening_and_preserves_energy(k):
    mol = imolecule()
    mol.sp = SimpleNamespace(rsym='oblate',J2c=2e-4,Jzc=-1e-4,Jcc=0.,Irep=3)
    sa = SimpleNamespace(sJ=5,SampInfo={},dist={'rotJ':{},
        'srotJz':{'cont':InitICDF(1,selected_state,[k],seed=12)},
        'srotGam':{'cont':InitICDF(1,uniform,[-np.pi,np.pi],seed=52)},
        'rotK':{'cont':InitICDF(1,RotorProjectionICDF,[.8],par0=[0,[1]],seed=41)}})
    errors, projections = [],[]
    for _ in range(6000):
        mol.SampleRigidRotorState(sa)
        errors.append(sa.SampInfo['rot']['energy_error'])
        projections.append(sa.srpar[-1][2])
    assert abs(np.mean(errors)) < 6*np.std(errors)/np.sqrt(len(errors))+1e-12
    if k:
        assert np.std(projections) > .1
    else:
        assert np.max(np.abs(projections)) < 1e-14


def test_revised_tutorials_do_not_use_wl(tmp_path):
    from icats.cli import _generate_tutorial
    for name in ('revision_nh3_h2o','thermal_atom_diatom_ar_no'):
        path = tmp_path/name
        _generate_tutorial(name,path)
        text = (path/'tutorial_input.txt').read_text()
        assert 'wang = False' in text
        assert 'orbital-sampling = geometric' in text
        assert not (path/'rd_tutorial_input'/'wang.pkl').exists()
