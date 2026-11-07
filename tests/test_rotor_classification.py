from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from icats.dist import IsotropicDistICDF
from icats.functions import iI, iX, Rabout
from icats.mc import InitICDF
from icats.molecules import imolecule


def molecule_from_geometry(x, mass=None):
    mol = imolecule()
    x = np.asarray(x, dtype=float)
    mol.sp = SimpleNamespace(x0=x.copy(), na=len(x), nd=3*len(x),
                             mass=np.ones(len(x)) if mass is None else np.asarray(mass),
                             el=['H']*len(x))
    return mol


def molecule_from_constants(constants):
    inertia = .5/np.asarray(constants)
    second = inertia.sum()/2-inertia
    assert np.all(second >= 0)
    positive = np.diag(np.sqrt(second/2))
    x = np.vstack((positive,-positive))
    # Test alignment from a non-principal initial frame.
    rotation = Rabout(.37, 0)@Rabout(.61, 2)
    return molecule_from_geometry(x@rotation.T+[1.,2.,3.])


@pytest.mark.parametrize('constants,kind', [
    ([1.,.8,.6], 'asym-prolate'),
    ([1.,.8001,.6], 'asym-oblate'),
    ([1.,1.,1.], 'spherical'),
    ([1.,1.-4e-6,1.-8e-6], 'spherical'),
    ([1.,1.-6e-6,1.-12e-6], 'asym-prolate'),
    ([1.,.6,.6], 'prolate'),
    ([1.,1.,.6], 'oblate'),
    ([1.,.61,.6], 'asym-prolate'),
    ([1.,.99,.6], 'asym-oblate'),
])
@pytest.mark.parametrize('scale', [1., 1e-4])
def test_classification_and_frame(constants,kind,scale):
    mol = molecule_from_constants(np.array(constants)*scale)
    with np.errstate(divide='raise',invalid='raise'):
        mol.StandardOrientat()
    assert mol.sp.rsym == kind
    np.testing.assert_allclose(np.linalg.det(mol.sp.S2B),1.,atol=1e-12)
    tensor = iI(iX(mol.sp.x0),iX(mol.sp.x0),mol.sp.mass)[2]
    np.testing.assert_allclose(tensor,np.diag(mol.sp.Ib),rtol=1e-11,atol=1e-9)
    if kind=='spherical':
        assert np.isnan(mol.sp.asymk)
        assert mol.sp.Jzc == mol.sp.Jcc == 0
        np.testing.assert_allclose(mol.sp.J2c,np.mean(constants)*scale)
    else:
        a,b,c = constants
        np.testing.assert_allclose(mol.sp.asymk,(2*b-a-c)/(a-c),atol=1e-9)
    np.testing.assert_allclose(mol.sp.rotor_spread,
                               (max(constants)-min(constants))/np.mean(constants),atol=1e-12)


@pytest.mark.parametrize('rep,order', [(1,[1,2,0]),(2,[2,0,1]),(3,[0,1,2])])
@pytest.mark.parametrize('constants', [[1.,.8,.6],[1.,.6,.6],[1.,1.,1.]])
def test_convention_defines_rotation_even_for_equal_constants(rep,order,constants):
    mol = imolecule()
    mol.sp = SimpleNamespace(rA=constants[0],rB=constants[1],rC=constants[2])
    inertia = .5/np.array(constants)
    mol.SetOrientatConvention(np.eye(3),inertia.copy(),rep)
    np.testing.assert_allclose(mol.sp.S2B@np.diag(inertia)@mol.sp.S2B.T,
                               np.diag(inertia[order]),atol=1e-12)
    np.testing.assert_allclose(mol.sp.Ib,inertia[order])


def test_linear_and_atom_unchanged():
    atom = molecule_from_geometry([[1.,2.,3.]])
    atom.StandardOrientat()
    assert atom.sp.rsym=='atom'
    linear = molecule_from_geometry([[0.,0.,-1.],[0.,0.,1.]])
    linear.StandardOrientat()
    assert linear.sp.rsym=='linear'
    assert np.isnan(linear.sp.asymk)
    assert np.isnan(linear.sp.rotor_spread)


def test_methane_reference_is_spherical():
    path = Path(__import__('icats').__file__).parent/'examples/meth_geom.xyz'
    x = np.loadtxt(path,skiprows=2,usecols=(1,2,3))
    mol = molecule_from_geometry(x,[12.,1.007825,1.007825,1.007825,1.007825])
    mol.StandardOrientat()
    assert mol.sp.rsym=='spherical'


def test_classified_spherical_top_samples_without_projection_sampler():
    mol = molecule_from_constants([1.,1.,1.])
    mol.StandardOrientat()
    cont = InitICDF(1,IsotropicDistICDF,
                    ['eul',np.array([2*np.pi,1,2*np.pi]),np.array([np.pi,0,np.pi])],seed=16)
    sample = SimpleNamespace(sJ=5,SampInfo={},dist={'rotJ':{},'arotJz':{'cont':cont}})
    mol.SampleRigidRotorState(sample)
    np.testing.assert_allclose(np.dot(sample.srpar[-1],sample.srpar[-1]),30)
    np.testing.assert_allclose(sample.srpar[1],30*mol.sp.J2c)


def test_kappa_zero_reaches_asymmetric_eigensolver(tmp_path):
    constants = np.array([1.,.8,.6])*1e-4
    mol = molecule_from_constants(constants)
    mol.ip.name = 'kappa_zero'
    mol.ip.diagdir = str(tmp_path)
    mol.ip.MaxR = 3
    mol.StandardOrientat()
    mol.AsymRigidRotorProjEnergies(3)
    assert mol.sp.rsym=='asym-prolate'
    for j in range(4):
        k = np.arange(-j,j+1)
        jp = np.diag(np.sqrt(j*(j+1)-k[:-1]*(k[:-1]+1)),-1)
        jx,jy,jz = (jp+jp.T)/2,(jp-jp.T)/(2j),np.diag(k)
        expected = np.linalg.eigvalsh(sum(c*(op@op) for c,op in zip(constants,(jx,jy,jz))))
        recovered = np.array(mol.sp.RR[j][1])+mol.sp.J2c*j*(j+1)
        np.testing.assert_allclose(recovered,expected,atol=1e-8)
    probabilities = mol.ThermalRotorProbabilities(500)
    assert np.isfinite(probabilities).all()
    np.testing.assert_allclose(probabilities.sum(),1.)


@pytest.mark.parametrize('kappa', [-.8,-.4,0.,.4,.8])
@pytest.mark.parametrize('epsilon', [1.,.1,.01])
def test_near_spherical_shape_keeps_closest_pair_axis(kappa,epsilon,tmp_path):
    d = np.array([1.,kappa,-1.])-kappa/3
    constants = (1.+.2*epsilon*d)*1e-4
    mol = molecule_from_constants(constants)
    mol.ip.name = 'axis_test'
    mol.ip.diagdir = str(tmp_path)
    mol.ip.MaxR = 5
    mol.StandardOrientat()
    assert mol.sp.rsym == ('asym-prolate' if kappa <= 0 else 'asym-oblate')
    assert mol.sp.Irep == (1 if kappa <= 0 else 3)
    np.testing.assert_allclose(.5/mol.sp.Ib[2],constants[0 if kappa <= 0 else 2])
    mol.AsymRigidRotorProjEnergies(5)
    for j in range(6):
        k = np.arange(-j,j+1)
        jp = np.diag(np.sqrt(j*(j+1)-k[:-1]*(k[:-1]+1)),-1)
        ops = [(jp+jp.T)/2,(jp-jp.T)/(2j),np.diag(k)]
        expected = np.linalg.eigvalsh(sum(c*(op@op) for c,op in zip(constants,ops)))
        recovered = np.array(mol.sp.RR[j][1])+mol.sp.J2c*j*(j+1)
        np.testing.assert_allclose(recovered,expected,atol=1e-8)
        assert len(mol.sp.RR[j][-1]) == 2*j+1
        np.testing.assert_allclose(np.square(mol.sp.RR[j][-1]),
                                   np.maximum(0,np.array(mol.sp.RR[j][5])-np.square(mol.sp.RR[j][-2])),
                                   atol=2e-10)


@pytest.mark.parametrize('constants,kind', [
    ([1.,.6001,.6], 'asym-prolate'),
    ([1.,.9999,.6], 'asym-oblate'),
    ([1.,.600001,.6], 'prolate'),
    ([1.,.999999,.6], 'oblate'),
])
def test_symmetric_model_requires_equal_pair_within_tolerance(constants,kind):
    mol = molecule_from_constants(constants)
    mol.StandardOrientat()
    assert mol.sp.rsym == kind
    assert (mol.sp.Jcc == 0) == ('asym' not in kind)
