from types import SimpleNamespace

import numpy as np
import pytest

from icats.constants import kboltz
from icats.functions import iI, iX
from icats.rotor_states import state_tables, residue_weights
from icats.ensemble_statistics import rotor_energy_mean
from test_rotor_classification import molecule_from_constants, molecule_from_geometry
from test_state_selection import molecule


def water(tmp_path, weights=(1., 3.), axis='b'):
    m = molecule_from_geometry([[0.,0.,0.], [.758,0.,.587], [-.758,0.,.587]], [16.,1.,1.])
    m.sp.el = ['O','H','H']
    m.sp.nm = 0
    m.ip.name = 'water'
    m.ip.diagdir = str(tmp_path)
    m.ip.filename = 'water'
    m.ip.MaxR = 4
    m.ip.MaxV = -1
    m.ip.Trot = 1e6  # Synthetic geometry uses small masses; populate all test levels.
    m.ip.rotor_axis = axis
    m.ip.rotor_k_weights = weights
    m.StandardOrientat()
    m.AsymRigidRotorProjEnergies(4)
    return m


@pytest.mark.parametrize('axis,rep', [('a',1),('b',2),('c',3)])
def test_axis_override_spectrum_and_cartesian_frame(tmp_path, axis, rep):
    m = water(tmp_path, axis=axis)
    assert m.sp.Irep == rep
    tensor = iI(iX(m.sp.x0),iX(m.sp.x0),m.sp.mass)[2]
    np.testing.assert_allclose(tensor, np.diag(m.sp.Ib), atol=1e-12)
    ref = water(tmp_path, axis='a')
    for j in range(5):
        np.testing.assert_allclose(np.array(m.sp.RR[j][1])+m.sp.J2c*j*(j+1),
                                   np.array(ref.sp.RR[j][1])+ref.sp.J2c*j*(j+1), atol=1e-10)


@pytest.mark.parametrize('weights', [(1.,0.),(0.,1.),(1.,3.)])
def test_water_joint_weights_and_energy(tmp_path, weights):
    m = water(tmp_path, weights)
    pj, ps = state_tables(m, m.ip.Trot)
    expected, actual, energies = [], [], []
    for j in range(5):
        p = m.sp.rotor_projection_probabilities[j]
        even = (np.arange(-j,j+1)%2 == 0).astype(float) @ p
        assert np.all(np.minimum(abs(even),abs(even-1)) < 1e-10)
        w = np.where(even > .5,weights[0],weights[1])
        e = m.sp.J2c*j*(j+1)+np.array(m.sp.RR[j][1])
        expected.extend((2*j+1)*w*np.exp(-e/(kboltz*m.ip.Trot)))
        actual.extend(pj[j]*ps[j])
        energies.extend(e)
    expected = np.array(expected); expected /= expected.sum()
    np.testing.assert_allclose(actual,expected,atol=1e-14)
    np.testing.assert_allclose(rotor_energy_mean(m,pj,m.ip.Trot),np.dot(expected,energies))
    # For water J=1: 1_01 and 1_10 are ortho, 1_11 is para.
    p1 = m.sp.rotor_projection_probabilities[1]
    np.testing.assert_allclose(np.array([0,1,0])@p1,[0,1,0],atol=1e-10)


def test_incompatible_rule_is_not_projection(tmp_path):
    m = water(tmp_path, (1.,0.,0.))
    with pytest.raises(ValueError,match='mixes differently weighted'):
        state_tables(m,m.ip.Trot)


@pytest.mark.parametrize('values', [(),(0,0),(-1,1),(np.nan,1),(np.inf,1),(1,2,3)])
def test_invalid_weights(values):
    with pytest.raises(ValueError):
        residue_weights(SimpleNamespace(rotor_k_weights=values))


def test_ammonia_preset_equivalence():
    for rule, weights in [('nh3-para',(0.,1.,1.)), ('nh3-ortho',(1.,0.,0.))]:
        old = molecule(rule=rule)
        new = molecule(); new.ip.rotor_k_weights = weights
        p,q = state_tables(old,100)
        pp,qq = state_tables(new,100)
        np.testing.assert_allclose(p,pp)
        for a,b in zip(q,qq):np.testing.assert_allclose(a,b)


def test_water_zero_temperature_and_sampling(tmp_path):
    m = water(tmp_path,(0.,1.))
    pj,ps = state_tables(m,0)
    assert pj[1] == 1
    assert ps[1][0] == 1
    m.ip.Trot = 0
    sa = SimpleNamespace(seed=42,SampInfo={})
    m.InitialDist(sa)
    for _ in range(100):
        m.SampleTotMolAngMom(sa)
        m.SampleRigidRotorState(sa)
        assert sa.sJ == 1
        np.testing.assert_allclose(np.dot(sa.srpar[-1],sa.srpar[-1]),2,atol=1e-12)


def test_ambiguous_or_non_symmetry_axis_rejected():
    for constants,axis in [([1.,1.,1.],'b'),([1.,1.,.6],'a')]:
        m = molecule_from_constants(constants)
        m.ip.rotor_axis = axis
        with pytest.raises(ValueError):m.StandardOrientat()


def test_input_parser_and_conflicting_rules(monkeypatch):
    m = molecule()
    m.sp.na = 4
    monkeypatch.setattr(m, 'GenNormalModesMat', lambda: None)
    m.ip.inpd = [('rotor-axis',['b']), ('rotor-k-weights',['1','3'])]
    m.GenerateInputData()
    assert m.ip.rotor_axis == 'b'
    assert m.ip.rotor_k_weights == (1.,3.)
    m.ip.rotor_k_rule = 'nh3-para'
    with pytest.raises(ValueError,match='Do not combine'):residue_weights(m.ip)


@pytest.mark.parametrize('width',[0.,.8])
def test_weighted_water_sampler_and_repeatability(tmp_path, width):
    m = water(tmp_path)
    m.ip.rotor_k_width = width
    pj, ps = state_tables(m,m.ip.Trot)
    def draw():
        sa = SimpleNamespace(seed=123,SampInfo={})
        m.InitialDist(sa)
        draws = []
        for _ in range(2000):
            m.SampleTotMolAngMom(sa)
            m.SampleRigidRotorState(sa)
            draws.append((sa.sJ,sa.srpar[1],*sa.srpar[-1]))
        return np.array(draws)
    data = draw()
    np.testing.assert_array_equal(data,draw())
    np.testing.assert_allclose(np.sum(data[:,2:]**2,axis=1),data[:,0]*(data[:,0]+1),atol=1e-11)
    observed = np.bincount(data[:,0].astype(int),minlength=len(pj))/len(data)
    assert np.max(abs(observed-pj)) < .04
    expected_energy = rotor_energy_mean(m,pj,m.ip.Trot)
    assert abs(data[:,1].mean()-expected_energy) < 5*data[:,1].std()/np.sqrt(len(data))
