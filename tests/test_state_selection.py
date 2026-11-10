from types import SimpleNamespace
import numpy as np
import pytest

from icats.constants import kboltz
from icats.rotor_states import state_tables, projection_icdf
from icats.molecules import imolecule
from icats.mc import InitICDF, ICDFscalar
from icats.rotor_sampling import RotorProjectionICDF
from icats.dist import RotorTotBoltzICDF, uniform
from icats.iscattering import icats


def molecule(jk=None, rule='all', allowed=None):
    m=imolecule()
    m.sp=SimpleNamespace(rsym='oblate',el=['N','H','H','H'],J2c=9.9441/219474.6313705,
                         Jzc=(6.2294-9.9441)/219474.6313705,Irep=3)
    m.ip.MaxR=15;m.ip.Trot=100.
    m.ip.rotor_state=jk;m.ip.rotor_k_rule=rule;m.ip.rotor_allowed_k=allowed
    return m


@pytest.mark.parametrize('jk',[(3,2),(5,4),(0,0)])
def test_fixed_tables(jk):
    m=molecule(jk)
    for t in (0,100,500):
        p,q=state_tables(m,t)
        assert p[jk[0]] == 1
        assert np.count_nonzero(p)==1
        assert np.isclose(q[jk[0]][jk[0]+jk[1]],1 if jk[1]==0 else .5)


@pytest.mark.parametrize('rule',['nh3-para','nh3-ortho'])
def test_thermal_joint_normalization(rule):
    m=molecule(rule=rule)
    p,q=state_tables(m,100)
    expected=[];actual=[]
    for j in range(len(p)):
        k=np.arange(-j,j+1)
        w=(k%3!=0 if rule=='nh3-para' else k%3==0).astype(float)
        expected.extend((2*j+1)*w*np.exp(-(m.sp.J2c*j*(j+1)+m.sp.Jzc*k*k)/(kboltz*100)))
        actual.extend(p[j]*q[j])
    expected=np.asarray(expected);expected/=expected.sum()
    np.testing.assert_allclose(actual,expected,atol=1e-15)


def test_ortho_inversion_unresolved_level_count():
    # JPL/LAMDA: j=3,K=3 has two inversion levels of weight 28 each;
    # j=3,K=0 has one of weight 28. Neglecting splitting gives ratio 2.
    m=molecule(rule='nh3-ortho')
    _,q=state_tables(m,100)
    ratio=(q[3][0]+q[3][6])/q[3][3]
    expected=2*np.exp(-9*m.sp.Jzc/(kboltz*100))
    assert np.isclose(ratio,expected)


def test_forbidden_and_unsupported_states():
    for m in (molecule((3,3),'nh3-para'),molecule((3,2),'nh3-ortho'),molecule((3,2),allowed=(1,))):
        with pytest.raises(ValueError,match='No rotational states'):state_tables(m,100)
    m=molecule((3,2));m.sp.rsym='asym-oblate'
    with pytest.raises(ValueError,match='symmetric top'):state_tables(m,100)
    m=molecule((3,2),'nh3-para');m.sp.el=['O','H','H']
    with pytest.raises(ValueError,match='require NH3'):state_tables(m,100)


def test_zero_temperature_restricted_ground_state():
    p,q=state_tables(molecule(rule='nh3-para'),0)
    assert p[1]==1
    np.testing.assert_array_equal(q[1],[.5,0,.5])


@pytest.mark.parametrize('width',[0.,.8])
def test_sampled_fixed_state_and_moments(width):
    m=molecule((3,2),'nh3-para');m.ip.rotor_k_width=width
    p,q=state_tables(m,0)
    sa=SimpleNamespace(SampInfo={},dist={
        'rotJ':{'cont':InitICDF(1,RotorTotBoltzICDF,[p],seed=33)},
        'srotJz':{'cont':InitICDF(1,projection_icdf,[q],par0=[3],seed=34)},
        'srotGam':{'cont':InitICDF(1,uniform,[-np.pi,np.pi],seed=35)},
        'rotK':{'cont':InitICDF(1,RotorProjectionICDF,[width],par0=[3,q[3]],seed=36)}})
    samples=[]
    for _ in range(5000):
        m.SampleTotMolAngMom(sa);m.SampleRigidRotorState(sa)
        assert sa.sJ==3 and abs(sa.SampInfo['rot']['jz'])==2
        samples.append(sa.srpar[-1].copy())
    v=np.array(samples)
    np.testing.assert_allclose(np.sum(v*v,axis=1),12,atol=1e-12)
    np.testing.assert_allclose(np.mean(v*v,axis=0),[4,4,4],atol=.22)
    if width==0:np.testing.assert_allclose(v[:,2]**2,4,atol=1e-12)


def test_flat_b_area_and_restart():
    sc=icats();sc.ip.orbital_sampling='flat-b';sc.ip.MinB=2.;sc.ip.MaxB=3.
    sc.ip.ImpactPhi=0.;sc.sp.rmass=2.
    def draw():
        sa=SimpleNamespace(sV=3.,SampInfo={},slog=[],dist={'flat_b':InitICDF(1,uniform,[2,3],seed=17)})
        result=[]
        for _ in range(2000):
            sc.SampleFixedImpactOrbital(sa)
            assert 2 <= sa.sb <=3
            assert np.isclose(np.linalg.norm(sa.scL),6*sa.sb)
            result.append(sa.SampInfo['orb']['area_weight_bohr2'])
        return result
    a=draw();assert a==draw()
    assert abs(np.mean(a)-5*np.pi)<.15
