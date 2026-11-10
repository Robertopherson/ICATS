import math
import unittest
from types import SimpleNamespace

from icats.orbital_channel import channel_impact


class OrbitalChannelTests(unittest.TestCase):
    def test_discrete_magnitudes_and_area(self):
        for ell in (0, 1, 40, 183):
            for k in (.5, 13.615392):
                b, area = channel_impact(ell,k,'sqrt')
                self.assertAlmostEqual((b*k)**2,ell*(ell+1))
                self.assertAlmostEqual(area,math.pi*((ell+1)**2-ell**2)/k**2)
                b2, area2 = channel_impact(ell,k,'langer')
                self.assertAlmostEqual(b2*k,ell+.5)
                self.assertEqual(area,area2)

    def test_continuous_shell_measure(self):
        for ell in (0,10):
            samples=[channel_impact(ell,2.,'uniform-b',(i+.5)/1000) for i in range(1000)]
            self.assertTrue(all(ell/2 <= b <= (ell+1)/2 for b,_ in samples))
            self.assertAlmostEqual(sum(w for _,w in samples)/1000,math.pi*(2*ell+1)/4)

    def test_invalid(self):
        for ell,k,rule in ((-.1,1,'sqrt'),(.5,1,'sqrt'),(0,0,'sqrt'),
                           (0,float('nan'),'sqrt'),(0,1,'unknown')):
            with self.assertRaises(ValueError):
                channel_impact(ell,k,rule)

    def test_cartesian_orbital_setup(self):
        import numpy as np
        from icats.iscattering import icats
        from icats.mc import InitICDF
        from icats.dist import uniform
        system=icats()
        system.sp.rmass=2.
        system.ip.ImpactPhi=0.
        system.ip.OrbitalChannel=0
        for rule in ('sqrt','langer','uniform-b'):
            system.ip.OrbitalChannelRule=rule
            sample=SimpleNamespace(sV=3.,SampInfo={},slog=[],
                dist={'orbital_channel':InitICDF(1,uniform,[0,1],seed=22)})
            system.SampleFixedImpactOrbital(sample)
            self.assertAlmostEqual(np.linalg.norm(sample.scL),6*sample.sb)
            self.assertEqual(sample.SampInfo['orb']['channel_L'],0)
            self.assertGreaterEqual(sample.SampInfo['orb']['area_weight_bohr2'],0)
        system.ip.OrbitalChannel=None
        system.ip.FixedB=.7
        sample=SimpleNamespace(sV=3.,SampInfo={},slog=[])
        system.SampleFixedImpactOrbital(sample)
        self.assertAlmostEqual(sample.snL,4.2)
        self.assertNotIn('channel_L',sample.SampInfo['orb'])


if __name__=='__main__':
    unittest.main()
