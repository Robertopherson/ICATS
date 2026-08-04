#!/usr/bin/env python3
from .constants import *
from .functions import *
from .dist import *
from .angmom import Pz2, Py2, Px2, Pz, Px, Py, WangTran, Wz2, Wy2, Wx2
from spherical import Wigner
import quaternionic as quat
from .mc import MCsample, ICDFsample, ICDFscalar, InitMC, InitICDF

INFO_LABEL_WIDTH = 22

def info_vec(label, vec, unit="", mag_label=None, mag=None, qm=None, fmt="{:12.5f}"):
    vals = ", ".join(fmt.format(float(v)) for v in vec)
    line = "{0:<{w}} = [ {1} ]".format(label, vals, w=INFO_LABEL_WIDTH)
    if unit:
        line += "  " + unit
    if mag_label is not None:
        if mag is None:
            mag = norm(vec)
        mag_name = "|{0}|".format(mag_label)
        line += ", {0:<7s} = {1:12.5f}".format(mag_name, float(mag))
    if qm is not None:
        line += "  (QM: {0:6.2f})".format(float(qm))
    return line + "\n"

def info_scalar(label, value, unit="", fmt="{:14.7f}"):
    unit_txt = ("  " + unit) if unit else ""
    return "{0:<{w}} = {1}{2}\n".format(label, fmt.format(float(value)), unit_txt, w=INFO_LABEL_WIDTH)

def info_angle_vec(label, vec, comment):
    vals = ", ".join("{:10.7f}".format(float(v) / pi) for v in vec)
    return "{0:<{w}} = [ {1} ]  pi rad   # {2}\n".format(label, vals, comment, w=INFO_LABEL_WIDTH)


def OrientationPDFSurface(angs, func, pars):
    """Sampling surface for user orientation PDFs.

    The user function is the physical density P(alpha,beta,gamma). ICATS
    applies the Euler measure sin(beta) here.
    """
    alpha, beta, gamma = [float(a) for a in angs]
    p = float(func(alpha, beta, gamma, *pars))
    if not np.isfinite(p):
        raise ValueError("orientation PDF returned a non-finite value")
    if p < 0.0:
        raise ValueError("orientation PDF returned a negative value")
    return p * sin(beta)


 

class imolecule:
    """
    Represents a molecule and provides methods to read input data and generate molecular information.

    Attributes:
        na (int): Number of atoms.
        nd (int): Number of atomic coordinates (3 times the number of atoms).
        el (list): List of element symbols.
        mass (numpy.ndarray): Array of atomic masses.
        log (list): List of log messages.
        MaxV (int): Maximum velocity.
        Tvib (float): Vibrational temperature.
        Trot (float): Rotational temperature.
        Tvel (float): Velocity temperature.
        ordist (int): Order distribution.

    Methods:
        ReadInput(self, fnam): Read input data from a file.
        GenerateInputData(self): Process and generate input data.
        NModesOut(self, filn): Generate normal modes output.
        GenNormalModesMat(self): Generate normal modes matrix.
    """
    class sample: 
      def __init__(self,ii,self1,**dic):
        sp = self1.sp 
        if 'seed' in dic.keys():
          self.seed = dic['seed']
        else:
          self.seed = int(ii*np.random.random()*709)
        self.svv = zeros(sp.xxe.shape)
        self.sxx = sp.xxe.copy()
        self.SampInfo = {}
        return

    def initsyspar(self):
        """Initialize a molecule parameters.

        This method initializes a new scattering sample, setting up necessary variables and data structures.
        """
        class msyspar: 
            def __init__(self):
             self.na = 0
             self.nd = 0
             self.el = []
             self.mass = array([])
             return
        self.sp = msyspar() 

    def initinppar(self):
        """Initialize the molecule simulation parameters.

        Args:
            dic (dict): parameters for initialization.
        """
        class minppar: 
            def __init__(self):
              self.mi = -1
              self.MaxV = 10
              self.Tvib = 0
              self.Trot = 0
              self.Tvel = 0
              self.VelPar = []
              self.ordist = "isotropic"
              self.orientation_frame = "scattering"
              self.orfilename = ""
              self.orfunction = ""
              self.orpars = []
              self.orthin = 25
              self.nfreeze = []
              self.rotpar = 'xyz'
              self.isotropic = True
              self.velfwhm = -1
              self.lowdin = False
              self.diagdir = "."
              return 
        self.ip = minppar()

    def __init__(self, **dic):
        """Initialize the imolecule object with default values."""
        self.initinppar()
        self.initsyspar()
        self.log = []
  
    def ReadInput(self, fnam):
        """Read input data from a file and generate molecular information.
        Args:
            fnam (str): Name of the input file."""
        self.ip.inpd = File2InputList(fnam)
        self.ip.filename = fnam
        self.ip.prefix = fnam.split(".")[0]
        self.ip.name = self.ip.prefix
        self.log.append("Reading input file : " + fnam + "\n")
        self.GenerateInputData()

    def GenerateInputData(self):
        """Process and generate molecular data from the input parameters."""
        sp = self.sp 
        ip = self.ip 
        for ky, val in ip.inpd:
            if ky == "xyz":
                sp.el, sp.x0 = ReadXYZ(val[0])
                self.log += ["   Input Orientation (Ang): \n"]
                self.log += ["   " + ln for ln in XYZlist(sp.el, sp.x0)]
                sp.x0 = sp.x0 * ang2au
                sp.na = len(sp.el)
                sp.nd = sp.na * 3
                sp.mass = el2Mass(sp.el) * amu2au
                sp.mass2 = np.repeat(sp.mass, 3)
            if ky == "trot":
                ip.Trot = float(val[0])
                if sp.na == 1:
                    self.log += ["   Ignoring rotational temperature for atom; using 0.0\n"]
                else:
                    self.log += ["   Temperature for Rotational States: " + val[0] + " \n"]
            if ky == "tvib":
                ip.Tvib = float(val[0])
                if sp.na == 1:
                    self.log += ["   Ignoring vibrational temperature for atom; using 0.0\n"]
                else:
                    self.log += ["   Temperature for Vibrational States: " + val[0] + " \n"]
            if ky == "vel":
                self.log +=['Molecular Velocity (m/s) centre: ' + val[0]+ '\n']
                self.log +=[' Full Width Half-Maximum (FWHM): ' + val[1] + '\n']
                self.log +=['  order of Velocify  v^n       : ' + val[2] + '\n']
                ip.VelPar = [float(val[0])*mps2au,float(val[1])*mps2au,int(val[2]) ]  # v0,fwhm, ^n 
            if ky == "hess":
                sp.HH = np.loadtxt(val[0])
            if ky == "name":
                self.name = val[0]
            if ky == "nfreeze":
                ip.nfreeze = [int(n) for n in val]
                sp.nfreeze = ip.nfreeze
            if ky == "ordist":
                if 'read' == val[0]:
                  self.log += ["Orientation distribution PDF file: " + val[1] + " function: " + val[2] + "\n"]
                  ip.ordist = "pdf"
                  ip.orfilename = val[1]
                  ip.orfunction = val[2]
                  ip.orpars = [float(v) for v in val[3:]]
                  ip.rotpar = "euler"
                elif 'fixed' == val[0]:
                  ip.ordist = "fixed"
                  ip.orpars = [float(v) for v in val[1:4]]
                  ip.rotpar = "euler"
            if ky == "orientation-mode":
                mode = val[0].lower()
                if mode in ("isotropic", "default"):
                    ip.ordist = "isotropic"
                    self.log += ["Molecular orientation mode: isotropic\n"]
                elif mode == "fixed":
                    if len(val) < 4:
                        raise ValueError("orientation-mode = fixed requires alpha beta gamma")
                    ip.ordist = "fixed"
                    ip.orpars = [float(v) for v in val[1:4]]
                    ip.rotpar = "euler"
                    self.log += ["Molecular orientation mode: fixed Euler angles\n"]
                elif mode == "pdf":
                    if len(val) < 3:
                        raise ValueError("orientation-mode = pdf requires file and function")
                    ip.ordist = "pdf"
                    ip.orfilename = val[1]
                    ip.orfunction = val[2]
                    ip.orpars = [float(v) for v in val[3:]]
                    ip.rotpar = "euler"
                    self.log += ["Molecular orientation mode: user PDF " + ip.orfilename + ":" + ip.orfunction + "\n"]
                else:
                    raise ValueError("Unknown orientation-mode: " + val[0])
            if ky == "orientation-frame":
                ip.orientation_frame = val[0].lower()
                if ip.orientation_frame != "scattering":
                    raise ValueError("Only orientation-frame = scattering is currently supported")
                self.log += ["Molecular orientation PDF frame: ICATS scattering frame\n"]
            if ky == "orientation-thin":
                ip.orthin = max(1, int(val[0]))
                self.log += ["Molecular orientation MC thinning: " + val[0] + "\n"]
            if ky == "rot-param":
                self.log += ["Rotation angle Parametrization : " + val[0] + "\n"]
                ip.rotpar = "euler" if val[0] == "eul" else val[0]
        if sp.na == 1:
            ip.Trot = 0.0
            ip.Tvib = 0.0
            self.log += ["   Atom detected: forced molecular Trot/Tvib to 0.0\n"]
        self.GenNormalModesMat()
        #self.StandardOrientat()
        #if self.nm > 0:
        #    self.NModesOut(self.name + "_nm.xyz")

    def NModesOut(self, filn):
        """
        Generate normal modes output and save it to a file.

        Args:
            filn (str): Output file name.
        """
        out = []
        sp = self.sp 
        for m in range(sp.ntr, sp.nd):
            vv = zeros(sp.nd)
            vv[m] = 4.0
            mo = reshape(matmul(sp.n2x.T, vv), (sp.na, 3))
            out += [ str(sp.na) + "\n",
                    "Mode " + str(m) + " Freq :" +
                    "{0:14.9f}".format(sp.w[m]) + " \n",]
            for i in range(sp.na):
                out += [sp.el[i] + " "+ "".join([
                            "{0:14.9f}".format(f) + " "
                            for f in (sp.x0[i, :] * au2ang).tolist()
                            + mo[i, :].tolist() ]) + "\n" ]
        open(filn, "w").writelines(out)

    def GenNormalModesMat(self):
        """
        Generate normal modes matrix and vibrational frequencies.

        If there is only one atom, it computes translational and rotational modes.
        If there are more atoms, it calculates the vibrational modes.

        Returns:
        list: Log messages.
        """
        sp = self.sp 
        ip = self.ip 
        if sp.na == 1:
            # Single atoms have only translational DOF in this model:
            # no internal vibrations and no rotational rigid-rotor space.
            sp.c2n = np.eye(sp.nd)
            sp.ntr = sp.nd
            sp.nm = 0
            sp.w = zeros(sp.ntr)
            self.log += ["   No harmonic frequencies :\n"]
        elif not hasattr(sp,'HH'):
            sp.c2n = GetRotTransVec(sp.x0, sp.mass, sp.el).T
            sp.ntr = sp.c2n.shape[1]
            sp.nm  = 0 
            sp.w = zeros(sp.ntr)
            self.log += ["   No harmonic frequencies :\n"]
        else:
            evl, sp.c2n = eigh(sp.HH)
            if evl[0] < -0.0001:
                sp.log.append("   WARNING: some negative EigenValues detected in Hessian Matrix...\n")
                print("WARNING: some negative EigenValues detected in Hessian Matrix...")
                print("MOLECULE ", sp.el)
                print(evl)
            sp.w = array([sqrt(max([1.0e-15, e])) for e in evl])
            self.log += ["   Frequencies (cm-1):\n"]
            self.log += ["   " + str(i).rjust(2) + " " + "{0:14.9f}".format(e * au2cm) + "\n"
                          for i, e in enumerate(sp.w) ]

    def MolecularVeloc(self,sa):
        """
        Calculate the center of mass velocity.

        Returns:
        numpy.ndarray: Center of mass velocity.
        """
        return COM(sa.svv, self.sp.mass)

    def MolecularPosition(self,sa):
        """
        Calculate the center of mass position.

        Returns:
        numpy.ndarray: Center of mass position.
        """
        return COM(sa.sxx, self.sp.mass)

    def CalcOrient(self,sa):
        """
        Calculate molecular orientation in terms of Euler angles.

        Returns:
        list: Log messages.
        """
        sp = self.sp 
        ip = self.ip 
        log = [f"{ip.name:<{INFO_LABEL_WIDTH}} = reconstructed orientation\n"]
        if sp.na > 1:
            R = EckartFrameTrans(sp.xxe, sa.sxx, sp.mass).T 
            alpha, beta, gamma = iR2q(R)
            wx, wy, wz = iR2xyz(R)
            if sp.na == 2:  # these are arbitrary for 2 atoms
              gamma, wz = 0.0, 0.0
            if 'ori' not in sa.SampInfo.keys():
             sa.SampInfo['ori'] = {}
            sa.SampInfo['ori']["seul"] = [alpha, beta, gamma]
            sa.SampInfo['ori']["sxyz"] = [wx, wy, wz]
            sa.SampInfo['ori']["salpha"] = alpha
            sa.SampInfo['ori']["sbeta"] = beta
            sa.SampInfo['ori']["sgamma"] = gamma
            # Backward-compatible names used by older histogram helpers.
            sa.SampInfo['ori']["sphi"] = alpha
            sa.SampInfo['ori']["stheta"] = beta
            sa.SampInfo['ori']["schi"] = gamma
            log += [info_angle_vec("alpha,beta,gamma", [alpha, beta, gamma], "molecular Euler")]
            log += [info_angle_vec("wx,wy,wz", [wx, wy, wz], "XYZ rotation angles")]
        else:
            log += [f"{'orientation space':<{INFO_LABEL_WIDTH}} = none\n"]
        return log

    def CalcInterEner(self,sa):
        """
        Calculate the internal energy of the system.

        Returns:
        list: Log messages.
        """
        sp = self.sp 
        ip = self.ip 
        log = [f"{ip.name:<{INFO_LABEL_WIDTH}} = vibrational analysis\n"]
        if sp.na == 1:
          log += [f"{'vibrational space':<{INFO_LABEL_WIDTH}} = none (atom)\n"]
          return log
        if sp.nm == 0:
          xx = sa.sxx.copy() - COM(sa.sxx, sp.mass).T
          vv = sa.svv.copy() - COM(sa.svv, sp.mass).T
          U = EckartFrameTrans(sp.xxe, xx, sp.mass)
          xx = matmul(U, xx.T).T
          vv = matmul(U, vv.T).T
          RR = sp.c2m[:, : sp.ntr]
          PP = np.eye(sp.nd)-matmul(RR, matmul(inv(matmul(RR.T, RR)), RR.T))
          dx = matmul(PP, reshape(xx - sp.xxe, (sp.nd)))
          dv = matmul(PP, reshape(vv, (sp.nd)))
          dp = sp.mass2 * dv
          ek = 0.5 * sum(sp.mass2 * dv**2)
          log += [f"{'vibrational modes':<{INFO_LABEL_WIDTH}} = unavailable; no Hessian/normal modes\n"]
          log += [f"{'internal residual':<{INFO_LABEL_WIDTH}} = projected outside translation/rotation\n"]
          log += [info_scalar("residual |dx|", norm(dx) * au2ang, "Ang", "{:14.7f}")]
          log += [info_scalar("residual |p|", norm(dp), "au", "{:14.7f}")]
          log += [info_scalar("residual kinetic", ek * au2ev, "eV", "{:14.7f}")]
          return log
        debug = True
        debug = False
        xx, vv = (
            sa.sxx.copy() - COM(sa.sxx, sp.mass).T,
            sa.svv.copy() - COM(sa.svv, sp.mass).T,)
        if debug:
            RR = sp.m2c[: sp.ntr, :]
            RR = mscale2(GetRotTransVec(
                sp.xxe, sp.mass, sp.el), sp.mass, -1)
            print("ZEROy? = ", norm(RR), diag(
                matmul(RR, sp.c2m[:, : sp.ntr])))
        # move molecule to its eckart frame
        U = EckartFrameTrans(sp.xxe, xx, sp.mass)
        # move to eckart frame:
        xx = matmul(U, xx.T).T
        vv = matmul(U, vv.T).T
        xe = sp.xxe
        # project velocity to space of internals coordinates:
        RR = sp.c2m[:, : sp.ntr]
        PP = np.eye(sp.nd)-matmul(RR, matmul(inv(matmul(RR.T, RR)), RR.T))
        xx = matmul(PP, reshape(xx - xe, (sp.nd)))
        pp = matmul(PP, reshape(vv, (sp.nd))) * sp.mass2
        # Convert to mass, frequency scaled Normal coordinates:
        Np = matmul(sp.p2n.T, pp)
        Nq = matmul(sp.x2n.T, xx)
        EE = 0.5 * (Nq**2 + Np**2)
        for i in range(sp.ntr, sp.nd):
            EE[i] = EE[i] * sp.w[i]
        if debug:
            for i in range(sp.nd):
                print(i,"{0:14.9f}".format(Nq[i]) + " "
                      + "{0:14.9f}".format(Np[i]) + " E = "
                      + "{0:14.9f}".format(EE[i]))
            quit()
        if 'vib' not in sa.SampInfo.keys(): 
          sa.SampInfo['vib'] = {}
        sa.SampInfo['vib']["sqp"] = [Nq,Np]
        sa.SampInfo['vib']["sQ"] = Nq
        sa.SampInfo['vib']["sP"] = Np
        sa.SampInfo['vib']["senergy"] = EE[sp.ntr:]
        if debug: 
          print('##################################')
          print('QQ = ', Nq)
          print('Q2 = ', sa.SampInfo['vib']["Q"])
          print('QE = ', Nq[6:]-sa.SampInfo['vib']["Q"][6:])
          print('PP = ', Np)
          print('P2 = ', sa.SampInfo['vib']["P"])
          print('PE = ', Np[6:]-sa.SampInfo['vib']["P"][6:])
          print('##################################')
        log += [" mode    freq     ~vstat       Q         P         QE        PE       EE (eV) \n"]
        debug = False
        # debug = True
        if debug:
            print(log[-1].strip())
            for i in range(sp.ntr):
                st = "NTR" + str(i).ljust(8)
                st += "{0:10.6f}".format(Nq[i])
                st += "{0:10.6f}".format(Np[i])
                st += "{0:10.6f}".format(au2ev * 0.5 * Nq[i] ** 2)
                st += "{0:10.6f}".format(au2ev * 0.5 * Np[i] ** 2)
                st += "{0:10.6f}".format(au2ev * 0.5 *
                                         (Np[i] ** 2 + Nq[i] ** 2))
                print(st)
        for i in range(sp.nm):
            e = (0.5* sp.w[sp.ntr + i]
                    * (Np[sp.ntr + i] ** 2 + Nq[sp.ntr + i] ** 2))
            st = "   " + str(i).ljust(2) + "  "+ "{0:6.1f}".format(sp.w[sp.ntr+i]*au2cm) + " "
            st += "  "+ "{0:8.4f}".format((e / sp.w[i + sp.ntr] - 0.5))
            st += "{0:10.6f}".format(Nq[sp.ntr + i])
            st += "{0:10.6f}".format(Np[sp.ntr + i])
            st += "{0:10.6f}".format(au2ev * 0.5 * sp.w[sp.ntr + i] * Nq[sp.ntr + i] ** 2)
            st += "{0:10.6f}".format(au2ev * 0.5 * sp.w[sp.ntr + i] * Np[sp.ntr + i] ** 2)
            st += "{0:10.6f}".format(e * au2ev)
            log += [st + "\n"]
            if debug:
                print(log[-1].strip())

        return log

    def InitializeSample(self,sa):
        """
        Initialize the sample by setting initial coordinates.

        This method sets the initial coordinates for the simulation.
        """
        sp = self.sp 
        sa.svv = zeros(sp.xxe.shape)
        sa.sxx = sp.xxe.copy()
        sa.sJ = 0.0
        sa.soR = np.eye(3)
        sa.siJ = np.zeros(3)
        sa.srpar = ['A', 0.0, 0.0, 0.0, 0.0, np.zeros(3)]
        sa.SampInfo = {
            'rot': {
                'J': 0,
                'ax': 'A',
                'jz': 0.0,
                'qjz': 0.0,
                'ang': [0.0, 0.0],
                'bet': 0.0,
                'gamm': 0.0,
                'vecJ': np.zeros(3),
                'energy': 0.0,
                'kenergy': 0.0,
                'jenergy': 0.0,
                'cenergy': 0.0,
                'svecJs': np.zeros(3),
                'svecJm': np.zeros(3),
                'svecJc': np.zeros(3),
                'svecJ0': np.zeros(3),
                'svecJ0s': np.zeros(3),
                'senergy': np.zeros(3),
                'sjz': 0.0,
                'sang': [0.0, 0.0],
                'sbet': 0.0,
                'sgamm': 0.0,
            }
        }
        return 

    def InitializeWorker(self,ii,**dic):
        """
        Initialize the sample by setting initial coordinates.

        This method sets the initial coordinates for the simulation.
        """
        sp = self.sp 
        sa = self.sample(ii,self) 
        if 'seed' in dic.keys():
          sa.seed = dic['seed']
        else:
          sa.seed = ii
        self.InitialDist(sa,seed=sa.seed)
        return sa


    def SampleTotMolAngMom(self,sa):
        J = int(ICDFscalar(sa.dist['rotJ']['cont']))
        if 'rot' not in sa.SampInfo.keys():
           sa.SampInfo['rot'] = {}
        sa.SampInfo['rot']['J'] = J
        sa.sJ = J
        return 
    # because the angular momentum and direction of the molecule should be conserved even as the molecule
    # is deformed (ie, vibrational motion), we can use the Jz and Energy from the equilibrium configuration
    # to define the angular momentum/angular velocity vector for the (potentially) deformed molecule.
    def SampleRigidRotorState(self,sa,**dic):
        """
        Sample a state for a rigid rotor.

        Returns:
        list: Log messages.
        """
        sp = self.sp 
        ip = self.ip 
        if 'rotJ' not in sa.dist.keys():
            log = ["    no rotational state\n"]
            return log
        if 'rot' not in sa.SampInfo.keys():
           sa.SampInfo['rot'] = {}
        #assumes we have sampled J
        J = abs(sa.sJ)
        eej = sp.J2c*J*(J+1)
    # If molecule is spherical top sample angles using isotropic distriobuition  
        if sp.rsym == 'spherical': 
           ax = 'A'
           alp,bet,gamm = ICDFsample(sa.dist['arotJz']['cont'])
           jz  = -1
           eek = 0.0
        # if symmetric top, then we pick projection analytically 
        if J == 0:
           ax = 'A'
           jz, bet, gamm, eek = [0.0]*4
           bet0,bets,sig = [0.0]*3
        elif sp.rsym == 'linear':
           ax = 'z'
           jz = 0.0
           bet = hpi
           gamm = ICDFsample(sa.dist['srotGam']['cont'])[0]
           eek = 0.0
        elif 'asym' not in sp.rsym:  # symmetric top 
           ax = 'z'
           jz = ICDFscalar(sa.dist['srotJz']['cont'], J)
           bet = np.arccos(jz/np.sqrt(J*(J+1)))
           gamm = ICDFsample(sa.dist['srotGam']['cont'])[0] 
           eek = sp.Jzc*jz**2
        else:   #if its asymmetric top : 
           id, eek,ax,jx2,jy2,jz2,jz,sig = ICDFsample(sa.dist['arotJz']['cont'],int(J))  # sample projection state
#           print('jzcomp = ', sqrt(jz2), jz)
           # sample polar coordinate IMPORTANT: it only samples the first quadrant 0-pi/2, need to sample the other symmetric halg :  
           if ax == 'A':  # if its an A irrep
              pars = [jx2,jy2,jz2] 
              alp,bet,gamm = ICDFsample(sa.dist['arotBing']['cont'],np.array(pars))   
              bet0,bets,sig = bet,sqrt(jz2), sqrt(jz2-jz**2) 
           else:
              pars = [jz,sig,J] 
              bet0 = np.arccos(jz/np.sqrt(J*(J+1)))
              bets = sig/(np.sqrt( J*(J+1) - jz**2  )+1e-8)
              bet = ICDFsample(sa.dist['arotBetB']['cont'],*pars)[0]
              jxy = [jx2,jy2] 
              gamm = ICDFsample(sa.dist['arotGam']['cont'],*jxy)[0]
        ee = eek + eej
        debug = False
        #gamm = 0.25*pi
        #debug = True
        if debug:
          print('AX = ',ax)
          print('GAM BET = ', gamm/pi, bet/pi, ' pi rad')
        # although we are rotating the body fixed coordinate, this corresponds to the alpha coordinate in the euler angles (since its rotating the J in the body fixed frame)
        vecJ = matmul(iq2R(np.array([0.0,bet,gamm])).T,z*np.sqrt(J*(J+1)))
        # classical energy: (different when beta has uncertainty) 
        cee = np.dot(np.array([sp.rA,sp.rB,sp.rC]),vecJ**2)
        #vecJ = matmul(iang2R(np.array([0.0,bet,gamm]),'eul').T,z*np.sqrt(J*(J+1)))
        if sp.Irep == 4:  # for strongly asymmetric tops we chose to sample along someother axis, so we need to rotate back... 
           if ax[-1] == 'x':
             vecJ = matmul(sp.xtozR.T,vecJ)
             if debug: 
              print('ROTATING VECTOR')
              print(matmul(sp.xtozR.T,np.eye(3)))
           elif ax[-1] == 'y': 
             vecJ = matmul(sp.ytozR.T,vecJ)
             if debug: 
              print('ROTATING VECTOR')
              print(matmul(sp.xtozR.T,np.eye(3)))
        # sets angular velocity from angular momentum:
        sa.srpar = [ax,ee,jz,bet,gamm,vecJ]
        sa.SampInfo['rot']['energy'] = ee 
        sa.SampInfo['rot']['kenergy'] = eek
        sa.SampInfo['rot']['jenergy'] = eej 
        sa.SampInfo['rot']['cenergy'] = cee 
        sa.SampInfo['rot']['ax']     = ax  
        sa.SampInfo['rot']['jz']     = jz
        if 'asym' in sp.rsym and J >0:
         sa.SampInfo['rot']['qjz']    = sqrt(jz2)
        else:
         sa.SampInfo['rot']['qjz']    = jz
        sa.SampInfo['rot']['ang']    = [bet,gamm]
        sa.SampInfo['rot']['bet']    = bet
        if 'asym' in sp.rsym and J >0:
          sa.SampInfo['rot']['asymd'] = [bet0,bets,sig]
          sa.SampInfo['rot']['idjz'] = id
        sa.SampInfo['rot']['gamm']   = gamm
        sa.SampInfo['rot']['vecJ']   = vecJ
        if debug: 
          print(ip.name,'VecJ = ', vecJ, 'J = ', J, 'Jz = ', jz )
          print(ip.name,'zero EE? = ', abs(np.dot(np.array([sp.rA,sp.rB,sp.rC]),vecJ**2)-ee)  )
          print(ip.name, ' quantum ee = ', ee,'classical e = ', np.dot(np.array([sp.rA,sp.rB,sp.rC]),vecJ**2))
        return []

    def SampleZVeloc(self,sa):
        """
        Sample z-velocity for the system.

        Returns:
        float: Z-velocity value.
        list: Log messages.
        """
        sp = self.sp 
        ip = self.ip 
        log = ["   :" + ip.name + " : \n"]
        if 'vel' not in sa.dist.keys():
            log += ["    No velocity sample \n"]
            return -1, log
        if 'cont' in sa.dist['vel'].keys():
          v = ICDFscalar(sa.dist['vel']['cont'])
        else:
          v = sa.dist['vel']['v']
        log += [ "    Molecular COM Velocity " + "{0:10.5f}".format(v)
            + ", Energy = " + "{0:10.5f}".format(0.5 * sum(sp.mass) * v**2 * au2ev) + " eV \n"]
        return v, log

    def SampleRotation(self,sa):
        sp = self.sp 
        ip = self.ip 
        if ip.ordist == 'pdf':
          for _ in range(max(1, int(ip.orthin))):
              ang = MCsample(sa.dist['ori']['cont'])
        elif ip.ordist == 'fixed':
          ang = np.array(sa.dist['ori']['fixed'], dtype=float)
        else:  
          ang = ICDFsample(sa.dist['ori']['cont']) 
#          print('name = ', self.ip.name, 'ang = ', ang)
        oR = iang2R(ang,ip.rotpar)
        sa.soR = oR 
        sa.sang = ang
        return

    def PrintAngMom(self,sa,mess):
        sp = self.sp 
        ip = self.ip 
        xx, vv =  sa.sxx.copy() - COM(sa.sxx, sp.mass).T,\
                  sa.svv.copy() - COM(sa.svv, sp.mass).T
        # move molecule to its eckart frame
        U = EckartFrameTrans(sp.xxe, xx, sp.mass)
        xx = matmul(U, xx.T).T
        vv = matmul(U, vv.T).T 
        RR = sp.c2m[:, : sp.ntr - 3]
        PP = matmul(RR, matmul(inv(matmul(RR.T, RR)), RR.T))
        vv = reshape(matmul(PP, reshape(vv, (sp.nd))), (sp.na, 3))
        # move to moment of Inertia frame:
        S2B, II, Is = iI(iX(xx), iX(xx), sp.mass)
        if abs(np.linalg.det(Is)) < 1e-12:
          iIs = pinv(Is)
        else:
          iIs = inv(Is)
        Lm = np.sum(cross(xx, vv).T * sp.mass, axis=1)
        Ls = matmul(U.T, Lm)
        log = [mess+ ip.name + "  Js : " + "".join(["{0:8.4f}".format(j) + " " for j in Ls]) +'\n'] 
        log += [mess+ ip.name + "  Jm : " + "".join(["{0:8.4f}".format(j) + " " for j in Lm]) +'\n'] 
        return log    

    def SetOrientat(self,sa,**dic):
        """
        Sample molecular orientation.

        Returns:
        list: Log messages.
        """
        sp = self.sp 
        ip = self.ip 
        debug = False
        # move molecule so that Eckart is in space frame 
        U = EckartFrameTrans(sp.xxe, sa.sxx, sp.mass)
        sa.sxx = matmul(U, sa.sxx.T).T
        sa.svv = matmul(U, sa.svv.T).T
        oR = sa.soR 
        # Molecule (eckart) to Sampled Rotation Matrix
        sa.svv = matmul(oR, sa.svv.T).T
        sa.sxx = matmul(oR, sa.sxx.T).T
        if debug:
          log += ["    Jvec was                    : "+ "".join(["{0:10.7f}".format(a) + " " for a in sa.srpar[-1] ]) + " \n"]
        # the sampled Jvec is used and also needs to be rotated
        sa.srpar[-1] = matmul(oR,sa.srpar[-1])
        sa.siJ = matmul(oR,sa.siJ)
        if debug:
          log += ["    Jvec is                     : "+ "".join(["{0:10.7f}".format(a) + " " for a in sa.srpar[-1] ]) + " \n"]
        if 'printlog' in dic.keys(): 
          log = [f"{ip.name:<{INFO_LABEL_WIDTH}} = sampled orientation\n"]
          if 'ori' not in sa.SampInfo.keys():
             sa.SampInfo['ori'] = {}
          xyz = list(iR2ang(sa.soR,'xyz'))
          eul = list(iR2ang(sa.soR,'eul'))
          if sp.na == 2:  # spin about a linear molecule axis is arbitrary
             xyz[2] = 0.0
             eul[2] = 0.0
          sa.SampInfo['ori']['xyz'] = xyz
          sa.SampInfo['ori']['eul'] = eul
          sa.SampInfo['ori']["alpha"] = sa.SampInfo['ori']['eul'][0]
          sa.SampInfo['ori']["beta"] = sa.SampInfo['ori']['eul'][1]
          sa.SampInfo['ori']["gamma"] = sa.SampInfo['ori']['eul'][2]
          # Backward-compatible names used by older histogram helpers.
          sa.SampInfo['ori']["phi"] = sa.SampInfo['ori']['eul'][0]
          sa.SampInfo['ori']["theta"] = sa.SampInfo['ori']['eul'][1]
          sa.SampInfo['ori']["chi"] = sa.SampInfo['ori']['eul'][2]
          log += [info_angle_vec("alpha,beta,gamma", eul, "molecular Euler")]
          log += [info_angle_vec("wx,wy,wz", xyz, "XYZ rotation angles")]
        return log

    def InertiaFrameTransform(self,sa):
        """
        Transform the coordinates to the inertia frame.

        Returns:
        numpy.ndarray: Transformed coordinates.
        """
        S2B, Ibm, Is = iI(iX(sa.sxx), iX(sa.sxx), sp.mass)
        RS = np.round(matmul(sp.Ro.T, S2B.T), 10)
        return RS

    # assumes mass frequency "dimensionless" scaled normal coordinates
    def SetHOVibrState(self, sa, Q, P):
        """
        Set harmonic oscillator vibrational state.

        Args:
        Q (numpy.ndarray): Normal coordinates for displacement.
        P (numpy.ndarray): Normal coordinates for momentum.
        """
        sp = self.sp 
        ip = self.ip 
        #x, p = matmul(sp.n2x.T, Q), matmul(sp.x2n, P) the one below turns out to be equivalent:
        x = matmul(sp.n2x.T, Q)
        p = matmul(sp.n2p.T, P)
        xx = reshape(x, (sp.na, 3))
        vv = (reshape(p, (sp.na, 3)).T /sp.mass).T
        xx -= COM(xx,sp.mass).T
        #vv -= COM(vv,sp.mass).T
        U = EckartFrameTrans(sp.xxe, sa.sxx, sp.mass)
        # move vectors to current (non) eckart frame:
        xx = matmul(U.T, xx.T).T
        vv = matmul(U.T, vv.T).T
        sa.sxx += xx  
        sa.svv += vv
        debug = False
        if debug:
          # Build velocity with Q=0 so geometry is pure equilibrium
          Qtest = np.zeros_like(Q)
          Ptest = np.random.randn(*P.shape)
          p_raw = (sp.n2p.T @ Ptest).reshape(sp.na,3)
          T_raw = (sp.mass[:,None]*p_raw).sum(0)
          R_raw = (sp.mass[:,None]*np.cross(sp.xxe, p_raw)).sum(0)
          print("Without /mass:  COM vel =", T_raw, "   Rot moment =", R_raw)
          p_div = (p_raw.T / sp.mass).T
          T_div = (sp.mass[:,None]*p_div).sum(0)
          R_div = (sp.mass[:,None]*np.cross(sp.xxe, p_div)).sum(0)
          print("With /mass:     COM vel =", T_div, "   Rot moment =", R_div)
          U = EckartFrameTrans(sp.xxe, sa.sxx, sp.mass)
          r_aligned = (U @ sp.xxe.T).T
          print("RMS after alignment:", np.sqrt(((r_aligned - sa.sxx)**2 * sp.mass[:,None]).sum()/sp.mass.sum()))
          quit()

    def CalcInstantRotJ(self,cJ):   
        sp = self.sp 
        #vector model 
        ww = solve(diag(sp.Ib),cJ)
        vv = cross(ww,sp.xxe)
        U = EckartFrameTrans(sp.xxe, xx, sp.mass)
        xx = matmul(U, xx.T).T
        # the actual instantaneous angular momentum then is: 
        Lm  = np.sum(cross(xx, vv).T * sp.mass, axis=1)  
        return Lm

    def SetAngularVeloc(self, sa, cJ):
        """
        Set angular velocity based on classical angular momentum.

        Args:
        cJ (numpy.ndarray): Classical angular momentum (in the eckart frame!!).
        """
        sp = self.sp 
        ip = self.ip 
        debug = False
        #debug = True
        xx =  sa.sxx.copy() - COM(sa.sxx, sp.mass).T
        #move to eckart...
        U = EckartFrameTrans(sp.xxe, xx, sp.mass)
        xx = matmul(U, xx.T).T
        # Two concievable approaches: Use the instantaneous moment of inertia...
        if False: 
         S2B, Ibm, Is = iI(iX(xx), iX(xx), sp.mass)
         ww = solve(Is,cJ)
         vv1 = cross(ww, xx)
        else:
        # or the rigid rotor with Equilibrium :
         #print('D = ')
         #print(diag(sp.Ib))
         #print('det = ', np.linalg.det(diag(sp.Ib)))
         if np.linalg.det(diag(sp.Ib)) == 0: 
          ww = svdsolve(diag(sp.Ib),cJ)
         else:
          ww = solve(diag(sp.Ib),cJ)
         vv1 = cross(ww,sp.xxe)
        sa.svren = sum(0.5*cJ**2/(sp.Ib+1e-10))
        #instantaneous angular momentum in eckart frame:  
        ccJ = np.sum(cross(xx, vv1).T * sp.mass, axis=1) 
        xx = matmul(U.T, xx.T).T
        vv1 = matmul(U.T, vv1.T).T 
        #instantaneous angular momentum in original frame:  
        LL = np.sum(cross(xx, vv1).T * sp.mass, axis=1)
        sa.svv += vv1
        if False: 
          xx, vv =  sa.sxx.copy() - COM(sa.sxx, sp.mass).T,\
                    sa.svv.copy() - COM(sa.svv, sp.mass).T
          Ls  = np.sum(cross(xx, vv).T * sp.mass, axis=1)
          print('Ls = ', Ls) 
          print('LL = ', LL) 
          # move molecule to its eckart frame
          U = EckartFrameTrans(sp.xxe, xx, sp.mass)
          xx = matmul(U, xx.T).T
          vv = matmul(U, vv.T).T
          L0  = np.sum(cross(sp.xxe, vv).T * sp.mass, axis=1)
          Lm  = np.sum(cross(xx, vv).T * sp.mass, axis=1)
          Lc  = Ls-L0
          print('Lm = ', Lm) 
          print('JJ = ', ccJ) 
        ax,ee,jz,bet,gamm, vecJ =  sa.srpar
        #the instantaneous angular momentum in original frame
        sa.siJ = np.sum(cross(sa.sxx, sa.svv).T * sp.mass, axis=1)
        sa.sren = ee 
        S2B, II, Is = iI(iX(xx), iX(xx), sp.mass)
        if abs(np.linalg.det(Is)) < 1e-12:
          iIs = pinv(Is)
        else:
          iIs = inv(Is)
        sa.scren = 0.5 * matmul(LL,matmul(iIs, LL))
        aax = {'x':0, 'y':1 ,'z':2, 'A':2}
        J = abs(sa.sJ)
        qJ = 0.5 * (-1 + sqrt(1 + 4.0 * norm(cJ) ** 2) )
        log = [f"{ip.name:<{INFO_LABEL_WIDTH}} = angular velocity from sampled rotor\n"]
        if J == 0 and norm(cJ) < 1.0e-10 and abs(sa.sren) < 1.0e-14 and abs(sa.svren) < 1.0e-14:
            log += [f"{'rotor state':<{INFO_LABEL_WIDTH}} = J = 0; angular momentum fixed at zero\n"]
            return log
        log += [info_scalar("quantum J", J, "", "{:14.1f}")]
        log += [info_scalar("quantum J_" + ax[-1], jz, "au", "{:14.3f}")]
        log += [info_scalar("equilibrium energy", sa.sren * au2ev, "eV", "{:14.5f}")]
        log += [info_vec("classical J", cJ, "au", "J", norm(cJ), qJ, "{:12.7f}")]
        log += [info_scalar("classical J_" + ax[-1], cJ[aax[ax[-1]]], "au", "{:14.3f}")]
        log += [info_scalar("classical energy", sa.svren * au2ev, "eV", "{:14.5f}")]
        log += [info_scalar("polar |J|", norm(cJ), "au", "{:14.6f}")]
        log += [info_scalar("polar beta", bet / pi, "pi rad", "{:14.6f}")]
        log += [info_scalar("polar gamma", gamm / pi, "pi rad", "{:14.6f}")]
        log += [info_vec("instantaneous J", cJ, "au", "J", norm(cJ), qJ, "{:12.7f}")]
        log += [info_scalar("instant. energy", sa.scren * au2ev, "eV", "{:14.5f}")]
        log += [info_vec("instant. I^-1", iIs[np.triu_indices_from(iIs)], "au", fmt="{:12.3e}")]
        return log   
 



    def CalcRotEner(self,sa):
        """
        Calculate the rotational energy of the system. s 

        Returns:
        list: Log messages.
        """
        sp = self.sp 
        ip = self.ip 
        log = [f"{ip.name:<{INFO_LABEL_WIDTH}} = {sp.rsym}, symmetry constant = {getattr(sp, 'asymk', 0.0):6.2f}\n"]
        asymk = getattr(sp, 'asymk', 0.0)
        if sp.na == 1:
         log += [f"{'rotational space':<{INFO_LABEL_WIDTH}} = none\n"]
         return log
        def inertia_tensor(r, masses):
           I = np.zeros((3,3))
           for (ri, m) in zip(r, masses):
               r2 = ri @ ri
               I += m * (r2 * np.eye(3) - np.outer(ri, ri))
           return I 

        debug = True
        debug = False
        xx, vv =  sa.sxx.copy() - COM(sa.sxx, sp.mass).T,\
                  sa.svv.copy() - COM(sa.svv, sp.mass).T
        xx0 = xx.copy()
        Ls  = np.sum(cross(xx, vv).T * sp.mass, axis=1)
        # move molecule to its eckart frame
        U = EckartFrameTrans(sp.xxe, xx, sp.mass)
        xx = matmul(U, xx.T).T
        vv = matmul(U, vv.T).T
        S2B, II, Is = iI(iX(xx), iX(xx), sp.mass)
        Lm  = np.sum(cross(xx, vv).T * sp.mass, axis=1)
        L0  = np.sum(cross(sp.xxe, vv).T * sp.mass, axis=1)
        E0  = 0.5*L0**2/(1.e-16+sp.Ib)
        Lc  = Lm-L0
        L0s = matmul(U.T,L0)
        # get projector
        #PP = matmul(RR, matmul(inv(matmul(RR.T, RR)), RR.T))
        #vv = reshape(matmul(PP, reshape(vv, (sp.nd))), (sp.na, 3))
        #vv = ProjectRTSpace(vv,sp.xxe,sp.xxe,sp.mass,'r')
        #vv = ProjectRTSpace(vv,xx,sp.xxe,sp.mass,'r')
        xx = matmul(S2B.T, xx.T).T
        vv = matmul(S2B.T, vv.T).T
        LL = np.sum(cross(xx, vv).T * sp.mass, axis=1)
        if any([e < 1e-12 for e in diag(II).tolist()]):
         ee = 0.5 * LL**2 / (1.0e-16 + diag(II))
        else:
         ee = 0.5 * LL**2 /  diag(II)
        EE = sum(ee)
        axx = ['x','y','z']
        rot_info = sa.SampInfo.get('rot', {})
        if 'ax' not in rot_info.keys():   
          ax = 2
        else:
          aa = rot_info['ax'][-1]
          if aa == 'z' or aa == 'A': 
           ax = 2 
          elif aa == 'x': 
           ax = 0 
          elif aa == 'y': 
           ax = 1
        nL = norm(L0)
        jz = L0[ax]
        if nL < 1e-5:
          bet, gamm = 0.0, 0.0
        else:
          bet = np.arccos(jz/nL)
          if ax == 2:
            gamm = np.arctan2(L0[1], -L0[0])
          elif ax == 1:
            gamm = np.arctan2(L0[0], -L0[2]) 
          elif ax == 0: 
            gamm = np.arctan2(L0[2], -L0[1]) 
        if 'rot' not in sa.SampInfo.keys(): 
         sa.SampInfo['rot'] = {} 
        sa.SampInfo['rot']["svecJs"]  = Ls  
        sa.SampInfo['rot']["svecJm"]  = Lm
        sa.SampInfo['rot']["svecJc"]  = Lc
        sa.SampInfo['rot']["svecJ0"]  = L0
        sa.SampInfo['rot']["svecJ0s"]  = L0s
        sa.SampInfo['rot']["senergy"] = ee
        sa.SampInfo['rot']["senergy_full"] = ee
        sa.SampInfo['rot']["senergy_vec"] = E0
        sa.SampInfo['rot']['sang']    = [bet,gamm]
        sa.SampInfo['rot']['sbet']    = bet
        sa.SampInfo['rot']['sgamm']   = gamm
        sa.SampInfo['rot']['sjz']     = jz
        nLs = norm(Ls) 
        nLm = norm(Lm)
        nLc = norm(Lc)
        nL0 = norm(L0)
        nL0s = norm(L0s)
        qLs  = 0.5*(-1+sqrt(1+4*nLs**2))
        qLm  = 0.5*(-1+sqrt(1+4*nLm**2))
        qLc  = 0.5*(-1+sqrt(1+4*nLc**2))
        qL0  = 0.5*(-1+sqrt(1+4*nL0**2))
        qL0s = 0.5*(-1+sqrt(1+4*nL0s**2))

        log += [info_vec("full J, space", Ls, "au", "J", nLs, qLs)]
        log += [info_vec("full J, Eckart", Lm, "au", "J", nLm, qLm)]
        log += [info_vec("vector J, space", L0s, "au", "J", nL0s, qL0s)]
        log += [info_vec("vector J, Eckart", L0, "au", "J", nL0, qL0)]
        log += [info_vec("vibrational J", Lc, "au", "J", nLc, qLc)]
        log += [info_vec("vector rot. energy", np.array(E0) * au2ev, "eV", "E", sum(E0) * au2ev)]
        log += [info_vec("full rot. energy", np.array(ee) * au2ev, "eV", "E", EE * au2ev)]
        log += [info_scalar("vector |J|", norm(L0), "au", "{:14.5f}")]
        log += [info_scalar("vector beta", bet / pi, "pi rad", "{:14.5f}")]
        log += [info_scalar("vector gamma", gamm / pi, "pi rad", "{:14.5f}")]
        log += [info_scalar("vector axis J" + axx[ax], jz, "au", "{:14.5f}")]
        return log



    def SampleHOVibrState(self,sa):
        """
        Generate a sample of harmonic oscillator vibrational states.

        This method generates a sample of harmonic oscillator vibrational states
        for a molecule, taking into account the vibrational frequencies and
        temperature. It returns the vibrational states and associated energies.

        Returns:
            list: Log messages containing information about generated states.
        """
        sp = self.sp 
        ip = self.ip 
        debug = True
        debug = False
        if debug:
            # samp = SampleMC(5000,array([0.01,-0.02]),HarmWigner,[0],idel=1.5)[0]
            samp = ICDFsample(sa.dist['wig']['cont'][0])[0]
            N = len(samp)
            for j in range(sp.nm):
                ee = 0.0
                Q, P = zeros(sp.nm), zeros(sp.nm)
                for i in range(N):
                    Q[j], P[j] = samp[i]
                    ee += 0.5 * sp.w[j + sp.ntr] * (Q[j] ** 2 + P[j] ** 2)
                print(j, " e = ",  au2cm * ee / float(N),  0.5 * au2cm * sp.w[j + sp.ntr])
            quit()
        log = ["  :" + ip.name + " : \n"]
        if 'nvib' not in sa.dist.keys():
            log += ["    no vibrational state\n"]
            return log
        vi = ICDFsample(sa.dist['nvib']['cont']).astype(int)
        #print('VI = ', vi, 'max = ', len(sa.dist['wig']['cont']))
        Q, P = zeros(sp.nd), zeros(sp.nd)
        sa.sven = zeros(sp.nm)
        log += ["    mode   freq  vstat    Q         P         QE        PE        EE (eV)      \n"]
        for j, ii in enumerate(vi, start=sp.ntr):
            i = int(ii)
            wi = ICDFsample(sa.dist['wig']['cont'][i])[0]
            if j - sp.ntr not in ip.nfreeze:
             Q[j], P[j] = wi
            sa.sven[j - sp.ntr] = 0.5 * sp.w[j] * (Q[j] ** 2 + P[j] ** 2)
            st = "    " + str(j - sp.ntr).ljust(3) + "  "+ \
                "{0:6.1f}".format(sp.w[j]*au2cm) + "  " + str(i).ljust(4)
            st += "{0:10.6f}".format(Q[j])
            st += "{0:10.6f}".format(P[j])
            st += "{0:10.6f}".format(au2ev * 0.5 * sp.w[+j] * Q[j] ** 2)
            st += "{0:10.6f}".format(au2ev * 0.5 * sp.w[+j] * P[j] ** 2)
            st += "{0:10.6f}".format(au2ev * 0.5 *
                                     sp.w[+j] * (P[j] ** 2 + Q[j] ** 2))
            log += [st + "\n"]
        if 'vib' not in sa.SampInfo.keys(): 
         sa.SampInfo['vib'] = {}
        sa.SampInfo['vib']['vi'] = vi 
        sa.SampInfo['vib']['qp'] = [Q,P]
        sa.SampInfo['vib']['P'] = P
        sa.SampInfo['vib']['Q'] = Q
        self.SetHOVibrState(sa, Q, P)
        return log

    def SetOrientatConvention(self, S2B, Ib, Rep):
        """
        Set the orientation convention for the molecule's rotational state.

        This method sets the orientation convention for the molecule's rotational state.
        The orientation convention is based on the shape of the molecule and its
        principal moments of inertia. It also adjusts the orientation transformation
        matrix S2B and the moment of inertia tensor Ib accordingly.

        Args:
            S2B (numpy.ndarray): The orientation transformation matrix.
            Ib (list): The principal moments of inertia [Ia, Ib, Ic].
            Rep (int): The representation code for the orientation convention.
        """
        # a b c are ordered with incresing 1/I eigenvalue (decreasing rotational constant)
        # the space-body frame rotation matrix S2B is fixed so that the figure axis z (the  symmetry axis) is pointing along Z
        # prolate (sausage) x y z (b c a), Ib~Ic > Ia
        sp = self.sp 
        if Rep == 1:  # b c a
            if sp.rB > sp.rC:
                sp.J2c = 0.5 * (sp.rB + sp.rC)
                sp.Jzc = sp.rA - 0.5 * (sp.rB + sp.rC)
                sp.Jcc = 0.25 * (sp.rB - sp.rC)
                Ib[0], Ib[1], Ib[2] = Ib[1], Ib[2], Ib[0]
                R = matmul(Rabout(hpi, 0), Rabout(hpi, 1))
            sp.S2B = matmul(R.T, S2B.T)
            sp.Irep = 1
        elif Rep == 2:  # c a b
            if sp.rC > sp.rA:
                sp.J2c = 0.5 * (sp.rC + sp.rA)
                sp.Jzc = sp.rB - 0.5 * (sp.rC + sp.rA)
                sp.Jcc = 0.25 * (sp.rC - sp.rA)
                Ib[0], Ib[1], Ib[2] = Ib[2], Ib[0], Ib[1]
                R = matmul(Rabout(hpi, 2), Rabout(hpi, 1))
            sp.S2B = matmul(R.T, S2B.T)
            sp.Irep = 2
        # oblate (frizbee) x y z ( a b c ), Ic > Ia~Ib
        elif Rep == 3:  # a b c
            sp.J2c = 0.5 * (sp.rA + sp.rB)
            sp.Jzc = sp.rC - 0.5 * (sp.rA + sp.rB)
            sp.Jcc = 0.25 * (sp.rA - sp.rB)
            sp.Irep = 3
            sp.S2B = S2B.T
            R = identity(3)
        # linear
        elif Rep == 0:
            sp.J2c, sp.Jzc, sp.Jcc = sp.rB, 0.0, 0.0
            Ib[0], Ib[1], Ib[2] = Ib[2], Ib[1], Ib[0]
            R = Rabout(hpi, 1)
            sp.S2B = matmul(R.T, S2B.T)
            sp.Irep = 0
        # stronly asymmetric top
        elif Rep == 4:
            sp.Irep = 4
            sp.J2c = (sp.rA + sp.rB +  sp.rC)/3.0
            sp.Jxc = sp.rA - sp.J2c  
            sp.Jyc = sp.rB - sp.J2c  
            sp.Jzc = sp.rC - sp.J2c  
            sp.Irep = 4
            sp.S2B = S2B.T
            R = identity(3)
        sp.Ib = Ib
        sp.Ro = R

    def StandardOrientat(self):
        """
        Standardize the orientation of the molecule.

        This method standardizes the orientation of the molecule based on its shape
        and principal moments of inertia. It also adjusts the orientation transformation
        matrix and other relevant properties.

        Returns:
            list: Log messages describing the standard orientation.
        """
        debug = True
        debug = False
        sp = self.sp 
        ip = self.ip 
        self.log.append("   #### Standard Orientation: \n")
        sp.x0 += -COM(sp.x0, sp.mass).T
        if sp.na == 1:
            self.log.append("   No standard orinetation: \n")
            sp.rsym = "atom"
            x0 = sp.x0
        else:
            S2B, Ibm, Is = iI(iX(sp.x0), iX(sp.x0), sp.mass)
            Ib = diag(Ibm).copy()
            sp.rA, sp.rB, sp.rC = 0.5 / (Ib + 1.0e-20)
            asymk = (2 * sp.rB - sp.rA - sp.rC) / (sp.rA - sp.rC)
            sp.asymk = asymk
            self.log.append( "   Asymmetry constant: " + "{0:14.9f}".format(asymk) + "\n")
            if abs(Ib[0]) < 0.001:
                sp.rsym = "linear"
                self.SetOrientatConvention(S2B, Ib, 0)
            elif asymk > 0.666:
            #elif asymk > 0.333:
                self.SetOrientatConvention(S2B, Ib, 3)
                if asymk - 0.98 > 0:
                    sp.rsym = "oblate"
                else:
                    sp.rsym = "asym-oblate"
            elif asymk < -0.666:
            #elif asymk < -0.333:
                self.SetOrientatConvention(S2B, Ib, 1)
                if asymk + 0.98 < 0:
                    sp.rsym = "prolate"
                else:
                    sp.rsym = "asym-prolate"
            else:
                if abs(asymk) < 0.02:
                    self.SetOrientatConvention(S2B, Ib, 2)
                    sp.rsym = "spherical"
                else:
                    self.SetOrientatConvention(S2B, Ib, 4)
                    sp.rsym = "asym-spherical"
            self.log.append("   Rotor type : " + sp.rsym + "\n")
            x0 = matmul(sp.S2B, sp.x0.T).T
        sp.xxe, sp.x0 = x0.copy(), x0.copy()
        # debugging standard orientation
        if debug:
          print('TYPE = ', sp.rsym, 'ASYM = ', asymk)
          S2B, II, Is = iI(iX(x0), iX(x0), sp.mass)
          print('II x0 = ')
          print(Is)
          print('S2B = ')
          print(S2B)
          R = Rabout(hpi, 1)
          x2 = matmul(R.T,x0.T).T
          U = EckartFrameTrans(sp.xxe, x2, sp.mass)
          print('MOL ',sp.el)
          print('U = ')
          print(U)
          print('R = ')
          print(R)
          print('prin = ', U[2,:])
          open('seed'+str(sp.na)+'.xyz','w').writelines(XYZlist(sp.el, x0) + XYZlist(sp.el, x2))
        self.log += ["   Standard Orientation :\n"]
        self.log += ["   " + ln for ln in XYZlist(sp.el, x0)]
        if hasattr(sp, "c2n"):
            if sp.na != 1 and hasattr(sp,'HH'):
             for i in range(sp.nd):
                 sp.c2n[:, i] = reshape(matmul(sp.S2B, reshape(sp.c2n[:, i], (sp.na, 3)).T).T,(sp.nd,))
             rt = GetRotTransVec(sp.x0, sp.mass, sp.el)
             sp.ntr = rt.shape[0]
             sp.nm = sp.nd - sp.ntr
             sp.c2n = sp.c2n[:, sp.ntr:]
             if ip.lowdin:
               PP = np.eye(sp.nd)-matmul(rt.T, rt)
               R = matmul(PP,sp.c2n) 
               # Lowdin orthonormalization in the space orthogonal to rt
               A,_,Vt = svd(R, full_matrices=False)
               sp.c2n = matmul(A,Vt)
             sp.c2n = np.concatenate([rt, sp.c2n.T]).T
             #sp.c2n = np.concatenate([rt, sp.c2n.T]).T
            # eigenvectors of mas-scalled hessian:
            sp.n2c = sp.c2n.T
            #transformation from mass-scaled to cartesian   
            sp.c2m, sp.m2c, sp.p2m, sp.m2p =  ScaleTransform2(sp.c2n, mass=sp.mass2)
            #transformation from mass and frequency scaled to cartesian
            sp.x2n, sp.n2x, sp.p2n, sp.n2p =  ScaleTransform2(sp.c2n, mass=sp.mass2,omega=sp.w)
            if debug:
             qout = [] 
             for n in range(sp.nm):
               Q, P = zeros(sp.nd), zeros(sp.nd) 
               for ii in np.arange(-3.0,3.0,0.05): 
                 Q[sp.ntr+n] = ii 
                 x = sp.xxe + reshape(matmul(sp.n2x.T, Q),(sp.na,3))
                 qout += XYZlist(sp.el,x)
             open('modesout'+ip.name+'.xyz','w').writelines(qout)
            ## normal coordinates
            #sp.c2m, sp.m2c = ScaleTransform(
            #    sp.c2n, sp.n2c, sp.w, sp.mass, +1, 0)
            ## classical momenta have the mass frequency reversed
            #sp.p2m, sp.m2p = ScaleTransform(
            #    sp.c2n, sp.n2c, sp.w, sp.mass, -1, 0)
            ## mass frequency scaled normal coordinates to
            #sp.x2n, sp.n2x = ScaleTransform(
            #    sp.c2n, sp.n2c, sp.w, sp.mass, +1, 1)
            #print('QZERO?', norm(sp.c2m2-sp.c2m), norm(sp.m2c2- sp.m2c),   norm(sp.p2m2-sp.p2m.T), norm(sp.m2p2- sp.m2p.T), norm(sp.x2n2-sp.x2n), norm(sp.n2x2- sp.n2x))
            #quit()

        debug = False
        #debug = True
        if debug:
            RR = GetRotTransVec(sp.x0, sp.mass, sp.el)
            print("ZERO?1 = ", norm(RR), diag(matmul(RR, sp.c2n[:, : sp.ntr])))
            print("ZERO?3 = ", norm(RR), diag(matmul(RR, sp.c2n[:, : sp.ntr])))
            RR = mscale2(GetRotTransVec(sp.x0, sp.mass, sp.el), sp.mass, -1)
            print("ZERO?5 = ", norm(RR), diag(matmul(RR, sp.c2m[:, : sp.ntr])))
            print("ZERO?6 = ", norm(RR), diag(matmul(RR, sp.c2m[:, : sp.ntr])))
        return

    def EstimateMaxR(self, T):
        """
        Estimate the maximum rotational state for a given temperature.

        This method estimates the maximum rotational state for a molecule at a
        given temperature based on the temperature and molecule's properties.

        Args:
            T (float): The temperature in Kelvin.

        Returns:
            int: The estimated maximum rotational state.
        """
        sp = self.sp 
        ip = self.ip 
        if T <= 0:
            return 1
        MaxR = 5
        ii = np.arange(5000)
        EE = sp.J2c * ii * (ii + 1)
        if ip.isotropic:
          rho = (2*ii+1)*Boltzmann(EE, T) 
        else:
          rho = Boltzmann(EE, T) 
        rho = rho/sum(rho) 
        for i in range(5000):
            if rho[i] < 0.001:
                MaxR = i 
                break 
        return int(MaxR)+1
  
    def AsymRigidRotorProjEnergies(self, MaxR):
        """
        Calculate the energies of the projection part of asymmetric rigid rotor states.

        This method calculates the projection energies of the rigid rotor states of the molecule
        up to a specified maximum state.

        Args:
            MaxR (int): The maximum rotational state.

        Returns:
            dict: List containing information about rotational states.
        """
        debug = True
        debug = False
        sp = self.sp 
        ip = self.ip 
        def D_matrix(j,ax,dd,wig):
            if ax == 'z':
             return np.eye(2*j+1)
            m_vals = np.arange(-j, j+1)
            dim = int(2*j + 1)
            D = zeros((dim, dim), dtype=complex)
            for i, mp in enumerate(m_vals):
                for k, m in enumerate(m_vals):
                    D[i, k] = dd[wig.Dindex(j,mp,m)]
            return D
        def PHamil(Jxc2,Jyc2,Jzc2):  
            if sp.rsym == 'asym-spherical':
              return sp.Jxc*Jxc2 + sp.Jyc*Jyc2 + sp.Jzc*Jzc2
            else:
              return sp.Jzc*Jzc2 + 2 * sp.Jcc*( Jxc2 - Jyc2 )
        # If spherical top :  just need Jcc  
        # If symmetric top :  Need Jcc and Jzc 
        # If asym  near oblate or near prolate, we stick to Z sampling.  We need: 
        # Jcc, discrete iCDF (from which to select eigenstate i)  for each eigenstate we need 
        #  <Jz> <Jz^2>  value and the sigma_Jz = sqrt(<Jz^2>- <Jz>^2) and <Jx^2> and <Jy^2> 
        # if asym is fairly spheroidal (neither prolate/oblate like), we also need: 
        #  the axis of rotation for each state (based on its irrep), and the rest.  
        #  
        if sp.rsym == 'oblate' or sp.rsym == 'prolate' or sp.rsym == 'spherical':
         return 
        if debug:
         np.set_printoptions(precision=4)
        RR = []
        if sp.rsym == 'asym-spherical':
          # y-> z, x-> y, z-> x   by uising Ry(bet=pi/2) Rz(gam=pi/2)
          # x-> z, y-> x, z-> y   by using Rx(-pi/2) Rz(-pi/2) in XYZ ->  Rz(alp=pi/2)Ry(bet=pi/2)Rz(gam=pi) for wigner needs to be in euler angles 
          if MaxR+1 > 500:
           print('Calculating Wigner Matrices...Will take a little while')
          wig = Wigner(MaxR+1)
          R = quat.array.from_euler_angles(pi/2,pi/2,pi)
          ddx = wig.D(R)
          R = quat.array.from_euler_angles(0.0,pi/2,pi/2)
          ddy = wig.D(R)
          ddz = 0 
          sp.ytozR = iq2R([0.0,hpi,hpi])
          sp.xtozR = iq2R([hpi,hpi,pi])
        evecs = []
        infos = []
        for j in range(MaxR+1):
            if debug:
              print('&&&&&&&&&&&&&&&&&&&& i = ', j)
            U, irreps,idx = WangTran(j)
            HC = PHamil(Wx2(j),Wy2(j),Wz2(j))
            HC = np.round(0.5*(HC +HC.conj().T),8) 
            P1 = { 'x':Px(j) , 'y':Py(j),  'z':Pz(j) }
            P2 = { 'x':Px2(j), 'y':Py2(j), 'z':Pz2(j)}
            if sp.rsym == 'asym-spherical':
               # get wigner matrices for rotating to different principal frames and the angular momentum operators in that new representation..
               # y-> z, x-> y, z-> x   by uising Ry(bet=pi/2) Rz(gam=pi/2)
               # x-> z, y-> x, z-> y   by using Rx(-pi/2) Rz(-pi/2) in XYZ ->  Rz(alp=pi/2)Ry(bet=pi/2)Rz(gam=pi) for wigner needs to be in euler angles 
               Dx = D_matrix(j,'x',ddx,wig)
               Dy = D_matrix(j,'y',ddy,wig)
               Dz = D_matrix(j,'z',ddz,wig)
               DD = { 'x':Dx    , 'y':Dy,     'z':Dz    }
               to, fm = {} , {}
               to['x'] = {'x':'z', 'y':'x' , 'z':'y'}
               to['y'] = {'x':'y', 'y':'z' , 'z':'x'}
               to['z'] = {'x':'x', 'y':'y' , 'z':'z'}
               fm['x'] = {'z':'x', 'x':'y' , 'y':'z'}
               fm['y'] = {'y':'x', 'z':'y' , 'x':'z'}
               fm['z'] = {'x':'x', 'y':'y' , 'z':'z'}
               Pn1, Pn2 = {} , {} # operator labelled in the new axis
               Po1, Po2 = {} , {} # operator labelled in the old axis 
               for ax in ['x','y','z']: # the original axis that was moved to the z coo 
                 Pn1[ax], Pn2[ax] = {}, {}
                 Po1[ax], Po2[ax] = {}, {}
                 if debug: 
                   print('##########  ORIGINAL AXIS : '+ax+' - > z' )
                 for ax1 in ['x','y','z']: # the new angular momentum operator 
                   Pn1[ax][ax1] = matmul( DD[ax] ,matmul(P1[fm[ax][ax1]],DD[ax].conj().T))
                   Pn2[ax][ax1] = matmul( DD[ax] ,matmul(P2[fm[ax][ax1]],DD[ax].conj().T))
                   Po1[ax][fm[ax][ax1]] = Pn1[ax][ax1] 
                   Po2[ax][fm[ax][ax1]] = Pn2[ax][ax1]
                   if debug: 
                     print('       NEW OPERATOR  for J2'+ax1, 'ERROR :', norm(Pn2[ax][ax1]-P2[ax1]) )
                     print('       NEW OPERATOR  for J'+ax1, 'ERROR :', norm(Pn1[ax][ax1]-P1[ax1]) )
               if debug:
                Pz2y, Pz2x, Pz2z =  matmul(Dy,matmul(Py2(j),Dy.conj().T)), matmul(Dx,matmul(Px2(j),Dx.conj().T)), Pz2(j)
                Pz1y, Pz1x, Pz1z =  matmul(Dy,matmul( Py(j),Dy.conj().T)), matmul(Dx,matmul( Px(j),Dx.conj().T)), Pz(j)

                Py2x, Pz2x, Px2x = matmul(Dx,matmul(Pz2(j),Dx.conj().T)), matmul(Dx,matmul(Px2(j),Dx.conj().T)), matmul(Dx,matmul(Py2(j),Dx.conj().T))   
                Px2y, Py2y, Pz2y = matmul(Dy,matmul(Pz2(j),Dy.conj().T)), matmul(Dy,matmul(Px2(j),Dy.conj().T)), matmul(Dy,matmul(Py2(j),Dy.conj().T))   
                pz1 = {'x':Pz1x , 'y':Pz1y , 'z':Pz1z , 'A':Pz1z}
                pz2 = {'x':Pz2x , 'y':Pz2y , 'z':Pz2z , 'A':Pz2z}
                print('DIS? ', norm(pz1['x']-Pn1['x']['z']), norm(pz2['x']-Pn2['x']['z']), norm(pz1['y']-Pn1['y']['z']), norm(pz2['y']-Pn2['y']['z']))
                #Hamiltonain on the new frame:
                HCr = {} 
                eig = []
                ei,ev = eigh(HC)
                for ax in ['x','y','z']:
                  HCr[ax] = PHamil(Po2[ax]['x'], Po2[ax]['y'], Po2[ax]['z'] ) 
                  eix,evx = eigh(HCr[ax])
                  eig.append(eix)
                  print(ax+' zero? = ', norm(eix-ei))
                print('Py2x = ',norm(Py2x-Po2['x']['z']))
                print('Pz2x = ',norm(Pz2x-Po2['x']['x']))
                print('Px2x = ',norm(Px2x-Po2['x']['y']))
                HCx = sp.Jzc * Py2x + 2 * sp.Jcc * ( Pz2x - Px2x )
                HCy = s/p.Jzc * Px2y + 2 * sp.Jcc * ( Py2y - Pz2y )
                print('HCXZ =', norm(HCr['x']-HCx))
                print('HCYZ =', norm(HCr['y']-HCy))
                print('ZETS? = ', norm(Py2x-Pn2['x']['y']) , norm(Pz2x-Pn2['x']['z']), norm(Px2x-Pn2['x']['x']) )
                print('DIFF D? =  DET = ', np.linalg.det(Dx), np.linalg.det(Dy))
                print('ROT zero? = ', norm(Py2(j)-Py2x), norm(Px2(j)-Px2x), norm(Pz2(j)-Pz2x), norm(Py2(j)-Py2y), norm(Px2(j)-Px2y), norm(Pz2(j)-Pz2y))
            blocks = []
            ee1 = []
            # diagonalise each wang sub-block and then rebuild 
            for i in range(len(idx)):
              if idx[i][1]-idx[i][0] > 0:
                el, ev = np.linalg.eigh(HC[idx[i][0]:idx[i][1],idx[i][0]:idx[i][1]] ) #+ 1e-10 * np.diag(np.arange(len(HC[idx[i][0]:idx[i][1],idx[i][0]:idx[i][1]]))) )
                blocks.append(ev)
                ee1 += el.tolist()
                if debug: 
                  print('####################')
                  print('  HC    = ')
                  print(HC[idx[i][0]:idx[i][1],idx[i][0]:idx[i][1]].real)
                  print(j,'Irrep   = ',irreps[idx[i][0]], idx[i][1]-idx[i][0]) 
                  print(' <x2>   = ')
                  print(Wx2(j)[idx[i][0]:idx[i][1],idx[i][0]:idx[i][1]].real)
                  print(eigh(Wx2(j)[idx[i][0]:idx[i][1],idx[i][0]:idx[i][1]].real)[0])
                  print(' <y2>   = ')
                  print(Wy2(j)[idx[i][0]:idx[i][1],idx[i][0]:idx[i][1]].real)
                  print(eigh(Wy2(j)[idx[i][0]:idx[i][1],idx[i][0]:idx[i][1]].real)[0])
                  print(' <z2>   = ')
                  print(Wz2(j)[idx[i][0]:idx[i][1],idx[i][0]:idx[i][1]].real)
                  print(eigh(Wz2(j)[idx[i][0]:idx[i][1],idx[i][0]:idx[i][1]].real)[0])
            #sort by eigenvalues and vectors:
            VV1 = block_diag(blocks)
            si = [u[1] for u in sorted(zip(ee1,np.arange(len(ee1))))]
            VV1 = VV1[:,si]
            ee1 = np.array([ee1[e] for e in si])
            sirreps = [irreps[s] for s in si]
            oo = [[i for i in range(len(ee1))],ee1,sirreps]
            ##Eigen in the Spherical Harmonic Representation:
            WS = matmul(U,VV1)
            #Eigen to Rotated Spherical Harmonic Representation:
            if sp.rsym == 'asym-spherical':
              WSR = {}
              WSR['x'] = matmul(Dx,WS)  
              WSR['y'] = matmul(Dy,WS)
              WSR['z'] = matmul(Dz,WS)
              WSR['A'] = WS
              if debug:
                HC2 = PHamil(Px2(j),Py2(j),Pz2(j))
                diav0 = diag(matmul(WS.conj().T,matmul(HC2,WS))).real
                for ax in ['x','y','z']: 
                  HCSRx = matmul(WSR[ax],matmul(diag(ee1),WSR[ax].conj().T))
                  print(ax+'ZOER0x = ', norm(HCSRx-HCr[ax]))
                  sm = 0.0
                  for i,e in enumerate(ee1):
                    diff = abs(norm(matmul(HCr[ax],WSR[ax][:,i])))- abs(e)
                    print(i,'vec  = ', diff)
                    sm += diff 
                  print('VECSERO  = ', sm)
                  print('diax2= ', norm(diag(matmul(WSR[ax].conj().T,matmul(HCr[ax],WSR[ax]))).real-diav0))
                  print(ax,' DIAG?') # should give the original z axis 
                  print(diag(matmul(DD[ax].conj().T,matmul(Po2[ax]['z'],DD[ax]))).real )
                  #print(matmul(WSR[ax].conj().T,matmul(Pn2[ax]['z'],WSR[ax])))
            if debug:
              HCS = matmul(WS,matmul(diag(ee1),WS.conj().T)) 
              HCS1 = matmul(VV1,matmul(diag(ee1),VV1.conj().T))
              print('ZERO1 = ',norm(HCS1-HC))
              HC2 = PHamil(Px2(j),Py2(j),Pz2(j))
              ei, vv = eigh(HC2)
              print('zero? = ', norm(array(sorted(ei))-array(sorted(ee1))),'sums = ', norm(HC2)- norm(HC), np.sum(np.sum(HC2))-np.sum(np.sum(HC)))
              print('zero? = ', norm(HC - matmul(VV1,matmul(diag(ee1),VV1.conj().T))))
              HC2 = np.round(0.5*(HC2 +HC2.conj().T),8) + 1e-10 * np.diag(np.arange(len(HC2)))
              print('ZERO2 = ', norm(HCS-HC2))
              sm = 0.0
              for i,e in enumerate(ee1):
                print(i,'vec  = ', norm(matmul(HC2,WS[:,i])), e)
                sm += abs(norm(matmul(HC2,WS[:,i]))-abs(e))
              print('VECSERO  = ', sm)
            def neghalve(j,vv):
#              pv = vv[:j]
              pv = np.hstack((vv[:j],[vv[j] * np.sqrt(0.5)]))
              pv = pv/(np.sqrt(np.dot(pv.conj(),pv))+1.e-15)
              pk1 = sum(pv.conj()*pv * -np.arange(j+1)[::-1])
              pk2 = sum(pv.conj()*pv * (-np.arange(j+1)[::-1])**2 )
              return pk1, pk2 
            def poshalve(j,vv):
#              pv = vv[j+1:]
              pv = np.hstack(([vv[j]*np.sqrt(0.5)], vv[j+1:]))
              pv = pv/(np.sqrt(np.dot(pv.conj(),pv))+1.e-15)
              pk1 = sum(pv.conj()*pv * np.arange(j+1))
              pk2 = sum(pv.conj()*pv * np.arange(j+1)**2 )
              return pk1, pk2 
            def fullproj(j,vv):
              pv = vv
              pv = pv/(np.sqrt(np.dot(pv.conj(),pv))+1.e-15)
              pk1 = sum(pv.conj()*pv * np.arange(-j,j+1))
              pk2 = sum(pv.conj()*pv * np.arange(-j,j+1)**2 )
              return pk1, pk2 
            ## expectation of the projection half-coo and squared
            k1, k2 = {}, {}
            ak1, ak2 = {}, {}
            if sp.rsym == 'asym-spherical':  #ak1 and ak2 are the half domain expectations about z only k1, k2 are the expectations on the eigenstates in the spherical harmonic basis (rotated if need be)
              for ax in ['x','y','z','A']:
               if ax == 'A':
                aa = 'z'
               else:
                aa = ax
               ak1[ax], ak2[ax] = [], []
               tk2, dk2 = [] , []
               for i in range(2*j+1):
                 pk1, pk2 = poshalve(j,vv = WSR[ax][:,i])
                 _, t2 = fullproj(j,vv = WSR[ax][:,i])
                 _, d2 = neghalve(j,vv = WSR[ax][:,i]) 
                 ak1[ax].append(pk1.real)
                 ak2[ax].append(pk2.real)
                 tk2.append(t2.real)
                 dk2.append(d2.real)
               k1[ax], k2[ax] = {}, {}
               for ax1 in ['x','y','z']: # angular mometnum in the new axis... note that Pn1[ax][x] = Px and so on.. I wasnt 100% it would be :D 
                 k1[ax][ax1] = diag(matmul(WSR[ax].conj().T,matmul(Pn1[aa][ax1],WSR[ax]))).real +1e-10
                 k2[ax][ax1] = diag(matmul(WSR[ax].conj().T,matmul(Pn2[aa][ax1],WSR[ax]))).real +1e-10
               np.set_printoptions(precision=4)
               if debug:
                 print(ax,' dK2 = ', dk2)
                 print(ax,'AK2 = ', ak2[ax])
                 print(ax,' K2 = ',  k2[ax]['z'], 'zero? = ', norm(k2[ax]['z']-tk2))
              akk1, akk2, asig = [], [], []
              kk1, kk2, sig = {'x':[],'y':[],'z':[]}, {'x':[],'y':[],'z':[]}, {'x':[],'y':[],'z':[]}
              for i in range(len(ee1)):
                ax = oo[-1][i][-1]
                for ax1 in ['x','y','z']:  # the full kx1, kx2, ky1, ky2 .... for the correct new frame ax
                  kk1[ax1].append(  k1[ax][ax1][i]  )   
                  kk2[ax1].append(  k2[ax][ax1][i]  )   
                  sig[ax1].append(sqrt( kk2[ax1][-1]- kk1[ax1][-1]**2 ))
                evecs.append('J'+str(j) + ' ev-' +str(i) + ' ax-' + ax +  ' : ' + "".join(["{0:8.4f}".format(ak) + " " for ak in  WSR[ax][:,i].tolist() ]) +'\n') 
                infos.append('J'+str(j) + ' ev-' +str(i) + ' ax-' + ax )
                infos[-1] += ' j2x: '+"{0:8.4f}".format(kk2['x'][0]) + ' j2y: '+ "{0:8.4f}".format(kk2['y'][0]) + ' j2z: '+"{0:8.4f}".format(kk2['z'][0]) +'\n'
                akk1.append(  ak1[ax][i] )   
                akk2.append(  ak2[ax][i] )   
                asig.append(sqrt( akk2[-1]- akk1[-1]**2 ))
              if debug:
                print('FINAL AK ')
                print(oo[1])
                print('   AK1 = ', akk1)
                print('   AK2 = ', akk2)
                print('   sig = ', asig)
            else:
              akk1, akk2, asig = [], [], []
              kk1, kk2, sig = {'x':[],'y':[],'z':[]}, {'x':[],'y':[],'z':[]}, {'x':[],'y':[],'z':[]}
              for i in range(2*j+1):
                pk1, pk2 = poshalve(j,WS[:,i])
                akk1.append(pk1.real)
                akk2.append(pk2.real)
              for ax in ['x','y','z']:
                kk1[ax] = diag(matmul(WS.conj().T,matmul( P1[ax],WS))).real +1e-10
                kk2[ax] = diag(matmul(WS.conj().T,matmul( P2[ax],WS))).real +1e-10
                sig[ax].append(sqrt( kk2[ax][-1]- kk1[ax][-1]**2 ))
            oo += [kk2['x'],kk2['y'],kk2['z'],akk1,asig]

            if debug:
              print('IRR = ', irreps) 
              print('irr = ', oo[1])
              print('aKK1 = ', akk1)
              print('aKK2 = ', akk2)
              print('KK2? = ', kk2['z'])
              print('aSIG = ', asig)
              for ax in ['x','y','z']: 
                print('KK1 = '+ax, kk1[ax])
                print('KK2 = '+ax, kk2[ax])
                print('SIG = '+ax, sig[ax])
            RR.append( oo )
        # if the axis was returned: we are returning the <J_i> in their NEW labels, so they need to be sampled like they were z and then at the end rotate
        # the IRREP will tell us which rotation was used, so we have to perform the reverse.
        # returns [ i, EE, IRREP, <J^2_x>, <J^2_y>, <J^2_z>, <J_z>+, sig^+_z ]
        sp.RR = RR
        asym_dir = os.path.join(self.ip.diagdir, "asym_rotor")
        os.makedirs(asym_dir, exist_ok=True)
        open(os.path.join(asym_dir, 'asym-rotor-'+self.ip.name+'-evecs.dat'),'w').writelines(evecs)
        open(os.path.join(asym_dir, 'asym-rotor-'+self.ip.name+'-info.dat'),'w').writelines(infos)
        return 

    def EstimateMaxVI(self):
        """
        Estimate the maximum value MaxV using nested loops and a condition based on Boltzmann probability.

        This function iterates over two loops to estimate the maximum value MaxV based on certain conditions.

        Returns:
        int: The estimated maximum value MaxV.
        """
        sp = self.sp 
        ip = self.ip 
        MaxV = 0
        i = 0 
        while True: 
          rho = VibBoltzmann(sp.w[sp.ntr:sp.ntr+sp.nm],i*np.ones(sp.nm),ip.Tvib)
          if not any(rho > 0.01) :
            MaxV = i
            break 
          i+=1 
        return MaxV

    def InitialDist(self,sa,**dic):
        """
        Generate an initial distribution of samples for the molecule.

        This method generates samples for various properties of the molecule, including vibrational states, rotational states,
        orientation, and velocity, based on the specified parameters.

        Args:
        Nsamp (int): The number of samples to generate.

        Returns:
        list: Log messages indicating the distribution generation process.
        """
        sp = self.sp 
        ip = self.ip 
        log = []
        log += [f"Generating Marcov Chain Parameters For Distribuitions (molecule :{ip.filename}\n"]
        if 'nsamp' in dic.keys(): 
          nsamp = dic['nsamp']
        else:
          nsamp = 0
        # generate vibrational dist:
        T = ip.Tvib
        sa.dist = {}
        if 'seed' in dic.keys():
          seed = dic['seed']
        else:
          seed = sa.seed
       
        if sp.nm > 0 and T >= 0:
            log += ["Generated vibrational state distro temperature " +   str(T) + "\n"]
            ip.MaxV = self.EstimateMaxVI()
            log += ["Maximum Excited Vibrational state "+str(ip.MaxV)+"\n"]
            for i in range(sp.nm):
              h  = np.array([exp(-sp.w[sp.ntr+i]*n/(kboltz*T)) for n in range(ip.MaxV+1)])
              h  = h/sum(h)
              log += [f'Mode {i} = '+"".join(["{0:6.3f}".format(v).rjust(8) for v in h])]
        if ip.MaxV >= 0:
            # generate some samples for vib states, assumes first sp.ntr vibrational frequencies are zero...
            sa.dist['nvib'] = {}
            sa.dist['nvib']['MaxV'] = ip.MaxV
            sa.dist['nvib']['cont'] = InitICDF(sp.nm,VibPartFuncICDF,
                                        [sp.w[sp.ntr:], T, False,ip.MaxV],seed=seed*607)
            if nsamp != 0:
               vvv = [ICDFsample(sa.dist['nvib']['cont']) for _ in range(nsamp)]
               for i in range(sp.nm):
                  vv = [samp[i] for samp in vvv]
                  hist_emit(vv, f"vi_mode{i}", stage="initial", scope=f"molecule_m{ip.mi}")
                  mx = int(max(vv))
                  hist, edg = np.histogram(vv, bins=mx+1)
                  hist = hist/sum(hist)
                  hedg = 0.5 * (edg[1] - edg[0])
                  hh  = np.array([exp(-sp.w[sp.ntr+i]*n/(kboltz*T)) for n in range(ip.MaxV+1)])
                  hh  = hh/sum(hh)
                  log += ['hist  = '+"".join(["{0:6.3f}".format(v).rjust(8) for v in hist[:ip.MaxV+1]])] 
                  log += ['exac  = '+"".join(["{0:6.3f}".format(v).rjust(8) for v in hh])]
               sa.dist['nvib']['cont'] = InitICDF(sp.nm,VibPartFuncICDF,
                                           [sp.w[sp.ntr:], T, False,ip.MaxV],seed=seed*439)
            log += [f"Generating distributions for leading-Wigner vibrational distributions up to excited state {ip.MaxV} \n"]
            # generate some samples from each leading-Wigner distribution for each vib state
            sa.dist['wig'] = {}
            sa.dist['wig']['MaxV'] = ip.MaxV
            sa.dist['wig']['cont'] = [] 
            for m in range(ip.MaxV + 1):
                sa.dist['wig']['cont'].append( InitICDF(1,HusimiFuncICDF,[m],seed=seed*983) )
            nsamp = 0
            #nsamp = 100000
            if nsamp != 0:
               for m in range(ip.MaxV + 1):
                   vv = [ICDFsample(sa.dist['wig']['cont'][m]) [0] for i in range(nsamp)]
                   open('husi'+str(m)+'.dat','w').writelines(["{0:6.3f}".format(v[0]).rjust(8)+',\n' for v in vv])
                   hist_emit(vv, f"iwig_state{m}", stage="initial", scope=f"molecule_m{ip.mi}")
                   ee, pp, qq, pq, pp2, qq2 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                   for v in vv:
                     ee += 0.5*(v[0]**2+v[1]**2)
                     pp2 += (v[0]**2)
                     qq2 += (v[1]**2)
                     pp += (v[0])
                     qq += (v[1])
                     pq += (v[1]*v[0])

                   ee, pp, qq, pq, pp2, qq2 = ee/float(nsamp), pp/float(nsamp), qq/float(nsamp), pq/float(nsamp), pp2/float(nsamp), qq2/float(nsamp)
                   log += ['Wigner expectation values for vibrational states = \n']
                   log += ['m = '+str( m)+ ' ee = '+str( ee)+ ' pp = '+str( pp)+' qq = '+str( qq)+' pq = '+str( pq)+' pp2 = '+str( pp2)+ ' qq2 = '+str( qq2)+'\n']
                   print('m = ', m, ' ee = ', ee, ' pp = ', pp, ' qq = ', qq, ' pq = ', pq, ' pp2 = ', pp2, ' qq2 = ', qq2)
               for m in range(ip.MaxV + 1):
                   sa.dist['wig']['cont'].append( InitICDF(1,HusimiFuncICDF,[m],seed=seed*661) )

        else:
            log += ["No vibrational state distro" ]
        # generate rotational dist:
        T = ip.Trot
        if sp.na > 1 and T >= 0:
            log += ["Generated rotational state distro temperature " + str(T) + "\n"]
            if T > 0:
                sa.dist['rotJ'] = {} 
                #ip.MaxR = self.EstimateMaxR(ip.Trot)
                MaxR = ip.MaxR
                sa.dist['rotJ']['MaxR'] = ip.MaxR
                log += ["Maximum rigid rotor state " + str(MaxR) + "\n"]
                # The field-free J manifold always carries the 2J+1
                # projection-state count. Orientation PDFs bias Euler angles,
                # not the underlying molecular J state count.
                func = IsoRotorTotBoltzICDF
                sa.dist['rotJ']['cont'] = InitICDF(1,func,[T,sp.J2c,ip.MaxR ],seed=seed*17)
                if nsamp != 0:
                  vv = np.round([ICDFsample(sa.dist['rotJ']['cont']) for _ in range(nsamp)])
                  hist_emit(vv, "j", stage="initial", scope=f"molecule_m{ip.mi}")
                  hist, edg = np.histogram(vv, bins=MaxR)
                  log += ['Total molecular angular momentum J Histogram  = \n']
                  log += [str(hist) +'\n']
                  log += [str(edg) + '\n']
                  sa.dist['rotJ']['cont'] = InitICDF(1,func,[T,sp.J2c,ip.MaxR ],seed=seed*283)
                # the molecule is a symmetric top, we sammple analytically from vector model... 
                #if True:

                if sp.rsym == 'spherical':  #spherical tops 
                  rn1,rn2 = array([tpi,1.0,tpi]), array([pi,0.0,pi])
                  sa.dist['arotJz'] = {}
                  sa.dist['arotJz']['cont'] = InitICDF(1,IsotropicDistICDF,[ip.rotpar,rn1,rn2],seed=seed*43)
                  if nsamp != 0:
                    vvv = np.round([ICDFsample(sa.dist['arotJz']['cont']) for _ in range(nsamp)])
                    aa = ['alpha','beta','gamma']
                    for i in range(3):
                      vv = [v[i] for v in vvv]
                      hist_emit(vv, f"iarotJz_{aa[i]}", stage="initial", scope=f"molecule_m{ip.mi}")  
                    sa.dist['arotJz']['cont'] = InitICDF(1,IsotropicDistICDF,[ip.rotpar,rn1,rn2],seed=seed*431)
                if sp.rsym == 'linear':
                   sa.dist['srotGam'] = {}
                   sa.dist['srotGam']['cont'] = InitICDF(1,uniform,[-pi,pi],seed=seed*433)
                   if nsamp != 0: 
                     vv = [ICDFsample(sa.dist['srotGam']['cont']) for _ in range(nsamp)]
                     hist, edg = np.histogram(vv, bins=20)
                     hist_emit(vv, "linear_rot_gamma", stage="initial", scope=f"molecule_m{ip.mi}")
                     log += ['Linear rotor transverse gamma histogram = \n']
                     log += [str(hist) +'\n']
                     log += [str(edg) + '\n']
                     sa.dist['srotGam']['cont'] = InitICDF(1,uniform,[-pi,pi],seed=seed*677)
                elif 'asym' not in sp.rsym: #  symmetric tops 
                   sa.dist['srotJz'] = {}
                   sa.dist['srotJz']['cont'] = InitICDF(1,SymRotorProjBoltzICDF,[T,sp.Jzc],par0=[0],seed=seed*73)
                   if nsamp != 0:
                     for iJ in range(MaxR):
                       vv = [ICDFsample(sa.dist['srotJz']['cont'],iJ) for _ in range(nsamp)]
                       hist_emit(vv, f"j{iJ}_jz", stage="initial", scope=f"molecule_m{ip.mi}")
                       hist, edg = np.histogram(vv, bins=21)
                       log += ['Rotational Jz Histogram (for J:'+str(iJ)+')  = \n']
                       log += [str(hist) +'\n']
                       log += [str(edg) + '\n']
                     sa.dist['srotJz']['cont'] = InitICDF(1,SymRotorProjBoltzICDF,[T,sp.Jzc],par0=[0],seed=seed*941)
                   sa.dist['srotGam'] = {}
                   sa.dist['srotGam']['cont'] = InitICDF(1,uniform,[-pi,pi],seed=seed*433)
                   if nsamp != 0: 
                     vv = [ICDFsample(sa.dist['srotGam']['cont']) for _ in range(nsamp)]
                     hist, edg = np.histogram(vv, bins=20)
                     log += ['Azimuthal phi angle Histogram  = \n']
                     log += [str(hist) +'\n']
                     log += [str(edg) + '\n']
                     sa.dist['srotGam']['cont'] = InitICDF(1,uniform,[-pi,pi],seed=seed*223)
                # the molecule is an asymmetric top, we sammple the eigenvalues 
                else:
                   #example dummy for J=10, ei eigenstate
                   Ji, ei = 0, 0
                   ki = 0 
                   sig3 = np.array([o[ei] for o in sp.RR[Ji][3:6]]) #[jx2,jy2,jz2]
                   jxy = [o[ei] for o in sp.RR[Ji][3:5]]
                   pars = [o[ei] for o in sp.RR[Ji][-2:]] + [Ji]
                   ll = sp.RR[Ji][1].tolist()
                   # asymmetric top list information:  ee,ax,jx2,jy2,jz2,jz,sig 
                   # for B irrep Wang eignestates:
                   sa.dist['arotJz'] = {}
                   sa.dist['arotJz']['cont'] = InitICDF(1,AsymRotorProjBoltzICDF,[T, sp.RR],par0=[Ji],seed=seed*127)
                   if nsamp != 0:
                     for iJ in range(MaxR):
                       vvv = [ICDFsample(sa.dist['arotJz']['cont'],iJ) for _ in range(nsamp)]
                       vv = [v[0] for v in vvv]
                       hist_emit(vv, f"j{iJ}_jz_en", stage="initial", scope=f"molecule_m{ip.mi}")
                       hist, edg = np.histogram(vv, bins=20,density=True)
                       log += ['Rotational Jz Energy Histogram (for J:'+str(iJ)+') = \n']
                       log += [str(edg) + '\n']
                       log += [str(hist) +'\n']
                       vv = [v[-2] for v in vvv]
                       hist_emit(vv, f"j{iJ}_jz", stage="initial", scope=f"molecule_m{ip.mi}")
                       hist, edg = np.histogram(vv, bins=20,density=True)
                       log += ['Rotational Jz Histogram (for J:'+str(iJ)+') = \n']
                       log += [str(edg) + '\n']
                       log += [str(hist) +'\n']
                       sa.dist['arotJz']['cont'] = InitICDF(1,AsymRotorProjBoltzICDF,[T, sp.RR],par0=[Ji],seed=seed*317)
                       vv = [v[0] for v in vvv]
                       hist_emit(vv, f"j{iJ}_jz_id", stage="initial", scope=f"molecule_m{ip.mi}")
                       hist, edg = np.histogram(vv, bins=20,density=True)
                       log += ['Rotational Jz Id Histogram (for J:'+str(iJ)+') = \n']
                       log += [str(edg) + '\n']
                       log += [str(hist) +'\n']
                       sa.dist['arotJz']['cont'] = InitICDF(1,AsymRotorProjBoltzICDF,[T, sp.RR],par0=[Ji],seed=seed*727)

                   sa.dist['arotBetB'] = {}
                   sa.dist['arotBetB']['cont'] = InitICDF(1,reject_GaussSin,[],par0=pars ,seed=seed*409)
                   #extract dummy Jx2 and Jy2
                   sa.dist['arotGam'] = {}
                   sa.dist['arotGam']['cont'] = InitICDF(1,reject_phi,[],par0=[*jxy],seed=seed*401)
                   #for A irrep Wang eigenstates
                   sa.dist['arotBing'] = {}
                   sa.dist['arotBing']['cont'] = InitICDF(1,reject_bingham, [],par0=[sig3],seed=seed*311)
                   if nsamp != 0:

                     #example for J=Ji, ei eigenstate
                     for iJ in range(MaxR):
                       idx = [o for o in sp.RR[iJ][0]]
                       pars = [o for o in sp.RR[iJ][-2:]] + [iJ]  #[jz,sig,J]
                       sigma3 = np.array([o for o in sp.RR[iJ][3:6]]) #[jx2,jy2,jz2]
                       axx = [o for o in sp.RR[iJ][2]]  
                       for i,jz in enumerate(pars[0]): 
                         sigz,ax, sig3 = pars[1][i], axx[i], [o[i] for o in sigma3]
                         log += ['Asymmetric vector model with proj. j-'+ax+'  '+str(jz)+' '+str(idx[i])+' \n']
                         if 'B' in ax:
                           vv = [ICDFsample(sa.dist['arotBetB']['cont'],*[jz,sigz,J]) for _ in range(nsamp)]
                           hist_emit(vv, f"iarotBetB_J{iJ}_{idx[i]}", stage="initial", scope=f"molecule_m{ip.mi}")
                           hist, edg = np.histogram(vv, bins=20)
                           bets = sp.RR[iJ][-1][i]/sqrt( iJ*(iJ+1) - sp.RR[iJ][-2][i]**2 )
                           bet0  = np.arccos(sp.RR[iJ][-2][i]/sqrt(iJ*(iJ+1)))
                           log += ['Rotational polar angle Histogram  = \n']
                           log += ['bet0     = '+ str(bet0) +'\n'] 
                           log += ['betsigma = '+ str(bets) + 'calculated : ' +"{0:8.4f}".format(np.std(vv)) + '\n']
                           log += [str(edg) + '\n']
                           log += [str(hist) +'\n']
                           jxy = [o[i] for o in sp.RR[iJ][3:5]]
                           vv = [ICDFsample(sa.dist['arotGam']['cont'],*jxy) for _ in range(nsamp)]
                           hist_emit(vv, f"iarotGam_J{iJ}_{idx[i]}", stage="initial", scope=f"molecule_m{ip.mi}")
                           hist, edg = np.histogram(vv, bins=20)
                           log += ['Rotational azimuthal angle Histogram  = \n']
                           log += ['J_x^2, J_y^2  = '+ str(jxy[0])+ '  ' + str(jxy[1])  +'\n' ] 
                           log += [str(edg) + '\n']
                           log += [str(hist) +'\n']
                         else:
                           vv = [ICDFsample(sa.dist['arotBing']['cont'],sig3) for _ in range(nsamp)]
                           vx,vy,vz = [[o[i] for o in vv] for i in range(3)]
                           hist_emit(vx, f"iarotAx_J{iJ}_{idx[i]}", stage="initial", scope=f"molecule_m{ip.mi}")
                           hist_emit(vy, f"iarotAy_J{iJ}_{idx[i]}", stage="initial", scope=f"molecule_m{ip.mi}")
                           hist_emit(vz, f"iarotAz_J{iJ}_{idx[i]}", stage="initial", scope=f"molecule_m{ip.mi}")
                           varx, vary,varz = np.var(vx), np.var(vy), np.var(vz)
                           log += ['Rotational Histogram  = \n']
                           log += ['J_x^2, J_y^2, J_z^2 0= '+ str(varx)+ '  ' + str(vary) + '  ' + str(varz)  +'\n' ] 
                           log += ['J_x^2, J_y^2, J_z^2 1= '+ str(sig3[0])+ '  ' + str(sig3[1]) + '  ' + str(sig3[2])  +'\n' ] 
                     sig3 = np.array([o[ei] for o in sp.RR[Ji][3:6]]) #[jx2,jy2,jz2]
                     jxy = [o[ei] for o in sp.RR[Ji][3:5]]
                     pars = [o[ei] for o in sp.RR[Ji][-2:]] + [Ji]
                     sa.dist['arotBetB']['cont'] = InitICDF(1,reject_GaussSin,[],par0=pars ,seed=seed*541)
                     sa.dist['arotGam']['cont'] = InitICDF(1,reject_phi,[],par0=[*jxy],seed=seed*211)
                     sa.dist['arotBing']['cont'] = InitICDF(1,reject_bingham, [],par0=[sig3],seed=seed*241)
            else:
                sa.rotsamp = [0 for i in range(nsamp)]
        else:
            log += ["No rotational state distro" ]
        # generate orientation 
        if sp.na > 1:
            log += ["Generated orientational distribuition type " +str(ip.ordist) + "\n"]
            # fix specific orientation
            if ip.ordist == 'fixed':
                sa.dist['ori'] = {}
                sa.dist['ori']['fixed'] = array(ip.orpars, dtype=float)
                if nsamp != 0:
                  vv = [sa.dist['ori']['fixed'] for _ in range(nsamp)]
                  valpha,vbeta,vgamma = [[o[i] for o in vv] for i in range(3)]
                  hist_emit(valpha, "ori_alpha", stage="initial", scope=f"molecule_m{ip.mi}")
                  hist_emit(vbeta,  "ori_beta", stage="initial", scope=f"molecule_m{ip.mi}")
                  hist_emit(vgamma, "ori_gamma", stage="initial", scope=f"molecule_m{ip.mi}")
            elif ip.ordist == 'pdf':  # polarization distribiuition
                orF = load_function_from_file(ip.orfilename,ip.orfunction)
                sa.dist['ori'] = {}
                sa.dist['ori']['cont'] = InitMC(
                    array([0.0, hpi, 0.0]),
                    OrientationPDFSurface,
                    [orF, ip.orpars],
                    domains=[[-pi, pi, True],[0.0, pi, False], [-pi, pi, True]],
                    idel=0.25,
                )
                if nsamp != 0: 
                  vv = []
                  for _ in range(nsamp):
                      for _ in range(max(1, int(ip.orthin))):
                          ang = MCsample(sa.dist['ori']['cont'])
                      vv.append(ang)
                  valpha,vbeta,vgamma = [[o[i] for o in vv] for i in range(3)]
                  hist_emit(valpha, "ori_alpha", stage="initial", scope=f"molecule_m{ip.mi}")
                  hist_emit(vbeta,  "ori_beta", stage="initial", scope=f"molecule_m{ip.mi}")
                  hist_emit(vgamma, "ori_gamma", stage="initial", scope=f"molecule_m{ip.mi}")
                  hist, edg = np.histogram(vv, bins=20)
                  log += ['Orientation Histogram  = \n']
                  log += [str(edg) + '\n']
                  log += [str(hist) +'\n']
                  sa.dist['ori']['cont'] = InitMC(
                      array([0.0, hpi, 0.0]),
                      OrientationPDFSurface,
                      [orF, ip.orpars],
                      domains=[[-pi, pi, True],[0.0, pi, False], [-pi, pi, True]],
                      idel=0.25,
                  )
            else:
              # xyz parametrisation has a uniform distribuition across all rotation axis (xyz jacobian is identity)
              # note the range of the beta angle is -hpi-hpi 
              rn1,rn2 = array([tpi,1.0,tpi]), array([pi,0.0,pi])
              sa.dist['ori'] = {}
              sa.dist['ori']['cont'] = InitICDF(1,IsotropicDistICDF,[ip.rotpar,rn1,rn2],seed=seed*932)
              if nsamp != 0: 
                vv = [ICDFsample(sa.dist['ori']['cont']) for _ in range(nsamp)]
                valpha,vbeta,vgamma = [[o[i] for o in vv] for i in range(3)]
                hist_emit(valpha, "ori_alpha", stage="initial", scope=f"molecule_m{ip.mi}")
                hist_emit(vbeta,  "ori_beta", stage="initial", scope=f"molecule_m{ip.mi}")
                hist_emit(vgamma, "ori_gamma", stage="initial", scope=f"molecule_m{ip.mi}")
                hist, edg = np.histogram(vv, bins=20)
                log += ['Orientation Histogram  = \n']
                log += [str(edg) + '\n']
                log += [str(hist) +'\n']
                sa.dist['ori']['cont'] = InitICDF(1,IsotropicDistICDF,[ip.rotpar,rn1,rn2],seed=seed*757)
        else:
            log += [" No orientational distribution " ]
        # generate velocity dist:
        if len(ip.VelPar) > 0:
            log += ["Generated molecular velocity temperature " + str(T) + "\n"]
            sa.dist['vel'] = {}
            v0, fwhm, n = ip.VelPar
            if fwhm > 0.0:
              sigma = fwhm/(2*sqrt(2.0*np.log(2.0)))
              A = 1.0/(2*sigma**2)
            else:
              sigma = 0 
              A = 1e10
            if v0 >= 0.0 and n > 0 and fwhm > 0.0: 
               #print(ip.name, '  ==', A, v0/mps2au, n, 'sigma = ', sigma/mps2au)
               cdf, v_grid = iCDFGenGauss(A,v0,n)
               sa.dist['vel']['cont'] = InitICDF(1, GenGaussSample, [cdf,v_grid],seed=seed*677)
               if nsamp != 0: 
                 vv = [ICDFsample(sa.dist['vel']['cont'])/mps2au for _ in range(nsamp)]
                 hist_emit(vv, "vel", stage="initial", scope=f"molecule_m{ip.mi}")
                 hist, edg = np.histogram(vv, bins=11)
                 log += ['Molecular Velocity Histogram  = \n']
                 log += [str(edg) + '\n']
                 log += [str(hist) +'\n']
                 sa.dist['vel']['cont'] = InitICDF(1, GenGaussSample, [cdf,v_grid],seed=seed*523) 
            else:
               if n == 0: 
                  sa.dist['vel']['cont'] = InitICDF(1, GaussianF2, [0,sigma],seed=seed*199)
                  if nsamp != 0: 
                    vv = [ICDFsample(sa.dist['vel']['cont'])/mps2au for _ in range(nsamp)]
                    hist_emit(vv, "vel", stage="initial", scope=f"molecule_m{ip.mi}")
                    sa.dist['vel']['cont'] = vv
                    hist, edg = np.histogram(vv, bins=11)
                    log += ['Molecular Velocity  Histogram  = \n']
                    log += [str(edg) + '\n']
                    log += [str(hist) +'\n']
                    sa.dist['vel']['cont'] = InitICDF(1, GaussianF2, [0,sigma],seed=seed*281)
               else: 
                  sa.dist['vel']['v'] = abs(v0)/mps2au
        return log
