#!/usr/bin/env python3
from .constants import *
from .functions import *
from .dist import *
from .mc import ICDFsample, InitICDF
from .molecules import imolecule
from . import wang
from .histograms import write_histogram_helpers_runtime
import pickle
from scipy.interpolate import CubicSpline
from joblib import Parallel, delayed
from tqdm import tqdm
import time

INFO_RULE = "=" * 67
INFO_LABEL_WIDTH = 22

def info_header(sample_id, stage):
    return [
        INFO_RULE + "\n",
        "Sample {0} | {1}\n".format(sample_id, stage),
        INFO_RULE + "\n",
    ]

def info_section(title):
    return ["\n[" + title + "]\n"]

def info_frame_marker(frame):
    return "{0:<{w}} = {1}\n".format("output frame", str(frame), w=INFO_LABEL_WIDTH)

def info_frame_transform(frame):
    if frame == "internal":
        return []
    if frame == "incoming-k-plus-z":
        return [
            "{0:<{w}} = {1}\n".format("frame transform", "Rx(pi): x,y,z -> x,-y,-z", w=INFO_LABEL_WIDTH),
            "{0:<{w}} = {1}\n".format("frame note", "sampled vectors below are reported after the output-frame transform", w=INFO_LABEL_WIDTH),
        ]
    return []

def info_scalar(label, value, unit="", fmt="{:14.7f}"):
    unit_txt = ("  " + unit) if unit else ""
    return "{0:<{w}} = {1}{2}\n".format(label, fmt.format(float(value)), unit_txt, w=INFO_LABEL_WIDTH)

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

def info_angle_vec(label, vec, comment):
    vals = ", ".join("{:10.4f}".format(float(v) / pi) for v in vec)
    return "{0:<{w}} = [ {1} ]  pi rad   # {2}\n".format(label, vals, comment, w=INFO_LABEL_WIDTH)

def info_matrix(label, mat, unit=""):
    rows = ["[ " + ", ".join("{:9.4f}".format(float(v)) for v in row) + " ]" for row in np.asarray(mat)]
    unit_txt = ("  " + unit) if unit else ""
    return "{0:<{w}} = {1}{2}\n".format(label, "; ".join(rows), unit_txt, w=INFO_LABEL_WIDTH)

class icats:
    """A class for simulating molecular scattering events."""
    class sample:   
      def __init__(self,id,scat,**dic): 
        sp = scat.sp  
        mol1 = scat.mol  
        self.slog = [] 
        if 'seed' in dic.keys():
          self.seed = dic['seed']
        else:
          self.seed = int((1+id)*np.random.random()*223)
        self.id = id
        self.mol = [mol1[0].InitializeWorker(id,seed=self.seed*17*(1+id)),mol1[1].InitializeWorker(id,seed=self.seed*11*(1+id))] 
        self.sxx = zeros(sp.shape) 
        self.svv = zeros(self.sxx.shape) 
        # Needed by pre-sampling histogram paths that call velocity sampling
        # before InitializeSample() populates per-sample state.
        self.SampInfo = {}
        self.slog = [] 
        return 

    def __init__(self, **dic):
        """Initialize the scattering simulation.

        Args:
            dic (dict): Additional parameters for initialization.
        """
        self.mol = [imolecule(), imolecule()]
        self.initsyspar() 
        self.initinppar()
        self.mol[0].ip.mi = 0
        self.mol[1].ip.mi = 1
        self.log = []
        self.slog = []
        self.sampls = {'cv': [], 'info': []}
        self.sii = 0
        self._logged_maxb_equiv = False
        return

    def initsyspar(self):
        """Initialize a system parameters.

        This method initializes a new scattering sample, setting up necessary variables and data structures.
        """
        class ssyspar: 
            def __init__(self):
             self.slog = []
             self.na = 0  
             self.nd = 0  
             self.mass = []  
             self.w0, self.w1 = 0.0, 0.0  
             self.mass2 = []
             self.mol = [] 
             self.beamang = pi/2.0
             self.chi = 0.0
             return
        self.sp = ssyspar() 

    def initinppar(self):
        """Initialize the system simulation parameters.

        Args:
            dic (dict): parameters for initialization.
        """
        class sinppar: 
            def __init__(self):
              self.nwork = 1
              self.MaxB = 0
              self.MaxJ = 0
              self.MaxL = 0
              self.FixedB = None
              self.ImpactPhi = None
              self.orbital_sampling = "geometric"
              self.velfwhm = 0.0
              self.relative_channel = None
              self.vib_mode = "sample"
              self.continues = False
              self.KeepInfo = False
              self.isotropic = True
              self.usewang = False
              self.wlmode = "default"
              self.wl_target = "auto"
              self.wl_target_user = None
              self.wl_ff = np.exp(0.10)
              self.wl_nstep_mult = 500
              self.wl_flatness = 0.90
              self.wl_wn_factor = 4.0
              self.wl_wn = None
              self.wl_j_range = None
              self.wl_l_cap = None
              self.wl_angular_sampler = "fast"
              self.wl_audit_angular_sampler = False
              self.audit_initial_sample = False
              self.audit_initial_energy_tol = 2.0e-2
              self.audit_initial_angular_tol = 0.0
              self.audit_initial_vib_tol = 0.0
              self.audit_initial_velocity_tol = 0.0
              self.wl_ff_user = None
              self.wl_nstep_user = None
              self.wl_flatness_user = None
              self.wl_wn_factor_user = None
              self.wl_wn_user = None
              self.wl_j_range_user = None
              self.wl_l_cap_user = None
              self.wl_tol = 1.000001
              self.wl_tol_user = None
              self.wl_max_iter = 0
              self.wl_log_every = 1
              self.seed_mode = "fixed"
              self.run_mode = "fresh"
              self.run_tag = None
              self.logfile_path = None
              self.progress = "normal"
              self.dry_run = False
              self.check_input = False
              self.save_frequency = 0
              self.output_format = "xyzvel"
              self.units_out = "ang-fs"
              self.output_frame = "internal"
              self.ostandard = True
              self.polarized_orientation = False
              self.plothist =False
              self.hist_initial = False
              self.hist_sampled = False
              self.hist_initial_user = False
              self.hist_sampled_user = False
              self.printout = [True,True,False,False] 
              self.pnsamp = 0
              return 
        self.ip = sinppar()
    def InitWang(self):
        ip = self.ip 
        class iWang: 
            def __init__(self): 
                if ip.wl_wn is not None:
                  self.wn = int(max(10, ip.wl_wn))
                else:
                  self.wn = int(max(10, np.round(ip.PeakJab*ip.wl_wn_factor)))
                #self.wn = int(np.round(ip.MaxJp/3))
                self.uu = np.zeros(self.wn) 
                self.hh = np.zeros(self.wn) 
                self.ff = float(ip.wl_ff)
                if self.ff <= 1.0:
                  self.ff = 1.000001
                #self.flatness = 0.85 
                self.flatness = float(ip.wl_flatness)
                self.nburn = 2
                self.nstep = int(self.wn*ip.wl_nstep_mult)
                #self.nstep = self.wn*400
        wg = iWang()
        return wg

    def InitializeWorker(self,id,**dic):  
        """  
        Initialize the sample by setting initial coordinates.  
  
        This method sets the initial coordinates for the simulation.  
        """  


        sa = self.sample(id,self)
        self.InitialDist(sa,seed=sa.seed)
        sa.sdat = self.PrepareSdat()
        return sa 


    def ReadInput(self, fnam):
        """Read input data from a file.

        Args:
            fnam (str): Input file name.
        """
        self.log.append("Reading scattering input file : " + fnam + "\n")
        ip = self.ip
        ip.inpd = File2InputList(fnam)
        ip.filename = fnam
        ip.prefix = fnam.split(".")[0]
        input_stem = os.path.splitext(os.path.basename(fnam))[0]
        run_tag = ""
        for ky, val in ip.inpd:
            if ky == "run-tag" and len(val) > 0:
                run_tag = val[0].strip()
        ip.rundir = "rd_" + (run_tag or input_stem)
        os.makedirs(ip.rundir, exist_ok=True)
        set_hist_base_dir(self._runpath("histograms"))
        self.GenerateInputData()
        if ip.hist_initial or ip.hist_sampled or ip.plothist or ip.pnsamp > 0:
            write_histogram_helpers_runtime(self._runpath("histograms"))
        self.seed = 661

    def _runpath(self, name):
        ip = self.ip
        return os.path.join(ip.rundir, name)

    def _write_costheta_convergence(self):
        """Write costheta convergence diagnostics to rd_<input>/convergence."""
        vals = self.sdat.get('orb', {}).get('cosLJab_thet', [])
        if vals is None or len(vals) == 0:
            return
        arr = np.asarray(vals, dtype=float)
        conv_dir = self._runpath("convergence")
        os.makedirs(conv_dir, exist_ok=True)

        # Final summary statistics.
        n = int(arr.size)
        mean = float(np.mean(arr))
        std = float(np.std(arr, ddof=1)) if n > 1 else 0.0
        sem = float(std / np.sqrt(n)) if n > 0 else 0.0
        mean_abs = float(np.mean(np.abs(arr)))
        summ = [
            "# cos(theta_{L,Jab}) convergence summary\n",
            f"n = {n}\n",
            f"mean = {mean:.8e}\n",
            f"std = {std:.8e}\n",
            f"sem = {sem:.8e}\n",
            f"mean_abs = {mean_abs:.8e}\n",
        ]
        with open(os.path.join(conv_dir, "costheta_summary.txt"), "w") as f:
            f.writelines(summ)

        # Cumulative trend for quick convergence inspection.
        csum = np.cumsum(arr)
        csum2 = np.cumsum(arr * arr)
        lines = ["sample_count\tmean\tstd\tsem\tmean_abs\n"]
        abs_csum = np.cumsum(np.abs(arr))
        for i in range(1, n + 1):
            mu = csum[i - 1] / i
            if i > 1:
                var = max((csum2[i - 1] - i * mu * mu) / (i - 1), 0.0)
                si = np.sqrt(var)
            else:
                si = 0.0
            sei = si / np.sqrt(i)
            mai = abs_csum[i - 1] / i
            lines.append(f"{i}\t{mu:.8e}\t{si:.8e}\t{sei:.8e}\t{mai:.8e}\n")
        with open(os.path.join(conv_dir, "costheta_cumulative.tsv"), "w") as f:
            f.writelines(lines)

    def GenerateInputData(self):
        """Generate input data from the provided parameters.

        Reads and processes input data from the input file, setting up simulation parameters and molecular properties.
        """
        mol = self.mol
        ip = self.ip
        sp = self.sp
        for ky, val in ip.inpd:
            if ky == "mol":
                self.log.append( "########################################################\n")
                self.log.append( "Reading molecular input : " + val[0] + " " + val[1] + " " + "\n")
                mi = int(val[0])
                mol[mi].ReadInput(val[1])
                self.log += mol[mi].log
                mol[mi].log = []
                if self.mol[mi].ip.ordist in ("pdf", "fixed"):
                  ip.polarized_orientation = True
                  self.log += ["Molecule-level fixed/PDF orientation detected; preserving space-fixed orientation frame\n"]
            if ky == "tvel":
                if float(val[0]) > 0.0:
                    self.log += ["Temperature for Intermolecular Velocity: "+ val[0]+ " \n"]
                    ip.Tvel = float(val[0])
                else:
                    if len(val) < 2:
                        val = [val[0], "0.0"]
                    self.log +=['Intermolecular Velocity (m/s) centre: '
                        + val[0]+ ' FWHM: ' + val[1] + '\n']
                    if len(val) > 1:
                        self.log +=['Full Width Half-Maximum (FWHM): '+ val[1] + '\n']
                        ip.velfwhm = float(val[1])*mps2au
                    ip.Tvel = float(val[0])*mps2au
            if ky == "relative-velocity":
                ip.Tvel = -abs(float(val[0])) * mps2au
                if len(val) > 1:
                    ip.velfwhm = abs(float(val[1])) * mps2au
                    self.log += ["Relative velocity FWHM: " + val[1] + " m/s\n"]
                else:
                    ip.velfwhm = 0.0
                self.log += ["Relative velocity centre: " + val[0] + " m/s\n"]
            if ky == "relative-velocity-fwhm":
                ip.velfwhm = abs(float(val[0])) * mps2au
                self.log += ["Relative velocity FWHM: " + val[0] + " m/s\n"]
            if ky == "collision-energy":
                ip.relative_channel = ("collision-energy", float(val[0]))
                ip.velfwhm = 0.0
                self.log += ["Fixed collision energy: " + val[0] + " eV\n"]
            if ky == "incoming-p0":
                ip.relative_channel = ("incoming-p0", float(val[0]))
                ip.velfwhm = 0.0
                self.log += ["Fixed incoming relative momentum p0: " + val[0] + " a.u.\n"]
            if ky == "incoming-k":
                ip.relative_channel = ("incoming-k", float(val[0]))
                ip.velfwhm = 0.0
                self.log += ["Fixed incoming wave number k: " + val[0] + " bohr^-1\n"]
            if ky == "fileout":
                self.log += ["Prefix output Prefix Name: " + val[0] + "\n"]
                ip.fileout = val[0]
            if ky == "dirout":
                raw_dir = val[0]
                if os.path.isabs(raw_dir) or raw_dir.startswith("rd_") or raw_dir.startswith("rd/"):
                    ip.dirout = raw_dir
                else:
                    ip.dirout = os.path.join(ip.rundir, raw_dir)
                self.log += ["Directory Output Name: " + ip.dirout + "\n"]
            if ky == "seed":
                self.log += ["RNG Seed: " + val[0] + "\n"]
                ip.seed = int(val[0])
                np.random.seed(ip.seed)
            if ky == "trot":
                self.log += ["Temperature for Rotational States: " + val[0] + " \n"]
                ip.Trot = float(val[0])
            if ky == "tvib":
                self.log += ["Temperature for Vibrational States: " + val[0] + " \n"]
                ip.Tvib = float(val[0])
            if ky == "maxb":
                self.log += ["Maximum impact parameter: " + val[0] + "\n"]
                ip.MaxB = float(val[0])*ang2au
            if ky == "fixed-b":
                ip.FixedB = float(val[0])*ang2au
                self.log += ["Fixed impact parameter: " + val[0] + " Ang\n"]
            if ky == "impact-phi":
                ip.ImpactPhi = float(val[0])
                self.log += ["Fixed impact-parameter azimuth phi: " + val[0] + " rad\n"]
            if ky == "output-frame":
                ip.output_frame = val[0].lower()
                if ip.output_frame not in ("internal", "incoming-k-plus-z"):
                    raise ValueError("output-frame must be internal or incoming-k-plus-z")
                self.log += ["Output/reporting frame convention: " + ip.output_frame + "\n"]
            if ky == "orbital-sampling":
                ip.orbital_sampling = val[0].lower()
                if ip.orbital_sampling not in ("geometric", "flat-l"):
                    raise ValueError("orbital-sampling must be geometric or flat-l")
                self.log += ["Orbital L sampling mode: " + ip.orbital_sampling + "\n"]
            if ky == "vib-mode":
                ip.vib_mode = val[0].lower()
                if ip.vib_mode not in ("sample", "rigid"):
                    raise ValueError("vib-mode must be sample or rigid")
                self.log += ["Vibrational mode: " + ip.vib_mode + "\n"]
            if ky == "maxj":
                self.log += ["Maximum total angular momentum J: " + val[0] + "\n"]
                ip.MaxJ = int(val[0])
            if ky == "maxl":
                self.log += ["Maximum orbital ang. momentum  L: " + val[0] + "\n"]
                ip.MaxL = int(val[0])
            if ky == "chi":
                self.log += ["Azimuthal scattering angule : " + val[0] + "\n"]
                sp.chi = float(val[0])
            if ky == "maxv":
                self.log += ["Maximum vibrational State: " + val[0] + "\n"]
                ip.MaxV = int(val[0])
            if ky == "nsamp":
                self.log += ["Number of generated samples: " + val[0] + "\n"]
                ip.Nsamp = int(val[0])
            if ky == "workers":
                self.log += ["Number of parallel workers : " + val[0] + "\n"]
                ip.nwork = int(val[0])
            if ky == "rz":
                self.log += ["Intermolecular Z-Distance: " + val[0] + "\n"]
                sp.Rz = float(val[0])*ang2au
            if ky == "beam-angle":
                self.log += ["Cross-Beam Angule        : " + val[0] + "\n"]
                sp.beamang = float(val[0])*pi/180.0
            if ky == "ordist":
                self.log += ["Orientation Distribuition Function: " + val[0] + "\n"]
                ip.ordist = val[0]
                ip.orpars = [float(v) for v in val[1:]]
            if ky == "printout":
                po  = self.ip.printout
                for i in range(4):
                  if len(val) > i:
                    po[i] = bool(int(val[i]))
                self.ip.printout = po
                if po[0]: 
                  self.log += ["Printout combined xyz/vel samples : " + str(po[0]) + "\n"]
                if po[1]: 
                  self.log += ["Printout sample info File   : " + str(po[1]) + "\n"]
                if po[2]: 
                  self.log += ["Printout directory of xyz   : " + str(po[2]) + "\n"]
                if po[3]: 
                  self.log += ["Printout info directory     : " + str(po[3]) + "\n"]
            if ky == "rot-param":
                self.log += ["Rotation angle Parametrization : " + val[0] + "\n"]
                self.log += [" (overrides molecular option)    \n"]
                ip.rotpar = "euler" if val[0] == "eul" else val[0]
                for i in range(2):
                    self.mol[i].ip.rotpar = ip.rotpar
                    if self.mol[i].ip.rotpar == 'xyz' and self.mol[i].ip.ordist == 'pdf':
                        self.log += ["for polarized distribuitions, switching to euler rotation parametrization (molecule " +str(i)+')' +'\n']
                        self.mol[i].ip.rotpar = 'euler'
            if ky == "phisample":
                self.log += ["Sample orbital azimuthal coordinate phi:" + val[0] + "\n"]
                ip.ostandard = bool(val[0]=='True')
            if ky == "plothist":
                self.log += ["Plotting Sample Histograms:" + val[0] + "\n"]
                ip.plothist = bool(val[0]=='True')
                if not ip.hist_sampled_user:
                    ip.hist_sampled = ip.plothist
            if ky == "hist_sampled":
                self.log += ["Sampled Histograms:" + val[0] + "\n"]
                ip.hist_sampled = bool(val[0]=='True')
                ip.hist_sampled_user = True
                ip.plothist = ip.hist_sampled
            if ky == "hist_initial":
                self.log += ["Initial Histograms:" + val[0] + "\n"]
                ip.hist_initial = bool(val[0]=='True')
                ip.hist_initial_user = True
            if ky == "wang":
                self.log += ["Generate WL Rejection Func:" + val[0] + "\n"]
                ip.usewang = val[0]=='True'
            if ky == "seed-mode":
                ip.seed_mode = val[0].lower()
                self.log += ["Seed mode: " + ip.seed_mode + "\n"]
            if ky == "run-mode":
                ip.run_mode = val[0].lower()
                self.log += ["Run mode: " + ip.run_mode + "\n"]
            if ky == "run-tag":
                ip.run_tag = val[0].strip()
                self.log += ["Run tag: " + ip.run_tag + "\n"]
            if ky == "wl-tol":
                ip.wl_tol_user = float(val[0])
                ip.wl_tol = ip.wl_tol_user
                self.log += ["Wang-Landau tolerance: " + val[0] + "\n"]
            if ky == "wl-max-iter":
                ip.wl_max_iter = int(val[0])
                self.log += ["Wang-Landau max iterations: " + val[0] + "\n"]
            if ky == "wl-log-every":
                ip.wl_log_every = max(1, int(val[0]))
                self.log += ["Wang-Landau log every: " + val[0] + "\n"]
            if ky == "wlmode":
                ip.wlmode = val[0].lower()
                self.log += ["Wang-Landau mode: " + ip.wlmode + "\n"]
            if ky == "wl-target":
                ip.wl_target_user = val[0].lower()
                if ip.wl_target_user not in ("linear-j", "flat-j"):
                    raise ValueError("wl-target must be linear-j or flat-j")
                self.log += ["Wang-Landau target: " + ip.wl_target_user + "\n"]
            if ky == "wl-ff":
                ip.wl_ff_user = float(val[0])
                self.log += ["Wang-Landau initial f: " + val[0] + "\n"]
            if ky == "wl-nstep":
                ip.wl_nstep_user = int(val[0])
                self.log += ["Wang-Landau nstep multiplier: " + val[0] + "\n"]
            if ky == "wl-flatness":
                ip.wl_flatness_user = float(val[0])
                self.log += ["Wang-Landau flatness: " + val[0] + "\n"]
            if ky == "wl-wn-factor":
                ip.wl_wn_factor_user = float(val[0])
                self.log += ["Wang-Landau wn factor: " + val[0] + "\n"]
            if ky == "wl-wn":
                ip.wl_wn_user = int(val[0])
                self.log += ["Wang-Landau wn bins: " + val[0] + "\n"]
            if ky == "wl-j-bins":
                ip.wl_wn_user = int(val[0])
                self.log += ["Wang-Landau J bins: " + val[0] + "\n"]
            if ky == "wl-j-range":
                ip.wl_j_range_user = float(val[0])
                self.log += ["Wang-Landau explicit J range: " + val[0] + "\n"]
            if ky == "wl-l-cap":
                ip.wl_l_cap_user = float(val[0])
                self.log += ["Wang-Landau orbital L cap: " + val[0] + "\n"]
            if ky == "wl-angular-sampler":
                ip.wl_angular_sampler = val[0].lower()
                if ip.wl_angular_sampler not in ("fast", "legacy"):
                    raise ValueError("wl-angular-sampler must be fast or legacy")
                self.log += ["Wang-Landau angular sampler: " + ip.wl_angular_sampler + "\n"]
            if ky == "wl-audit-angular-sampler":
                ip.wl_audit_angular_sampler = bool(val[0]=='True')
                self.log += ["Wang-Landau audit angular sampler: " + val[0] + "\n"]
            if ky == "audit-initial-sample":
                ip.audit_initial_sample = bool(val[0]=='True')
                self.log += ["Audit initial sample: " + val[0] + "\n"]
            if ky == "audit-initial-energy-tol":
                ip.audit_initial_energy_tol = float(val[0])
                self.log += ["Audit initial energy tolerance (eV): " + val[0] + "\n"]
            if ky == "audit-initial-angular-tol":
                ip.audit_initial_angular_tol = float(val[0])
                self.log += ["Audit initial angular tolerance: " + val[0] + "\n"]
            if ky == "audit-initial-vib-tol":
                ip.audit_initial_vib_tol = float(val[0])
                self.log += ["Audit initial vibrational coordinate tolerance: " + val[0] + "\n"]
            if ky == "audit-initial-velocity-tol":
                ip.audit_initial_velocity_tol = float(val[0])
                self.log += ["Audit initial relative-velocity tolerance (m/s): " + val[0] + "\n"]
            if ky == "keepinfo":
                self.log += ["Keeping Sample Info:" + val[0] + "\n"]
                ip.KeepInfo = bool(val[0]=='True') 
            if ky == "plotinit":
                self.log += ["Plotting Sample Histograms:" + val[0] + "\n"]
                ip.pnsamp = int(val[0])
                if not ip.hist_initial_user:
                    ip.hist_initial = ip.pnsamp > 0
            if ky == "continue":
                self.log += ["Continuation of Sampling  :" + val[0] + "\n"]
                ip.continues = bool(val[0]=='True')
            if ky == "progress":
                ip.progress = val[0].lower()
                self.log += ["Progress mode: " + ip.progress + "\n"]
            if ky == "dry-run":
                ip.dry_run = bool(val[0]=='True')
                self.log += ["Dry run: " + val[0] + "\n"]
            if ky == "check-input":
                ip.check_input = bool(val[0]=='True')
                self.log += ["Check input only: " + val[0] + "\n"]
            if ky == "save-frequency":
                ip.save_frequency = max(0, int(val[0]))
                self.log += ["Save frequency: " + val[0] + "\n"]
            if ky == "output-format":
                ip.output_format = val[0].lower()
                self.log += ["Output format: " + ip.output_format + "\n"]
            if ky == "units-out":
                ip.units_out = val[0].lower()
                self.log += ["Output units: " + ip.units_out + "\n"]

        # Keep legacy/input ergonomics: allow maxj as the primary orbital cap
        # when maxl is not explicitly provided.
        if ip.MaxL <= 0 and ip.MaxJ > 0:
          ip.MaxL = int(ip.MaxJ)
          self.log += ["Using maxj as orbital cap maxl: " + str(ip.MaxL) + "\n"]

        if ip.FixedB is not None and ip.usewang:
          raise ValueError("fixed-b is currently supported only with wang = False")

        if ip.wl_target_user is not None and not ip.usewang:
          raise ValueError("wl-target is only meaningful with wang = True")

        if not any(str(k).lower() == "orbital-sampling" for k, _ in ip.inpd):
          mode_label = "fixed-b" if ip.FixedB is not None else "geometric"
          self.log += ["Orbital L sampling mode: " + mode_label + "\n"]

        # Wang-Landau profile defaults (KIS):
        # default = current behavior, fast = cheaper, accurate = heavier.
        if ip.wlmode == "fast":
          ip.wl_ff = np.exp(0.20)
          ip.wl_nstep_mult = 300
          ip.wl_flatness = 0.85
          ip.wl_wn_factor = 3.0
          ip.wl_wn = None
        elif ip.wlmode == "accurate":
          # User requested ff near 1 and larger nstep.
          ip.wl_ff = 1.000001
          ip.wl_nstep_mult = 750
          ip.wl_flatness = 0.92
          ip.wl_tol = 1.000001
          ip.wl_wn_factor = 4.0
          ip.wl_wn = None
        else:
          ip.wl_ff = np.exp(0.10)
          ip.wl_nstep_mult = 500
          ip.wl_flatness = 0.88
          ip.wl_tol = 1.000005
          ip.wl_wn_factor = 4.0
          ip.wl_wn = None

        # Explicit per-parameter overrides:
        if ip.wl_ff_user is not None:
          ip.wl_ff = ip.wl_ff_user
        if ip.wl_nstep_user is not None:
          ip.wl_nstep_mult = ip.wl_nstep_user
        if ip.wl_flatness_user is not None:
          ip.wl_flatness = ip.wl_flatness_user
        if ip.wl_wn_factor_user is not None:
          ip.wl_wn_factor = ip.wl_wn_factor_user
        if ip.wl_wn_user is not None:
          ip.wl_wn = ip.wl_wn_user
        if ip.wl_j_range_user is not None:
          ip.wl_j_range = ip.wl_j_range_user
        if ip.wl_l_cap_user is not None:
          ip.wl_l_cap = ip.wl_l_cap_user

        if ip.usewang:
          if ip.wl_target_user is None:
            ip.wl_target = "flat-j" if ip.orbital_sampling == "flat-l" else "linear-j"
            self.log += ["Wang-Landau target: " + ip.wl_target + " (automatic)\n"]
          else:
            ip.wl_target = ip.wl_target_user
          if ip.wl_target == "flat-j" and ip.orbital_sampling != "flat-l":
            raise ValueError("wl-target = flat-j needs orbital-sampling = flat-l")
          if ip.wl_target == "linear-j" and ip.orbital_sampling != "geometric":
            raise ValueError("wl-target = linear-j needs orbital-sampling = geometric")
        if ip.wl_tol_user is not None:
          ip.wl_tol = ip.wl_tol_user

        if ip.run_mode == "continue":
          ip.continues = True
        elif ip.run_mode == "fresh" or ip.run_mode == "rebuild-wang":
          ip.continues = False

        if ip.seed_mode == "time":
          ip.seed = int(time.time())
          np.random.seed(ip.seed)
          self.log += ["Seed mode time -> RNG Seed: " + str(ip.seed) + "\n"]
        elif ip.seed_mode == "fixed" and hasattr(ip, "seed"):
          np.random.seed(ip.seed)
        # per-worker is handled by InitializeWorker/sample seeds derived from worker ids

        self.log += ["WL resolved settings: mode=" + str(ip.wlmode)
                     + " ff=" + str(ip.wl_ff)
                     + " nstep_mult=" + str(ip.wl_nstep_mult)
                     + " flatness=" + str(ip.wl_flatness)
                     + " wn_factor=" + str(ip.wl_wn_factor)
                     + " wn=" + str(ip.wl_wn)
                     + " j_range=" + str(ip.wl_j_range)
                     + " l_cap=" + str(ip.wl_l_cap) + "\n"]
 
        sp.na = mol[0].sp.na + mol[1].sp.na
        sp.el = mol[0].sp.el + mol[1].sp.el
        sp.mass = array(mol[0].sp.mass.tolist() + mol[1].sp.mass.tolist())
        sp.mass2 = np.repeat(sp.mass, 3)
        sp.rmass = (sum(mol[0].sp.mass) * sum(mol[1].sp.mass) /
                   (sum(mol[0].sp.mass) + sum(mol[1].sp.mass)) )
        sp.shape = (sp.na,3)
        self.log += ["Reduced Mass           : " + "{0:10.3e}".format(sp.rmass) + "\n"]
        sp.nd = sp.na * 3
        w0, w1 = sum(mol[0].sp.mass), sum(mol[1].sp.mass)
        sp.w0, sp.w1 = w0 / (w1 + w0), w1 / (w1 + w0)
        self.ResolveRelativeChannel()
        self.consistentT()
        if all(m.sp.na == 1 for m in mol):
            ip.Trot = 0.0
            ip.Tvib = 0.0
            self.log += ["Atom-only system detected: forcing system Trot/Tvib to 0.0\n"]
        diag_root = self._runpath("diagnostics")
        nmodes_dir = os.path.join(diag_root, "nmodes")
        ref_dir = os.path.join(diag_root, "reference")
        asym_dir = os.path.join(diag_root, "asym_rotor")
        os.makedirs(nmodes_dir, exist_ok=True)
        os.makedirs(ref_dir, exist_ok=True)
        os.makedirs(asym_dir, exist_ok=True)
        for i in range(2):
          mol[i].ip.diagdir = diag_root
          mol[i].StandardOrientat()
          if mol[i].sp.nm > 0:
            mol[i].NModesOut(os.path.join(nmodes_dir, mol[i].ip.name + "_nm.xyz"))
        self.AsymRigidRotorProjEnergies()
        sa = self.InitializeWorker(0)
        if hasattr(sp,'Rz'):
          self.SetInterZDist(sa,sp.Rz*2)
        else:
          self.SetInterZDist(sa,10.0)
        self.Mol2Image(sa)
        xyz = XYZlist(sp.el, sa.sxx * au2ang,mess="Reference Geometry, double Rz (Ang)")
        ref_path = os.path.join(ref_dir, ip.fileout.split('.')[0] + '_reference.xyz')
        open(ref_path,'w').writelines(xyz)

    # overwrite temperatures if necessary:
    def consistentT(self):
        """Ensure temperature consistency within the system.

        This method enforces consistency in temperatures by overwriting molecular temperatures with system temperatures if provided.
        """
        mol = self.mol
        for i in range(2):
            if mol[i].sp.na == 1:
                mol[i].ip.Trot = 0.0
                mol[i].ip.Tvib = 0.0
                self.log.append("Molecule " + str(i) + " is an atom: forcing molecular Trot/Tvib to 0.0\n")
                continue
            if hasattr(self.ip, "Trot"):
                self.log.append("Overwriting molecular Trot to system Trot.." + str(self.ip.Trot)+ "\n")
                mol[i].ip.Trot = self.ip.Trot
            if hasattr(self.ip, "Tvib"):
                self.log.append("Overwriting molecular Tvib to system Tvib.." + str(self.ip.Tvib)+ "\n")
                mol[i].ip.Tvib = self.ip.Tvib

    def ResolveRelativeChannel(self):
        """Convert direct relative-channel inputs to a fixed relative speed."""
        ip = self.ip
        sp = self.sp
        if ip.relative_channel is None:
            return
        kind, raw = ip.relative_channel
        if kind == "collision-energy":
            E = raw * ev2au
            if E <= 0.0:
                raise ValueError("collision-energy must be positive")
            V = sqrt(2.0 * E / sp.rmass)
            p0 = sp.rmass * V
        elif kind in ("incoming-p0", "incoming-k"):
            p0 = abs(raw)
            if p0 <= 0.0:
                raise ValueError(kind + " must be positive")
            V = p0 / sp.rmass
            E = 0.5 * p0 * p0 / sp.rmass
        else:
            raise ValueError("Unknown relative-channel input: " + str(kind))
        ip.Tvel = -V
        self.log += [
            "Resolved direct relative channel: v_rel = "
            + "{0:.7f}".format(V * au2mps)
            + " m/s, E_coll = "
            + "{0:.7f}".format(E * au2ev)
            + " eV, p0 = k = "
            + "{0:.7f}".format(p0)
            + " a.u.\n"
        ]
        return

    def AsymRigidRotorProjEnergies(self):
        """Calculate rigid rotor energies for molecules.

        This method calculates the rigid rotor energies for both molecules based on their estimated maximum rotational radii.
        """
        mol = self.mol
        for i in range(2):
            if mol[i].sp.na > 1:
                mol[i].ip.MaxR = mol[i].EstimateMaxR(self.ip.Trot)
                mol[i].AsymRigidRotorProjEnergies(mol[i].ip.MaxR)

    def saveworkers(self,wks):
        """Save distribution data to files.

        This method saves various distribution data to files, including vibrational, rotational, and orientation distributions.
        """
        ip = self.ip
        with open(self._runpath('work_sys_'+os.path.basename(ip.filename)+'.pkl'),'wb') as f:
           pickle.dump(wks,f)

    def savedata(self):
        ip = self.ip
        with open(self._runpath('dat_'+os.path.basename(ip.filename)+'.pkl'),'wb') as f:
           pickle.dump(self.sdat ,f)

    def loaddata(self):           
        ip = self.ip
        pnam = self._runpath('dat_'+os.path.basename(ip.filename)+'.pkl')
        if os.path.exists(pnam):
           with open(pnam,'rb') as f:
             self.sdat=  pickle.load(f)
           print('Read ... '+ str(len( self.sdat['vel']['ivel'])) + ' samples')

    def saveinfo(self):
        ip = self.ip
        with open(self._runpath('info_'+os.path.basename(ip.filename)+'.pkl'),'wb') as f:
           pickle.dump([self.sampls['cv'],self.sampls['info']] ,f)

    def loadinfo(self):           
        ip = self.ip
        pnam = self._runpath('info_'+os.path.basename(ip.filename)+'.pkl')
        if os.path.exists(pnam):
           with open(pnam,'rb') as f:
             self.sampls['cv'], self.sampls['info'] =   pickle.load(f)

    def loadworkers(self):
        """Load distribution data from files.

        This method loads distribution data from previously saved files, allowing for resuming simulations.
        """
        idd = str(id)
        ip = self.ip
        pnam = self._runpath('work_sys_'+os.path.basename(ip.filename)+'.pkl')
        if os.path.exists(pnam):
           with open(pnam,'rb') as f:
              wks = pickle.load(f)
        else: 
           return False
        return wks

    def _log_maxb_equivalence(self):
        """Log how maxb maps to rough L/J scales for current velocity settings."""
        ip = self.ip
        sp = self.sp
        if ip.MaxB <= 0:
            return

        def _emit(msg):
            self.log += [msg + "\n"]

        bmax_au = float(ip.MaxB)
        bmax_ang = bmax_au * au2ang
        _emit(f"maxb estimate basis: b_max = {bmax_ang:.6f} Ang ({bmax_au:.6e} a.u.)")
        _emit("Assuming JAB = 0: J ~= L ~= mu * v_rel * b_max")

        v_cases = []
        if hasattr(ip, "Tvel"):
            if ip.Tvel > 0.0:
                # Most probable relative speed for MB-like profile used in sampler.
                v_mp = np.sqrt(2.0 * kboltz * ip.Tvel / sp.rmass)
                v_cases.append(("Tvel>0 (MB most-probable v)", v_mp))
            elif ip.Tvel < 0.0:
                v0 = abs(ip.Tvel)
                if ip.velfwhm > 0.0:
                    sig = ip.velfwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
                    v_cases.append(("Tvel<0 gaussian center v", v0))
                    v_cases.append(("Tvel<0 gaussian center+sigma v", v0 + sig))
                else:
                    v_cases.append(("Tvel<0 fixed v", v0))

        if len(v_cases) == 0:
            # Fallback: derive an estimate from per-molecule beam velocity settings.
            # Molecule inputs usually provide `vel = center(m/s) fwhm(m/s) n`.
            try:
                vp0 = getattr(self.mol[0].ip, "VelPar", None)
                vp1 = getattr(self.mol[1].ip, "VelPar", None)
                if vp0 is not None and vp1 is not None and len(vp0) >= 2 and len(vp1) >= 2:
                    v10 = float(vp0[0])
                    v20 = float(vp1[0])
                    f10 = abs(float(vp0[1]))
                    f20 = abs(float(vp1[1]))
                    s10 = f10 / (2.0 * np.sqrt(2.0 * np.log(2.0)))
                    s20 = f20 / (2.0 * np.sqrt(2.0 * np.log(2.0)))

                    def _vrel_from_mol(v1, v2):
                        # Same geometry used by GetInterMolZVelocFromMolV.
                        vv2 = -v2 * z
                        vv1 = matmul(Rabout(sp.beamang, 0), -z) * v1
                        return norm(vv2 - vv1)

                    v_cases.append(("molecular-beam center v", _vrel_from_mol(v10, v20)))
                    if s10 > 0.0 or s20 > 0.0:
                        v_cases.append(("molecular-beam center+sigma v", _vrel_from_mol(v10 + s10, v20 + s20)))
            except Exception:
                v_cases = v_cases

        if len(v_cases) == 0:
            _emit(
                "maxb -> L/J note: could not infer a representative relative speed; "
                "equivalent L/J from maxb is trajectory-dependent."
            )
            return

        for label, vrel in v_cases:
            nL = sp.rmass * vrel * bmax_au
            qL = 0.5 * (-1.0 + np.sqrt(1.0 + 4.0 * nL * nL))
            qJ_max = int(np.floor(qL + 0.5))
            _emit(
                f"  [{label}] v_rel={vrel*au2mps:.3f} m/s -> "
                f"approximate maxJ~{qJ_max} (J~=L~{qL:.3f})"
            )
         

    # Generates all the necessary Nsamp distribuition samples for orientational, rotational, vibrational, velocity,
    def InitialDist(self,sa,**dic):
        """Generate initial distribution samples for different states.

        This method generates initial distribution samples for various states, including vibrational, rotational, and velocity distributions.
        """
        ip = self.ip
        sp = self.sp 
        mol = self.mol 
        # initialize sample counter and log:
        sa.slog += ["Generating " + str(ip.Nsamp) + " Samples from distribuition \n"]
        if 'nsamp' in dic.keys():
          nsamp = dic['nsamp']
        elif ip.pnsamp != 0:
          nsamp = ip.pnsamp
        else:
          nsamp = 0 
        if not ip.hist_initial:
          nsamp = 0
        if 'seed' in dic.keys():
          seed = dic['seed']
        else: 
          seed = abs((sa.id+1)*1151)
        ip.MaxJab = 0
        ip.PeakJab = 0
        # generates molecular vibrational and rotational distribuitions.
        for i in range(2):
            # overwrites molecular temperatures if system is provided:
            mol[i].ip.Tvib, mol[i].ip.Trot = ip.Tvib, ip.Trot
            mol[i].InitialDist(sa.mol[i],seed=seed*(1+i)*10000,nsamp=nsamp)
            ip.MaxJab += int(getattr(mol[i].ip, 'MaxR', 0))
            # approximate peak J only for multi-atom molecules with rotational space
            if mol[i].sp.na > 1 and hasattr(mol[i].sp, 'J2c') and mol[i].sp.J2c > 0:
                mol[i].ip.PeakJab = int(np.ceil(np.sqrt((kboltz * ip.Trot) / (2 * mol[i].sp.J2c)) - 0.5))
            else:
                mol[i].ip.PeakJab = 0
            ip.PeakJab += mol[i].ip.PeakJab
            sa.slog += mol[i].log
            mol[i].log = []
            #print("Molecule " + str(i) + " done ...")
        if sa.id == 0 and not self._logged_maxb_equiv:
            self._log_maxb_equivalence()
            self._logged_maxb_equiv = True
        sa.dist = {}
        # impact parameter / orbital angular-momentum sampling:
        if ip.FixedB is not None:
          sa.slog += ["Using fixed impact parameter b= "
                      + str(ip.FixedB * au2ang) + " Ang\n"]
          sa.dist['phi'] = {}
          sa.dist['phi']['cont'] = InitICDF(1,uniform,[0,tpi])
        elif ip.MaxL == 0:
          sa.dist['b'] = {} 
          sa.dist['b']['MaxB'] = ip.MaxB
          sa.slog += ["Generated impact parater with maximum b= " + str(ip.MaxB) + "\n"]
          sa.dist['b']['cont'] = InitICDF(1,IPICDF,[ip.MaxB],seed=seed*887)
          if nsamp != 0:
            vv = [ICDFsample(sa.dist['b']['cont']) for _ in range(nsamp)]
            sa.dist['b']['samp'] = vv
            hist_emit(vv, "b", stage="initial", scope="system")
            hist, edg = np.histogram(vv, bins=11)
            sa.slog += ['IP Histogram  = \n']
            sa.slog += [str(edg) + '\n']
            sa.slog += [str(hist) +'\n']
          sa.dist['phi'] = {}
          sa.dist['phi']['cont'] = InitICDF(1,uniform,[0,tpi])
        elif ip.MaxL > 0:
          ip.MaxLp = int(ip.MaxL*1.2)
          ip.MaxJ = ip.MaxL + int(ip.PeakJab)
          ip.MaxJp = int(ip.MaxJ*1.2)
          sa.dist['J'] = {} 
          sa.dist['J']['MaxJ'] = ip.MaxJ
          sa.slog += ["Generated total angular momentum MaxJ = " + str(ip.MaxJ) + "\n"]
          #self.dist['J']['cont'] = InitICDF(1,JcrossICDF,[self.isotropic,self.MaxJ])
          sa.dist['J']['cont'] = InitICDF(1,JcrossICDFc,[ip.isotropic,ip.MaxJp],seed=seed*31)
          sa.dist['L'] = {}
          sa.dist['L']['MaxL'] = ip.MaxL
          sa.slog += ["Generated orbital angular moment MaxL = " + str(ip.MaxL) + "\n"]
          sa.dist['L']['cont'] = InitICDF(1,JcrossICDFc,[ip.isotropic,ip.MaxLp],seed=seed*757)
          #self.dist['L']['cont'] = InitICDF(1,JcrossICDF,[self.isotropic,self.MaxL])
          if nsamp != 0:
            vv = [ICDFsample(sa.dist['J']['cont']) for _ in range(nsamp)]
            sa.dist['J']['samp'] = vv
            hist_emit(vv, "J", stage="initial", scope="system")
            hist, edg = np.histogram(vv, bins=ip.MaxJ+1)
            sa.slog += ['IP Histogram  = \n']
            sa.slog += [str(edg) + '\n']
            sa.slog += [str(hist) +'\n']
            if getattr(ip, "orbital_sampling", "geometric") == "flat-l":
              flat_l_cont = InitICDF(1, uniform, [0, 1.0], seed=seed * 73)
              vv = [float(ICDFsample(flat_l_cont)) * ip.MaxL for _ in range(nsamp)]
            else:
              vv = [ICDFsample(sa.dist['L']['cont']) for _ in range(nsamp)]
            sa.dist['L']['samp'] = vv
            hist_emit(vv, "L", stage="initial", scope="system")
            hist, edg = np.histogram(vv, bins=ip.MaxL+1)
            sa.slog += ['L Histogram  = \n']
            sa.slog += [str(edg) + '\n']
            sa.slog += [str(hist) + '\n']
          sa.dist['perchi'] = {} 
          if ip.isotropic:
              sa.dist['perchi']['cont'] = InitICDF(1,uniform,[0,tpi],seed=seed*23) 
          else:  
              raise NotImplementedError(
                  "System-level anisotropic chi sampling is not currently enabled. "
                  "Use molecule-level orientation-mode = pdf for polarized molecular orientations."
              )
        else:
         quit('MaxB or MaxJ needs to be set')
        # generate Intermolecular velocity dist:
        if hasattr(ip, "Tvel"):
            T = ip.Tvel
        else:
            T = 0
        if abs(T) > 0:
            sa.slog += ["Generated intermolecular velocity temperature " + str(T) + "\n"]
            sa.dist['vel'] = {}
            if T > 0.0: 
               A = sp.rmass / (2 * kboltz * T)
               sa.dist['vel']['cont'] = InitICDF(1,MBiCDF, [A],seed=seed*997)
               if nsamp != 0: 
                 vv = [ICDFsample(sa.dist['vel']['cont'])/mps2au for _ in range(nsamp)]
                 hist_emit(vv, "vel", stage="initial", scope="system")
                 sa.dist['vel']['cont'] = vv
                 hist, edg = np.histogram(vv, bins=11)
                 sa.slog += ['Intermolecular Velocity Histogram  = \n']
                 sa.slog += [str(edg) + '\n']
                 sa.slog += [str(hist) +'\n']
                 sa.dist['vel']['cont'] = InitICDF(1, MBiCDF, [A], seed=seed*419)
            elif T < 0.0:
               v0 = abs(T)
               if ip.velfwhm > 0.0: 
                  sigma = ip.velfwhm/(2*sqrt(2.0*np.log(2.0)))              
                  A = 1.0 / (2.0 * sigma ** 2)
                  sa.dist['vel']['cont'] = InitICDF(1, GaussianICDF, [v0,A],seed=seed*383) 
                  if nsamp != 0:                                                     
                    vv = [ICDFsample(sa.dist['vel']['cont'])/mps2au for _ in range(nsamp)]     
                    hist_emit(vv, "vel", stage="initial", scope="system")
                    hist, edg = np.histogram(vv, bins=11)                            
                    sa.slog += ['Intermolecular Velocity Histogram  = \n']                   
                    sa.slog += [str(edg) + '\n']                                         
                    sa.slog += [str(hist) +'\n']                                         
                    sa.dist['vel']['cont'] = InitICDF(1, GaussianICDF, [v0,A],seed=seed*151) 
               else:
                  sa.dist['vel']['v'] = v0
        elif nsamp != 0:
           vv = []
           for _ in range(nsamp):
             self.SampleInterMolZVeloc(sa) 
             vv.append(sa.sV/mps2au) 
           hist_emit(vv, "vel", stage="initial", scope="system")
           hist, edg = np.histogram(vv, bins=11)
           sa.slog += ['Inter-Molecular Velocity Histogram  = \n']
           sa.slog += [str(edg) + '\n']
           sa.slog += [str(hist) +'\n']   
        # finally some generic RNG for arbitrary purposes...  
        sa.dist['gen'] = {}
        sa.dist['gen']['cont'] = InitICDF(1,uniform,[0,1.0],seed=seed*73) 
        if 'phi' in sa.dist and nsamp != 0:
          vv = [ICDFsample(sa.dist['phi']['cont']) for _ in range(nsamp)]
          sa.dist['phi']['samp'] = vv
          hist_emit(vv, "phi", stage="initial", scope="system")
        if 'printlog' in dic.keys(): 
         open(ip.fileout.split(".")[0] +"_dist.log", "w").writelines(log)

    def CalcInterMolMomentum(self,sa):
        """Calculate intermolecular momentum and energy.

        This method calculates the intermolecular momentum and energy, including rotational and radial components.
        """
        
        ip = self.ip
        sp = self.sp
        mol = self.mol 
        msa = sa.mol
        msp = [mol[0].sp,mol[1].sp]
        sa.slog += info_section("intermolecular")
        sa.slog.append("{0:<{w}} = {1} x {2}\n".format("molecules", mol[0].ip.name, mol[1].ip.name, w=INFO_LABEL_WIDTH))
        #self.Mol2Image()
        com, vcom = COM(sa.sxx, sp.mass), COM(sa.svv, sp.mass)
        xx, vv = zeros((2, 3)), zeros((2, 3))
        ms = array([sum(msp[0].mass), sum(msp[1].mass)])
        for i in range(2):
            xx[i, :] = (mol[i].MolecularPosition(msa[i]) - com).T
            vv[i, :] = (mol[i].MolecularVeloc(msa[i]) - vcom).T 
        KE = np.sum(0.5*(vv.T**2*ms),axis=0).tolist()
        # difference jacobi:
        rr = xx[0, :] - xx[1, :]
        vv_rel = vv[0, :] - vv[1, :]
        # projection onto radial and angular parts
        Pr = np.outer(rr, rr) / norm(rr) ** 2
        vr = matmul(Pr, vv.T).T
        vp = vv - vr
        # angular and radial momentums:
        LL = np.sum(np.cross(xx, vv).T * ms, axis=1)
        np.set_printoptions(precision=6)
        #LL = np.sum(np.cross(xx, vp).T * ms, axis=1)
        debug = True 
        debug = False 
        if debug: 
         pp1 = (ms*vv.T).T
         vv2 = vv[0,:]-vv[1,:]
         pp2 = sp.rmass * vv2
         LL2 = np.cross(rr,pp1)[0]
         LL3 = np.cross(rr,pp2)
         if not np.allclose(LL,LL2) or not np.allclose(LL,LL3):
          print('WARNING, inconsistent COM angular momentum:')
          print('another J way1 = ', np.allclose(LL,LL2))
          print('another J way2 = ', np.allclose(LL,LL3))
        PP = vr * sp.rmass
        PP = PP[1, :] - PP[0, :]
        # angular and radial energy:
        II = sp.rmass * norm(rr) ** 2
        RotE = np.dot(LL, LL) / (2 * II)
        RadE = 0.5 * norm(PP) ** 2 / sp.rmass
        # imact parameter
        b = (norm(LL)) / (sp.rmass*norm(vv_rel))
        b0 = sum(0.5 * norm(LL) / (norm(vv, axis=1) * ms))
        phi = np.arctan2(-LL[0],LL[1])
        iJa = msa[0].SampInfo['rot']["svecJs"] 
        iJb = msa[1].SampInfo['rot']["svecJs"] 
        iJab = iJa+iJb 
        Ja = msa[0].SampInfo['rot']["svecJ0s"] 
        Jb = msa[1].SampInfo['rot']["svecJ0s"] 
        Jab = Ja+Jb 
        cJ = Jab + LL
        nJab = norm(Jab)
        qJab = 0.5 * (-1 + sqrt(1 + 4.0 * norm(Jab) ** 2))
        qJa  = 0.5 * (-1 + sqrt(1 + 4.0 * norm(Ja) ** 2) )
        qJb  = 0.5 * (-1 + sqrt(1 + 4.0 * norm(Jb) ** 2) )
        qL   = 0.5 * (-1 + sqrt(1 + 4.0 * norm(LL) ** 2) )
        qR   = 0.5 * (-1 + sqrt(1 + 4.0 * norm(PP) ** 2) )
        qJ   = 0.5 * (-1 + sqrt(1 + 4.0 * norm(cJ) ** 2) )
        niJab, niJa, niJb = norm(iJab), norm(iJa), norm(iJb)
        qiJab = 0.5 * (-1 + sqrt(1 + 4.0 * norm(iJab) ** 2))
        qiJa  = 0.5 * (-1 + sqrt(1 + 4.0 * norm(iJa) ** 2) ) 
        qiJb  = 0.5 * (-1 + sqrt(1 + 4.0 * norm(iJb) ** 2) ) 
        #print('##############')
        #print('sQL = ', qL)
        #print('oQL = ', sa.SampInfo['orb']['iL'] )
        #print('sQJ = ', qJ)
        #print('oQJ = ', sa.SampInfo['orb']['iJ'] )
        if 'orb' not in sa.SampInfo.keys():
            sa.SampInfo['orb'] = {}
        sa.SampInfo['orb']['sJab'] = Jab 
        sa.SampInfo['orb']["scL"] = LL
        sa.SampInfo['orb']["scJ"] = cJ
        sa.SampInfo['orb']['snJab'] = norm(Jab)
        sa.SampInfo['orb']["sncL"] = norm(LL)
        sa.SampInfo['orb']["sncJ"] = norm(cJ)
        sa.SampInfo['orb']['sJab'] = qJab 
        sa.SampInfo['orb']["sL"] = qL
        sa.SampInfo['orb']["sJ"] = qJ
        sa.SampInfo['orb']["sII"] = II
        sa.SampInfo['orb']["RoE"] = RotE
        sa.SampInfo['orb']["RaE"] = RadE
        sa.SampInfo['orb']["senergy"] = KE
        sa.SampInfo['orb']["sb"] = b
        sa.SampInfo['orb']["sphi"] = phi
 
        sa.slog += [info_scalar("angular energy", RotE * au2ev, "eV", "{:14.5e}")]
        sa.slog += [info_scalar("radial energy", RadE * au2ev, "eV", "{:14.5e}")]
        sa.slog += [info_scalar("total energy", (RotE + RadE) * au2ev, "eV", "{:14.5f}")]
        sa.slog += [info_vec("Ja, full", iJa, "au", "Ja", niJa, qiJa)]
        sa.slog += [info_vec("Jb, full", iJb, "au", "Jb", niJb, qiJb)]
        sa.slog += [info_vec("Jab, full", iJab, "au", "Jab", niJab, qiJab)]
        sa.slog += [info_vec("Ja, vector model", Ja, "au", "Ja", norm(Ja), qJa)]
        sa.slog += [info_vec("Jb, vector model", Jb, "au", "Jb", norm(Jb), qJb)]
        sa.slog += [info_vec("Jab, vector model", Jab, "au", "Jab", nJab, qJab)]
        sa.slog += [info_vec("L", LL, "au", "L", norm(LL), qL)]
        sa.slog += [info_vec("P_R", PP, "au", "P_R", norm(PP), qR)]
        sa.slog += [info_vec("J = L + Jab", cJ, "au", "J", norm(cJ), qJ)]
        sa.slog += [info_scalar("moment of inertia", II, "au", "{:14.3e}")]
        sa.slog += [info_scalar("b", float(b) * au2ang, "Ang", "{:14.5f}")]
        sa.slog += [info_scalar("phi", float(phi / pi), "pi rad", "{:14.5f}")]
        sa.slog += [info_vec("COM 1", mol[0].MolecularPosition(msa[0]).flatten() * au2ang, "Ang")]
        sa.slog += [info_vec("COM 2", mol[1].MolecularPosition(msa[1]).flatten() * au2ang, "Ang")]
        return 
   
         
 
    def CalcJacobiCoordinates(self,sa,**dic):
        """This takes the principal axis (the symmetry axis if its a symmetric top or near top) of each molecule and 
        treats it like a Jacobi distance (R) vector.
        We then calculate the Jacobi COM-to-COM distance and the system Euler
        angles (phi, beta, chi). The polar angle beta was historically stored
        under the key theta; that alias is retained for compatibility.
        returns information into the log file

        Parameters
        ----------
        phi, beta, chi : floats   (radians)
        """
        #copy to sxx 
        mol = self.mol
        ip = self.ip
        sp = self.sp
        self.Mol2Image(sa)
        sa.sxx -= COM(sa.sxx, sp.mass).T
        sa.svv -= COM(sa.svv, sp.mass).T
        sxx = sa.sxx.copy()
        svv = sa.svv.copy()
        mol = self.mol
        msa = sa.mol
        msp = [mol[0].sp,mol[1].sp] 
        # difference jacobi:
        r = mol[0].MolecularPosition(msa[0]).flatten() - mol[1].MolecularPosition(msa[1]).flatten()
        nr = norm(r)
        vv_rel = mol[0].MolecularVeloc(msa[0]).flatten() - mol[1].MolecularVeloc(msa[1]).flatten()
        rx, ry, rz = r
        phi   = np.arctan2(ry, rx)              # (-pi, pi]
        beta = np.arccos(rz / nr)               #   [0, pi]
        R1Z = Rabout(-phi,2)                    # rotate phi about Z 
        R2N = Rabout(-beta,1)                   # rotate beta about line of nodes
        R12 = matmul(R2N, R1Z)                  # takes r -> [0,0,|r|]
        sxx = matmul(sxx,R12.T)
        sx1, sx2 = sxx[:msp[0].na,:], sxx[msp[0].na:,:] 
        c1, c2 = COM(sx1,msp[0].mass).flatten(), COM(sx2,msp[1].mass).flatten()
        sx1, sx2 = sx1 - c1.T, sx2 - c2.T
        # Get Principal Eckart vector of molecule 0:
        vv = EckartFrameTrans(msp[0].xxe, sx1, msp[0].mass)[2,:]
        vv2 = EckartFrameTrans(msp[1].xxe, sx2, msp[1].mass)[2,:]
        vx,vy,_ = vv
        chi = np.arctan2(vy, vx)                # (-pi, pi]
        R1z = Rabout(-chi,2)                    # rotate phi about z  
        sxx = matmul(sxx,R1z.T)
        sx1 = matmul(sx1,R1z.T)
        vtest = EckartFrameTrans(msp[0].xxe, sx1, msp[0].mass)[2,:]
        vtest2 = EckartFrameTrans(msp[0].xxe, sxx[:msp[0].na,:], msp[0].mass)[2,:]
        # matrix in ZYZ which makes system to standard isotropic frame:
        RJ = matmul(Rabout(-chi,2),matmul(Rabout(-beta,1),Rabout(-phi,2)))
        Rr = matmul(Rabout(beta,1),RJ)
        rii = matmul(Rr,vv_rel)
        rii = matmul(Rr,z)
        c1, c2 = COM(sxx[:msp[0].na,:],msp[0].mass).flatten(), COM(sxx[msp[0].na:,:],msp[1].mass).flatten()
        sxx[:msp[0].na,:] = sxx[:msp[0].na,:]-c1
        sxx[msp[0].na:,:] = sxx[msp[0].na:,:]-c2
        Ra1, Ra2 = EckartFrameTrans(msp[0].xxe, sxx[:msp[0].na,:], msp[0].mass), EckartFrameTrans(msp[1].xxe, sxx[msp[0].na:,:], msp[1].mass) 
        sBFa1 = iR2q(Ra1.T)
        sBFa2 = iR2q(Ra2.T)
        if msp[0].na == 2:
          sBFa1 = [sBFa1[0], sBFa1[1], 0.0]
        if msp[1].na == 2:
          sBFa2 = [sBFa2[0], sBFa2[1], 0.0]
        sxx[:msp[0].na,:] = matmul(sxx[:msp[0].na,:],Ra1.T)+c1
        sxx[msp[0].na:,:] = matmul(sxx[msp[0].na:,:],Ra2.T)+c2
        if abs(sBFa1[1]) > 1e-6 and abs(sBFa2[1]) > 1e-6: 
          u1,u2, u3 = vv, r, vv2
          uc12, uc23 = np.cross(u1,u2), np.cross(u2,u3)
          ucc12c23, udc12c23 = np.cross(uc12,uc23), np.dot(uc12,uc23)
          dphi = np.arctan2(np.dot(u2,ucc12c23),norm(u2)*udc12c23)
        else:
          dphi = 0.0
        sa.slog += info_section("system angles")
        sa.slog += [f"{'frame note':<{INFO_LABEL_WIDTH}} = SFF is the lab Cartesian frame; Jacobi BF is the collision frame\n"]
        sa.slog += [info_angle_vec("system Euler", [phi, beta, chi], "SFF -> Jacobi BF; phi, beta, chi")]
        sa.slog += [info_angle_vec("mol 1 BF Euler", sBFa1, "molecule 1 in Jacobi BF; alpha, beta, gamma")]
        sa.slog += [info_angle_vec("mol 2 BF Euler", sBFa2, "molecule 2 in Jacobi BF; alpha, beta, gamma")]
        sa.slog += [info_scalar("v1-v2 dihedral", dphi / pi, "pi rad", "{:14.4f}")]
        sa.slog += [info_scalar("Jacobi R", float(norm(r)) * au2ang, "Ang", "{:14.5f}")]
        sa.slog += [info_matrix("SFF->Jacobi BF", RJ.T)]
        sa.slog += [info_matrix("Jacobi BF->mol1 BF", Ra1.T)]
        sa.slog += [info_matrix("Jacobi BF->mol2 BF", Ra2.T)]
        sa.SampInfo['2bJac'] = {} 
        sa.SampInfo["2bJac"]['R'] = r
        sa.SampInfo['2bJac']['dphi'] = dphi 
        sa.SampInfo['2bJac']['phi'] = phi 
        sa.SampInfo['2bJac']['beta'] = beta
        sa.SampInfo['2bJac']['theta'] = beta
        sa.SampInfo['2bJac']['chi'] = chi
        sa.SampInfo['2bJac']['Rr'] = Rr 
        sa.SampInfo['2bJac']['alpha1'],sa.SampInfo['2bJac']['beta1'], sa.SampInfo['2bJac']['gamma1'] = sBFa1
        sa.SampInfo['2bJac']['alpha2'],sa.SampInfo['2bJac']['beta2'], sa.SampInfo['2bJac']['gamma2'] = sBFa2
        return  

    

    def InitializeSample(self,sa,ii):
        """Initialize a scattering sample.

        This method initializes a new scattering sample, setting up necessary variables and data structures.
        """
        sp = self.sp
        sa.sii = ii 
        sa.slog = info_header(ii, "generation")
        sa.slog += [info_frame_marker(self.ip.output_frame)]
        sa.slog += info_frame_transform(self.ip.output_frame)
        sa._output_frame_applied = False
        sa.SampInfo = {}
        sa.svv = zeros(sp.shape)  
        sa.sxx = zeros(sp.shape)
        mol = self.mol
        msa = sa.mol
        sa.smkin = [0.0,0.0]
        for i in range(2):
            mol[i].InitializeSample(msa[i])

    def GenerateSample(self,sa,rang):
        """Generate a scattering sample.

        This method generates a complete scattering sample, including vibrational, rotational, and orientational states, as well as intermolecular parameters.
        """
        sp = self.sp
        ip = self.ip
        wn = len(sp.td) 
        debug = False
        #if sa.id ==0:
         #debug = True
        xyz = []
        slog = []
        sa.ixyz = []
        sa.ivxyz = []
        if ip.progress == "quiet":
          hidebar = True
        elif sa.id == 0:
          hidebar = False
        else:
          hidebar = True
        ff = 1.0
        ntrial = 0
        naccepted = 0
        nreject_jcap = 0
        nreject_wl = 0
        nreject_wl_range = 0
        live_log_path = getattr(ip, "logfile_path", None)
        def append_live_sampling_log(force=False):
          if live_log_path is None:
            return
          if not force and naccepted % 50 != 0:
            return
          acc_rate = float(naccepted) / max(1.0, float(ntrial))
          line = (
              "Sampling worker {wid}: accepted={acc} trials={tr} "
              "trial_acceptance={ar:.5f} J_cap_rejects={jc} "
              "WL_rejects={wr} WL_range_rejects={rr}\n"
          ).format(
              wid=sa.id,
              acc=naccepted,
              tr=ntrial,
              ar=acc_rate,
              jc=nreject_jcap,
              wr=nreject_wl,
              rr=nreject_wl_range,
          )
          try:
            with open(live_log_path, "a") as lf:
              lf.write(line)
          except OSError:
            pass
        pbar = tqdm(range(rang[0],rang[1]), disable=hidebar, dynamic_ncols=True)
        for ii in pbar:
          log = []
          log.append("Generating Sample number " + str(ii + 1) + "\n")
          # Initialize sample
          debug and print('Initialize...',ii)  
          while True:
            ntrial += 1
            self.InitializeSample(sa,ii)
            log = []
            if ip.FixedB is not None:
              log += self.SampleRigidRotorState0(sa)
              log += self.SampleOrientat0(sa)
              self.SampleInterMolZVeloc(sa)
              log += self.SampleFixedImpactOrbital(sa)
              log += self.SampleJ(sa)
              break
            else:
              log  = self.SampleOrbitalL(sa,cap=ip.MaxL*ff)
              log += self.SampleRigidRotorState0(sa)
              log += self.SampleOrientat0(sa)
              log += self.SampleJ(sa)
              # Old-style gate: cap very large J before WL acceptance.
              if ip.MaxL > 0 and sa.sJ >= ip.MaxL*ff:
                nreject_jcap += 1
                continue
              if len(sp.td) == 0:
                break
              iiL = int(np.floor(wn*abs(sa.sJ)/ip.MaxJ))
              if iiL < len(sp.td):
                if np.random.rand() <= sp.td[iiL]:
                  break
                nreject_wl += 1
              else:
                nreject_wl_range += 1
                continue
          naccepted += 1
          if not hidebar and ip.progress == "normal" and (naccepted == 1 or naccepted % 50 == 0):
            acc_rate = float(naccepted) / max(1.0, float(ntrial))
            pbar.set_postfix_str(
                f"trial_acc={acc_rate:.3f} "
                f"jcap={nreject_jcap} wl={nreject_wl} wl_range={nreject_wl_range}"
            )
          append_live_sampling_log()
          sa.slog += log
          if ip.isotropic and ip.ostandard and ip.ImpactPhi is None and not ip.polarized_orientation:
             self.SetStandardOrientation(sa)
          debug and print('Sample HOVib... ')  
          if ip.vib_mode == "rigid":
            sa.slog += info_section("vibration")
            sa.slog += [f"{'vib-mode':<{INFO_LABEL_WIDTH}} = rigid; harmonic oscillator sampling skipped\n"]
          else:
            self.SampleHOVibrState(sa)
          # sample intermolecular DOF
          debug and print('Sample InterMolZ...')
          self.SetInterZDist(sa,sp.Rz) # set z distance ...
          if ip.FixedB is None:
            self.SampleInterMolZVeloc(sa) # need to get the magnitude of intermol v before getting impact parameter
          debug and print('Sample Impact Param...')  
          if not hasattr(sa, "sb"):
            sa.sb = sa.snL/(sp.rmass*sa.sV)
            sa.sphi = np.arctan2(-sa.scL[0],sa.scL[1])
          self.SetImpactParam(sa,sa.sb, sa.sphi)
          self.SetInterMolZVeloc(sa)
          self.ApplyOutputFrameConvention(sa)
          self.Mol2Image(sa)
          self.StoreOrbitalInfoLog(sa)
          debug and print('Calc Jacobi ...')  
          self.CalcJacobiCoordinates(sa)
          # summarize energy from generated sample
          self.SummarizeLogEnergy(sa,False)
          self.CaptureInitialAuditState(sa, "generation")
          # summarize energy from calculated/analysed sample:
          self.AnalyseSample(sa)
          self.CaptureInitialAuditState(sa, "analysis")
          self.AuditInitialSample(sa)

          ixyz, ivxyz = self.ListXYZOut(sa,mess="Sample " + str(sa.sii))
          if ip.printout[0]:
            sa.ixyz += ixyz
            sa.ivxyz += ivxyz
          if ip.printout[1]:
            slog += sa.slog
          if ip.printout[2]:
            self.ImageXYZOut(ixyz,ivxyz,sa.sii,sa=sa)
          if ip.printout[3]:
            open(ip.dirout + "/" + ip.fileout +"_" + str(sa.sii) + ".info", "w").writelines(sa.slog)
          if sa.id == 0:
            self.AddInfoToSamples(sa.sdat,sa.SampInfo)
          else:
            self.AddInfoToSamples(sa.sdat,sa.SampInfo)
          if (ii+1)%1000 == 0 and sa.id == 0 and ip.progress == "verbose":
            print('AVERAGE COSTHETA = ', np.mean(sa.sdat['orb']['cosLJab_thet']))
          if ip.save_frequency > 0 and ((ii+1) % ip.save_frequency == 0):
            if not os.path.isdir(ip.dirout):
              os.system("mkdir " + ip.dirout)
            cnam = ip.dirout + "/" + ip.fileout + "_checkpoint_w" + str(sa.id) + "_n" + str(ii+1) + ".pkl"
            with open(cnam, "wb") as cf:
              pickle.dump(sa.sdat, cf)
        if ntrial > 0:
          acc_rate = float(naccepted) / float(ntrial)
          slog += [
              "\n[sampling diagnostics]\n",
              info_scalar("accepted samples", naccepted, "", "{:14.0f}"),
              info_scalar("total trials", ntrial, "", "{:14.0f}"),
              info_scalar("trial acceptance", acc_rate, "", "{:14.5f}"),
              info_scalar("J-cap rejects", nreject_jcap, "", "{:14.0f}"),
              info_scalar("WL rejects", nreject_wl, "", "{:14.0f}"),
              info_scalar("WL range rejects", nreject_wl_range, "", "{:14.0f}"),
          ]
          append_live_sampling_log(force=True)
        sa.slog = slog 
        return sa

    def SampleJ(self,sa):
        msa = sa.mol
        Jab = msa[0].srpar[-1] + msa[1].srpar[-1]
        cJ = sa.scL + Jab 
        sa.scJ =  cJ
        nJ = norm(cJ) 
        sa.snJ =  nJ
        sa.sJ = 0.5 * (-1 + sqrt(1 + 4.0 * norm(nJ) ** 2) )
        log = [info_scalar("total J", sa.sJ, "", "{:14.2f}")]
        sa.SampInfo['orb']['iJ'] = sa.sJ
        sa.SampInfo['orb']['inJ'] = sa.snJ
        sa.SampInfo['orb']['icJ'] = sa.scJ
        return log 

    def GenerateWang(self,wks):
        def SampleWangAngularMomenta(self, sa, ii, cap, fast=True, audit=False):
          if audit:
            import copy
            import random
            rng_state = np.random.get_state()
            py_rng_state = random.getstate()
            old_sa = copy.deepcopy(sa)
            new_sa = copy.deepcopy(sa)
            old_vals = SampleWangAngularMomenta(self, old_sa, ii, cap, fast=False, audit=False)
            np.random.set_state(rng_state)
            random.setstate(py_rng_state)
            new_vals = SampleWangAngularMomenta(self, new_sa, ii, cap, fast=True, audit=False)
            dJ = abs(old_vals[0] - new_vals[0])
            dL = abs(old_vals[1] - new_vals[1])
            dJab = norm(old_vals[2] - new_vals[2])
            if dJ > 1.0e-9 or dL > 1.0e-9 or dJab > 1.0e-9:
              raise ValueError(
                  "Fast Wang-Landau angular sampler mismatch: "
                  + "dJ=" + str(dJ) + " dL=" + str(dL) + " dJab=" + str(dJab)
              )
            np.random.set_state(rng_state)
            random.setstate(py_rng_state)

          self.InitializeSample(sa, ii)
          if not fast:
            self.SampleOrbitalL(sa,cap=cap)
            self.SampleRigidRotorState0(sa)
            self.SampleOrientat0(sa)
            self.SampleJ(sa)
            return sa.sJ, sa.sL, sa.sjab

          self.SampleOrbitalL(sa,cap=cap)
          Jab = np.zeros(3)
          for imol in range(2):
            moli = self.mol[imol]
            msai = sa.mol[imol]
            if moli.sp.na > 1 and 'rotJ' in msai.dist:
              moli.SampleTotMolAngMom(msai)
              moli.SampleRigidRotorState(msai)
              if 'ori' in msai.dist:
                moli.SampleRotation(msai)
                Ji = matmul(msai.soR, msai.srpar[-1])
              else:
                Ji = msai.srpar[-1]
            else:
              Ji = np.zeros(3)
            Jab += Ji
          sa.sjab = Jab
          sa.snjab = norm(Jab)
          sa.scJ = sa.scL + Jab
          sa.snJ = norm(sa.scJ)
          sa.sJ = 0.5 * (-1 + sqrt(1 + 4.0 * norm(sa.snJ) ** 2) )
          return sa.sJ, sa.sL, Jab

        def GetiL(self,sa,rnge,cap):
          cost = []
          iiLL = np.empty(rnge[1] - rnge[0], dtype=np.int32)
          for k, s in enumerate(range(rnge[0],rnge[1])):
            SampleWangAngularMomenta(
                self, sa, s, cap,
                fast=ip.wl_angular_sampler == "fast",
                audit=ip.wl_audit_angular_sampler and n == 0 and k < 20,
            )
            if sa.snL*sa.snjab > 0.0: 
              costhet = np.dot(sa.scL,sa.sjab)/(sa.snL*sa.snjab)
            else:
              costhet = -2
            sa.costhet = costhet
            if abs(sa.costhet) < 1.0:
             cost.append(sa.costhet)
            iiL = int(np.floor(wg.wn*abs(sa.sJ)/wg.maxr))
            iiLL[k] = iiL
          return iiLL
        def flatten(xss):
          return [x for xs in xss for x in xs] 
        wg = self.InitWang()
        n = 0
        iL = 4
        ip = self.ip 
        sp = self.sp 
        ff_init = float(wg.ff)
        ff_tol = float(ip.wl_tol)
        log_denom = np.log(max(ff_init, ff_tol + 1e-15)) - np.log(max(ff_tol, 1e-15))
        if abs(log_denom) < 1.0e-15:
          log_denom = 1.0
        wg.maxr = float(ip.wl_j_range) if ip.wl_j_range is not None else float(ip.PeakJab*4)
        cap = float(ip.wl_l_cap) if ip.wl_l_cap is not None else float(ip.PeakJab*5)
        wl_bar = None
        wl_bar_total = 1000
        if ip.progress == "normal":
          wl_bar = tqdm(
              total=wl_bar_total,
              desc="WL",
              unit="‰",
              dynamic_ncols=True,
              bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}] {postfix}",
          )
        while wg.ff > ip.wl_tol: 
          if ip.wl_max_iter > 0 and n >= ip.wl_max_iter:
            break
          ac = 0
          wg.hh = np.zeros(wg.wn)
          edges = np.linspace(0, wg.nstep, ip.nwork + 1, dtype=int)
          iiLL_parts = Parallel(n_jobs=ip.nwork)(delayed(GetiL)(self, wks[i], [int(edges[i]), int(edges[i+1])], cap) for i in range(ip.nwork))
          iiLL = np.concatenate(iiLL_parts)
          for iiL in iiLL: 
            if iiL <= wg.wn-1:
              accept_prob = min(1,np.exp(wg.uu[iL] - wg.uu[iiL]))         
              if np.random.rand() < accept_prob:
                  iL = iiL
                  ac += 1
            else:
              continue
            wg.uu[iL] += np.log(wg.ff)
            wg.hh[iL] += 1
          hmin = float(wg.hh.min())
          hcrit = float(wg.flatness * wg.hh.mean())
          acc_rate = float(ac) / float(wg.nstep)
          flat_ratio = hmin / (hcrit + 1.0e-15)
          log_gap = np.log(max(wg.ff, ff_tol + 1e-15)) - np.log(max(ff_tol, 1e-15))
          log_rem = log_gap / log_denom
          log_rem = float(np.clip(log_rem, 0.0, 1.0))
          prog = 1.0 - log_rem
          prog = float(np.clip(prog, 0.0, 1.0))
          if ip.progress != "quiet" and (n % ip.wl_log_every == 0):
            if ip.progress == "verbose":
              print(n, 'Hmin = ', hmin, ' Criterion: ', hcrit, 'ac = ', acc_rate)
              print(f'   WL status: flat={flat_ratio:.4f} (>=1.0 target), ff={wg.ff:.8f}, progress={100.0*prog:.2f}%, log_gap={log_gap:.3e}, rem_log={100.0*log_rem:.2f}%')
            elif wl_bar is not None:
              new_n = int(round(prog * wl_bar_total))
              wl_bar.update(max(0, new_n - wl_bar.n))
              wl_bar.set_postfix_str(
                  f"iter={n} flat={flat_ratio:.4f} ff={wg.ff:.8f} "
                  f"acc={acc_rate:.3f} log_gap={log_gap:.3e}"
              )
          if hmin > hcrit and n > wg.nburn: 
            wg.ff = np.sqrt(wg.ff) 
            if ip.progress != "quiet":
              log_gap2 = np.log(max(wg.ff, ff_tol + 1e-15)) - np.log(max(ff_tol, 1e-15))
              log_rem2 = log_gap2 / log_denom
              log_rem2 = float(np.clip(log_rem2, 0.0, 1.0))
              prog2 = 1.0 - log_rem2
              prog2 = float(np.clip(prog2, 0.0, 1.0))
              msg = f'WL update: flatness reached, reducing ff -> {wg.ff:.8f} ({100.0*prog2:.2f}% to tol, rem_log={100.0*log_rem2:.2f}%)'
              if wl_bar is not None:
                tqdm.write(msg)
              else:
                print('   ' + msg)
            if ip.progress == "verbose":
              print('UU = ') 
              print(wg.uu - wg.uu.min()) 
          elif n <= wg.nburn and n%2 == 0:
            wg.hh[:] = 0. 
          if n%10 == 0: 
            wg.hh[:] = 0. 
          wg.uu -= wg.uu.min()
          n += 1 
        if ip.progress != "quiet":
          log_gapf = np.log(max(wg.ff, ff_tol + 1e-15)) - np.log(max(ff_tol, 1e-15))
          log_remf = log_gapf / log_denom
          log_remf = float(np.clip(log_remf, 0.0, 1.0))
          progf = 1.0 - log_remf
          progf = float(np.clip(progf, 0.0, 1.0))
          if wl_bar is not None:
            wl_bar.update(max(0, wl_bar_total - wl_bar.n))
            wl_bar.set_postfix_str(
                f"iter={n} ff={wg.ff:.8f} tol={ip.wl_tol:.8f} "
                f"log_gap={log_gapf:.3e}"
            )
            wl_bar.close()
          print(f'WL complete: iter={n}, ff={wg.ff:.8f}, tol={ip.wl_tol:.8f}, progress={100.0*progf:.2f}%, log_gap={log_gapf:.3e}, rem_log={100.0*log_remf:.2f}%')
        if ip.progress == "verbose":
          print('U = ', wg.uu)
        dlj = float(wg.maxr)/float(wg.wn)
        if ip.progress != "quiet":
          print(f'Convergence reached: WL profile built (wn={wg.wn}, dlj={dlj:.4f}, maxr={wg.maxr})')
        jj = dlj*0.5+ np.linspace(0,wg.maxr-dlj*1,wg.wn)
        uu_lin = np.exp(wg.uu - wg.uu.min()) 
        uu_lin = uu_lin / 0.5
        if ip.progress == "verbose":
          print('MaxJP = ', ip.MaxJp) 
        if False:
         mxj  = wg.uu.tolist().index(wg.uu.max())
         print('MXJ = ', mxj+1) 
         self.ip.MaxJ = int(jj[mxj+1])
         print('Updated MaxJ = ', self.ip.MaxJ) 
        x0, y0 = jj[0], uu_lin[0]  
        x1, y1 = jj[1], uu_lin[1] 
        if ip.progress == "verbose":
          print('x0 = ', x0, ' y0 =', y0)
          print('x1 = ', x1, ' y1 =', y1)
        # some suitable positive value 
        u0 =  (y1 * x0**2 - y0 * x1**2) / ( x0**2 - x1**2 ) 
        u0 = u0*0.9 + y0*0.1
        jj = np.array([-jj[1],-jj[0],0.0] + jj.tolist())
        uu = uu_lin.tolist() 
        uu = np.array([uu[1],uu[0], u0] + uu  )

        if ip.progress == "verbose":
          print('jlens = ', len(jj))
          print('ulens = ', len(uu))
          print('jj = ', jj)
          print('uu = ', uu) 
        wg.iwl = CubicSpline(jj, uu, bc_type="natural", extrapolate=True)
        if ip.progress == "verbose":
          print('iwl= ', [wg.iwl(j) for j in range(ip.PeakJab*4)])
        #print('iwl= ', [wg.iwl(j) for j in range(ip.MaxJp)])
        nn = min([wg.iwl(float(j)) for j in range(0,3)])
        if ip.progress == "verbose":
          print('nn = ', nn)
        wg.iwl = CubicSpline(jj, uu/nn, bc_type="natural", extrapolate=True)
        if False:
          sp.iwld = np.array([wg.iwl(j) for j in range(ip.MaxJ)])
          if ip.progress == "verbose":
            print('IWLD = ', sp.iwld)
          sp.td = np.array([1+2*J for J in range(ip.MaxJ)])/sp.iwld 
          sp.td = sp.td/sp.td.max()
        else:
          wl_range = int(np.ceil(wg.maxr))
          sp.iwld = np.array([wg.iwl(j) for j in range(wl_range)])
          if ip.progress == "verbose":
            print('IWLD = ', sp.iwld)
          if ip.wl_target == "linear-j":
            target_j = np.array([1+2*J for J in range(wl_range)], dtype=float)
          elif ip.wl_target == "flat-j":
            target_j = np.ones(wl_range, dtype=float)
          else:
            raise ValueError("Unknown wl-target: " + str(ip.wl_target))
          sp.td = target_j/sp.iwld
          tail0 = min(len(sp.td) - 1, max(0, int(len(sp.td)*0.75)))
          mntd = np.mean(sp.td[tail0:])
          td2 = [mntd for _ in range(wl_range,ip.MaxJ)]
          sp.td = np.array(sp.td.tolist() +td2)
          sp.td = sp.td/sp.td.max()

        if ip.progress == "verbose":
          print('TD = ', sp.td)
        uu_lin = np.exp(wg.uu - wg.uu.min()) 
        sp.uu = uu_lin 
        wang.save(self._runpath('wang.pkl'), sp.uu, sp.iwld, sp.td, wang.metadata_from_input(ip))
        wl_dir = self._runpath("histograms/wl")
        os.makedirs(wl_dir, exist_ok=True)
        write_wl_plot_script(
            sp.td,
            J_range=(0, ip.MaxJ),
            script_path=os.path.join(wl_dir, "wl_td_plot.py"),
            outfile="wl_td_plot.png",
            title="WL sampling correction vs. J",
            ylabel="target weight (2J+1)/Omega(J)",
        )
        write_wl_plot_script(
            sp.iwld,
            J_range=(0, len(sp.iwld)),
            script_path=os.path.join(wl_dir, "wl_wl_plot.py"),
            outfile="wl_wl_plot.png",
            title="Estimated sampled J density",
            ylabel="Omega(J), normalized",
        )





    def SummarizeLogEnergy(self,sa,FromSample):
        """Summarize energy-related information.

        This method summarizes energy-related information, including vibrational, rotational, and velocity contributions.
        """
        ven1, ven2, ren1, ren2, ken1, ken2 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        msa = sa.mol
        if FromSample:
            if 'senergy' in msa[0].SampInfo.get('vib', {}).keys():
                ven1 = sum(msa[0].SampInfo['vib']['senergy']) * au2ev 
            if 'senergy' in msa[1].SampInfo.get('vib', {}).keys():
                ven2 = sum(msa[1].SampInfo['vib']['senergy']) * au2ev 
            if 'senergy' in msa[0].SampInfo.get('rot', {}).keys():
                ren1 = sum(msa[0].SampInfo['rot']['senergy']) * au2ev 
            if 'senergy' in msa[1].SampInfo.get('rot', {}).keys():
                ren2 = sum(msa[1].SampInfo['rot']['senergy']) * au2ev 
            if 'senergy' in sa.SampInfo['orb'].keys():
                ken1 = sa.SampInfo['orb']['senergy'][0] * au2ev
                ken2 = sa.SampInfo['orb']['senergy'][1] * au2ev
            stage = "analysis"
        else:
            if hasattr(msa[0], "sven"):
                ven1 = sum(msa[0].sven) * au2ev
            if hasattr(msa[1], "sven"):
                ven2 = sum(msa[1].sven) * au2ev
            if hasattr(msa[0], "sren"):
                ren1 = msa[0].scren * au2ev
            if hasattr(msa[1], "sren"):
                ren2 = msa[1].scren * au2ev
            if hasattr(sa, "smkin"):
                ken1 = sa.smkin[0] * au2ev
                ken2 = sa.smkin[1] * au2ev
            stage = "generation"
        sa.slog += info_section("energy summary")
        sa.slog.append(
            "{0:<{w}} = {1:>16s} {2:>14s} {3:>14s}\n".format(
                "component",
                self.mol[0].ip.name + "/eV",
                self.mol[1].ip.name + "/eV",
                "total/eV",
                w=INFO_LABEL_WIDTH,
            )
        )
        rows = [
            ("vibrational", ven1, ven2, ven1 + ven2),
            ("rotational", ren1, ren2, ren1 + ren2),
            ("velocity", ken1, ken2, ken1 + ken2),
            ("total", ken1 + ven1 + ren1, ken2 + ven2 + ren2, ken1 + ken2 + ven1 + ven2 + ren1 + ren2),
        ]
        for name, v1, v2, vt in rows:
            sa.slog.append("{0:<{w}} = {1:16.4f} {2:14.4f} {3:14.4f}\n".format(name, v1, v2, vt, w=INFO_LABEL_WIDTH))
        summary = {
            "vib": ven1 + ven2,
            "rot": ren1 + ren2,
            "vel": ken1 + ken2,
            "total": ken1 + ken2 + ven1 + ven2 + ren1 + ren2,
        }
        if FromSample:
            sa.audit_energy_analysis = summary
        else:
            sa.audit_energy_generation = summary
        return summary

    def CaptureInitialAuditState(self, sa, stage):
        """Store scalar generation/analysis quantities for the optional t=0 audit."""
        orb = sa.SampInfo.get("orb", {})
        def _arr(value):
            if value is None:
                return None
            return np.asarray(value, dtype=float).copy()

        def _scalar(value):
            if value is None:
                return None
            try:
                return float(value)
            except (TypeError, ValueError):
                return None

        def _relative_velocity_mps():
            try:
                v0 = self.mol[0].MolecularVeloc(sa.mol[0]).flatten()
                v1 = self.mol[1].MolecularVeloc(sa.mol[1]).flatten()
                return float(norm(v0 - v1) * au2mps)
            except Exception:
                return None

        vectors = {}
        scalars = {}
        angles = {}
        molecule_vectors = {}
        molecule_scalars = {}
        vibrational = {}

        jac = sa.SampInfo.get("2bJac", {})
        if "R" in jac:
            scalars["Jacobi_R"] = float(norm(jac["R"]))
        for key in ("phi", "beta", "theta", "chi", "dphi"):
            if key in jac:
                angles["Jacobi_" + key] = float(jac[key])

        if stage == "generation":
            angular = {
                "L": orb.get("ncL"),
                "J": orb.get("ncJ"),
                "Jab": orb.get("nJab"),
            }
            vectors["L"] = _arr(orb.get("cL"))
            vectors["J"] = _arr(orb.get("cJ"))
            vectors["Jab"] = _arr(orb.get("Jab"))
            scalars["b"] = _scalar(orb.get("b"))
            if self.ip.ImpactPhi is not None:
                angles["impact_phi"] = _scalar(orb.get("sampled_impact_phi", orb.get("phi")))
            else:
                angles["impact_phi"] = _scalar(orb.get("phi"))
            scalars["relative_velocity_mps"] = _scalar(sa.SampInfo.get("vel", {}).get("ivel"))
            energy = getattr(sa, "audit_energy_generation", {})
        else:
            angular = {
                "L": orb.get("sncL"),
                "J": orb.get("sncJ"),
                "Jab": orb.get("snJab"),
            }
            vectors["L"] = _arr(orb.get("scL"))
            vectors["J"] = _arr(orb.get("scJ"))
            if vectors["L"] is not None and vectors["J"] is not None:
                vectors["Jab"] = vectors["J"] - vectors["L"]
            scalars["b"] = _scalar(orb.get("sb"))
            angles["impact_phi"] = _scalar(orb.get("sphi"))
            scalars["relative_velocity_mps"] = _relative_velocity_mps()
            energy = dict(getattr(sa, "audit_energy_analysis", {}))
            vec_rot = []
            for msa in sa.mol:
                rot = msa.SampInfo.get("rot", {})
                if "senergy_vec" in rot:
                    vec_rot.append(sum(rot["senergy_vec"]) * au2ev)
                elif "senergy" in rot:
                    vec_rot.append(sum(rot["senergy"]) * au2ev)
                else:
                    vec_rot.append(0.0)
            if energy:
                energy["rot"] = sum(vec_rot)
                energy["total"] = energy.get("vib", 0.0) + energy["rot"] + energy.get("vel", 0.0)

        for idx, msa in enumerate(sa.mol):
            label = "m" + str(idx)
            rot = msa.SampInfo.get("rot", {})
            vib = msa.SampInfo.get("vib", {})
            if stage == "generation":
                molecule_vectors[label + "_J_vector"] = _arr(msa.srpar[-1])
                molecule_scalars[label + "_J"] = _scalar(norm(msa.srpar[-1]))
            else:
                molecule_vectors[label + "_J_vector"] = _arr(rot.get("svecJ0s"))
                if rot.get("svecJ0s") is not None:
                    molecule_scalars[label + "_J"] = _scalar(norm(rot.get("svecJ0s")))
            if stage == "generation":
                if "Q" in vib:
                    vibrational[label + "_Q"] = _arr(vib.get("Q"))
                if "P" in vib:
                    vibrational[label + "_P"] = _arr(vib.get("P"))
            else:
                if "sQ" in vib:
                    vibrational[label + "_Q"] = _arr(vib.get("sQ"))
                if "sP" in vib:
                    vibrational[label + "_P"] = _arr(vib.get("sP"))

        setattr(
            sa,
            "audit_initial_" + stage,
            {
                "energy": energy,
                "angular": angular,
                "vectors": vectors,
                "scalars": scalars,
                "angles": angles,
                "molecule_vectors": molecule_vectors,
                "molecule_scalars": molecule_scalars,
                "vibrational": vibrational,
            },
        )

    def AuditInitialSample(self, sa):
        """Check that generated sample bookkeeping survives analysis of its coordinates."""
        ip = self.ip
        if not ip.audit_initial_sample:
            return
        gen = getattr(sa, "audit_initial_generation", None)
        ana = getattr(sa, "audit_initial_analysis", None)
        if gen is None or ana is None:
            raise ValueError("Initial sample audit requested, but audit state is incomplete.")

        failures = []
        def _angle_diff(a, b):
            return abs(float(np.arctan2(np.sin(float(a) - float(b)), np.cos(float(a) - float(b)))))

        def _vec_diff(a, b):
            aa = np.asarray(a, dtype=float)
            bb = np.asarray(b, dtype=float)
            if aa.shape != bb.shape:
                return None
            return float(norm(aa - bb))

        def _audit_scalar(group_name, key, gv, av, tol, unit="", enforce=True):
            if gv is None or av is None:
                return
            diff = abs(float(gv) - float(av))
            status = "OK" if diff <= tol else ("FAIL" if enforce else "INFO")
            if status == "FAIL":
                failures.append(f"{group_name} {key}: generation={float(gv):.8g}, analysis={float(av):.8g}, diff={diff:.3g}")
            unit_txt = (" " + unit) if unit else ""
            sa.slog.append(
                "Audit {0:<10s} {1:>18s}: generation {2:13.6g}, analysis {3:13.6g}, diff {4:10.3e}{5} [{6}]\n".format(
                    group_name, key, float(gv), float(av), diff, unit_txt, status
                )
            )

        def _audit_vector(group_name, key, gv, av, tol, unit="", enforce=True):
            if gv is None or av is None:
                return
            diff = _vec_diff(gv, av)
            if diff is None:
                if enforce:
                    failures.append(f"{group_name} {key}: vector shapes differ")
                status = "FAIL" if enforce else "INFO"
                sa.slog.append("Audit {0:<10s} {1:>18s}: vector shapes differ [{2}]\n".format(group_name, key, status))
                return
            status = "OK" if diff <= tol else ("FAIL" if enforce else "INFO")
            if status == "FAIL":
                failures.append(f"{group_name} {key}: vector norm diff={diff:.3g}")
            unit_txt = (" " + unit) if unit else ""
            sa.slog.append(
                "Audit {0:<10s} {1:>18s}: vector-norm diff {2:10.3e}{3} [{4}]\n".format(
                    group_name, key, diff, unit_txt, status
                )
            )

        def _audit_component_vector(group_name, key, gv, av, tol, unit="", enforce=True):
            if gv is None or av is None:
                return
            aa = np.asarray(gv, dtype=float)
            bb = np.asarray(av, dtype=float)
            if aa.shape != bb.shape:
                if enforce:
                    failures.append(f"{group_name} {key}: vector shapes differ")
                status = "FAIL" if enforce else "INFO"
                sa.slog.append("Audit {0:<10s} {1:>18s}: vector shapes differ [{2}]\n".format(group_name, key, status))
                return
            if aa.size == 0:
                return
            delta = aa - bb
            rms = float(np.sqrt(np.mean(delta * delta)))
            max_abs = float(np.max(np.abs(delta)))
            status = "OK" if rms <= tol else ("FAIL" if enforce else "INFO")
            if status == "FAIL":
                failures.append(f"{group_name} {key}: component rms diff={rms:.3g}")
            unit_txt = (" " + unit) if unit else ""
            sa.slog.append(
                "Audit {0:<10s} {1:>18s}: component-rms diff {2:10.3e}{3}, max {4:10.3e}{3}, ndof {5:d} [{6}]\n".format(
                    group_name, key, rms, unit_txt, max_abs, int(aa.size), status
                )
            )

        sa.slog.append("########## Initial Sample Audit ############################### \n")
        for key in ("vib", "rot", "vel", "total"):
            gv = gen["energy"].get(key)
            av = ana["energy"].get(key)
            if gv is None or av is None:
                continue
            diff = abs(float(gv) - float(av))
            status = "OK" if diff <= ip.audit_initial_energy_tol else "FAIL"
            if status == "FAIL":
                failures.append(f"energy {key}: generation={gv:.8g}, analysis={av:.8g}, diff={diff:.3g} eV")
            sa.slog.append(
                "Audit energy {0:>5s}: generation {1:12.6f} eV, analysis {2:12.6f} eV, diff {3:10.3e} eV [{4}]\n".format(
                    key, float(gv), float(av), diff, status
                )
            )

        angular_tol = getattr(ip, "audit_initial_angular_tol", 0.0)
        angular_enforce = angular_tol > 0.0
        for key in ("L", "Jab", "J"):
            _audit_scalar("angular", key, gen["angular"].get(key), ana["angular"].get(key), angular_tol, "au", enforce=angular_enforce)
            _audit_vector("vector", key, gen["vectors"].get(key), ana["vectors"].get(key), angular_tol, "au", enforce=angular_enforce)
        for key in sorted(set(gen["molecule_scalars"]) & set(ana["molecule_scalars"])):
            _audit_scalar("mol scalar", key, gen["molecule_scalars"].get(key), ana["molecule_scalars"].get(key), angular_tol, "au", enforce=angular_enforce)
        for key in sorted(set(gen["molecule_vectors"]) & set(ana["molecule_vectors"])):
            _audit_vector("mol vector", key, gen["molecule_vectors"].get(key), ana["molecule_vectors"].get(key), angular_tol, "au", enforce=angular_enforce)
        for key in sorted(set(gen["scalars"]) & set(ana["scalars"])):
            if key == "relative_velocity_mps":
                vel_tol = getattr(ip, "audit_initial_velocity_tol", 0.0)
                _audit_scalar("scalar", key, gen["scalars"].get(key), ana["scalars"].get(key), vel_tol, "m/s", enforce=(vel_tol > 0.0))
            else:
                unit = "Ang" if key in ("b", "Jacobi_R") else "au"
                _audit_scalar("scalar", key, gen["scalars"].get(key), ana["scalars"].get(key), angular_tol, unit, enforce=angular_enforce)
        vib_tol = getattr(ip, "audit_initial_vib_tol", 0.0)
        for key in sorted(set(gen["vibrational"]) & set(ana["vibrational"])):
            _audit_component_vector(
                "vib",
                key,
                gen["vibrational"].get(key),
                ana["vibrational"].get(key),
                vib_tol,
                enforce=(vib_tol > 0.0),
            )
        for key in sorted(set(gen["angles"]) & set(ana["angles"])):
            diff = _angle_diff(gen["angles"][key], ana["angles"][key])
            status = "OK" if diff <= angular_tol else ("FAIL" if angular_enforce else "INFO")
            if status == "FAIL":
                failures.append(f"angle {key}: diff={diff:.3g} rad")
            sa.slog.append(
                "Audit {0:<10s} {1:>18s}: circular diff {2:10.3e} rad [{3}]\n".format(
                    "angle", key, diff, status
                )
            )
        if failures:
            raise ValueError("Initial sample audit failed for sample " + str(sa.sii) + ": " + "; ".join(failures))
        sa.slog.append("Initial sample audit: OK\n")

    def Mol2Image(self,sa):
        """Convert molecular coordinates and velocities to the image frame.

        This method transforms molecular coordinates and velocities to the image frame for calculations.
        """
        sa.sxx = np.concatenate([sa.mol[0].sxx, sa.mol[1].sxx])
        sa.svv = np.concatenate([sa.mol[0].svv, sa.mol[1].svv])

    def Image2Mol(self,sa):
        mol = sa.mol
        mol[0].sxx = sa.sxx[:mol[0].na,:]
        mol[1].sxx = sa.sxx[mol[0].na:,:]
        mol[0].svv = sa.svv[:mol[0].na,:]
        mol[1].svv = sa.svv[mol[0].na:,:]

    def ListXYZOut(self,sa, **dic): 
        """Generate output files for image coordinates and velocities.

        This method generates output files for image coordinates and velocities, storing them in the specified directory.
        """
        ip = self.ip
        sp = self.sp
        if "mess" in dic.keys():
            message = dic["mess"]
        else:
            message = " "
        if ip.units_out == "au":
          xscale = 1.0
          vscale = 1.0
          xunit = "au"
          vunit = "au"
        else:
          xscale = au2ang
          vscale = au2ang / au2fmt
          xunit = "Ang"
          vunit = "Ang/fmts"
        if 'logout' in dic.keys():
          sa.slog += [" Sample " + str(sa.id) + " Coordinates (" + xunit + ") : \n"]
        xyz = XYZlist(sp.el, sa.sxx * xscale,mess=message + " Coordinate (" + xunit + ")")
        debug = True
        debug = False
        if debug:
          vvo = dic['vvo'] 
          xyz.append('F ' + ''.join(["{0:14.7f}".format(f*0.5+sa.sxx[0,i]*au2ang) for i,f in enumerate(vvo.tolist()) ])+'\n')  
          print('XYZ = ')
          print(xyz)
          xyz[0] = '8\n'
        vxyz = XYZlist(sp.el, sa.svv * vscale,mess=message + " Velocities (" + vunit + ")")
        if not os.path.isdir(ip.dirout):
            os.system("mkdir " + ip.dirout)
        if 'logout' in dic.keys():
          sa.slog += xyz
          sa.slog += [" Sample " +str(sa.id) + " Velocities (" + vunit + ") : \n"]
          sa.slog += vxyz[2:]
        return xyz, vxyz

    def ImageXYZOut(self, xyz,vxyz,sii,sa=None):
        ip = self.ip
        fmt = ip.output_format
        if fmt in ("xyzvel", "both"):
          open(ip.dirout + "/" + ip.fileout +"_" + str(sii) + ".xyz", "w").writelines(xyz)
          open(ip.dirout + "/" + ip.fileout +"_" + str(sii) + ".vel", "w").writelines(vxyz)
        if fmt in ("npz", "both") and sa is not None:
          np.savez(
              ip.dirout + "/" + ip.fileout + "_" + str(sii) + ".npz",
              x=sa.sxx.copy(),
              v=sa.svv.copy(),
              units=ip.units_out,
          )
        return 

    def SetInterMolZVeloc(self,sa):
        sa.mol[0].svv -= sa.svel[0] 
        sa.mol[1].svv -= sa.svel[1]
        
    def SampleInterMolZVeloc(self,sa):
        """Sample intermolecular z-velocity for the scattering event.

        This method samples the intermolecular z-velocity based on molecular velocities or direct input parameters.
        """
        ip = self.ip
        sp = self.sp
        mol = self.mol
        msa = sa.mol 
        msp = [mol[0].sp,mol[1].sp]
        # if the user chose intermolecular energy/velocity directly
        if 'vel' in sa.dist.keys():
            if 'cont' in sa.dist['vel']:
                V = abs(ICDFsample(sa.dist['vel']['cont']))
            else:
                V = abs(sa.dist['vel']['v'])
            sa.sV  = V 
            vv = self.GetInterMolZVeloc(V)
            sa.svel = vv
        else:
            # if user chose molecular temperatures/velocities
            v1, log1 = mol[0].SampleZVeloc(msa[0])
            v2, log2 = mol[1].SampleZVeloc(msa[1])
            sa.slog += log1 + log2
            vv = self.GetInterMolZVelocFromMolV(v1, v2,sp.beamang)
            sa.sV = norm(vv[0,:]-vv[1,:])
            sa.svel = vv
        mm = array([sum(msp[0].mass), sum(msp[1].mass)])
        sa.smkin = np.sum(0.5 * vv.T**2 * mm, axis=0).tolist()
        vn = norm(vv)
        if 'vel' not in sa.SampInfo.keys():
            sa.SampInfo['vel'] = {}
        sa.SampInfo['vel']['ivel'] = sa.sV*au2mps 
        sa.SampInfo['vel']['velen'] = [k*au2ev for k in sa.smkin]
        sa.slog += info_section("intermolecular")
        sa.slog += [info_scalar("relative velocity", sa.sV * au2mps, "m/s")]
        sa.slog += [info_scalar("collision energy", sum(sa.smkin) * au2ev, "eV")]
        sa.slog += [info_scalar(mol[0].ip.name + " z velocity", vv[0, 2] * au2mps, "m/s")]
        sa.slog += [info_scalar(mol[1].ip.name + " z velocity", vv[1, 2] * au2mps, "m/s")]
        sa.slog += [info_scalar(mol[0].ip.name + " kinetic", sa.smkin[0] * au2ev, "eV")]
        sa.slog += [info_scalar(mol[1].ip.name + " kinetic", sa.smkin[1] * au2ev, "eV")]
        sa.slog += [f"{'orbital sampling':<{INFO_LABEL_WIDTH}} = {getattr(self.ip, 'orbital_sampling', 'geometric')}\n"]
        return


    def SampleOrbitalL(self,sa,**dic):
        cap = dic.get('cap', None)
        mode = getattr(self.ip, "orbital_sampling", "geometric")
        if mode == "flat-l":
          mxL = float(cap if cap is not None else self.ip.MaxL)
          if mxL <= 0:
            raise ValueError("flat-l orbital sampling requires a positive maxl/cap")
          L = float(ICDFsample(sa.dist['gen']['cont'])) * mxL
        else:
          L = ICDFsample(sa.dist['L']['cont'])
          if cap is not None:
            mxL = cap
            while L > mxL:
              L = ICDFsample(sa.dist['L']['cont'])
        sa.sL = L
        nL = np.sqrt(L*(L+1))
        sa.snL = nL
        phi = self.SampleImpactPhi(sa)
        cL = matmul(Rabout(phi,2),y)
        cL = nL*cL/norm(cL)
        sa.scL = cL
        log = [" - Orbital Angular Q.N. L = : " + "{0:3.2f}".format(L)+"\n"]
        if 'orb' not in sa.SampInfo.keys():
            sa.SampInfo['orb'] = {}
        sa.SampInfo['orb']['sampling'] = mode
        sa.SampInfo['orb']['iL'] = sa.sL 
        sa.SampInfo['orb']['inL'] = sa.snL 
        sa.SampInfo['orb']['icL'] = sa.scL 
        sa.SampInfo['orb']['sampled_impact_phi'] = phi
        return log

    def SampleImpactPhi(self, sa):
        ip = self.ip
        if ip.ImpactPhi is not None:
            return float(ip.ImpactPhi)
        return float(ICDFsample(sa.dist['gen']['cont'])) * tpi

    def SampleFixedImpactOrbital(self, sa):
        ip = self.ip
        sp = self.sp
        if ip.FixedB is None:
            raise ValueError("SampleFixedImpactOrbital called without fixed-b")
        if not hasattr(sa, "sV"):
            raise ValueError("fixed-b orbital setup requires relative velocity first")
        b = float(ip.FixedB)
        phi = self.SampleImpactPhi(sa)
        nL = sp.rmass * sa.sV * b
        L = 0.5 * (-1.0 + sqrt(1.0 + 4.0 * nL ** 2))
        cL = matmul(Rabout(phi, 2), y)
        cL = nL * cL / norm(cL)
        sa.sb = b
        sa.sphi = phi
        sa.sL = L
        sa.snL = nL
        sa.scL = cL
        if 'orb' not in sa.SampInfo.keys():
            sa.SampInfo['orb'] = {}
        sa.SampInfo['orb']['iL'] = sa.sL
        sa.SampInfo['orb']['inL'] = sa.snL
        sa.SampInfo['orb']['icL'] = sa.scL
        sa.SampInfo['orb']['sampled_impact_phi'] = phi
        sa.SampInfo['orb']['sampling'] = "fixed-b"
        sa.SampInfo['orb']['fixed_b'] = b
        return [
            *info_section("intermolecular"),
            f"{'orbital sampling':<{INFO_LABEL_WIDTH}} = fixed-b\n",
            info_scalar("fixed b", b * au2ang, "Ang", "{:14.5f}"),
            info_scalar("implied L", L, "", "{:14.2f}"),
        ]

    def SetStandardOrientation(self,sa):
         msa = sa.mol
         sa.slog += ['   *** Setting to Standard Azimuthal Orientation (Chi=0) *** \n']
         L = sa.scL
         if norm(L) >= 1e-8:
          uL = L/norm(L)
          chii = np.arctan2(uL[0], uL[1])
          Rs = Rabout(chii,2)
          Jab = msa[0].srpar[-1] + msa[1].srpar[-1]
          Jab = matmul(Rs,Jab) 
          L = matmul(Rs,L) 
          cJ = L + Jab
          if abs(norm(cJ)-sa.snJ) > 0.01: # this is dumb, but I cant be asked.... 
            quit('WORRY!')
          sa.scL = matmul(Rs,sa.scL)
          sa.scJ = matmul(Rs,sa.scJ)
          for mi in range(2):
            msa[mi].siJ = matmul(Rs,msa[mi].siJ)
            msa[mi].srpar[-1] = matmul(Rs,msa[mi].srpar[-1])
            msa[mi].sxx = matmul(msa[mi].sxx,Rs.T)
            msa[mi].svv = matmul(msa[mi].svv,Rs.T)
          sa.slog += ['       delChi = '+"{0:10.5f}".format(chii/pi)+' pi rad  \n']
         else:
          sa.slog += ['       delChi = '+"{0:10.5f}".format(0.0/pi)+' pi rad  \n']
        

    def OutputFrameRotation(self):
        frame = getattr(self.ip, "output_frame", "internal")
        if frame == "internal":
            return np.eye(3)
        if frame == "incoming-k-plus-z":
            return np.array([[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]])
        raise ValueError("Unknown output-frame convention: " + str(frame))

    def ApplyOutputFrameConvention(self, sa):
        """Rotate generated Cartesian/vector state into the requested output frame."""
        frame = getattr(self.ip, "output_frame", "internal")
        if frame == "internal" or getattr(sa, "_output_frame_applied", False):
            return
        R = self.OutputFrameRotation()

        def rot_vec(v):
            return matmul(R, np.asarray(v, dtype=float))

        def rot_rows(a):
            return matmul(np.asarray(a, dtype=float), R.T)

        def rotate_dict_vectors(dic, keys):
            if not isinstance(dic, dict):
                return
            for key in keys:
                val = dic.get(key)
                if val is None:
                    continue
                arr = np.asarray(val, dtype=float)
                if arr.shape == (3,):
                    dic[key] = rot_vec(arr)

        for msa in sa.mol:
            msa.sxx = rot_rows(msa.sxx)
            msa.svv = rot_rows(msa.svv)
            if hasattr(msa, "siJ"):
                msa.siJ = rot_vec(msa.siJ)
            if hasattr(msa, "srpar") and len(msa.srpar) > 0:
                msa.srpar[-1] = rot_vec(msa.srpar[-1])
            rotate_dict_vectors(
                msa.SampInfo.get("rot", {}),
                ("vecJ", "svecJs", "svecJm", "svecJc", "svecJ0", "svecJ0s"),
            )

        for attr in ("scL", "scJ", "sjab"):
            if hasattr(sa, attr):
                setattr(sa, attr, rot_vec(getattr(sa, attr)))
        if hasattr(sa, "svel"):
            sa.svel = rot_rows(sa.svel)

        rotate_dict_vectors(
            sa.SampInfo.get("orb", {}),
            ("icL", "icJ", "cL", "cJ", "Jab", "scL", "scJ", "sJab"),
        )
        sa._output_frame_applied = True
        sa.slog += info_section("frame convention")
        sa.slog += [info_frame_marker(frame)]
        sa.slog += info_frame_transform(frame)


    def StoreOrbitalInfoLog(self,sa):
        mol = self.mol
        sp = self.sp
        msa = sa.mol
        cL = sa.scL 
        cJ = sa.scJ
        nJ = sa.snJ
        nL = sa.snL
        ncL = norm(cL)
        ncJ = norm(cJ)
        #print('oooooooooooooooooo')
        #print('NL = ', ncL,nL)
        #print('NJ = ', ncJ,nJ)
        Ja, Jb = msa[0].srpar[-1], msa[1].srpar[-1] 
        Jab = Ja + Jb
        #print('jab = ', sa.snjab, ' = ', norm(Jab))
        iJa, iJb = msa[0].siJ, msa[1].siJ 
        iJab = iJa+iJb 

        nJa = norm(Ja)
        nJb = norm(Jb)
        nJab = norm(Jab)
        niJab = norm(iJab)
        niJa, niJb = norm(iJa), norm(iJb) 
        qJ    = 0.5 * (-1 + sqrt(1 + 4.0 * nJ ** 2)   ) 
        qL    = 0.5 * (-1 + sqrt(1 + 4.0 * nL ** 2)   ) 
        qJab  = 0.5 * (-1 + sqrt(1 + 4.0 * nJab ** 2) ) 
        qJa   = 0.5 * (-1 + sqrt(1 + 4.0 * nJa ** 2)  ) 
        qJb   = 0.5 * (-1 + sqrt(1 + 4.0 * nJb ** 2)  ) 
        qiJab = 0.5 * (-1 + sqrt(1 + 4.0 * niJab ** 2)) 
        qiJa  = 0.5 * (-1 + sqrt(1 + 4.0 * niJa ** 2) ) 
        qiJb  = 0.5 * (-1 + sqrt(1 + 4.0 * niJb ** 2) ) 
        if nL > 1e-5 and nJab > 1e-5:
          #L and Jab should be uncorrelated 
          cosLJab_thet = np.dot(cL,Jab)/(nL*nJab)
          # gen spherica polar coordinate angles between L and Jab (should be isotropic)
          _, R = rot_match_vec(reshape(cL,(1,3)),reshape(z,(1,3))) 
          _,LJab_be,LJab_al =  xyz2polar(matmul(R,Jab))
        if max(niJa, niJb, niJab, nJa, nJb, nJab) < 1.0e-10:
          sa.slog += [f"{'molecular J':<{INFO_LABEL_WIDTH}} = Ja, Jb, and Jab are zero in full and vector-model generation\n"]
        else:
          sa.slog += [info_vec("Ja, full", iJa, "au", "Ja", niJa, qiJa)]
          sa.slog += [info_vec("Jb, full", iJb, "au", "Jb", niJb, qiJb)]
          sa.slog += [info_vec("Jab, full", iJab, "au", "Jab", niJab, qiJab)]
          sa.slog += [info_vec("Ja, vector model", Ja, "au", "Ja", nJa, qJa)]
          sa.slog += [info_vec("Jb, vector model", Jb, "au", "Jb", nJb, qJb)]
          sa.slog += [info_vec("Jab, vector model", Jab, "au", "Jab", nJab, qJab)]
        sa.slog += [info_vec("L", cL, "au", "L", ncL, qL)]
        sa.slog += [info_vec("J = L + Jab", cJ, "au", "J", ncJ, qJ)]
        # semiclassical mapping
        b = nL/(sp.rmass*sa.sV)
        phi = np.arctan2(-cL[0],cL[1])
        if 'orb' not in sa.SampInfo.keys():
            sa.SampInfo['orb'] = {}
        sa.SampInfo['orb']['cL'] = cL 
        sa.SampInfo['orb']['cJ'] = cJ
        sa.SampInfo['orb']['Jab'] = Jab
        sa.SampInfo['orb']["ncL"] = nL
        sa.SampInfo['orb']["ncJ"] = nJ
        sa.SampInfo['orb']["nJab"] = nJab
        if nL > 1e-5 and nJab > 1e-5:
          sa.SampInfo['orb']["LJab_al"] = LJab_al
          sa.SampInfo['orb']["LJab_be"] = LJab_be
          sa.SampInfo['orb']["cosLJab_thet"] = cosLJab_thet
        sa.SampInfo['orb']['b'] = b
        sa.SampInfo['orb']['phi'] = phi
        sa.sb = b 
        sa.sphi = phi
        sa.slog += [info_scalar("b", float(b) * au2ang, "Ang", "{:14.5f}")]
        sa.slog += [info_scalar("impact phi", phi / pi, "pi rad", "{:14.5f}")]
        return   

    # in comes lab fixed z-velocity magnitude, out goes z-velocity
    # in centre of mass
    def GetInterMolZVelocFromMolV(self, V1, V2,ang):
        """Get intermolecular z-velocity from molecular velocities.

        This method calculates the intermolecular z-velocity from the molecular velocities of two molecules.

        Args:
            V1 (float): Z-velocity of the first molecule.
            V2 (float): Z-velocity of the second molecule.

        Returns:
            numpy.ndarray: Inter-molecular z-velocity.
        """
        v2 = -V2*z  
        v1 = matmul(Rabout(ang,0),-z)*V1
        VV = norm(v2-v1)
        return self.GetInterMolZVeloc(VV) 

    def GetInterMolZVeloc(self, V):
        """Calculate intermolecular z-velocity from given velocity magnitude.

        This method calculates the intermolecular z-velocity from a given velocity magnitude.

        Args:
            V (float): Inter-molecular velocity magnitude.

        Returns:
            numpy.ndarray: Inter-molecular z-velocity for both molecules.
        """
        sp = self.sp 
        vv = zeros((2, 3))
        vv[0, :] = (z * V)  * sp.w1
        vv[1, :] = -(z * V) * sp.w0
        return vv

    def SetInterZDist(self,sa, Rz):
        """Set the intermolecular z-coordinate distance.

        This method sets the intermolecular z-coordinate distance for both molecules.

        Args:
            Rz (float): Inter-molecular distance along the z-axis.
        """
        mol = sa.mol
        sp = self.sp 
        mol[0].sxx += (z * Rz) * sp.w1
        mol[1].sxx -= (z * Rz) * sp.w0
        sa.slog += info_section("intermolecular")
        sa.slog += [info_scalar("Rz", Rz * au2ang, "Ang", "{:14.5f}")]

    def SetImpactParam(self,sa, b, phi):
        """Set the impact parameter for scattering.

        This method sets the impact parameter and azimuthal angle (phi) for the scattering process.

        Args:
            b (float): Impact parameter.
            phi (float): Azimuthal angle.
        """
        mol = sa.mol
        sp = self.sp 
        if hasattr(sa, "svel") and hasattr(sa, "scL"):
            p_rel = -sp.rmass * (sa.svel[0] - sa.svel[1])
            p2 = float(np.dot(p_rel, p_rel))
            if p2 > 1.0e-20:
                d = np.cross(p_rel, sa.scL) / p2
            else:
                d = b * (y * sin(phi) + x * cos(phi))
        else:
            d = b * (y * sin(phi) + x * cos(phi))
        mol[0].sxx += d * sp.w1
        mol[1].sxx -= d * sp.w0

    def SampleRigidRotorState0(self,sa,**dic):
        """Sample the rigid rotor state for each molecule.

        This method samples the rigid rotor state for each molecule, considering their respective temperatures.
        """
        debug = False
        log = info_section("rotation")
        mol = self.mol
        msa = sa.mol 
        for i in range(2):
            if mol[i].sp.na > 1 and 'rotJ' in msa[i].dist:
                mol[i].SampleTotMolAngMom(msa[i])
                log += mol[i].SampleRigidRotorState(msa[i])
            else:
                log += [f"{mol[i].ip.name:<{INFO_LABEL_WIDTH}} = no rotational state\n"]
        return log

    # sets angular velocity from angular momentum once :
    def SetAngularVelocity(self,sa,**dic):
        slog = []
        mol = self.mol
        msa = sa.mol
        for i in range(2):
            if mol[i].sp.na > 1:
                slog += mol[i].SetAngularVeloc(msa[i],msa[i].srpar[-1])
        if 'printlog' in dic.keys():
          return slog

    def SampleOrientat0(self,sa):
        """Sample molecular orientation.

        This method samples molecular orientation in three-dimensional space for each molecule.
        """
        debug = True
        debug = False
        mol = self.mol
        msa = sa.mol 
        log = info_section("orientation")
        for i in range(2):
            if mol[i].sp.na > 1 and 'ori' in msa[i].dist:
                mol[i].SampleRotation(msa[i])
            else:
                msa[i].soR = np.eye(3)
                log += [f"{mol[i].ip.name:<{INFO_LABEL_WIDTH}} = no orientational state\n"]
        j0, j1 = matmul(msa[0].soR,msa[0].srpar[-1]), matmul(msa[1].soR,msa[1].srpar[-1])
        jab = j0+j1 
        log += self.SetAngularVelocity(sa,printlog=True)
        for i in range(2):
            if mol[i].sp.na > 1:
                log += mol[i].SetOrientat(msa[i],printlog=True)
            else:
                log += [f"{mol[i].ip.name:<{INFO_LABEL_WIDTH}} = orientation unchanged for atom\n"]
        sa.sjab = jab
        sa.snjab = norm(jab)
        return log

    def SampleHOVibrState(self,sa):
        """Sample harmonic oscillator vibrational states.

        This method samples the harmonic oscillator vibrational states for each molecule.
        """
        sa.slog += info_section("vibration")
        mol = self.mol
        msa = sa.mol 
        sa.slog += mol[0].SampleHOVibrState(msa[0])
        sa.slog += mol[1].SampleHOVibrState(msa[1])

    def CalcRotEner(self,sa):
        """Calculate rotational energies.

        This method calculates rotational energies for each molecule and the overall system.
        """
        debug = False
        sa.slog += info_section("rotation")
        mol = self.mol
        msa = sa.mol 
        sa.slog += mol[0].CalcRotEner(msa[0])
        sa.slog += mol[1].CalcRotEner(msa[1])
        if 'rot' not in sa.SampInfo.keys(): 
           sa.SampInfo['rot'] = {}
        sa.SampInfo['rot']['m0'] = msa[0].SampInfo.get('rot', {}).copy()
        sa.SampInfo['rot']['m1'] = msa[1].SampInfo.get('rot', {}).copy()

        if debug:
          sa.slog += [" -Ja         = "+''.join(["{0:10.5f}".format(j)+' ' for j in  msa[0].SampInfo['rot']['svecJs']])+" \n"]
          sa.slog += [" -J b        = "+''.join(["{0:10.5f}".format(j)+' ' for j in  msa[1].SampInfo['rot']['svecJs']])+" \n"]
          sa.slog += ["             = "+''.join(["{0:10.5f}".format(j)+' ' for j in  msa[0].SampInfo['rot']['svecJs']+msa[1].SampInfo['rot']['svecJs']])+" \n"]

    def CalcInterEner(self,sa):
        """Calculate vibrational energies.

        This method calculates vibrational energies for each molecule and the overall system.
        """
        sa.slog += info_section("vibration")
        mol = self.mol
        msa = sa.mol 
        sa.slog += mol[0].CalcInterEner(msa[0])
        sa.slog += mol[1].CalcInterEner(msa[1])
        if 'vib' not in sa.SampInfo.keys(): 
           sa.SampInfo['vib'] = {}
        sa.SampInfo['vib']['m0'] = msa[0].SampInfo.get('vib', {}).copy()
        sa.SampInfo['vib']['m1'] = msa[1].SampInfo.get('vib', {}).copy()

    def CalcOrient(self,sa):
        """Calculate molecular orientation.

        This method calculates the molecular orientation and Euler angles for each molecule.
        """
        sa.slog += info_section("orientation")
        mol = self.mol
        msa = sa.mol 
        sa.slog += mol[0].CalcOrient(msa[0])
        sa.slog += mol[1].CalcOrient(msa[1])
        if 'ori' not in sa.SampInfo.keys(): 
           sa.SampInfo['ori'] = {}
        sa.SampInfo['ori']['m0'] = msa[0].SampInfo.get('ori', {}).copy()
        sa.SampInfo['ori']['m1'] = msa[1].SampInfo.get('ori', {}).copy()

    def GenSamples(self, **dic):
        """Generate multiple scattering samples.

        This method generates multiple scattering samples, allowing for the customization of the number of samples and other parameters.

        Args:
            dic (dict): Additional parameters for sample generation.
        """
        ip = self.ip
        sp = self.sp
        if "N" in dic.keys():
            N = dic["N"]
        else:
            N = ip.Nsamp
        if ip.progress != "quiet":
          print("Generating " + str(N) + " Samples")
        if ip.check_input or ip.dry_run:
          self.log += ["Input check/dry-run requested; skipping sample generation.\n"]
          return
        start_offset = 0
        if ip.continues:
          wks =self.loadworkers()
        else:
          wks = False 
        if wks != False: 
          self.loaddata() 
          if ip.KeepInfo: 
            self.loadinfo()
          try:
            start_offset = int(sum(len(sa.sdat['vel']['ivel']) for sa in wks))
          except Exception:
            start_offset = 0
        else:
          ip.continues = False
          wks = [ self.InitializeWorker(i) for i in range(ip.nwork)  ] 
        if ip.progress == "verbose":
          print('WKS = ', wks)
        if not ip.usewang:
         sp.td = []
        else:
          wang_path = self._runpath('wang.pkl')
          if ip.run_mode == "rebuild-wang" and os.path.exists(wang_path):
            raise ValueError(
                "run-mode = rebuild-wang was requested, but "
                + wang_path
                + " already exists.\nMove or rename the existing wang.pkl first; "
                + "icats will not overwrite it automatically."
            )
          if os.path.exists(wang_path):
            sp.uu, sp.iwld, sp.td, warning = wang.load_validated(
                wang_path, wang.metadata_from_input(ip)
            )
            if warning:
              self.log += [warning]
          else:
             self.GenerateWang(wks)
        self.sdat = self.PrepareSdat()
        base = int(N // ip.nwork) if ip.nwork > 0 else 0
        rem = int(N % ip.nwork) if ip.nwork > 0 else 0
        ranges = []
        cursor = int(start_offset)
        for i in range(ip.nwork):
          size = base + (1 if i < rem else 0)
          ranges.append([cursor, cursor + size])
          cursor += size
        wks = Parallel(n_jobs=ip.nwork)(
            delayed(self.GenerateSample)(wks[i], ranges[i]) for i in range(ip.nwork)
        )
        #for i,sa in enumerate(wks): 
        #  print(i, 'sa = ', sa.sdat)
        if ip.printout[1]:
         slog = [] 
         for sa in wks: 
          slog += sa.slog
         open(ip.fileout +"_full.info", "w").writelines(slog)
        if ip.printout[0]:
         slog = []
         vslog = []
         for sa in wks: 
          slog += sa.ixyz
          vslog += sa.ivxyz
         open(ip.fileout +"_full.xyz", "w").writelines(slog)
         open(ip.fileout +"_full.vel", "w").writelines(vslog)
        self.sdat = self.MergeSdats(wks)
        self._write_costheta_convergence()
        self.saveworkers(wks)
        self.savedata()
        if ip.KeepInfo: 
           self.saveinfo()
        if ip.hist_sampled:
            self.PlotSamples()
#        _ = plot_wl_weights(self.ww)
 

    def MergeSdats(self,wks): 
      fsdat = wks[0].sdat.copy()
      debug = True
      debug = False
      for sa in wks[1:]:
        sdat =sa.sdat 
        for kys in sdat.keys():
          debug and print('kys = ', kys)
          for ky in sdat[kys].keys(): 
            debug and print('ky = ', ky)
            if 'm0' == ky or 'm1' == ky:
              if 'rot' == kys:
                fsdat[kys][ky]['J'] += sdat[kys][ky]['J']
                for J in set(sdat[kys][ky]['J']):
                  if J not in fsdat[kys][ky].keys():
                    fsdat[kys][ky][J] = {}
                  for k1 in sdat[kys][ky][J].keys():
                    if 'jz' == k1 or 'sjz' == k1 or 'qjz' == k1:
                      if k1 not in fsdat[kys][ky][J].keys():
                        fsdat[kys][ky][J][k1] = []
                      fsdat[kys][ky][J][k1] += sdat[kys][ky][J][k1] 
                    else: 
                      if k1 not in fsdat[kys][ky][J].keys():
                        fsdat[kys][ky][J][k1] = {}
                      for k2 in sdat[kys][ky][J][k1].keys():
                        if k2 not in fsdat[kys][ky][J][k1].keys():
                          fsdat[kys][ky][J][k1][k2] = []
                        fsdat[kys][ky][J][k1][k2] += sdat[kys][ky][J][k1][k2] 
              elif 'vib' == kys:
                fsdat[kys][ky]['vi'] += sdat[kys][ky]['vi']
                for vi in range(self.mol[int(ky[-1])].ip.MaxV+1):
                 for k in ['Q','sQ','P','sP']:
                   fsdat[kys][ky][k][vi] += sdat[kys][ky][k][vi]
              else:
                for k in sdat[kys][ky].keys():
                  fsdat[kys][ky][k] += sdat[kys][ky][k]
            else:
              fsdat[kys][ky] += sdat[kys][ky]
      return fsdat

    def PrepareSdat(self): 
      debug = True
      debug = False
      sdat = {}
      ky = 'vel'
      sdat[ky] = {}
      #plot intermolecular velocities and energies: 
      sdat[ky]['ivel'], sdat[ky]['velen'] = [], []
      ky = 'orb'
      sdat[ky] = {}
      # plot lengths of total angular momentum, orbital angular momentum, Jab 
      #sdat[ky]['sncJ'], sdat[ky]['sncL'], sdat[ky]['snJab'] = [], [], [] 
      sdat[ky]['iJ'], sdat[ky]['iL'] = [], []
      sdat[ky]['sJ'], sdat[ky]['sL'], sdat[ky]['sJab'] = [], [], [] 
      sdat[ky]['LJab_be'], sdat[ky]['LJab_al'], sdat[ky]["cosLJab_thet"] = [], [], []
      # plot orbital cylibdrical coordinates 
      sdat[ky]['sb'], sdat[ky]['sphi'] = [], [] 
      # plot intermolecular COM kinetic-energy shares separately; the raw
      # senergy record is a two-component vector and should not be flattened
      # into one mixed histogram.
      sdat[ky]['senergy_m0_ev'], sdat[ky]['senergy_m1_ev'], sdat[ky]['senergy_total_ev'] = [], [], []
      #plot rotational body-fixed coordinates: 
      ky = '2bJac'
      sdat[ky] = {}
      sdat[ky]['phi'], sdat[ky]['beta'], sdat[ky]['theta'], sdat[ky]['chi'] = [], [], [], []
      sdat[ky]['alpha1'], sdat[ky]['beta1'], sdat[ky]['gamma1'] = [], [], []
      sdat[ky]['alpha2'], sdat[ky]['beta2'], sdat[ky]['gamma2'] = [], [], []
      #plot molecular orientations 
      ky = 'ori'
      sdat[ky] = {}
      sdat[ky]['m0'], sdat[ky]['m1'] = {}, {}
      for m in sdat[ky].keys():
        sdat[ky][m]['salpha'], sdat[ky][m]['sbeta'], sdat[ky][m]['sgamma'] = [], [] , []
        sdat[ky][m]['alpha'], sdat[ky][m]['beta'], sdat[ky][m]['gamma'] = [], [] , []
        sdat[ky][m]['sphi'], sdat[ky][m]['stheta'], sdat[ky][m]['schi'] = [], [] , []
        sdat[ky][m]['phi'], sdat[ky][m]['theta'], sdat[ky][m]['chi'] = [], [] , []
      #molecular angular momentum 
      ky = 'rot'
      sdat[ky] = {}
      sdat[ky]['m0'], sdat[ky]['m1'] = {}, {}
      # each projection depends on J
      # each vector model distribuition depends on J and its projection pz and axis ax
      for m in sdat[ky].keys():
        #molecular angular momentum J, projection and classical vector model 
        sdat[ky][m]['J'] = []
      #molecular vibrational coordinates 
      ky = 'vib'
      sdat[ky] = {}
      sdat[ky]['m0'], sdat[ky]['m1'] = {}, {}
      for im,m in enumerate(sdat[ky].keys()):
        #print('m = ', m, im, ' s = ', self.mol[im].nm)
        # vibrational state, and position and coordinate modes.. since the modes are standardized, we only need to plot one for each vibrational state. 
        sdat[ky][m] = {} 
        sdat[ky][m]['vi'] = [[] for _ in range(self.mol[im].sp.nm)] 
        sdat[ky][m]['Q'],  sdat[ky][m]['P'] = [[] for i in range(self.mol[im].ip.MaxV+1) ], [[] for i in range(self.mol[im].ip.MaxV+1) ]
        sdat[ky][m]['sQ'], sdat[ky][m]['sP'] = [[] for i in range(self.mol[im].ip.MaxV+1)], [[] for i in range(self.mol[im].ip.MaxV+1)]
        sdat[ky][m]['senergy'] = [[] for i in range(self.mol[im].ip.MaxV+1)]
         # vibrational energy  
      for si, info in enumerate(self.sampls['info']):
        if debug:
          print('######',si)
          for kys in info.keys():
             print('########### ', kys)
             for ky in info[kys].keys():
                 print('        ###  ', ky)
                 if 'm0' == ky or 'm1' == ky:
                    for kk in info[kys][ky]:
                       print('          #', kk)
      return sdat

    def AddInfoToSamples(self,sdat,info,**dic):
      debug = False 
      if 'debug' in dic.keys():
        debug = True
        print('INFO = ', info)
      for kys in info.keys():
        debug  and print('kys=',kys)
        for ky in info[kys].keys():
          debug  and print('ky=',ky)
          if 'm0' == ky or 'm1' == ky:
            mi = int(ky[-1])
            if 'rot' == kys:
              debug and print('YES')
              sjz = info[kys][ky]['sjz'] 
              jz = info[kys][ky]['jz'] 
              qjz = info[kys][ky]['qjz'] 
              J   = info[kys][ky]['J'] 
              sdat[kys][ky]['J'].append(J)
              debug  and print('ky=',ky, 'append = ', len(sdat[kys][ky]['J']))
              if J not in sdat[kys][ky].keys(): 
                sdat[kys][ky][J] = {}
                sdat[kys][ky][J]['jz'] = []
                sdat[kys][ky][J]['qjz'] = []
                sdat[kys][ky][J]['sjz'] = []
              sdat[kys][ky][J]['jz'].append(jz)
              sdat[kys][ky][J]['sjz'].append(sjz)
              sdat[kys][ky][J]['qjz'].append(qjz)
              if 'idjz' in info[kys][ky].keys():
                id = info[kys][ky]['idjz']
                if id not in sdat[kys][ky][J].keys():
                  sdat[kys][ky][J][id] = {'sbet':[],'sgamm':[], 'bet':[],'gamm':[]}
                sdat[kys][ky][J][id]['sbet'].append(info[kys][ky]['sbet'])
                sdat[kys][ky][J][id]['sgamm'].append(info[kys][ky]['sgamm'])
                sdat[kys][ky][J][id]['bet'].append(info[kys][ky]['bet'])
                sdat[kys][ky][J][id]['gamm'].append(info[kys][ky]['gamm'])
            elif 'vib' == kys:
              if 'vi' not in info[kys][ky].keys():
                 continue
              vi = info[kys][ky]['vi']
              for i,v in enumerate(vi):
               sdat[kys][ky]['vi'][i].append(v)
               for k in ['Q','sQ','P','sP']:
                 sdat[kys][ky][k][v] += [info[kys][ky][k][self.mol[mi].sp.ntr+i]]
               sdat[kys][ky]['senergy'][v] += [info[kys][ky]['senergy'][i]]
            else:
              for k in info[kys][ky].keys():
                if k in sdat[kys][ky].keys():
                  sdat[kys][ky][k].append(info[kys][ky][k])
          elif kys == 'orb' and ky == 'senergy':
            vals = np.asarray(info[kys][ky], dtype=float).ravel()
            if len(vals) > 0:
              sdat[kys]['senergy_m0_ev'].append(float(vals[0] * au2ev))
            if len(vals) > 1:
              sdat[kys]['senergy_m1_ev'].append(float(vals[1] * au2ev))
            if len(vals) > 0:
              sdat[kys]['senergy_total_ev'].append(float(np.sum(vals) * au2ev))
          elif ky in sdat[kys].keys():
            sdat[kys][ky].append(info[kys][ky]) 
      return 

    def PlotSamples(self): 
      sdat = self.sdat
      print('generating histograms...')
      debug = False
      #debug = True
      hist, edg = {}, {}
      for kys in sdat.keys():
        if debug:
         print('###################', kys)
        hist[kys], edg[kys] = {}, {}
        for ky in sdat[kys].keys(): 
          if debug:
           print('   ###################', ky)
          hist[kys][ky], edg[kys][ky] = {}, {}
          if 'm0' == ky or 'm1' == ky:
            if 'rot' == kys:
#              hist[kys][ky]['J'], edg[kys][ky]['J'] = np.histogram(sdat[kys][ky]['J'],bins='auto')
              hist_emit(sdat[kys][ky]['J'], "J",
                        stage="sampled", scope=f"molecule_{ky}")
              for J in set(sdat[kys][ky]['J']): 
                if debug:
                 print('      ###################', J)
                hist[kys][ky][J], edg[kys][ky][J] = {}, {}
                for k1 in sdat[kys][ky][J].keys(): 
                  if debug:
                   print('         ###################', k1, ' ln  = ',len(sdat[kys][ky][J][k1]))
                  if 'jz' == k1 or 'sjz' == k1 or 'qjz' == k1:
#                    hist[kys][ky][J][k1], edg[kys][ky][J][k1] = np.histogram(sdat[kys][ky][J][k1],bins='auto')
                    hist_emit(sdat[kys][ky][J][k1], f"J{J}_{k1}",
                              stage="sampled", scope=f"molecule_{ky}")
                  else: 
                    hist[kys][ky][J][k1], edg[kys][ky][J][k1] = {}, {}
                    for k2 in sdat[kys][ky][J][k1].keys(): 
                      if debug:
                       print('            ###################', k2, ' ln  = ',len(sdat[kys][ky][J][k1][k2]))
#                      hist[kys][ky][J][k1][k2], edg[kys][ky][J][k1][k2] = np.histogram(sdat[kys][ky][J][k1][k2],bins='auto')
                      hist_emit(sdat[kys][ky][J][k1][k2], f"J{J}_i{k1}_{k2}",
                                stage="sampled", scope=f"molecule_{ky}")
            elif 'vib' == kys: 
              for mi, vlist in enumerate(sdat[kys][ky]['vi']):
                hist_emit(vlist, f"vi_mode{mi}",
                          stage="sampled", scope=f"molecule_{ky}")
              for vi in range(self.mol[int(ky[-1])].ip.MaxV+1):
               for k in ['Q','sQ','P','sP','senergy']:
                 hist_emit(sdat[kys][ky][k][vi], f"v{vi}_{k}",
                           stage="sampled", scope=f"molecule_{ky}")
            else:
              for k in sdat[kys][ky].keys():
                if debug:
                 print('      ###################', k, ' ln  = ',len(sdat[kys][ky][k]))
#                hist[kys][ky][k], edg[kys][ky][k] = np.histogram(sdat[kys][ky][k],bins='auto')
                hist_emit(sdat[kys][ky][k], str(k),
                          stage="sampled", scope=f"molecule_{ky}")
          else:
            if debug:
              print('                      ', ky, ' ln  = ',len(sdat[kys][ky]))
#            hist[kys][ky], edg[kys][ky] = np.histogram(sdat[kys][ky],bins='auto')
            hist_emit(sdat[kys][ky], f"{kys}_{ky}",
                      stage="sampled", scope="system")
          
          

    def AnalyseSample(self,sa):
        """
        Analyze a sample by performing various calculations and generating a summary.

        This method calculates various physical properties of the sample and appends the results to a log.

        It consists of the following steps:
        1. Calculate Euler orientation of the molecular frame from position and velocity (self.CalcOrient()).
        2. Calculate rotational information from position and velocity (self.CalcRotEner()).
        3. Calculate vibrational information from position and velocity (self.CalcInterEner()).
        4. Calculate intermolecular information from position and velocity (self.CalcInterMolMomentum()).
        5. Summarize and log energy information (self.SummarizeLogEnergy()).
        6. Generate output including final coordinates and velocity.

        The results are stored in the 'slog', 'sampls', and 'SampInfo' attributes of the object.

        Returns:
        None
        """
        ip = self.ip
        sp = self.sp
        sa.slog += info_header(sa.sii, "analysis")
        sa.slog += [info_frame_marker(ip.output_frame)]
        sa.slog += info_frame_transform(ip.output_frame)
        # Calculates Euler orientation of molecular frame from xx and vv
        self.CalcOrient(sa)
        # Calculates Rotational information from xx and vv
        self.CalcRotEner(sa)
        # Calculates Vibrational information from xx and vv
        self.CalcInterEner(sa)
        # Calculate InterMolecular information from xx and vv
        self.CalcInterMolMomentum(sa)
        # Calculate Jacobi coordinates 
        self.CalcJacobiCoordinates(sa)
        # write up some energy information
        self.SummarizeLogEnergy(sa,True)
        # generate output
        sa.slog += info_section("coordinates and velocities")
        if ip.KeepInfo:
          self.sampls['cv'].append([sa.sxx.copy(), sa.svv.copy()])
          self.sampls['info'].append(sa.SampInfo.copy())
        else: 
          self.sampls['cv'].append([sa.sxx.copy(), sa.svv.copy()]) 

    def ReadSamples(self,sa, filx, filv):
        """
        Read sample data from files, analyze each sample, and write results to an output file.

        Parameters:
        filx (str): The filename containing position data (XYZs).
        filv (str): The filename containing velocity data.

        This method reads data from 'filx' and 'filv', assigns the data to the appropriate properties
        of the 'mol' object, and then analyzes each sample using the 'AnalyseSample' method. The results
        are stored in 'self.slog' and written to an output file named '[self.fileout]_samples.info'.

        Note:
        - The 'ReadXYZs' function is used to read position and velocity data from the specified files.
        - The 'AnalyseSample' method is called for each sample to perform the analysis.

        Returns:
        None
        """
        el, xyzs, xmess = ReadXYZs(filx)
        el, vels, vmess = ReadXYZs(filv)
        mol = sa.mol
        n0 = self.mol[0].sp.na
        for i, x in enumerate(xyzs):
            v = vels[i]
            m = xmess[i] + ' & ' + vmess[i]
            sa.sii = i
            sa.SampInfo = {}
            # Ensure per-molecule sample defaults exist for analysis-only paths.
            self.mol[0].InitializeSample(mol[0])
            self.mol[1].InitializeSample(mol[1])
            mol[0].sxx = x[:n0, :]*ang2au 
            mol[1].sxx = x[n0:, :]*ang2au 
            mol[0].svv = v[:n0, :]*ang2au/fmt2au
            mol[1].svv = v[n0:, :]*ang2au/fmt2au
            sa.slog += [ '###### Sample Name ' + m + '\n']
            self.AnalyseSample(sa)


iscatter = icats


# Example Usage:
if __name__ == "__main__":
    # Create an instance of the icats class
    sc = icats()

    # Read input data from a file
    sc.ReadInput("input.txt")

    # Generate scattering samples
    sc.GenSamples(N=10)
