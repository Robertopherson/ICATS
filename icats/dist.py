#!/usr/bin/env python3
from .constants import *
from .functions import *
       
def VibPartFunc(n, w, T, izpe):
    """
    Calculate the vibrational partition function.

    Parameters:
    n (array): Array of quantum numbers for vibrational states.
    w (float): Vibrational frequency.
    T (float): Temperature.
    izpe (Bool): Flag whether or not to include ZPE.

    Returns:
    float: The vibrational partition function.
    """
    nn = array([int(np.round(na)) for na in n])
    for n in nn:
        if n < 0:
            return 0.0
    if izpe:
     ee = w * (nn + 0.5)
    else:
     ee = w * nn 
    return prod(exp(-ee / (kboltz * T)))

def VibPartFuncICDF(rng,w,T, izpe,VMax): 
    """
    The inverse CDF of the vibrational partition function.

    Parameters:
    rng (array): Array of RNG objects 
    w (float): Vibrational frequency vector.
    T (float): Temperature.

    Returns:
    float: The inverse of the  vibrational partition function, n.
    float: The boltzmann weighting  for that n  
    """
    nn = np.array([100000,100000] )
    while any(nn > VMax): 
      u = np.array([ rn.random() for rn in rng]) 
      nn = np.floor( -np.log( 1- u ) * kboltz * T / w  ) 
    if izpe:
     ee = w * (nn + 0.5)
    else:
     ee = w * nn 
    ee = prod(exp(-ee / (kboltz * T)))
    return nn, ee  


def iCDFGenGauss(A,v0,n):
    """
    Returns inverse CDF of a gaussian x^n exp(-A(x-x0)^2) for x > 0
    
    Parameters
    ----------
    A, v0, n : floats           # PDF parameters
    n      : int (≈ 1e4–1e5)    # grid size for tabulation
    tail_prob : float           # probability mass to leave outside the domain
    
    Returns
    -------
    v_grid : (n,) ndarray       # velocity grid
    cdf    : (n,) ndarray       # cumulative distribution (0→1)
    icdf(u): callable           # vectorised inverse CDF, u in (0,1)
    """
    nn=50000
    tail_prob=1e-8
    # --- 1. finite domain ----------------------------------------------------
    # An equivalent Gaussian width is sigma = 1/sqrt(2A).
    sigma = 1.0 / np.sqrt(2.0 * A)
    v_min = 0.0
    #  erf^{-1}(1-ϵ) ≈ √(ln(2/ϵ)) ; choose enough sigmas to cover tail_prob
    n_sig = np.sqrt(np.log(2.0 / tail_prob))
    v_max = v0 + n_sig * sigma
    # --- 2. tabulate PDF and CDF --------------------------------------------
    v_grid = np.linspace(v_min, v_max, nn)
    pdf    = v_grid**n * np.exp(-A * (v_grid - v0)**2)
    # normalise numerically
    dv     = v_grid[1] - v_grid[0]
    pdf   /= pdf.sum() * dv
    # cumulative integral (simple trapezoidal cumtrapz inline)
    cdf = np.cumsum(pdf) * dv
    # ensure last point is exactly 1.0
    cdf[-1] = 1.0
    # --- 3. build the iCDF via interpolation --------------------------------
    def icdf(u):
        """
        Vectorised inverse CDF using np.interp (monotone → O(log n) search).
        u must be in (0,1).  Returns v with same shape as u.
        """
        return np.interp(u, cdf, v_grid) 
    return cdf, v_grid 

def GenGaussSample(rng,cdf, v_grid):
    x = np.interp(rng[0].random(), cdf, v_grid)
    return x, 0.0

def MaxwellBoltzmann(v, T, m):
    """
    Calculate the Maxwell-Boltzmann distribution.

    Parameters:
    v (float): Velocity.
    T (float): Temperature.
    m (float): Mass.

    Returns:
    float: The Maxwell-Boltzmann distribution value.
    """
    #print('v = ', v, 'm = ', m, 'T = ', T, 'kb = ', kboltz, ' RETURN = ', v**2 * exp(-m * v**2 / (2 * kboltz * T)) )
    return exp(-m * v**2 / (2 * kboltz * T))


def MBiCDF(rng,A):
    """
    Draw 'size' samples from f(v) ∝ v^2 exp(-A v^2), v ≥ 0.
    A = m/(2 k_B T);  sigma = 1/sqrt(2A)
    """
    sigma = 1.0 / np.sqrt(2.0 * A) 
    vxyz  = np.random.normal(0.0, sigma, (1, 3))   # 3 Gaussians
    return np.linalg.norm(vxyz, axis=1), 0.0 


def RadialMaxwellBoltzmann(v, T, m):
    """
    Calculate the Maxwell-Boltzmann distribution in intermolecular coordinate V_R.

    Parameters:
    v (float): Velocity along R.
    T (float): Temperature.
    m (float): Mass.

    Returns:
    float: The Maxwell-Boltzmann distribution value.
    """
    return v**2 * exp(-m * v**2 / (2 * kboltz * T))

def GaussianICDF(rng,x0,A):
    sigma = 1/np.sqrt(2*A)
    u1 = rng[0].random()
    u2 = rng[0].random()
    r = np.sqrt(-2*np.log(u1)) * sigma
    theta  = 2*np.pi*u2 
    x = r * np.cos(theta) + x0
    return x, GaussianF2(x,x0,sigma)


def GaussianF(x, x0, fwhm):
    """
    Gaussian function.

    Parameters:
    x (float): Input value.
    x0 (float): Mean of Gaussian
    fwhm (float): Full-width half maximum

    Returns:
    float: Gaussian function value.
    """
    sigma = fwhm/(2*sqrt(2.0*log(2.0)))
    a     = 1.0/(2*sigma**2)
    return exp(-a * (x-x0)**2)

def GaussianF2(x, x0, sigma):
    """
    Gaussian function.

    Parameters:
    x (float): Input value.
    x0 (float): Mean of Gaussian
    sigma (float): standard deviation  

    Returns:
    float: Gaussian function value.
    """
    a     = 1.0/(2*sigma**2)
    return exp(-a * (x-x0)**2)

def IPICDF(rng,MaxB):
    u = rng[0].random()
    return np.array([np.sqrt(u)*MaxB]), u

def JcrossICDFc(rng, iso, MaxJ):
    """ICDF for continuous 2J+1 weighting on [0, Jmax]."""
    u = rng[0].random()
    if iso:
      J = sqrt(u * ((MaxJ + 1)**2 - 1) + 1) - 1
    else:
      J = u*MaxJ 
    return J, J 


def JcrossICDF(rng,iso,MaxJ):
    rho = np.arange(0,MaxJ)
    if iso:
      prob = 2*rho +1 
      prob = prob / prob.sum()
      J = rng[0].choice(rho, p=prob)
      return J, 2*J+1 
    else: 
      J = rng[0].choice(rho)
      return J, J
 
def IPDist(b, MaxB):
    """
    Calculate the Impact Parameter Distribution.

    Parameters:
    b (float): Impact parameter.
    MaxB (float): Maximum impact parameter.

    Returns:
    float: The Impact Parameter Distribution value.
    """
    if abs(b) > MaxB:
        return 0.0
    else:
        return abs(b)

def HarmWigner(QP, n):
    """
    Calculate the Wigner distribution for harmonic oscillators.

    Parameters:
    QP (tuple): Tuple containing Q and P coordinates.
    n (int): Quantum number.

    Returns:
    float: Wigner distribution value.
    """
    Q, P = QP
    def facfac_loop(n):
        yield 0, 1.
        r = 1.
        for m in range(1, n + 1):
            r *= float(n - m + 1) / m**2
            yield m, r
    def ana_laguerre(n, x):
        total = 0.
        for m, r in facfac_loop(n):
            entry = (-1.)**m * r * x**m
            total += entry
        return total
    if n == 0:
        return exp(-Q**2) * exp(-P**2)
    else:
        rhosquare = 2.0 * (P**2 + Q**2)
        W = (-1.0)**n * ana_laguerre(n, rhosquare) * exp(-rhosquare / 2.0)
        return max([W, 0.0])

def HarmHusimi(QP, n): 
    """
    Calculate the Husimi distribution for harmonic oscillators.

    Parameters:
    QP (tuple): Tuple containing Q and P coordinates.
    n (int): Quantum number.

    Returns:
    float: Husimi distribution value.
    """
    Q, P = QP
    if n == 0:
        return exp(-Q**2) * exp(-P**2)
    else:
        rhosquare = 2.0 * (P**2 + Q**2)
        W = (rhosquare**2)**n * exp(-rhosquare / 2.0)
        return W

def gamma_mt(n,rng):
      """Generate one sample from Gamma(n+1, 1) using Marsaglia & Tsang's method."""
      k = n+1
      d = k - 1.0 / 3.0
      c = 1.0 / np.sqrt(9.0 * d)
      while True:
          x = rng.normal()
          v = (1 + c * x)**3
          if v <= 0:
              continue
          u = rng.random()
          if u < 1 - 0.0331 * (x**4):
              return d * v
          if np.log(u) < 0.5 * x**2 + d * (1 - v + np.log(v)):
              return d * v

def HusimiFuncICDF(rng,n):  
    """
    The inverse CDF of the Husimi distribution.

    Parameters:
    rng (array): Array of RNG objects (ONE ONE ELEMENT)
    n (integer): vibrational state 

    Returns:
    float: The inverse of the  ICDF ( the momenta and positions).
    float: The Husimi distributiion value for that vibrational state 
    """
    r = np.sqrt(gamma_mt(n,rng[0]))
    theta = rng[0].random()*tpi
    Q = r * np.cos(theta)
    P = r * np.sin(theta)
    return [[Q,P]], HarmHusimi([Q,P],n) 

def IsoRotorTotBoltzICDF(rng,T,A,Jmax):
    abet = A/(kboltz*T)
    J = 1000000 
    while J > Jmax:
     J =  np.sqrt( - np.log(1-rng[0].random()) / abet ) - 0.5
    sigma = 1/(np.sqrt(2*A/(kboltz*T))+1e-10 )
    return np.round(J), (2J+1)*GaussianF2(J+0.5,0.0,sigma)

def AniRotorTotBoltzICDF(rng,T,A,Jmax):
    sigma = 1/(np.sqrt(2*A/(kboltz*T))+1e-10)
    J = 1000000 
    while J > Jmax:
     u1 = rng[0].random()
     u2 = rng[0].random()
     r = np.sqrt(-2*np.log(u1)) * sigma
     theta = 2*np.pi*u2
     J = abs(r * np.cos(theta)) - 0.5
    return np.round(J), GaussianF2(J+0.5,0.0,sigma)

def SymRotorProjBoltzICDF(rng,T,A,J):
    if A > 0.0:   #for prolates 
      sigma = 1/np.sqrt(2*A/(kboltz*T))
      u1 = rng[0].random()
      u2 = rng[0].random()
      r = np.sqrt(-2*np.log(u1)) * sigma
      theta = 2*np.pi*u2
      K = np.round(max(min(r * np.cos(theta),J),-J))
      return K, GaussianF2(K,0.0,sigma)
    else:  # for oblates 
      beta = 1.0/(kboltz*T)
      KK = np.arange(-J, J+1)
      rho = np.exp(-beta * A * KK**2)   
      prob = rho / rho.sum()
      K = rng[0].choice(KK, p=prob)
      return K, np.exp(-beta * A * K**2)
 

def reject_bingham(rng,sigma):
    g = rng[0].multivariate_normal(mean=np.zeros(3), cov=np.diag(sigma), size=None)
    # 2) Normalize each to unit length
    u = g / norm(g)
    ang = iR2q(rot_match_vec(z.reshape((1,3)), u.reshape((1,3)) )[1])
    return ang, 1

def reject_GaussSin(rng,jz,sig,J): 
    # node jz is actually the positive half jz, so need to include +/- at the beta leve
    if J == 0:
     return np.array([0.0]), 0.0
    bet0 = np.arccos(jz/np.sqrt(J*(J+1)))
    bets = sig/(np.sqrt( J*(J+1) - jz**2  )+1e-8)
    if bets < 1e-5: 
      return np.array([bet0]), 1.0
    while True: 
      bet = rng[0].normal(bet0, bets)
      if 0 <= bet <= np.pi:
          if rng[0].random() <= np.sin(bet): 
            if rng[0].random() >= 0.5:
              bet = pi-bet
              bet0 = pi-bet0
            return np.array([bet]), GaussianF2(bet,bet0,bets)*sin(bet) 

def reject_phi(rng, Jx2, Jy2):    #used as  arotGam key
    while True:
        phi = rng[0].random() * 2 * np.pi - np.pi
        ff = Jx2 * np.cos(phi)**2 + Jy2 * np.sin(phi)**2
        if rng[0].random() <=  ff / max(Jx2,Jy2):
            return np.array([phi]), ff 

def uniform(rng,min,max):
    return np.array([min + rng[0].random() * ( max-min )]), 1

def AsymRotorProjBoltzICDF(rng,T,rr,J):
    if J >= len(rr):
      print('WARNING J = ',J, 'rr = ', len(rr))
    ee1 = rr[J][1]
    pdf = np.exp(-ee1/(T*kboltz))
    pdf = pdf/pdf.sum()
    #Build discrete inverse CDF (step-function CDF)
    cdf = np.concatenate(([0.0], np.cumsum(pdf)))  # length = len(ee1)+1
    u = rng[0].random()
    # find index i such that cdf[i] < u <= cdf[i+1] ⇒ energy = ee1[i]
    ii = np.searchsorted(cdf, u, side='right') - 1
    #ee,ax,jx2,jy2,jz2,jz,sig = [o[ii] for o in oo[J]]
    return [o[ii] for o in rr[J]], 1


def Boltzmann(E, T): 
    """
    Calculate the Boltzmann distribution.

    Parameters:
    E (float): Energy (au).
    T (float): Temperature.

    Returns:
    float: The Boltzmann distribution value (numerator).
    """
    return exp(-E / (kboltz * T))

def VibBoltzmann(w,ni,T):
    if T <= 0:
        ni = np.asarray(ni)
        return np.where(ni == 0, 1.0, 0.0)
    bet = 1/(kboltz*T) 
    Z = 1/(1-np.exp(-bet*w))
    p = np.exp(-bet*w*ni) / Z 
    return p 


def REBoltzmann(xi, EE, JJ, T):
    """
    Calculate the Rotational Energy Boltzmann distribution.

    Parameters:
    xi (float): Quantum number.
    EE (list): List of energy levels.
    JJ (list): List of angular momentum quantum numbers.
    T (float): Temperature.

    Returns:
    float: The Rotational Energy Boltzmann distribution value.
    """
    ii = int(np.round(xi))
    if ii < 0:
        return 0.0
    elif ii >= len(EE):
        return 0.0
    else:
        E, sigma = EE[ii], JJ[ii] * 2 + 1
        return sigma * Boltzmann(E, T)

def AniChiPerpICDF3(rng, j, jab):
    j2, jab2, jjab = j*j, jab*jab, j*jab
    def L(chi): return np.sqrt(j2 + jab2 - 2*jjab*np.cos(chi))
    def f(chi): return (2*L(chi) + 1) * jjab * np.sin(chi) / (1e-10+L(chi))
    if j <= 0.0 or jab == 0.0:
        chi = rng[0].random()*np.pi
        return np.array([chi]), f(chi)

    cc = np.linspace( np.pi, 5000)
    fvals = f(cc)
    mx = fvals.max()
    while True:
        chi = rng[0].random()*np.pi     
        if rng[0].random() <= f(chi)/mx:
            return np.array([chi]), f(chi)

def AniChiPerpICDFX(rng):
    chi = np.arccos(1-2*rng[0].random())
    return np.array([chi]), sin(chi)

def AniChiPerpICDF(rng, j, jab):
    if j <= 0.0 or jab == 0.0:
        return np.array([0.0]), 0.0  # trivial
    j2, jab2, jjab = j*j, jab*jab, j*jab
    def LL(chi): return np.sqrt(j2 + jab2 - 2*jjab*np.cos(chi))
    def f(chi): return (2*LL(chi) + 1) * jjab * np.sin(chi) / (1e-10+LL(chi))
    Lmin = abs(j - jab)
    Lmax = j + jab
    u = rng[0].random()
    L = np.sqrt(u*(Lmax**2 - Lmin**2) + Lmin**2)
    coschi = (j**2 + jab**2 - L**2)/(2*j*jab)
    coschi = np.clip(coschi, -1.0, 1.0)
    chi = np.arccos(coschi)
    # randomize the sign symmetry if needed
    if rng[0].random() < 0.5:
        chi = 2*np.pi - chi  # or just use -chi in your convention
    return np.array([chi]), f(chi)



def AniChiPerpICDF4(rng, J, J_AB):
    """
    Draw χ from  P(χ) ∝ ((L + 0.5)/k) * sin(χ) / L .
    Constant 1/k cancels in normalisation but is left explicit.
    """

    j2, jab2, jjab = J*J, J_AB*J_AB, J*J_AB
    def LL(chi): return np.sqrt(j2 + jab2 - 2*jjab*np.cos(chi))
    def f(chi): return (2*LL(chi) + 1) * jjab * np.sin(chi) / (1e-10+LL(chi))
    JJ = J * J_AB
    L_min = abs(J - J_AB)
    
    # --- envelope: biggest value of (L+0.5)/L is at L=L_min
    M = (L_min + 0.5) / L_min            # ensures p_accept ≤ 1
    i = 0  
    while True:
        i += 1
        # 1. proposal from ½ sin χ
        chi = np.arccos(1 - 2*rng[0].random())

        # 2. compute L(χ)
        L = np.sqrt(J**2 + J_AB**2 - 2 * JJ * np.cos(chi))

        # 3. acceptance probability
        p_accept = ((L + 0.5) / L) / M   # 1/k cancels out
        if i == 100000:
          return np.array([0.0]), 0.0
        if rng[0].random() < p_accept:
            return np.array([chi]), f(chi) 



#def AniChiPerpICDF3(rng, j, jab):
#    if j <= 0.0 or jab == 0.0:
#        return np.array([0.0]), 0.0
#    j2, jab2, jjab = j*j, jab*jab, j*jab
#    def L(chi): return np.sqrt(j2 + jab2 - 2*jjab*np.cos(chi))
#    def f(chi): return (2*L(chi) + 1) * jjab * np.sin(chi) / (1e-10+L(chi))
#    cc = np.linspace( np.pi, 5000)
#    fvals = f(cc)
#    mx = fvals.max()
#    while True:
#        chi = rng[0].random()*np.pi     
#        if rng[0].random() <= f(chi)/mx:
#            return np.array([chi]), f(chi)
#

def AniChiPerpICDF2(rng,j,jab): 
    j2 = j**2                    
    jab2 = jab**2                
    jjab = j*jab                 
    def L(chi):                  
        return np.sqrt(j2 + jab2 - 2*jjab*np.cos(chi))  
    def f(chi):                  
        return (2*L(chi) + 1) * jjab * np.sin(chi) / L(chi)  
    def g(l,chi):  
        return np.cos(chi)*(2*l + 1)*l**2 - jjab*np.sin(chi)**2  
    if j <= 0.0 or jab == 0.0:  #chi can be anything   
      return np.array([0.0]), 0.0  
    if jab+j < 1.5:    
      cc = np.linspace(0,np.pi,5000)                                                      
      chi = cc[np.argmax(f(cc))]                                                            
    else:  # for most values of j and jab the peak is at near pi/2                          
      cc = np.pi/2.0 + np.linspace(-0.1,0.1,50)                                             
      chi = cc[np.argmax(f(cc))]                                                            
    h = 1e-7                                                                                
    max_iter = 50                                                                   
    for _ in range(max_iter):                                                       
       l = L(chi)                                                                   
       gg = g(l,chi)                                                                
       chi += gg*h                                                                  
    mx = f(chi)                                                                     
    while True:                                                                     
       chi = rng[0].random()*np.pi                                                
       if rng[0].random() <= f(chi)/mx:                                             
         return np.array([chi]), f(chi)  


def EulerSurface(angs):
    """
    Calculate the Euler surface.

    Parameters:
    angs (list): List of angles.

    Returns:
    float: The Euler surface value.
    """
    if angs[1] < 0.0 or angs[1] > pi:
#        print('NEYOND!',1)
        return 0.0
    for i in [0, 2]:
        if angs[i] < -pi or angs[i] > pi:
 #           print('NEYOND!',i)
            return 0.0
    return sin(angs[1])

def EulerPolarSurface(angs, alpha, beta):
    """
    Calculate the Euler Polar surface.

    Parameters:
    angs (list): List of angles.
    alpha (float): Parameter.
    beta (float): Parameter.

    Returns:
    float: The Euler Polar surface value.
    """
    if angs[1] < 0.0 or angs[1] > pi:
        return 0.0
    for i in [0, 2]:
        if angs[i] < -pi or angs[i] > pi:
            return 0.0
    return sin(angs[1]) * 0.5 * (1 - alpha * beta * cos(angs[1]))

def IsotropicDistICDF(rng,rotpar,rn1,rn2):
    d1 = rng[0].random(3)*rn1-rn2
    d1[1] = np.arccos(1-2*d1[1])
    if rotpar == 'xyz':
       d1[1] = d1[1]-np.pi/2.0
    return d1, 1  

