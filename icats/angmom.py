from numpy import ndarray, array, zeros, linspace, matmul, sqrt, diag, arange, complex128, cos, sin, exp, eye
from numpy.linalg import eigh, norm
from math import lgamma, pi, factorial

def absM(M):
    el, ev = eigh(M)
    return matmul(ev,(matmul(diag(abs(el)),ev.conj().T)))

def J2(j: int) -> ndarray:
    """
    Calculate the square of the total angular momentum operator for a given quantum number.

    Parameters:
    j (int): Quantum number.

    Returns:
    ndarray: The squared total angular momentum operator.
    """
    return Jz2(j) + Jx2(j) + Jy2(j) 

def Jz2(j: int) -> ndarray:
    """
    Calculate the square of the z-component of the angular momentum operator for a given quantum number.

    Parameters:
    j (int): Quantum number.

    Returns:
    ndarray: The squared z-component of the angular momentum operator.
    """
    return matmul(J_z(j), J_z(j))

def Jx2(j: int) -> ndarray:
    """
    Calculate the square of the x-component of the angular momentum operator for a given quantum number.

    Parameters:
    j (int): Quantum number.

    Returns:
    ndarray: The squared x-component of the angular momentum operator.
    """
    return matmul(J_x(j), J_x(j))

def Jy2(j: int) -> ndarray:
    """
    Calculate the square of the y-component of the angular momentum operator for a given quantum number.

    Parameters:
    j (int): Quantum number.

    Returns:
    ndarray: The squared y-component of the angular momentum operator.
    """
    return matmul(J_y(j), J_y(j))

def J_x(j: int) -> ndarray:
    """
    Calculate the x-component of the angular momentum operator for a given quantum number.

    Parameters:
    j (int): Quantum number.

    Returns:
    ndarray: The x-component of the angular momentum operator.
    """
    return (J_plus(j) + J_minus(j)) / 2

def J_y(j: int) -> ndarray:
    """
    Calculate the y-component of the angular momentum operator for a given quantum number.

    Parameters:
    j (int): Quantum number.

    Returns:
    ndarray: The y-component of the angular momentum operator.
    """
    return (J_plus(j) - J_minus(j)) / (2j)

def J_z(j: int) -> ndarray:
    """
    Calculate the z-component of the angular momentum operator for a given quantum number.

    Parameters:
    j (int): Quantum number.

    Returns:
    ndarray: The z-component of the angular momentum operator.
    """
    return diag([-j + i for i in range(int(2 * j + 1))])

def Jx(j: int) -> ndarray:
    """
    Alias for J_x(j).

    Parameters:
    j (int): Quantum number.

    Returns:
    ndarray: The x-component of the angular momentum operator.
    """
    return J_x(j)

def Jy(j: int) -> ndarray:
    """
    Alias for J_y(j).

    Parameters:
    j (int): Quantum number.

    Returns:
    ndarray: The y-component of the angular momentum operator.
    """
    return J_y(j)

def Jz(j: int) -> ndarray:
    """
    Alias for J_z(j).

    Parameters:
    j (int): Quantum number.

    Returns:
    ndarray: The z-component of the angular momentum operator.
    """
    return J_z(j)

def P2(j: int) -> ndarray:
    """
    Calculate the square of the total angular momentum operator for a given quantum number.

    Parameters:
    j (int): Quantum number.

    Returns:
    ndarray: The squared total angular momentum operator.
    """
    return J2(j) 

def Pz2(j: int) -> ndarray:
    """
    Calculate the square of the z-component of the angular momentum operator for a given quantum number.

    Parameters:
    j (int): Quantum number.

    Returns:
    ndarray: The squared z-component of the angular momentum operator.
    """
    return matmul(P_z(j), P_z(j))

def Px2(j: int) -> ndarray:
    """
    Calculate the square of the x-component of the angular momentum operator for a given quantum number.

    Parameters:
    j (int): Quantum number.

    Returns:
    ndarray: The squared x-component of the angular momentum operator.
    """
    return matmul(P_x(j), P_x(j))

def Py2(j: int) -> ndarray:
    """
    Calculate the square of the y-component of the angular momentum operator for a given quantum number.

    Parameters:
    j (int): Quantum number.

    Returns:
    ndarray: The squared y-component of the angular momentum operator.
    """
    return matmul(P_y(j), P_y(j))

def Pp2(j: int) -> ndarray:
    """
    Calculate the square of the raising operator for a given quantum number.

    Parameters:
    j (int): Quantum number.

    Returns:
    ndarray: The squared raising operator.
    """
    return matmul(P_plus(j), P_plus(j))

def Pm2(j: int) -> ndarray:
    """
    Calculate the square of the lowering operator for a given quantum number.

    Parameters:
    j (int): Quantum number.

    Returns:
    ndarray: The squared lowering operator.
    """
    return matmul(P_minus(j), P_minus(j))

def P_x(j: int) -> ndarray:
    """
    Calculate the x-component of the angular momentum operator for a given quantum number.

    Parameters:
    j (int): Quantum number.

    Returns:
    ndarray: The x-component of the angular momentum operator.
    """
    return (P_plus(j) + P_minus(j)) / 2

def P_y(j: int) -> ndarray:
    """
    Calculate the y-component of the angular momentum operator for a given quantum number.

    Parameters:
    j (int): Quantum number.

    Returns:
    ndarray: The y-component of the angular momentum operator.
    """
    return (P_plus(j) - P_minus(j)) / (2j)

def P_z(j: int) -> ndarray:
    """
    Calculate the z-component of the angular momentum operator for a given quantum number.

    Parameters:
    j (int): Quantum number.

    Returns:
    ndarray: The z-component of the angular momentum operator.
    """
    return diag([-j + i for i in range(int(2 * j + 1))])

def Px(j: int) -> ndarray:
    """
    Alias for P_x(j).

    Parameters:
    j (int): Quantum number.

    Returns:
    ndarray: The x-component of the angular momentum operator.
    """
    return P_x(j)

def Py(j: int) -> ndarray:
    """
    Alias for P_y(j).

    Parameters:
    j (int): Quantum number.

    Returns:
    ndarray: The y-component of the angular momentum operator.
    """
    return P_y(j)

def Pz(j: int) -> ndarray:
    """
    Alias for P_z(j).

    Parameters:
    j (int): Quantum number.

    Returns:
    ndarray: The z-component of the angular momentum operator.
    """
    return P_z(j)

def P_minus(j: int) -> ndarray:
    """
    Calculate the lowering operator for a given quantum number.

    Parameters:
    j (int): Quantum number.

    Returns:
    ndarray: The lowering operator.
    """
    if int(2 * j + 1) != 2 * j + 1:
        raise ValueError(f"j must be a half-integer. Found: {j}")
    dim = int(2 * j + 1)
    mat = zeros((dim, dim))
    m_prime_list = linspace(-j, j, dim)
    m_list = linspace(-j, j, dim)
    for row, m_prime in enumerate(m_prime_list):
        for col, m in enumerate(m_list):
            mat[row, col] = J_plus_component(j, m_prime, j, m)
    return mat

def P_plus(j: int) -> ndarray:
    """
    Calculate the raising operator for a given quantum number.

    Parameters:
    j (int): Quantum number.

    Returns:
    ndarray: The raising operator.
    """
    if int(2 * j + 1) != 2 * j + 1:
        raise ValueError(f"j must be a half-integer. Found: {j}")
    dim = int(2 * j + 1)
    mat = zeros((dim, dim))
    m_prime_list = linspace(-j, j, dim)
    m_list = linspace(-j, j, dim)
    for row, m_prime in enumerate(m_prime_list):
        for col, m in enumerate(m_list):
            mat[row, col] = J_minus_component(j, m_prime, j, m)
    return mat

def J_plus(j: int) -> ndarray:
    """
    Calculate the raising operator for a given quantum number.

    Parameters:
    j (int): Quantum number.

    Returns:
    ndarray: The raising operator.
    """
    if int(2 * j + 1) != 2 * j + 1:
        raise ValueError(f"j must be a half-integer. Found: {j}")
    dim = int(2 * j + 1)
    mat = zeros((dim, dim))
    m_prime_list = linspace(-j, j, dim)
    m_list = linspace(-j, j, dim)
    for row, m_prime in enumerate(m_prime_list):
        for col, m in enumerate(m_list):
            mat[row, col] = J_plus_component(j, m_prime, j, m)
    return mat

def J_minus(j: int) -> ndarray:
    """
    Calculate the lowering operator for a given quantum number.

    Parameters:
    j (int): Quantum number.

    Returns:
    ndarray: The lowering operator.
    """
    if int(2 * j + 1) != 2 * j + 1:
        raise ValueError(f"j must be a half-integer. Found: {j}")
    dim = int(2 * j + 1)
    mat = zeros((dim, dim))
    m_prime_list = linspace(-j, j, dim)
    m_list = linspace(-j, j, dim)
    for row, m_prime in enumerate(m_prime_list):
        for col, m in enumerate(m_list):
            mat[row, col] = J_minus_component(j, m_prime, j, m)
    return mat

def J_plus_component(j_prime: int, m_prime: int, j: int, m: int) -> float:
    """
    Get the matrix element of the raising operator.

    Parameters:
    j_prime (int): Quantum number of the final state.
    m_prime (int): Magnetic quantum number of the final state.
    j (int): Quantum number of the initial state.
    m (int): Magnetic quantum number of the initial state.

    Returns:
    float: The matrix element of the raising operator.
    """
    if (j_prime != j) or (m_prime != m + 1):
        return 0
    return J_plus_coefficient(j, m)

def J_minus_component(j_prime: int, m_prime: int, j: int, m: int) -> float:
    """
    Get the matrix element of the lowering operator.

    Parameters:
    j_prime (int): Quantum number of the final state.
    m_prime (int): Magnetic quantum number of the final state.
    j (int): Quantum number of the initial state.
    m (int): Magnetic quantum number of the initial state.

    Returns:
    float: The matrix element of the lowering operator.
    """
    if (j_prime != j) or (m_prime != m - 1):
        return 0
    return J_minus_coefficient(j, m)

def J_plus_coefficient(j: int, m: int) -> float:
    """
    Calculate the coefficient for the raising operator.

    Parameters:
    j (int): Quantum number.
    m (int): Magnetic quantum number.

    Returns:
    float: The coefficient for the raising operator.
    """
    return sqrt((j - m) * (j + m + 1))

def J_minus_coefficient(j: int, m: int) -> float:
    """
    Calculate the coefficient for the lowering operator.

    Parameters:
    j (int): Quantum number.
    m (int): Magnetic quantum number.

    Returns:
    float: The coefficient for the lowering operator.
    """
    return sqrt((j + m) * (j - m + 1))



def WangTran(J):
  isqrt = sqrt(2) * 0.5
  nJ = 2 * J + 1
  U = zeros((nJ, nJ))
  k = 0
  nd = [0,0,0,0]
  if J == 0:
   return array([[1.0]]), ['A'], [[0,1]]
  def Ep(U,J,k,nd):
    # E+ Even K, Even combination 
    for j in range(2, J+1, 2):
      U[k, J + j] =  isqrt
      U[k, J - j] =  isqrt
      k += 1
      nd += 1
    return k, nd
  def Em(U,J,k,nd):
    # E- Even K, Odd combination 
    for j in range(2, J+1, 2):
      U[k, J + j] =  isqrt
      U[k, J - j] = -isqrt
      k += 1
      nd += 1
    return k, nd
  def Op(U,J,k,nd):
    # O+ Odd K, Even combination 
    for j in range(1, J+1, 2):
      U[k, J + j] =  isqrt
      U[k, J - j] =  isqrt
      k += 1
      nd += 1
    return k, nd
  def Om(U,J,k,nd):
    # O- Odd K, Odd combination 
    for j in range(1, J+1, 2):
      U[k, J + j] =  isqrt
      U[k, J - j] = -isqrt
      k += 1
      nd += 1
    return k, nd
  idx = [[0]]
  ####### A symmetry:
  if J%2 == 0: # Even J
    U[k, J] = 1  # A irrep (K = 0) 
    nd[0] = 1
    k += 1
    
    k, nd[0] = Ep(U,J,k,nd[0]) 
  else:   # Odd J  
    k, nd[0] = Em(U,J,k,nd[0])
  idx[-1].append(k)
  idx.append([k]) 
  ####### Bc (Bz) symmetry:  
  if J%2 == 0: # Even J 
    k, nd[1] = Em(U,J,k,nd[1]) 
  else:   # Odd J  
    U[k, J] = 1  # Bz irrep (K = 0) 
    nd[1] = 1
    k += 1

    k, nd[1] = Ep(U,J,k,nd[1]) 
  idx[-1].append(k)
  idx.append([k]) 
  ####### Bb (By) symmetry:  
  if J%2 == 0: # Even J 
    k, nd[2] = Om(U,J,k,nd[2]) 
  else:   # Odd J  
    k, nd[2] = Op(U,J,k,nd[2]) 
  idx[-1].append(k)
  idx.append([k]) 
  ####### Ba (Bx) symmetry:  
  if J%2 == 0: # Even J 
    k, nd[3] = Op(U,J,k,nd[3]) 
  else:   # Odd J  
    k, nd[3] = Om(U,J,k,nd[3]) 
  idx[-1].append(k)
  irreps = ['A']*nd[0] + ['B_z']*nd[1] + ['B_y']*nd[2] + ['B_x']*nd[3]
  return U.T, irreps, idx

def W2(j: int) -> ndarray:
    """
    Calculate the square of the total angular momentum operator for a given quantum  number (Wang Representation).

    Parameters:
    j (int): Quantum  number (Wang Representation).

    Returns:
    ndarray: The squared total angular momentum operator.
    """
    U, *dum = WangTran(j)
    return matmul(U.conj().T,matmul(P2(j),U))

def Wz2(j: int) -> ndarray:
    """
    Calculate the square of the z-component of the angular momentum operator for a given quantum  number (Wang Representation).

    Parameters:
    j (int): Quantum  number (Wang Representation).

    Returns:
    ndarray: The squared z-component of the angular momentum operator.
    """
    U, *dum = WangTran(j)
    return matmul(U.conj().T,matmul(Pz2(j),U))

def Wx2(j: int) -> ndarray:
    """
    Calculate the square of the x-component of the angular momentum operator for a given quantum  number (Wang Representation).

    Parameters:
    j (int): Quantum  number (Wang Representation).

    Returns:
    ndarray: The squared x-component of the angular momentum operator.
    """
    U, *dum = WangTran(j)
    return matmul(U.conj().T,matmul(Px2(j),U))

def Wy2(j: int) -> ndarray:
    """
    Calculate the square of the y-component of the angular momentum operator for a given quantum  number (Wang Representation).

    Parameters:
    j (int): Quantum  number (Wang Representation).

    Returns:
    ndarray: The squared y-component of the angular momentum operator.
    """
    U, *dum = WangTran(j)
    return matmul(U.conj().T,matmul(Py2(j),U))

def Wp2(j: int) -> ndarray:
    """
    Calculate the square of the raising operator for a given quantum  number (Wang Representation).

    Parameters:
    j (int): Quantum  number (Wang Representation).

    Returns:
    ndarray: The squared raising operator.
    """
    U, *dum = WangTran(j)
    return matmul(U.conj().T,matmul(Pp2(j),U))

def Wm2(j: int) -> ndarray:
    """
    Calculate the square of the lowering operator for a given quantum  number (Wang Representation).

    Parameters:
    j (int): Quantum  number (Wang Representation).

    Returns:
    ndarray: The squared lowering operator.
    """
    U, *dum = WangTran(j)
    return matmul(U.conj().T,matmul(Pm2(j),U))

def Wx(j: int) -> ndarray:
    """
    Alias for P_x(j).

    Parameters:
    j (int): Quantum  number (Wang Representation).

    Returns:
    ndarray: The x-component of the angular momentum operator.
    """
    U, *dum = WangTran(j)
    return matmul(U.conj().T,matmul(Px(j),U))

def Wy(j: int) -> ndarray:
    """
    Alias for P_y(j).

    Parameters:
    j (int): Quantum  number (Wang Representation).

    Returns:
    ndarray: The y-component of the angular momentum operator.
    """
    U, *dum = WangTran(j)
    return matmul(U.conj().T,matmul(Py(j),U))

def Wz(j: int) -> ndarray:
    """
    Alias for P_z(j).

    Parameters:
    j (int): Quantum  number (Wang Representation).

    Returns:
    ndarray: The z-component of the angular momentum operator.
    """
    U, *dum = WangTran(j)
    return matmul(U.conj().T,matmul(Pz(j),U))

def W_minus(j: int) -> ndarray:
    """
    Calculate the lowering operator for a given quantum  number (Wang Representation).

    Parameters:
    j (int): Quantum  number (Wang Representation).

    Returns:
    ndarray: The lowering operator.
    """
    U, *dum = WangTran(j)
    return matmul(U.conj().T,matmul(P_minus(j),U))

def W_plus(j: int) -> ndarray:
    """
    Calculate the raising operator for a given quantum  number (Wang Representation).

    Parameters:
    j (int): Quantum  number (Wang Representation).

    Returns:
    ndarray: The raising operator.
    """
    U, *dum = WangTran(j)
    return matmul(U.conj().T,matmul(P_plus(j),U))

def aWx(j: int) -> ndarray:
    """
    Alias for P_x(j).

    Parameters:
    j (int): Quantum  number (Wang Representation).

    Returns:
    ndarray: The x-component of the angular momentum operator.
    """
    U, *dum = WangTran(j)
    return matmul(U.conj().T,matmul(aPx(j),U))

def aWy(j: int) -> ndarray:
    """
    Alias for P_y(j).

    Parameters:
    j (int): Quantum  number (Wang Representation).

    Returns:
    ndarray: The y-component of the angular momentum operator.
    """
    U, *dum = WangTran(j)
    return matmul(U.conj().T,matmul(aPy(j),U))

def aWz(j: int) -> ndarray:
    """
    Alias for P_z(j).

    Parameters:
    j (int): Quantum  number (Wang Representation).

    Returns:
    ndarray: The z-component of the angular momentum operator.
    """
    U, *dum = WangTran(j)
    return matmul(U.conj().T,matmul(aPz(j),U))

def aW_minus(j: int) -> ndarray:
    """
    Calculate the lowering operator for a given quantum  number (Wang Representation).

    Parameters:
    j (int): Quantum  number (Wang Representation).

    Returns:
    ndarray: The lowering operator.
    """
    U, *dum = WangTran(j)
    return matmul(U.conj().T,matmul(P_minus(j),U))

def aW_plus(j: int) -> ndarray:
    """
    Calculate the raising operator for a given quantum  number (Wang Representation).

    Parameters:
    j (int): Quantum  number (Wang Representation).

    Returns:
    ndarray: The raising operator.
    """
    U, *dum = WangTran(j)
    return matmul(U.conj().T,matmul(P_plus(j),U))

def aP_x(j: int) -> ndarray:
    """
    Calculate the x-component of the angular momentum operator for a given quantum number.

    Parameters:
    j (int): Quantum number.

    Returns:
    ndarray: The x-component of the angular momentum operator.
    """
    return absM((P_plus(j) + P_minus(j)) / 2)

def aP_y(j: int) -> ndarray:
    """
    Calculate the y-component of the angular momentum operator for a given quantum number.

    Parameters:
    j (int): Quantum number.

    Returns:
    ndarray: The y-component of the angular momentum operator.
    """
    return absM((P_plus(j) - P_minus(j)) / (2j))

def aP_z(j: int) -> ndarray:
    """
    Calculate the z-component of the angular momentum operator for a given quantum number.

    Parameters:
    j (int): Quantum number.

    Returns:
    ndarray: The z-component of the angular momentum operator.
    """
    return absM(diag([-j + i for i in range(int(2 * j + 1))]))

def aPx(j: int) -> ndarray:
    """
    Alias for P_x(j).

    Parameters:
    j (int): Quantum number.

    Returns:
    ndarray: The x-component of the angular momentum operator.
    """
    return aP_x(j)

def aPy(j: int) -> ndarray:
    """
    Alias for P_y(j).

    Parameters:
    j (int): Quantum number.

    Returns:
    ndarray: The y-component of the angular momentum operator.
    """
    return aP_y(j)

def aPz(j: int) -> ndarray:
    """
    Alias for P_z(j).

    Parameters:
    j (int): Quantum number.

    Returns:
    ndarray: The z-component of the angular momentum operator.
    """
    return aP_z(j)

def wigner_small_d_exact(j, mp, m, beta):
    """Compute the exact Wigner small-d matrix element (Varshalovich et al. definition)."""
    prefactor = sqrt(factorial(j + mp) * factorial(j - mp) * factorial(j + m) * factorial(j - m))
    s_min = max(0, mp - m)
    s_max = min(j + mp, j - m)
    d = 0.0
    for s in range(s_min, s_max + 1):
        denom = factorial(j + mp - s) * factorial(s) * factorial(j - m - s) * factorial(m - mp + s)
        d += ((-1)**(s + m - mp) * prefactor / denom *
              (cos(beta / 2))**(2*j + mp - m - 2*s) *
              (sin(beta / 2))**(2*s + m - mp))
    return d

def wigner_D_matrix2(j, alpha, beta, gamma):
    dim = int(2*j + 1)
    D = zeros((dim, dim), dtype=complex)
    m_vals = arange(-j, j+1)
    for i, mp in enumerate(m_vals):
        for k, m in enumerate(m_vals):
            d = wigner_small_d_exact(j, mp, m, beta)
            D[i, k] = exp(-1j * mp * alpha) * d * exp(-1j * m * gamma)
    return D

def wigner_d_halfpi_high_precision2(j):
    dim = int(2*j + 1)
    beta = pi/2.0
    D = zeros((dim, dim), dtype=complex)
    m_vals = arange(-j, j+1)
    for i, mp in enumerate(m_vals):
        for k, m in enumerate(m_vals):
            d = wigner_small_d_exact(j, mp, m, beta)
            D[i, k] =  d 
    return D
def wigner_D_halfpi_matrix(j):
    dim = int(2*j + 1)
    D = zeros((dim, dim), dtype=complex)
    m_vals = arange(-j, j+1)
    dd = wigner_d_halfpi_high_precision(j)
    alpha, gamma = 0.0, 0.0
    for i, mp in enumerate(m_vals):
        for k, m in enumerate(m_vals):
            D[i, k] = exp(-1j * mp * alpha) * dd[i,k] * exp(-1j * m * gamma)
    return D.T

# either: 
# y-> z, x-> y, z-> x   by uising Ry(bet=pi/2) Rz(gam=pi/2)
# x-> z, y-> x, z-> y   by using Rx(-pi/2) Rz(-pi/2) in XYZ ->  Rz(alp=pi/2)Ry(bet=pi/2)Rz(gam=pi) for wigner needs to be in euler angles 
def wigner_swap_principal_axis(j,ax):
    dd = wigner_d_halfpi(j)
    if ax != 'z':
      m_vals = arange(-j, j+1)
      dim = int(2*j + 1)
      D = zeros((dim, dim), dtype=complex)
      for i, mp in enumerate(m_vals):
          for k, m in enumerate(m_vals):
             if ax == 'x':
                 D[i, k] = exp(-1j * mp * pi * 0.5 ) * dd[i,k] * exp(-1j * m * pi )
             else:
                 D[i, k] = exp(-1j * mp * pi * 0.0 ) * dd[i,k] * exp(-1j * m * pi * 0.5 )
      return D
    else:
      return eye(2*j+1)


