from fpylll.algorithms.bkz2 import BKZReduction
from fpylll import BKZ as BKZ_FPYLLL, GSO, IntegerMatrix

from memoized import memoized

def mat_floor(B, precision):
  f(x) = floor(x * 2^precision)
  return matrix(ZZ, B.apply_map(f))

def BKZ_reduce(A, precision = "mpfr", transformation = False, block_size=40):
  A = IntegerMatrix.from_matrix(A)
  M = IntegerMatrix.identity(A.nrows)

  gso = GSO.Mat(A, float_type=precision,
                U=IntegerMatrix.identity(M.nrows, int_type=M.int_type),
                UinvT=IntegerMatrix.identity(M.nrows, int_type=M.int_type))
  gso.update_gso()

  bkz = BKZReduction(gso)

  for b in range(7, block_size + 1):
    par = BKZ_FPYLLL.Param(b,
		  strategies=BKZ_FPYLLL.DEFAULT_STRATEGY,
		  max_loops=8,
		  flags=BKZ_FPYLLL.MAX_LOOPS
	  )
    bkz(par)

  H = gso.B.to_matrix(matrix(ZZ, gso.B.nrows, gso.B.ncols))
  M = gso.U.to_matrix(matrix(ZZ, gso.U.nrows, gso.U.ncols))
  if transformation:
    return H,M
  else:
    return H

def BKZ_reduce_approx(B, prec, fplll_precision = "mpfr", transformation = False, block_size=40):
  A = mat_floor(B, prec)
  A = BKZ_reduce(A, precision=fplll_precision, transformation = transformation, block_size=block_size)
  return matrix(B.base_ring(), A / 2^prec)

def isometry(r1, r2):
  '''
  Matrix that represents an isometry used to make the log-unit lattice full-rank. Taken from [BR20, $3.1].
  
  [BR20] Bernard O., Roux-Langlois A. - Twisted-PHS - Using the product formula to Solve Approx-SVP in Ideal Lattices (2020).
  '''
  r1 = ZZ(r1)
  r2 = ZZ(r2)
  v = r1 + r2 - 1
  M1 = matrix.circulant([-1,1]+[0]*(v-1))[:-1]
  M1 = matrix(QQ, M1)
  #print(f"M1:\n{M1}")
  assert M1.nrows() == v
  assert M1.ncols() == v+1
  M2 = matrix(QQ, [1/2,1/2])
  M2 = block_diagonal_matrix([M2]*r2)
  M2 = matrix(QQ, M2)
  #print(f"M2:\n{M2}")
  assert M2.nrows() == r2, f"{M2.nrows()} != {r2}"
  assert M2.ncols() == 2*r2
  I1 = identity_matrix(QQ, r1)
  M3 = block_diagonal_matrix([I1, M2])
  M3 = matrix(QQ, M3)
  #print(f"M3 = {M3}")
  assert M3.nrows() == r1 + r2
  assert M3.ncols() == r1 + 2*r2
  M = M1 * M3
  return matrix(QQbar, M).gram_schmidt(orthonormal=True)[0]

def dual(B):
  #return (B^(-1)).transpose() # for full rank matrices only
  B_t = B.transpose()
  B_dual = B_t * (B * B_t)^(-1)
  return B_dual.transpose()

@memoized
def log_unit(U, fld=QQ):
  '''
  Basis matrix of the log-unit lattice.
  '''
  if len(U) == 0:
    return matrix(ZZ)
  N = len(U[0].approxlog)
  logU = matrix(fld, [[u.approxlog[i] for i in range(N)] for u in U])
  return logU

@memoized
def log_unit_ZZ(U, precision, fld=QQ):
  '''
  Computes basis matrix of a lattice floor(Log(U)*2^precision).
  '''
  if len(U) == 0:
    return matrix(ZZ)
  logU = log_unit(U, fld)
  return mat_floor(logU, precision)

@memoized
def log_unit_isom(U, r1, r2, fld=QQ):
  '''
  Computes a basis matrix of a full-rank lattice isometric to log unit lattice.
  '''
  logU = log_unit(U, fld)
  M = matrix(fld, isometry(r1, r2))
  logU_ism = logU * M.transpose()
  return logU_ism

@memoized
def log_unit_isom_ZZ(U, r1, r2, precision, fld=QQ):
  '''
  Computes a basis matrix of a full-rank lattice isometric to log unit lattice.
  The result of computation is converted to integers.
  '''
  logU_isom = log_unit_isom(U, r1, r2, fld)
  return mat_floor(logU_isom, precision)

@memoized
def log_unit_dual(U, fld=QQ):
  '''
  Computes a basis matrix of the dual lattice of log unit lattice.
  '''
  logU = log_unit(U, fld)
  return dual(logU)

@memoized
def log_unit_isom_dual(U, r1, r2, fld=QQ):
  '''
  Computes a basis matrix of the lattice dual to a full-rank lattice isometric to the lattice of log unit lattice.
  '''
  logU_isom = log_unit_isom(U, r1, r2, fld)
  return dual(logU_isom)

@memoized
def log_unit_isom_dual_ZZ(U, r1, r2, precision, fld=QQ):
  '''
  Computes a basis matrix of a full-rank lattice isometric to dual lattice of the log unit lattice.
  The result of computation is converted to integers.
  '''
  logU_dual_isom = log_unit_isom_dual(U, r1, r2, fld)
  return mat_floor(logU_dual_isom, precision)