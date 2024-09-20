import goodprime
import goodprimecheat
import centermod
import mult
import units

from memoized import memoized

def evaluate(n,N,q,s,f):
  if n == 0: return f[0]
  n -= 1
  N /= 2
  g = [(f[j] + s[n] * f[j + N]) % q for j in range(N)]
  return evaluate(n,N,q,s,g)
  
def symbol2(d,f,denom,t):
  n = len(d)
  N = 2^n
  q,s,sinv = goodprimecheat.sqrtprime(d,t)
  v = evaluate(n,N,q,s,f)
  return v.jacobi(q) / denom.jacobi(q)

def symbols2(d,f,denom,low,high):
  n = len(d)
  N = 2^n
  if high < low + 1:
    raise Exception('%s must be below %s' % (low,high))
  if high == low + 1:
    return (symbol2(d,f,denom,low),)
  mid = low + 1
  while mid + (mid - low) < high:
    mid = mid + (mid - low)
  q0 = goodprimecheat.sqrtprime_product(d,low,mid)
  q1 = goodprimecheat.sqrtprime_product(d,mid,high)
  f0 = centermod.vector(f,q0)
  denom0 = denom % q0
  f1 = centermod.vector(f,q1)
  denom1 = denom % q1
  return symbols2(d,f0,denom0,low,mid) + symbols2(d,f1,denom1,mid,high)

def symbol(d,f,denom,t):
  n = len(d)
  N = 2^n
  q,s,sinv = goodprime.sqrtprime(d,t)
  v = evaluate(n,N,q,s,f)
  #print "symbol",t, v, q, denom
  return v.jacobi(q) / denom.jacobi(q)

def symbols(d,f,denom,low,high):
  n = len(d)
  N = 2^n
  if high < low + 1:
    raise Exception('%s must be below %s' % (low,high))
  if high == low + 1:
    return (symbol(d,f,denom,low),)
  mid = low + 1
  while mid + (mid - low) < high:
    mid = mid + (mid - low)
  q0 = goodprime.sqrtprime_product(d,low,mid)
  q1 = goodprime.sqrtprime_product(d,mid,high)
  f0 = centermod.vector(f,q0)
  denom0 = denom % q0
  f1 = centermod.vector(f,q1)
  denom1 = denom % q1
  return symbols(d,f0,denom0,low,mid) + symbols(d,f1,denom1,mid,high)

def is_2nth_root(a, denom, q):
  res = None
  denom0 = denom
  if mod(q,2*denom) == 1: # fast check
    while denom0 > 1:
      denom0 /= 2
      # for q = 1 mod (2*denom) the sqrt is always computable
      a = mod(a, q).sqrt()
    res = ZZ(a).jacobi(q)
  else:
    if mod(a, q) == 0:
      res = 0
    else:
      try:
        rt = mod(a, q).nth_root(2*denom)
        res = 1
      except: #ValueError:
        res = -1
  return ZZ(res)

# Evaluation of characters corresponding to primes in qs at power products of small field elements (compact representation of a field element).
#
# We assume that the list of base elements is separated into two parts in such way that second part is sigma-conjugates of the first part.
#
# Input: 
# -- Aq is the (precomputed) result of evaluation of characters at the first part of base field elements (Aqd is the same for denominators of the first part).
# -- Aqc is the (precomputed) result of evaluation of characters at the second part of base field elements (Aqcd is the same for denominators of the second part).
# -- qs list of primes corresponding to characters.
# -- vec is the integer vector of exponents at first part of base elements (for a field element in compact representation).
# -- vecc is the integer vector of exponents at second part of base elements (for a field element in compact representation).
# -- denom is a common denominator for exponents in vec and vecc (for support of rational exponents).
def altfinal(Aq, Aqd, Aqc, Aqcd, qs, vec, vecc, denom):
  chi = []
  for i in range(len(qs)):
    #denom0 = denom
    q = qs[i]
    k = GF(q)
    newq = prod([ZZ(k((Aq[v][i]/Aqd[v]))).powermod(vec[v],q) for v in range(len(vec))]) % q
    newq *= prod([ZZ(k((Aqc[v][i]/Aqcd[v]))).powermod(vecc[v],q) for v in range(len(vecc))]) 
    newq = newq % q
    chi += [is_2nth_root(newq, denom, q)]
  return tuple(chi)

# The method is for the case when we have only one list of base field elements (e.g. we do not split into two parts).
def altfinal2(Aq, Aqd, qs, vec, denom):
  chi = []
  for i in range(len(qs)):
    #denom0 = denom
    q = qs[i]
    k = GF(q)
    newq = prod([ZZ(k((Aq[v][i]/Aqd[v]))).powermod(vec[v],q) for v in range(len(vec))]) % q
    chi += [is_2nth_root(newq, denom, q)]
    #print(f"newq = {newq}, q = {q}, chi = {chi[-1]}")
  return tuple(chi)

def altfinalzeta(Aq, Aqd, Aqc, Aqcd, qs, vec, vecc, denom):
  chi = []
  for i in range(len(qs)):
    #denom0 = denom
    q = qs[i]
    k = GF(q)
    newq = prod([ZZ(k((Aq[v][i]/Aqd[v]))).powermod(vec[v],q) for v in range(len(vec))]) % q
    newq *= prod([ZZ(k((Aqc[v][i]/Aqcd[v]))).powermod(vecc[v],q) for v in range(len(vecc))]) 
    newq = newq % q
    newq = mod(newq,q)^((q-1)/(2*denom))
    chi += [log(newq, primroot(q, 2*denom))]
  return vector(Zmod(2*denom), chi)

def altsymbol(d, f, denom, t, cheat=False):
  n = len(d)
  N = 2^n
  if not cheat:
    q,s,sinv = goodprime.sqrtprime(d,t)
  else:
    q,s,sinv = goodprimecheat.sqrtprime(d,t)
  v = evaluate(n,N,q,s,f)
  return v, q

def altsymbols(d, f, denom, low, high, cheat=False):
  n = len(d)
  N = 2^n
  if high < low + 1:
    raise Exception('%s must be below %s' % (low,high))
  if high == low + 1:
    return (altsymbol(d,f,denom,low,cheat=cheat),)
  mid = low + 1
  while mid + (mid - low) < high:
    mid = mid + (mid - low)
  if not cheat:
    q0 = goodprime.sqrtprime_product(d,low,mid)
    q1 = goodprime.sqrtprime_product(d,mid,high)
  else:
    q0 = goodprimecheat.sqrtprime_product(d,low,mid)
    q1 = goodprimecheat.sqrtprime_product(d,mid,high)
	#
	# centermod is defined in centermod.pyx
	# centermod.vector(v,q) takes the coos of vector v mod q in [-q/2, q/2) (not sure about the end points)
	#
  f0 = centermod.vector(f,q0)
  denom0 = denom % q0
  f1 = centermod.vector(f,q1)
  denom1 = denom % q1
  return altsymbols(d,f0,denom0,low,mid,cheat=cheat) + altsymbols(d,f1,denom1,mid,high,cheat=cheat)

@memoized
def primroot(q, r):
  return GF(q).primitive_element()^((q-1)/r)
