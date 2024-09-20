# Arithmetic using representation of field elements by power products.

import field
import ideal
import fb
import units
import char
from memoized import memoized

@memoized
def field_compact(d):
  n = len(d)
  N = 2^n
  K = field.field(d)
  class KC:
    gens = d
    def __init__(f,*args):
      if len(args) == 1:
        g = args[0]
        if g.__class__ == K:
          f.elements = (g,)
          f.powers = (1,)
          return
        if n >= 1 and  g.__class__.gens == field_compact(d[:-1]).gens:
          f.elements = tuple(K(g.elements[i]) for i in range(len(g.elements))) 
          f.powers = tuple(g.powers)
          return
        if n >= 2 and g.__class__.gens == field_compact(d[:-2] + d[-1:]).gens:
          f.elements = tuple(K(g.elements[i]) for i in range(len(g.elements)))
          f.powers = tuple(g.powers)
          return
        if n >= 2 and g.__class__.gens == field_compact(d[:-2] + (d[-2]*d[-1],)).gens:
          f.elements = tuple(K(g.elements[i]) for i in range(len(g.elements)))
          f.powers = tuple(g.powers)
          return
      if len(args) == 2: # XXX: trusting caller to provide suitable values
        f.elements,f.powers = args
        return
      raise Exception('not known how to initialize field_compact(%s)(%s)' % (str(d),args))

    def __repr__(f):
      return 'field_compact.field_compact(%s)(%s,%s)' % (d,f.elements,f.powers)

    def __eq__(f,g):
      if g.__class__ != f.__class__: raise Exception('do not know how to compare %s to %s' % (f,g))
      f.remove_dup_elements(inplace=True)
      g.remove_dup_elements(inplace=True)
      if len(f.elements) != len(g.elements): # We assume that remove_dup_elements() joined duplicate entries.
        return False
      #if f.absnorm() != g.absnorm(): # FIXME: this can be expensive, make compact representation of norms and uncomment
      #  return False
      for i in range(len(f.elements)):
        if i >= len(g.elements):
          return False
        if f.elements[i] != g.elements[i]:
          try:
            j = g.elements.index(f.elements[i])
            if g.powers[j] != f.powers[i]:
              return False
          except ValueError:
            return False
        elif f.powers[i] != g.powers[i]:
          return False
      return True

    def __mul__(f,g):
      if f.__class__ != g.__class__:
        g = KC(g)
      helements = list(f.elements)
      hpowers = list(f.powers)
      for i in range(len(g.elements)):
        if i < len(f.elements) and g.elements[i] == f.elements[i]:
          hpowers[i] += g.powers[i]
        else:
          try:
            j = f.elements.index(g.elements[i])
            hpowers[j] += g.powers[i]
          except ValueError:
            helements.append(g.elements[i])
            hpowers.append(g.powers[i])
      h = KC(tuple(helements),tuple(hpowers))
      return h.trim()

    def __truediv__(f,g):
      g_inv_powers = [-g.powers[i] for i in range(len(g.powers))]
      g_inv = KC(g.elements, tuple(g_inv_powers))
      h = f * g_inv
      # TODO: make reduction/shortening
      return h

    def __pow__(f,e):
      if e == 0: return KC(K.one())
      if e == 1: return f
      g = vector(f.powers) * e
      return KC(tuple(f.elements), tuple(g))

    def sqrt(f, ideal_sqrt=False, compact_repr=False):
      '''
      Square root computation assuming that the powers in input are integers.
      
      When compact_repr=True the method uses compact representation trying to make the sizes elements in the result smaller.
      
      For ideal_sqrt=True the method tries to shorten the elements in the result by computing modulo units.
      '''
      
      helements = [0]*len(f.elements)
      hpowers = [0]*len(f.elements)
      el = K.one()
      for i in range(len(f.elements)):
        assert type(f.powers[i]) == Integer or (type(f.powers[i]) == Rational and f.powers[i].denominator()==1), f"Powers of type {f.powers[i].__class__} are not supported!" # call compute square roots first.
        e1,e0 = divmod(ZZ(f.powers[i]), 2) # f.powers[i] = e0 + 2*e1
        hpowers[i] = e1
        helements[i] = f.elements[i]
        if e0 == 1:
           # We assume that the length of f.elements is small as well as the sizes of its entries.
           # So, we perform slow compact representation computation after this loop.
           # This also may enable reduction due to cancellation.
           # The size of the el after the loop should be quasi-polynomial.
          el *= f.elements[i]
      h = KC(tuple(helements),tuple(hpowers))
      if compact_repr:
        elc = el.compact_repr(2)
        assert elc.powers[0].is_power_of(2)
        if elc.powers[0] == 1:
          el = elc.elements[0]
          if len(elc.elements) > 1:
            h *= KC(tuple(elc.elements[1:]), tuple(vector(elc.powers[1:]) / 2))
        else:
          el = K.one()
          h *= KC(tuple(elc.elements), tuple(vector(elc.powers) / 2))

      if ideal_sqrt:
        el_sqrt = ideal.idealsqrtshorten(len(d), 2^len(d), d, el)
      else:
        el_sqrt = el.sqrt()
      # there is no need to compute a compact representation of el_sqrt since it is already reduced.
      h = KC(el_sqrt) * h
      return h.trim(inplace=True, mod_units=ideal_sqrt)

    def conj(f): # conjugate last entry of d
      helements = [f.elements[i].conj() for i in range(len(f.elements))]
      hpowers = f.powers
      return KC(tuple(helements),tuple(hpowers))

    def symbols2_log(f,low,high):
      if len(f.elements) == 0:
        return 0
      r = vector(GF(2), f.elements[0].symbols2_log(low,high)) * ZZ(f.powers[0])
      for i in range(1,len(f.elements)):
        r += vector(GF(2), f.elements[i].symbols2_log(low,high)) *  ZZ(f.powers[i])
      return vector(ZZ, r)

    def symbols_log(f, low, high):
      if len(f.elements) == 0:
        return 0
      r = vector(GF(2), f.elements[0].symbols_log(low,high)) * f.powers[0]
      for i in range(1,len(f.elements)):
        r += vector(GF(2), f.elements[i].symbols_log(low,high)) *  f.powers[i]
      return vector(ZZ, r)

    def altsymbols_log(f, low, high, cheat=False):
      if len(f.elements) == 0:
        return 0
      K = field.field(d)
      f_s = []
      f_d = []
      qs = []
      for i in range(len(f.elements)):
        a = f.elements[i]
        assert a != K.zero()
        f_d += [a.denom]
        aq = a.get_qs(low,high, cheat=cheat)
        f_s += [[aq[j][0] for j in range(len(aq))]]
        if qs == []:
          qs = [aq[j][1] for j in range(len(aq))]
      
      aq = []
      aqd = []
      vec = []
      for i in range(len(f.powers)):
        if (f.powers[i] != 0):
           vec += [f.powers[i]]
           aq += [f_s[i]]
           aqd += [f_d[i]]
      vec = vector(vec)
      denom = vec.denominator()

      s = char.altfinal2(aq, aqd, qs, denom*vec, denom)
      if 0 in s:
        raise Exception('zero symbol')
      return tuple(0 if si == 1 else 1 for si in s)

    def bit_size(f):
      r = 0
      for i in range(len(f.elements)):
        r += f.elements[i].bit_size() + ZZ(f.powers[i]).nbits()
      return r

    def absnorm_old(f):
      r = 1
      s_r = 1
      for i in range(len(f.elements)):
        if f.powers[i].denominator() == 2:
          s_r *= f.elements[i].absnorm() ^ f.powers[i].numerator()
        else:
          r *= f.elements[i].absnorm() ^ f.powers[i]
      return r * s_r ^ (1/2)

    def absnorm(f, known_factors = []):
      #print(f"\nf.elements.absnorm() = {[f.elements[i].absnorm().factor() for i in range(len(f.elements))]}")
      #print(f"powers = {f.powers}")
      nrm = 1
      r = 1
      denom = 1
      l = 1
      for i in range(len(f.elements)):
        if f.powers[i].denominator() == 1:
          nrm *= f.elements[i].absnorm() ^ ZZ(f.powers[i])
        else:
          l = lcm(denom, ZZ(f.powers[i].denominator()))
          r = r^ZZ(l/denom) * f.elements[i].absnorm() ^ ZZ(f.powers[i].numerator() * l / f.powers[i].denominator())
          denom = l
      #print(f"r (not reduced) = {r.factor()}")
      # reduce part with non-integer powers using known factors of norm
      for i in range(len(known_factors)):
        v = r.valuation(known_factors[i])
        if v == 0:
          continue
        assert Mod(v, denom) == 0
        nrm *= ZZ(known_factors[i]) ^ ZZ(v / denom)
        r = r / ZZ(known_factors[i]) ^ ZZ(v)
        #print(f"r (new) = {r}, p = {known_factors[i]}, v = {v}, type(v) = {type(v)}, type(known_factors[i] ^ v) = {known_factors[i] ^ v}, type(p) = {type(known_factors[i])}, denom = {denom}")
        #print(f"r.factor() = {r.factor()}")
      #print(f"nrm = {nrm.factor()}")
      #print(f"r = {r.factor()}")
      #print(f"denom = {denom}")
      #print(f"QQ(nrm) * QQ(r^(1/denom)) = {(QQ(nrm) * QQ(r^(1/denom))).factor()}")
      #if r != 1:
      #  print(f"r (remaining) = {r}")

      res = QQ(nrm) * QQ(r^(1/denom))
      #assert f.absnorm_old() == res, f"{f.absnorm_old()} != {res}"
      #if f.absnorm_old() != res:
      #  print(f"{f.absnorm_old()} != {res}")
      #  exit(1)
      return res

    def norm_size(f):
      return abs(f.absnorm().numerator()) + abs(f.absnorm(denominator()))

    def approxlog(f, prec=None):
      if prec == None:
        prec = units.logbits
      R = RealField(prec)
      v = vector(R, f.elements[0].approxlog(prec=prec)) * f.powers[0]
      for i in range(1, len(f.elements)):
        v += f.powers[i] * vector(R, f.elements[i].approxlog(prec=prec))
      return v

    def lognorm(f, p=2, prec=None):
      if len(f.elements) == 0:
        return 0
      v = f.approxlog(prec=prec)
      return v.norm(p=p)
    
    def canonical_embedding_sage(f, prec=None):
      if prec == None:
        prec = units.logbits
      R = RealField(prec)
      v = vector(R, [1]*N)
      for i in range(len(f.elements)):
        u = f.elements[i].canonical_embedding_sage(prec=prec)
        u = vector(R, u)
        u = u.apply_map(lambda x: x^f.powers[i])
        for j in range(len(v)):
          v[j] *= u[j]
      return v
    
    def approxlog_sage(f, prec=None):
      if prec == None:
        prec = units.logbits
      R = RealField(prec)
      v = vector(R, f.elements[0].approxlog_sage(prec=prec)) * f.powers[0]
      for i in range(1, len(f.elements)):
        v += f.powers[i] * vector(R, f.elements[i].approxlog_sage(prec=prec))
      return v

    def lognorm_sage(f, p=2, prec=None):
      if len(f.elements) == 0:
        return 0
      v = f.approxlog_sage(prec=prec)
      return v.norm(p=p)

    def valuation(f, P):
      v = 0
      for el,e in zip(f.elements, f.powers):
        if e == 0:
          continue
        v += el.valuation(P) * e
      return v

    def factor(f, S): # S is a list of prime ideals
      h = KC(K.one())
      r = vector([0]*len(S))
      for el,e in zip(f.elements, f.powers):
        nrm = el.absnorm()
        el_v = vector([0]*len(S))
        for i in range(len(S)):
          el_v[i] = el.valuation(S[i]) * e
          nrm /= S[i].absolute_norm() ^ el_v[i]
        if abs(nrm) != 1:
          h *= el
        else:
          r += el_v
      return (h,tuple(r))

    def evaluate(f, ideal_sqrt=False, mod_units=False, rollback=True):
      # evaluates product. Very slow, it works only for small elements
      el = K.one()
      el_r = K.one()
      for i in range(len(f.elements)):
        if f.powers[i].denominator() == 1:
          el *= f.elements[i] ^ ZZ(f.powers[i])
        elif f.powers[i].denominator() == 2:
          el_r *= f.elements[i] ^ f.powers[i].numerator()
        else:
          raise NotImplementedError(f"Extraction of {f.powers[i].denominator()}-roots is not implemented!")
        if mod_units:
          el = el.reduce_mod_units(rollback=rollback)
          el_r = el_r.reduce_mod_units(rollback=rollback)
      if not ideal_sqrt:
        el_r = el_r.sqrt()
      else:
        el_r = ideal.idealsqrtshorten(len(d), 2^len(d), d, el_r)
      return el * el_r

    def compute_roots_old(f, ideal_sqrt=False): # expand rational exponents
      el = KC.one()
      el_r = KC.one()
      for i in range(len(f.elements)):
        if f.powers[i].denominator() == 1:
          el *= KC(f.elements[i]) ^ ZZ(f.powers[i])
        elif f.powers[i].denominator() == 2:
          if ideal_sqrt:
            try:
              rt = ideal.idealsqrtshorten(len(d), 2^len(d), d, f.elements[i])
              el *= KC(rt) ^ f.powers[i].numerator()
            except:
              el_r *= KC(f.elements[i]) ^ f.powers[i].numerator()
          else:
            try:
              rt = f.elements[i].sqrt()
              el *= KC(rt) ^ f.powers[i].numerator()
            except:
              el_r *= KC(f.elements[i]) ^ f.powers[i].numerator()
        else:
          raise NotImplementedError(f"Extraction of {f.powers[i].denominator()}-roots is not implemented!")
      el_r = el_r.sqrt(ideal_sqrt=ideal_sqrt)
      return el * el_r

    def compute_roots(f, ideal_sqrt=False, compact_repr=False):
      # expand rational exponents
      el = KC.one()
      el_r = KC.one()
      denom = 1
      for i in range(len(f.elements)):
        if f.powers[i].denominator() == 1:
          el *= KC(f.elements[i]) ^ ZZ(f.powers[i])
        elif Mod(f.powers[i].denominator(), 2) == 0:
          found = False
          if f.powers[i].denominator() == 2:
            if ideal_sqrt:
              try:
                rt = ideal.idealsqrtshorten(len(d), 2^len(d), d, f.elements[i])
                el *= KC(rt) ^ f.powers[i].numerator()
                found = True
              except:
                pass
            else:
              try:
                rt = f.elements[i].sqrt()
                el *= KC(rt) ^ f.powers[i].numerator()
                found = True
              except:
                pass
          if not found:
            l = lcm(denom, f.powers[i].denominator())
            el_r = el_r ^ (l / denom) * KC(f.elements[i]) ^ (f.powers[i].numerator() * l / f.powers[i].denominator())
            denom = l
        else:
          raise NotImplementedError(f"Extraction of {f.powers[i].denominator()}-roots is not implemented!")
      r = el
      if denom != 1:
        while Mod(denom, 2) == 0:
          el_r = el_r.sqrt(ideal_sqrt=ideal_sqrt, compact_repr=compact_repr)
          denom /= 2
      assert denom == 1, f"Extraction of {denom}-roots is not implemented!"
      r = r * el_r
      return r

    def to_sage(f, K = None, names='a', quotient=False):
      helements = tuple(f.elements[i].to_sage(K, names, quotient) for i in range(len(f.elements)))
      return (helements, f.powers)

    def trim(f, inplace=False, mod_units=False):
      helements = []
      hpowers = []
      for i in range(len(f.elements)):
        if f.elements[i] == K.one() or f.powers[i] == 0:
          continue
        if mod_units and f.elements[i].absnorm() in [1,-1]:
          ff = f.elements[i].reduce_mod_units(rollback=True)
          if ff not in [K.one(), -K.one()]: # TODO: check also full torsion for imaginary fields
            helements.append(ff)
            hpowers.append(f.powers[i])
          continue
        helements.append(f.elements[i])
        hpowers.append(f.powers[i])
      if len(helements) == 0:
        return KC.one()
      if inplace:
        f.elements = helements
        f.powers = hpowers
        return f
      return KC(tuple(helements), tuple(hpowers))

    def remove_dup_elements(f, inplace=False, mod_units=False):
      '''
      Joins duplicate elements. If mod_units = True it also try to join elements that differ by a unit.
      '''

      helements = list(f.elements)
      hpowers = list(f.powers)
  
      # TODO: make flag to mark already processed elements
      for i in range(len(f.elements)):
        if helements[i] == K.one():
          continue
        for j in range(i+1, len(helements)):
          if helements[i] == helements[j]:
            hpowers[i] += hpowers[j]
            helements[j] = K.one()
          elif mod_units and (helements[j].absnorm() / helements[i].absnorm() in [1,-1]):
            fji = (helements[j] / helements[i]).reduce_mod_units(rollback=True)
            # We can't just exclude element fji since there are elements h of K s.t. 
            # N(h) = +-1, h is not integral, and h is not a unit in O_K.
            if fji.bit_size() < helements[j].bit_size():
              hpowers[i] += hpowers[j]  
              helements[j] = fji

      if inplace:
        f.elements = tuple(helements)
        f.powers = tuple(hpowers)
        h = f
      else:
        h = KC(tuple(helements), tuple(hpowers))

      return h.trim(inplace=True, mod_units=mod_units)
    
    def shorten(f, prec=None):
      N = 2^len(d)
      print("-> approxlog ... ", flush=True)

      if min(d) > 0:
        v = f.approxlog(prec=prec)
      else:
        # OPTIMIZE: the method doesn't reduce well with fast approxlog for some reason
        v = f.approxlog_sage(prec=prec)
      print("-> units ... ", flush=True)
      S = units.generators_mod_torsion(d)
      print("-> unitloginverse ... ", flush=True)
      if min(d) > 0:
        # solving matrix equation is much faster if ideal.unitloginverse(d) is not precomputed
      #  e = v * ideal.unitloginverse(d)
        M = matrix(QQ,[(1,)*N] + [u.approxlog for u in S])
        e = M.solve_left(v)
      else:
        #M = matrix(QQ,[(1,)*N] + [u.approxlog for u in S])
        M = matrix(QQ,[(1,)*N] + [u.element.approxlog_sage(prec=prec) for u in S])
        e = M.solve_left(v)
      #e = [-ZZ(round(ej)) for ej in e[1:]]
      print("-> round ... ", flush=True)
      e = [-ZZ(ej.round()) for ej in e[1:]]
      u = KC(tuple([u.element for u in S]), tuple(e))
      print("-> mul ... ", flush=True) 
      r = f*u
      print("-> done ... ", flush=True) 
      return r

    def compactify_elements(f):
      '''
      Computes compact representation for the base elements.
      
      Transforms the product f = prod f_i^a_i to the form f = prod h_i^(a_i * 2^i), where h_i are small elements of the field.
      '''
      return prod([f.elements[i].compact_repr(2)^f.powers[i] for i in range(len(f.elements))], KC.one())

    
    def compactify(f, base=False):
      '''
      Computes compact representation described in [1, 4C].
      
      Transforms the product f = prod f_i^a_i to the form f = prod h_i^(2^i), where h_i are small elements of the field.
      
      Note 1: this procedure is slow and useable only for small degree fields.

      Note 2: if f is S-unit then it is possible that h_i is not an S-unit due to cancellation.

      [1] Biasse J.-F., Vredendaal v. C. - Fast multiquadratic S-unit computation and application to the calculation of class groups (2019)
      '''

      if base:
        f = f.compactify_elements()
      
      # processing of negative powers
      felements = list(f.elements)
      fpowers = list(f.powers)
      for i in range(len(felements)):
        if fpowers[i] < 0:
          felements[i] = K.one() / felements[i]
          fpowers[i] = abs(fpowers[i])
      f = KC(tuple(felements), tuple(fpowers))

      denom = 1
      for e in f.powers:
        denom = lcm(denom, QQ(e).denominator())
      assert denom.is_power_of(2) # unimplemented otherwise
      v = vector(f.powers) * denom
      v = vector(ZZ, v)

      # Computing the representation of f as f = prod f_i^(2^i).
      # We assume that the lengths of f.elements are small and the list f.elements consists of small elements.
      # So, the result of the following computations have quasi-polynomial size.
      s = {}
      for i in range(len(f.elements)):
        b = v[i].digits(2)
        for j in range(len(b)):
          if b[j] == 0:
            continue
          if 2^j in s:
            s[2^j] *= f.elements[i]
          else:
            s[2^j] = f.elements[i]
      
      #assert prod([g^e for e,g in s.items()], K.one()) == f.evaluate()

      # Now we compute a compact representation of quasi-polynomial size elements to make them of polynomial size.
      hpowers = []
      helements = []
      i = 0
      while i < floor(log(max(s.keys()))/log(2)+1):
        e = 2^i
        #print(f"s = {s}")
        if e not in s:
          i += 1
          continue
        g = s[2^i]
        gr = g.compact_repr(2)
        if gr.powers[0] == 1:
          helements.append(gr.elements[0])
          hpowers.append(e)
        else:
          e0 = e*gr.powers[0]
          if e0 in s:
            s[e0] *= gr.elements[0]
          else:
            s[e0] = gr.elements[0]
        for j in range(1,len(gr.elements)):
          e0 = e*gr.powers[j]
          if e0 in s:
            s[e0] *= gr.elements[j]
          else:
            s[e0] = gr.elements[j]
        i += 1

      hpowers = vector(hpowers) / denom
      return KC(tuple(helements), tuple(hpowers))

    def rdiv(f, g, bits=True, norms=False, log_units=False, log_sunits=False, l=2, prec=None, FB=None, mod_units=True):
      ''' 
      Given two elements f = prod f_i^a_i and g = prod g_i^b_i of the field divides elements f_i by g_j when 
      the size of f_i/g_j^sign(b_i) reduces with respect to given metric measuring a size of elements.
      
      The "size" is measured using the following metrics:
      1. ||ApproxLog(.)||_l, norm of logarithmic complex embeddings of the element (log_units = True).
      2. ||Log_S(.)||_l, norm of logarithmic S-embeddings of the element. Requires factor base (FB) to compute valuations (log_sunits = True, FB != None).
      3. bit_size(.), bit size of element's internal representation (bits = True).
      4. N(.), absolute norm of the element (norms=True).
      
      If mod_units = True the method reduces all intermediate results modulo units.
      
      Returns triple (f/s, g/s, s) where s s.t. size(f/s) <= size(f) and s divides (exactly) g.
      '''

      assert bits or norms or log_units or (log_sunits and FB!=None), "enable at least one metric"

      h = KC.one()
      s = KC.one()
      gpowers = list(g.powers)

      for i in range(len(f.elements)):
        el = f.elements[i]

        if vector(gpowers).is_zero():
          break

        j = 0
        while j < len(g.elements):
          el_old = el
          if gpowers[j] > 0:
            if log_units:
              l1 = vector(el.approxlog(prec=prec))
              l2 = vector(g.elements[j].approxlog(prec=prec))
              if RR((l1-l2).norm(l)) >= RR(l1.norm(l)):
                j += 1
                continue
            if norms:
              n1 = el.absnorm()
              n2 = g.absnorm()
              if n1/n2 >= n1:
                j += 1
                continue
            if log_sunits:
              v1 = vector([el.valuation(FB[k]) for k in range(len(FB))])
              v2 = vector([g.elements[j].valuation(FB[k]) for k in range(len(FB))])
              if (v1-v2).norm(l) >= v1.norm(l):
                j += 1
                continue
            el = el / g.elements[j]
            if mod_units:
              el = el.reduce_mod_units(rollback=True)
            if bits and (el.bit_size() >= el_old.bit_size()):
              el = el_old
              j += 1
              continue
            gpowers[j] -= f.powers[i]
            s *= KC(g.elements[j]) ^ f.powers[i] 
            j = 0
          elif gpowers[j] < 0:
            if log_units:
              l1 = vector(el.approxlog(prec=prec))
              l2 = vector(g.elements[j].approxlog(prec=prec))
              if RR((l1+l2).norm(l)) >= RR(l1.norm(l)):
                j += 1
                continue
            if norms:
              n1 = el.absnorm()
              n2 = g.absnorm()
              if n1*n2 >= n1:
                j += 1
                continue
            if log_sunits:
              v1 = vector([el.valuation(FB[k]) for k in range(len(FB))])
              v2 = vector([g.elements[j].valuation(FB[k]) for k in range(len(FB))])
              if RR((v1+v2).norm(l)) >= RR(v1.norm(l)):
                j += 1
                continue
            el = el * g.elements[j]
            if mod_units:
              el = el.reduce_mod_units(rollback=True)
            if bits and (el.bit_size() >= el_old.bit_size()):
              el = el_old
              j += 1
              continue
            s /= KC(g.elements[j]) ^ f.powers[i]
            gpowers[j] += f.powers[i]
            j = 0
          else:
            j += 1

        h *= KC(el)^ZZ(f.powers[i])

      h = h.trim(mod_units=mod_units)
      r = KC(g.elements, tuple(gpowers)).trim(mod_units=mod_units)
      return (h, r, s)

    def pseudo_gcd(f, g):
      f.remove_dup_elements(inplace=True)
      g.remove_dup_elements(inplace=True)
      r = KC.one()
      for i in range(len(f.elements)):
        if f.powers[i] == 0:
          continue
        for j in range(len(g.elements)):
          if g.powers[j] == 0:
            continue
          if f.elements[i] == g.elements[j]:
          #if abs(f.elements[i].absnorm() / g.elements[j].absnorm()) == 1:
            r *= KC(f.elements[i]) ^ min(f.powers[i], g.powers[j])
            continue
      return r.trim()

    def reduce_mod_units(f, rollback=True):
      # TODO: try to reduce all elements at once, not just per element
      helements = []
      hpowers = []
      for i in range(len(f.elements)):
#        if f.elements[i].absnorm() in [1,-1]: # doesn't work, since we need to check that element is integral
#          continue
        el = f.elements[i].reduce_mod_units(rollback=rollback)
        if rollback and (el.bit_size() >= f.elements[i].bit_size()):
          helements.append(f.elements[i])
          hpowers.append(f.powers[i])
        else:
          helements.append(el)
          hpowers.append(f.powers[i])
      return KC(tuple(helements), tuple(hpowers)).remove_dup_elements(inplace=True, mod_units=True)

    def reduce_lin_dep(f, logbits=None, mod_units=False, mod_sunits=False, SU=None, rollback=False):

      if mod_sunits:
        assert SU != None
      #M = units.shortening_matrix(len(d), 2^len(d), f.elements)

      #U = units.generators_mod_torsion(d)

      # shortening matrix that is used to reduce vectors. Taken from units.shortening_matrix
      S = f.elements
      #S = f.elements + [u.element for u in U]

      if logbits==None:
        logbits = units.logbits*2

      #c = 2^(logbits+n-1)
      #c = 2^(logbits+len(S)-1)
      #M = matrix(ZZ,[[ZZ(i==j) for i in range(len(S))] + [ZZ(sign(f.powers[j]) * v*2^(logbits+n-1)) for v in S[j].approxlog()] for j in range(len(S))])
      #M = matrix(ZZ,[[ZZ(i==j) for i in range(len(S))] + [ZZ(v*2^(logbits+n-1)) for v in S[j].approxlog()] for j in range(len(S))])
      #print("ApproxLogs:")
      #print([S[i].approxlog() for i in range(len(S))])

      #M = matrix(ZZ,[[ZZ(i==j) for i in range(len(S))] + [ZZ(v*c) for v in S[j].approxlog()] for j in range(len(S))])
      #M = matrix(QQ,[[QQ(i==j) for i in range(len(S))] + [ZZ(v*c) for v in S[j].approxlog()] for j in range(len(S))])
      #M = matrix(ZZ,[[ZZ(i==j) for i in range(len(S))] + [ZZ(round(v*c)) for v in S[j].approxlog()] for j in range(len(S))])
      
      t = walltime()
      
      #M = matrix(QQ,[[QQ(i==j) for i in range(len(S))] + [QQ(sign(f.powers[j])*v) for v in S[j].approxlog(prec=logbits)] for j in range(len(S))])
      M = matrix(QQ,[[QQ(i==j) for i in range(len(S))] + [QQ(v) for v in S[j].approxlog(prec=logbits)] for j in range(len(S))])
      
      #print(f"M = {M}")
      M = M.LLL()

      #print(f"M (without zeroes) = {[Mi[:len(S)] for Mi in M if not Mi[len(S):].is_zero()]}")

      H = [Mi[len(S):] for Mi in M]

      # obtaing transformation matrix
      #T = [Mi[:len(S)] for Mi in M if not Mi[len(S):].is_zero()]
      T = [Mi[:len(S)] for Mi in M]

      #H,T = M.LLL(transformation=True)
      #H,T = M.hermite_form(transformation=True)

      T = matrix(ZZ, T)

      print(f"T (transform. matrix):\n{T}")

      MM = matrix(QQ,[[QQ(v) for v in S[j].approxlog(prec=logbits)] for j in range(len(S))])
      #MM = matrix(QQ,[[QQ(sign(f.powers[j])*v) for v in S[j].approxlog(prec=logbits)] for j in range(len(S))])
      assert T * MM == matrix(H), f"M*MM = {T * MM}, H = {H}"

      v = vector(f.powers)
      #v = vector([abs(f.powers[i]) for i in range(len(f.powers))])

      T_inv = T^(-1)
      print(f"T^(-1) = {T_inv}")

      powers = tuple(v * T_inv)

      assert len(powers) == T.nrows()

      if mod_sunits:
        s = SU.one()
      
      r = KC.one()

      #helements = []
      #hpowers = []
      for i in range(T.nrows()):
        if powers[i] == 0:
          continue
        #if H[i].is_zero():
        #  continue
        if mod_sunits:
          s_el = SU.one()
        el = K.one()
        for j in range(T.ncols()):
          if mod_sunits:
            ex,s0 = S[j].pow_mod_sunits(T[i,j], SU)
            s_el *= s0
            el,s0 = el.mul_mod_sunits(ex, SU, mod_units=mod_units)
            s_el *= s0
          elif mod_units:
            el *= (S[j] ^ T[i,j]).reduce_mod_units(rollback=True)
          else:
            el *= S[j] ^ T[i,j]

        #helements.append(el)
        #hpowers.append(powers[i])
        r *= KC(el) ^ powers[i]

        # We expect that all values in r.elements are independent, after LLL-reduction.
        # So, bit sizes can't decrease when we add new element el. Thus we perform early rollback.
        if rollback and r.bit_size() >= f.bit_size():
          print(f"[reduce_lin_dep(rollback)|{walltime(t)} sec.]")
          if mod_sunits:
            return (f, SU.one())
          return f
        
        if mod_sunits:
          s *= s_el ^ powers[i]
        if H[i].is_zero():
          #assert el == K.one(), f"{el.to_sage(K.sage())} != 1"
          print(f"H_i = 0, el = {el.to_sage(K.sage())}")

      #assert len(hpowers) == len(helements)

      #r = KC(tuple(helements), tuple(hpowers)).trim(inplace=True, mod_units=mod_units)

      print(f"[reduce_lin_dep(elements: {len(S)} => {len(r.elements)}, bits: {f.bit_size()} => {r.bit_size()}])|{walltime(t)} sec.]")

      #if abs(r.absnorm()) != abs(f.absnorm()):
      # if attempts <= 0:
      #   print(f"[info] Failed to reduce")
      #   return f
      # print(f"[info] Norms of reduced and the original element are not equal ({r.absnorm()} != {f.absnorm()})! Not enough precision. Trying again with logbits = {logbits*2}")
      # return f.reduce_LLL(logbits=logbits*2, attempts=attempts-1)

      #assert r.absnorm() == f.absnorm(), f"{r.absnorm()} != {f.absnorm()}"
      #assert r.evaluate() == f.evaluate(), f"{r.evaluate().to_sage(K.sage())} != {r.evaluate().to_sage(K.sage())}"
      if mod_sunits:
        return (r,s)
      return r

    def reduce_norms(f, times=1000):
      helements = list(f.elements)
      hpowers = list(f.powers)
      for k in range(times):
        for i in range(len(helements)):
          if helements[i].absnorm() in [1,-1]:
            continue
          nrm_i = helements[i].absnorm()
          for j in range(i+1, len(helements)):
            nrm_j = helements[j].absnorm()
            # if helements[j].absnorm() in [1,-1]:
            #   continue

            # if nrm_i / nrm_j in [1,-1]:
            #   hpowers[i] += hpowers[j]
            #   helements[j] = K.one()
            #   hpowers[j] = 0
            #   continue

            # if nrm_i * nrm_j in [1,-1]:
            #   hpowers[i] -= hpowers[j]
            #   helements[j] = K.one()
            #   hpowers[j] = 0
            #   continue

            nrm1 = nrm_i * nrm_j
            if abs(nrm_i).numerator() + abs(nrm_i).denominator() > abs(nrm1).numerator() + abs(nrm1).denominator():
              helements[i] = (helements[i] * helements[j]).reduce_mod_units()
              #helements[i] = (helements[i] * helements[j])
              hpowers[j] = hpowers[j] - hpowers[i]
              continue
            
            nrm2 = nrm_i / nrm_j
            if abs(nrm_i).numerator() + abs(nrm_i).denominator() > abs(nrm2).numerator() + abs(nrm2).denominator():
              helements[i] = (helements[i] / helements[j]).reduce_mod_units() # slow, since division uses Sage artithmetic.
              #helements[i] = (helements[i] / helements[j]) # slow, since division uses Sage artithmetic.
              hpowers[j] = hpowers[j] + hpowers[i]
              continue

            nrm3 = nrm_j / nrm_i
            if abs(nrm_j).numerator() + abs(nrm_j).denominator() > abs(nrm3).numerator() + abs(nrm3).denominator():
              helements[j] = (helements[j] / helements[i]).reduce_mod_units() # slow, since division uses Sage artithmetic.
              #helements[j] = (helements[j] / helements[i]) # slow, since division uses Sage artithmetic.
              hpowers[i] = hpowers[i] + hpowers[j]
              continue 
            
            nrm4 = nrm_j * nrm_i
            if abs(nrm_j).numerator() + abs(nrm_j).denominator() > abs(nrm4).numerator() + abs(nrm4).denominator():
              helements[j] = (helements[j] * helements[i]).reduce_mod_units() # slow, since division uses Sage artithmetic.
              #helements[j] = (helements[j] * helements[i]) # slow, since division uses Sage artithmetic.
              hpowers[i] = hpowers[i] - hpowers[j]
              continue
      #return KC(helements, tuple(hpowers)).trim(mod_units=True)
      return KC(helements, tuple(hpowers)).trim()

    def reduce_base_mod_sunits(f, SU, mod_units=True):
      '''
      Given element f = prod_i f_i^e_i of the field K reduces each component f_i modulo S-units.

      Returns pair (f/s, s) where s is an S-unit.
      '''
      h = KC.one()
      s = KC.one()
      for i in range(len(f.elements)):
        hi,si = f.elements[i].reduce_mod_sunits(SU, mod_units=mod_units)
        h *= KC(hi) ^ f.powers[i]
        s *= si ^ f.powers[i]
      return (h,s)

    def reduce_all_mod_sunits(f, SU, mod_units=True):
      '''
      Given element f = prod_i f_i^e_i of the field K reduces (mod S-units) all product at once.

      Returns pair (f/s, s) where s is an S-unit and size(f/s) <= size(s).
      Where size is the bit_size of elements.
      '''

      A = SU.relations_matrix()
      FB = SU.factor_base(sage=True)

      v = vector([f.valuation(FB[i]) for i in range(len(FB))])

      # f.absnorm() can be slow for non-reduced f.
      # if fb.is_smooth(f.absnorm(), FB0): 
      #   return (K.one(), SU.from_prime_prod(v))

      w = clgp.approx_CVP(d, list(v))
      assert len(w) == len(SU.gens())

      wA = vector(w) * A

      # reduction failed, resulting vector has larger norm
      if RR((vector(v)-wA).norm()) >= RR(vector(v).norm()):
        return (f, SU.one())
 
      s = prod([SU.gens()[i]^w[i] for i in range(len(w))], KC.one()).compute_roots(ideal_sqrt=mod_units)

      h,r,s = f.rdiv(s, bits=False, log_sunits=True, FB = FB)

      # TODO: reduce linear dependency using LLL, otherwise this is not very usefull procedure.

      return (h,s)

    def reduce_mod_sunits(f, SU, mod_units=True):
      '''
      Reduces element f = prod f_i^e_i modulo S-units in the following way:
      1) reduces base elements f_i individualy to remove S-unit components.
      2) reduces base elements f_j (j=1, ...) all together to remove dependence by S-unit between them.

      Returns pair (f/s, s) where s is an S-unit and size(f/s) <= size(s).
      Where size is the bit_size of elements.
      '''

      h,s = f.reduce_base_mod_sunits(SU, mod_units=mod_units)
      #h,s0 = h.reduce_all_mod_sunits(SU, mod_units=mod_units) # rdiv doesn't provide good reduction in this case
      #s *= s0

      # TODO: other reductions and verbose mode
      return (h,s)

    @staticmethod
    def random(bits = 8, exp_bound = 5, elems = 4, exp_neg = False):
      felements = tuple(K.random(bits) for i in range(elems))
      if exp_neg:
        fpowers = tuple(ZZ.random_element(-exp_bound,exp_bound) for j in range(elems))
      else:
        fpowers = tuple(ZZ.random_element(1,exp_bound) for j in range(elems))
      return KC(felements, fpowers)
    @staticmethod
    def base_field():
      return K
    @staticmethod
    def one():
      return KC(K.one())
    @staticmethod
    def zero():
      return KC(K.zero())
    @staticmethod
    def from_ZZ(a):
      return KC((K.from_ZZ(a),),(1,))

  return KC
