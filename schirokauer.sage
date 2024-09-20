# Implementation of Schirokauer maps for l = 2.

import field

from memoized import memoized

def coeffs_mod(K, a, l):
    num = vector(ZZ, vector(Zmod(l), a.numer.c))
    denom = ZZ(mod(a.denom, l))
    assert denom != 0, "division by zero"
    return K(tuple(num), denom)

def coeffs_div(K, a, l):
    num = vector(ZZ, a.numer.c) / l
    return K(tuple(num), a.denom)

@memoized
def ramification_index(d, l):
    K = field.field(d)
    if mod(K.discriminant(), l) != 0:
        return 1
    # multiquadratic fields are abelian, so the index is the same 
    # for all prime ideals in the decompostion
    e = K.sage().factor(l)[0][1]
    return e

@memoized
def residue_degree(d, l):
    K = field.field(d)
    # multiquadratic fields are abelian, so the residue degree is the
    # same for all prime ideals in the decompostion
    #f = K.sage().factor(l)[0][0].residue_field().degree()
    P = K.sage().factor(l)
    e = P[0][1]
    f = ZZ(K.degree() / (len(P)*e))
    return f

@memoized
def exponent(d, l, i):
    # Returns an integer eps such that x^eps = 1 mod l^i O_K for any non zero x in O_K\l O_K.
    #
    # Note: this value is not minimal in general.
    #
    # The value is determined by using Proposition 4.2.12 from
    # * Cohen H. - Advanced topics in computational number theory (2000)

    assert l == 2, "unimplemented for l != 2"

    f = residue_degree(d, l)
    e = ramification_index(d, l)
    #print(f"d = {d}, e = {e}, f = {f}")
    i *= e # l^i O_K = prod_j L_j^(i*e) for prime ideal L_j | l O_K

    eps = l^f - 1
    if i == 1:
        return eps

    if i == 2:
        eps *= l
        return eps
    
    if i == 3:
        if e == 1 and f==1:
            eps *= l
        else:
            eps *= l^2
        return eps
    
    if i == 4:
        if (e == 1 and f==1) or e >= 2:
            eps *= l^2
        else:
            eps *= l^3
        return eps

    # See [Coh2000, Ch.4, Exercise 20]
    if e >= l-1 and e < (l-1)*l^floor(log(i)/log(l)):
        s = ZZ(floor( i - l^floor(log(e/(l-1))/log(l)) )) + ZZ(floor( log(e/(l-1))/log(l) ))
        #print(f"s = {s}")
        eps *= l^s
    elif e >= (l-1)*l^floor(log(i)/log(l)):
        s = floor(log(i)/log(l))
        eps *= l^s
    else:
        raise Exception("Unimplemented!")
    return eps

@memoized
def mod_exp_1(d, l):
    e = ramification_index(d, l)
    f = residue_degree(d, l)
    ex = (e-1)*f
    return ex

@memoized
def mod_exp_2(d, l):
    ex = mod_exp_1(d, l)
    return ex-2

def schirokauer_map(d, a, l, i=1):
    K = field.field(d)
    if a == K.zero():
        #TODO: write tests for this case
        return vector(Zmod(l^i), [0]*2^(len(d)))

    t = walltime()
    eps = exponent(d, l, i)
    print(f"[schirokauer_map:exponent:{walltime(t)}", flush=True)
    
    #print(f"d = {d}, eps = {eps}, i = {i}")

    num = K(a.numer.c,1)
    denom = a.denom

    if denom == 1:
        num = coeffs_mod(K, num, l^(2*i))

    if mod(a.denom, l) == 0:
        #TODO: write tests for this case
        assert (a.absnorm(), l) != 0, "Schirokauer map is undefined"
        # for l = 2 an integral element of the field can have l in denominator,
        # since for d_i = 1 mod 4 we have O_K = Z[(1+sqrt(d_1))/2, (1+sqrt(d_2))/2, ...]
        num.denom = l^denom.valuation(l)
        denom = ZZ(denom / l^a.denom.valuation(l))

    denom = K.from_ZZ(denom)

    t=walltime()
    num = num^eps-1 # TODO: optimize by exponentiation mod l^(2 i)
    print(f"[schirokauer_map:eps:{walltime(t)}", flush=True)

    # OK = K.sage().ring_of_integers()
    # for j in range(1,10+1):
    #    print(f"a.num^eps-1 mod l^{j}*OK = {num.to_sage(K.sage()).mod(l^j*OK)}")

    # The element a.numer^eps-1 belongs to the ring Z[sqrt(d_1), ..., sqrt(d_n)].
    # This is only a subring of O_K, so we transform it with respect to O_K-basis.
    # Otherwise, it is possible that coefficients of a.numer^eps-1 are not divisible by l.
    t = walltime()
    X_a_num = vector(Zmod(l^(2*i)), pari.nfalgtobasis(K.sage().pari_nf(), num.to_sage(K.sage())))
    print(f"[schirokauer_map:modl2i:{walltime(t)}", flush=True)
    X_a_num = vector(ZZ, X_a_num)
    X_a_num /= l^i
    X_a_num = vector(Zmod(l^i), X_a_num)
    #print(f"X_a_num = {X_a_num}")

    denom = coeffs_mod(K, denom^eps-1, l^(2*i))
    X_a_denom = coeffs_div(K, denom, l^i).numer.c
    X_a_denom = vector(Zmod(l^i), X_a_denom)
    #print(f"X_a_denom = {X_a_denom}")
    #print(f"X_a_num - X_a_denom = {X_a_num - X_a_denom}")

    return X_a_num - X_a_denom
