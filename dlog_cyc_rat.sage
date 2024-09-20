# Computing discrete logarithm using computation of square root in cyclic subgroups of CL_K decompose I using decomposition of I^2.

import field
import ideal
import relations
import fb
import trees
import verify
import clgp
import dlog_quad
import ideal_sqrt_cyc
import smooth_ideal_sqrt_sat
import field_compact
import idealprod
import units
import sunits
from memoized import memoized
import argparse

from fpylll import GSO, IntegerMatrix, LLL

import cProfile
import pstats
from pstats import SortKey

profile = cProfile.Profile()

food = trees.get_food() # load food from trees

parser = argparse.ArgumentParser(description='Computing class group discrete logarithm for a random ideal.')
parser.add_argument("--quadchars", dest="quadchars", action=argparse.BooleanOptionalAction, default=False, help="Use quadratic characters in ideal square root to avoid obstructions.")

args = parser.parse_args()

verify.set(verify.NONE)

proof.all(False)

# Folder to save results of DLOG computations for the field and its subfields. Set it to None if you don't want to store results
SAVE_ALL_RESULTS = "experiments/dlogs"

pari.allocatemem(600*1024^3)
print(f"Pari stack size: {pari.stacksize()}")

def apply_aut_rec(K, K_gens, f, mu, i):
    n = len(K_gens)
    if i >= n:
        return f # f in QQ
    res = K.zero()
    g = K(K_gens[i])
    c = f.list()
    res = apply_aut_rec(K, K_gens, c[0], mu, i+1) + apply_aut_rec(K, K_gens, c[1], mu, i+1) * (-1)^mu[i] * g
    return res

def apply_aut(K, f, mu):
    K_gens = K.gens()
    return apply_aut_rec(K, K_gens, f, mu, 0)

# Heuristic check for equality I == g * prod P_i^(e_i)], i = 1, ..., #FB.
def check_sqrt(K, FB, I, e, g):
    nrm = [K.ideal(FB[i].prime, K(FB[i].elts).absolute_norm()^e[i]) for i in range(len(e))]
    if I.q != g.absnorm() * nrm:
        return False

    g = g.to_sage(K)
    #return True
    for i in range(len(e)):
        # TODO: add check for e[i] < 0
        if e[i] > 0:
            P = K.ideal(FB[i].prime, K(FB[i].elts))
            el = I.random_element()
            if not (el/g in P):
            #if not (el in g*P):
                print(f"-> check for e[{i}] = {e[i]} [fail]")    
                return False
            else:
                print(f"-> check for e[{i}] = {e[i]} [ok]") 
        elif e[i] < 0:
            P = K.ideal(FB[i].prime, K(FB[i].elts))
            el = I.random_element()
            P_inv = clgp.prime_ideal_inverse(P)
            if not (el/g in P_inv):
            #if not (el in g*P_inv):
                print(f"-> check for e[{i}] = {e[i]} [fail]")    
                return False
            else:
                print(f"-> check for e[{i}] = {e[i]} [ok]")
    return True

# Heuristic check for equality I == g * prod P_i^(e_i)], i = 1, ..., #FB.
def ideal_eq_test(K, FB, I, e, g):
    assert len(FB) == len(e), f"Wrong exponent vector e, len(e) = {len(e)} != {len(FB)}"
    known_factors = [I.q] + [FB[i].prime for i in range(len(FB))]
    nrm = prod([K.ideal(FB[i].prime, K(FB[i].elts)).absolute_norm()^ZZ(e[i]) for i in range(len(e))])
    if abs(I.q) != abs(g.absnorm(known_factors=known_factors) * nrm): # compare upto units
        print("\t[fail (norms)]")
        return False
    #else:
    #    print("\t[ok]")
    #    return True

    I_sage = I.to_sage(K)
    for i in range(len(e)):
        P = K.ideal(FB[i].prime, K(FB[i].elts))
        if g.valuation(P) + e[i] != I_sage.valuation(P):
            print("\t[fail (valuations)]")
            return False
    return True

# Heuristic check for equality I^2 == g * prod P_i^(e_i), i = 1, ..., #FB.
def ideal_sq_eq_test(K, FB, I, e, g):
    assert len(FB) == len(e), f"Wrong exponent vector e, len(e) = {len(e)} != {len(FB)}"
    nrm = prod([K.ideal(FB[i].prime, K(FB[i].elts)).absolute_norm()^ZZ(e[i]) for i in range(len(e))])
    known_factors = [I.q] + [FB[i].prime for i in range(len(FB))]
    #print(f"[ideal_sq_eq_test] I = {I}\n\tN(I^2) = {(I.q^2).factor()}\n\tN(g J) = ({g.absnorm().factor()}) * ({(nrm).factor()})")
    if abs(I.q)^2 != abs(g.absnorm(known_factors = known_factors) * nrm): # compare upto units
        #print("\t[fail (norms)]")
        return False
    
    I_sage = I.to_sage(K)
    for i in range(len(e)):
        P = K.ideal(FB[i].prime, K(FB[i].elts))
        if g.valuation(P) + e[i] != I_sage.valuation(P) * 2:
            #print("\t[fail (valuations)]")
            return False
    return True

def save_result(d, G, I):
    if SAVE_ALL_RESULTS != None:
        fld = "_".join([str(d[i]) for i in range(len(d))])
        fn = f"{SAVE_ALL_RESULTS}/{fld}_dlog"
        with open(fn, 'w') as f: f.write(str(G))
        with open(f"{fn}_input", 'w') as f: f.write(str(I))\

def clean_results():
    trees.clear_dir("experiments/dlogs")

def cl_dlog(d, I, d_parent = ()):
    '''
    For an ideal I returns an ideal product g * prod_i(P[i]^e[i]) = I,
    where g is a field element in compact representation and P[i] are primes from the
    factor base (precomputed during class group computation)
    '''

    print(f"Init data for d = {d} ... ", end="", flush=True)
    t = walltime()
    print("[fields]", end=" ", flush=True)
    K = field.field(d).sage(food)
    KC = field_compact.field_compact(d)
    print("[clgp]", end=" ", flush=True)
    CL = clgp.clgp(d, d_parent=d_parent)

    print("[matrices]", end=" ", flush=True)
    M,B,U,V,U_inv,V_inv = clgp.get_matrices(d)

    #FB = clgp.get_FB(d, d_parent)
    FB = CL.factor_base()
    print("[fb]", end=" ", flush=True)

    IPS = idealprod.idealprods(d, food=food, FB=FB)
    print("[idealprod]", end=" ", flush=True)
    print(f"[ok] {walltime(t)} sec.", flush=True)

    d = tuple(d)
    if len(d) == 1:
        # should return smooth ideal J and an element g, s.t. I*J = g O_K
        g, e = dlog_quad.dlog_quad(d[0], I, d_parent=d_parent, food=food)
        G = IPS(KC(g), -vector(e))
        save_result(d, G, I)
        return [G]

    B_t = [ZZ(B[i,i]) for i in range(B.ncols())]
    print(f"Computing DLOG for ideal I = {I} in the group {clgp.strip_oz(B_t)} ... ", flush=True)
    comp_time = walltime()

    if set(B_t) == {1}:
        e = [0]*M.ncols()
        e_g = clgp.primes_to_gens(e, V, B, strip_zeroes=True)
        #print(f"h_K = 1, computation of dlog is trivial")
        g = I.generator()
        print(f"Finished computing DLOG for ideal I = {I} in the group {clgp.strip_oz(B_t)} ... {e_g}, {walltime(comp_time)} sec.")
        print(f"-> I = < {g.to_sage(K)} >\n")
        #assert K.ideal(g.to_sage(K)) == I.to_sage(K), f"Wrong computation for trivial class group (d = {d})"
        if verify.level() > verify.LIGHT or (verify.level() == verify.LIGHT and len(d) <= 5):
            print(f"-> Sage's result: {I.to_sage(K).ideal_class_log(proof=False)}")
            print()
        G = IPS(KC(g), e)
        save_result(d, G, I)
        return [G]

    d_s  = d[:-1]
    d_t = d[:-2] + d[-1:]
    d_st = d[:-2] + (d[-2]*d[-1],)

    I_s = ideal.ideals(d_s)(I.q, I.s[:-1])
    I_t = ideal.ideals(d_t)(I.q, I.s[:-2] + I.s[-1:])
    I_st = ideal.ideals(d_st)(I.q, I.s[:-2] + (I.s[-2]*I.s[-1],))

    # loading precomputed ramification trees with factor base
    CL_s = clgp.clgp(d_s, d_parent=d)
    CL_t = clgp.clgp(d_t, d_parent=d)
    CL_st = clgp.clgp(d_st, d_parent=d)

    FB_s = CL_s.factor_base()
    FB_t = CL_t.factor_base()
    FB_st = CL_st.factor_base()

    #K_s = field.field(d_s).sage(food)
    #K_t = field.field(d_t).sage(food)
    #K_st = field.field(d_st).sage(food)

    # computing DLOGs in subfields
    print(f"Lift dlog into subfield {d_s}")
    dlog_s  = cl_dlog(d_s, I_s, d_parent=d)

    print(f"Lift dlog into subfield {d_t}")
    dlog_t  = cl_dlog(d_t, I_t, d_parent=d)

    print(f"Lift dlog into subfield {d_st}")
    dlog_st = cl_dlog(d_st, I_st, d_parent=d)

    print(f"{len(dlog_s)*len(dlog_t)*len(dlog_st)} triples lifted from subfields for the field {d}")

    if len(dlog_s) == 0 or len(dlog_t) == 0 or len(dlog_st) == 0:
        print(f"Failed to find dlog for the ideal I = {I}")
        return []

    #K_st = NumberField([x^2 - d_st[i] for i in range(len(d_st))], names=list(sorted(food.keys()))[:len(d_st)])
    #FB_st_sage = [K_st.ideal(FB_st[i].prime, K_st(FB_st[i].elts)) for i in range(len(FB_st))]
    #FB_st_conj = [K_st.ideal(FB_st[i].prime, apply_aut(K_st, K_st(FB_st[i].elts), [0]*(len(d_st)-1) + [1])) for i in range(len(FB_st))]

    I_sage = I.to_sage(K)
    #FB_sage = [FB[i].to_sage(food) for i in range(len(FB))]

    res = []

    for H_s in dlog_s:
        # Check correstness of DLOG computation in subfield
        #assert ideal_eq_test(K_s, FB_s, I_s, H_s.powers, H_s.element), f"incorrect dlog computation for the field {d_s}"

        for H_t in dlog_t:
            # Check correstness of DLOG computation in subfield
            #assert ideal_eq_test(K_t, FB_t, I_t, H_t.powers, H_t.element), f"incorrect dlog computation for the field {d_t}"

            for H_st in dlog_st:
                # Check correstness of DLOG computation in subfield
                #assert ideal_eq_test(K_st, FB_st, I_st, H_st.powers, H_st.element), f"incorrect dlog computation for the field {d_st}"

                print(f"Applying norm equation... ", end="", flush=True)
                t = walltime()

                e_s_lift = clgp.lift_e(H_s.powers, FB_s)
                e_t_lift = clgp.lift_e(H_t.powers, FB_t)

                # We assume that the automorphism sigma is applied to the prime ideals in trees for the subfield field K_st.
                e_st_lift = clgp.lift_e(H_st.powers, FB_st)

                assert(len(e_s_lift) == len(e_t_lift))
                assert(len(e_t_lift) == len(e_st_lift))

                # applying norm equation
                e = [e_s_lift[i] + e_t_lift[i] - e_st_lift[i] for i in range(len(e_s_lift))]
                #e = [e_st_lift[i] - e_s_lift[i] - e_t_lift[i] for i in range(len(e_s_lift))]

                g_s = KC(H_s.element)
                g_t = KC(H_t.element)
                g_st = KC(H_st.element)

                g = g_s * g_t / g_st.conj()
                
                #print(f"g = {g.to_sage(K)}")

                G = IPS(g, vector(e))
                #print(f"e = {e}")
                print(f"{walltime(t)} sec.", flush=True)

                #print(f"I^2 = {G.to_sage(K)}")
                #print(f"-> G.element: {G.element.to_sage(K)}")
                #print(f"-> G.powers: {G.powers}")
                #print(f"-> G.element (norm factor): {[G.element.elements[i].absnorm().factor() ^ G.element.powers[i] for i in range(len(G.element.elements))]}")

                print(f"\n-> loading S-units for {d} ... ", end="", flush=True)
                t = walltime()
                SU = sunits.sunits(d, d_parent=d_parent)
                SU.preload()
                print(f"{walltime(t)} sec.", flush=True)

                clgp.check_clgp(d, d_parent, food=food)
                
                print(f"Reducing exponents ... ", end="", flush=True)
                t = walltime()
                G = G.reduce_powers_mod_sunits(SU)
                print(f"{walltime(t)} sec.", flush=True)

                #if len(d) <= 3:
                #     g = G.element
                #     print(f"g (after reducing exponents) = {g.to_sage(K)}")

                #     g = g.compactify()
                #     print(f"g (compactified) = {g.to_sage(K)}")
                #     g = g.compute_roots(ideal_sqrt=True, compact_repr=True).compactify(base=True)
                #     print(f"g (expanded roots) = {g.to_sage(K)}")

                #     G = IPS(g, tuple(G.powers))

                e_g = clgp.primes_to_gens(e, V, B, strip_zeroes=True, reduce=True)
                print(f"Computing square root for {e_g} in class group {clgp.strip_oz(B_t)} of the field {d}... ", flush=True)

                G = G.sqrt_rat(SU, quadchars=args.quadchars)
                
                print(f"-> G.element: {G.element.to_sage(K)}")
                print(f"-> G.powers: {G.powers}")
                
                # print(f"Reducing exponents (2) ... ", end="", flush=True)
                # t = walltime()
                # G0 = G
                # G = G.reduce_powers_mod_sunits(SU)
                
                # assert G0.evaluate() == G.evaluate()
                # print(f"{walltime(t)} sec.", flush=True)
            
                # print(f"-> G.element: {G.element.to_sage(K)}")
                # print(f"-> G.powers: {G.powers}")
                
                #G.element = G.element.compute_roots(ideal_sqrt=True)
                #G.element = G.element.compactify()
                
                profile.create_stats()
                pstats.Stats(profile).strip_dirs().sort_stats(SortKey.CUMULATIVE).print_stats(int(50))
                profile.enable()

                e_g = clgp.primes_to_gens(vector(G.powers), V, B, strip_zeroes=True, reduce=True)
                print(f"Finished computing DLOG for ideal I = {I} in the group {clgp.strip_oz(B_t)} ... {e_g}, {walltime(comp_time)} sec.\n", flush=True)
                
                save_result(d, G, I)
                
                if verify.level() > verify.LIGHT or (verify.level() == verify.LIGHT and len(d) <= 5):
                    tt = walltime()
                    print(f"-> Sage's result: {I_sage.ideal_class_log(proof=False)}, {walltime(tt)} sec.")
                    #print(f"-> is_principal(I)?: {I_sage.is_principal(proof=False)}")
                    #print(f"-> is_principal(I^2)?: {(I_sage^2).is_principal(proof=False)}")
                    print()
                
                return [G]
    return res

trees_food = trees.get_food()
if food != None:
    assert food == trees_food, f"Run trees generation and relation computation first! Trees are generated for {trees_food} != {food}."

print(f"food = {food}", flush=True)
#d = tuple(sorted(food.values(), key=lambda di: abs(di)))
d = trees.get_d()
#d = d[:-1]
print(f"d = {d}", flush=True)
d_parent = ()
bound = prod(di.abs() for di in d)
#I = ideal.random(d, bound=bound)
I = ideal.random(d)

# The ideal can be fixed as in the following examples.
#I = ideal.ideals((5, 13, 17, 29, 37, 41))(9032339,(2001020, 3845142, 269749, 1568238, 3382698, 560462))

print(f"Target ideal: {I}", flush=True)
#K = NumberField([x^2 - d[i] for i in range(len(d))], names=list(sorted(food.keys()))[:len(d)])
#I_hnf = I.hnf(K)
#print(f"Target ideal (hnf): {I_hnf}")
#I_sage = I.to_sage(K)
#print(f"I_sage.norm = {I_sage.absolute_norm()}")
#assert I_sage.absolute_norm() == I.q

FB = clgp.get_FB(d, d_parent=d_parent)

assert not fb.is_smooth(I.q, FB), "This code is for non-smooth ideals. For smooth ideal run dlog_smooth.sage"

for i in range(len(FB)):
    assert gcd(FB[i].prime, I.q) == 1 # FIXME: change ideal_test_eq

clean_results()

profile.enable()
t = walltime()
dls = cl_dlog(d, I, d_parent = d_parent)
print(f"Computation time: {walltime(t)} sec.", flush=True)
profile.disable()
