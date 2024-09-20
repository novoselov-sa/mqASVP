# This script uses precomputed dlogs in experiments/dlogs.
# So, you should run the dlog_cyc_rat.sage script first.

import idealprod
import trees
import sunits
import ideal
import ring
import field
import field_compact
from pathlib import Path
import relations
import fb
import clgp
import random
import sys
import numpy
import argparse

from fpylll import FPLLL

parser = argparse.ArgumentParser(description='Finding a short element of an ideal of n-multiquadratic field given precomputed input ideal, a solution of discrete logarithm problem (stored in dlogs/ folder), and class group (folders "relations" and "trees").')
parser.add_argument("n", type=Integer, help="n for n-quadratic field")
parser.add_argument("-lp", type=Integer, dest="log_emb_precision", default=10000, help="precision for computing of log-embeddings")
parser.add_argument("-p", type=Integer, dest="emb_precision", default=50000, help="precision for computing of the cannonical embedding")
parser.add_argument("-f", type=Integer, dest="fplll_precision", default=1000, help="FPLLL precision")
parser.add_argument("-fi", type=Integer, dest="fplll_precision_inc", default=None, help="increment FPLLL's precision by this value when the result of CVP computation doesn't fit the expected bound (default: None, disabled)")
parser.add_argument("-c", "--check", dest="check", action=argparse.BooleanOptionalAction, default=False, help="evaluate output element explicitly and check that it belongs to the input ideal (extemely slow)")
parser.add_argument("-ds", "--disable-stages", dest="disabled_stages", type=Integer, nargs='+', default=[], help='disable specified stages')
parser.add_argument("--smooth", dest="smooth", action=argparse.BooleanOptionalAction, default=False, help="Assume that input ideal is smooth. Target ideal is selected as an ideal of the biggest norm in the factor base.")

args = parser.parse_args()

precision = args.log_emb_precision
R = RealField(precision)
print(f"log-embeddings precision: {precision}")
print(f"canonical embedding precision: {args.emb_precision}")

# This is for Babai's algorithm in Approx-CVP.
# Try to increase this value if Approx-CVP solver doesn't reduce target vector enough.
FPLLL.set_precision(args.fplll_precision)
print(f"FPLLL precision: {FPLLL.get_precision()}")

proof.all(False)

food = trees.get_food()

d_top = tuple(sorted(food.values(), key=lambda di: abs(di)))

n = args.n
d_parent = d_top[:n+1]
d = d_top[:n]

m = 2^len(d)
D = prod(d).abs()

K_sage = field.field(d).sage()

file = relations.convert_letters(d, seed=1, food=food)
FB = fb.load_primes(tuple(d), file, food=food)
IPS = idealprod.idealprods(d, food=food, FB=FB)
K = field.field(d)
KC = field_compact.field_compact(d)

print(f"d_top = {d_top}")
print(f"d_parent = {d_parent}")
print(f"d = {d}")

if args.smooth:
    i_max = 0
    for i in range(1,len(FB)):
        if FB[i].prime > FB[i_max].prime:
            i_max = i
    I = FB[i_max]
    v = [0] * len(FB)
    v[i_max] = 1
    G = IPS(KC.one(), tuple(v))
else:
    fn = "_".join([str(d[i]) for i in range(len(d))])
    I = eval(preparse(Path(f"experiments/dlogs/{fn}_dlog_input").read_text()))

    G = eval(preparse(Path(f"experiments/dlogs/{fn}_dlog").read_text()))

print(f"Target ideal:\n{I}")
print(f"\nDLOG:\n{G}")

print(f"\nlognorm_2: {G.element.lognorm(2)}")
#print(f"lognorm_2 (sage): {G.element.lognorm_sage(2, prec=precision)}") # slow for deg K >= 64

assert len(G.powers) == len(FB)

print("\nLoading S-units ...")
SU = sunits.sunits(d, d_parent=d_parent)
SU.preload()

assert len(G.powers) == SU.relations_matrix().ncols()

tt = walltime()

found = False

if 1 not in args.disabled_stages:
    print("\n\nstage 1 ...\n")
    G_red = G.short_element(SU, fprec_increment=args.fplll_precision_inc)
    print(f"powers without drift:\n{G_red.powers}")

    if min(G_red.powers) >= 0:
        s0 = 0
        drift0 = vector(ZZ, [s0]*len(G_red.powers))
        G_red0 = G_red
        found = True
        print(f"found good drift:\n{drift0}")
        print(f"new powers:\n{G_red0.powers}")
        print(f"new norm: {RR(vector(ZZ, G_red0.powers).norm(2))}")

    if not found:
        drift = []
        for i in range(len(G_red.powers)):
            s0 = abs(min(G_red.powers))
            if G_red.powers[i] < 0:
                drift.append(-G_red.powers[i])
            else:
                drift.append(0)
        drift = vector(ZZ, drift)
        #print(f"stage 1 drift 0: {drift}")
        G_red = G.short_element(SU, drift=drift, fprec_increment=args.fplll_precision_inc)
        #print(f"stage 1 powers 0: {G_red.powers}")
        if min(G_red.powers) >= 0:
            s0 = floor(numpy.mean(drift)) #max(drift)
            drift0 = drift
            G_red0 = G_red
            found = True
            print(f"found good drift:\n{drift0}")
            print(f"new powers:\n{G_red0.powers}")
            print(f"new_norm: {RR(vector(ZZ, G_red0.powers).norm(2))}")
        else:
            for j in range(m):
                for i in range(len(G_red.powers)):
                    if G_red.powers[i] < 0:
                        drift[i] = -G_red.powers[i]
                #print(f"stage 1 drift {j+1}: {drift}")
                G_red = G.short_element(SU, drift=drift, fprec_increment=args.fplll_precision_inc)
                #print(f"stage 1 powers {j+1}: {G_red.powers}")
                if min(G_red.powers) >= 0:
                    s0 = floor(numpy.mean(drift)) #max(drift)
                    drift0 = drift
                    G_red0 = G_red
                    found = True
                    print(f"found good drift:\n{drift0}")
                    print(f"new powers:\n{G_red0.powers}")
                    print(f"new_norm: {RR(vector(ZZ, G_red0.powers).norm(2))}")
                    break

if 2 not in args.disabled_stages:
    print("\nstage 2 ...")
    bound = ceil(RR(sqrt(m * log(D)/log(2))))+1
    for i in range(m):
        if i < 2:
            s = i
        elif i == 2:
            s = bound
        else:
            s = ZZ.random_element(0, bound)

        print(f"drift (base) = {s}", flush=True)
        drift = vector(ZZ, [s]*len(G.powers))
        G_red = G.short_element(SU, drift=drift, fprec_increment=args.fplll_precision_inc)

        if min(G_red.powers) >= 0 and (not found or vector(ZZ, G_red.powers).norm(2) < vector(ZZ, G_red0.powers).norm(2)):
            s0 = s
            drift0 = drift
            G_red0 = G_red
            found = True
            print(f"found good drift:\n{drift0}")
            print(f"new powers:\n{G_red0.powers}")
            print(f"new_norm: {RR(vector(ZZ, G_red0.powers).norm(2))}")

        for i in range(5):
            if s == 0:
                v = vector(ZZ, random.choices([0,1], k=len(G.powers)))
            else:
                v = vector(ZZ, random.choices([0,1,-1], k=len(G.powers)))
            drift = vector(ZZ, [s]*len(G.powers)) + v
            G_red = G.short_element(SU, drift=drift, fprec_increment=args.fplll_precision_inc)
            if min(G_red.powers) >= 0:
                if not found or vector(ZZ, G_red.powers).norm(2) < vector(ZZ, G_red0.powers).norm(2):
                    s0 = s
                    drift0 = drift
                    G_red0 = G_red
                    found = True
                    print(f"found good drift:\n{drift0}")
                    print(f"new powers:\n{G_red0.powers}")
                    print(f"new_norm: {RR(vector(ZZ, G_red0.powers).norm(2))}")

if 3 not in args.disabled_stages:
    assert found, "failed to find s0"

    print("\nstage 3 (finding the best neighbor) ...")
    print(f"s0 = {s0}")
    print(f"drift0 = {drift0}")

    for i in range(m):
        v = vector(random.choices(range(ceil(0.9*s0), floor(1.1*s0)+1), k=len(G.powers)))
        if s0 == 0:
            b0 = vector(ZZ, random.choices([0,1], k=len(G.powers)))
        else:
            b0 = vector(ZZ, random.choices([0,1,-1], k=len(G.powers)))
        drift = v + b0
        G_red = G.short_element(SU, drift=drift, fprec_increment=args.fplll_precision_inc)
        if min(G_red.powers) >= 0:
            if vector(ZZ, G_red.powers).norm(2) < vector(ZZ, G_red0.powers).norm(2):
                s0 = s
                drift0 = drift
                G_red0 = G_red
                found = True
                print(f"found good drift:\n{drift0}")
                print(f"new powers:\n{G_red0.powers}")
                print(f"new_norm: {RR(vector(ZZ, G_red0.powers).norm(2))}")

assert found
G_red = G_red0

print(f"\nfinal drift:\n{drift0}")

print(f"\nDLOG (reduced):\n{G_red}\n")

assert min(G_red.powers) >= 0, "wrong drift"

lognorm_inf = G_red.element.lognorm_sage(infinity, prec=precision)

print(f"lognorm_2: {(G_red.element.lognorm(2, prec=precision))}")
print(f"lognorm_inf: {(G_red.element.lognorm(infinity, prec=precision))}")
#print(f"lognorm_2 (sage): {(G_red.element.lognorm_sage(2, prec=precision))}")
#print(f"lognorm_inf (sage): {lognorm_inf}")

#assert G_red.absnorm() == G.absnorm() # this condition doesn't hold since we added a drift.

I_vol = sqrt(abs(K.discriminant())) * I.absolute_norm()
print(f"Vol(I) = {I_vol}")

gh = sqrt(m/(2*pi*exp(1))) * I_vol^(1/m)

print(f"Gaussian Heuristic = {RR(gh)}")

print("Shortening ...")
el_s = G_red.element.shorten(prec=precision)

print(f"\ncomputation time: {walltime()-tt} sec.")

print(f"\nel (shortened) = {el_s}\n")
print(f"-> lognorm_2: {(el_s.lognorm(2, prec=precision))}")
print(f"-> lognorm_inf: {lognorm_inf}")

print(f"-> lognorm_2 (sage): {(el_s.lognorm_sage(2, prec=precision))}")
lognorm_inf = el_s.lognorm_sage(infinity, prec=precision)
print(f"-> lognorm_inf (sage): {lognorm_inf}")
print(f"-> log(af_gh) (estim.): {R(lognorm_inf + log(sqrt(m)) - log(gh))}")

R2 = RealField(args.emb_precision)
el_s_emb = vector(R2, el_s.canonical_embedding_sage(prec=args.emb_precision))
#print(f"el_s_emb = {el_s_emb}")
el_inf = el_s_emb.norm(infinity)
el_2 = el_s_emb.norm(2)
print(f"\n|el|_inf = {el_inf}")
print(f"\n|el|_2 = {el_2}")
print(f"\nlog(af_gh): {R(log(el_2/gh))}")

if args.check:
    el = el_s
    # trying to reduce the element and simplify computations of the roots.
    el = el.compactify(base=True).trim()
    #print(f"el (compactified) = {el}")
    el = el.compute_roots(ideal_sqrt=True).evaluate(mod_units=True)
    el = el.reduce_mod_units()
    print(f"\nel (explicit) = {el}\n")
    print(f"lognorm_2: {(el.lognorm(2))}")
    print(f"lognorm_2 (sage): {(el.lognorm_sage(2, prec=precision))}")
    emb = el.canonical_embedding_sage(prec=args.emb_precision)
    el_inf = vector(R2, emb).norm(infinity)
    el_2 = vector(R2, emb).norm(2)
    print(f"\n|el|_inf = {el_inf}")
    print(f"\n|el|_2 = {el_2}")

    print(f"\naf_gh: {R(el_2 / gh)}")
    print(f"\nln(af_gh): {R(log(el_2 / gh))}")

    lambda_1_2_bound_up = R(sqrt(K.degree()) * I_vol^(1/K.degree()))
    lambda_1_2_bound_lo = R(I.absolute_norm()^(1/K.degree()))
    approx_factor_up_bound = R(el_2 / lambda_1_2_bound_lo)
    approx_factor_lo_bound = R(el_2 / lambda_1_2_bound_up)
    print(f"\naf_sup = {approx_factor_up_bound}")
    print(f"\naf_inf = {approx_factor_lo_bound}")
    assert el.to_sage(K_sage) in I.to_sage(K_sage)
