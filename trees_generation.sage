import nprofile
import trees
import prime_decomp
import polynomial_ring
import field
import sys
import argparse

food = None

parser = argparse.ArgumentParser(description='Generation of trees describing splitting of primes over subfields of multiquadratic field.')
parser.add_argument("mquad_field", type=Integer, nargs='+', help='list of d_1, ..., d_n describing multiquadratic field')
parser.add_argument("--bound", type=Integer, dest="bound", help="bound for primes")
parser.add_argument("--gm-add", type=Integer, dest="gm_add", default=120, help="value to add to maximum of GM-bounds for quadratic subfields")
parser.add_argument("--primes", type=Integer, dest="primes", nargs="+", help='additional primes for building the factor base', default=[])
parser.add_argument("--hard-primes", action=argparse.BooleanOptionalAction, default=True, help="use primes p s.t. p | [O_K : Z[theta_K]] or not")
parser.add_argument("--ramified-primes", action=argparse.BooleanOptionalAction, default=True, help="use primes p s.t. p | disc(K) or not")
parser.add_argument("--split-primes-only", action=argparse.BooleanOptionalAction, default=False, help="use only primes that split completely")

args = parser.parse_args()

if food == None:
    l = "abcdeghijklmnopqrstuvwxyz"
    d = tuple(sorted(args.mquad_field, key=lambda di: abs(di)))
    food = {}
    for i in range(len(d)):
        food[l[i]] = d[i]
else:
    print(f"Warning: input field is hardcoded ({food}), ignoring command line arguments.")

# Additional primes to factor base. For DLOG computation add here primes from the factorization of norm of target ideal.
PRIMES = args.primes

trees.USE_RAMIFIED_PRIMES = args.ramified_primes
trees.USE_EXCEPTIONAL_PRIMES = args.hard_primes
trees.USE_SPLIT_PRIMES_ONLY = args.split_primes_only

d = tuple(sorted(food.values(), key=lambda di: abs(di)))
gens = list(sorted(food.keys()))
working_folder="trees/"

# set this variable to True if you want to compute norm bound using Grenié-Molteni algorithm.
COMPUTE_GM_NBOUND = (len(d) <= 5)
COMPUTE_GM_NBOUND_QUAD = True

print("Generating trees for field: {})".format(d))

norm_bound = trees.compute_bound(d)
print("Bach bound =", norm_bound)

# compute norm bound using Grenié-Molteni algorithm
if COMPUTE_GM_NBOUND: # computation may be very slow 
    gp.read("norm_bounds/GM_bounds.gp")
    K = NumberField([x^2 - di for di in food.values()])
    GM_bound = ZZ(gp.GRHoptimize(K.pari_nf()))
    print("GM bound =", GM_bound)

# compute maximal bound for quadratic subfields.
if COMPUTE_GM_NBOUND_QUAD:
    gp.read("norm_bounds/GM_bounds.gp")
    I = [1]
    for di in d:
        I += [di*i for i in I]
    GM_bound_quad = 0
    for di in I:
        if di == 1:
            continue
        K = QuadraticField(di)
        b = ZZ(gp.GRHoptimize(K.pari_nf()))
        print(f"GM bound for D = {di}: {b}")
        GM_bound_quad = max(b, GM_bound_quad)
    print("GM bound (quadr.subf.) =", GM_bound_quad)

# comment out the following line if you want to use Bach bound
if args.bound == None:
    norm_bound = GM_bound_quad + args.gm_add
else:
    norm_bound = args.bound
# norm_bound = GM_bound # bound computed using algorithm of Grenié-Molteni

print("Used bound =", norm_bound)

print("Computing trees ...")
t = walltime()
K = field.field(d)
K.set_names(food)
tr = trees.build_tree(K, {}, norm_bound=norm_bound, append_primes=PRIMES)
print(f"-> done in {walltime(t)} sec.")

print("Removing existing trees ...")
trees.clean_trees(folder=working_folder)
print("-> done")

print("Saving new trees ...")
trees.save_tree(tr, K, folder=working_folder)
trees.save_food(food, folder=working_folder)
trees.save_field(d, folder=working_folder)
print("-> done")

nprofile.output([trees, prime_decomp, polynomial_ring])
