# This is a script for verification of class group computation.

import field
import units
import argparse
import trees
import clgp

prec = 1000
RF = RealField(prec)

parser = argparse.ArgumentParser(description='Verification of class group computation for multiquadratic field.')
parser.add_argument("mquad_field", type=Integer, nargs='*', help='list of integers d_1, ..., d_n describing multiquadratic field. Loaded from precomputed data if ommited.')
parser.add_argument("--prec", type=Integer, dest="prec", default=1000, help="precision")
parser.add_argument("-r", type=RF, dest="regulator", help="regulator, if not given the approximation is computed")
parser.add_argument("--class-number", type=Integer, dest="class_number", help="class number")
parser.add_argument("--class-group", type=Integer, nargs='+', dest="class_group", help="class group structure")
parser.add_argument("--log-zeta-res", type=RF, required=True, help="residue of the zeta function of O at 1. Use zeta_log_residue from Hecke to compute.")

args = parser.parse_args()

if args.mquad_field != []:
    d = tuple(args.mquad_field)
else:
    d = trees.get_d()

print(f"Field: {d}")

print(f"Precision: {args.prec}")

r = args.regulator
if r == None:
    print("Computation of regulator approximation ...")
    r = units.approxregulator(d, prec=args.prec)

print("Regulator: ", r)

n = len(d)
K = field.field(d)

print(f"Discriminant: {K.discriminant()}")

z_log = args.log_zeta_res
print(f"Euler product (log): {z_log}")
z = exp(z_log)
print(f"Euler product: {z}")

print(f"Field init")
L = NumberField([x^2-di for di in d]).absolute_field(names="a")

print(f"Computing number of roots of unity")
w = units.torsion(d)[0]
#assert w == L.number_of_roots_of_unity()
print(f"Num. of roots of unity: {w}")

r1,r2 = K.signature()
print(f"r1 = {r1}, r2 = {r2}")

if args.class_number == None and args.class_group == None:
    cl = clgp.structure(d) # load from precomputed data
    print(f"Class group = {cl}")
    h = prod(cl)
elif args.class_number == None:
    cl = args.class_group
    print(f"Class group = {cl}")
    h = prod(cl)
else:
    h = args.class_number

print(f"Class number = {h}")

pr1 = RF(w * sqrt(abs(K.discriminant())) / (2^r1 * (2 * pi)^r2) * z)
print(f"h(K)*R(K) from class number formula: {pr1}")
pr2 = h*r
print(f"Given product h(K)'*R(K)': {pr2}")
df = abs(log(pr2) - log(pr1))
print(f"log(h'(K)*R'(K) / h(K)*R(K)): {df}")

# since we know regulator, we do not need epsilon in zeta_log_residue.
if df < 1/2*log(2.0):
    print("OK")
else:
   #print("FAIL")
   if pr1/pr2 < 1:
       #print(f"FAIL: we are in an overgroup of class group. Approximate index: {floor(pr2/pr1)}")
       print(f"FAIL: h'/h_k = {round(pr2/pr1)}")
   else:
       #print(f"FAIL: we are in a subgroup of class group. Approximate index: {floor(pr1/pr2)}")
       print(f"FAIL: h'(K)/h(K): 1/{round(pr1/pr2)}")
