import clgp
import trees
import relations

import argparse
parser = argparse.ArgumentParser(description='Norms of rows in relation matrices.')
parser.add_argument("mquad_field", type=Integer, nargs='*', help='list of integers d_1, ..., d_n describing multiquadratic field. Loaded from precomputed data if ommited.')
parser.add_argument("--gamma", action=argparse.BooleanOptionalAction, default=False, help="print data for subfield of type (d_1, ..., d_{n-2}, d_n)")
parser.add_argument("--mu", action=argparse.BooleanOptionalAction, default=False, help="print data for subfield of type (d_1, ..., d_{n-1}), default")
parser.add_argument("--gamma-mu", action=argparse.BooleanOptionalAction, dest="gamma_mu", default=False, help="print data for subfield of type (d_1, ..., d_{n-1}*d_{n})")
parser.add_argument("-m", type=Integer, dest="degree", help='degree of the subfield')

args = parser.parse_args()

d_parent = trees.get_d()

print(f"d_parent = {d_parent}")

if args.mquad_field != []:
    d = d_parent
else:
    d = trees.get_d()

subfield = False
if args.degree != None:
    if args.degree > 2^(len(d)) / 2 or args.degree < 2 or not args.degree.is_power_of(2):
        print(f"The degree of subfield should be power of two in the range 2 .. {2^(len(d)) / 2}")
    ms = floor(log(args.degree)/log(2))
    d = d[:ms+1]
    subfield = True

subfield = subfield or args.gamma or args.gamma_mu or args.mu

if subfield:  
    if args.gamma:
        d = d[:-2] + (d[-1],)    
    elif args.gamma_mu:
        d = d[:-2] + (d[-1] * d[-2],)
    else:
        d = d[:-1]

print(f"d = {d}")

N = 2^len(d)
print(f"N = {N}")
bound = sqrt(N)
print(f"sqrt(N) = {RR(sqrt(N))}")
print(f"2*sqrt(N) = {RR(2*sqrt(N))}")
print(f"sqrt(N*log_2(N)) = {RR(sqrt(N*log(N)/log(2)))}")

print(pari.allocatemem(60*1024^3))

#M,B,U,V,U_inv,V_inv = clgp.get_matrices(d)
M = relations.load_matrix(d)

print("CL:", clgp.structure(d))

bad = []
mx_2 = 0
mx_inf = 0
for r in M.rows():
    if RR(r.norm(2)) > bound:
        print(f"r = {r}")
        bad.append([RR(r.norm(2)), r.norm(infinity)])
    if RR(r.norm(2)) > mx_2:
        mx_2 = RR(r.norm(2))
    if RR(r.norm(infinity)) > mx_inf:
        mx_inf = r.norm(infinity)
    print([RR(r.norm(2)), r.norm(infinity)])

print(f"max |r|_2 = {mx_2}")
print(f"max |r|_inf = {mx_inf}")
print(f"max |r|_2 <= sqrt(N): {mx_2 <= RR(bound)}")
print(f"max |r|_2 <= 2*sqrt(N): {mx_2 <= RR(2*bound)}")
print(f"max |r|_2 <= sqrt(N log_2(N)): {mx_2 <= RR(sqrt(N*log(N)/log(2)))}")

