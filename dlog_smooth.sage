import ideal
import relations
import fb
import trees
import clgp

food = trees.get_food()

d = tuple(sorted(food.values(), key=lambda di: abs(di)))

proof.all(False)

#pari.allocatemem(1024*1024*1024)
#print(f"Pari stack size: {pari.stacksize()}")

trees_food = trees.get_food()
if food != None:
    assert food == trees_food, f"Run trees generation and relation computation first! Trees are generated for {trees_food} != {food}."

bound = prod(di.abs() for di in d)
I = ideal.random(d, bound=bound)
#I = ideal.ideals([5, 229, 257])(669289,(303651, 6337, 202003))

print(f"Target ideal: {I}")
K = NumberField([x^2 - d[i] for i in range(len(d))], names=list(sorted(food.keys()))[:len(d)])
I_hnf = I.hnf(K)
print(f"Target ideal (hnf): {I_hnf}")
I_sage = I.to_sage(K)
print(f"I_sage.norm = {I_sage.absolute_norm()}")
assert I_sage.absolute_norm() == I.q

file = relations.convert_letters(d, seed=1, food=food)
FB = fb.load_primes(tuple(d), file, food=food)
M = relations.load_matrix(d)
B,U,V=M.smith_form()
V_inv = V^(-1)

assert fb.is_smooth(I.q, FB), f"Ideal is not smooth. Run class group generation including all factors of in the norm of ideal!"

print(f"CL_K: {[B[i,i] for i in range(B.nrows())]}")

w = walltime()
#print(f"FB = {FB}")
e = fb.factor(K, I_sage, FB)
print(f"e = {e}")

e_g = clgp.primes_to_gens(e, V, B)
print(f"cl_dlog: {e_g}, time: {walltime(w)}", flush=True)

print("Checking using Sage ...") # Unreachable for n >= 6
if len(d) >= 6:
    print("-> DISABLED")
    exit(0)

w = walltime()
cl_sage = I_sage.ideal_class_log(proof=False)
cl_sage = [ZZ(i) for i in cl_sage]
print(f"-> [sage] cl_dlog = {cl_sage}, time: {walltime(w)}")
print(f"-> [sage] I.factor() = {I_sage.factor()}")

e_rev = clgp.gens_to_primes(e_g, V_inv)
print(f"e_rev = {e_rev}")
I2 = fb.prod_ideal(K, FB, e_rev)
CL_K = K.class_group(proof=False)
assert CL_K(I_sage) == CL_K(I2)

#print([CL_K.gens()[i].order() for i in range(len(CL_K.gens()))])
I3 = prod([CL_K.gens()[i]^cl_sage[i] for i in range(len(CL_K.gens()))])
assert CL_K(I_sage) == CL_K(I3)
