# Computation of square roots for S-smooth ideal I using saturation.
# Requires precomputed class group.

import field
import field_compact
import ideal
import clgp
import fb
import units
import sunits
import saturation

# Finding a set of elements s such that s*h is a square, where h is a field 
# element and s is from set S.
# The r determines a number of characters to use (r + 64).
def find_squares_mult(d, h, S, r):
    low,high = 0,r + 64
    h_s = h.symbols2_log(low,high)
    A = matrix(GF(2), len(S) + 1, len(h_s))
    A[0] = h_s
    for i in range(len(S)):
        A[i+1] = S[i].symbols2_log(low,high)
    #print(f"A = {A}")
    V = A.left_kernel().basis_matrix().rows()
    #print(f"V = {V}")
    K = field.field(d)
    KC = field_compact.field_compact(d)
    squares = []
    for v in V:
        if v[0] == 0:
            continue
        assert v[0] == 1
        sq = KC(K.one())
        for j in range(1,len(v)):
            sq *= S[j-1] ^ ZZ(v[j])
        squares.append(sq)
    return squares

# Generators of prime ideals (or their inverses) in FB according to coefficients in alpha.
# Parameters:
#   FB is a set of prime ideals of the field K = field.field(d).sage() (Sage).
#   alpha is an integer vector.
# Returns a set with coefficients in the field.field(d).
def prime_ideals_generators(d, alpha, FB):
    res = [[]]*len(alpha)
    for i in range(len(alpha)):
        if alpha[i] >= 0:
            res[i] = FB[i].gens()
        elif alpha[i] < 0:
            res[i] = clgp.prime_ideal_inverse(FB[i]).gens()
        res[i] = [field.field(d).from_sage(res[i][j]) for j in range(len(res[i]))]
    return res

# Random elements of prime ideals (or their inverses) in FB according to coefficients in alpha.
# Parameters:
#   FB is a set of prime ideals of the field K = field.field(d).sage() (Sage).
#   alpha is an integer vector.
# Returns a set with coefficients in the field.field(d).
def prime_ideals_elements(d, alpha, FB, num = 2):
    res = [[]]*len(alpha)
    for i in range(len(alpha)):
        if alpha[i] >= 0:
            res[i] = [FB[i].random_element() for j in range(num)]
        elif alpha[i] < 0:
            res[i] = [clgp.prime_ideal_inverse(FB[i]).random_element() for j in range(num)]
        res[i] = [field.field(d).from_sage(res[i][j]) for j in range(len(res[i]))]
    return res

# Generates random vector e of length l with elements in [0, bound] and s.t. sum(e) == bound.
def random_vec_bounded(l, bound):
    e = [0]*l
    b = bound
    while b != 0:
        v = ZZ.random_element(0, b+1)
        i = ZZ.random_element(0, l)
        e[i] += v
        b -= v
    return e

# Selecting upto num random generators from smooth ideal h prod P[i]^alpha[i], where P[i] are prime ideals from FB.
# Parameters:
#   K is a number field (Sage)
#   FB is a factor base (set of prime ideals from K)
#   alpha is a vector of integers. The function returns generator of P[i]^(-1) if alpha[i] < 0 and generator of P[i] otherwise.  
#   h is an element from field.field(d)
# Returns values in power products representation of elements from field.field(d).
def smooth_ideal_generators_random(d, h, alpha, FB, num):
    gamma = prime_ideals_generators(d, alpha, FB)

    #print(f"gamma = {[[gamma[i][j].to_sage() for j in range(len(gamma[i]))] for i in range(len(gamma))]}")

    KC = field_compact.field_compact(d)
    T = []
    for i in range(num):
        elems = [h]
        pows = [1]
        for j in range(len(FB)):
            if alpha[j] == 0:
                continue
            # Generating set of prod P[j]^alpha[j] consist of elements prod_k gamma[j][k]^delta_k s.t.
            # beta_k >= 0 and sum_k delta_k = alpha[j]. So, we select elements from generating set with such property.
            delta = random_vec_bounded(len(gamma[j]), abs(alpha[j]))
            assert sum(delta) == abs(alpha[j])

            for k in range(len(gamma[j])):
                if delta[k] != 0:
                    elems.append(gamma[j][k])
                    pows.append(delta[k])
        #print(f"elems={elems}")
        #print(f"pows={pows}")
        gen = KC(elems, pows)
        #print(f"gen = {gen.to_sage()}")
        
        if not gen in T: 
            T.append(gen)
    return T

# Selecting generators of the power of an ideal given by set of generators. 
def ideal_power_generators(gens, e, num = infinity):
    res = []
    if e == 0 or gens == []:
        return res

    res += [gens[0]^e]
    for j in range(e):
        res += [gens[0]^ZZ(j) * g for g in ideal_power_generators(gens[1:], e-j, num)]
        if len(res) > num:
            return res
    return res

# Finding generators of the product of ideals given by generators.
def ideal_product_generators(ideals, num = infinity):
    res = []
    if ideals == []:
        return res
    I = ideals[0]
    if len(ideals) == 1:
        return I
    for j in range(len(I)):
        res += [I[j] * g for g in ideal_product_generators(ideals[1:], num)]
        if len(res) > num:
            return res
    return res

# Selecting upto num generators from smooth ideal h prod P[i]^alpha[i], where P[i] are prime ideals from FB.
# Parameters:
#   K is a number field (Sage)
#   FB is a factor base (set of prime ideals from K)
#   alpha is a vector of integers. The function returns generator of P[i]^(-1) if alpha[i] < 0 and generator of P[i] otherwise.  
#   h is an element from field.field(d)
# Returns values in power products representation of elements from field.field(d).
def smooth_ideal_generators(d, h, alpha, FB, num, random_elements=False):
    K = field.field(d)
    KC = field_compact.field_compact(d)

    gamma = prime_ideals_generators(d, alpha, FB)

    #print(f"gamma = {[[gamma[i][j].to_sage() for j in range(len(gamma[i]))] for i in range(len(gamma))]}")

    ideals = []
    for i in range(len(alpha)):
        if alpha[i] == 0:
            continue
        gens = [KC(g) for g in gamma[i]]
        #print(f"gens = {gens}")
        I = ideal_power_generators(gens, abs(alpha[i]), num)
        ideals.append(I)
    res = ideal_product_generators(ideals, num)
    if h != K.one():
        res = [KC(h) * g for g in res]
    # if len(res) < num and num != infinity and random_elements:
    #     gamma2 = prime_ideals_elements(d, alpha, FB)
    #     print(f"gamma2 = {[[gamma2[i][j].to_sage() for j in range(len(gamma2[i]))] for i in range(len(gamma2))]}")
    #     ideals = []
    #     for i in range(len(alpha)):
    #         if alpha[i] == 0:
    #             continue
    #         gens = [KC(g) for g in gamma2[i]]
    #         #print(f"gens = {gens}")
    #         I = ideal_power_generators(gens, abs(alpha[i]), num)
    #         ideals.append(I)
    #     res += ideal_product_generators(ideals, num)
    return res

# Computation of square root of smooth ideal I^2 = h prod P[i]^alpha[i] where P[i] are prime ideals from the factor base.
# Requires precomputed relations (run testrelations.sage with the same 'food').
def smooth_ideal_sqrt(d, I, h, alpha, d_parent=(), food=None):
    M,B,U,V,U_inv,V_inv = clgp.get_matrices(d)

    clgp.check_clgp(d, d_parent, food=food)

    # Write the input ideal J in terms of class group generators
    e_g = clgp.primes_to_gens(alpha, V, B, reduce=True)
    #print(f"-> Input ideal exponents in terms of Cl_K gens, e_g = {e_g}")

    # check that ideal is written correctly in terms of generators of class group.
    assert clgp.check_prods_p2g(d, alpha, e_g, d_parent, food=food), "Wrong conversion prod P_i^e_i =>  prod g_i^et_i"

    FB = clgp.get_FB(d, d_parent = d_parent, food = food)
    K = field.field(d).sage(food)
    FB_sage = [FB[i].to_sage(food) for i in range(len(FB))]

    assert len(alpha) == len(FB)
    II = I.to_sage(K)^2
    assert II == K.ideal(h.to_sage(K)) * fb.prod_ideal(K, FB, alpha)

    if set(alpha) == {0}:
        h_sqrt = ideal.idealsqrtshorten(len(d), 2^len(d), d, h)
        assert h_sqrt.to_sage(K) in I.to_sage(K)
        return (h_sqrt, alpha)

    #h_sqrt = ideal.idealsqrtshorten(len(d), 2^len(d), d, h)
    #print(f"h_sqrt = {h_sqrt}")

    num = 2^(len(d)) # heuristic number of elements of I^2 s.t. after recognizing the square roots we obtain generators of I
    num = 2^(5*len(d))

    print(f"alpha = {alpha}")

    # Choosing generating set for I^2 = h * prod_i P ^ e_i
    #T = smooth_ideal_generators_random(d, h, alpha, FB_sage, num)
    T = smooth_ideal_generators(d, h, alpha, FB_sage, num)
    print(f"T = {[t.evaluate().to_sage(K) for t in T]}")

    for t in T:
        print(f"t = {t.to_sage(K)}, in I^2? {t.evaluate().to_sage(K) in II}, sq? {t.evaluate().to_sage(K).is_square()}, N(t) = {t.absnorm().factor()}, {[t.valuation(P) for P in FB_sage]}")
        assert t.evaluate().to_sage(K) in II

    II_gens = [t.evaluate().to_sage(K) for t in T]
    I2 = K.ideal(II_gens)
    assert I2 == II

    print(f"I^2 = {II}")

    U = units.generators(d)
    for u in U:
        print(f"u = {u.element.to_sage(K)}")
        assert u.element.to_sage(K).absolute_norm() in [1,-1]
    KC = field_compact.field_compact(d)
    U = [KC(-field.field(d).one()), KC.one()] + [KC(u.element) for u in U]
    #U = [KC(u.element) for u in U]

    #print("h-u squares:", [sq.to_sage() for sq in saturation.find_squares(d, [KC(h)] + U, len(FB))])
    #print("->", KC(h).factor(FB_sage))

    #print("h-u-T squares:", [sq.to_sage() for sq in saturation.find_squares(d, [KC(h)] + U + [t / KC(h) for t in T], len(FB))])

    #print("u-T squares:", [sq.to_sage() for sq in saturation.find_squares(d, U + [t / KC(h) for t in T], len(FB))])

    #print("h-T squares:", [sq.to_sage() for sq in saturation.find_squares(d, [KC(h)] + [t / KC(h) for t in T], len(FB))])

    #E = [KC(field.field(d).from_sage(II.random_element())) for i in range(10)]
    
    # Add some random elements of I^2.
    E = []
    if False:
        for i in range(10):
            #beta = [ZZ.random_element(0,2) for i in range(len(T))]
            beta = [ZZ.random_element(0,2) for i in range(len(FB_sage))]
            # select coprime elements
            #for j in range(len(beta)):
            #    if alpha[j] == 0:
            #        beta[j] = 0
            for j in range(len(beta)):
                #if alpha[j] < 0:
                if alpha[j] < 0:
                    beta[j] = 0
            gamma = prime_ideals_generators(d, beta, FB_sage)
            for j in range(len(beta)):
                if beta[j] == 0:
                    continue
                k = ZZ.random_element(0, len(T))
                l = ZZ.random_element(0, len(gamma[j]))
                E.append(T[k] * gamma[j][l])

    #II_gens = [t.evaluate().to_sage(K) for t in T+E]
    #assert K.ideal(II_gens) == II, "Wrong random elements selected"

        for e in E:
            print(f"e = {e.to_sage(K)}")

    # list of S-units (SU) and the list of powers (SUF), s.t. s O_K = prod P_i^e_i for e in SUF, s in S.
    SU,SUF = sunits.load(d) # SUF is currently broken.

    SU0 = []
    print("Checking S-units ...")
    for i in range(len(SU)):
        su = SU[i]
        e = SUF[i]
        #print(f"su = {su.to_sage(K)}, e = {e}, {[su.valuation(P) for P in FB_sage]}, N(su) = {su.absnorm().factor()}")
        print(f"su = {su.to_sage(K)}, N(su) = {su.absnorm().factor()}")
        assert fb.is_smooth(su.absnorm(), FB)
        try:
            su = su.compute_roots(ideal_sqrt=True)
            #su = KC(su.evaluate(ideal_sqrt=True))
        except Exception as e:
            print(f"-> broken S-unit: {e}")
            continue
        assert SU[i].absnorm() == su.absnorm()
        SU0.append(su)
        print(f"su (new) = {su.to_sage(K)}, N(su) = {su.absnorm().factor()}")
        #suO = K.ideal(su.evaluate().to_sage(K))
        #print(f"suO = {suO}")
        #JJ = fb.prod_ideal(K, FB, e)
        #assert suO == JJ
    SU = SU0

    print(f"S-units left: {len(SU)}")

    #T = T + U
    #v = saturation.find_squares(d, T + U, len(FB))
    #v = saturation.find_squares(d, T + U + E, len(FB))
    #v = saturation.find_squares(d, T + U + E + SU, len(FB))
    v = saturation.find_squares(d, T + U + E + SU, len(FB))
    print(f"squares = {[sq.to_sage() for sq in v]}")

    #v = list(set(v)) # removing duplicates

    for sq in v:
        #print(f"{sq.evaluate().to_sage(K)}, is square? {sq.evaluate().to_sage(K).is_square()}")
        assert sq.evaluate().to_sage(K).is_square()

    #b = [sq.sqrt() for sq in v]
    b = [sq^(1/2) for sq in v]
    #b = [bi for bi in b if abs(bi.absnorm()) != 1]

    if True:
        print("\nLooking for (ideal) squares among generators of I^2 ...")
        print(f"T = {[t.to_sage(K) for t in T]}")
        for t in T:
            try:
                t_r = t.sqrt(ideal_sqrt = True)
                print(f"t = {t.to_sage(K)} = {[t.valuation(P) for P in FB_sage]}")
                print(f"sqrt(t) = {t_r.to_sage(K)} = {(t^(1/2)).to_sage(K)}, {[t_r.valuation(P) for P in FB_sage]}")
                print(f"N(sqrt(t)) = {t_r.absnorm().factor()}")
                assert ((t_r^2) / t).absnorm() in [1,-1]
                if t_r not in b:
                    #b.append(t_r)
                    b.append(t^(1/2))
                    print("-> new")
                    print(f"sqrt(t) in T: {t_r in T}")
                    print(f"t^1/2 in T: {t^(1/2) in T}")
                    #print(f"t.symbols2_log = {t.symbols2_log(0,len(FB)+64)}")
                    #print(f"t.symbols_log = {t.symbols_log(0,len(FB)+64)}")
                else:
                    print("-> old")
            except Exception as e:
                #print(e)
                continue

    b0 = []
    for bi in b:
         if abs(bi.absnorm()) == 1:
             continue
    #     if abs(bi.absnorm().denominator()) != 1:
    #         continue
    #     print(f"bi = {bi.absnorm().factor()}")
    #     if Mod(abs(bi.absnorm()), I.q) != 0:
    #         continue
         if bi in T:
             continue
         if fb.is_smooth(bi.absnorm(), FB):
             continue
         if bi.evaluate().to_sage(K) in I2:
             continue
         b0.append(bi)
    b = b0

    print(f"square roots (filtered): {[bi.to_sage(K) for bi in b]}")
    #print(f"square roots: {[bi.evaluate().to_sage(K) for bi in b]}")

    print(f"I = {I.to_sage(K)}, principal? {I.to_sage(K).is_principal(proof=False)}")
    for r in b:
        r0 = r.evaluate().to_sage(K)
        print(f"r = {r0}")
        print(f"N(r) = {r0.absolute_norm().factor()}")
        print(f"r in I: {r0 in I.to_sage(K)}")
        print(f"r in I^2: {r0 in II}")
        print(f"r.is_square? {r0.is_square()}")
        print(f"r in T+U? {r in T + U}")

    print(f"N(I) = ", I.to_sage(K).absolute_norm().factor())
    print(f"N(I^2) = ", II.absolute_norm().factor())
    print(f"N(h) = {h.to_sage(K).absolute_norm().factor()}")
    print(f"h = {h.to_sage(K)}, square? {h.to_sage(K).is_square()}, powers: {[h.valuation(P) for P in FB_sage]}")

    I_gens = [t.evaluate().to_sage(K) for t in b] + II_gens
    assert K.ideal(I_gens) == I.to_sage(K), f"Wrong ideal I generators, field: {d}"

    # Obtaining (partial) common factor
    g = b[0]
    for i in range(1,len(b)):
        g = b[i].pseudo_gcd(g)
    print(f"gcd(b) = {g.to_sage(K)}")

    for i in range(0,len(T)):
        g = T[i].pseudo_gcd(g)
    print(f"gcd(b + T) = {g.to_sage(K)}")

    for i in range(len(g.elements)):
        if g.elements[i] == h:
            e0 = g.powers[i]
            break
    print(f"e0 = {e0}")
    assert e0 != 0
    #print(f"g = {g.evaluate().to_sage(K)}")
    #print(f"N(g) = {g.evaluate().to_sage(K).absolute_norm().factor()}")
    print(f"g = {g.to_sage(K)}")
    print(f"N(g) = {g.absnorm().factor()}")

    # Obtaining prime decomposition of remaining smooth part
    e_sqrt = [infinity]*len(alpha)
    for i in range(len(FB_sage)):
        for t in b+T:
            e_sqrt[i] = min(e_sqrt[i], (t/g).valuation(FB_sage[i]))
    print(f"e_sqrt = {e_sqrt}")

    gO = K.ideal(g.evaluate().to_sage(K)) 
    print(f"gO = {gO}")
    I2 = gO * fb.prod_ideal(K, FB, e_sqrt)
    print(f"I2 = {I2}")
    assert I.to_sage(K) == I2

    #TODO: reduce e_sqrt in Log-S-unit lattice.

    return (g.evaluate(ideal_sqrt=True), e_sqrt)
