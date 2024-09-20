# Computation of square roots for S-smooth ideal I using saturation.
# Requires precomputed class group.

import field
import field_compact

def find_squares(d, T, r):
    '''
    Recognizing squares among products of elements from a set T.
    '''
    low,high = 0,r + 64
    # The procedure works with symbols_log if T doesn't contain units.
    t0_s = T[0].symbols2_log(low,high) 
    #t0_s = T[0].symbols_log(low,high)
    A = matrix(GF(2), len(T), len(t0_s))
    A[0] = t0_s
    for i in range(1,len(T)):
        A[i] = T[i].symbols2_log(low,high)
        #A[i] = T[i].symbols_log(low,high)
    V = A.left_kernel().basis_matrix().rows()
    K = field.field(d)
    KC = field_compact.field_compact(d)
    squares = []
    for v in V:
        sq = KC(K.one())
        for j in range(len(v)):
            sq *= T[j] ^ ZZ(v[j])
        squares.append(sq)
    return squares

def find_square_exp(h, T, r = 0):
    '''
    Given a list T of field elements in and element of the field h the method finds binary vector e such that h*T[i]^e[i] is a square.

    Important: this method doesn't support rational exponents in compact representation of field elements.
    '''
    low,high = 0,r + 64
    assert len(T) != 0
    v = vector(GF(2), h.symbols2_log(low, high))
    
    A = matrix(GF(2), len(T), len(v))
    for i in range(len(T)):
        A[i] = T[i].symbols2_log(low,high)
    e = A.solve_left(v)
    return e

def find_square_rat_exp(h, T, r = 0, quadchars=False):
    '''
    Given a list T of field elements in and element of the field h the method finds binary vector e such that h*T[i]^e[i] is a square.

    This method supports rational exponents in compact representation of field elements.
    '''
    low,high = 0,r + 64
    assert len(T) != 0
    v = vector(GF(2), h.altsymbols_log(low, high))
    
    A = matrix(GF(2), len(T), len(v))
    for i in range(len(T)):
        A[i] = T[i].altsymbols_log(low,high)
    
    if quadchars:
        K = field.field(h.gens)
        denom = vector(QQ, h.powers).denominator()
        for t in T:
            denom = lcm(denom, vector(QQ, t.powers).denominator())
        v0 = (h^denom).symbols2_log(low, high)
        v = vector(GF(2), list(v)+list(v0))
        
        A0 = matrix(GF(2), len(T), len(v0))
        for i in range(len(T)):
            A0[i] = (T[i]^denom).symbols2_log(low, high)

        v1 = K.from_ZZ(-1).symbols2_log(0, high)
        A0 = A0.stack(matrix(GF(2), v1))
        A = A.stack(vector(GF(2), [0]*A.ncols())) # for -1
        A = A.augment(A0)

    e = A.solve_left(v)
    if quadchars:
        e = e[:-1] # exclude -1
    return e

def altsymbols_log(T, r):
    low,high = 0,r + 64
    t0_s = T[0].altsymbols_log(low,high) 
    #t0_s = T[0].symbols2_log(low,high)
    A = matrix(GF(2), len(T), len(t0_s))
    A[0] = t0_s
    for i in range(1,len(T)):
        A[i] = T[i].altsymbols_log(low,high)
        #A[i] = T[i].symbols2_log(low,high)
    return A

def find_squares_rat(d, T, r):
    '''
    Recognizing squares among products of elements from a set T.
    '''
    A = altsymbols_log(T, r)
    V = A.left_kernel().basis_matrix().rows()
    K = field.field(d)
    KC = field_compact.field_compact(d)
    squares = []
    for v in V:
        sq = KC(K.one())
        for j in range(len(v)):
            sq *= T[j] ^ ZZ(v[j])
        squares.append(sq)
    return squares
