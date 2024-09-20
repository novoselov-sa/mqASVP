# This is a hack to run unittests from parent directory
import os
import sys
sys.path.append(os.path.dirname(__file__) + "/../")

import unittest
import field
import field_compact
import smooth_ideal_sqrt_sat
import fb

class TestSmoothIdealSqrtSatMethods(unittest.TestCase):

  def setUp(self):
    self.fields = [
      (5, 13),
      (-11, -31),
      (-3, -7, -11),
      (-11, -31, -19),
      (-7, -11, -19, -23),
      (-11, -19, -31, -47),
      (-3, 5*17),
      (7*11, 23*43),
      (3*5, 7*13, 3*13)
    ]

  def test_find_squares_mult(self):
    for d in self.fields:
      KC = field_compact.field_compact(d)
      K = KC.base_field()
      r = 10
      S = [KC.random() for i in range(r)]
      h = K.random()
      squares = smooth_ideal_sqrt_sat.find_squares_mult(d, h, S, r)
      for sq in squares:
        sq_sage = sq.evaluate().to_sage(K.sage()) * h.to_sage(K.sage())
        self.assertTrue(sq_sage.is_square(), f"element should be a square: {sq_sage}")

  def test_smooth_ideal_generators(self):
    for d in self.fields:
      K = field.field(d)
      KC = field_compact.field_compact(d)

      FB = []
      for p in primes(10):
        for P,e in K.sage().factor(p):
          FB.append(P)
  
      h = K.random()
      alpha = [ZZ.random_element(-10, 10) for i in range(len(FB))]
      I = K.sage().ideal(h.to_sage(K.sage()))
      I = prod([FB[i]^alpha[i] for i in range(len(alpha))], I)

      gens = smooth_ideal_sqrt_sat.smooth_ideal_generators(d, h, alpha, FB, 10)
      for a in gens:
        self.assertIn(a.evaluate().to_sage(K.sage()), I)

if __name__ == '__main__':
  unittest.main()
