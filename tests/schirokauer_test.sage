# This is a hack to run unittests from parent directory
import os
import sys
sys.path.append(os.path.dirname(__file__) + "/../")

import unittest
import field
import field_compact
import smooth_ideal_sqrt_sat
import fb
import schirokauer

class TestSchirokauerMethods(unittest.TestCase):

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
  
  def test_schirokauer_map(self):
    l = 2
    for d in self.fields:
      # 2 is unramified if d_i = 1 mod 4 for all i
      #if set([Mod(di, 4) for di in d]) != {1}:
      #  continue
      K = field.field(d)

      OK = K.sage().ring_of_integers()
      while True:
        a = K.sage()(OK.random_element())
        b = K.sage()(OK.random_element())
        if Mod(a.absolute_norm(), l) != 0 and Mod(b.absolute_norm(), l) != 0:
          break
      a = K.from_sage(a)
      b = K.from_sage(b)
      X_a = schirokauer.schirokauer_map(d, a, l, 1)
      X_b = schirokauer.schirokauer_map(d, b, l, 1)
      X_ab = schirokauer.schirokauer_map(d, a*b, l, 1)

      X_a_b = X_a + X_b
      self.assertEqual(X_ab, X_a_b, f"K=field.field({d}); a = {a}; b = {b}; ab = {a*b}")
      

      #print(f"d = {d}")
      for i in range(5):
        found = False
        while not found:
          a = K.random(bits=8)
          b = K.random(bits=8)
          a_num = K(a.numer, 1)
          b_num = K(b.numer, 1)
          if Mod(a_num.absnorm(), l) == 0 or Mod(b_num.absnorm(), l) == 0 or Mod(a.denom, l) == 0 or Mod(b.denom, l) == 0:
            continue
          found = True
        
        #a = a_num
        #b = b_num

        for j in range(1,5):
          X_a = schirokauer.schirokauer_map(d, a, l, j)
          X_b = schirokauer.schirokauer_map(d, b, l, j)
          X_ab = schirokauer.schirokauer_map(d, a*b, l, j)

          X_a_b = X_a + X_b
          self.assertEqual(X_ab, X_a_b, f"K=field.field({d}); a = {a}; b = {b}; ab = {a*b}")

if __name__ == '__main__':
  unittest.main()
