# This is a hack to run unittests from parent directory
import os
import sys
sys.path.append(os.path.dirname(__file__) + "/../")

import unittest
import field
import field_compact
import saturation
import ring

class TestSaturationMethods(unittest.TestCase):

  def setUp(self):
    self.fields = [
      (5, 13),
      (-3, 35),
      (-11, -31),
      (-3, -7, -11),
      (-11, -31, -19),
      (-7, -11, -19, -23),
      (-11, -19, -31, -47),
      (-3, 5*17),
      (7*11, 23*43),
      (3*5, 7*13, 3*13)
    ]

    with open("tests/saturation_test_data.txt", "r") as f:
      data = f.read()
      self.T = eval(preparse(data))

  def test_find_squares(self):
    for d in self.fields:
      KC = field_compact.field_compact(d)
      K = KC.base_field()
      r = 10
      t = [KC.random() for i in range(r)]
      squares = saturation.find_squares(d, t, r)
      for sq in squares:
        sq_sage = sq.evaluate().to_sage(K.sage())
        self.assertTrue(sq_sage.is_square(), f"element should be a square: {sq_sage}")

  def test_find_squares_rat(self):
    d = (-3, 35)
    KC = field_compact.field_compact(d)
    K = KC.base_field()
    r = 10

    T = self.T #self.T[0:16]
    
    T += [KC.from_ZZ(2), KC.from_ZZ(3)]

    #print(f"len(T) = {len(T)}")
    #print(f"{[t.evaluate().to_sage(K.sage()) for t in T]}")
    squares = saturation.find_squares_rat(d, T, r)
    for sq in squares:
        sq_sage = sq.evaluate().to_sage(K.sage())
        self.assertTrue(sq_sage.is_square(), f"element should be a square: {sq_sage}")

if __name__ == '__main__':
  unittest.main()
