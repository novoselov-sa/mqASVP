# This is a hack to run unittests from parent directory
import os
import sys
sys.path.append(os.path.dirname(__file__) + "/../")

import unittest
import field
import units
import ring
import random as rnd

class TestFieldMethods(unittest.TestCase):
  def setUp(self):
    self.fields = [
      (5, 13),
      (-3, 35),
      (-11, -31),
      (-3, -7, -11),
      (-11, -19, -31),
      (-7, -11, -19, -23),
      (-11, -19, -31, -47),
      (-3, 5*17),
      (7*11, 23*43),
      (3*5, 7*13, 3*13)
    ]

  def test_absolute_polynomial(self):
    for d in self.fields:
      K = field.field(d)
      pol = K.absolute_polynomial().to_sage().change_ring(QQ)
      K2 = NumberField([x^2-d[i] for i in range(len(d))], names='a')
      self.assertEqual(K2.absolute_polynomial(), pol)

  def test_field_discriminant(self):
    self.assertEqual(field.field((2, 3, 5)).discriminant(), 3317760000)
    self.assertEqual(field.field((-2, 3, 5)).discriminant(), 3317760000)
    self.assertEqual(field.field((3, 5, 11)).discriminant(), 189747360000)
    self.assertEqual(field.field((5, 13, 17)).discriminant(), 1490902050625)
    self.assertEqual(field.field((-5, -13, 17)).discriminant(), 381670924960000)
    self.assertEqual(field.field((-3, -7, -11)).discriminant(), 2847396321)
    self.assertEqual(field.field((-3*5, -7, -11)).discriminant(), 1779622700625)

  def test_pow(self):
    for d in self.fields:
      K = field.field(d)
      self.assertEqual(K.zero()^0, K.one())
      v = ZZ.random_element(1000)
      a = ZZ.random_element(100)
      self.assertEqual(K.from_ZZ(v)^a, K.from_ZZ(v^a))

  def test_sqrt(self):
    for d in self.fields:
      K = field.field(d)
      a = K.random()
      s = a^2
      self.assertIn(s.sqrt(), [-a, a])
  
  def test_add_int(self):
    for d in self.fields:
      K = field.field(d)
      self.assertEqual(K.zero() + 1, K.one())
      a = ZZ.random_element(1000)
      b = ZZ.random_element(1000)
      self.assertEqual(K.from_ZZ(a) + b, K.from_ZZ(a+b))
      self.assertEqual(K.from_ZZ(a) - b, K.from_ZZ(a-b))

  def test_mul_int(self):
    for d in self.fields:
      K = field.field(d)
      self.assertEqual(K.zero()*2, K.zero())
      self.assertEqual(K.one()*0, K.zero())
      a = ZZ.random_element(1000)
      b = ZZ.random_element(1000)
      self.assertEqual(K.from_ZZ(a) * b, K.from_ZZ(a*b))
  
  def test_fixed_sqrt_mult(self):
    di = [2, -3, 5, 7*11, -13*17, 19*23*47]
    for p in primes(50):
      for a in di:
         for b in di:
            a_s = field.fixed_sqrt(a, p)
            b_s = field.fixed_sqrt(b, p)
            ab_s = field.fixed_sqrt(a*b, p)
            self.assertEqual(a_s*b_s, ab_s, f"choice of square roots should be compatible with multiplication, a = {a}, b = {b}")

  def test_fixed_sqrt_mult_reverse(self):
    di = reversed([2, -3, 5, 7*11, -13*17, 19*23*47])
    for p in primes(50):
      for a in di:
        for b in di:
            a_s = field.fixed_sqrt(a, p)
            b_s = field.fixed_sqrt(b, p)
            ab_s = field.fixed_sqrt(a*b, p)
            self.assertEqual(a_s*b_s, ab_s, f"choice of square roots should be compatible with multiplication, a = {a}, b = {b}")

  def test_from_sage(self):
    for d in self.fields:
      K = field.field(d)

      # Testing conversion to/from number field NumberField([x^2 - d_1, ..., x^2 - d_n])
      a0 = K.random()
      a1 = a0.to_sage(K.sage())
      a1 = K.from_sage(a1)
      self.assertEqual(a0, a1)
  
  def test_reduce_mod_units(self):
    for d in self.fields:
      K = field.field(d)
      U = units.generators_mod_torsion(d)
      u = prod(rnd.choices(U, k=3))

      a0 = K.random()
      a1 = a0 * u.element

      a1_red = a1.reduce_mod_units()

      I1 = K.sage().ideal(a1.to_sage(K.sage()))
      I2 = K.sage().ideal(a1_red.to_sage(K.sage()))

      self.assertEqual(I1, I2)
  
  def test_div(self):
    for d in self.fields:
      K = field.field(d)
      a = K.random()
      b = K.random()
      c = a / b
      self.assertEqual(a, c * b)

  def test_approxlog(self):
    for d in self.fields:
      K = field.field(d)
      a = K.random()
      b = K.random()
      c = a * b
      self.assertLessEqual(vector(c.approxlog()).norm(), vector(a.approxlog()).norm() + vector(b.approxlog()).norm())
      
  def test_signature(self):
    for d in self.fields:
      K = field.field(d)
      self.assertEqual(K.sage().signature(), K.signature())
  
  def test_set_names(self):
    l = "abcdeghijklmnopqrstuvwxyz"
    for d in self.fields:
      K = field.field(d)
      dp = prod(d).prime_factors()
      food = dict([(l[i], dp[i]) for i in range(len(dp))])
      K.set_names(food)
      self.assertNotEqual(K.names(), tuple(), f"d = {d}")
    
    d = (-2, 7, 11*13*17, -17*47)
    K = field.field(d)
    dp = prod(d).prime_factors()
    food = dict([(l[i], dp[i]) for i in range(len(dp))])
    K.set_names(food)
    self.assertEqual(K.names(), ("a", "b", "cde", "eg"))
    
    K.set_names(names=('a', 'b', 'c', 'd'))
    self.assertEqual(K.names(), ("a", "b", "c", "d"))
  
  def test_compact_repr(self):
    for d in self.fields:
      if len(d) >= 4: #slow fields
        continue
      K = field.field(d)
      a0 = K.random()^5
      a1 = a0.compact_repr().evaluate()
      self.assertEqual(a0, a1)
      
      a0 = K.random()^-5
      a1 = a0.compact_repr().evaluate()
      self.assertEqual(a0, a1)

    # test for fixing of "not enough precision error" in idealred
    a0 = field.field((515290009,))(ring.ring((515290009,))((-323097494626260472025871537445757735701422228011507144111403942294467805879, 14233369683024994509981043570936964669348334604104152978928747802969334)),1)
    a1 = a0.compact_repr().evaluate()
    self.assertEqual(a0, a1)

if __name__ == '__main__':
  unittest.main()
