#!/usr/bin/env python

#  Codac tests
# ----------------------------------------------------------------------------
#  \date       2024
#  \author     Damien Massé
#  \copyright  Copyright 2024 Codac Team
#  \license    GNU Lesser General Public License (LGPL)


import unittest
from codac import *

class TestIntFullPivLU(unittest.TestCase):


# bindings not implemented
   def tests_IntFullPivLU_1(self):
      M = Matrix([
      [ 1, -4, 3, 7 ],
      [ 2, 1, -4, 6 ],
      [ 5, 2, 1 , 9 ],
      [ -1, 0, 3, 2 ]
      ])

      LUdec = IntFullPivLU(M)
      self.assertTrue(LUdec.isInjective()==BoolInterval.TRUE)
      self.assertTrue(LUdec.isSurjective()==BoolInterval.TRUE)
      self.assertTrue(LUdec.rank()==Interval(4))
      self.assertTrue((LUdec.determinant()+602).mag()<=1e-10)
      self.assertTrue((LUdec.ReconstructedMatrix()-M).norm().ub()<=1e-10)
      I1 = LUdec.solve(IntervalMatrix.eye(4,4));
      self.assertTrue((I1*M-Matrix.eye(4,4)).norm().ub()<1e-10)

   def tests_IntFullPivLU_2(self):
      M = Matrix([
      [ 1, -4, 6, 7 ],
      [ 2, 1, 3, 6 ],
      [ 5, 2, 8 , 9 ],
      [ -1, 0, -2, 2 ]
      ])

      LUdec = IntFullPivLU(M)
      self.assertTrue(LUdec.isInjective()==BoolInterval.UNKNOWN)
      self.assertTrue(LUdec.isSurjective()==BoolInterval.UNKNOWN)
      self.assertTrue(LUdec.rank()==Interval([3,4]))
      self.assertTrue((LUdec.determinant()).mag()<=1e-10)
      self.assertTrue((LUdec.ReconstructedMatrix()-M).norm().ub()<=1e-10)
      K = LUdec.kernel()
      self.assertTrue(K.cols()==1)
      self.assertTrue((M*K).norm().ub()<1e-10)
      Im = LUdec.image(M)
      self.assertTrue(Im.cols()==3)


if __name__ ==  '__main__':
  unittest.main()
