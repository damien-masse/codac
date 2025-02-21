#!/usr/bin/env python

#  Codac tests
# ----------------------------------------------------------------------------
#  \date       2024
#  \author     Simon Rohou
#  \copyright  Copyright 2024 Codac Team
#  \license    GNU Lesser General Public License (LGPL)

import unittest
from codac import *
import sys

class TestTransformations(unittest.TestCase):

  def tests_affine_transformation(self):

    t = ScalarVar()
    f1 = AnalyticFunction([t], [ cos(2*t)+0.1*cos(10*t), 2*sin(t)+0.1*sin(10*t) ])

    a = PI + PI/3 ; b = 2.5
    T = Vector([ -3.2, 2.5 ])
    f2 = AnalyticFunction([t], [
      b*cos(a)*f1(t)[0]-b*sin(a)*f1(t)[1] + T[0] + 0.05*cos(100*t),
      b*sin(a)*f1(t)[0]+b*cos(a)*f1(t)[1] + T[1] + 0.05*sin(100*t)
    ])

    src = AnalyticTraj(f1,[-1,3]).sampled(0.01)
    dst = AnalyticTraj(f2,[-1,3]).sampled(0.01)

    tr = affine_transformation(src, dst)
    src2 = SampledVectorTraj()
    for ti,src_i in src:
      src2.set(ti, tr*src_i)

    scale = tr.linear().col(0).norm()
    rotation_matrix = tr.linear() / scale
    rotation_angle = atan2(rotation_matrix(1, 0), rotation_matrix(0, 0))
    translation = tr.translation()

    self.assertTrue(Approx(scale,1e-3) == b)
    self.assertTrue(Approx(rotation_angle,1e-3) == (a-2*PI))
    self.assertTrue(Approx(Vector(translation),1e-3) == T)
      
if __name__ ==  '__main__':
  unittest.main()