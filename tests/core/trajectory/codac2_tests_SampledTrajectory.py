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

class TestSampledTrajectory(unittest.TestCase):

  def tests_SampledTrajectory(self):

    x = SampledTrajectory({
      0.25: [-0.5,0.5],
      1.: [0,0],
      2.: [1,0],
      3.: [1,1],
      4.: [-1,1],
      5.: [-1,-1],
      6.: [2,-1]
    })

    self.assertTrue(x.tdomain() == Interval(0.25,6))
    self.assertTrue(x.size() == 2)
    self.assertTrue(x.nb_samples() == 7)
    self.assertTrue(not x.is_empty())
    self.assertTrue(x.codomain() == IntervalVector([[-1,2],[-1,1]]))
    self.assertTrue(x(0.25) == Vector([-0.5,0.5]))
    self.assertTrue(x(1.) == Vector([0,0]))
    self.assertTrue(x(6.) == Vector([2,-1]))
    # Interpolations:
    self.assertTrue(Approx(x(0.75)) == Vector([-1/6.,1/6.]))
    self.assertTrue(x(1.5) == Vector([0.5,0]))
    self.assertTrue(x(5.5) == Vector([0.5,-1]))
    self.assertTrue(x(Interval(1,4)) == IntervalVector([[-1,1],[0,1]]))

    x_sampled = x.sampled(0.1)
    self.assertTrue(x_sampled.tdomain() == Interval(0.25,6))
    self.assertTrue(x_sampled.size() == 2)
    self.assertTrue(x_sampled.nb_samples() > 5*x.nb_samples()) # approx
    self.assertTrue(not x_sampled.is_empty())
    self.assertTrue(x_sampled.codomain() == IntervalVector([[-1,2],[-1,1]]))
    self.assertTrue(x_sampled(0.25) == Vector([-0.5,0.5]))
    self.assertTrue(x_sampled(1.) == Vector([0,0]))
    self.assertTrue(x_sampled(6.) == Vector([2,-1]))
    # Interpolations:
    self.assertTrue(Approx(x_sampled(0.75)) == Vector([-1/6.,1/6.]))
    self.assertTrue(Approx(x_sampled(1.5)) == Vector([0.5,0]))
    self.assertTrue(Approx(x_sampled(5.5)) == Vector([0.5,-1]))
    self.assertTrue(x_sampled(Interval(1,4)) == IntervalVector([[-1,1],[0,1]]))

    #DefaultView.set_window_properties([75,75],[700,700])
    #DefaultView.draw_trajectory(x, Color.blue())
    #DefaultView.draw_trajectory(x_sampled, Color.red())

if __name__ ==  '__main__':
  unittest.main()