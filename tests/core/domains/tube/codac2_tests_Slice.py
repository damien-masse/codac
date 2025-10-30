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
import math

class TestSlice(unittest.TestCase):

  def test_slice_inversion(self):

    tdomain = create_tdomain([0,1]);
    x = SlicedTube(tdomain, Interval(0,10));
    sx = x.first_slice();
    x.set([2,3],0);
    x.set([5,6],1);
    self.assertTrue(sx.invert([4,6], 0) == Interval.empty());
    self.assertTrue(sx.invert([2.5,6], 0) == 0.);
    self.assertTrue(sx.invert([0,1], 1) == Interval.empty());
    self.assertTrue(sx.invert([2.5,6], 1) == 1.);
    self.assertTrue(sx.invert([2.5,6], [0.2,0.5]) == Interval(0.2,0.5));
    self.assertTrue(sx.invert([2.5,6], 0.2) == Interval(0.2));

  def test_slice_inversion_unbounded_derivative(self):

    tdomain = create_tdomain([0,1]);
    x = SlicedTube(tdomain, Interval(0,10));
    v = SlicedTube(tdomain, Interval());
    sx = x.first_slice();
    sv = v.first_slice();
    self.assertTrue(sx.invert(5., sv, sx.t0_tf()) == Interval(0,1));
    self.assertTrue(sx.invert(15., sv, sx.t0_tf()) == Interval.empty());

if __name__ ==  '__main__':
  unittest.main()