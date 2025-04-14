#!/usr/bin/env python

#  Codac tests
# ----------------------------------------------------------------------------
#  \date       2024
#  \author     Simon Rohou
#  \copyright  Copyright 2024 Codac Team
#  \license    GNU Lesser General Public License (LGPL)

import unittest
from codac import *

class TestGeometry(unittest.TestCase):

  def tests_Edge(self):

    self.assertTrue(Edge([[0,0],[10,0]]).intersects(Edge([[4,0],[6,0]])) == BoolInterval.TRUE)
    self.assertTrue(Edge([[0,0],[10,0]]).intersects(Edge([[5,0],[15,0]])) == BoolInterval.TRUE)
    self.assertTrue(Edge([[0,0],[10,0]]).intersects(Edge([[5,-5],[5,10]])) == BoolInterval.TRUE)
    self.assertTrue(Edge([[0,0],[10,0]]).intersects(Edge([[15,-5],[15,10]])) == BoolInterval.FALSE)
    self.assertTrue(Edge([[0,0],[10,10]]).intersects(Edge([[1,1],[9,9]])) == BoolInterval.TRUE)
    self.assertTrue(Edge([[1,1],[9,9]]).intersects(Edge([[0,0],[10,10]])) == BoolInterval.TRUE)

    self.assertTrue(Edge(Vector([next_float(-oo),2]),Vector([0,2])).intersects(Edge(Vector([0,0]),Vector([0,1]))) == BoolInterval.FALSE)
    self.assertTrue(Edge(Vector([next_float(-oo),2]),Vector([0,2])).intersects(Edge(Vector([0,0]),Vector([0,2]))) == BoolInterval.TRUE)

if __name__ ==  '__main__':
  unittest.main()