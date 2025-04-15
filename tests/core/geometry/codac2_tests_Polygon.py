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

  def tests_Polygon(self):

    p1 = Polygon([[3,-1],[3,4],[5,6],[-1,1]])
    self.assertTrue(p1.contains([3.1,3]) == BoolInterval.FALSE)
    self.assertTrue(p1.contains([2.9,3]) == BoolInterval.TRUE)
    self.assertTrue(p1.contains([3,3]) == BoolInterval.TRUE)
    self.assertTrue(p1.contains([0,1]) == BoolInterval.TRUE)
    self.assertTrue(p1.contains([4,1]) == BoolInterval.FALSE)
    self.assertTrue(p1.contains([2,4]) == BoolInterval.FALSE)
    self.assertTrue(p1.contains([2.8,4]) == BoolInterval.TRUE)
    self.assertTrue(p1.contains([3,4]) == BoolInterval.TRUE)
    self.assertTrue(p1.contains([4,4]) == BoolInterval.FALSE)
    self.assertTrue(p1.contains([5,6]) == BoolInterval.TRUE)
    self.assertTrue(p1.contains([6,6]) == BoolInterval.FALSE)

    transect = Edge(Vector([next_float(-oo),3]), Vector([3,3]))
    e1,e2 = Edge(Vector([5,6]),Vector([-1,1])), Edge(Vector([3,-1]),Vector([3,4]))
    self.assertTrue(transect.intersects(e1) == BoolInterval.TRUE)
    self.assertTrue(transect.intersects(e2) == BoolInterval.TRUE)

    p2 = Polygon([[0,0],[0,1],[1,1],[1,0]])
    self.assertTrue(p2.contains([0,0]) == BoolInterval.TRUE)
    self.assertTrue(p2.contains([1,1]) == BoolInterval.TRUE)
    self.assertTrue(p2.contains([0,2]) == BoolInterval.FALSE)
    self.assertTrue(p2.contains([2,0]) == BoolInterval.FALSE)
    self.assertTrue(p2.contains([0.5,1]) == BoolInterval.TRUE)
    self.assertTrue(p2.contains([1,0.5]) == BoolInterval.TRUE)

    transect = Edge(Vector([next_float(-oo),2]), Vector([0,2]))
    e1 = Edge(Vector([0,0]),Vector([0,1]))
    e2 = Edge(Vector([0,1]),Vector([1,1]))
    e3 = Edge(Vector([1,1]),Vector([1,0]))
    e4 = Edge(Vector([1,0]),Vector([0,0]))

    self.assertTrue(transect.intersects(e1) == BoolInterval.FALSE)
    self.assertTrue(transect.intersects(e2) == BoolInterval.FALSE)
    self.assertTrue(transect.intersects(e3) == BoolInterval.FALSE)
    self.assertTrue(transect.intersects(e4) == BoolInterval.FALSE)

    p3 = Polygon([[0,1],[1,0],[0,0]])
    self.assertTrue(p3.contains([1,1]) == BoolInterval.FALSE)

if __name__ ==  '__main__':
  unittest.main()