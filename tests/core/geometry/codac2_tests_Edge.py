#!/usr/bin/env python

#  Codac tests
# ----------------------------------------------------------------------------
#  \date       2024
#  \author     Simon Rohou
#  \copyright  Copyright 2024 Codac Team
#  \license    GNU Lesser General Public License (LGPL)

import unittest
from codac import *

class TestEdge(unittest.TestCase):

  def tests_intersects(self):

    self.assertTrue(Edge([[0,0],[10,0]]).intersects(Edge([[4,0],[6,0]])) == BoolInterval.TRUE)
    self.assertTrue(Edge([[0,0],[10,0]]).intersects(Edge([[5,0],[15,0]])) == BoolInterval.TRUE)
    self.assertTrue(Edge([[0,0],[10,0]]).intersects(Edge([[5,-5],[5,10]])) == BoolInterval.TRUE)
    self.assertTrue(Edge([[0,0],[10,0]]).intersects(Edge([[15,-5],[15,10]])) == BoolInterval.FALSE)
    self.assertTrue(Edge([[0,0],[10,10]]).intersects(Edge([[1,1],[9,9]])) == BoolInterval.TRUE)
    self.assertTrue(Edge([[1,1],[9,9]]).intersects(Edge([[0,0],[10,10]])) == BoolInterval.TRUE)

    self.assertTrue(Edge(Vector([next_float(-oo),2]),Vector([0,2])).intersects(Edge(Vector([0,0]),Vector([0,1]))) == BoolInterval.FALSE)
    self.assertTrue(Edge(Vector([next_float(-oo),2]),Vector([0,2])).intersects(Edge(Vector([0,0]),Vector([0,2]))) == BoolInterval.TRUE)

  def tests_intersects_bis(self):

    p1,p2,p3 = Vector([0.,0.]), Vector([1.5,3.]), Vector([0.,2.])
    e1 = Edge(Vector([-10.,-10.]), Vector([10.,10.]))
    e2 = Edge(Vector([1.,1.]), Vector([10.,10.]))
    e3 = Edge(Vector([0.,0.]), Vector([2.,4.]))

    self.assertTrue(e1.contains(p1) != BoolInterval.FALSE)
    self.assertTrue(e2.contains(p1) == BoolInterval.FALSE)
    self.assertTrue(e3.contains(p2) != BoolInterval.FALSE)
    self.assertTrue(e3.contains(p3) == BoolInterval.FALSE)

  def tests_proj_intersection(self):

    self.assertTrue(proj_intersection([[0,0],[1,1]], [[3,3],[4,4]])
      == IntervalVector(2))

    self.assertTrue(proj_intersection([[0,0],[1,0]], [[0,2],[1,2]])
      == IntervalVector.empty(2))

    self.assertTrue(proj_intersection([[2,0],[6,4]],[[6,5],[5,6]])
      == Approx(IntervalVector([6.5,4.5])))

    self.assertTrue(proj_intersection([[-0.7,0.7],[0,1]], [[-1,0],[-0.7,-0.7]])
      == Approx(IntervalVector([-1.2069,0.4828]), 1e-4))

    self.assertTrue(proj_intersection([[1,2],[2,3]], [[8,4],[7,4]])
      == IntervalVector([3,4]))

    self.assertTrue(proj_intersection([[4,11],[8,11]], [[12,12],[12,16]])
      == IntervalVector([12,11]))

    self.assertTrue(proj_intersection([[1,2],[1,2]], [[8,4],[7,4]])
      == IntervalVector(2))

    self.assertTrue(proj_intersection([[0,1],[0,2]], [[2,1],[2,2]])
      == IntervalVector.empty(2))

    self.assertTrue(proj_intersection([[-4,4],[0,0]], [[0,0],[4,4]])
      == IntervalVector([0,0]))

    self.assertTrue(proj_intersection([[-4,4],[0,0]], [[4,4],[0,0]])
      == IntervalVector([0,0]))

    self.assertTrue(proj_intersection([[0,4],[4,0]], [[4,4],[0,0]])
      == IntervalVector([2,2]))


if __name__ ==  '__main__':
  unittest.main()