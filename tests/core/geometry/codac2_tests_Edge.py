#!/usr/bin/env python

#  Codac tests
# ----------------------------------------------------------------------------
#  \date       2024
#  \author     Simon Rohou
#  \copyright  Copyright 2024 Codac Team
#  \license    GNU Lesser General Public License (LGPL)

import unittest
from codac import *

def hull(x):
  if x == -oo:
    return Interval(-oo, next_float(-oo))
  if x == oo:
    return Interval(previous_float(oo), oo)
  else:
    return INterval(x)


class TestEdge(unittest.TestCase):

  def tests_intersects(self):

    self.assertTrue(Edge([[0,0],[10,0]]).intersects(Edge([[4,0],[6,0]])) == BoolInterval.TRUE)
    self.assertTrue(Edge([[0,0],[10,0]]).intersects(Edge([[5,0],[15,0]])) == BoolInterval.TRUE)
    self.assertTrue(Edge([[0,0],[10,0]]).intersects(Edge([[5,-5],[5,10]])) == BoolInterval.TRUE)
    self.assertTrue(Edge([[0,0],[10,0]]).intersects(Edge([[15,-5],[15,10]])) == BoolInterval.FALSE)
    self.assertTrue(Edge([[0,0],[10,10]]).intersects(Edge([[1,1],[9,9]])) == BoolInterval.TRUE)
    self.assertTrue(Edge([[1,1],[9,9]]).intersects(Edge([[0,0],[10,10]])) == BoolInterval.TRUE)

    self.assertTrue(Edge(IntervalVector([hull(-oo),2]),Vector([0,2])).intersects(Edge(Vector([0,0]),Vector([0,1]))) == BoolInterval.FALSE)
    self.assertTrue(Edge(IntervalVector([hull(-oo),2]),Vector([0,2])).intersects(Edge(Vector([0,0]),Vector([0,2]))) == BoolInterval.TRUE)

  def tests_contains(self):

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

  def tests_intersection(self):

    self.assertTrue((Edge([0,0],[0,0]) & Edge([0,0],[5,0])) == IntervalVector([0,0])) # degenerate line, horizontal edge line
    self.assertTrue((Edge([0,0],[0,0]) & Edge([0,0],[0,5])) == IntervalVector([0,0])) # degenerate line, vertical edge line

    self.assertTrue((Edge([0,0],[0,1]) & Edge([0,0],[0,5.0])) == IntervalVector([0,[0,1.]])) # vertical edge line (colinear)
    self.assertTrue((Edge([0,0],[0,1]) & Edge([0,0],[0,0.5])) == IntervalVector([0,[0,0.5]])) # vertical edge line (colinear)
    self.assertTrue((Edge([0,0],[1,0]) & Edge([0,0],[5.0,0])) == IntervalVector([[0,1.],0])) # horizontal edge line (colinear)
    self.assertTrue((Edge([0,0],[1,0]) & Edge([0,0],[0.5,0])) == IntervalVector([[0,0.5],0])) # horizontal edge line (colinear)

    self.assertTrue((Edge([7,4],[7,5]) & Edge([7,6],[7,8])) == IntervalVector.empty(2)) # vertical edge line (colinear), no intersection
    self.assertTrue((Edge([4,7],[5,7]) & Edge([6,7],[8,7])) == IntervalVector.empty(2)) # horizontal edge line (colinear), no intersection

    self.assertTrue((Edge([0,0],[3,0]) & Edge([1,2],[1,-999])) == IntervalVector([1,0])) # perpendicular intersection
    self.assertTrue((Edge([0.5,-1],[0.5,5]) & Edge([-3,2],[5,2])) == IntervalVector([0.5,2])) # perpendicular intersection

    self.assertTrue((Edge([0,0],[3,0]) & Edge([1,-10],[1,-999])) == IntervalVector.empty(2)) # perpendicular lines, no intersection
    self.assertTrue((Edge([0.5,-1],[0.5,5]) & Edge([-3,-2],[5,-2])) == IntervalVector.empty(2)) # perpendicular lines, no intersection

    self.assertTrue((Edge([8,4],[9,2]) & Edge([7,3],[9,3])) == IntervalVector([8.5,3])) # perpendicular oblique lines
    self.assertTrue((Edge([8,4],[9,2]) & Edge([8.5,8],[8.5,0])) == IntervalVector([8.5,3])) # perpendicular oblique lines

    self.assertTrue((Edge([8,4],[9,2]) & Edge([7.,3],[7.5,3])) == IntervalVector.empty(2)) # secant oblique lines, no intersection
    self.assertTrue((Edge([8,4],[9,2]) & Edge([8.5,8],[8.5,7])) == IntervalVector.empty(2)) # secant oblique lines, no intersection
    self.assertTrue((Edge([6,-1],[8,1]) & Edge([7.5,0],[9,0])) == IntervalVector.empty(2)) # secant oblique lines, no intersection
    self.assertTrue((Edge([6,-1],[8,1]) & Edge([6.5,0.5],[6.5,2])) == IntervalVector.empty(2)) # secant oblique lines, no intersection
  
    # Other tests
    self.assertTrue((Edge([8,14],[6,13]) & Edge([1,1],[1,14])) == IntervalVector.empty(2))
    self.assertTrue((Edge([hull(-oo),0],[1,0]) & Edge([0,1],[0,-1])) == IntervalVector.zero(2))
    self.assertTrue((Edge([0,0],[1,2]) & Edge([0.5,0],[0.5,2])) == IntervalVector([0.5,1]))
    self.assertTrue((Edge([-1,-1],[0,0]) & Edge([9,3],[0,0])) == IntervalVector([0,0]))
    self.assertTrue((Edge([0,0],[10,0]) & Edge([2,0],[8,0])) == IntervalVector([[2,8],[0]])) # colinear case
    self.assertTrue((Edge([0,0],[10,10]) & Edge([2,2],[8,8])) == IntervalVector([[2,8],[2,8]])) # colinear case
    self.assertTrue((Edge([0,0],[10,0]) & Edge([2,1],[8,1])) == IntervalVector.empty(2)) # parallel case

    self.assertTrue((Edge([2,0],[6,4]) & Edge([6,5],[5,6])) == IntervalVector.empty(2))
    self.assertTrue(proj_intersection(Edge([2,0],[6,4]), Edge([6,5],[5,6])) == IntervalVector([6.5,4.5]))

if __name__ ==  '__main__':
  unittest.main()