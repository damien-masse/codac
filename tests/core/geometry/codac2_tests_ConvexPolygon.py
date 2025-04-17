#!/usr/bin/env python

#  Codac tests
# ----------------------------------------------------------------------------
#  \date       2024
#  \author     Simon Rohou
#  \copyright  Copyright 2024 Codac Team
#  \license    GNU Lesser General Public License (LGPL)

import unittest
from codac import *

class TestConvexPolygon(unittest.TestCase):

  def tests_ConvexPolygon_degenerate_cases(self):

    p1 = ConvexPolygon([[0,4]])
    self.assertTrue(len(p1.edges()) == 1)
    self.assertTrue(p1.edges()[0][0] == IntervalVector([0,4]))
    self.assertTrue(p1.edges()[0][1] == IntervalVector([0,4]))

    p1 = ConvexPolygon([[0,4],[2,8]])
    self.assertTrue(len(p1.edges()) == 2)
    self.assertTrue(p1.edges()[0][0] == IntervalVector([2,8]))
    self.assertTrue(p1.edges()[0][1] == IntervalVector([0,4]))
    self.assertTrue(p1.edges()[1][0] == IntervalVector([0,4]))
    self.assertTrue(p1.edges()[1][1] == IntervalVector([2,8]))
    self.assertTrue(p1 == ConvexPolygon([[2,8],[0,4]]))

  def tests_ConvexPolygon_intersection(self):

    p1 = ConvexPolygon([[0,4],[4,8],[7,0]])
    p2 = ConvexPolygon([[-1,2],[-1,3],[5,3],[5,2]])
    print("1 ", p1 & p2)
    self.assertTrue((p1 & p2) == ConvexPolygon([
        [3.5,2],[5,2],[5,3],[1.75,3]
      ]))

    p1 = ConvexPolygon([[3,2],[1,6],[6,5]])
    p2 = ConvexPolygon(IntervalVector([[0,4],[0,4]]))
    print("2 ", p1 & p2)
    self.assertTrue((p1 & p2) == ConvexPolygon([
        [3,2],[2,4],[4,4],[4,3]
      ]))
  
    p1 = ConvexPolygon([[1,2],[3,4],[5,1],[2,1]])
    p2 = ConvexPolygon(IntervalVector([[2,6],[2,6]]))
    print("3 ", p1 & p2)
    self.assertTrue(Approx(p1 & p2) == ConvexPolygon([
        [2,2],[2,3],[3,4],[4.+1./3.,2.]
      ]))

    # Polygons intersections, test 3 (big box)
    p1 = ConvexPolygon([[1,2],[3,4],[5,1],[2,1]])
    p2 = ConvexPolygon(IntervalVector([[-10,10],[-10,10]]))
    print("4 ", p1 & p2)
    self.assertTrue((p1 & p2) == p1) # same polygon

    # Polygons intersections, test 4 (inner box)
    p1 = ConvexPolygon([[1,2],[3,4],[5,1],[2,1]])
    p2 = ConvexPolygon(IntervalVector([[2.8,3],[2.8,3]]))
    print("5 ", p1 & p2)
    self.assertTrue((p1 & p2) == p2) # same box

    # Polygons intersections, test 5
    p1 = ConvexPolygon([[2,1],[3,1],[4,2],[4,3],[3,4],[2,4],[1,3],[1,2]])
    p2 = ConvexPolygon(IntervalVector([[1,4],[1,4]]))
    print("6 ", p1 & p2)
    self.assertTrue((p1 & p2) == p1) # same polygon

    # Polygons intersections, test 6 (shifted polygon points declaration)
    p1 = ConvexPolygon([[3,4],[2,4],[1,3],[1,2],[2,1],[3,1],[4,2],[4,3]])
    p2 = ConvexPolygon(IntervalVector([[1,4],[1,4]]))
    print("7 ", p1 & p2)
    self.assertTrue((p1 & p2) == p1) # same polygon

    # Polygons intersections, test 7 (degenerate case)
    p1 = ConvexPolygon([[4000,200]])
    p2 = ConvexPolygon(IntervalVector([4000,200]))
    print("8 ", p1 & p2)
    self.assertTrue((p1 & p2) == p1) # same polygon

    # Polygons intersections, test 8
    p1 = ConvexPolygon([[1,1],[2,4],[7,5],[6,2]])
    p2 = ConvexPolygon(IntervalVector([[2,6],[1,5]]))
    print("9 ", p1 & p2)
    q = p1 & p2
    self.assertTrue(Approx(q) == Polygon([[2,1.2],[6,2],[6,4.8],[2,4]]))
    self.assertTrue(len(q.edges()) == 4)
    self.assertTrue(len(q.unsorted_vertices()) == 4)

    # Polygons intersections, test 9
    p1 = ConvexPolygon([[1,1],[2,4],[7,5],[6,2]])
    p2 = ConvexPolygon(IntervalVector([[3,5],[1,5]]))
    print("10 ", p1 & p2)
    q = p1 & p2
    self.assertTrue(Approx(q) == Polygon([[3,1.4],[5,1.8],[5,4.6],[3,4.2]]))
    self.assertTrue(len(q.edges()) == 4)
    self.assertTrue(len(q.unsorted_vertices()) == 4)

    # Polygons intersections, test 10 (degenerated box)
    p1 = ConvexPolygon([[1,1],[2,4],[7,5],[6,2]])
    p2 = ConvexPolygon(IntervalVector([[4],[1,5]]))
    print("11 ", p1 & p2)
    q = p1 & p2
    self.assertTrue(Approx(q) == Polygon([[4,4.4],[4,1.6]]))
    self.assertTrue(len(q.unsorted_vertices()) == 2)

    # Polygons intersections, test 10 (degenerated polygon)
    p1 = ConvexPolygon([[1,1],[2,4],[7,5],[6,2]])
    p2 = ConvexPolygon([[4,1],[4,5]])
    print("12 ", p1 & p2)
    q = p1 & p2
    self.assertTrue(Approx(q) == Polygon([[4,4.4],[4,1.6]]))
    self.assertTrue(len(q.unsorted_vertices()) == 2)

if __name__ ==  '__main__':
  unittest.main()