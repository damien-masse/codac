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

  def tests_geometry(self):
    
    # Align points
    self.assertTrue(orientation([0.,0.],[5.,0.],[8.,0.]) == OrientationInterval.COLINEAR)
    self.assertTrue(orientation([0.,0.],[5.,0.],[-2.,0.]) == OrientationInterval.COLINEAR)
    self.assertTrue(orientation([0.,0.],[-2.,0.],[5.,0.]) == OrientationInterval.COLINEAR)
    self.assertTrue(orientation([0.,0.],[0.,0.],[5.,0.]) == OrientationInterval.COLINEAR)
    self.assertTrue(orientation([0.,0.],[0.,0.],[0.,0.]) == OrientationInterval.COLINEAR)

    # Clockwise points
    self.assertTrue(orientation([0.,0.],[0.,-5.],[-2.,-5.]) == OrientationInterval.CLOCKWISE)
    self.assertTrue(orientation([-999.,0.1],[0.,0.],[-999.,-0.1]) == OrientationInterval.CLOCKWISE)

    # Anticlockwise points
    self.assertTrue(orientation([0.,0.],[0.,-5.],[2.,-5.]) == OrientationInterval.COUNTERCLOCKWISE)
    self.assertTrue(orientation([999.,0.1],[0.,0.],[999.,-0.1]) == OrientationInterval.COUNTERCLOCKWISE)

    # Alignments
    p1 = IntervalVector([0.,0.])
    p2 = IntervalVector([0.,1.])
    p3 = IntervalVector([0.,10.])
    p4 = IntervalVector([1.,10.])
    p5 = IntervalVector([0.,9.+1e-8*Interval(-1.,1.)])
    p6 = IntervalVector([0.+1e-8*Interval(-1.,1.),9.])
    p7 = IntervalVector([0.+1e-8*Interval(-1.,1.),9.+1e-8*Interval(-1.,1.)])

    self.assertTrue(aligned(p1,p2,p3) == BoolInterval.TRUE)
    self.assertTrue(aligned(p1,p2,p4) == BoolInterval.FALSE)
    self.assertTrue(aligned(p1,p1,p4) == BoolInterval.TRUE)
    self.assertTrue(aligned(p1,p2,p5) == BoolInterval.TRUE)
    self.assertTrue(aligned(p1,p2,p6) == BoolInterval.UNKNOWN)
    self.assertTrue(aligned(p1,p2,p7) == BoolInterval.UNKNOWN)
    self.assertTrue(aligned([0,0],[1,1],[10,10]) == BoolInterval.TRUE)
    self.assertTrue(aligned([[-oo,0],[-oo,0]],[1,1],[10,10]) == BoolInterval.UNKNOWN)
    self.assertTrue(aligned(IntervalVector(2),[1,1],[10,10]) == BoolInterval.UNKNOWN)

if __name__ ==  '__main__':
  unittest.main()