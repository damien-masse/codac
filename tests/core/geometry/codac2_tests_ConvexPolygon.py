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

  def tests_ConvexPolygon(self):

    p1 = ConvexPolygon([[0,4],[4,8],[7,0]])
    p2 = ConvexPolygon([[-1,2],[-1,3],[5,3],[5,2]])

    self.assertTrue((p1 & p2) == ConvexPolygon([
        [3.5,2],[5,2],[5,3],[1.75,3]
      ]))

if __name__ ==  '__main__':
  unittest.main()