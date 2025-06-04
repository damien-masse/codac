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

class TestTDomain(unittest.TestCase):

  def test_tdomain(self):

    tdomain = create_tdomain()
    tdomain.sample(Interval(0,1), 0.5)
    self.assertTrue(tdomain.nb_tslices() == 4)
    self.assertTrue(tdomain.t0_tf() == Interval(-oo,oo))
    
    #vector<TSlice> vector_tslices{
    #  make_move_iterator(tdomain.begin()),
    #  make_move_iterator(tdomain.end()) }
    #  
    #self.assertTrue(vector_tslices.size() == 4)
    #self.assertTrue(vector_tslices[0] == Interval(-oo,0))
    #self.assertTrue(vector_tslices[1] == Interval(0,0.5))
    #self.assertTrue(vector_tslices[2] == Interval(0.5,1))
    #self.assertTrue(vector_tslices[3] == Interval(1,oo))

    #self.assertTrue(*tdomain.iterator_tslice(-10.) == Interval(-oo,0))
    #self.assertTrue(*tdomain.iterator_tslice(-120.) == Interval(-oo,0))
    #self.assertTrue(*tdomain.iterator_tslice(0.2) == Interval(0,0.5))
    #self.assertTrue(*tdomain.iterator_tslice(5540.2) == Interval(1,oo))

    #self.assertTrue(tdomain.nb_tubes() == 0)
    #x = SlicedTube(tdomain, IntervalVector(2))
    #self.assertTrue(tdomain.nb_tubes() == 1)

    #{ # new scope
    #  SlicedTube v(tdomain, IntervalVector(3))
    #  self.assertTrue(tdomain.nb_tubes() == 2)
    #} # end of scope, removing the tube
    
    #self.assertTrue(tdomain.nb_tubes() == 1)
  # }

  # SECTION("Test degenerated TDomain")
  # {
  #   tdomain = create_tdomain(Interval(1))
  #   self.assertTrue(tdomain.nb_tslices() == 1)
  #   self.assertTrue(tdomain.t0_tf() == Interval(1))
  #   self.assertTrue(tdomain.nb_tubes() == 0)
    
  #   vector<TSlice> vector_tslices{
  #     make_move_iterator(tdomain.begin()),
  #     make_move_iterator(tdomain.end()) }
      
  #   self.assertTrue(vector_tslices.size() == 1)
  #   self.assertTrue(vector_tslices[0] == Interval(1))
  # }

  # SECTION("Test TDomain with gates")
  # {
  #   tdomain = create_tdomain(Interval(0,1), 0.5, true)
  #   self.assertTrue(tdomain.nb_tslices() == 5)
  #   self.assertTrue(tdomain.t0_tf() == Interval(0,1))
  #   self.assertTrue(tdomain.nb_tubes() == 0)
    
  #   vector<TSlice> vector_tslices{
  #     make_move_iterator(tdomain.begin()),
  #     make_move_iterator(tdomain.end()) }

  #   self.assertTrue(vector_tslices.size() == 5)
  #   self.assertTrue(vector_tslices[0] == Interval(0))
  #   self.assertTrue(vector_tslices[1] == Interval(0,0.5))
  #   self.assertTrue(vector_tslices[2] == Interval(0.5))
  #   self.assertTrue(vector_tslices[3] == Interval(0.5,1))
  #   self.assertTrue(vector_tslices[4] == Interval(1,1))

  #   self.assertTrue(*tdomain.iterator_tslice(0.) == Interval(0))
  #   self.assertTrue(*tdomain.iterator_tslice(0.1) == Interval(0,0.5))
  #   self.assertTrue(*tdomain.iterator_tslice(0.5) == Interval(0.5))
  #   self.assertTrue(*tdomain.iterator_tslice(0.6) == Interval(0.5,1))
  #   self.assertTrue(*tdomain.iterator_tslice(1.) == Interval(1))
  # }

  # SECTION("Test TDomain with sampling")
  # {
  #   tdomain = create_tdomain()
  #   tdomain.sample(1.)
  #   self.assertTrue(tdomain.nb_tslices() == 2)
  #   tdomain.sample(10.,false) # no gate
  #   self.assertTrue(tdomain.nb_tslices() == 3)
  #   tdomain.sample(10.,true) # second sampling with gate
  #   self.assertTrue(tdomain.nb_tslices() == 4)
  #   tdomain.sample(10.,true) # no more action
  #   self.assertTrue(tdomain.nb_tslices() == 4)
    
  #   vector<TSlice> vector_tslices{
  #     make_move_iterator(tdomain.begin()),
  #     make_move_iterator(tdomain.end()) }

  #   self.assertTrue(vector_tslices.size() == 4)
  #   self.assertTrue(vector_tslices[0] == Interval(-oo,1))
  #   self.assertTrue(vector_tslices[1] == Interval(1,10))
  #   self.assertTrue(vector_tslices[2] == Interval(10))
  #   self.assertTrue(vector_tslices[3] == Interval(10,oo))
  # }

  # SECTION("Test unbounded TDomain")
  # {
  #   tdomain = create_tdomain()
  #   self.assertTrue(tdomain.nb_tslices() == 1)
  #   self.assertTrue(tdomain.t0_tf() == Interval(-oo,oo))
  #   vector<TSlice> vector_tslices{
  #     make_move_iterator(tdomain.begin()),
  #     make_move_iterator(tdomain.end()) }
  #   self.assertTrue(vector_tslices.size() == 1)
  #   self.assertTrue(vector_tslices[0] == Interval(-oo,oo))
  # }
      
if __name__ ==  '__main__':
  unittest.main()