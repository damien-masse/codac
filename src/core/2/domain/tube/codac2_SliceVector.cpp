/** 
 *  \file
 *  
 * ----------------------------------------------------------------------------
 *  \date       2022
 *  \author     Simon Rohou
 *  \copyright  Copyright 2022 Codac Team
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */

#include "codac2_TubeVector.h"
#include "codac2_SliceVector.h"

using namespace std;

namespace codac2
{
  SliceVector::SliceVector(size_t n, const TubeVector& tube_vector, list<TSlice>::iterator it_tslice) :
    _tubevector(tube_vector), _it_tslice(it_tslice),
    _codomain(IntervalVector(n))
  {

  }

  SliceVector::SliceVector(const SliceVector& s) :
    _tubevector(s._tubevector), _it_tslice(s._it_tslice), _codomain(s._codomain)
  {
    
  }

  SliceVector::~SliceVector()
  {

  }

  const TubeVector& SliceVector::tube_vector() const
  {
    return _tubevector;
  }

  size_t SliceVector::size() const
  {
    return _codomain.size();
  }

  bool SliceVector::contains(const TrajectoryVector& value) const
  {
    return true;
  }

  const Interval& SliceVector::t0_tf() const
  {
    return _it_tslice->t0_tf();
  }

  const SliceVector* SliceVector::prev_slice() const
  {
    if(&_tubevector.first_slice() == this)
      return nullptr;
    return &prev(_it_tslice)->slices().at(&_tubevector);
  }

  SliceVector* SliceVector::prev_slice()
  {
    return const_cast<SliceVector*>(
      static_cast<const SliceVector&>(*this).prev_slice());
  }

  const SliceVector* SliceVector::next_slice() const
  {
    if(&_tubevector.last_slice() == this)
      return nullptr;
    return &next(_it_tslice)->slices().at(&_tubevector);
  }

  SliceVector* SliceVector::next_slice()
  {
    return const_cast<SliceVector*>(
      static_cast<const SliceVector&>(*this).next_slice());
  }

  const IntervalVector& SliceVector::codomain() const
  {
    return _codomain;
  }

  IntervalVector SliceVector::input_gate() const
  {
    IntervalVector gate = codomain();
    if(prev_slice())
      gate &= prev_slice()->codomain();
    return gate;
  }

  IntervalVector SliceVector::output_gate() const
  {
    IntervalVector gate = codomain();
    if(next_slice())
      gate &= next_slice()->codomain();
    return gate;
  }

  void SliceVector::set(const IntervalVector& codomain)
  {
    assert((size_t)codomain.size() == size());
    _codomain = codomain;
  }

  ostream& operator<<(ostream& os, const SliceVector& x)
  {
    os << x.t0_tf()
       << "↦" << x.codomain()
       << flush;
    return os;
  }
} // namespace codac