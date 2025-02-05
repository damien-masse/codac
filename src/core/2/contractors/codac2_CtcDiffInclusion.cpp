/** 
 *  CtcDiffInclusion class
 * ----------------------------------------------------------------------------
 *  \date       2022
 *  \author     
 *  \copyright  Copyright 2022 Codac Team
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */

#include "codac2_CtcDiffInclusion.h"

using namespace std;
using namespace ibex;

namespace codac2
{
  CtcDiffInclusion::CtcDiffInclusion(const TFunction& f)
    : _f(f)
  {

  }

  const TFunction& CtcDiffInclusion::f() const
  {
    return _f;
  }

  void CtcDiffInclusion::contract(Tube<IntervalVector>& x, const Tube<IntervalVector>& u, TimePropag t_propa)
  {
    // Verifying that x and u share exactly the same tdomain and slicing:
    assert(x.tdomain() == u.tdomain());
    // Verifying that the provided tubes are consistent with the function
    assert((size_t)_f.nb_var() == 2);
    assert((size_t)_f.image_dim() == x.size());

    for(auto& sx : x) // sx is a SliceVector of the TubeVector x
    {
      if(sx.is_gate()) // the slice may be on a degenerated temporal domain, i.e. a gate
        continue;

      // su is a SliceVector of the TubeVector u:
      const shared_ptr<Slice<IntervalVector>> su = static_pointer_cast<Slice<IntervalVector>>(sx.tslice().slices().at(&u));
      
      //const double dt = sx.t0_tf().diam();


      // sx.set(su.codomain());
      // cout << sx << " " << su << endl;

      // ...

      if(t_propa & TimePropag::FORWARD)
      {
        // Computations related to forward propagation
        // ...
      }

      if(t_propa & TimePropag::BACKWARD)
      {
        // Computations related to backward propagation
        // ...
      }
    }
  }

  void CtcDiffInclusion::contract(Slice<IntervalVector>& x, const Slice<IntervalVector>& u, TimePropag t_propa)
  {
    // Verifying that x and u share exactly the same tdomain
    assert(&x.tslice() == &u.tslice());
    // Verifying that the provided slices are consistent with the function
    assert((size_t)_f.nb_var() == 2);
    assert((size_t)_f.image_dim() == x.size());

    //const double dt = x.t0_tf().diam();

    // ...

    if(t_propa & TimePropag::FORWARD)
    {
      // Computations related to forward propagation
      // ...
    }

    if(t_propa & TimePropag::BACKWARD)
    {
      // Computations related to backward propagation
      // ...
    }
  }
}