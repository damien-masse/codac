/** 
 *  \file codac2_CtcDeriv.h
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Simon Rohou
 *  \copyright  Copyright 2025 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include "codac2_Ctc.h"
#include "codac2_TimePropag.h"
#include "codac2_ConvexPolygon.h"
#include "codac2_TDomain.h"

namespace codac2
{
  template<class T>
  class Slice;

  template<typename T>
  class SlicedTube;
  
  class CtcDeriv
  {
    public:

      CtcDeriv(const TimePropag& time_propag = TimePropag::FWD_BWD);

      void restrict_tdomain(const Interval& tdomain)
      {
        _tdomain = tdomain;
      }

      template<typename T>
      inline void contract(Slice<T>& x, const Slice<T>& v) const
        requires std::is_same_v<T,Interval> || std::is_same_v<T,IntervalVector>
      {
        assert(x.size() == v.size());
        assert(x.t0_tf() == v.t0_tf());
        
        T input = x.input_gate();
        T output = x.output_gate();
        T envelope = x.codomain();

        if constexpr(std::is_same_v<T,Interval>)
          contract(x.t0_tf(), envelope, input, output, v.codomain(), _time_propag);

        else if constexpr(std::is_same_v<T,IntervalVector>)
          for(Index i = 0 ; i < x.size() ; i++)
            contract(x.t0_tf(), envelope[i], input[i], output[i], v.codomain()[i], _time_propag);

        auto x_next = x.next_slice();
        if(x_next && x_next->is_gate())
          x_next->set(output, false);

        auto x_prev = x.prev_slice();
        if(x_prev && x_prev->is_gate())
          x_prev->set(input, false);

        x.set(envelope, false);
      }

      template<typename T>
      inline void contract(SlicedTube<T>& x, const SlicedTube<T>& v) const
        requires std::is_same_v<T,Interval> || std::is_same_v<T,IntervalVector>
      {
        assert_release(TDomain::are_same(x.tdomain(), v.tdomain()));
        Interval t = x.tdomain()->t0_tf() & _tdomain;
        auto it_beg = x.tdomain()->sample(t.lb(),false);
        auto it_end = x.tdomain()->sample(t.ub(),false);

        if((_time_propag & TimePropag::FWD) == TimePropag::FWD)
        {
          for(auto it = it_beg ; it != it_end ; it++)
          {
            auto sx = x.slice(it);
            if(!sx->is_gate())
              this->contract(*sx, *v.slice(it));
          }
        }

        if((_time_propag & TimePropag::BWD) == TimePropag::BWD)
        {
          for(auto it = std::prev(it_end) ; it != std::prev(it_beg) ; it--)
          {
            auto sx = x.slice(it);
            if(!sx->is_gate())
              this->contract(*sx, *v.slice(it));
          }
        }
      }

      static ConvexPolygon polygon_slice(
        const Interval& t, const Interval& envelope,
        const Interval& input, const Interval& proj_input,
        const Interval& output, const Interval& proj_output,
        const Interval& v);

    protected:

      static void contract(
        const Interval& t, Interval& envelope,
        Interval& input, Interval& output,
        const Interval& v, const TimePropag& time_propag);

      const TimePropag _time_propag;
      Interval _tdomain;
  };
}