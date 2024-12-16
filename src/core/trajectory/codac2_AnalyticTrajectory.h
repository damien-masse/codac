/** 
 *  \file codac2_AnalyticTrajectory.h
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include "codac2_AnalyticFunction.h"
#include "codac2_TrajectoryBase.h"

namespace codac2
{
  template<typename O, typename S = typename O::Scalar>
  class AnalyticTrajectory : public TrajectoryBase<S>, public AnalyticFunction<O>
  {
    public:

      AnalyticTrajectory(const AnalyticFunction<O>& f, const Interval& tdomain)
        : TrajectoryBase<S>(), AnalyticFunction<O>(f), _tdomain(tdomain)
      {
        assert_release(f.args().total_size() == 1 && "domain of f must be 1d");
      }

      virtual Index size() const
      {
        return AnalyticFunction<O>::output_size();
      }

      virtual bool is_empty() const
      {
        return _tdomain.is_empty();
      }

      virtual Interval tdomain() const
      {
        return _tdomain;
      }

      virtual void truncate_tdomain(const Interval& new_tdomain)
      {
        assert_release(this->tdomain().is_superset(new_tdomain));
        _tdomain &= new_tdomain;
      }

      virtual typename Wrapper<S>::Domain codomain() const
      {
        return AnalyticFunction<O>::eval(this->_tdomain);
      }

      virtual S operator()(double t) const
      {
        if(!this->tdomain().contains(t))
          return this->nan_value();
        return AnalyticFunction<O>::eval(t).mid();
      }

      virtual typename Wrapper<S>::Domain operator()(const Interval& t) const
      {
        if(!this->tdomain().is_superset(t))
          return this->nan_value();
        return AnalyticFunction<O>::eval(t);
      }

    protected:

      Interval _tdomain;
  };
}