/** 
 *  \file codac2_SampledTrajectory.h
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <map>
#include "codac2_TrajectoryBase.h"
#include "codac2_analytic_variables.h"

namespace codac2
{
  template<typename T>
  class SampledTrajectory : public TrajectoryBase<T>, public std::map<double,T>
  {
    public:

      SampledTrajectory()
        : TrajectoryBase<T>(), std::map<double,T>()
      { }

      SampledTrajectory(const std::map<double,T>& m)
        : TrajectoryBase<T>(), std::map<double,T>(m)
      { }

      // size is not the std::map<double,T>::size() !
      virtual Index size() const
      {
        if constexpr(std::is_same_v<typename ValueType<T>::Type,ScalarType>)
          return 1;

        else
        {
          if(this->empty())
            return 0;
          else
            return this->begin()->second.size();
        }
      }

      size_t nb_samples() const
      {
        return std::map<double,T>::size();
      }

      virtual bool is_empty() const
      {
        return std::map<double,T>::empty();
      }

      virtual Interval tdomain() const
      {
        if(this->empty())
          return Interval::empty();
        else
          return { this->begin()->first, std::prev(this->end())->first };
      }

      virtual void truncate_tdomain(const Interval& new_tdomain)
      {
        assert_release(this->tdomain().is_superset(new_tdomain));

        // Interpolation on the limits of the new tdomain
        T y_lb = (*this)(new_tdomain.lb());
        T y_ub = (*this)(new_tdomain.ub());

        auto it = this->begin();
        while(it != this->end())
        {
          if(!new_tdomain.contains(it->first))
            it = this->erase(it);
          else
            ++it;
        }

        (*this)[new_tdomain.lb()] = y_lb; // clean truncation
        (*this)[new_tdomain.ub()] = y_ub;
      }

      virtual typename Wrapper<T>::Domain codomain() const
      {
        typename Wrapper<T>::Domain hull(this->begin()->second);
        for(const auto& [t,v] : *this)
          hull |= v;
        return hull;
      }

      virtual T operator()(double t) const
      {
        if(!this->tdomain().contains(t))
          return this->nan_value();

        auto it_lower = this->lower_bound(t);
        if(it_lower->first == t)
          return it_lower->second;

        auto it_upper = it_lower;
        it_lower--;

        // Linear interpolation
        return it_lower->second +
               (t - it_lower->first) * (it_upper->second - it_lower->second) /
               (it_upper->first - it_lower->first);
      }

      virtual typename Wrapper<T>::Domain operator()(const Interval& t) const
      {
        // We obtain the output dimension by an evalution...
        typename Wrapper<T>::Domain hull(this->begin()->second);

        if(!this->tdomain().is_superset(t))
          return hull.init(Interval(-oo,oo));

        else
        {
          hull = (*this)(t.lb());
          for(auto it = this->lower_bound(t.lb()) ; it != this->upper_bound(t.ub()) ; it++)
            hull |= it->second;
          hull |= (*this)(t.ub());
          return hull;
        }
      }

      virtual SampledTrajectory<T> sampled(double dt) const
      {
        return sampled(dt, true);
      }

      SampledTrajectory<T> sampled(double dt, bool keep_original_values) const
      {
        assert(dt > 0.);
        assert(!is_empty());

        auto straj = TrajectoryBase<T>::sampled(dt);

        if(keep_original_values)
        {
          // Appending values from the initial map:
          for(const auto& [ti,xi] : *this)
            straj[ti] = xi;
        }
        
        return straj;
      }
  };
}