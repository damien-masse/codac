/** 
 *  \file
 *  AbstractDomain class
 * ----------------------------------------------------------------------------
 *  \date       2020
 *  \author     Simon Rohou
 *  \copyright  Copyright 2020 Simon Rohou
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */

#ifndef __TUBEX_ABSTRACT_DOMAIN_H__
#define __TUBEX_ABSTRACT_DOMAIN_H__

#include "ibex_Interval.h"
#include "ibex_IntervalVector.h"
#include "tubex_Slice.h"
#include "tubex_Tube.h"
#include "tubex_TubeVector.h"

namespace tubex
{
  enum class DomainType { INTERVAL, INTERVAL_VECTOR, SLICE, TUBE, TUBE_VECTOR };

  enum class DomainRelation
  {
    IN = 0x01,
    OUT = 0x02
  };

  inline int operator&(DomainRelation a, DomainRelation b)
  { return static_cast<int>(static_cast<int>(a) & static_cast<int>(b)); }

  inline DomainRelation operator|(DomainRelation a, DomainRelation b)
  { return static_cast<DomainRelation>(static_cast<int>(a) | static_cast<int>(b)); }
  

  class AbstractContractor;

  class AbstractDomain
  {
    public:

      AbstractDomain(const AbstractDomain& ad);
      explicit AbstractDomain(ibex::Interval& i);
      explicit AbstractDomain(ibex::IntervalVector& iv);
      explicit AbstractDomain(tubex::Slice& s);
      explicit AbstractDomain(tubex::Tube& t);
      explicit AbstractDomain(tubex::TubeVector& tv);
      ~AbstractDomain();

      DomainType type() const;

      double volume() const;
      bool is_empty() const;
      
      bool operator==(const AbstractDomain& x) const;
      bool is_component_of(const AbstractDomain* x) const;
      bool is_prev_slice_of(const AbstractDomain* x) const;
      bool is_next_slice_of(const AbstractDomain* x) const;
      bool is_slice_of(const AbstractDomain* x) const;
      
      friend std::ostream& operator<<(std::ostream& str, const AbstractDomain& x);

      std::vector<AbstractContractor*> m_v_ctc;

    //protected:

      const DomainType m_type;

      ibex::Interval& m_i;
      ibex::IntervalVector& m_iv;
      tubex::Slice& m_s;
      tubex::Tube& m_t;
      tubex::TubeVector& m_tv;

      static ibex::Interval m_dump_i;
      static ibex::IntervalVector m_dump_iv;
      static tubex::Slice m_dump_s;
      static tubex::Tube m_dump_t;
      static tubex::TubeVector m_dump_tv;
  };

  struct DomainParams
  {
    AbstractDomain *ad;
    DomainRelation rel = DomainRelation::IN | DomainRelation::OUT;
    double volume = 0.;
  };
}

#endif