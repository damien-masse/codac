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

#ifndef __CODAC2_TUBEVECTOR_H__
#define __CODAC2_TUBEVECTOR_H__

#include <list>
#include <memory>
#include "codac2_TDomain.h"
#include "codac2_SliceVector.h"

#include "codac2_TubeAbstract_const.h"

namespace codac2
{
  class TDomain;
  //class SliceVector;

  class TubeVector : public TubeVector_const
  {
    public:

      explicit TubeVector(size_t n, TDomain& tdomain);
      size_t size() const;

      size_t nb_slices() const;
      const std::shared_ptr<SliceVector>& first_slice() const;
      std::shared_ptr<SliceVector>& first_slice();
      const std::shared_ptr<SliceVector>& last_slice() const;
      std::shared_ptr<SliceVector>& last_slice();

      //TubeVectorComponent operator[](size_t index);
      //const TubeVectorComponent operator[](size_t index) const;

      bool contains(const TrajectoryVector& value) const;

      TDomain& tdomain() const;
      Interval t0_tf() const;
      IntervalVector codomain() const;
      IntervalVector operator()(double t) const;
      IntervalVector operator()(const Interval& t) const;
      void set(const IntervalVector& codomain);

      friend std::ostream& operator<<(std::ostream& os, const TubeVector& x);


    protected:

      //friend class TubeVectorComponent;
      TDomain& _tdomain;


    public:

      using base_container = std::list<std::shared_ptr<TSlice>>;

      struct iterator : public base_container::iterator
      {
        using iterator_category = typename base_container::iterator::iterator_category;
        using difference_type   = typename base_container::iterator::difference_type;

        using value_type        = SliceVector;
        using pointer           = std::shared_ptr<SliceVector>;
        using reference         = SliceVector&;

        public:
          
          iterator(const TubeVector& tube_vector, base_container::iterator it) : 
            base_container::iterator(it), _tube_vector(tube_vector) { }

          reference operator*()
          {
            return *((*this)->get()->_slices.at(&_tube_vector));
          }

        protected:

          const TubeVector& _tube_vector;
      };

      iterator begin() { return iterator(*this, _tdomain._tslices.begin()); }
      iterator end()   { return iterator(*this, _tdomain._tslices.end()); }


      struct const_iterator : public base_container::const_iterator
      {
        using iterator_category = typename base_container::const_iterator::iterator_category;
        using difference_type   = typename base_container::const_iterator::difference_type;

        using value_type        = SliceVector;
        using pointer           = const std::shared_ptr<SliceVector>;
        using reference         = const SliceVector&;

        public:
          
          const_iterator(const TubeVector& tube_vector, base_container::const_iterator it) : 
            base_container::const_iterator(it), _tube_vector(tube_vector) { }

          reference operator*() const
          {
            return *((*this)->get()->_slices.at(&_tube_vector));
          }

        protected:

          const TubeVector& _tube_vector;
      };

      const_iterator begin() const { return const_iterator(*this, _tdomain._tslices.cbegin()); }
      const_iterator end() const   { return const_iterator(*this, _tdomain._tslices.cend()); }
  };
} // namespace codac

#endif