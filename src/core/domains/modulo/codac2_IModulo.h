/** 
 *  \file codac2_IModulo.h
 * 
 *  class representing ([A]+k)[B] with k in Z. [B] is supposed to
 *  be non-empty, small (e.g. punctual interval) and positive
 *  
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Damien Massé
 *  \copyright  Copyright 2025 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <list>
#include <array>
#include "codac2_Interval.h"
#include "codac2_Interval_operations.h"
#include "codac2_Domain.h"
#include "codac2_assert.h"


namespace codac2
{

  /**
   * \class IModulo
   * \brief IModulo is a representation of ([A]+k)[B] from two intervals
   * [A] and [B]. [B] is non-empty, positive, and with a small radius.
   * We consider mid([A]) in ] -0.5, 0.5 ] and rad([A])<0.5
   * In general the interface use [A][B] (the "real" values) and not the
   * internal representation [A] (the ratios).
   *
   *
   * Notes: "special" case full : [A]=[-oo,oo], [B]=[0,+oo]
   * "empty" : A=empty,  [B] = whatever
   */
  class IModulo 
  {
    public:

      /**
       * \brief Creates \f$[-\infty,\infty] + k \mathbb{R}+\f$
       */
      IModulo();

      /**
       * \brief Creates domain \f$(a/b+k)b\f$, \f$a\in\mathbb{R}\f$
       * and \f$b>0\f$
       * 
       * \param a single value contained in the resulting interval
       * \param b single (positive) value for the modulo
       */
      IModulo(double a, double b);

      /**
       * \brief Creates interval \f$A + kb\f$ 
       * 
       * \param A an interval
       * \param b single (positive) value for the modulo
       */
      IModulo(const Interval &a, double b);

      /**
       * \brief Creates interval \f$A + kB\f$ 
       * 
       * \param A an interval
       * \param B positive interval
       */
      IModulo(const Interval &a, const Interval &b);

      /**
       * \brief Creates an interval from another one
       * 
       * \param x the interval to be copied
       */
      IModulo(const IModulo& x);

      /**
       * \brief Creates an interval using ratio and modulo.
       *
       * \param R ratio
       * \param B modulo
       */
      static IModulo from_ratio(const Interval &R, const Interval &B);

      /**
       * \brief Sets the value of this interval to 
       * $f\[-\infty,\infty],]0,\infty]\f$
       * 
       * \note This function is used for template purposes.
       * 
       * \param x the value for re-initialization
       * \return a reference to this
       */
      IModulo& init();

      /**
       * \brief Sets the value of this interval to x
       * 
       * \note This function is used for template purposes.
       * 
       * \param x the value for re-initialization
       * \return a reference to this
       */
      IModulo& init(const IModulo& x);

      /**
       * \brief Sets this to x, keeping the modulo part
       * 
       * \param x real to be intervalized
       * \return a reference to this
       */
      IModulo& operator=(double x);

      /**
       * \brief Sets this to x, keeping the modulo part
       * 
       * \param x interval to be copied
       * \return a reference to this
       */
      IModulo& operator=(const Interval& x);

      /**
       * \brief Sets this to x
       * 
       * \param x IModulo to be copied
       * \return a reference to this
       */
      IModulo& operator=(const IModulo& x);

      /**
       * \brief Comparison (equality) between two imodulos
       * 
       * \param x imodulo to be compared with
       * \return true iff this and x are exactly the same intervals
       */
      bool operator==(const IModulo& x) const;

      /**
       * \brief Comparison (non equality) between two imodulos 
       * 
       * \param x interval to be compared with
       * \return true iff this and x are not exactly the same intervals
       */
      bool operator!=(const IModulo& x) const;


      /**
       * \brief Returns [A][B] (from \f$(A+k)B\f$)
       * 
       * \return the interval [A][B]
       */
      Interval get_range() const;

      /**
       * \brief Returns [A] (from \f$(A+k)B\f$)
       * 
       * \return the interval [A]
       */
      const Interval& get_ratio() const;

      /**
       * \brief Returns B (from \f$A+kB\f$)
       * 
       * \return the interval B
       */
      const Interval& get_modulo() const;

      /**
       * \brief Returns the lower bound of [A][B]
       * 
       * \return lower bound
       */
      double lb() const;

      /**
       * \brief Returns the upper bound of [A][B]
       * 
       * \return upper bound
       */
      double ub() const;
      
      /**
       * \brief Returns the midpoint of [A][B]
       * 
       * \return midpoint
       */
      double mid() const;

      /**
       * \brief Returns the radius of [A]B.ub()
       * 
       * \return the radius, or 0 if this is empty
       */
      double rad() const;

      /**
       * \brief Returns the diameter of [A]B.ub()
       * 
       * \return the diameter, or 0 if this is empty
       */
      double diam() const;

      /**
       * \brief Returns the diameter of A
       * 
       * \note This function is used for template purposes. See diam()
       * 
       * \return the diameter, or 0 if this is empty
       */
      double volume() const;

      /**
       * \brief Returns the dimension of A (which is always 1)
       * 
       * \note This function is used for template purposes.
       * 
       * \return 1
       */
      Index size() const;
    
      /**
       * \brief Sets this to the empty set
       */
      void set_empty();

      /**
       * \brief Tests if this is empty
       * 
       * \return true in case of empty set
       */
      bool is_empty() const;

      /**
       * \brief Sets this to \f$\mathbb{R}\f$
       */
      void set_full();

      /**
       * \brief Tests if this is \f$\mathbb{R}\f$
       * 
       * \return true in case of empty set
       */
      bool is_full() const;

      /**
       * \brief Tests if this (may) contains x
       * 
       * \note if x is an "open bound", i.e., infinity, returns
       * true iff this is \f$\mathbb{R}\f$
       * 
       * \param x real value
       * \return true if this contains x
       */
      bool contains(const double& x) const;

      /**
       * \brief Tests if the interior of this contains x
       * 
       * \param x real value
       * \return true if the interior of this contains x
       */
      bool interior_contains(const double& x) const;

      /**
       * \brief Tests if the representation interval of this is degenerated,
       * that is, in the form of \f$[a,a]\f$
       *
       * \note An empty interval is considered here as degenerated
       * 
       * \return true if this is degenerated
       */
      bool is_degenerated() const;

      /**
       * \brief Tests if this and x "intersect"
       * that is : \f$\exists k \in Z, this.A - x.A \cap k (this.B\cap A.B)\f$
       * not empty this property is symmetric.
       * 
       * \param x the other interval
       * \return true if the intervals intersect
       */
      bool intersects(const IModulo &x) const;

      /**
       * \brief Tests if this and x do not "intersect"
       * 
       * \param x the other interval
       * \return true if the intervals are disjoint
       */
      bool is_disjoint(const IModulo& x) const;

      /**
       * \brief Tests if this and x "intersect" and their intersection
       * has a non-null volume
       * this.A and x.A must be non-degenerated, and 
       * \f$\exists k \in Z, this.A - x.A \cap k this.B\f$ has a non-null volume
       * this property is symmetric.
       * 
       * \param x the other interval
       * \return true
       */
      bool overlaps(const IModulo &x) const;

      /**
       * \brief Tests if this is a "subset" of x, that is, x "contains"
       * all elements of this
       *
       * \param x the other interval (same modulo)
       * \return true iff this interval is a subset of x
       */
      bool is_subset(const IModulo& x) const;

      /**
       * \brief Tests if this is a "subset" of x, but x is not a "subset"
       * of this
       *
       * \note In particular, \f$[-\infty,\infty]\f$ is not a strict
       * subset of \f$[-\infty,\infty]\f$ and \f$\varnothing\f$ is not 
       * a strict subset of \f$\varnothing\f$, although in both cases,
       * the first is inside the interior of the second.
       * 
       * \param x the other interval
       * \return true iff this interval is a strict subset of x
       */
      bool is_strict_subset(const IModulo& x) const;

      /**
       * \brief Tests if this is in the interior of x, that is,
       * the interior of x contains all elements of x
       * 
       * \note In particular, \f$[-\infty,\infty]\f$ is in the interior
       * of \f$[-\infty,\infty]\f$ and \f$\varnothing\f$ is in the 
       * interior of \f$\varnothing\f$.
       * 
       * \param x the other interval
       * \return true iff this interval is in the interior of x
       */
      bool is_interior_subset(const IModulo& x) const;

      /**
       * \brief Tests if this is a superset of x
       *
       * \note If this is empty, always returns true
       * 
       * \param x the other interval
       * \return true iff this interval is a superset of x
       */
      bool is_superset(const IModulo& x) const;

      /**
       * \brief Tests if this is a superset of x and different from x
       *
       * \see is_strict_subset()
       * 
       * \param x the other interval
       * \return true iff this interval is a strict superset of x
       */
      bool is_strict_superset(const IModulo& x) const;

      /**
       * \brief Adds [-rad,+rad] to this
       * 
       * \param rad the radius of the inflation
       * \return a reference to this
       */
      IModulo& inflate(const double& rad);

      /**
       * \brief Tests if this can be bisected into two non-degenerated intervals
       * 
       * \note Examples of non bisectable intervals are [0,next_float(0)] or [DBL_MAX,+oo).
       * 
       * \return true iff this can be bisected
       */
      bool is_bisectable() const;

      /**
       * \brief Bisects this into two subintervals
       * 
       * \pre is_bisectable() must be true
       * \pre 0 < ratio < 1
       * \param ratio optional proportion of split (0.5 corresponds to middle)
       * \return a pair of resulting subintervals
       */
      std::pair<IModulo,IModulo> bisect(float ratio = 0.49) const;

      /**
       * \brief Computes the ``complementary'' of this. Unlike intervals
       * the result is always an IModulo.
       *
       * The complementary satisfies the following properties:
       * 1) \f$\forall x in \reals, x \in_B this or x \in_B comp(this)\f$
       * 2) this and comp(this) do not overlap
       * Two intervals are possible for this.A, we select the one with
       * mid closest to 0.
       * 
       * \return the complementary
       */
      IModulo complementary() const;

      /**
       * \brief Computes the result of \f$[this]\[y]\f$
       * can be seen as the intersection with the complementary
       * (hence intersection, and not difference of the modulo) 
       * 
       * \param y IModulo to ``remove'' from this
       * \param compactness optional boolean to obtain or not disjoint intervals
       * \return a vector of resulting diff intervals
       */
      std::vector<IModulo> diff(const IModulo& y, bool compactness = true) const;

      /**
       * \brief Computes the result of \f$[this] \cap_B [y]\f$
       * basically, the elements of [this] which are related to [y]
       * 
       * \param y IModulo to ``intersect'' to this
       * \param compactness optional boolean to obtain or not disjoint intervals
       * \return a vector of resulting diff intervals
       */
      std::vector<IModulo> meet(const IModulo& y, bool compactness = true) const;

      /**
       * \brief Computes the result of \f$[this] \cup_B [y]\f$
       * 
       * \param y IModulo to ``join'' to this
       * \param compactness optional boolean to obtain or not disjoint intervals
       * \return a vector of resulting diff intervals
       */
      std::vector<IModulo> join(const IModulo& y, bool compactness = true) const;

      /**
       * \brief Self union of this and x
       * 
       * \param x the other interval
       * \return a reference to this
       */
      IModulo& operator|=(const IModulo& x);

      /**
       * \brief Self intersection of this and x
       * 
       * \param x the other interval
       * \return a reference to this
       */
      IModulo& operator&=(const IModulo& x);

      /**
       * \brief Self addition of this and an interval x
       * 
       * \param x the real to add
       * \return a reference to this
       */
      IModulo& operator+=(const Interval& x);

      /**
       * \brief Self addition of this and x
       * 
       * \param x the other interval (same modulo)
       * \return a reference to this
       */
      IModulo& operator+=(const IModulo& x);

      /**
       * \brief negation of this
       * 
       * \return (- this)
       */
      IModulo operator-() const;

      /**
       * \brief Self substraction of this and a real x
       * 
       * \param x the real to substract
       * \return a reference to this
       */
      IModulo& operator-=(const Interval& x);

      /**
       * \brief Self substraction of this and x
       * 
       * \param x the other interval
       * \return a reference to this
       */
      IModulo& operator-=(const IModulo& x);

      /**
       * \brief Self multiplication of this and a real x not null
       * the modulo is also changed
       * 
       * \param x the real to multiply
       * \return a reference to this
       */
      IModulo& operator*=(double x);

      /**
       * \brief Self multiplication of this and x (constant sign, bounded)
       * the modulo is also changed
       * 
       * \param x the other interval
       * \return a reference to this
       */
      IModulo& operator*=(const Interval& x);

      /**
       * \brief Self multiplication of this and a real x not null
       * the modulo is also changed
       * 
       * \param x the real to multiply
       * \return a reference to this
       */
      IModulo& operator/=(double x);

      /**
       * \brief Self multiplication of this and x (constant sign, bounded)
       * the modulo is also changed
       * 
       * \param x the other interval
       * \return a reference to this
       */
      IModulo& operator/=(const Interval& x);

      /**
       * \brief Provides an empty interval
       * 
       * \return an empty set
       */
      static IModulo empty(const Interval &b);

      /**
       * \brief Provides an empty interval
       * 
       * \return an empty set
       */
      static IModulo empty();

      /**
       * \brief Provides an interval for \f$[0]\f$
       * 
       * \return an interval containing \f$0\f$
       */
      static IModulo zero(const Interval &b);

      friend IModulo operator|(const IModulo &x, const IModulo &y);
      friend IModulo operator&(const IModulo &x, const IModulo &y);
      friend IModulo operator+(const IModulo &x, const IModulo &y);
      friend IModulo operator+(const Interval &x, const IModulo &y);
      friend IModulo operator+(const IModulo &x, const Interval &y);
      friend IModulo operator-(const IModulo &x, const IModulo &y);
      friend IModulo operator-(const Interval &x, const IModulo &y);
      friend IModulo operator-(const IModulo &x, const Interval &y);
      friend IModulo operator*(const IModulo &x, double y);
      friend IModulo operator*(const IModulo &x, const Interval& y);
      friend IModulo operator*(double x, const IModulo &y);
      friend IModulo operator*(const Interval& x, const IModulo &y);
      friend IModulo operator/(const IModulo &x, double y);
      friend IModulo operator/(const IModulo &x, const Interval &y);

      
      friend std::ostream& operator<<(std::ostream& os, const IModulo& x);

#if 0
      /***** angle specific functions (TODO) *****/

      /**
       * \brief Provides an interval for \f$[\frac{\pi}{2}]\f$
       * 
       * \return an interval containing \f$\frac{\pi}{2}\f$
       */
      static Interval half_pi();

      /**
       * \brief Provides an interval for \f$[\pi]\f$
       * 
       * \return an interval containing \f$\pi\f$
       */
      static Interval pi();
#endif

      protected:
         /**
          * \brief Contains testing, with the value of k 
          * given (i.e. k is a value, integral if possible,
          * s.t. x + k[B] \in [A]).
          *
          * \param x the value
          * \param k the (returned) integral value if possible, or a fp if not
          * \return true if k is integral,false if not */
         bool contains_k(double x, double &k) const;

         /**
          * \brief Interior contains testing, with the value of k 
          * given (i.e. k is a value, integral if possible,
          * s.t. x + k[B] \in [A]).
          *
          * \param x the value
          * \param k the (returned) integral value if possible, or a fp if not
          * \return true if k is integral,false if not */
         bool interior_contains_k(double x, double &k) const;


#if 0
      _dec_friend_interval2_unary_op(sqr)
      _dec_friend_interval2_unary_op(sqrt)
      _dec_friend_interval2_unary_op(exp)
      _dec_friend_interval2_unary_op(log)
      _dec_friend_interval2_unary_op(cos)
      _dec_friend_interval2_unary_op(sin)
      _dec_friend_interval2_unary_op(tan)
      _dec_friend_interval2_unary_op(acos)
      _dec_friend_interval2_unary_op(asin)
      _dec_friend_interval2_unary_op(atan)
      _dec_friend_interval2_unary_op(cosh)
      _dec_friend_interval2_unary_op(sinh)
      _dec_friend_interval2_unary_op(tanh)
      _dec_friend_interval2_unary_op(acosh)
      _dec_friend_interval2_unary_op(asinh)
      _dec_friend_interval2_unary_op(atanh)
      _dec_friend_interval2_unary_op(abs)
      _dec_friend_interval2_unary_op(sign)
      _dec_friend_interval2_unary_op(integer)
      _dec_friend_interval2_unary_op(floor)
      _dec_friend_interval2_unary_op(ceil)

      #define _dec_friend_interval2_binary_op(f) \
        friend Interval f(const Interval&, const Interval&); \

      _dec_friend_interval2_binary_op(max)
      _dec_friend_interval2_binary_op(min)
      _dec_friend_interval2_binary_op(atan2)

      friend Interval pow(const Interval&, int);
      friend Interval pow(const Interval&, double);

      _dec_friend_interval2_binary_op(pow)

      friend Interval root(const Interval&, int);
      friend Interval chi(const Interval&, const Interval&, const Interval&);

      friend struct AbsOp;
      friend struct AcosOp;
      friend struct AddOp;
      friend struct ChiOp;
      friend struct DivOp;
      friend struct MulOp;
      friend struct SubOp;
      friend struct AsinOp;
      friend struct AtanOp;
      friend struct Atan2Op;
      friend struct CosOp;
      friend struct CoshOp;
      friend struct DetOp;
      friend struct ExpOp;
      friend struct LogOp;
      friend struct PowOp;
      friend struct SinOp;
      friend struct SinhOp;
      friend struct SqrOp;
      friend struct SqrtOp;
      friend struct TanOp;
      friend struct TanhOp;
#endif
     private:
        Interval A,B;
  
        void recenter();
        static IModulo join_close(const Interval &A1, const Interval &A2,
		const Interval &B);
  };
      
       IModulo operator|(const IModulo &x, const IModulo &y);
       IModulo operator&(const IModulo &x, const IModulo &y);
       IModulo operator+(const IModulo &x, const IModulo &y);
       IModulo operator+(const Interval &x, const IModulo &y);
       IModulo operator+(const IModulo &x, const Interval &y);
       IModulo operator-(const IModulo &x, const IModulo &y);
       IModulo operator-(const Interval &x, const IModulo &y);
       IModulo operator-(const IModulo &x, const Interval &y);
       IModulo operator*(const IModulo &x, double y);
       IModulo operator*(const IModulo &x, const Interval& y);
       IModulo operator*(double x, const IModulo &y);
       IModulo operator*(const Interval& x, const IModulo &y);
       IModulo operator/(const IModulo &x, double y);
       IModulo operator/(const IModulo &x, const Interval &y);

      
       std::ostream& operator<<(std::ostream& os, const IModulo& x);
}
