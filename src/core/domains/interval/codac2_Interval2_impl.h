/** 
 *  \file codac2_Interval2_impl.h
 * 
 *  This class reuses several functions developed for ibex::Interval2.
 *  See ibex::Interval2 (IBEX lib, main author: Gilles Chabert)
 *    https://ibex-lib.readthedocs.io
 *  
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Gilles Chabert, Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

// Inline functions

namespace codac2
{
  inline Interval2::Interval2()
    : gaol::interval(-oo,oo)
  { }

  inline Interval2::Interval2(double a)
    : gaol::interval(a)
  {
    if(a == -oo || a == oo)
      set_empty();
  }

  inline Interval2::Interval2(double a, double b)
    : gaol::interval(a,b)
  { }

  inline Interval2::Interval2(const Interval2& x)
    : gaol::interval(x)
  { }

  inline Interval2::Interval2(const std::array<double,1>& array)
    : Interval2(array[0])
  { }

  inline Interval2::Interval2(const std::array<double,2>& array)
    : Interval2(array[0],array[1])
  { }

  inline Interval2::Interval2(std::initializer_list<double> l)
    : Interval2()
  {
    init_from_list(l);
  }

  inline Interval2& Interval2::init(const Interval2& x)
  {
    *this = x;
    return *this;
  }

  inline Interval2& Interval2::init_from_list(const std::list<double>& l)
  {
    if(l.size() == 1)
      *this = Interval2(*l.begin());

    else if(l.size() == 2)
      *this = Interval2(*l.begin(),*std::prev(l.end()));

    else
    {
      assert_release("'Interval' can only be defined by one or two 'double' values.");
    }

    return *this;
  }

  inline Interval2& Interval2::operator=(double x)
  {
    if(x == -oo || x == oo)
      set_empty();
    else
      gaol::interval::operator=(x);

    return *this;
  }

  inline Interval2& Interval2::operator=(const Interval2& x)
  {
    gaol::interval::operator=(x);
    return *this;
  }

  inline bool Interval2::operator==(const Interval2& x) const
  {
    return (is_empty() && x.is_empty()) || (lb() == x.lb() && ub() == x.ub());
  }
  
  inline bool Interval2::operator!=(const Interval2& x) const
  {
    return !(*this == x);
  }

  inline double Interval2::lb() const
  {
    return gaol::interval::left();
  }

  inline double Interval2::ub() const
  {
    return gaol::interval::right();
  }

  inline double Interval2::mid() const
  {
    return gaol::interval::midpoint();
  }

  inline double Interval2::mag() const
  {
    return gaol::interval::mag();
  }

  inline double Interval2::mig() const
  {
    return gaol::interval::mig();
  }

  inline double Interval2::rand() const
  {
    if(is_empty())
      return std::numeric_limits<double>::quiet_NaN();

    double a = std::max<double>(next_float(-oo),lb());
    double b = std::min<double>(previous_float(oo),ub());
    double r = a + (((double)std::rand())/(double)RAND_MAX)*(b-a);
    // The above operation may result in a floating point outside the bounds,
    // due to floating-point errors. Such possible error is corrected below:
    return std::max<double>(lb(),std::min<double>(r,ub()));
  }

  inline double Interval2::rad() const
  {
    if(is_empty())
      return std::numeric_limits<double>::quiet_NaN();

    else if(is_unbounded())
      return oo;

    else
    {
      double t = mid();
      double t1 = (t-*this).ub();
      double t2 = (*this-t).ub();
      return (t1>t2) ? t1 : t2;
    }
  }

  inline double Interval2::diam() const
  {
    if(is_empty())
      return std::numeric_limits<double>::quiet_NaN();

    else
      return gaol::interval::width();
  }

  inline double Interval2::volume() const
  {
    return diam();
  }

  inline Index Interval2::size() const
  {
    return 1;
  }

  inline void Interval2::set_empty()
  {
    *this = Interval2::empty();
  }
  
  inline bool Interval2::is_empty() const
  {
    return gaol::interval::is_empty();
  }

  inline bool Interval2::contains(const double& x) const
  {
    return gaol::interval::set_contains(x);
  }

  inline bool Interval2::interior_contains(const double& x) const
  {
    return !is_empty() && x > lb() && x < ub();
  }

  inline bool Interval2::is_unbounded() const
  {
    return !gaol::interval::is_finite();
  }
  
  inline bool Interval2::is_degenerated() const
  {
    return is_empty() || gaol::interval::is_a_double();
  }

  inline bool Interval2::intersects(const Interval2 &x) const
  {
    return !is_empty() && !x.is_empty() && lb() <= x.ub() && ub() >= x.lb();
  }
  
  inline bool Interval2::is_disjoint(const Interval2& x) const
  {
    return is_empty() || x.is_empty() || lb() > x.ub() || ub() < x.lb();
  }
  
  inline bool Interval2::overlaps(const Interval2& x) const
  {
    return !is_empty() && !x.is_empty() && ub() > x.lb() && x.ub() > lb();
  }
  
  inline bool Interval2::is_subset(const Interval2& x) const
  {
    return is_empty() || (!x.is_empty() && x.lb() <= lb() && x.ub() >= ub());
  }
  
  inline bool Interval2::is_strict_subset(const Interval2& x) const
  {
    return !x.is_empty() && (is_empty() || (x.lb() < lb() && x.ub() >= ub()) || (x.ub() > ub() && x.lb() <= lb()));
  }
  
  inline bool Interval2::is_interior_subset(const Interval2& x) const
  {
    return is_empty() || (!x.is_empty() && (x.lb() == -oo || x.lb() < lb()) && (x.ub() == oo || x.ub() > ub()));
  }
  
  inline bool Interval2::is_strict_interior_subset(const Interval2& x) const
  {
    return !x.is_empty() && (is_empty() || (
         (x.lb() < lb() && (x.ub() == oo  || x.ub() > ub()))
      || (x.ub() > ub() && (x.lb() == -oo || x.lb() < lb()))
    ));
  }

  inline bool Interval2::is_superset(const Interval2& x) const
  {
    return x.is_subset(*this);
  }
  
  inline bool Interval2::is_strict_superset(const Interval2& x) const
  {
    return x.is_strict_subset(*this);
  }
  
  inline Interval2& Interval2::inflate(const double& rad)
  {
    (*this) += Interval2(-rad,rad);
    return *this;
  }

  inline bool Interval2::is_bisectable() const
  {
    if(is_empty())
      return false;
    double m = mid();
    return lb() < m && m < ub();
  }

  inline std::pair<Interval2,Interval2> Interval2::bisect(float ratio) const
  {
    assert_release(is_bisectable());
    assert_release(Interval2(0,1).interior_contains(ratio));

    if(lb() == -oo)
    {
      if(ub() == oo)
        return { Interval2(-oo,0), Interval2(0,oo) };
      else
        return { Interval2(-oo,-std::numeric_limits<double>::max()), Interval2(-std::numeric_limits<double>::max(),ub()) };
    }

    else if(ub() == oo)
      return { Interval2(lb(),std::numeric_limits<double>::max()), Interval2(std::numeric_limits<double>::max(),oo) };

    else
    {
      double m;

      if(ratio == 0.5)
        m = mid();

      else
      {
        m = lb() + ratio*diam();
        if(m >= ub())
          m = next_float(lb());

        assert(m < ub());
      }

      return { Interval2(lb(),m), Interval2(m,ub()) };
    }
  }

  inline std::vector<Interval2> Interval2::complementary(bool compactness) const
  {
    if(is_empty() || (compactness && is_degenerated()))
      return { {-oo,oo} };

    std::vector<Interval2> l;

    if(lb() > -oo)
      l.push_back({-oo,lb()});

    if(ub() < oo)
      l.push_back({ub(),oo});

    return l;
  }

  inline std::vector<Interval2> Interval2::diff(const Interval2& y, bool compactness) const
  {
    if(compactness && is_degenerated())
    {
      if(is_empty() || y.contains(lb()))
        return {};
      else
        return { *this };
    }

    std::vector<Interval2> l;
    for(const auto& li : y.complementary(compactness))
    {
      Interval2 inter = li & *this;
      if(!inter.is_degenerated())
        l.push_back(inter);
    }

    return l;
  }

  inline Interval2 operator&(const Interval2& x, const Interval2& y)
  {
    if(x.is_empty() || y.is_empty() || x.ub() < y.lb())
      return Interval2::empty();
    
    else
      return gaol::operator&(x,y);
  }

  inline Interval2 operator|(const Interval2& x, const Interval2& y)
  {
    return gaol::operator|(x,y);
  }

  inline const Interval2& operator+(const Interval2& x)
  {
    return x;
  }

  inline Interval2 operator+(const Interval2& x, double y)
  {
    if(y == -oo || y == oo)
      return Interval2::empty();

    else
      return gaol::operator+(x,y);
  }

  inline Interval2 operator+(double x, const Interval2& y)
  {
    if(x == -oo || x == oo)
      return Interval2::empty();

    else
      return gaol::operator+(x,y);
  }

  inline Interval2 operator+(const Interval2& x, const Interval2& y)
  {
    return gaol::operator+(x,y);
  }

  inline Interval2 operator-(const Interval2& x, double y)
  {
    if(y == -oo || y == oo)
      return Interval2::empty();

    else
      return gaol::operator-(x, y);
  }

  inline Interval2 operator-(double x, const Interval2& y)
  {
    if(x == -oo || x == oo)
      return Interval2::empty();

    else
      return gaol::operator-(x, y);
  }

  inline Interval2 operator-(const Interval2& x, const Interval2& y)
  {
    return gaol::operator-(x, y);
  }

  inline Interval2 operator*(const Interval2& x, double y)
  {
    if(y == -oo || y == oo)
      return Interval2::empty();

    else
      return gaol::operator*(x,y);
  }

  inline Interval2 operator*(double x, const Interval2& y)
  {
    if(x == -oo || x == oo)
      return Interval2::empty();

    else
      return gaol::operator*(x,y);
  }

  inline Interval2 operator*(const Interval2& x, const Interval2& y)
  {
    return gaol::operator*(x,y);
  }

  inline Interval2 operator/(const Interval2& x, double y)
  {
    if(y == -oo || y == oo)
      return Interval2::empty();

    else
      return gaol::operator/(x,y);
  }

  inline Interval2 operator/(double x, const Interval2& y)
  {
    if(x == -oo || x == oo)
      return Interval2::empty();

    else
      return gaol::operator/(x,y);
  }

  inline Interval2 operator/(const Interval2& x, const Interval2& y)
  {
    return gaol::operator/(x,y);
  }

  inline Interval2& Interval2::operator|=(const Interval2& x)
  {
    gaol::interval::operator|=(x);
    return *this;
  }

  inline Interval2& Interval2::operator&=(const Interval2& x)
  {
    gaol::interval::operator&=(x);
    return *this;
  }

  inline Interval2& Interval2::operator+=(double x)
  {
    if(x == -oo || x == oo)
      set_empty();
    else
      gaol::interval::operator+=(x);
    return *this;
  }

  inline Interval2& Interval2::operator+=(const Interval2& x)
  {
    gaol::interval::operator+=(x);
    return *this;
  }

  inline Interval2 Interval2::operator-() const
  {
    return 0.-*this;
  }

  inline Interval2& Interval2::operator-=(double x)
  {
    if(x == -oo || x == oo)
      set_empty();
    else
      gaol::interval::operator-=(x);
    return *this;
  }

  inline Interval2& Interval2::operator-=(const Interval2& x)
  {
    gaol::interval::operator-=(x);
    return *this;
  }

  inline Interval2& Interval2::operator*=(double x)
  {
    if(x == -oo || x == oo)
      set_empty();
    else
      gaol::interval::operator*=(x);
    return *this;
  }

  inline Interval2& Interval2::operator*=(const Interval2& x)
  {
    gaol::interval::operator*=(x);
    return *this;
  }

  inline Interval2& Interval2::operator/=(double x)
  {
    if(x == -oo || x == oo)
      set_empty();
    else
      gaol::interval::operator/=(x);
    return *this;
  }

  inline Interval2& Interval2::operator/=(const Interval2& x)
  {
    gaol::interval::operator/=(x);
    return *this;
  }

  inline Interval2 Interval2::empty()
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
  
  inline Interval2 Interval2::zero()
  {
    return gaol::interval::zero();
  }
  
  inline Interval2 Interval2::one()
  {
    return gaol::interval::one();
  }

  inline Interval2 Interval2::half_pi()
  {
    return gaol::interval::half_pi();
  }

  inline Interval2 Interval2::pi()
  {
    return gaol::interval::pi();
  }

  inline Interval2 Interval2::two_pi()
  {
    return gaol::interval::two_pi();
  }

  inline std::ostream& operator<<(std::ostream& os, const Interval2& x)
  {
    gaol::interval::precision(os.precision());
    gaol::operator<<(os,x);
    return os;
  }

  inline Interval2::Interval2(const gaol::interval& x)
    : gaol::interval(x)
  { }

  inline Interval2 operator""_i(long double x)
  {
    return Interval2(x);
  }
}