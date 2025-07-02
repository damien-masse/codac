/** 
 *  \file codac2_Modulo.cpp
 * 
 *  class representing [A] + k[B] with k \in Z. [B] is supposed to
 *  be small (e.g. punctual interval)
 *  
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Damien Massé
 *  \copyright  Copyright 2025 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */


#include <list>
#include <array>
#include <cmath>
#include "codac2_Interval.h"
#include "codac2_Domain.h"
#include "codac2_assert.h"
#include "codac2_IModulo.h"

using namespace codac2;

/* utility function : test of inclusion 
   x + k b \in [A]b with k \in Z and b in [B]
   (ie [A] - x/[B] contains k with k in Z
   2 cases : k exists, returns true and the value of k (smallest absolute value)
             k does not exist, returns false and k (non-integral) 
	     with x + kb \in [A]b (that is, x/[B] + k intersects [A])
   A must be non-empty and B>0
*/
static bool contains_k_util(const Interval &A, const Interval &B, double x, double &k) {
   Interval dK = A - x/B;
   Interval iK = integer(dK);
   if (iK.is_empty()) { k=dK.mid(); return false; }
   k = (iK.lb()>0.0 ? iK.lb() : iK.ub());
   return true;
}

bool IModulo::contains_k(double x, double &k) const {
   if (this->is_empty() || this->is_full()) { k = 0.0; return this->is_full(); }
   return contains_k_util(this->A,this->B,x,k);
}

static bool interior_contains_k_util(const Interval &A, const Interval &B, double x, double &k) {
  Interval dK = A-x/B;
  if (A.is_degenerated()) { k=dK.mid(); return false; }
  if (dK.lb()<0.0) {
     if (dK.ub()>0.0) { k=0; return true; }
     k = std::ceil(dK.ub())-1.0;
     return dK.interior_contains(k);
  } else {
     k = std::floor(dK.lb())+1.0;
     return dK.interior_contains(k);
  }
}

bool IModulo::interior_contains_k(double x, double &k) const {
  if (this->is_empty() || this->is_full()) { k = 0.0; return this->is_full(); }
  return interior_contains_k_util(this->A,this->B,x,k);
}

/* utility function :
   from two intervals A and B, compute
   A+kB such that (A+kB).mid() in [-B.lb()/2, B.ub()/2] 
   exceptions :
       A must not be empty, B must satisfy B.lb>0.0
       if A.diam() >= B.lb() (or (A+kB).diam() >= B.lb()), then
         A = [-oo,oo]
         (we may be more precise, but it is computationally costly
         and most probably useless)
*/      
static void center_interval(Interval &A, const Interval &B) {
   if (A.diam()>=1.0) {  A.init();  return; }
   double mid=A.mid();
   double kint;
   double mfrac=std::modf(mid,&kint);
   if (mfrac<=-0.5) kint+=1.0;
   else if (mfrac>0.5) kint-=1.0;
   A -= kint;
   if (A.diam()>=0.5) { A.init(); } /* I'm not sure it can happen */
}


/* recenter operation */
void IModulo::recenter() {
   if (this->is_empty()) return;
   if (this->is_full()) return;
   center_interval(this->A,this->B);
}

IModulo IModulo::join_close(const Interval &A1, const Interval &A2,
    const Interval &B) {
    double mA1 = A1.mid();
    double mA2 = A2.mid();
    if (mA1-mA2 > 0.5*B.ub()) { /* either A1-B, or A2+B ? */
       Interval E1 = A1 | (A2+B);
       if (E1.mid()>0.5*B.ub()) {
          return IModulo((A1-B) | A2, B);
       } 
       return IModulo(E1,B);
    } else if (mA1-mA2 < -0.5*B.ub()) {
       Interval E1 = (A1+B) | A2;
       if (E1.mid()>0.5*B.ub()) {
          return IModulo(A1 | (A2-B), B);
       } 
       return IModulo(E1,B);
    } else {
       return IModulo(A1|A2, B);
    }
}

IModulo::IModulo() : A(), B(Interval(0.0,oo)) {
}

IModulo::IModulo(double a, double b) : A(Interval(a)/b), B(b) {
  assert_release(b>0.0);
  this->recenter();
}

IModulo::IModulo(const Interval &a, double b) : A(a/b), B(b) {
   assert_release(b>0.0);
   this->recenter();
}

IModulo::IModulo(const Interval &a, const Interval &b) : A(a/b), B(b) {
   assert_release(b.lb()>0.0);
   this->recenter();
}

IModulo::IModulo(const IModulo& x) : A(x.A), B(x.B) {
}

IModulo& IModulo::init() {
   this->A.init();
   this->B = Interval(0.0,oo);
   return (*this);
}

IModulo& IModulo::init(const IModulo& x) {
   this->A = x.A;
   this->B = x.B;
   return (*this);
}

IModulo& IModulo::operator=(double x) {
   this->A = x;
   this->recenter();
   return (*this);
}

IModulo& IModulo::operator=(const Interval &x) {
   this->A = x;
   this->recenter();
   return (*this);
}

IModulo& IModulo::operator=(const IModulo &x) {
   this->A = x.A;
   this->B = x.B;
   return (*this);
}

bool IModulo::operator==(const IModulo& x) const {
   return (this->B == x.B && this->A == x.A);
}

bool IModulo::operator!=(const IModulo& x) const {
   return (this->B != x.B || this->A != x.A);
}

Interval IModulo::get_range() const {
   return this->A*this->B;
}

const Interval& IModulo::get_ratio() const {
   return this->A;
}

const Interval& IModulo::get_modulo() const {
   return this->B;
}

double IModulo::lb() const { return this->get_range().lb(); }
double IModulo::ub() const { return this->get_range().ub(); }
double IModulo::mid() const { return this->get_range().mid(); }
double IModulo::rad() const { return (this->A*this->B.ub()).rad(); }
double IModulo::diam() const { return (this->A*this->B.ub()).diam(); }
double IModulo::volume() const { return this.diam(); }
Index IModulo::size() const { return 1; }

void IModulo::set_empty() { this->A.set_empty(); }
bool IModulo::is_empty() const { return this->A.is_empty(); }
void IModulo::set_full() { this->A.init(); }
bool IModulo::is_full() const { return this->A.is_unbounded(); }

bool IModulo::contains(const double& x) const {
  double dummy;
  return this->contains_k(x,dummy);
}

bool IModulo::interior_contains(const double& x) const {
  double dummy;
  return this->interior_contains_k(x,dummy);
}

bool IModulo::is_degenerated() const {
  return this->A.is_degenerated();
}

bool IModulo::intersects(const IModulo& x) const {
  if (this->is_empty() || x.is_empty()) return false;
  if (this->is_full() || x.is_full()) return true;
  Interval B2 = this->B & x.B;
  if (B2.is_empty()) return false;
  double dummy;
  return contains_k_util(x.A-this->A,B2,0.0,dummy);
}

bool IModulo::is_disjoint(const IModulo& x) const {
  return !this->intersects(x);
}

bool IModulo::overlaps(const IModulo& x) const {
  if (this->is_empty() || x.is_empty()) return false;
  if (this->is_full()) return !x.is_degenerated();
  if (x.is_full()) return !this->is_degenerated();
  Interval B2 = this->B & x.B;
  if (B2.is_empty()) return false;
  if (x.is_degenerated() || this->is_degenerated()) return false;
  double dummy;
  return interior_contains_k_util(x.A-this->A,B2,0.0,dummy);
}

bool IModulo::is_subset(const IModulo& x) const {
  if (this->is_empty() || x.is_full()) return true;
  if (this->is_full() || x.is_empty()) return false;
  if (!this->B.is_subset(x.B)) return false;
  double k;
  if (!contains_k_util(x.A,this->B,this->A.lb(),k)) return false;
  double y = (x.A-k*this->B).ub();
  if (y>=this->A.ub()) return true;
  return ((x.A-(k-1)*this->B).lb()<=y);
      /* we know that (x.A-(k-1)*this->B).ub() > this->A.ub() because
              a) (x.A-k*this->B).ub() >= this->A.lb()
              b) this->A.lb + this->B.lb() > this->A.ub() (this not full)
              c) (x.A-(k-1)*this->B).ub() >= (x.A-k*this->B).ub()+this->B.lb()
      */
}

bool IModulo::is_strict_subset(const IModulo& x) const {
  return (this->is_subset(x) && !x.is_subset(*this));
}

bool IModulo::is_interior_subset(const IModulo& x) const {
  if (this->is_empty() || x.is_full()) return true;
  if (this->is_full() || x.is_empty()) return false;
  if (!this->B.is_subset(x.B)) return false;
  double k;
  if (x.is_degenerated()) return false;
  if (!interior_contains_k_util(x.A,this->B,this->A.lb(),k)) return false;
  double y = (x.A-k*this->B).ub();
  if (y>this->A.ub()) return true;
  return ((x.A-(k-1)*this->B).lb()<y);
}

bool IModulo::is_superset(const IModulo& x) const {
  return (x.is_subset(*this));
}

bool IModulo::is_strict_superset(const IModulo& x) const {
  return (x.is_strict_subset(*this));
}

IModulo& IModulo::inflate(const double& rad) {
  if (this->is_empty() || this->is_full()) return (*this);
  this->A.inflate(rad);
  if (this->A.diam()>=B.lb()) this->set_full();
  return (*this);
}

bool IModulo::is_bisectable() const {
  if (this->is_empty()) return false;
  if (this->is_full()) return (this->B.lb()>0);
  return (this->A.is_bisectable());
}

std::pair<IModulo,IModulo> IModulo::bisect(float ratio) const {
  assert_release(is_bisectable());
  assert_release(Interval(0,1).interior_contains(ratio));

  std::pair<Interval,Interval> apair;
  if (this->is_full()) {
     double m = -0.5*B.lb()+ratio*B.mid();
     apair.first = Interval(-0.5*B.lb(),m);
     apair.second = Interval(m,0.5*B.ub());
  } else {
     apair=this->A.bisect(ratio);
  }
  return { IModulo(apair.first,this->B), IModulo(apair.second,this->B) };
}

IModulo IModulo::complementary() const {
  if (this->is_empty()) return IModulo(Interval(),this->B);
  if (this->is_full()) return IModulo::empty(this->B);
  /* using the sign of A.mid */
  if (this->A.mid()>B.rad()) {
     Interval A2 ((this->A.ub()-B).lb(),A.lb());
     return IModulo(A2,this->B); /* new mid :
				   A.mid()-B.ub()/2 which is > -B.lb()/2 */
  } else {
     Interval A2 (A.ub(),(this->A.lb()+B).ub());
     return IModulo(A2,this->B); /* no centering needed in the exact case */
  }
}

std::vector<IModulo> IModulo::diff(const IModulo& y, bool compactness) const {
  if (this->is_empty() || y.is_full()) return { };
  if (y.is_empty()) return { (*this) };
  if (this->is_full()) return { y.complementary() }; /* we know y.B.lb>0 */
  Interval B2 = this->B & y.B;
  if (B2.is_empty()) return { };
  double k;
  if (contains_k_util(y.A,B2,this->A.lb(),k)) {
      double z = (y.A-k*B2).ub();
      if (z>=this->A.ub()) return { };
      double z2 = (y.A-(k-1)*B2).lb();
      if (z2<=z) return { };
      if (z2>this->A.ub()) z2=this->A.ub();
      return { IModulo(Interval(z,z2),B2) };
  } else {
      k = std::floor(k);
      Interval C = y.A-k*B2;
      std::vector<Interval> res = this->A.diff(C, compactness);
      std::vector<IModulo> ret;
      for (const auto &i : res) {
          ret.push_back(IModulo(i,B2));
      }
      return ret;
  }
}

std::vector<IModulo> IModulo::meet(const IModulo& y, bool compactness) const {
  if (this->is_empty() || y.is_empty()) return { };
  if (y.is_full()) return { (*this) };
  if (this->is_full()) return { y };
  Interval B2 = this->B & y.B;
  if (B2.is_empty()) return { };
  double k;
  if (contains_k_util(y.A,B2,this->A.lb(),k)) {
      Interval C1 = (y.A-k*B2) & this->A;
      Interval C2 = (y.A-(k-1)*B2) & this->A;
      if (C2.is_empty() || C1.intersects(C2)) {
         return { IModulo(C1 | C2, B2) };
      }
      if (compactness) {
         return { IModulo::join_close(C1,C2,B2) };
      } else {
         return { IModulo(C1,B2), IModulo(C2,B2) };
      }
  } else {
      k = std::floor(k);
      Interval C = (y.A-k*y.B) & this->A;
      if (C.is_empty()) return { };
      else {
	 return { IModulo(C,this->B) };
      }
  }
}

/* meet with compactness=true */
IModulo codac2::operator&(const IModulo &x, const IModulo &y) {
  std::vector<IModulo> z = x.meet(y,true);
  if (z.empty()) return IModulo::empty(x.B);
  return z[0];
}

/* meet with compactness=true */
IModulo& IModulo::operator&=(const IModulo& x) {
  std::vector<IModulo> z = this->meet(x,true);
  if (z.empty()) this->set_empty();
  else *this = z[0];
  return (*this);
}

std::vector<IModulo> IModulo::join(const IModulo& y, bool compactness) const {
  if (this->is_empty() && y.is_empty()) return { };
  if (this->is_full()) return { (*this) };
  if (y.is_full()) return { y };
  Interval B2 = this->B | y.B;
  if (this->A.diam()>=B2.lb() || y.A.diam()>=B2.lb()) 
	return { IModulo(Interval(),B2) };
  double k;
  Interval C;
  if (!contains_k_util(y.A,B2,this->A.lb(),k)) {
      k = std::floor(k);
      C = (y.A-k*B2);
      if (C.lb() > this->A.ub()) {
        if (compactness) {
           return { IModulo::join_close(this->A,y.A,B2) };
        } else {
           return { IModulo(this->A,B2), IModulo(y.A,B2) };
        }
      }
      C |= this->A;
  } else {
      C = (y.A-k*y.B) | this->A;
  }
  return { IModulo(C, this->B) };
}

/* join with compactness=true */
IModulo codac2::operator|(const IModulo& x, const IModulo& y) {
  std::vector<IModulo> z = x.join(y,true);
  if (z.empty()) return IModulo::empty(x.B);
  return z[0];
}

/* meet with compactness=true */
IModulo& IModulo::operator|=(const IModulo& x) {
  std::vector<IModulo> z = this->join(x,true);
  if (z.empty()) this->set_empty();
  else *this = z[0];
  return (*this);
}

/* Self addition of this and x */
IModulo& IModulo::operator+=(const Interval& x) {
  if (this->is_empty() || this->is_full()) return (*this);
  if (x.is_empty()) { this->set_empty(); return (*this); }
  this->A += x;
  this->recenter();
  return (*this);
}

/* Self addition of this and IModulo with same modulo */
IModulo& IModulo::operator+=(const IModulo& x) {
  if (this->is_empty() || this->is_full()) return (*this);
  if (x.is_empty() || x.is_full()) { *this=x; return (*this); }
  this->B &= x.B;
  if (this->B.is_empty()) { this->set_empty(); return (*this); }
  this->A += x.A;
  this->recenter();
  return (*this);
}

IModulo operator+(const IModulo &x, const IModulo &y) {
   IModulo r (x);
   return (r+=y);
}
IModulo operator+(const IModulo &x, const Interval &y) {
   IModulo r (x);
   return (r+=y);
}
IModulo operator+(const Interval &x, const IModulo &y) {
   IModulo r (y);
   return (r+=x);
}
/* negation */
IModulo IModulo::operator-() const {
  if (this->is_empty() || this->is_full()) return (*this);
  return IModulo(this->A,this->B);
}

/* Self substraction of this and x */
IModulo& IModulo::operator-=(const Interval& x) {
  if (this->is_empty() || this->is_full()) return (*this);
  if (x.is_empty()) { this->set_empty(); return (*this); }
  this->A -= x;
  this->recenter();
  return (*this);
}

/* Self addition of this and IModulo with same modulo */
IModulo& IModulo::operator-=(const IModulo& x) {
  if (this->is_empty() || this->is_full()) return (*this);
  if (x.is_empty() || x.is_full()) { *this=x; return (*this); }
  this->B &= x.B;
  if (this->B.is_empty()) { this->set_empty(); return (*this); }
  this->A -= x.A;
  this->recenter();
  return (*this);
}

IModulo operator-(const IModulo &x, const IModulo &y) {
   IModulo r (x);
   return (r-=y);
}
IModulo operator-(const IModulo &x, const Interval &y) {
   IModulo r (x);
   return (r-=y);
}
IModulo operator-(const Interval &x, const IModulo &y) {
   IModulo r (-y);
   return (r+=x);
}

IModulo& IModulo::operator*=(double x) {
  assert_release(x!=0.0);
  if (this->is_empty() || this->is_full()) return (*this);
  this->B *= std::abs(x);
  this->A *= x;
  return (*this);
}

IModulo& IModulo::operator*=(const Interval &x) {
  assert_release(!x.contains(0.0) && !x.is_unbounded());
  if (this->is_empty() || this->is_full()) return (*this);
  this->B *= abs(x);
  this->A *= x;
  this->recenter();
  return (*this);
}

IModulo operator*(const IModulo &x, double y) {
   IModulo r (x);
   return (r*=y);
}
IModulo operator*(const IModulo &x, const Interval &y) {
   IModulo r (x);
   return (r*=y);
}
IModulo operator*(double x, const IModulo &y) {
   IModulo r (y);
   return (r*=x);
}
IModulo operator*(const Interval &x, const IModulo &y) {
   IModulo r (y);
   return (r*=x);
}

IModulo& IModulo::operator/=(double x) {
  assert_release(x!=0.0);
  if (this->is_empty() || this->is_full()) return (*this);
  this->B /= std::abs(x);
  this->A /= x;
  return (*this);
}

IModulo& IModulo::operator/=(const Interval &x) {
  assert_release(!x.contains(0.0) && !x.is_unbounded());
  if (this->is_empty() || this->is_full()) return (*this);
  this->B /= abs(x);
  this->A /= x;
  this->recenter();
  return (*this);
}

IModulo operator/(const IModulo &x, double y) {
   IModulo r (x);
   return (r/=y);
}
IModulo operator/(const IModulo &x, const Interval &y) {
   IModulo r (x);
   return (r/=y);
}

IModulo IModulo::empty() {
  return IModulo(Interval::empty(), 0.0);
}

IModulo IModulo::empty(const Interval &b) {
  assert_release(b.lb()>0.0);
  return IModulo(Interval::empty(), b);
}

IModulo IModulo::zero(const Interval &b) {
  assert_release(b.lb()>0.0);
  return IModulo(Interval::zero(),b);
}

std::ostream& codac2::operator<<(std::ostream& os, const IModulo& x) {
   if (x.is_empty() || x.is_full()) {
     return (os << x.A);
   }
   return (os << x.A << "+k" << x.B);
}


