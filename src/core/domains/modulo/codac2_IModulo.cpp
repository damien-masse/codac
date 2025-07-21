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
   (ie [A] - x/[B] contains k with k in Z )
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
   from A, compute A+k such that (A+k).mid() in ]-0.5,0.5]
   exceptions :
       A must not be empty
       if A.diam() >= 1.0 
         A = [-oo,oo]
*/      
static void center_interval(Interval &A) {
   if (A.diam()>=1.0) {  A.init();  return; }
   double mid=A.mid();
   if (mid<=0.5 && mid>-0.5) return;
   double kint;
   double mfrac=std::modf(mid,&kint);
   if (mfrac<=-0.5) kint-=1.0;
   else if (mfrac>0.5) kint+=1.0;
   A -= kint;
   if (A.diam()>=1.0) { A.init(); } /* I'm not sure it can happen */
}


/* recenter operation */
void IModulo::recenter() {
   if (this->is_empty()) return;
   if (this->is_full()) return;
   center_interval(this->A);
}

IModulo IModulo::join_close(const Interval &A1, const Interval &A2,
    const Interval &B) {
    double mA1 = A1.mid();
    double mA2 = A2.mid();
    if (mA1-mA2 > 0.5) { /* either A1-1, or A2+1 */
       if (A1.lb()<=-A2.ub()) /* A1 | A2+1 */
          return IModulo::from_ratio(A1 | (A2+1.0),B);
       return IModulo::from_ratio((A1-1.0) | A2,B);
    } else if (mA1-mA2 < -0.5) {
       if (A2.lb()<=-A1.ub()) /* A2 | A1+1 */
          return IModulo::from_ratio(A2 | (A1+1.0),B);
       return IModulo::from_ratio((A2-1.0) | A1,B);
    } else {
       return IModulo::from_ratio(A1|A2, B);
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

IModulo IModulo::from_ratio(const Interval &A, const Interval &B) {
    IModulo build;
    build.A=A;
    build.B=B;
    build.recenter();
    return build;
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
double IModulo::volume() const { return this->diam(); }
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
  if (this->B.is_disjoint(x.B)) return false;
  if (this->A.lb()>x.A.ub()) {
     return (this->A.intersects(x.A+1.0));
  } else if (this->A.ub()<x.A.lb()) {
     return (this->A.intersects(x.A-1.0));
  } else return true;
}

bool IModulo::is_disjoint(const IModulo& x) const {
  return !this->intersects(x);
}

bool IModulo::overlaps(const IModulo& x) const {
  if (this->is_empty() || x.is_empty()) return false;
  if (this->is_full()) return !x.is_degenerated();
  if (x.is_full()) return !this->is_degenerated();
  if (this->B.is_disjoint(x.B)) return false;
  if (this->A.lb()>=x.A.ub()) {
     return (this->A.overlaps(x.A+1.0));
  } else if (this->A.ub()<=x.A.lb()) {
     return (this->A.overlaps(x.A-1.0));
  } else return true;
}

bool IModulo::is_subset(const IModulo& x) const {
  if (this->is_empty() || x.is_full()) return true;
  if (this->is_full() || x.is_empty()) return false;
  if (!this->B.is_subset(x.B)) return false;
  if (this->A.is_subset(x.A)) return true;
  if (this->A.is_subset(x.A+1.0)) return true;
  return (this->A.is_subset(x.A-1.0));
}

bool IModulo::is_strict_subset(const IModulo& x) const {
  if (this->is_empty() || x.is_full()) return true;
  if (this->is_full() || x.is_empty()) return false;
  if (!this->B.is_subset(x.B)) return false;
  if (this->A.is_strict_subset(x.A)) return true;
  if (this->A.is_subset(x.A+1.0)) return true;
  return (this->A.is_subset(x.A-1.0));
}

bool IModulo::is_interior_subset(const IModulo& x) const {
  if (this->is_empty() || x.is_full()) return true;
  if (this->is_full() || x.is_empty()) return false;
  if (!this->B.is_subset(x.B)) return false;
  if (this->A.is_interior_subset(x.A)) return true;
  if (this->A.is_interior_subset(x.A+1.0)) return true;
  return (this->A.is_interior_subset(x.A-1.0));
}

bool IModulo::is_superset(const IModulo& x) const {
  return (x.is_subset(*this));
}

bool IModulo::is_strict_superset(const IModulo& x) const {
  return (x.is_strict_subset(*this));
}

IModulo& IModulo::inflate(const double& rad) {
  if (this->is_empty() || this->is_full()) return (*this);
  this->A.inflate((rad/this->B).ub());
  if (this->A.diam()>=1.0) this->set_full();
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
     apair=Interval(-0.5,0.5).bisect(ratio);
  } else {
     apair=this->A.bisect(ratio);
  }
  return { IModulo::from_ratio(apair.first,this->B),
	   IModulo::from_ratio(apair.second,this->B) };
}

IModulo IModulo::complementary() const {
  if (this->is_empty()) return IModulo(Interval(),this->B);
  if (this->is_full()) return IModulo::empty(this->B);
  /* using the sign of A.mid */
  if (this->A.mid()>0.0) {
     Interval A2 ((this->A-1.0).ub(),A.lb());
     return IModulo::from_ratio(A2,this->B); /* new mid : A.mid()-0.5 */
  } else {
     Interval A2 (A.ub(),(this->A+1.0).lb());
     return IModulo::from_ratio(A2,this->B); /* new mid : A.mid()+0.5 */
  }
}

std::vector<IModulo> IModulo::diff(const IModulo& y, bool compactness) const {
  if (this->is_empty() || y.is_full()) return { };
  if (y.is_empty()) return { (*this) };
  if (this->is_full()) return { y.complementary() }; /* we know y.B.lb>0 */
  Interval B2 = this->B & y.B;
  if (B2.is_empty()) return { };
  if (y.A.lb()<=this->A.lb()) {
     if (y.A.ub()<=this->A.lb()) {
       Interval yA = y.A+1.0; /* Warning: use "inner" sum */
       if (yA.lb()<this->A.ub()) {
          if (yA.lb()<=this->A.lb()) /* A subset yA (Va) */
             return { };
          else 
          if (yA.ub()<this->A.ub()) { /* A superset yA (VIa) */
             if (compactness) return { IModulo::from_ratio(this->A,B2) };
             return { IModulo::from_ratio(Interval(this->A.lb(),yA.lb()),B2),
                      IModulo::from_ratio(Interval(yA.ub(),this->A.ub()),B2) };
          }
          else /* yA intersects A ( IVa ) */
          return { IModulo::from_ratio(Interval(this->A.lb(),yA.lb()),B2) };
       } else /* disjoint (Ia) */
          return { IModulo::from_ratio(this->A,B2) };
     } 
     else if (y.A.ub()<this->A.ub()) {
       Interval yA = y.A+1.0; /* Warning: use "inner" sum */
       if (yA.lb()<=this->A.ub()) /* double intersection (IIIa) */
           return { IModulo::from_ratio(Interval(y.A.ub(),yA.lb()),B2) };
       else /* y.A intersects A (IIa) */
           return { IModulo::from_ratio(Interval(y.A.ub(),this->A.ub()),B2) };
     } else  /* A subset y.A (Vc) */
        return { };
   } else {
     if (y.A.lb()<this->A.ub()) {
        if (y.A.ub()<this->A.ub()) { /* A superset y.A (VIc) */
           if (compactness) return { IModulo::from_ratio(this->A,B2) };
           return { IModulo::from_ratio(Interval(this->A.lb(),y.A.lb()),B2),
                    IModulo::from_ratio(Interval(y.A.ub(),this->A.ub()),B2) };
        } else {
           Interval yA = y.A-1.0; /* Warning: use "inner" sum */
           if (yA.ub()<=this->A.lb()) /* y.A intersects A (IVb) */
              return { IModulo::from_ratio(Interval(this->A.lb(),y.A.lb()),B2) };
           else { /* double intersection (IIIb) */
              return { IModulo::from_ratio(Interval(yA.ub(),y.A.lb()),B2) };
           }
        }
     } else {
        Interval yA = y.A-1.0; /* Warning: use "inner" sum */
        if (yA.ub()<=this->A.lb()) /* disjoint (Ib) */
           return { IModulo::from_ratio(this->A,B2) };
        else if (yA.ub()<this->A.ub()) {
           if (yA.lb()<=this->A.lb()) { /* yA intersects A (IIb) */
               return { IModulo::from_ratio(Interval(yA.ub(),this->A.ub()),B2) };
           } else { /* A superset yA (VIb) */
              if (compactness) return { IModulo::from_ratio(this->A,B2) };
              return { IModulo::from_ratio(Interval(this->A.lb(),y.A.lb()),B2),
                       IModulo::from_ratio(Interval(y.A.ub(),this->A.ub()),B2) };
           }
        } else /* A subset yA (Vb) */
           return { };
     }   
   }
}

std::vector<IModulo> IModulo::meet(const IModulo& y, bool compactness) const {
  if (this->is_empty() || y.is_empty()) return { };
  if (y.is_full()) return { (*this) };
  if (this->is_full()) return { y };
  Interval B2 = this->B & y.B;
  if (B2.is_empty()) return { };

  if (y.A.lb()<=this->A.lb()) {
     if (y.A.ub()<this->A.lb()) {
       Interval yA = y.A+1.0;
       if (yA.lb()<=this->A.ub()) {
          if (yA.lb()<=this->A.lb()) /* A subset yA (Va) */
             return { IModulo::from_ratio(this->A,B2) };
          else 
          if (yA.ub()<this->A.ub()) { /* A superset yA (VIa) */
             return { IModulo::from_ratio(y.A,B2) };
          }
          else /* yA intersects A ( IVa ) */
          return { IModulo::from_ratio(Interval(yA.lb(),this->A.ub()),B2) };
       } else /* disjoint (Ia) */
          return { };
     } 
     else if (y.A.ub()<=this->A.ub()) {
       Interval yA = y.A+1.0; 
       if (yA.lb()<=this->A.ub()) { /* double intersection (IIIa) */
           if (compactness) return { IModulo::from_ratio(this->A,B2) };
           return { IModulo::from_ratio(Interval(this->A.lb(),y.A.ub()),B2),
                    IModulo::from_ratio(Interval(yA.lb(),this->A.ub()),B2) };
       }
       else /* y.A intersects A (IIa) */
           return { IModulo::from_ratio(Interval(this->A.lb(),y.A.ub()),B2) };
     } else  /* A subset y.A (Vc) */
        return { IModulo::from_ratio(this->A,B2) };
   } else {
     if (y.A.lb()<this->A.ub()) {
        if (y.A.ub()<this->A.ub()) { /* A superset y.A (VIc) */
           return { IModulo::from_ratio(y.A,B2) };
        } else {
           Interval yA = y.A-1.0; 
           if (yA.ub()<this->A.lb()) /* y.A intersects A (IVb) */
              return { IModulo::from_ratio(Interval(y.A.lb(),this->A.ub()),B2) };
           else { /* double intersection (IIIb) */
              if (compactness) return { IModulo::from_ratio(this->A,B2) };
              return { IModulo::from_ratio(Interval(this->A.lb(),yA.ub()),B2),
                       IModulo::from_ratio(Interval(y.A.lb(),this->A.ub()),B2) };
           }
        }
     } else {
        Interval yA = y.A-1.0;
        if (yA.ub()<this->A.lb()) /* disjoint (Ib) */
           return { };
        else if (yA.ub()<this->A.ub()) {
           if (yA.lb()<this->A.lb()) { /* yA intersects A (IIb) */
               return { IModulo::from_ratio(Interval(this->A.lb(),yA.ub()),B2) };
           } else { /* A superset yA (VIb) */
               return { IModulo::from_ratio(y.A,B2) };
           }
        } else /* A subset yA (Vb) */
           return { IModulo::from_ratio(this->A,B2) };
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

  if (y.A.lb()<=this->A.lb()) {
     if (y.A.ub()<this->A.lb()) {
       Interval yA = y.A+1.0;
       if (yA.lb()<=this->A.ub()) {
          if (yA.lb()<=this->A.lb()) /* A subset yA (Va) */
             return { IModulo::from_ratio((this->A-1.0) | y.A,B2) };
          else 
          if (yA.ub()<this->A.ub()) { /* A superset yA (VIa) */
             return { IModulo::from_ratio(this->A,B2) };
          }
          else /* yA intersects A ( IVa ) */
          return { IModulo::from_ratio(Interval(this->A.lb(),yA.ub()),B2) };
       } else {/* disjoint (Ia) */
          if (compactness) return { IModulo::join_close(this->A,y.A,B2) };
          return { IModulo::from_ratio(this->A,B2),
                   IModulo::from_ratio(y.A,B2) };
       }
     } 
     else if (y.A.ub()<=this->A.ub()) {
       Interval yA = y.A+1.0; 
       if (yA.lb()<=this->A.ub()) /* double intersection (IIIa) */
           return { IModulo(Interval(),B2) };
       else /* y.A intersects A (IIa) */
           return { IModulo::from_ratio(Interval(y.A.lb(),this->A.ub()),B2) };
     } else  /* A subset y.A (Vc) */
        return { IModulo::from_ratio(y.A,B2) };
   } else {
     if (y.A.lb()<this->A.ub()) {
        if (y.A.ub()<this->A.ub()) { /* A superset y.A (VIc) */
           return { IModulo::from_ratio(this->A,B2) };
        } else {
           Interval yA = y.A-1.0; 
           if (yA.ub()<this->A.lb()) /* y.A intersects A (IVb) */
              return { IModulo::from_ratio(Interval(this->A.lb(),y.A.ub()),B2) };
           else  /* double intersection (IIIb) */
              return { IModulo(Interval(),B2) };
        }
     } else {
        Interval yA = y.A-1.0;
        if (yA.ub()<this->A.lb()) { /* disjoint (Ib) */
          if (compactness) return { IModulo::join_close(this->A,y.A,B2) };
          return { IModulo::from_ratio(this->A,B2),
                   IModulo::from_ratio(y.A,B2) };
        } else if (yA.ub()<this->A.ub()) {
           if (yA.lb()<this->A.lb()) { /* yA intersects A (IIb) */
               return { IModulo::from_ratio(Interval(yA.lb(),this->A.ub()),B2) };
           } else { /* A superset yA (VIb) */
               return { IModulo::from_ratio(this->A,B2) };
           }
        } else /* A subset yA (Vb) */
           return { IModulo::from_ratio((this->A+1.0) & y.A,B2) };
     }   
  }
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
  this->A += x/this->B;
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
  if (x<0.0) { this->A = -this->A; this->recenter(); }
  return (*this);
}

IModulo& IModulo::operator*=(const Interval &x) {
  assert_release(!x.contains(0.0) && !x.is_unbounded());
  if (this->is_empty() || this->is_full()) return (*this);
  this->B *= abs(x);
  if (x.ub()<0.0) { this->A = -this->A; this->recenter(); }
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
  if (x<0.0) { this->A = -this->A; this->recenter(); }
  return (*this);
}

IModulo& IModulo::operator/=(const Interval &x) {
  assert_release(!x.contains(0.0) && !x.is_unbounded());
  if (this->is_empty() || this->is_full()) return (*this);
  this->B /= abs(x);
  if (x.ub()<0.0) { this->A = -this->A; this->recenter(); }
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
   return (os << (x.A*x.B) << "+k" << x.B);
}


