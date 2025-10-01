/**
 * \file codac2_Facet.h  
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Damien Massé
 *  \copyright  Copyright 2025
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */


#include <utility> /* uses swap */

#include "codac2_Index.h"
#include "codac2_Vector.h"
#include "codac2_Row.h"
#include "codac2_Interval.h"
#include "codac2_IntervalVector.h"
#include "codac2_IntervalRow.h"
#include "codac2_BoolInterval.h"
#include "codac2_Facet.h"

namespace codac2 {

polytope_inclrel Facet::relation_Box(const IntervalVector &b) const {
      if (b.is_empty()) return inclrel_includes | inclrel_disjoint;
      IntervalVector a(b);
      /* check the vertex that maximizes row */
      for (Index i=0;i<row.size();i++) {
          if (row[i]>0) a[i]=a[i].ub(); else a[i]=a[i].lb();
      }
      Interval maxv = row.dot(a)-rhs;
      polytope_inclrel r1=0;
      if (maxv.ub()<=0.0) {
         if (!eqcst) return inclrel_includes | inclrel_intersects;
         else if (maxv.ub()<0.0) return inclrel_notinclude | inclrel_disjoint;
         r1 = inclrel_includes;
      } else if (maxv.lb()<=0.0) {
         r1 = inclrel_mayinclude;
      } else {
         r1 = inclrel_notinclude;
      }
      /* check the vertex that minimizes row */
      for (Index i=0;i<row.size();i++) {
          if (row[i]>0) a[i]=b[i].lb(); else a[i]=b[i].ub();
      }
      Interval minv = row.dot(a)-rhs;
      if (minv.lb()>0.0) return inclrel_notinclude | inclrel_disjoint;
      if (!eqcst) {
         if (minv.ub()>0.0) return r1 | inclrel_mayintersect;
         return r1 | inclrel_intersects;
      } else {
         if (r1[INCLUDES]) { /* maxv.ub == 0 */
              if (minv.lb()>=0.0) return 
				inclrel_includes | inclrel_intersects;
              return (maxv.lb()==0.0 ? inclrel_intersects : inclrel_mayintersect)
                  |  (minv.ub()<0.0 ? inclrel_notinclude : inclrel_mayinclude);
         }
         if (minv.ub()<0.0) r1 = inclrel_notinclude;
         if (maxv.lb()>=0.0 && minv.ub()<=0.0) 
			return r1 | inclrel_intersects;
         else return r1 | inclrel_mayintersect;
      }
} 

void Facet::contract_Box(IntervalVector &b) const {
   /* use MulOp:bwd */
   IntervalRow x1(this->row);
   MulOp::bwd(Interval(-oo,this->rhs),x1,b);
   if (b[0].is_empty()) b.set_empty();
}

void Facet::contract_out_Box(IntervalVector &b) const {
   /* use MulOp:bwd */
   IntervalRow x1(this->row);
   MulOp::bwd(Interval(this->rhs,oo),x1,b);
   if (b[0].is_empty()) b.set_empty();
}

std::ostream& operator<<(std::ostream& os, const Facet& f) {
   os << f.row << (f.eqcst ? "=" : "<=" ) << f.rhs;
   return os;
}

}
