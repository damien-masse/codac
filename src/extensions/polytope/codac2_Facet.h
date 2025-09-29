/**
 * \file codac2_Facet.h  
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Damien Massé
 *  \copyright  Copyright 2025
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */

#pragma once

#include <utility> /* uses swap */

#include "codac2_Index.h"
#include "codac2_Vector.h"
#include "codac2_Row.h"
#include "codac2_Interval.h"
#include "codac2_IntervalVector.h"
#include "codac2_IntervalRow.h"
#include "codac2_BoolInterval.h"

namespace codac2 {

/** \brief Potential inclusion relation between sets 
 * INCLUDES : A is included in B
 * MAYINCLUDE : A may be included in B
 * NOTINCLUDE : A\B is not empty
 * INTERSECTS : A inter B is non empty
 * MAYINTERSECT : A inter B may be non empty
 * DISJOINT : A inter B is empty
 */
enum INCLREL {
  INCLUDES,
  MAYINCLUDE,
  NOTINCLUDE,
  INTERSECTS,
  MAYINTERSECT,
  DISJOINT,
  SIZEINCLREL
};
using polytope_inclrel = std::bitset<SIZEINCLREL>;
const polytope_inclrel inclrel_includes = 1<<INCLUDES;
const polytope_inclrel inclrel_mayinclude = 1<<MAYINCLUDE;
const polytope_inclrel inclrel_notinclude = 1<<NOTINCLUDE;
const polytope_inclrel inclrel_intersects = 1<<INTERSECTS;
const polytope_inclrel inclrel_mayintersect = 1<<MAYINTERSECT;
const polytope_inclrel inclrel_disjoint = 1<<DISJOINT;


struct Facet {
   Row row;
   double rhs;
   bool eqcst;

   Facet(const Row &row, double rhs, bool eqcst) :
		 row(row), rhs(rhs), eqcst(eqcst)
   {}

   Facet() {}

   void swap(Facet &f) {
      this->row.swap(f.row);
      std::swap(this->rhs,f.rhs);
      std::swap(this->eqcst,f.eqcst);
   }

   /* test only includes and intersects */
   polytope_inclrel relation_Box(const IntervalVector &b) const {
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

};

}
