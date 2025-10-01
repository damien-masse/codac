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

   inline void swap(Facet &f) {
      this->row.swap(f.row);
      std::swap(this->rhs,f.rhs);
      std::swap(this->eqcst,f.eqcst);
   }

   /* test only includes and intersects */
   polytope_inclrel relation_Box(const IntervalVector &b) const;

   void contract_Box(IntervalVector &b) const;
   void contract_out_Box(IntervalVector &b) const;
   friend std::ostream& operator<<(std::ostream& os, const Facet& f); 


};

}
