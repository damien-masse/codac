/**
 * \file codac2_Facet.cpp
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Damien Massé
 *  \copyright  Copyright 2025
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */


#include <utility> /* uses swap */
#include <stdexcept>

#include "codac2_Index.h"
#include "codac2_Vector.h"
#include "codac2_Row.h"
#include "codac2_Interval.h"
#include "codac2_IntervalVector.h"
#include "codac2_IntervalRow.h"
#include "codac2_BoolInterval.h"
#include "codac2_Facet.h"

namespace codac2 {

namespace Facet_ {

polytope_inclrel relation_Box(const Facet &f, const IntervalVector &b) {
      const Row &row = f.first.row;
      const double &rhs = f.second.rhs;
      const bool &eqcst = f.second.eqcst;
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

void contract_Box(const Facet &f, IntervalVector &b) {
   const Row &row = f.first.row;
   const double &rhs = f.second.rhs;
   /* use MulOp:bwd */
   IntervalRow x1(row);
   MulOp::bwd(Interval(-oo,rhs),x1,b);
   if (b[0].is_empty()) b.set_empty();
}

void contract_out_Box(const Facet &f, IntervalVector &b) {
   const Row &row = f.first.row;
   const double &rhs = f.second.rhs;
   /* use MulOp:bwd */
   IntervalRow x1(row);
   MulOp::bwd(Interval(rhs,oo),x1,b);
   if (b[0].is_empty()) b.set_empty();
}

Interval bound_linear_form(const Facet &f,
		Row &row2, const IntervalVector &b) {
   const Row &row = f.first.row;
   const double &rhs = f.second.rhs;
   const Index &bdim = f.first.bdim;
   const bool &eqcst = f.second.eqcst;
   if (bdim==-1) return row2.dot(b);
   Index abdim = bdim/(2*row.size());
   bool neg = (abdim>=row.size());
   if (neg) { abdim -= row.size(); }
   if (!eqcst && 
       ((!neg && row2[abdim]<=0.0) || (neg && row2[abdim]>=0))) 
	return row2.dot(b);
   Interval quotient = Interval(row2[abdim])/Interval(row[abdim]);
   Interval res=quotient*rhs;
   IntervalRow rem=row2-quotient*row;
   rem[abdim]=0.0;
   return res + rem.dot(b);
}

} /* end of namespace Facet_ */

std::ostream& operator<<(std::ostream& os, const Facet& f) {
   os << f.first.row << (f.second.eqcst ? "=" : "<=" ) << f.second.rhs;
   return os;
}

CollectFacets::CollectFacets(const Matrix &mat, const Vector &rhsvect,
		const std::vector<Index> &eqSet) {
   assert(mat.rows()==rhsvect.size());
   for (Index i=0;i<mat.rows();i++) {
      std::pair<Index,bool> a =this->insert_facet(mat.row(i),rhsvect[i],false);
      if (!a.second) 
	throw std::domain_error("CollectFacets with duplicate rows");
   }
   for (Index i: eqSet) { 
      _eqFacets.push_back(i);
      (_allFacets[i])->second.eqcst=true;
   }
}

std::pair<Index,bool> CollectFacets::result_insertion
	(Index id, std::pair<mapIterator,bool> res,
	double rhs, bool eqcst, ACTION_DUPLICATE act) {
   if (res.second) {
     _allFacets.push_back(res.first);
     if (eqcst) _eqFacets.push_back(id-1);
     return { id, true };
   } else {
      double old_rhs=res.first->second.rhs;
      if (act==CHANGE_RHS ||
	   (act==MAX_RHS && old_rhs<rhs) || (act==MIN_RHS && old_rhs>rhs)) {
           bool old_eqcst = res.first->second.eqcst;
           res.first->second.rhs=rhs;
           res.first->second.eqcst=eqcst;
           if (eqcst && !old_eqcst) 
		_eqFacets.push_back(res.first->second.Id-1);
           else if (!eqcst && old_eqcst) 
		remove_in_eqFacets(res.first->second.Id-1);
           return { res.first->second.Id, false };
      }
   }
   return { 0, false };
}

void CollectFacets::remove_in_eqFacets(Index id) {
   for (Index i=0; i<(Index)_eqFacets.size(); i++) {
      if (_eqFacets[i]==id) {
         _eqFacets[i]=_eqFacets.back();
         _eqFacets.pop_back();
         return;
      }
   }
}

std::pair<Index,bool> CollectFacets::insert_facet(const Row &row, 
		double rhs, bool eqcst, ACTION_DUPLICATE act) {
   Index id = _allFacets.size()+1;
   std::pair<mapIterator,bool> res =
	_map.emplace(std::move(Facet_::make(row,rhs,eqcst,id)));
   return result_insertion(id,res,rhs,eqcst,act);
}
std::pair<Index,bool> CollectFacets::insert_facet(Row &&row, double rhs, bool eqcst, ACTION_DUPLICATE act) {
   Index id = _allFacets.size()+1;
   std::pair<mapIterator,bool> res =
	_map.emplace(std::move(Facet_::make(row,rhs,eqcst,id)));
   return result_insertion(id, res,rhs,eqcst,act);
}
std::pair<Index,bool> CollectFacets::insert_facet(const Facet &fct, ACTION_DUPLICATE act) {
   Index id = _allFacets.size()+1;
   std::pair<mapIterator,bool> res =
	_map.emplace(std::pair(fct.first, FacetRhs(fct.second,id)));
   return result_insertion(id, res,fct.second.rhs,fct.second.eqcst,act);
}

std::pair<CollectFacets::mapIterator,bool>
   CollectFacets::change_eqFacet(Index id, Row &&nrow, double nrhs,
		ACTION_DUPLICATE act) {
   mapIterator fcIt = _allFacets[_eqFacets[id]];
   auto nde = _map.extract(fcIt);
   nde.key().change_row(nrow);
   nde.mapped().rhs=nrhs;
   const auto result_insert = _map.insert(std::move(nde));
   if (!result_insert.inserted) {
      double old_rhs=result_insert.position->second.rhs;
      if (act==CHANGE_RHS ||
	   (act==MAX_RHS && old_rhs<nrhs) || (act==MIN_RHS && old_rhs>nrhs)) {
           bool old_eqcst = result_insert.position->second.eqcst;
           result_insert.position->second.rhs=nrhs;
           result_insert.position->second.eqcst=true;
	   result_insert.position->second.Id=_eqFacets[id]+1;
           _allFacets[_eqFacets[id]] = result_insert.position;
           if (old_eqcst) 
		remove_in_eqFacets(result_insert.position->second.Id-1);
	   _allFacets[result_insert.position->second.Id-1] = endFacet();
           return { result_insert.position, false };
      } else {
         _allFacets[_eqFacets[id]] = this->_map.end();
         _eqFacets[id] = _eqFacets.back();
         _eqFacets.pop_back();
         return { this->_map.end(), false };
      }
   }
   _allFacets[_eqFacets[id]] = result_insert.position;
   return { result_insert.position, true };
} 

std::pair<CollectFacets::mapIterator,bool>
   CollectFacets::change_ineqFacet(Index id, Row &&nrow, double nrhs,
		ACTION_DUPLICATE act) {
   mapIterator fcIt = _allFacets[id-1];
   auto nde = _map.extract(fcIt);
   nde.key().change_row(nrow);
   nde.mapped().rhs=nrhs;
   const auto result_insert = _map.insert(std::move(nde));
   if (!result_insert.inserted) {
      double old_rhs=result_insert.position->second.rhs;
      if (act==CHANGE_RHS ||
	   (act==MAX_RHS && old_rhs<nrhs) || (act==MIN_RHS && old_rhs>nrhs)) {
         bool old_eqcst = result_insert.position->second.eqcst;
         result_insert.position->second.rhs=nrhs;
         result_insert.position->second.eqcst=true;
         if (old_eqcst) 
		remove_in_eqFacets(result_insert.position->second.Id-1);
	 _allFacets[result_insert.position->second.Id-1] = endFacet();
	 result_insert.position->second.Id=id;
         _allFacets[id-1] = result_insert.position;
         return { result_insert.position, false };
      } else {
         _allFacets[id-1] = this->_map.end();
         return { this->_map.end(), false };
      }
   }
   _allFacets[id-1] = result_insert.position;
   return { result_insert.position, true };
} 

CollectFacets::mapIterator CollectFacets::change_ineqFacet_rhs(mapIterator it,
			 double nrhs) {
   it->second.rhs=nrhs;
   return it;
}

CollectFacets::mapIterator 
	CollectFacets::dissociate_eqFacet(Index id, double nbound,
		ACTION_DUPLICATE act) {
   Index aId = _eqFacets[id]; 
		/* start with removing id in _eqFacets */
   _eqFacets[id] = _eqFacets.back();
   _eqFacets.pop_back();
   mapIterator fcIt = _allFacets[aId];
   mapIterator ret = this->_map.end();
   double rhs = fcIt->second.rhs;
   if (nbound<=rhs) {
       fcIt->second.eqcst=false; 
       if (nbound>-oo) {
          std::pair<Index,bool> rinsert =
		 this->insert_facet(-fcIt->first.row,-nbound,false,act);
          if (rinsert.first!=0) ret = _allFacets[rinsert.first-1];
       } else ret = fcIt;
   } else {
       fcIt->second.eqcst=false; 
       if (nbound<+oo) {
          std::pair<Index,bool> rinsert = this->insert_facet(-fcIt->first.row,
				-fcIt->second.rhs,false,act);
          fcIt->second.rhs = nbound;  
          if (rinsert.first!=0) ret = _allFacets[rinsert.first-1];
       } else { /* we need to extract and reinsert the facet */
          auto nde = _map.extract(fcIt);
          nde.key().negate_row();
          double nrhs = -nde.mapped().rhs;
          nde.mapped().rhs= nrhs;
          nde.mapped().eqcst= false;
          const auto result_insert = _map.insert(std::move(nde));
          if (!result_insert.inserted) {
             double old_rhs=result_insert.position->second.rhs;
             if (act==CHANGE_RHS ||
	        (act==MAX_RHS && old_rhs<nrhs) || 
		(act==MIN_RHS && old_rhs>nrhs)) {
                bool old_eqcst = result_insert.position->second.eqcst;
                result_insert.position->second.rhs=nrhs;
                result_insert.position->second.eqcst=false;
                if (old_eqcst) 
	   	   remove_in_eqFacets(result_insert.position->second.Id-1);
	        _allFacets[result_insert.position->second.Id-1] = endFacet();
	        result_insert.position->second.Id=aId+1;
	        _allFacets[aId] = result_insert.position;
                ret=result_insert.position;
             } else {
                _allFacets[aId] = endFacet();
                ret=endFacet();
             }
         } else {
            _allFacets[aId] = result_insert.position;
            ret=result_insert.position;
         }
      }
   }
   return ret;
}

}
