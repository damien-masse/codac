/**
 * \file codac2_Polytope.cpp
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Damien Massé
 *  \copyright  Copyright 2025
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */

#include <vector>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <ctime>
#include <cmath>
#include <bitset>
#include <memory>

#include "codac2_Index.h"
#include "codac2_Matrix.h"
#include "codac2_Vector.h"
#include "codac2_Row.h"
#include "codac2_IntervalRow.h"
#include "codac2_IntervalVector.h"
#include "codac2_IntervalMatrix.h"
#include "codac2_inversion.h"
#include "codac2_clp.h"
#include "codac2_Polytope.h"
#include "codac2_Facet.h"

using namespace codac2;


namespace codac2 {

Polytope::Polytope() : _dim(-1), _empty(true)
{
}

Polytope::Polytope(Index dim) : Polytope(dim,false)
{
}

Polytope::Polytope(Index dim, bool empty) : _dim(dim), _empty(empty),
		_box(dim, empty ? Interval::empty() : Interval()),
		_facets(std::make_shared<CollectFacets>(dim)),
		_box_updated(true)
{
}

#if 0
void Polytope::build_clpForm_from_box() {
  _nbEqcsts=0;
  _nbineqCsts=0;
  _dim=_box.size();
  std::vector<Facet> built;
  for (Index i=0;i<_box.size();i++) {
    Row r = Row::zero(_box.size());
    r[i] = 1.0;    
    if (_box[i].is_degenerated()) {
      _nbEqcsts++;
      built.emplace_back(r,_box[i].mid(),true);
    }
    else {
     if (_box[i].ub()<oo) {
       _nbineqCsts++;
       built.emplace_back(r,_box[i].ub(),false);
     }
     if (_box[i].lb()>-oo) {
       _nbineqCsts++;
       built.emplace_back(r,_box[i].ub(),false);
     }
    }
  }
  _clpForm=std::make_unique<LPclp>(_box.size(),built);
}
#endif

Polytope::Polytope(const IntervalVector &box) : _dim(box.size()),
		_empty(box.is_empty()), _box(box), 
		_facets(std::make_shared<CollectFacets>(dim)), 
		_box_updated(true)
{
}

Polytope::Polytope(const Polytope &P) : _dim(P._dim), _empty(P._empty),
        _box(P._box),
	_facets(P._empty ? std::make_shared<CollectFacets>(dim) :
		           std::make_shared<CollectFacets>(P._facets)),
	_box_updated(P._box_updated) 
{
}

Polytope::Polytope(const std::vector<Vector> &vertices) : 
  Polytope()
{
   if (vertices.empty()) return;
   _dim=vertices[0].size();
   _empty= false;
   /* build the V2F form */
   _DDbuildV2F = std::make_unique<DDbuildV2F>(1,vertices[0]);
   for (Index i=1;i<(Index)vertices.size();i++) {
       _DDbuildV2F->add_vertex(i+1,vertices[i]);
   }
   /* get the CollectFacets */
   _facets=_DDbuildV2F->getFacets();
   _box = _facets->extractBox();
   _facets->include_vertices(vertices,_box,true);
   _box_updated=true;
   _buildV2F_updated=true; /* keep buildV2F ? */
}   

Polytope::Polytope(const std::vector<IntervalVector> &vertices,
		const CollectFacets &facetsform) :  Polytope()
{
  _dim = facetsform.getDim();
  if (_dim==-1) return;
  _facets=std::make_shared<CollectFacets>(facetsform);
  if (vertices.empty()) {
     _box = InvervalVector::constant(_dim,Interval::empty());
     _empty=true;
     return;
  }
  _facets=_DDbuildV2F->getFacets();
  _box = _facets->extractBox();
  _facets->include_vertices(vertices,_box,true);
  _box_updated=true;
}


Polytope::Polytope(const IntervalVector &box, 
		const std::vector<std::pair<Row,double>> &facets,
		bool minimize) : Polytope(box) {
   for (const std::pair<Row,double> &cst : facets) {
      this->add_constraint(cst);
   }
   if (minimize) this->minimize_constraints();
}


Interval Polytope::fast_bound(const FacetBase &base) const {
   if (base.isNull()) return Interval::zero();
   Index gdim = base.gtDim();
   Interval res;
   if (base.row[gdim]>0.0) {
      res= Interval(_box[gdim].ub())*base.row[gdim];
   } else {
      res= Interval(-_box[gdim].lb())*base.row[gdim];
   }
   if (!base.isCoord()) {
      Row bcopy = base.row;
      bcopy[gdim]=0.0;
      res += bcopy.dot(_box);
      std::pair<CollectFacets::mapIterator,CollectFacets::mapIterator>
         pIt = _facets->equal_range(base);
      if (pIt.first!=pIt.end) { /* we found an exact key */
         return Interval(pIt.first->second.rhs);
      } else {
         if (pIt.first!=_facets->get_map().begin()) --pIt.first;
         Interval a1 = bound_linear_form(*(pIt.first), base.row, _box);
         if (pIt.second!=_facets->endFacet()) {
             Interval a2 = bound_linear_form(*(pIt.second), base.row, _box);
             if (a2.ub()<a1.ub()) a1=a2;
         }
         if (a1.ub()<res.ub()) res=a1;
      }
      /* check the equalities ? */
      for (Index i=0;i<_facets->nbeqfcts();i++) {
         CollectFacets::mapIterator it = _facets.get_eqFacet(i);
         Interval a = bound_linear_form(*(pIt.first), base.row, _box);
         if (a.ub()<res.ub()) res=a;
      }
   }
   /* check the vertices if _buildF2V is up-to-date */
   if (_DDbuildF2V && _buildF2V_updated) {
       const std::forward_list<DDvertex> &vertices=_DDbuildF2V->get_vertices();
       for (const DDvertex &vt : vertices) {
            IntervalVector vtB = _DDbuildF2V->compute_vertex(vt);
            Interval a = base.row.dot(vtB);
            if (a.ub()<res.ub()) res=a;
       }
   }
   return res;
}

#if 0

Polytope::Polytope(const IntervalMatrix &M, const IntervalMatrix &rM,
	const IntervalVector &V) :
	Polytope(V.size(),V.is_empty()) {
   if (V.is_empty()) return;
   /* rM x \in V */
   _box = M*V;
   std::vector<Facet> built;
   for (Index i=0;i<_dim;i++) {
      IntervalRow rowI = rM.row(i);
      if (rowI.is_unbounded()) continue;
      Row rowIm = rowI.mid();
      Interval rhs = V[i] - (rowI-rowIm)*_box;
      if (rhs.is_degenerated()) {
        _nbEqcsts++;
         built.emplace_back(rowIm,rhs.mid(),true);
      } else {
         if (rhs.ub()<oo) {
           _nbineqCsts++;
            built.emplace_back(rowIm,rhs.lb(),false);
         }
         if (rhs.lb()>-oo) {
           _nbineqCsts++;
            built.emplace_back(-rowIm,-rhs.lb(),false);
         }

     }
  }
}

Polytope::Polytope(const IntervalMatrix &M, const IntervalVector &V) :
	Polytope(M,inverse_enclosure(M),V) {
}

Polytope::~Polytope() {
}

void Polytope::minimize_constraints() { /* TODO */ }

const IntervalVector& Polytope::update_box() {
   if (_dim==-1 || _empty) return this->_box;
   if (!_clpForm) return this->_box;
   Row obj = Row::zero(_dim);
   for (Index i=0;i<_dim;i++) {
       obj[i]=1.0;
       _clpForm->setObjective(obj);
       LPclp::lp_result_stat ret =  _clpForm->solve(false,4);
       if (ret[LPclp::EMPTY]) { _empty=true; _box.set_empty(); return _box; }
       if (ret[LPclp::BOUNDED]) _box[i] = min(_box[i],_clpForm->getValobj());
       obj[i]=-1.0;
       _clpForm->setObjective(obj);
       ret =  _clpForm->solve(false,4);
       if (ret[LPclp::EMPTY]) { _empty=true; _box.set_empty(); return _box; }
       if (ret[LPclp::BOUNDED]) _box[i] = max(_box[i],-_clpForm->getValobj());
   }
   _clpForm->set_bbox(_box);
   return this->_box;
}

polytope_inclrel Polytope::relation_Box(const IntervalVector& p) const {
   assert_release(_dim>=0);
   if (p.is_empty()) return (1<<INCLUDES) | (1<<DISJOINT);
   if (_empty) return (1<<NOTINCLUDE) | (1<<DISJOINT);
   if (!_clpForm) return (1<<INCLUDES) | (1<<INTERSECTS);
   /* to get INCLUDES, must have INCLUDES for all 
      to get NOTINCLUDE, must have NOTINCLUDE once
      to get INTERSECTS, must have INCLUDES for all but one, and INTERSECTS
                         for this one 
      to get DISJOINT, must have DISJOINT once */
   bool notinclude=false;
   bool allinclude=true;
   bool intersects=true;
   const std::vector<Facet>& facets = _clpForm->getFacets();
   for (const Facet &facet : facets) {
       polytope_inclrel res2 = facet.relation_Box(p);
       if (res2[DISJOINT]) return (1<<NOTINCLUDE) | (1<<DISJOINT);
       if (res2[INCLUDES]) continue;
       if (allinclude && res2[INTERSECTS]) { allinclude=false; }
       else { allinclude=false; intersects=false; }
       if (res2[NOTINCLUDE]) notinclude=true;
   }
   return (allinclude ? (1<<INCLUDES | 1<<INTERSECTS)
	: ((intersects ? 1<<INTERSECTS : 1<<MAYINTERSECT) |
           (notinclude ? 1<<NOTINCLUDE : 1<<MAYINCLUDE)));
}

BoolInterval Polytope::contains(const IntervalVector& p) const {
   polytope_inclrel r = this->relation_Box(p);
   return (r[INCLUDES] ? BoolInterval::TRUE : 
		(r[NOTINCLUDE] ? BoolInterval::FALSE : BoolInterval::UNKNOWN));
}

BoolInterval Polytope::intersects(const IntervalVector& p) const {
   polytope_inclrel r = this->relation_Box(p);
   return (r[INTERSECTS] ? BoolInterval::TRUE : 
		(r[DISJOINT] ? BoolInterval::FALSE : BoolInterval::UNKNOWN));
}

const Interval& Polytope::operator[](Index i) const {
   assert_release(i<=_dim);
   return _box[i];
}

Vector Polytope::mid() const {
   assert_release(_dim>=0);
   if (_empty) return 
	Vector::Constant(_dim,std::numeric_limits<double>::quiet_NaN());
   if (!_clpForm) return _box.mid();
   LPclp::lp_result_stat stat = _clpForm->solve(true,0);
   if (stat[LPclp::NOTEMPTY] || stat[LPclp::NOTEMPTY_APPROX])
       return _clpForm->getFeasiblePoint().mid();
   return Vector::Constant(_dim,std::numeric_limits<double>::quiet_NaN());
}

double Polytope::distance_cst(const Facet &fc) const {
   assert_release(_dim==fc.row.size());
   if (_empty) return -oo;
   if (!_clpForm) return (fc.row.dot(_box)-fc.rhs).ub();
   _clpForm->setObjective(fc.row);
   LPclp::lp_result_stat stat = _clpForm->solve(false,0);
   if (stat[LPclp::EMPTY]) return -oo;
   return (_clpForm->getValobj()-fc.rhs).ub();
}

bool Polytope::box_is_included(const IntervalVector& x) const {
   return _box.is_subset(x);
}

void Polytope::set_empty() {
   _empty=true;
   _box.set_empty();
   _nbEqcsts=_nbineqCsts=0;
   _clpForm.reset(nullptr);
}

void Polytope::clear() {
   assert_release(_dim>=0);
   _empty=false;
   _box = IntervalVector::Zero(_dim);
   this->build_clpForm_from_box();
}

#endif

bool Polytope::add_constraint(const std::pair<Row,double>& facet, 
	double tolerance) {
   if (_empty) return false;
   FacetBase base = FacetBase(facet.first);
   if (base.isNull()) {
      if (facet.second>=tolerance) return false;
      this->set_empty(); return true;
   }
   Interval act = this->fast_bound(base);
   if (act.ub()<=facet.second.tolerance) return false;
   if (base.isCoord()) {
      Index gdim = base.gtDim();
      Interval val = Interval(facet.second)/row[gdim];
      if (row[gdim]>0.0) {
         _box[gdim] &= Interval(-oo,val.ub());
         if (_box[gdim].is_empty()) { 
	   this->set_empty(); 
           return true;
	 }
      } else {
         _box[gdim] &= Interval(-val.ub(),+oo);
         if (_box[gdim].is_empty()) { 
	   this->set_empty(); 
           return true;
	 }
      }
      /////////////// update ?
      return true; 
   }
   auto res= _facets->insert_facet(facet.first, facet.second, 
		false, CollectFacets::MIN_RHS);
   //////// update ?
   return res.second;
}

}
#endif

