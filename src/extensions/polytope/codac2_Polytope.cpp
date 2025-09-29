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

Polytope::Polytope() : _dim(-1), _nbEqcsts(0), _nbineqCsts(0), _empty(true)
{
}

Polytope::Polytope(Index dim) : Polytope(dim,false)
{
}

Polytope::Polytope(Index dim, bool empty) : _dim(dim), _nbEqcsts(0), _nbineqCsts(0), _empty(empty)
{
}

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

Polytope::Polytope(const IntervalVector &box) : Polytope(box.size())
{
  if (box.is_empty()) { _empty=true; return; }
  _box=box;
  this->build_clpForm_from_box();
}

Polytope::Polytope(const Polytope &P) : _dim(P._dim), _nbEqcsts(P._nbEqcsts),
	_nbineqCsts(P._nbineqCsts), _empty(P._empty),
        _box(P._box) {
   if (P._clpForm) this->_clpForm = std::make_unique<LPclp>(*P._clpForm);

}

Polytope::Polytope(const std::vector<Vector> &vertices) : 
  Polytope()
{
  /* note : for now we just take the bounding box of the vertices,
     so that's definitely not good */
  if (vertices.empty()) return;
  _empty= false;
  _box = *(vertices.begin());
  for (auto v : vertices) {
     _box |= v;
  }
  this->build_clpForm_from_box();
}   

Polytope::Polytope(const std::vector<Vector> &vertices,
		const std::vector<Row> &facetforms) 
		: _dim(facetforms[0].size()),
	_nbEqcsts(0), _nbineqCsts(facetforms.size()),
         _empty(vertices.empty()) {
  if (vertices.empty()) { _nbineqCsts=0; return; }
  _box = *(vertices.begin());
  for (auto v : vertices) {
     _box |= v;
  }
  std::vector<Facet> facets;
  for (const auto &row : facetforms) {
      IntervalRow rI(row); /* use intervals to ensure upward bound */
      double bnd = -oo;
      for (const auto &vtx : vertices) {
          bnd = std::max(bnd,rI.dot(vtx).ub());
      }     
      facets.emplace_back(row,bnd,false);
  }
  _clpForm=std::make_unique<LPclp>(_dim,facets,_box);
}


Polytope::Polytope(const IntervalVector &box,
	const std::vector<Facet> &facets, bool minimize) : Polytope(box.size()) {
   _box=box;
   _clpForm = std::make_unique<LPclp>(_box.size(), facets);
   if (minimize) this->minimize_constraints();
}

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

}

