/**
 * \file codac2_Polytope_vertex.h vertex (box + cone)
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Damien Massé
 *  \copyright  Copyright 2026
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */

#pragma once

#include <vector>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <cmath>
#include <bitset>
#include <memory>
#include <tuple>
#include <stdexcept>
#include <forward_list>

#include "codac2_Index.h"
#include "codac2_Vector.h"
#include "codac2_Row.h"
#include "codac2_Facet.h"

namespace codac2 {

/** \brief multiply a vector by a power of 2 (without rounding)
 *  such that 2^(-5) <= norm(vec) <= 2^5
 *  \param vec the vector (modified)
 *  \return the multiplying factor */
inline double reduce_vector_without_rounding(IntervalVector &vec) {
      int fx;
      if (frexp(vec[0].mag(),&fx)==0.0) fx=INT_MIN;
      for (Index i=1;i<vec.size();i++) {
         int fx2;
         if (frexp(vec[i].mag(),&fx2)==0.0) continue;
         if (fx2>fx) fx=fx2;
      }
      assert_release (fx!=INT_MIN);
      if (std::abs(fx)<5) return 1.0;
      double mult=exp2(-fx);
      for (Index i=0;i<vec.size();i++) {
//           vec[i] = Interval(ldexp(vec[i].lb(),-fx),ldexp(vec[i].ub(),-fx));
         vec[i] *= mult;
      }
      return mult;
}

/** \brief multiply a row by a power of 2 (without rounding)
 *  such that 2^(-5) <= norm(vec) <= 2^5
 *  \param vec the row (modified)
 *  \return the multiplying factor */
inline double reduce_row_without_rounding(Row &vec) {
      int fx;
      if (frexp(vec[0],&fx)==0.0) fx=INT_MIN;
      for (Index i=1;i<vec.size();i++) {
         int fx2;
         if (frexp(vec[i],&fx2)==0.0) continue;
         if (fx2>fx) fx=fx2;
      }
      assert_release (fx!=INT_MIN);
      if (std::abs(fx)<5) return 1.0;
      double mult=exp2(-fx);
      for (Index i=0;i<vec.size();i++) {
//           vec[i] = ldexp(vec[i],-fx);
         vec[i] *= mult;
      }
      return mult;
}

/** \brief a vertex has a point (IntervalVector) and a admissible cone for
 *  ``extension''. The cone is defined using an IntervalVector (one point
 *  is guaranteed to be "ok", i.e    vertex + lambda (c+B_n)
 *  A vertex is also associated to a set of constraints :
 *      * "positive" constraints
 *      * "zero" constraints (for rays), intersected with the vertex
 *      * (zero or one) "negative" constraint
 *      * (rays have one negative constraint)
 *  positive and zero constraints are stored in a set
 *  a "ray" level gives the number of zero constraints
 */
struct DDvertex_base {
   IntervalVector vertex;
   Interval lambda;

   /** \brief build the initial vertex,
    *  vertex = (1,0...}
    *  @param dim the dimension of the space (=vertex.size-1)
    */
   DDvertex_base(Index dim);

   /** \brief build a new ray
    *  vertex = ray
    *  @param ray the ray (as a dim vector, built from the existing lines)
    */
   DDvertex_base(const IntervalVector &ray);

   /** \brief build a new intersection
    *  @param v1 vertex1 (outside)
    *  @param v2 vertex2 (inside)
    */
   DDvertex_base(const DDvertex_base &v1, const DDvertex_base &v2);

   /** use a line to make the vertex/ray satisfy a constraint 
    *  @param line the line used
    *  @param facet the constraint which must be satisfied (facet X >= 0)
    *  @param prod the product bfacet.line (must be >0), to avoid 
    *  recomputing it
    */ 
   void add_constraint_line(const IntervalVector &line,
	const IntervalRow &facet, const Interval &prod);

   /** \brief is the vertex a ray (first component is 0) */
   bool is_ray() const;

   /** \brief with a constraint, compute lambda the gap between the 
    * vertex and the constraint 
    * warning : lambda is in fact v[0]*lambda when it is a vertex */
   void compute_lambda(const IntervalVector &facet);

  /** \brief "extend" if possible a vertex to satisfy a constraint
   *   - compute the "closest" vertex satisfying the constraint,
   *     using interior Point 
   *  return also an double with the "contribution" of the interior point
   *  (-1 means not possible)
   * Note: lambda is already computed , and is supposed to have .lb<0 */
   std::pair<double,IntervalVector>
		compute_closest(const IntervalRow &facet,
			const IntervalVector &intPoint) const;
};

inline DDvertex_base::DDvertex_base(Index dim) :
	vertex(IntervalVector::Zero(dim+1))
{
   this->vertex[0]=1.0;
}

inline DDvertex_base::DDvertex_base(const IntervalVector &ray):
	vertex(ray)
{
    reduce_vector_without_rounding(this->vertex);
}

inline DDvertex_base::DDvertex_base(const DDvertex_base &v1,
	 const DDvertex_base &v2) {
   /* vertex, cone, level, lambda, proj */
    /* vertex : max v1, min v2 */
    this->vertex = -v2.lambda.lb()*v1.vertex + v1.lambda.lb()*v2.vertex;
    reduce_vector_without_rounding(this->vertex);
} 

inline void DDvertex_base::add_constraint_line(const IntervalVector &line,
	const IntervalRow &facet, const Interval &prod) {
   bool signprod=prod.lb()>0.0;
   /* first : vertex. we want facet.(vertex+alpha line)>=0 
      ie alpha >= -(facet.vertex)/prod  (provided prod>0)
      we consider ub of the result (facet.vertex.lb, and
           prod.lb if facet.vertex.lb<0, prod.ub if facet.vertex.lb>0)
    */
   Interval alpha = -(facet.dot(this->vertex)/prod);
   this->vertex += (signprod ? alpha.ub() : alpha.lb())*line;
   reduce_vector_without_rounding(this->vertex);
}


inline bool DDvertex_base::is_ray() const { return this->vertex[0].ub()==0.0; }

inline void DDvertex_base::compute_lambda(const IntervalVector &facet) {
    this->lambda= facet.dot(this->vertex);
}

inline std::pair<double, IntervalVector>
	 DDvertex_base::compute_closest(const IntervalRow &facet,
			const IntervalVector &intPoint) const {
   if (intPoint.is_empty()) return { -1.0, this->vertex };
   if (this->is_ray() && intPoint[0].ub()>0.0) return { -1.0, this->vertex };
   Interval facetPoint = facet.dot(intPoint);
   if (facetPoint.ub()>=0.0) /* failed */
         return { -1.0, this->vertex };
      /* new vertex = vertex - mu IntPoint
	 facet . (vertex - mu intPoint) >= 0 means
             mu >= facet.vertex / facetPoint
         furthermore, if no ray, we want (new vertex)[0]>=0.0
	 so vertex[0] - mu intPoint[0] >= 0.0
             so mu <= vertex[0]/intPoint[0] */
   double mu = (this->lambda / facetPoint).ub();
   IntervalVector ret = this->vertex - mu*intPoint;
   if (ret[0].lb()<0.0 || (ret[0].lb()==0.0 && this->vertex[0].ub()>0.0))
       return { -1.0, this->vertex };
   reduce_vector_without_rounding(ret);
   return { mu, ret };
}

}
