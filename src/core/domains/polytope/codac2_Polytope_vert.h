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
   IntervalVector cone;
   Index level;
   Interval lambda;

   /** \brief build the initial vertex,
    *  vertex = (1,0...}, cone = (0,...,0)
    *  level = 0
    *  @param dim the dimension of the space (=vertex.size-1)
    */
   DDvertex_base(Index dim);

   /** \brief build a new ray
    *  vertex = ray, cone = (0,...,0)
    *  @param level the level of the ray
    *  @param ray the ray (as a dim vector, built from the existing lines)
    */
   DDvertex_base(Index level, const IntervalVector &ray, 
			const IntervalMatrix *A);

   /** \brief build a new intersection
    *  vertex = ray, cone = (0,...,0)
    *  @param level the level of the ray
    *  @param ray the ray (as a dim vector, built from the existing lines)
    */
   DDvertex_base(const DDvertex_base &v1, const DDvertex_base &v2, 
		const IntervalRow &facet);

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

  /** \brief intersects the vertex with a constraint :
   *   - compute the "closest" vertex satisfying the constraint in the cone
   *   - compute the new "direction" for the cone
   *  return also an integer with the state of the result
   *   0 : returned a new vertex, no change in direction
   *   1 : direction restricted
   *  -1 : can't satisfy the constraint.
   * Note: lambda is already computed , and is supposed to have .lb<0 */
   std::tuple<int,IntervalVector,IntervalVector>
		compute_closest(const IntervalRow &facet) const;


   /** \brief compute the unit vector v maximising the direction of a
    *  constraint rho
    *  and satisfying A^T v = 0 (with A a matrix of orthogonal unit vectors)
    *  If A is not defined, returns d_rho/||d_rho||
    *  otherwise, returns (Id - AA^T)rho / ||(Id - AA^T) rho||
    *  (also, gives the norm ||(Id - AA^T) rho|| )
    */
   static IntervalVector max_unit_vector (const IntervalRow &rho,
			 const IntervalMatrix *A, Interval &norm);

   private:
      const IntervalMatrix *proj; /* projection matrix for rays */

};

inline IntervalVector DDvertex_base::max_unit_vector(const IntervalRow &rho,
			const IntervalMatrix *A, Interval &norm) {
    IntervalVector ret(rho.transpose());
    if (A==nullptr) {
         ret[0]=0.0;
    } else { 
         ret = (*A)*ret;
	 std::cout << "ret " << ret << "\n";
    }
    norm = ret.norm();
    std::cout << "norm " << norm << "\n";
    return ret/norm;
}

inline DDvertex_base::DDvertex_base(Index dim) :
	vertex(IntervalVector::Zero(dim+1)), cone(IntervalVector::Zero(dim+1)),
	level(0), proj(nullptr)
{
   this->vertex[0]=1.0;
}

inline DDvertex_base::DDvertex_base(Index level, const IntervalVector &ray,
			const IntervalMatrix *proj) :
	vertex(ray),
        cone(IntervalVector::Zero(ray.size())),
	level(level), proj(proj)
{
}

inline DDvertex_base::DDvertex_base(const DDvertex_base &v1,
	 const DDvertex_base &v2, const IntervalRow &facet) {
   /* vertex, cone, level, lambda, proj */
    /* vertex : max v1, min v2 */
    this->vertex = -v2.lambda.lb()*v1.vertex + v1.lambda.lb()*v2.vertex;
    reduce_vector_without_rounding(this->vertex);
    /* level : min v1,v2 */
    std::cout << v1.cone << " " << v2.cone << "\n";
    std::cout << v1.level << " " << v2.level << "\n";
    std::cout << facet << "\n";
    this->level=std::min(v1.level,v2.level);
    this->proj = (v1.level<=v2.level ? v1.proj : v2.proj);
    Interval norm_bfacet;
    IntervalVector vmax = max_unit_vector(facet,this->proj,norm_bfacet);
    Interval rhoc;
    IntervalVector mcone;
    if (v1.level<v2.level) { /* note : we known v2.lambda.lb()<0.0 */
       /* a new cone can be build with v1.cone + alpha * v2.vertex
          with 0 <= alpha <= (norm_bfacet - rhoc)/lambda
          if rhoc <= norm_bfacet, attempt with v1 cone */
       if (v1.cone.is_empty()) { this->cone = v1.cone; return; }
       rhoc = facet.dot(v1.cone); 
       Interval calc = -rhoc+norm_bfacet;
       if (calc.ub()<0.0) {
          Interval alpha = calc/v2.lambda.lb();
          this->cone = v1.cone + alpha.lb()*v2.vertex;
          std::cout << "mix cone : " << this->cone << "\n";
          return;
       }
       mcone = v1.cone;
    } else if (v1.level>v2.level) {
       /* note : we known v1.lambda.lb()>0.0 */
       /* a new cone can be build with v2.cone + alpha * v1.vertex
          with alpha >= (norm_bfacet -rhoc)/lambda and alpha>=0.0
          if rhoc >= norm_bfacet, alpha = 0 => don't move. */
       if (v2.cone.is_empty()) { this->cone = v2.cone; return; }
       rhoc = facet.dot(v2.cone); 
       Interval calc = -rhoc+norm_bfacet;
       if (calc.ub()<=0.0) { 
           this->cone = v2.cone; 
           std::cout << "mix cone : " << this->cone << "\n";
           return;
       }
       Interval alpha = calc/v1.lambda.lb();
       this->cone = v2.cone + alpha.lb()*v1.vertex;
       std::cout << "mix cone : " << this->cone << "\n";
       return;
    } else {
       Interval r1 = facet.dot(v1.cone);
       Interval r2 = facet.dot(v2.cone);
    std::cout << "rho1 : " << r1 << "\n";
    std::cout << "rho2 : " << r2 << "\n";
       if (r1.is_empty()) {
         mcone=v2.cone;
         rhoc = r2;
       } else if (r2.is_empty()) {
         mcone=v1.cone;
         rhoc = r1;
       } else {
         Interval alpha = (norm_bfacet-r2)/(r1-r2);
         if (alpha.lb()<0.0) {
            mcone=v2.cone; rhoc=r2;
         } else if (alpha.ub()>1.0) {
            mcone=v1.cone; rhoc=r1;
         } else {
            mcone=(alpha*v1.cone+(1.0-alpha)*v2.cone);
            rhoc=facet.dot(mcone);
         }
       }
    } 
    std::cout << "rhoc : " << rhoc << "\n";
    if (rhoc.is_empty()) {
        this->cone = IntervalVector::empty(this->vertex.size());
        return;
    }
    if (rhoc.lb()>=norm_bfacet.ub())  {
       this->cone=mcone;
    } else if (rhoc.lb()<=-norm_bfacet.lb()) {
       this->cone = IntervalVector::empty(this->vertex.size());
    } else {
       Interval ratio = rhoc/norm_bfacet;
       Interval mu = -1.0 + 2.0/(1.0+ratio);
       Interval beta = mu+1;
       this->cone = beta.ub()*mcone+mu.ub()*vmax;
    }
    std::cout << "mix cone : " << this->cone << "\n";
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
   /* the cone :  bfacet is the base of the facet (without constant part)
      we want bfacet.(ncone+B_n) >= 0
      i.e. with bfacet.ncone >= ||bfacet||
      we consider ncone=lambda.cone + mu.line, lambda>=1
      then lambda bfacet.cone + mu bfacet.line >= ||bfacet||
      with lambda=1, we get:
      mu = (||bfacet|| - bfacet.cone)/prod (once again, prod>0)
      once again, we consider the ub the result */
   Interval norm_bfacet;
   std::cout << "part1a : " << facet << "\n";
   if (proj) std::cout << *proj << "\n"; else std::cout << "nullptr\n";
   max_unit_vector(facet,proj,norm_bfacet); /* check : nullptr ? */
   std::cout << "part2 : " << norm_bfacet << "\n";
   Interval mu = (norm_bfacet - facet.dot(this->cone))/prod;
   this->cone += (signprod ? mu.ub() : mu.lb())*line;
   std::cout << "cone : " << this->cone << "\n";
}


inline bool DDvertex_base::is_ray() const { return this->level>0; }

inline void DDvertex_base::compute_lambda(const IntervalVector &facet) {
    this->lambda= facet.dot(this->vertex);
}

inline std::tuple<int, IntervalVector, IntervalVector>
	 DDvertex_base::compute_closest(const IntervalRow &facet) const {
   if (this->cone.is_empty()) 
      return { -1, this->vertex, this->cone };
   Interval norm_bfacet;
   IntervalVector vmax = max_unit_vector(facet,proj,norm_bfacet);
   /* for newv = vertex + alpha (cone + b) we want
    * facet newv >= 0  and alpha as small as possible
    * we get alpha >=  (-lambda)/(facet (cone+b)) with b maximising facet.b
    */
   Interval rhoc = facet.dot(cone);
   std::cout << " try_extend : " << facet << " " << cone << "\n";
   std::cout << rhoc << " " << norm_bfacet << "\n";
   if ((rhoc+norm_bfacet).lb()<=0.0) return { -1, this->vertex, this->cone };
   IntervalVector newv = this->vertex + 
			(-lambda/(rhoc+norm_bfacet)).ub() * (cone+vmax);
      /* FIXME : can we use vertex.mid ? as lambda is used to handle any
	   value in this->vertex */
   /* for newc = beta cone + mu b we want
    * facet.(newc-b)>= 0 , so beta facet.cone + (mu-1) facet.b >= 0
    * and beta >= mu+1 (ensure inclusion of new cone)
    * we consider ratio = facet.c / facet.b (with b maximising facet.b)
    * we get mu >= -1 + 2/(1+ratio)
    * and beta >= mu+1 
    */
   if (rhoc.lb()>=norm_bfacet.ub()) 
	return { 0, newv, this->cone };
   Interval ratio = rhoc/norm_bfacet;
   Interval mu = -1.0 + 2.0/(1.0+ratio);
   Interval beta = mu+1;
   IntervalVector newc = beta.ub()*this->cone+mu.ub()*vmax;
      /* FIXME : can we use cone.mid ? as rhoc is used to handle any
	   value in this->cone */
   std::cout << "newc : " << newc << "\n";
   return { 1, newv, newc };
}

}
