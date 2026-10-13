/**
 * \file codac2_Polytope_dd.h implementation of the dd algorithms
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
#include <iomanip>
#include <iostream>
#include <ctime>
#include <cmath>
#include <bitset>
#include <memory>
#include <map>
#include <forward_list>
#include <set>
#include <unordered_set>
#include <stdexcept>

#include "codac2_Index.h"
#include "codac2_Vector.h"
#include "codac2_Row.h"
#include "codac2_Facet.h"
#include "codac2_PosetMax.h"
#include "codac2_Polytope_vert.h"


namespace codac2 {


/** \brief structure for vertices in BuildF2V */
struct DDvertex : DDvertex_base {
   Index Id;

   std::unordered_set<DDvertex*> links; /* bonne méthode ? */
   PosetElem facets;
	/* note : double redirection here, a ``real'' facet may have several
           Ids */
   
   enum ddvertexstatus {
     VERTEXNEW, /* created */
     VERTEXLT, /* < (GT,LT) => new vertex */
     VERTEXON, /* <,=,> (GT,LE) => for now, equivalent to LT */
     VERTEXGE, /* = and possible small >
		  keep the vertex. no change with (GT,GE) */
     VERTEXGT, /* vertex GT, untreated and not in stack */
     VERTEXREM /* vertex GT treated, to be removed */
   };
   ddvertexstatus status;

   /* constructor */
   DDvertex(Index id, Index dim) : DDvertex_base(dim), Id(id), status(VERTEXNEW) { }
   DDvertex(Index id, const IntervalVector &v, Index vdim, const IntervalMatrix *proj) : DDvertex_base(vdim,v,proj), Id(id), status(VERTEXNEW) { }
   DDvertex(Index id, const DDvertex &v1, DDvertex &v2, const IntervalRow &facet) : DDvertex_base(v1,v2,facet), Id(id), facets(v1.facets & v2.facets),
    status(VERTEXNEW) { }  /* links are added after */

   /* should not be copyable */
   DDvertex(const DDvertex &v) = delete;
   DDvertex &operator=(const DDvertex &v) = delete;

   /* possibility of extension along the cone ? */
   bool try_extend(const IntervalRow &facet);

   inline void reduce_vector() {
      reduce_vector_without_rounding(vertex);
   }

   inline bool is_ray() {
      return (this->vertex[0].lb()<=0.0);
   }

   inline bool add_link(DDvertex *a) {
      return links.insert(a).second;
   }

//   void addLnkVertex(const std::forward_list<DDvertex>::iterator &V1);

   inline bool remove_link(DDvertex *a) {
      return links.erase(a)==1;
   }

   inline void change_link(DDvertex *a, DDvertex *b) {
      links.erase(a);
      links.insert(b);
   }

   inline bool add_facet(Index a) {
      return facets.add(a);
   }

   inline const PosetElem& get_facets() const {
      return this->facets;
   }
 
};

/** \brief Structure to build an enclosing set of IntervalVector
 *  from a set of facets */
class DDbuildF2V {
  public:
     /* create the build with a collection of facet, but
        use only the equality facets */
     DDbuildF2V(Index dim, const IntervalVector &box,
 		std::shared_ptr<CollectFacets> facetsptr,
		bool include_box=true);
     
     /* add an inequality facet in the build */
     /* return 1 : changed, 0 : not changed, -1 : empty */
     int add_facet(Index id);
     int add_facet(CollectFacets::mapCIterator itFacet);

     int add_constraint_box(const IntervalVector &box);
     int add_bound_var(Index c, bool mx, double rhs);

     IntervalVector compute_vertex(const IntervalVector &vect) const;
     Index get_dim() const;
     Index get_fdim() const;

     bool is_empty() const;
     const IntervalMatrix &get_M_EQ() const;

     /* return a set of the Id of (detected as) redundant facets ,
        and remove these facets of the representation */
     std::set<Index> redundantFacets();

     /* update reffacets following a renumbering */
     void update_renumber(const std::vector<Index> &ref);

     /* return the bounding box of all vertices */
     IntervalVector build_bbox() const;

     const std::vector<IntervalVector> &get_lines() const;
     const std::forward_list<std::shared_ptr<DDvertex>> &get_vertices() const;
     const std::vector<Index> &get_reffacets() const;

     friend std::ostream& operator<<(std::ostream& os, const DDbuildF2V& build);

     static double tolerance;

  private:
     int add_facet(Index idFacet, const IntervalRow &facet);
     int add_facet_internal(Index idFacet, const IntervalRow &facet);
     int add_facet_new_dimension(Index idFacet, const IntervalRow &facet, Index bestline);
     static bool add_link(DDvertex *V1, DDvertex *V2);
     void add_zero(const IntervalRow &facet);
     static bool find_path(DDvertex *v1, DDvertex *v2, 
			const PosetElem &fcts);

     DDvertex* addDDvertex(Index fdim);
     DDvertex* addDDvertex(const IntervalVector &v, Index vdim); 
     DDvertex* addVertexSon(const DDvertex &p1,
		 DDvertex &p2, const IntervalRow &facet, Index idFacet); 
/*     void removeVertex(Index Id, bool removeLnks); */

     Index dim;
     Index fdim; /* dimension - nb independant equalities constraints */
     Index vdim; /* actual dimension without the lines */
     IntervalVector bbox;

     std::shared_ptr<CollectFacets> facetsptr;

     /* if fdim<dim :
        constraint a X <= b becomes (-b 0) + (a M_EQ) 
        vertex : (c Y) becomes X = M_EQ (1 Y/c) 
        ray : (0 Y) becomes a ray X = M_EQ (0 Y) ?
        lines : (0 Y) becomes a line X = M_EQ (0 Y) ?
     */
     IntervalMatrix M_EQ;

     std::vector<IntervalVector> lines; /* lines (for unbounded polyhedron) */
               /* rays can be handled as classical vertices */
     std::vector<IntervalMatrix> projection_rays; 
	/* projection matrices for the extension cones of rays 
	   each matrix is orthogonal */

     std::forward_list<std::shared_ptr<DDvertex>> vertices;
     std::vector<Index> reffacets; /* indirection for facets */

     bool empty, flat;
     Index nbIn;
};

inline Index DDbuildF2V::get_dim() const { return this->dim; }
inline Index DDbuildF2V::get_fdim() const { return this->fdim; }
inline bool DDbuildF2V::is_empty() const { return this->empty; }
inline const IntervalMatrix &DDbuildF2V::get_M_EQ() const { return this->M_EQ; }
inline const std::vector<IntervalVector> &DDbuildF2V::get_lines() const 
	{ return this->lines; }
inline const std::forward_list<std::shared_ptr<DDvertex>> &DDbuildF2V::get_vertices() const 
	{ return this->vertices; }
inline const std::vector<Index> &DDbuildF2V::get_reffacets() const 
	{ return this->reffacets; }


/** \brief structure for facets in BuildF2V */
struct DDfacet {
   CollectFacets::mapIterator facetIt;
   double lambda;
   
   enum ddfacetstatus {
     FACETNEW,
     FACETIN,
     FACETON,
     FACETOUT,
     FACETREM
   };
   ddfacetstatus status;

   std::vector<std::forward_list<DDfacet>::iterator> links;
   PosetElem vtx;
 
   DDfacet(CollectFacets::mapIterator facetIt,
		std::vector<std::forward_list<DDfacet>::iterator> &nlinks,
		PosetElem &nvtx) : facetIt(facetIt), status(FACETNEW) {
         this->links.swap(nlinks);
         this->vtx.swap(nvtx);
   }

   bool add_link(std::forward_list<DDfacet>::iterator a, bool check=false) {
      if (check) {
         for (auto b : links) {
            if (a==b) return false;
         } 
      }
      this->links.push_back(a);
      return true;
   }

   void remove_link(const std::forward_list<DDfacet>::iterator &a) {
      for (Index i=0;i<(Index)links.size();i++) {
         if (links[i]==a) { 
	    links[i]=links[links.size()-1]; 
	    links.pop_back();
            return;
         }
      }
   }

   void change_link(const std::forward_list<DDfacet>::iterator &a,
		const std::forward_list<DDfacet>::iterator &b) {
      for (Index i=0;i<(Index)links.size();i++) {
         if (links[i]==a) { links[i]=b; return; }
      }
      /* should not happen */
      links.push_back(b);
   }

   void add_vertex(Index a) {
      this->vtx.add(a);
   }

};

/** \brief construction of an (approximate) set of facets from a
 *  a set of vertices */
class DDbuildV2F {
  public:
     DDbuildV2F(Index idVertex, const Vector &vertex);
   
     int add_vertex(Index idVertex, const Vector& vertex);
     int add_vertex(Index idVertex, const IntervalVector& vertex);

     std::shared_ptr<CollectFacets> getFacets();

     friend std::ostream& operator<<(std::ostream& os, const DDbuildV2F& build);


  private:
     /* can NOT return facets.end() */
     std::forward_list<DDfacet>::iterator addDDfacet(
		const CollectFacets::mapIterator &it,
		std::vector<std::forward_list<DDfacet>::iterator> &links,
		PosetElem &vtx);
     /* can return facets.end() if insertion failed */
     std::forward_list<DDfacet>::iterator addDDfacet(
		std::pair<Row,double> &&row_rhs,
		std::vector<std::forward_list<DDfacet>::iterator> &links,
		PosetElem &vtx);
     /* can return facets.end() if insertion failed */
     std::forward_list<DDfacet>::iterator addDDfacet(
 		const std::pair<Row,double> &row_rhs,
		std::vector<std::forward_list<DDfacet>::iterator> &links,
		PosetElem &vtx);
     /* can return facets.end() if insertion failed */
     std::forward_list<DDfacet>::iterator
	 addFacetSon(const std::forward_list<DDfacet>::iterator &It1, 
		std::forward_list<DDfacet>::iterator &It2, Index idVertex); 
     void removeFacet(std::forward_list<DDfacet>::iterator It,
			 bool removeLnks);

     inline std::pair<Row,double> reduce_facet(const Row &row, double rhs) {
        std::pair<Row,double> ret(row,rhs);
        double mult = reduce_row_without_rounding(ret.first);
        ret.second *= mult;
        return ret;
     }
     
     std::shared_ptr<CollectFacets> facetsptr;
     std::forward_list<DDfacet> facets;
     PosetElem idVertices; /* list of vertices entered, useful as
                long as there are eqfacets */
     Index nbIn=0;
     double tolerance =1e-8;
     
};

inline std::shared_ptr<CollectFacets> DDbuildV2F::getFacets() { return facetsptr; }

}

