/**
 * \file codac2_Polytope_dd.cpp implementation of the dd algorithm
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Damien Massé
 *  \copyright  Copyright 2026
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
#include <list>
#include <queue>
#include <utility>
#include <tuple>
#include <unordered_map>

#include "codac2_Index.h"
#include "codac2_Vector.h"
#include "codac2_Row.h"
#include "codac2_Facet.h" 
#include "codac2_Polytope_dd.h"
#include "codac2_IntvFullPivLU.h"
#include "codac2_PosetMax.h"

#define DEBUG_CODAC2_DD

namespace codac2 {


/********************/
/*** DDbuildF2V *****/
/********************/

double DDbuildF2V::tolerance=1e-5; /* check if needed */

DDbuildF2V::DDbuildF2V(Index dim, const IntervalVector &box,
	std::shared_ptr<CollectFacets> facetsptr, bool include_box) : dim(dim),
		fdim(dim), bbox(box), facetsptr(facetsptr),
		empty(false), flat(false),
		nbIn(0) {
   /* counting the number of equalities */
   std::cout << "counting equalities\n";
   int nb_eq_dim= (facetsptr==nullptr ? 0 : facetsptr->nbeqfcts());
   if (include_box) {
      for (int i=0;i<dim;i++) {
         if (bbox[i].is_degenerated()) nb_eq_dim++;
      }
   }
   if (nb_eq_dim>0) {
      /* managing the equalities */
      /* we first check the existence of a solution, and get a finite
         number of idependant equalities */
      Matrix MatEQ = Matrix::Zero(nb_eq_dim,dim);
      IntervalVector RhsEQ = 
			IntervalVector::Zero(nb_eq_dim);
      int base=0; 
      if (include_box) {
        for (int i=0;i<dim;i++) {
           if (!bbox[i].is_degenerated()) continue;
           MatEQ(base,i)=1.0;
           RhsEQ[base]=box[i];
           base++;
        }
      }
      if (facetsptr!=nullptr) {
        for (Index i=0;i<(Index)facetsptr->nbeqfcts();i++) {
            MatEQ.row(base+i) = facetsptr->get_eqFacet(i)->first.row;
            RhsEQ[base+i] = facetsptr->get_eqFacet(i)->second.rhs;
        }
      }
      IntvFullPivLU LUdec(MatEQ);
      if ((LUdec.is_surjective())!=BoolInterval::TRUE) {
          /* first we check if emptiness is guaranteed */
          IntervalMatrix eqsolve = box;
          LUdec.solve(RhsEQ, eqsolve);
          if (eqsolve.is_empty()) { empty=true; return; }  
          if (LUdec.rank().lb()>0) {
            /* otherwise, the best option is probably to generate a smaller
             * subset of equalities, and discard the others */
            Matrix MatEQ2 = LUdec.coimage(MatEQ);
            IntervalVector RhsEQ2 = LUdec.coimage(RhsEQ);
            IntvFullPivLU LUdec2(MatEQ2); /* we can hope LUdec2 is surjective */
            if (LUdec2.is_surjective()==BoolInterval::TRUE) {
              this->fdim = this->dim - MatEQ2.rows();
              this->M_EQ = IntervalMatrix::Zero(this->dim,this->fdim+1);
              this->M_EQ.col(0) = LUdec2.solve(RhsEQ2);
              this->M_EQ.rightCols(this->fdim) = LUdec2.kernel();
              this->flat=true;
            }
          } /* if we cannot have a surjective matrix, we ignore the 
               equalities */
      } else 
      {
         this->fdim = this->dim - nb_eq_dim;
         this->M_EQ = IntervalMatrix::Zero(this->dim,this->fdim+1);
         this->M_EQ.col(0) = LUdec.solve(RhsEQ);
         this->M_EQ.rightCols(this->fdim) = LUdec.kernel();
	 this->flat=true;
      }
#ifdef DEBUG_CODAC2_DD
      std::cout << "M_EQ : " << M_EQ << std::endl;
#endif
   } 
   std::cout << "initial vertice\n";
   /* initial vertice :
      (1 0 0 ... ), + rays */
   this->vdim=0;
   for (int i=0;i<fdim;i++) {
       IntervalVector line = IntervalVector::Zero(this->fdim+1);
       line[i+1]=1.0;
       lines.push_back(line);
   }
   this->addDDvertex(fdim);
   /* initial projection for rays */
   this->projection_rays.push_back(IntervalMatrix::Identity(this->fdim+1,this->fdim+1));
   this->projection_rays[0](0,0)=0.0;
   std::cout << "before add_constraint_box\n";

   if (include_box) this->add_constraint_box(bbox);
}

int DDbuildF2V::add_bound_var(Index c, bool mx, double rhs) {
   if (empty) return -1;
   if (mx && rhs==+oo) return 0;
   if (!mx && rhs==-oo) return 0;
   IntervalRow facet = IntervalRow::Zero(this->fdim+1);
   facet[0] = -rhs;
   if (flat) {
      facet += this->M_EQ.row(c);
   } else {
      facet[c+1] = 1.0;
   }
   if (!mx) facet=-facet;
   Index idFacet = -c-1-(mx ? 0 : this->dim);
   return this->add_facet_internal(idFacet,facet);
}

int DDbuildF2V::add_constraint_box(const IntervalVector &box) {
   int ret=0;
   for (Index i=0;i<box.size();i++) {
      if (box[i].is_degenerated()) continue;
      int ret1=this->add_bound_var(i,true,box[i].ub());
      if (ret1==-1) return -1;
      if (ret1==1) ret=1;
   }
   for (Index i=0;i<box.size();i++) {
      if (box[i].is_degenerated()) continue;
      int ret1=this->add_bound_var(i,false,box[i].lb());
      if (ret1==-1) return -1;
      if (ret1==1) ret=1;
   }
   return ret;
}

DDvertex* DDbuildF2V::addDDvertex(Index dim) {
    nbIn++;
    vertices.push_front(std::make_shared<DDvertex>(nbIn,dim));
    return vertices.front().get();
}
DDvertex* DDbuildF2V::addDDvertex(const IntervalVector &v, Index vdim) {
    nbIn++;
    IntervalMatrix *proj=nullptr;
    if (vdim>0) {
       assert(projection_rays.size()>vdim);
       proj=&(projection_rays[vdim]);
    }
    vertices.push_front(std::make_shared<DDvertex>(nbIn,v,vdim,proj));
    return vertices.front().get();
}

inline std::vector<Index> fctSon(Index fct,const std::vector<Index> fcts1,
		const std::vector<Index> fcts2) {
    std::vector<Index> nfcts;
    Index a1=0; Index a2=0;
    bool idFacetAdded=(fct<0);
    while (a1<(Index)fcts1.size() && a2<(Index)fcts2.size()) {
       if (fcts1[a1]<fcts2[a2]) {
          a1++; if (a1==(Index)fcts1.size()) break;
       } else if (fcts2[a2]<fcts1[a1]) {
          a2++; if (a2==(Index)fcts2.size()) break;
       } else {
         if (fcts1[a1]>fct && !idFacetAdded) {
            nfcts.push_back(fct);
            idFacetAdded=true;
         }
         nfcts.push_back(fcts1[a1]);
         a1++; a2++;
       }
    }
    if (!idFacetAdded) nfcts.push_back(fct);

    return nfcts;
}

/* p1 has lambda>rhs, p2 has lambda<rhs ;
   p1 is outside the new constraint, p2 is inside.
   For inequalities, we need to get "as much" p1 as possible, 
   that is maximize |lambda2| and minimize |lambda1| */
DDvertex* DDbuildF2V::addVertexSon(const DDvertex &p1,DDvertex &p2,
			const IntervalRow &facet,Index idFacet) {
    nbIn++;
    std::cout << "add vertex son : " << nbIn << "\n";
    std::cout << p1.vertex << " " << p2.vertex << " " << p1.lambda << " " << p2.lambda << "\n";
    std::shared_ptr<DDvertex> ret = 
			std::make_shared<DDvertex>(nbIn,p1,p2,facet);
    vertices.push_front(ret);
    std::cout << "result : " << ret->vertex << "\n";
    ret->facets.add(idFacet);
    return ret.get();
}

#if 0
Index DDbuildF2V::addVertexSon(const DDvertex &p1, DDvertex &p2, double rhs,
        Index idFacet) {
    Interval lambdaP1=p1.lambda.lb();
    Interval lambdaP2=p2.lambda.lb();
    IntervalVector newV = ((lambdaP1-rhs)*p2.vertex
          + (rhs-lambdaP2)*p1.vertex)/(lambdaP1-lambdaP2);
#if 0
    Interval ratio2=1.0/(1.0-p1.lambda/p2.lambda); /* maximize ratio */
    ratio2 &= Interval(0.0,1.0);
    IntervalVector newV = ratio2.ub()*p1.vertex + 
		(Interval(1.0)-ratio2.ub())*p2.vertex;
#endif
    /* computation of facets */
    std::vector<Index> nfcts;
    Index a1=0; Index a2=0;
    bool idFacetAdded=false;
    while (a1<(Index)p1.fcts.size() && a2<(Index)p2.fcts.size()) {
       if (p1.fcts[a1]<p2.fcts[a2]) {
          a1++; if (a1==(Index)p1.fcts.size()) break;
       } else if (p2.fcts[a2]<p1.fcts[a1]) {
          a2++; if (a2==(Index)p2.fcts.size()) break;
       } else {
         if (p1.fcts[a1]>idFacet && !idFacetAdded) {
            nfcts.push_back(idFacet);
            idFacetAdded=true;
         }
         nfcts.push_back(p1.fcts[a1]);
         a1++; a2++;
       }
    }
    if (!idFacetAdded) nfcts.push_back(idFacet);
    std::set<Index> lnks;
    lnks.insert(p2.Id);
    Index newId=this->addDDvertex(newV,lnks,nfcts);
    p2.replaceLnks(p1.Id,newId);
    return newId;
}
#endif

IntervalVector DDbuildF2V::compute_vertex(const IntervalVector &vect) const {
   /* FIXME : does not handle "hybrid" vertices correctly */
   if (!flat) {
     if (vect[0].lb()<=0.0) return vect.tail(this->fdim);
     return vect.tail(this->fdim)/vect[0];
   }
   if (vect[0].lb()<=0.0) return this->M_EQ*vect;
   return this->M_EQ*(vect/vect[0]);
}


#if 0  /* useful ? */
struct CmpDD {
   bool operator()(const std::forward_list<DDvertex>::iterator &a,
		const std::forward_list<DDvertex>::iterator &b) {
       return (b->vertex[0].lb()*a->lambda.lb()<=a->vertex[0].lb()*b->lambda.lb());
   }
};
#endif

bool DDbuildF2V::add_link(DDvertex *V1,
                      DDvertex *V2) {
      return (V1->add_link(V2) && V2->add_link(V1));
}

 /* add facet by removing one line , hence increasing the dimension
    of the vertices/rays set */
int DDbuildF2V::add_facet_new_dimension(Index idFacet, const IntervalRow &facet,
                       Index bestline) {
     std::cout << "add_facet_new_dimension : " << bestline << "\n";
     Index codeFacet = reffacets.size();
     reffacets.push_back(idFacet);
     /* swap the line so the the removed line is the last one */
     if (bestline!=(Index)lines.size()-1) {
       lines[bestline].swap(lines.back());
     }
     IntervalVector& bstline = lines.back();
     Interval ubest = facet.dot(bstline);
     /* each line l is modified using the "best" line L 
      * l' = alpha l + beta L
      * s.t facet.l' = 0
      *  (ie alpha = ubest and beta = -facet.l 
      */
     for (Index i=0;i<(Index)lines.size()-1;i++) {
        Interval u = facet.dot(lines[i]);
        lines[i] = ubest*lines[i] - u*bstline;
        reduce_vector_without_rounding(lines[i]);
     }
     std::cout << "after lines modified\n";
     /* vertices management.
        + change each existing vertex to satisfy the constraint and
          add the (positive) constraint
        + create a new direction of level+1, add all the previously existing
          (positive) constraints */
     for (std::shared_ptr<DDvertex> &v : vertices) {
         v->add_facet(codeFacet);
         std::cout << "add constraint line " << v->vertex << "\n";
         v->add_constraint_line(bstline, facet, ubest);
         std::cout << "after add constraint line " << v->vertex << "\n";
     }
     std::forward_list<std::shared_ptr<DDvertex>>::iterator itVertex = vertices.begin();
     /* create the new ray, with the new dimension */
     this->vdim++;
     this->add_zero(facet); /* FIXME : what's the meaning of this? */
     std::cout << "after add_zero\n";
     std::cout << "addDDvertex : " << bstline << " " << vdim << "\n";
     DDvertex* newRay = this->addDDvertex(ubest.lb()<0.0 ? bstline : (-bstline),this->vdim);
     std::cout << "new vertex : " << newRay->Id << " " << newRay->vertex << "\n";
     /* add the zero constraint for the new ray */
     for (;itVertex!=vertices.end();++itVertex) {
         this->add_link(itVertex->get(),newRay);
     }
     newRay->add_facet(-1); /* FIXME : not sure */
     for (Index a=0;a<codeFacet;a++) {
		 /* add all facets except the last one */
         newRay->add_facet(a);
     }
     this->lines.pop_back();
     return 1;
}

void DDbuildF2V::add_zero(const IntervalRow &facet) {
     const IntervalMatrix &M = projection_rays.back();
     std::cout << "add_zero : " << M << "\n" << facet << "\n";
     Interval denom = facet * (M*facet.transpose());
     projection_rays.push_back
	(M * (IntervalMatrix::Identity(fdim+1,fdim+1)
		- (1.0/denom) * (facet.transpose()*facet*M)));
     std::cout << "result : " << projection_rays.back() << "\n";
}


bool DDvertex::try_extend(const IntervalRow &facet) {
    auto ret = this->compute_closest(facet);
    if (std::get<0>(ret) ==-1) return false;
    this->vertex=std::get<1>(ret);
    this->cone=std::get<2>(ret);
    return true;
}

bool DDbuildF2V::find_path(DDvertex *v1, DDvertex* v2,
		const PosetElem &fcts) {
    std::unordered_set<DDvertex *> explored;
    std::vector<DDvertex *> stack;
    stack.push_back(v1);
    explored.insert(v1);
    while (!stack.empty()) {
        DDvertex *v = stack.back();
        stack.pop_back();
        if (v->links.find(v2)!=v->links.end()) return true;
        for (auto l : v->links) {
           if (explored.contains(l)) continue;
           if (fcts.comp(l->facets)>0) continue;
           stack.push_back(l);
           explored.insert(l);
        }
    }
    return false;
}

/* add_facet with id and facet IntervalRow : 
   equalities have already been taken into account , emptiness also ,
   return 1: changed    0 : not changed   -1 : empty */
int DDbuildF2V::add_facet_internal(Index idFacet, const IntervalRow &facet) {
   std::cout << "add_facet_internal" << facet << "\n";
   /*  consider lines */
   Index bestline=-1;
   double maxmig=0.0;
   for (Index i=0;i<(Index)lines.size();i++) {
      Interval u = facet.dot(lines[i]);
      if (u.contains(0.0)) continue;
      double mg=u.mig();
      if (mg>maxmig) {
         bestline=i;
         maxmig=mg;
      }
   }  
   if (bestline!=-1) {
      return add_facet_new_dimension(idFacet,facet,bestline);
   }
   /* end of line treatments , the following 
      is when no line is available */
   int nbGT=0;
   int nbGE=0;
   int nbTot=0;
   bool notempty=false;
   std::vector<std::shared_ptr<DDvertex>> onvector;
   std::vector<std::shared_ptr<DDvertex>> gtvector;
   for (auto &vt_ptr : vertices) {
      DDvertex &vt = (*vt_ptr);
      vt.lambda = facet.dot(vt.vertex);
      notempty |= (vt.lambda.lb()<=0.0);
      nbTot++;
      if (vt.lambda.ub()< -this->tolerance) {
         vt.status = DDvertex::VERTEXLT;
         continue;
      }
      if (vt.lambda.lb()>=0.0) {
         nbGE++;
         if (vt.lambda.lb()>this->tolerance) {
            nbGT++;
            vt.status=DDvertex::VERTEXGT;
	    gtvector.push_back(vt_ptr);
            continue;
         }
         vt.status=DDvertex::VERTEXGE;
         continue;
      }
      onvector.push_back(vt_ptr);
   }
   if (!notempty) { this->empty=true; vertices.clear(); return -1; }
   if (nbGT==0) return 0;
   for (auto &vt_ptr : onvector) {
      DDvertex &vt = (*vt_ptr);
      if (vt.try_extend(facet)) { 
         nbGE++;
         vt.status = DDvertex::VERTEXGE;
         continue;
      }
      vt.status = DDvertex::VERTEXON;
   }

   std::vector<DDvertex *> newlist;

   Index codeFacet = reffacets.size();
   reffacets.push_back(idFacet);
   for (std::shared_ptr<DDvertex> &actVertex : gtvector) {
#ifdef DEBUG_CODAC2_DD
       std::cout << " elimination sommet : " << actVertex->Id << ", " << actVertex->links.size() << "liens \n";
#endif
       actVertex->status = DDvertex::VERTEXREM;
       for (DDvertex *lnk_ptr : actVertex->links) {
	   std::cout << "lien second : " << lnk_ptr->Id << " " << lnk_ptr->vertex << "\n";
           DDvertex &destlink = *lnk_ptr;
           if (destlink.status==DDvertex::VERTEXGT
			|| destlink.status==DDvertex::VERTEXREM) continue;
           if (destlink.status==DDvertex::VERTEXGE) {
		if (destlink.add_facet(codeFacet)) newlist.push_back(lnk_ptr);
		/* note: ce serait pas mieux de ne pas ajouter tout de suite
			codeFacet ? */
		destlink.remove_link(actVertex.get());
                continue;
           }  
#ifdef DEBUG_CODAC2_DD
           std::cout << " lien : " << destlink.Id << "\n";
#endif
           assert_release(destlink.status!=DDvertex::VERTEXNEW);
           /* in the following we are supposing that we can't extend
              the vertex */
	   /* creation of a new vertex. Note: even with a ON, we consider
              the full vertex. The goal is the create a new vertex+cone */
	   DDvertex* newVx =
		this->addVertexSon(*actVertex,destlink,facet,codeFacet);
#if 0
           IntervalVector partV=destlink.vertex;
	   IntervalVector newV = -destlink.lambda.lb()*actVertex->vertex
		      + actVertex->lambda.lb()*partV;
           assert(!newV.is_empty());
           std::forward_list<DDvertex>::iterator newVx = 
		     this->addDDvertex(newV);
#endif
#ifdef DEBUG_CODAC2_DD
//       std::cout << " creation new V " << (*newVx) << "id : " << newVx->Id << " coef " << (-destlink.lambda.lb()/actVertex->lambda.lb()) << "\n";
#endif
	   destlink.change_link(actVertex.get(),newVx);
	   newVx->add_link(lnk_ptr);
           newlist.push_back(newVx);
       }
   }
   /* travail sur les liens */
   PosetElemMap<std::pair<DDvertex*, DDvertex*>> mapFcts;
   for (auto it1=newlist.begin();it1!=newlist.end();++it1) {
       for (auto it2 = (it1+1);it2!=newlist.end();++it2) {
	  PosetElem p((*it1)->get_facets(),(*it2)->get_facets());
          if (p.nb_elems<this->vdim-1) continue;
          mapFcts.insert(std::make_pair(
		std::move(p),
		std::make_pair(*it1,*it2)));
       }
   }
   for (auto mapIt = mapFcts.begin();
	mapIt!=mapFcts.end();++mapIt) {
        if (!DDbuildF2V::find_path(mapIt->second.first,
			mapIt->second.second,mapIt->first)) {
	    /* note : check also if there is a direct link */
	    DDbuildF2V::add_link(mapIt->second.first,
			mapIt->second.second);
        }
   }
   vertices.remove_if([](const std::shared_ptr<DDvertex> &v) 
			{ return v->status==DDvertex::VERTEXREM; });
   return 1;
}


int DDbuildF2V::add_facet(Index id) {
   return this->add_facet((*facetsptr)[id]);
}


int DDbuildF2V::add_facet(CollectFacets::mapCIterator itFacet) {
   Index idFacet = itFacet->second.Id;
   if (empty) return 0;
   if (itFacet->second.eqcst)
	throw std::domain_error("add_facet with an equality facet.");
   /* the first step is to take the equalities constraints into account */
   IntervalRow facet = IntervalRow::Zero(this->fdim+1);
   facet[0] = -itFacet->second.rhs;
   if (flat) {
      facet += itFacet->first.row * this->M_EQ;
   } else {
      facet.tail(this->fdim) = itFacet->first.row;
   }
   return this->add_facet_internal(idFacet,facet);
}

/* detection of "irredundant facets" */
std::set<Index> DDbuildF2V::redundantFacets() {
   std::unordered_map<Index,PosetElem> mapping;
   std::set<Index> resultat; /* set of existing facets */
   for (Index i = -2*dim;i<0;i++) {
      resultat.insert(i);
   }
   if (facetsptr) {
     for (const auto &u : facetsptr->get_map()) {
        resultat.insert(u.second.Id);
     } 
   }
   for (const std::shared_ptr<DDvertex> &vt : vertices) { 
      /* to help adding elements at the end, and since id are in 
         decreasing order, we add -vt->id */
      vt->facets.reverse_relation(-vt->Id,mapping);
   }
   PosetElemMaxSet<Index> facetElems;
   for (std::pair<const Index,PosetElem> &p : mapping) {
       facetElems.addElement(p); /* swap p.second */
   }
   /* to get the result of facetElems, we store booleans for each facets */
   std::vector<bool> redundant(reffacets.size(),true);
   for (const std::pair<Index,PosetElem> &elm : facetElems) {
	if (elm.first<0) continue;
        redundant[elm.first]=false;
        resultat.erase(reffacets[elm.first]); 
		/* note : this approach (fill resultat then remove
		   irredundant facets) is ok when one facet may "appear"
                   several times */
   }
   /* renumbering reffacets */
   std::vector<Index> new_reffacets;
   Index nb_irredundant=0;
   for (unsigned int i=0;i<reffacets.size();i++) {
       if (redundant[i]) { reffacets[i]=-1; continue; }
       new_reffacets.push_back(reffacets[i]);
       reffacets[i]=nb_irredundant;
       nb_irredundant++;
   } 
   /* renumbering the index of vertices */
   for (std::shared_ptr<DDvertex> &vt : vertices) { 
      vt->facets = vt->facets.renumber(reffacets);
   }
   reffacets = std::move(new_reffacets);
   return resultat;
}

/* update reffacets following a renumbering */
void DDbuildF2V::update_renumber(const std::vector<Index> &ret) {
   for (Index &a : reffacets) {
      if (a<=0) continue;
      if (ret[a-1]==-1) {
         a=0;
      } else { 
         a=ret[a-1];
      }
   }
}

/* build the bbox of all vertices */
IntervalVector DDbuildF2V::build_bbox() const {
   IntervalVector b = IntervalVector::empty(dim);
   for (const auto &a : vertices) {
      b |= this->compute_vertex(a->vertex);
   }
   return b;
}


/* normalisation for printing of vertice... good idea ? */
IntervalVector normalise(const IntervalVector &v) {
   IntervalVector res(v);
   if (res[0].lb()<=0.0) return res;
   return res/res[0].lb();
}

std::ostream& operator<<(std::ostream& os, const DDbuildF2V& build) {
   os << "DDbuildF2V (" << /*build.vertices.size() <<*/ " vtx) : " << std::endl;
   os << "   " << build.lines.size() << " lines : " << std::endl;
   for (const auto &l : build.lines) {
     os << " " << l << std::endl;
   }
   Index nbvx=0, nbrays=0;
   for (const auto &vt : build.vertices) {
      if (vt->vertex[0].lb()>0.0) nbvx++; else nbrays++;
      os << " " << vt->Id << " : " << normalise(vt->vertex) << 
		" cone " << vt->cone << std::endl;
      os << " fct : (";
      os << vt->facets; /* TODO : use reffacets ! */
#if 0
      for (const auto &a : vt->facets) {
         os << a << ":";
	 if (a>=0) os << build.reffacets[a];
	 else os << "oo";
	 os << ",";
      }
#endif
      os << ")" << std::endl;
      os << " lnks: (";
      for (const auto &a : vt->links) {
         os << a->Id << ",";
      }
      os << ")" << std::endl;
   }
   os << "endDDbuildF2V ( "<< nbrays << " rays, " << nbvx << " vertices, " << (nbrays+nbvx) << " total)" << std::endl;
   return os;
}

/*****************************/
/*** DDbuildF2V::EqFacet *****/
/*****************************/
/* management of the equality constraints for DDbuildF2V */
#if 0
void DDbuildF2V::EqFacet::adapt_eqconstraint(const IntervalVector &bbox,
                                Row &cst, Interval &rhs) const {
           Interval alpha(cst[dim]); cst[dim]=0.0;
           IntervalRow newcst = cst + alpha*this->eqcst;
           cst = newcst.mid();
           rhs += (newcst-cst).dot(bbox)-alpha*valcst;
}
void DDbuildF2V::EqFacet::adapt_ineqconstraint(IntervalRow &cst, double &rhs)
	 const {
           Interval alpha(cst[dim]); cst[dim]=0.0;
           cst = cst + alpha*this->eqcst;
           rhs = (rhs-alpha*valcst).ub();
}
DDbuildF2V::EqFacet 
	DDbuildF2V::EqFacet::build_EqFacet(const IntervalVector &bbox,
                        const Row &cst, Interval &rhs) {
            /* identification of the dimension */
            Index dim = -1;
            double val=0;
            for (Index i=0;i<cst.size();i++) {
                double v2 = std::abs(cst[i]);
                if (v2>val) {
                   dim=i;
                   val = v2;
                }
            }
            if (dim==-1) throw 
			std::domain_error("zero vector in build_EqFacet");
            IntervalRow newcst(cst);
            val=cst[dim];
            newcst[dim]=0.0;
            newcst = newcst/val;
            rhs = rhs/val - (newcst-newcst.mid()).dot(bbox);
            return EqFacet(dim,-newcst.mid(),rhs);
}

void DDbuildF2V::EqFacet::adapt_box(IntervalVector &bbox) const {
    Interval eval = this->eqcst.dot(bbox) + valcst;
    if ((bbox[this->dim] & eval).is_empty()) {
       bbox.set_empty();
    }
    bbox[this->dim]=0;
}
#endif

/*****************************/
/*** DDbuildV2F::EqFacet *****/
/*****************************/

DDbuildV2F::DDbuildV2F(Index idVertex, const Vector &vertex) :
	facetsptr(std::make_shared<CollectFacets>()) {
   Row r = Row::Zero(vertex.size());
   for (Index i=0;i<vertex.size();i++) {
      r[i]=1.0;
      facetsptr->insert_facet(r,vertex[i],true);
      r[i]=0.0;
   }
   idVertices.add(idVertex);
}

std::forward_list<DDfacet>::iterator DDbuildV2F::addFacetSon(
		const std::forward_list<DDfacet>::iterator &It1, 
                std::forward_list<DDfacet>::iterator &It2, Index idVertex) {
    const CollectFacets::mapIterator &p1 = It1->facetIt;
    const CollectFacets::mapIterator &p2 = It2->facetIt;
#ifdef DEBUG_CODAC2_DD
    std::cout << " p1.facet.row " <<  p1->first.row << std::endl;
    std::cout << " p2.facet.row " <<  p2->first.row << std::endl;
#endif
    Row newV = - It2->lambda*p1->first.row + It1->lambda*p2->first.row;
    double newRhs = - It2->lambda*p1->second.rhs + It1->lambda*p2->second.rhs;
#ifdef DEBUG_CODAC2_DD
    std::cout << " newV " << newV << std::endl;
    std::cout << " newRhs " << newRhs << std::endl;
#endif
    /* computation of vertices */
    PosetElem nvtx(It1->vtx,It2->vtx);
    nvtx.add(idVertex);
    std::vector<std::forward_list<DDfacet>::iterator> lnks;
    lnks.push_back(It2);
    std::forward_list<DDfacet>::iterator
	 newF=this->addDDfacet(reduce_facet(newV,newRhs),lnks,nvtx); 
    if (newF!=facets.end()) It2->change_link(It1,newF);
	else It2->remove_link(It1);
    return newF;
}

int DDbuildV2F::add_vertex(Index idVertex, const Vector& vertex) {
    return this->add_vertex(idVertex,IntervalVector(vertex));
}

int DDbuildV2F::add_vertex(Index idVertex, const IntervalVector& vertex) {
   int cnt=0;
   /* first check equalities */
   if (facetsptr->nbeqfcts()>0) {
      Index eqViolated=-1;
      CollectFacets::mapIterator eqViolatedIt= facetsptr->endFacet();
      Interval delta;
      for (Index i=facetsptr->nbeqfcts()-1;i>=0;i--) {
          CollectFacets::mapIterator fcIt = facetsptr->get_eqFacet(i);
          Interval calc = fcIt->first.row.dot(vertex)-fcIt->second.rhs;
          if (calc.mig()<this->tolerance) continue;
          if (eqViolated==-1) { 
	     eqViolated=i; 
	     delta=calc; 
 	     eqViolatedIt = fcIt;
	     continue;
          } else {
             double correction = -(calc/delta).mid();
 	     bool success;
             std::tie(fcIt,success) = facetsptr->change_eqFacet(i,
		  fcIt->first.row + correction*eqViolatedIt->first.row,
		  fcIt->second.rhs + correction*eqViolatedIt->second.rhs);
	     /* if success is false, the eqFacet index may have changed */
             if (!success && 
		 eqViolated==facetsptr->nbeqfcts()) eqViolated=i; 
          }
      }
      if (eqViolated!=-1) {
         if (facets.empty()) {
            this->nbIn=0;
            std::vector<std::forward_list<DDfacet>::iterator> lnk;
            PosetElem vtx; 
            vtx.add(idVertex);
            Interval nrhsI = eqViolatedIt->first.row.dot(vertex);
            double nrhs = delta.lb()>0.0 ? nrhsI.ub() : nrhsI.lb();
	    CollectFacets::mapIterator newIt =
		facetsptr->dissociate_eqFacet(eqViolated,nrhs);
	    assert_release(newIt!=facetsptr->endFacet());
            if (delta.lb()>0.0) {
               auto a = this->addDDfacet(eqViolatedIt,lnk,vtx);
               vtx = idVertices;
               lnk = { a };
               auto b = this->addDDfacet(newIt,lnk,vtx);
	       a->add_link(b,false);
            } else {
               auto a = this->addDDfacet(newIt,lnk,vtx);
               vtx = idVertices;
               lnk = { a };
               auto b = this->addDDfacet(eqViolatedIt,lnk,vtx);
	       a->add_link(b,false);
           }
           cnt=2;
         } else {
           /* add the remaining inequality, linked with all existing facets */
           std::vector<std::forward_list<DDfacet>::iterator> lnk;
           for (std::forward_list<DDfacet>::iterator it =
			facets.begin(); it!=facets.end(); ++it) {
                 it->status=DDfacet::FACETOUT; /* to remove FACETNEW */
		 lnk.push_back(it);
           } 
           Row saveEqrow = eqViolatedIt->first.row;
           double saveRhs = eqViolatedIt->second.rhs;
	   CollectFacets::mapIterator nIt =
		facetsptr->dissociate_eqFacet(eqViolated,
			(delta.lb()>0.0 ? +oo : -oo));
//	   std::cout << "dissociated : " << (delta.lb()>0.0 ? +oo : -oo) <<
//                  " : " << (nIt == facetsptr->endFacet()) << "\n";
           eqViolatedIt=nIt; /* in the following, the case nIt==facets.end()
			        is not treated, as it should not happen */
           PosetElem vtx = idVertices;
           
           std::forward_list<DDfacet>::iterator u = 
		(eqViolatedIt != facetsptr->endFacet() ? 
			this->addDDfacet(eqViolatedIt,lnk,vtx) :
			facets.end());
           /* change all existing facets to take the new vertex into account */
           /* no facet should be removed, theoretically, but we
	    * nevertheless consider this possibility */
           std::forward_list<DDfacet>::iterator it_f=facets.before_begin();
           std::forward_list<DDfacet>::iterator nxt = facets.begin();
           while (nxt!=facets.end()) {
               DDfacet &f = (*nxt);
               if (f.status!=DDfacet::FACETOUT) { it_f=nxt; ++nxt; continue; }
               Interval q = f.facetIt->second.rhs-
				f.facetIt->first.row.dot(vertex);
	       double correction = (q/delta).mid();
               if (u!=facets.end()) f.add_link(u);
               f.add_vertex(idVertex);
 	       bool success;
               std::tie(f.facetIt, success) =
		   facetsptr->change_ineqFacet(f.facetIt->second.Id,
			    f.facetIt->first.row + correction*saveEqrow,
			    f.facetIt->second.rhs + correction*saveRhs);
	       if (!success) { this->removeFacet(nxt,true); 
			       facets.erase_after(it_f); 
			       nxt=it_f; ++nxt; continue; }
               it_f=nxt; ++nxt;
           }
           cnt=(u!=facets.end() ? 1 : 0);
         }
         if (facetsptr->nbeqfcts()>0) idVertices.add(idVertex);
         return cnt;
      } else idVertices.add(idVertex); /* and follows */
   }
   /* no equalities removed */
   for (DDfacet &d : facets) {
       if (d.status==DDfacet::FACETREM) continue;
       Interval calc = d.facetIt->first.row.dot(vertex)-d.facetIt->second.rhs;
       d.lambda = calc.lb();
       if (calc.mig()<this->tolerance) { d.status=DDfacet::FACETON; }
       else if (calc.lb()>0.0) { d.status=DDfacet::FACETOUT; }
       else { d.status=DDfacet::FACETIN; }
   } 
   for (std::forward_list<DDfacet>::iterator it = facets.begin();
                  it!=facets.end();++it) {
     DDfacet &d = (*it);
     if (d.status==DDfacet::FACETON) {
        d.add_vertex(idVertex);
        continue;
     }
     if (d.status==DDfacet::FACETIN || d.status==DDfacet::FACETNEW) continue;
     if (d.links.empty()) { /* not supposed to happen */
        /* very specific case of a 1-dimension polyhedra */
        d.vtx.clear();
        d.vtx.add(idVertex);
        facetsptr->change_ineqFacet_rhs(d.facetIt, 
			d.facetIt->first.row.dot(vertex).ub());
        continue;
     }
     for (std::forward_list<DDfacet>::iterator idL : d.links) {
        if (idL==facets.end()) continue; /* can it happen ? */
        DDfacet &d2 = (*idL);
        if (d2.status==DDfacet::FACETON) {
           d2.remove_link(it);
           continue;
        }
        if (d2.status!=DDfacet::FACETIN) { 
	   continue; /* do nothing for out */
        }
#ifdef DEBUG_CODAC2_DD
        std::cout << "    addFacet " << d.facetIt->second.Id << " " << d2.facetIt->second.Id << " " << d.status << " " << d2.status << std::endl;
        std::cout << "    (" << d.lambda << " " << d2.lambda << ")    " << std::endl;
#endif
        [[maybe_unused]]  std::forward_list<DDfacet>::iterator 
		newIt = this->addFacetSon(it,idL,idVertex);
#ifdef DEBUG_CODAC2_DD
        if (newIt==facets.end()) 
	  std::cout << " ... FAILED! " << std::endl;
        else 
          std::cout << " ... done " << newIt->facetIt->second.Id << std::endl;
#endif
        cnt++;
     }
     facetsptr->removeFacetById(it->facetIt->second.Id);
     this->removeFacet(it,false);
   }
   /*construction of new links for the eq facets */
   PosetElemMaxSet<std::pair<std::forward_list<DDfacet>::iterator,
			std::forward_list<DDfacet>::iterator>> maxset;
   for (std::forward_list<DDfacet>::iterator it = facets.begin();
                  it!=facets.end();++it) {
      DDfacet &d = (*it);
      if (d.status==DDfacet::FACETON || d.status==DDfacet::FACETNEW) {
         for (std::forward_list<DDfacet>::iterator it2 = next(it);
		it2!=facets.end(); ++it2) {
            DDfacet &d2 = (*it2);
            if (d2.status==DDfacet::FACETON || d2.status==DDfacet::FACETNEW) {
	       maxset.addElement(
			std::make_pair(std::pair(it,it2),d.vtx & d2.vtx));
            }
         }
      }
   }
   for (const auto &pm : maxset) {
//      std::cout << " maxset : " << pm.a << " " << pm.b << std::endl;
      pm.first.first->add_link(pm.first.second);
      pm.first.second->add_link(pm.first.first);
   }
   facets.remove_if([](const DDfacet &v) 
			{ return v.status==DDfacet::FACETREM; });
   facetsptr->renumber();
   return cnt;
}

std::forward_list<DDfacet>::iterator 
	DDbuildV2F::addDDfacet(const CollectFacets::mapIterator &it,
		 std::vector<std::forward_list<DDfacet>::iterator> &links,
		 PosetElem &vtx) {
    facets.emplace_front(it,links,vtx);
    return facets.begin();
}

std::forward_list<DDfacet>::iterator 
	DDbuildV2F::addDDfacet(const std::pair<Row, double> &row_rhs,
		 std::vector<std::forward_list<DDfacet>::iterator> &links,
		 PosetElem &vtx) {
    std::pair<Index,bool> result =
	 facetsptr->insert_facet(row_rhs.first,row_rhs.second,false);
    if (!result.second) {
       return facets.end();
    }
    facets.emplace_front((*facetsptr)[result.first-1],links,vtx);
    return facets.begin();
}

std::forward_list<DDfacet>::iterator 
	DDbuildV2F::addDDfacet(std::pair<Row,double> &&row_rhs,
		 std::vector<std::forward_list<DDfacet>::iterator> &links,
		 PosetElem &vtx) {
    std::pair<Index,bool> result  =
	 facetsptr->insert_facet(row_rhs.first,row_rhs.second,false);
    if (!result.second) {
       return facets.end();
    }
    facets.emplace_front((*facetsptr)[result.first-1],links,vtx);
    return facets.begin();
}

void DDbuildV2F::removeFacet(std::forward_list<DDfacet>::iterator It,
		 bool removeLnks) {
   DDfacet &d = (*It);
   if (removeLnks) {
      for (std::forward_list<DDfacet>::iterator it : d.links) {
         it->remove_link(It);
      }
   }
   d.facetIt=facetsptr->endFacet();
   d.status=DDfacet::FACETREM;
}

std::ostream& operator<<(std::ostream& os, const DDbuildV2F& build) {
   os << "DDbuildV2F (" << build.facetsptr->nbeqfcts() << " eq, " <<
		build.facetsptr->nbfcts() << " total fcts) : " << std::endl;
   for (Index i=0;i<build.facetsptr->nbeqfcts();++i) {
      const Facet &eqf = (*build.facetsptr->get_eqFacet(i));
      os << eqf.first.get_row() << " = " << eqf.second.get_rhs() << std::endl;
   }
   for (const auto &ifacet : build.facets) {
      auto &ieqf = *(ifacet.facetIt);
      os << ieqf.second.get_Id() << " : " << ieqf.first.get_row() << " <= " << ieqf.second.get_rhs() << std::endl;
      os << " vtx: (";
      os << ifacet.vtx;
      os << ")" << std::endl;
      os << " lnks: (";
      for (std::forward_list<DDfacet>::iterator a : ifacet.links) {
         os << a->facetIt->second.get_Id() << ",";
      }
      os << ")" << std::endl;
   }
   os << "endDDbuildV2F" << std::endl;
   return os;
}

}

