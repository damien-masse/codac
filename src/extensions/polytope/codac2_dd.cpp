/**
 * \file codac2_dd.h implementation of the dd algorithm
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
#include <list>

#include "codac2_Index.h"
#include "codac2_Vector.h"
#include "codac2_Row.h"
#include "codac2_Facet.h" 
#include "codac2_dd.h"
#include "codac2_IntvFulPivLU.h"

#define DEBUG_CODAC2_DD

namespace codac2 {

/*************************/
/****** PairMaxSet *******/
/*************************/
template <typename T>
struct PairMaxSets<T> {
   T a,b;
   std::vector<Index> elems;

   PairMaxSets() {}

   PairMaxSets(T a, T b, const std::vector<Index> e1, const std::vector<Index> e2) : a(a), b(b) {
       Index i=0,j=0;
       while(i<(Index)e1.size() && j<(Index)e2.size()) {
         if (e1[i]==e2[j]) { 
           elems.push_back(e1[i]);
           i++; j++;
         } else if (e1[i]<e2[j]) i++;
         else j++;
       }
   }

   void swap(PairMaxSets &p2) {
      std::swap(this->a,p2.a);
      std::swap(this->b,p2.b);
      this->elems.swap(p2.elems);
   }

   int comp(const PairMaxSets &pm) const { 
		/* -1 if elems included in pm.elems; 0 if ==,
		   +1 elems includes pm.elems, +2 if not comparable */
       const std::vector<Index> &e1=elems;
       const std::vector<Index> &e2=pm.elems;
       std::vector<long int>::size_type i=0,j=0;
       bool leq=(e1.size()<=e2.size()), geq=(e1.size()>=e2.size()); 
			/* leq=F,geq=T => +1
                           leq=F,geq=F => 2
                           leq=T,geq=F => -1
                           leq=T,geq=T => 0 */
       while(i<e1.size() && j<e2.size()) {
         if (e1[i]<e2[j]) { 
            if (!geq) return 2;
            leq=false; i++;
         } else if (e1[i]>e2[j]) {
            if (!leq) return 2;
            geq=false; j++;
         } else { i++; j++; }
       }
       if (i<e1.size()) {
          if (!geq) return 2;
          leq=false;
       } 
       if (j<e2.size() || !geq) {
          if (!leq) return 2;
          return -1;
       } 
       return (leq ? 0 : 1);
   }
};


/********************/
/*** DDbuildF2V *****/
/********************/

DDbuildF2V::DDbuildF2V(Index dim, const IntervalVector &box,
	const std::vector<Facet> &eqconstraints) : dim(dim),
		fdim(dim), 
		bbox(box), empty(false), nbInVertex(0), nbInLinks(0) {
   if (eqconstraints.size()!=0) {
      /* managing the equalities */
      /* we first check the existence of a solution, and get a finite
         number of idependant equalities */
      Matrix MatEQ = Matrix(eqconstraints.size(),dim);
      IntervalVector RhsEQ = IntervalVector(eqconstraints.size());
      for (Index i=0;i<eqconstraints.size();i++) {
          MatEQ.row(i) = eqconstraints[i].row;
          RhsEQ[i] = eqconstraints[i].rhs;
      }
      IntvFullPivLU LUdec(MatEQ);
      if (LUdec.is_injective()!=BoolInterval::TRUE) {
         /* TODO : check emptiness */
         /* IntervalMatrix eqsolve = LUdec.solve(RhsEQ);
         if (eqsolve.is_empty()) { empty=true; return; }  */
			/* FIXME : NOT CORRECT, as solve DOES NOT
			   guarantee the emptiness of solutions. We need
			   to use e.g. a bounding box */
      }
      
      

   }
   for (const Facet &f : eqconstraints) {
       assert(f.eqcst);
       Row cst = f.row;
       Interval rhs(f.rhs);
       for (const EqFacet &e : eqfacets) {
           e.adapt_eqconstraint(bbox,cst,rhs);
       }
       try {
           eqfacets.push_back(EqFacet::build_EqFacet(bbox,cst,rhs));
       } catch (const std::domain_error &error) {
           if (rhs.contains(0.0)) continue;
           empty=true;
           eqfacets.clear();
           return;
       }
   } 
   for (const EqFacet &e : eqfacets) {
       e.adapt_box(bbox);
       if (bbox.is_empty()) {
           empty=true;
           eqfacets.clear();
           return;
       }
   }
   /* we order dim_facets in decreasing order, to build fcts in the
      correct order later */
   for (Index i=bbox.size()-1; i>=0; i--) {
       bool found=false;
       for (const EqFacet &e : eqfacets) {
           if (e.dim ==i) { found=true; break; }
       }
       if (!found) dim_facets.push_back(i);
   }
   /* we use "negative" number for the box constraints :
      one for each negative bounds (-s to -1), one for all positive bounds 
	(-s-1).
      idem for links, which means we use specific id for the vertices */
   std::vector<Index> fcts;
   for (Index i : dim_facets) {
      fcts.push_back(-i-1);
   }
   IntervalVector l = bbox.lb();
   this->initDDvertex(-bbox.size()-1,l,fcts);
   Interval sum = bbox.sum();
   for (Index i : dim_facets) {
      l=box.lb();
      l[i]=l[i]+sum.ub()-sum.lb();
      fcts.push_back(-bbox.size()-1);
      for (Index j : dim_facets) {
         if (j==i) continue;
         fcts.push_back(-j-1);
      }
      this->initDDvertex(-i-1,l,fcts);
   }
   for (Index i : dim_facets) {
      for (Index j : dim_facets) {
         if (j==i) continue;
         fcts.push(back(-j-1);
      }
      this->initDDlinks(vertices.find(-bbox.size()-1),
			vertices.find(-i-1),fcts));
      for (Index j : dim_facets) {
         if (j<=i) continue;
         fcts.push_back(-bbox.size()-1);
         for (Index k : dim_facets) {
           if (k!=i) continue;
           if (k!=j) continue;
           fcts.push_back(-k-1);
         }
         this->initDDlinks(vertices.find(-i-1),
			   vertices.find(-j-1),fcts);
      }
   }  
}

Index DDbuildF2V::initDDvertex(Index id, const IntervalVector &v,
		 std::vector<Index> &fcts) {
    vertices.emplace(id,
         DDvertex(v,id,fcts));
    return id;
}
Index DDbuildF2V::addDDvertex(const IntervalVector &v, std::vector<Index> &fcts) {
    nbIn++;
    vertices.emplace(nbIn, DDvertex(v,nbIn,fcts));
    return nbIn;
}
Index DDbuildF2V::addDDvertex(IntervalVector &&v, std::vector<Index> &fcts) {
    nbIn++;
    vertices.emplace(nbIn, DDvertex(v,nbIn,fcts));
    return nbIn;
}

/* p1 has lambda>rhs, p2 has lambda<rhs ;
   p1 is outside the new constraint, p2 is inside.
   For inequalities, we need to get "as much" p1 as possible, 
   that is maximize |lambda2| and minimize |lambda1| */
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



int DDbuildF2V::add_facet(Index idFacet, const Facet& fct) {
   if (empty) return 0;
   IntervalRow rowI = fct.row;
   double rhs = fct.rhs;
   /* treat eqfacets to remove some dimensions */
   for (const EqFacet &e : eqfacets) {
      e.adapt_ineqconstraint(rowI,rhs);
   }
   ineqfacets.emplace(idFacet,std::pair<IntervalRow,double>(rowI,rhs));
   int cnt=0;
   /* compute the status of each vertices */
   for (std::pair<const Index, DDvertex> &v : vertices) {
      DDvertex &d = v.second;
      d.lambda = rowI*d.vertex;
      if (d.lambda.ub()<rhs) { d.status=DDvertex::VERTEXLT; }
      else if (d.lambda.lb()<rhs) { d.status=DDvertex::VERTEXON; }
      else if (d.lambda.lb()<=rhs+tolerance) { d.status=DDvertex::VERTEXGE; }
      else { d.status=DDvertex::VERTEXGT;  }
   }
   for (std::map<Index,DDvertex>::iterator it = vertices.rbegin();
                  it!=vertices.rend();++it) {
     DDvertex &d = it->second;
     std::cout << "     vertex: " << d.Id << " " << d.status << "\n";
     if (d.status==DDvertex::VERTEXLT || 
	 d.status==DDvertex::VERTEXON ||
         d.status==DDvertex::VERTEXNEW)
	 continue;
     /* check the links towards LT and ON, if the facet is a valid
        intersection */
     if (d.status==DDvertex::VERTEXGT) {
        bool oktoremove=true;
        for (auto &l : d.links) {
            if (l==total_links.end()) continue;
            auto &u = l->second.otherVertex(it);
            if (u==vertices.end()) {
               total_links.remove(l);
               l=total_links.end();
               continue;
            }
            if (u->second.status!=VERTEXON 
		&& u->second.status!=VERTEXLT) continue;
        }
     }
     if (d.status==DDvertex::VERTEXGE) {
        d.addFct(idFacet);
        continue;
     }
     for (auto &l : d.links) {
        if (l==total_links.end()) continue;
        auto &u = l->second.otherVertex(it);
        auto &d2 = u->second;
        if (d2.status==DDvertex::VERTEXGE || d2.status==DDvertex::VERTEXGT) {
           l->second.remove(it);
           continue;
        }
        if (d2.status!=DDvertex::VERTEXLT && d2.status!=DDvertex::VERTEXON) { 
	   continue; /* do nothing */
        }
#ifdef DEBUG_CODAC2_DD
        std::cout << "    addVertex " << d.Id << " " << d2.Id << " " << d.status << " " << d2.status << std::endl;
        std::cout << "    (" << d.lambda << " " << d2.lambda << ")    " << std::endl;
#endif
#ifdef DEBUG_CODAC2_DD
        Index newId = this->addVertexSon(d,d2,rhs,idFacet);
        std::cout << " ... done " << newId <<std::endl;
#else
        this->addVertexSon(d,d2,idFacet);
#endif
        cnt++;
     }
     this->removeVertex(d.Id,false);
   }
   /*construction of new links for the eq vertices */
#ifdef DEBUG_CODAC2_DD
        std::cout << "   (" << idFacet << ")" << "    construction of PairMaxSet  " << cnt << std::endl;
#endif
   std::vector<PairMaxSets> maxset;
   for (std::map<Index,DDvertex>::iterator it = vertices.begin();
                  it!=vertices.end();++it) {
      DDvertex &d = it->second;
      if (d.status==DDvertex::VERTEXGE || d.status==DDvertex::VERTEXNEW) {
         for (std::map<Index,DDvertex>::iterator it2 = next(it);
                it2!=vertices.end(); ++it2) {
            DDvertex &d2 = it2->second;
            if (d2.status==DDvertex::VERTEXGE 
		|| d2.status==DDvertex::VERTEXNEW) {
               PairMaxSets actuel(d.Id,d2.Id,d.fcts,d2.fcts);
//               if (actuel.elems.size()<fct.row.size()-2) continue;
               bool placed=false;
               unsigned long int k=0, l=0;
               while (k<maxset.size()) {
                  int u=actuel.comp(maxset[k]);
                  if (u==2) { 
		    if (k!=l) maxset[k].swap(maxset[l]);
	    	    k++; l++;
		    continue;
		  }
                  if (u==-1 || u==0) { placed=true; break;  }
			/* note : if u==-1 or 0, 
			   placed should be false so k==l */
                  if (u==1) {
                     if (!placed) { 
			maxset[k]=actuel; 
			placed=true;
			l++;
                     }
                     k++;
                  } 
               }
               if (l<k) maxset.resize(l);
	       else if (!placed) {
                   maxset.push_back(actuel);
               }
            }
         }
      }
   }
   for (const PairMaxSets &pm : maxset) {
#ifdef DEBUG_CODAC2_DD
      std::cout << " maxset : " << pm.a << " " << pm.b << " " << pm.elems.size() <<  std::endl;
#endif
      vertices.at(pm.a).addLnk(pm.b);
      vertices.at(pm.b).addLnk(pm.a);
   }
   return cnt;
}

/* TODO : make a specific computation */
bool DDbuildF2V::vertex_compatible(const std::vector<Index> &fcts,
                const IntervalVector &v) {
   IntervalMatrix M = IntervalMatrix::zero(fcts.size(), v.size());
   for (unsigned long int i=0;i<fcts.size();i++) {
       M.row(i) = fcts;
   }
   M.row(fcts.size())=v;
   IntvFullPivLU decLU(M);
   return decLU.rank.lb()+eqfacets.size()==v.size();
}

void DDbuildF2V::removeVertex(Index Id, bool removeLnks) {
   DDvertex &d = vertices.at(Id);
   if (removeLnks) {
      for (Index l : d.links) {
         auto it = vertices.find(l);
         if (it==vertices.end()) continue;
         it->second.removeLnk(Id);
      }
   }
   this->vertices.erase(Id);
}

IntervalVector DDbuildF2V::recompute_vertex(const IntervalVector &vect) const {
   IntervalVector v=vect;
   for (Index i = eqfacets.size()-1; i>=0; i--) {
       const EqFacet &ef = eqfacets[i];
       v[ef.dim]=ef.eqcst.dot(v)+ef.valcst;
   }
   return v;
}

std::ostream& operator<<(std::ostream& os, const DDbuildF2V& build) {
   os << "DDbuildF2V (" << build.vertices.size() << " vtx) : " << std::endl;
   for (const auto &ivt : build.vertices) {
      auto &vt = ivt.second;
      os << vt.Id << " : " << vt.vertex << std::endl;
      os << " fct: (";
      for (Index a : vt.fcts) {
         os << a << ",";
      }
      os << ")" << std::endl;
      os << " lnks: (";
      for (Index a : vt.links) {
         os << a << ",";
      }
      os << ")" << std::endl;
   }
   os << "endDDbuildF2V" << std::endl;
   return os;
}

/*****************************/
/*** DDbuildF2V::EqFacet *****/
/*****************************/
/* management of the equality constraints for DDbuildF2V */
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


/*****************************/
/*** DDbuildV2F::EqFacet *****/
/*****************************/

DDbuildV2F::DDbuildV2F(Index idVertex, const Vector &vertex) {
   Row r = Row::Zero(vertex.size());
   for (Index i=0;i<vertex.size();i++) {
      r[i]=1.0;
      eqfacets.emplace_back(r,vertex[i],true);
      r[i]=0.0;
   }
   idVertices.push_back(idVertex);
}

Index DDbuildV2F::addFacetSon(const DDfacet &p1, DDfacet &p2,
        Index idVertex) {
#ifdef DEBUG_CODAC2_DD
    std::cout << " p1.facet.row " <<  p1.facet.row << std::endl;
    std::cout << " p2.facet.row " <<  p2.facet.row << std::endl;
#endif
    Row newV = - p2.lambda*p1.facet.row + p1.lambda*p2.facet.row;
    double newRhs = - p2.lambda*p1.facet.rhs + p1.lambda*p2.facet.rhs;
#ifdef DEBUG_CODAC2_DD
    std::cout << " newV " << newV << std::endl;
    std::cout << " newRhs " << newRhs << std::endl;
#endif
    /* computation of vertices */
    std::vector<Index> nvtx;
    Index a1=0; Index a2=0;
    bool idVertexAdded=false;
    while (a1<(Index)p1.vtx.size() && a2<(Index)p2.vtx.size()) {
       if (p1.vtx[a1]<p2.vtx[a2]) {
          a1++; if (a1==(Index)p1.vtx.size()) break;
       } else if (p2.vtx[a2]<p1.vtx[a1]) {
          a2++; if (a2==(Index)p2.vtx.size()) break;
       } else {
         if (p1.vtx[a1]>idVertex && !idVertexAdded) {
            nvtx.push_back(idVertex);
            idVertexAdded=true;
         }
         nvtx.push_back(p1.vtx[a1]);
         a1++; a2++;
       }
    }
    if (!idVertexAdded) nvtx.push_back(idVertex);
    std::set<Index> lnks;
    lnks.insert(p2.Id);
    Index newId=this->addDDfacet(Facet(newV,newRhs,false),lnks,nvtx); 
    p2.replaceLnks(p1.Id,newId);
    return newId;
}

int DDbuildV2F::add_vertex(Index idVertex, const Vector& vertex) {
    return this->add_vertex(idVertex,IntervalVector(vertex));
}

int DDbuildV2F::add_vertex(Index idVertex, const IntervalVector& vertex) {
   int cnt=0;
   /* first check equalities */
   if (!eqfacets.empty()) {
      Index eqViolated=-1;
      Interval delta;
      for (Index i=eqfacets.size()-1;i>=0;i--) {
          Interval calc = eqfacets[i].row.dot(vertex)-eqfacets[i].rhs;
          if (calc.mig()<this->tolerance) continue;
          if (eqViolated==-1) { eqViolated=i; delta=calc; continue; }
          else {
             double correction = -(calc/delta).mid();
             eqfacets[i].row += correction*eqfacets[eqViolated].row;
             eqfacets[i].rhs += correction*eqfacets[eqViolated].rhs;
          }
      }
      if (eqViolated!=-1) {
         if (facets.empty()) {
            this->nbIn=0;
            std::set<Index> lnk = { 2 };
            std::vector<Index> vtx = { idVertex };
            if (delta.lb()>0.0) {
               this->addDDfacet(Facet(eqfacets[eqViolated].row,
			eqfacets[eqViolated].row.dot(vertex).ub(),false),
		     lnk,vtx);
               vtx = idVertices;
               lnk = { 1 };
               this->addDDfacet(Facet(-eqfacets[eqViolated].row,
			-eqfacets[eqViolated].rhs,false),
		     lnk,vtx);
            } else {
               this->addDDfacet(Facet(-eqfacets[eqViolated].row,
			(-eqfacets[eqViolated].row.dot(vertex)).ub(),false),
		     lnk,vtx);
               vtx = idVertices;
               lnk = { 1 };
               this->addDDfacet(Facet(eqfacets[eqViolated].row,
			eqfacets[eqViolated].rhs,false),
		     lnk,vtx);
           }
           cnt=2;
         } else {
           /* add the remaining inequality, linked with all existing facets */
           std::set<Index> lnk;
           for (std::pair<const Index, DDfacet> &f : facets) {
                 f.second.status=DDfacet::FACETOUT; /* to remove FACETNEW */
		 lnk.insert(f.first);
           }
           Index u;
           std::vector<Index> vtx = idVertices;
           if (delta.lb()>0.0) {
              u = this->addDDfacet(Facet(-eqfacets[eqViolated].row,
				-eqfacets[eqViolated].rhs,false),lnk,vtx);
           } else {
              u = this->addDDfacet(Facet(eqfacets[eqViolated].row,
			eqfacets[eqViolated].rhs,false),lnk,vtx);
           }
           /* change all existing facets to take the new vertex into account */
    
           for (std::pair<const Index, DDfacet> &f : facets) {
               if (f.second.status!=DDfacet::FACETOUT) continue;
               Interval q = f.second.facet.rhs-
				f.second.facet.row.dot(vertex);
	       double correction = (q/delta).mid();
               f.second.facet.row += correction*eqfacets[eqViolated].row;
               f.second.facet.rhs += correction*eqfacets[eqViolated].rhs;
               f.second.addLnk(u);
               f.second.addVtx(idVertex);
           }
           cnt=1;
         }
         if (eqViolated!=(Index)eqfacets.size()-1) {
            eqfacets[eqViolated]=eqfacets[eqfacets.size()-1];
         } 
         eqfacets.pop_back();
         if (!eqfacets.empty()) idVertices.push_back(idVertex);
         return cnt;
      } else idVertices.push_back(idVertex); /* and follows */
   }
   /* no equalities removed */
   for (std::pair<const Index, DDfacet> &f : facets) {
       DDfacet &d = f.second;
       Interval calc = d.facet.row.dot(vertex)-d.facet.rhs;
       d.lambda = calc.lb();
       if (calc.mig()<this->tolerance) { d.status=DDfacet::FACETON; }
       else if (calc.lb()>0.0) { d.status=DDfacet::FACETOUT; }
       else { d.status=DDfacet::FACETIN; }
   } 
   Index endF = nbIn;
   for (std::map<Index,DDfacet>::iterator it = facets.begin();
                  it!=facets.end() && it->first<=endF;) {
     DDfacet &d = it->second;
     ++it;
     if (d.status==DDfacet::FACETON) {
        d.addVtx(idVertex);
        continue;
     }
     if (d.status==DDfacet::FACETIN || d.status==DDfacet::FACETNEW) continue;
     if (d.links.empty()) { /* not supposed to happen */
        /* very specific case of a 1-dimension polyhedra */
        d.vtx.clear();
        d.vtx.push_back(idVertex);
        d.facet.rhs = d.facet.row.dot(vertex).ub();
        continue;
     }
     for (Index idL : d.links) {
        auto it2 = facets.find(idL); 
        if (it2==facets.end()) continue;
        DDfacet &d2 = it2->second;
        if (d2.status==DDfacet::FACETON) {
           d2.removeLnk(d.Id);
           continue;
        }
        if (d2.status!=DDfacet::FACETIN) { 
	   continue; /* do nothing for out */
        }
#ifdef DEBUG_CODAC2_DD
        std::cout << "    addFacet " << d.Id << " " << d2.Id << " " << d.status << " " << d2.status << std::endl;
        std::cout << "    (" << d.lambda << " " << d2.lambda << ")    " << std::endl;
#endif
#ifdef DEBUG_CODAC2_DD
        Index newid = this->addFacetSon(d,d2,idVertex);
        std::cout << " ... done " << newid <<std::endl;
#else
        this->addFacetSon(d,d2,idVertex);
#endif
        cnt++;
     }
     this->removeFacet(d.Id,false);
     
   }
   /*construction of new links for the eq facets */
   std::vector<PairMaxSets> maxset;
   for (std::map<Index,DDfacet>::iterator it = facets.begin();
                  it!=facets.end();++it) {
      DDfacet &d = it->second;
      if (d.status==DDfacet::FACETON || d.status==DDfacet::FACETNEW) {
         for (std::map<Index,DDfacet>::iterator it2 = next(it);
		it2!=facets.end(); ++it2) {
            DDfacet &d2 = it2->second;
            if (d2.status==DDfacet::FACETON || d2.status==DDfacet::FACETNEW) {
               PairMaxSets actuel(d.Id,d2.Id,d.vtx,d2.vtx);
               bool placed=false;
               unsigned long int k=0, l=0;
               while (k<maxset.size()) {
                  int u=actuel.comp(maxset[k]);
                  if (u==2) { 
		    if (k!=l) maxset[k].swap(maxset[l]);
	    	    k++; l++;
		    continue;
		  }
                  if (u==-1 || u==0) { placed=true; break; }
			/* note : if u==-1 or 0, 
			   placed should be false so k==l */
                  if (u==1) {
                     if (!placed) { 
			maxset[k]=actuel;
			placed=true;
			l++;
                     } 
                     k++;
                  } 
               }
               if (l<k) maxset.resize(l);
               else if (!placed) {
                   maxset.push_back(actuel);
               }
            }
         }
      }
   }
   for (const PairMaxSets &pm : maxset) {
//      std::cout << " maxset : " << pm.a << " " << pm.b << std::endl;
      facets.at(pm.a).addLnk(pm.b);
      facets.at(pm.b).addLnk(pm.a);
   }
   return cnt;
}

Index DDbuildV2F::addDDfacet(const Facet &f, std::set<Index> &links,
		std::vector<Index> &vtx) {
    
    nbIn++;
    facets.emplace(nbIn, 
	 DDfacet(f,nbIn,links,vtx));
    return nbIn;
}

Index DDbuildV2F::addDDfacet(Facet &&f, std::set<Index> &links,
		std::vector<Index> &vtx) {
    
    nbIn++;
    facets.emplace(nbIn, 
	 DDfacet(f,nbIn,links,vtx));
    return nbIn;
}

void DDbuildV2F::removeFacet(Index Id, bool removeLnks) {
   DDfacet &d = facets.at(Id);
   if (removeLnks) {
      for (Index l : d.links) {
         auto it = facets.find(l);
         if (it==facets.end()) continue;
         it->second.removeLnk(Id);
      }
   }
   this->facets.erase(Id);
}

std::ostream& operator<<(std::ostream& os, const DDbuildV2F& build) {
   os << "DDbuildV2F (" << build.eqfacets.size() << " eq, " <<
		build.facets.size() << " ineq) : " << std::endl;
   for (const Facet &eqf : build.eqfacets) {
      os << eqf.row << " = " << eqf.rhs << std::endl;
   }
   for (const auto &ifacet : build.facets) {
      auto &ieqf = ifacet.second;
      os << ieqf.Id << " : " << ieqf.facet.row << " <= " << ieqf.facet.rhs << std::endl;
      os << " vtx: (";
      for (Index a : ieqf.vtx) {
         os << a << ",";
      }
      os << ")" << std::endl;
      os << " lnks: (";
      for (Index a : ieqf.links) {
         os << a << ",";
      }
      os << ")" << std::endl;
   }
   os << "endDDbuildV2F" << std::endl;
   return os;
}

}

