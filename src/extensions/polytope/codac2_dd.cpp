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
#include <queue>
#include <utility>

#include "codac2_Index.h"
#include "codac2_Vector.h"
#include "codac2_Row.h"
#include "codac2_Facet.h" 
#include "codac2_dd.h"
#include "codac2_IntvFullPivLU.h"

#define DEBUG_CODAC2_DD

namespace codac2 {

/*************************/
/****** PairMaxSet *******/
/*************************/
template <typename T>
struct PairMaxSets {
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

   PairMaxSets(T a, T b, const std::vector<Index> e1, const std::vector<Index> e2, const Index rejectF) : a(a), b(b) {
       Index i=0,j=0;
       while(i<(Index)e1.size() && j<(Index)e2.size()) {
         if (e1[i]==e2[j]) { 
           if (e1[i]!=rejectF) 
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
	const std::vector<Facet> &eqconstraints, bool include_box) : dim(dim),
		fdim(dim), 
		bbox(box), empty(false), flat(false),
		nbIn(0) {
   int nb_eq_dim= eqconstraints.size();
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
      for (Index i=0;i<(Index)eqconstraints.size();i++) {
          MatEQ.row(base+i) = eqconstraints[i].row;
          RhsEQ[base+i] = eqconstraints[i].rhs;
      }
      IntvFullPivLU LUdec(MatEQ);
#if 0
      if (LUdec.is_surjective)!=BoolInterval::TRUE) {
         /* TODO : check emptiness */
         /* IntervalMatrix eqsolve = LUdec.solve(RhsEQ);
         if (eqsolve.is_empty()) { empty=true; return; }  */
			/* FIXME : NOT CORRECT, as solve DOES NOT
			   guarantee the emptiness of solutions. We need
			   to use e.g. a bounding box... wait for 
			   pull request */
      } else 
#endif
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
   /* initial vertice :
      (1 0 0 ... ), + rays */
   IntervalVector v = IntervalVector::Zero(this->fdim+1);
   v[0]=1.0;
   this->addDDvertex(v);
   for (int i=0;i<fdim;i++) {
       IntervalVector line = IntervalVector::Zero(this->fdim+1);
       line[i+1]=1.0;
       lines.push_back(line);
   }
//   if (include_box) this->add_constraint_box(bbox);
}

std::forward_list<DDvertex>::iterator DDbuildF2V::addDDvertex(const IntervalVector &v) {
    nbIn++;
    vertices.emplace_front(DDvertex(v,nbIn));
    return vertices.begin();
}
std::forward_list<DDvertex>::iterator  DDbuildF2V::addDDvertex(IntervalVector &&v) {
    nbIn++;
    vertices.emplace_front(DDvertex(v,nbIn));
    return vertices.begin();
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

#if 0
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
#endif

IntervalVector DDbuildF2V::compute_vertex(const IntervalVector &vect) const {
   if (!flat) {
     if (vect[0]==0.0) return vect.tail(this->fdim);
     return vect.tail(this->fdim)/vect[0];
   }
   if (vect[0]==0.0) return this->M_EQ*vect;
   return this->M_EQ*(vect/vect[0]);
}


struct CmpDD {
   bool operator()(const std::forward_list<DDvertex>::iterator &a,
		const std::forward_list<DDvertex>::iterator &b) {
       return (b->vertex[0].lb()*a->lambda.lb()<=a->vertex[0].lb()*b->lambda.lb());
   }
};

bool DDbuildF2V::addDDlink(const std::forward_list<DDvertex>::iterator& V1,
                      const std::forward_list<DDvertex>::iterator& V2) {
      if (V1->addLnk(V2,true)) return V2->addLnk(V1,false);
      return false;
//    return (V1->addLnk(V2,true) && V2->addLnk(V1,false));
}

int DDbuildF2V::add_facet(Index idFacet, const Facet& fct) {
   if (empty) return 0;
   /* the first step is to take the equalities constraints into account */
   IntervalRow facet = IntervalRow::Zero(this->fdim+1);
   facet[0] = -fct.rhs;
   if (flat) {
      facet += fct.row * this->M_EQ;
   } else {
      facet.tail(this->fdim) = fct.row;
   }
   /* then consider lines */
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
      Index codeFacet = reffacets.size();
      reffacets.push_back(idFacet);
      if (bestline!=(Index)lines.size()-1) {
	lines[bestline].swap(lines.back());
      }
      IntervalVector& bstline = lines.back();
      Interval ubest = facet.dot(bstline);
      /* traitement des lines */
      for (Index i=0;i<(Index)lines.size()-1;i++) {
         Interval u = facet.dot(lines[i]);
         lines[i] = ubest*lines[i] - u*bstline;
      }
      /* traitement des sommets */
      for (DDvertex &v : vertices) {
         double lower =  facet.dot(v.vertex).lb();
         if (ubest.lb()>0.0) {
            if (lower>0.0) 
              v.vertex = ubest.ub()*v.vertex - lower*bstline;
            else
	      v.vertex = ubest.lb()*v.vertex - lower*bstline;
	 } else  {
            if (lower>0.0) 
              v.vertex = -ubest.lb()*v.vertex + lower*bstline;
            else
	      v.vertex = -ubest.ub()*v.vertex + lower*bstline;
         }
	 v.reduce_vector();
         v.addFct(codeFacet);
      }
      /* ajout du nouveau sommet et des liens */
      std::forward_list<DDvertex>::iterator itVertex = vertices.begin();
      std::forward_list<DDvertex>::iterator newVertex = 
		this->addDDvertex(ubest.lb()<0.0 ? bstline : (-bstline));
      for (;itVertex!=vertices.end();++itVertex) {
          this->addDDlink(itVertex,newVertex);
      }
      for (Index a=0;a<codeFacet;a++) {
          newVertex->addFct(a);
      }
      /* ce nouveau sommet a une face spéciale, "négation" de la contrainte,
         de code 0 */
      newVertex->addFct(reffacets.size());
      reffacets.push_back(0);
      /* suppression de la ligne */
      this->lines.pop_back();
      return 1;
   }
   /* end of line treatments , the following 
      is when no line is available */
   int nbGT=0;
   bool notempty=false;
   for (DDvertex &vt : vertices) {
      vt.lambda = facet.dot(vt.vertex);
      notempty |= (vt.lambda.lb()<=0.0);
      if (vt.lambda.ub()<0.0) {
         vt.status = DDvertex::VERTEXLT;
         continue;
      }
      if (vt.lambda.lb()<0.0) {
         vt.status = DDvertex::VERTEXON;
         continue;
      }
      if (vt.lambda.lb()<this->tolerance) {
         vt.status = DDvertex::VERTEXGE;
         continue;
      }
      nbGT++;
      vt.status=DDvertex::VERTEXGT;
   }
   if (!notempty) { this->empty=true; vertices.clear(); return -1; }
   if (nbGT==0) return 0;

   std::forward_list<DDvertex>::iterator itVertex = vertices.begin();
   std::priority_queue<std::forward_list<DDvertex>::iterator,
		std::vector<std::forward_list<DDvertex>::iterator>,
		CmpDD> stack;
   Index refFacet=-1;
   while (itVertex!=vertices.end()) {
//       bool inBoundary;
       std::forward_list<DDvertex>::iterator actVertex;
       if (stack.size()>0) {
          actVertex = stack.top();
          stack.pop();
          if (actVertex->status!=DDvertex::VERTEXSTACK) continue;
//          inBoundary=true;
       } else {
          actVertex=itVertex;
          ++itVertex;
          if (actVertex->status!=DDvertex::VERTEXGT) continue;
          refFacet = reffacets.size(); /* new facet identification */
          reffacets.push_back(idFacet);
//          inBoundary=false;
       }
#ifdef DEBUG_CODAC2_DD
       std::cout << " elimination sommet : " << actVertex->Id << ", " << actVertex->links.size() << "liens \n";
#endif
       actVertex->status = DDvertex::VERTEXREM;
       std::vector<std::pair<std::forward_list<DDvertex>::iterator,
			std::vector<Index>>>
		adjacentVertices;
       for (std::forward_list<DDvertex>::iterator lnk : actVertex->links) {
           DDvertex &destlink = *lnk;
#ifdef DEBUG_CODAC2_DD
       std::cout << " lien : " << destlink.Id << "\n";
#endif
           if (destlink.status==DDvertex::VERTEXGT) {
               destlink.status=DDvertex::VERTEXSTACK;
               stack.push(lnk);
           }
           if (destlink.status==DDvertex::VERTEXGE ||
               destlink.status==DDvertex::VERTEXSTACK) {
	       std::cout << " GE ou STACK " << destlink.lambda << "\n";
               /* add the facet to the vertex, if it does not already exists */
	       adjacentVertices.push_back(std::pair(lnk,fctSon(-1,lnk->fcts,actVertex->fcts)));
               destlink.addFct(refFacet);
	       destlink.removeLnk(actVertex);
	       continue;
           }
           assert(destlink.status!=DDvertex::VERTEXREM);
			 /* should not happen */
           IntervalVector partV=destlink.vertex;
           if (destlink.status==DDvertex::VERTEXON) {
             IntervalRow x1(facet);
             MulOp::bwd(Interval(-oo,0.0),x1,partV);
             if (partV.is_empty()) {
	       std::cout << " Sommet conserve " << facet << " " << destlink.vertex << " " << destlink.lambda << "\n";
#if 0
              IntervalVector partV2=destlink.vertex;
              IntervalRow x2(facet);
              MulOp::bwd(Interval(0.0,+oo),x2,partV2);
	       std::cout << " reste " << partV2 << "\n";
#endif
	       adjacentVertices.push_back(std::pair(lnk,fctSon(-1,lnk->fcts,actVertex->fcts)));
               destlink.addFct(refFacet);
	       destlink.removeLnk(actVertex);
	       continue;
             }
           }
           if (-destlink.lambda.lb()/actVertex->lambda.lb()<tolerance) {
	       std::cout << " Extension sommet " << -destlink.lambda.lb()/actVertex->lambda.lb() << "\n";
	      /* extend the vertex, and add fct */
	      adjacentVertices.push_back(std::pair(lnk,fctSon(-1,lnk->fcts,actVertex->fcts)));
              destlink.vertex |=
                 (partV - 
		(destlink.lambda.lb()/actVertex->lambda)*actVertex->vertex);
              destlink.addFct(refFacet);
	      destlink.removeLnk(actVertex);
	      continue;
           }
	   /* creation of a new vertex */
	   IntervalVector newV = -destlink.lambda.lb()*actVertex->vertex
		      + actVertex->lambda.lb()*partV;
           assert(!newV.is_empty());
#if 0
           if (newV.is_empty()) {
	      std::cout << "  newV : " << newV << "\n";
              std::cout << destlink.lambda << " " << actVertex->vertex << " " << actVertex->lambda << " " << partV << "\n";
	   }
#endif
           std::forward_list<DDvertex>::iterator newVx = 
		     this->addDDvertex(newV);
#ifdef DEBUG_CODAC2_DD
       std::cout << " creation new V " << newV << "id : " << newVx->Id << "\n";
#endif
	   newVx->lambda=0.0;
	   destlink.changeLnk(actVertex,newVx);
	   newVx->addLnk(lnk,false);
	   newVx->fcts = fctSon(-1,destlink.fcts,actVertex->fcts);
	   adjacentVertices.push_back(std::pair(newVx,newVx->fcts));
           newVx->addFct(refFacet);
           newVx->status=DDvertex::VERTEXGE;
       }
#ifdef DEBUG_CODAC2_DD
       std::cout << " fin traitement liens : " << adjacentVertices.size() << "\n";
#endif
       std::vector<PairMaxSets<std::forward_list<DDvertex>::iterator>> maxset;
       for (Index i=0;i<=(Index)adjacentVertices.size()-1;i++) {
#ifdef DEBUG_CODAC2_DD
         std::cout << " sommet : " << adjacentVertices[i].first->Id << " : ";
         for (Index u : adjacentVertices[i].second) {
            std::cout << u << ",";
         }        
         std::cout << std::endl;
#endif
         for (Index j=i+1;j<(Index)adjacentVertices.size();j++) {
              PairMaxSets<std::forward_list<DDvertex>::iterator> 
			actuel(adjacentVertices[i].first,
			       adjacentVertices[j].first,
			       adjacentVertices[i].second,
			       adjacentVertices[j].second,refFacet);
	      if (actuel.elems.size()<fdim-3-lines.size()) continue;
	      bool placed=false;
              unsigned long int k=0,l=0;
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
#ifdef DEBUG_CODAC2_DD
       std::cout << " boucle maxset : " << actVertex->Id << "\n";
#endif
       for (const PairMaxSets<std::forward_list<DDvertex>::iterator>
			 &pm : maxset) {
#ifdef DEBUG_CODAC2_DD
	      if (pm.elems.size()<fdim-3-lines.size()) {
           std::cout << " maxset strange : " << pm.a->Id << " " << pm.b->Id << " " << pm.elems.size() <<  std::endl;
            }
		std::cout << pm.a->Id << " " << pm.b->Id << "\n";
#endif
	    bool res = this->addDDlink(pm.a,pm.b);
#ifdef DEBUG_CODAC2_DD
            std::cout << (res ? "   added " : "   already here ") << "\n";
#endif

       }
   }
       
   vertices.remove_if([](const DDvertex &v) 
			{ return v.status==DDvertex::VERTEXREM; });
   return 1;
}
   

/* TODO : make a specific computation */
#if 0
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
#endif

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
      if (vt.vertex[0]!=0) nbvx++; else nbrays++;
      os << " " << vt.Id << " : " << normalise(vt.vertex) << std::endl;
      os << " fct : (";
      for (const auto &a : vt.fcts) {
         os << a << ":" << build.reffacets[a] << ",";
      }
      os << ")" << std::endl;
      os << " lnks: (";
      for (const auto &a : vt.links) {
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
   std::vector<PairMaxSets<Index>> maxset;
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
   for (const PairMaxSets<Index> &pm : maxset) {
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

