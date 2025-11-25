/**
 * Linear programming problem, use coin-or clp
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
#include <CoinBuild.hpp>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wregister"
#include <ClpModel.hpp>
#include <ClpSimplex.hpp>
#pragma GCC diagnostic pop


#include "codac2_BoolInterval.h"
#include "codac2_Facet.h"
#include "codac2_clp.h"

#undef CODAC2_EXTENSION_POLYTOPE_CLP_DEBUG

namespace codac2 {

/**********  Problem building funtions ***********/

LPclp::LPclp (Index dim, std::shared_ptr<CollectFacets> &facets,
      	  const Row &objvect) :
   nbRows(facets->nbfcts()), nbCols(dim),
   Afacets(facets), objvect(objvect),
   cststat(nbRows,0),
   rowBasis(nbRows+nbCols,false), built(false), model(nullptr),
   primalSol(nbCols), dualSol(nbRows+nbCols), 
   primalRay(nbCols), dualRay(nbRows+nbCols),
   bbox(nbCols)
{
}

LPclp::LPclp (Index dim, std::shared_ptr<CollectFacets> &facets) :
   nbRows(facets->nbfcts()), nbCols(dim),
   Afacets(facets), objvect(Row::Zero(nbCols)),
   cststat(nbRows,0),
   rowBasis(nbRows+nbCols,false), built(false), model(nullptr),
   primalSol(nbCols), dualSol(nbRows+nbCols), 
   primalRay(nbCols), dualRay(nbRows+nbCols),
   bbox(nbCols)
{
}

LPclp::LPclp (Index dim, std::shared_ptr<CollectFacets> &facets,
		const IntervalVector &box) :
   nbRows(facets->nbfcts()), nbCols(dim),
   Afacets(facets), objvect(Row::Zero(nbCols)),
   cststat(nbRows,0),
   rowBasis(nbRows+nbCols,false), built(false), model(nullptr),
   primalSol(nbCols), dualSol(nbRows+nbCols),
   primalRay(nbCols), dualRay(nbRows+nbCols),
   bbox(box)
{
}

LPclp::LPclp(const Matrix &mat, const Row &objvect,
                const Vector &rhsvect, const std::vector<Index> &eqSet) : 
   nbRows(mat.rows()), nbCols(mat.cols()), 
   Afacets(std::make_shared<CollectFacets>(mat,rhsvect,eqSet)),
   objvect(objvect),
   cststat(nbRows,0), 
   rowBasis(nbRows+nbCols,false), built(false), model(nullptr),
   primalSol(nbCols), dualSol(nbRows+nbCols),
   primalRay(nbCols), dualRay(nbRows+nbCols),
   bbox(nbCols)
{
}

LPclp::LPclp(const Matrix &mat, 
                const Vector &rhsvect, const std::vector<Index> &eqSet) : 
   nbRows(mat.rows()), nbCols(mat.cols()), 
   Afacets(std::make_shared<CollectFacets>(mat,rhsvect,eqSet)),
   objvect(Row::Zero(nbCols)),
   cststat(nbRows,0), 
   rowBasis(nbRows+nbCols,false), built(false), model(nullptr),
   primalSol(nbCols), dualSol(nbRows+nbCols),
   primalRay(nbCols), dualRay(nbRows+nbCols),
   bbox(nbCols)
{
}

LPclp::LPclp(Index dim) :
   nbRows(0), nbCols(dim), 
   objvect(Row::Zero(nbCols)),
   built(false), model(nullptr),
   primalSol(nbCols), primalRay(nbCols),
   bbox(nbCols)
{
}

#if 0   
/* no copy constructor it shares the facet pointer */
LPclp::LPclp(const LPclp &P) :
    nbRows(P.nbRows), nbCols(P.nbCols),
    Afacets(P.Afacets), objvect(P.objvect), cststat(P.cststat),
    rowBasis(nbRows+nbCols,false), built(false), model(nullptr),
    primalSol(nbCols), dualSol(nbRows+nbCols),
    primalRay(nbCols), dualRay(nbRows+nbCols),
    bbox(P.bbox), timeout(P.timeout), maxIteration(P.maxIteration),
    tolerance(P.tolerance)
{
}
#endif

LPclp::~LPclp() {
   delete(model);
}

void LPclp::addFacetToCoinBuild(CoinBuild &buildObject,
	     const Facet &facet, lp_cststatus status,
             int *row2Index, double *row2Vals) const {
    double bndUp = ((status[REMOVED] || status[INACTIVE]
	|| status[REDUNDANT]) ?  COIN_DBL_MAX : facet.second.get_rhs());
    double bndLow = (!facet.second.is_eqcst() ||
	  (status[REMOVED] || status[INACTIVE] || status[REDUNDANT]) ?
                        -COIN_DBL_MAX : facet.second.get_rhs());
    if (built_emptytest) {
      for (Index j=0;j<nbCols;j++) { 
	/* direct copy, needed if emptytest is true */
           row2Vals[j]=facet.first.get_row()[j];
      }
      row2Vals[nbCols]=(facet.second.is_eqcst() ? 0.0 : 1.0); 
      buildObject.addRow(nbCols+1, row2Index, row2Vals, bndLow, bndUp);
   } else {
      buildObject.addRow(nbCols, row2Index, 
		facet.first.get_row().data(), bndLow, bndUp);
   }
}


Index LPclp::addConstraint(const Facet &facet) {
   std::pair<Index,bool> a = this->Afacets->insert_facet(facet);
   if (a.first==0) return 0;
   this->cststat.push_back(0);
   this->dualSol.resize(nbRows+nbCols+1);
   this->dualRay.resize(nbRows+nbCols+1);
   if (built) {
     CoinBuild buildObject;
     Index nbCols2=nbCols+(built_emptytest ? 1 : 0);
     int* row2Index = new int[nbCols2];
     for(Index i=0;i<nbCols2;i++){row2Index[i]=i;}
     double* row2Vals = built_emptytest ? (new double[nbCols2]) : nullptr;
     this->addFacetToCoinBuild(buildObject,facet,0,row2Index,row2Vals);
     delete row2Index;
     if (built_emptytest) delete row2Vals;
     model->addRows(buildObject);
     this->status=1<<CHANGED;
   }
   return (nbRows++);
}

Index LPclp::addConstraint(const Row &vst, double rhs, bool isEq) {
   return this->addConstraint(Facet_::make(vst,rhs,isEq));
}

Index LPclp::updateConstraint(Index id) {
   const Facet &fct = (*(*this->Afacets)[id]);
   if (id<nbRows) {
       if (!built) return id;
       cststat[id][INACTIVE]=false;
       model->setRowUpper(id,fct.second.get_rhs());
       if (fct.second.is_eqcst()) 
  	   model->setRowLower(id,fct.second.get_rhs());
       this->status=1<<CHANGED;
       return id;
   }
   assert_release(id==nbRows);
   nbRows++;
   this->dualSol.resize(nbRows+nbCols);
   this->dualRay.resize(nbRows+nbCols);
   this->cststat.push_back(0);
   if (built) {
     CoinBuild buildObject;
     Index nbCols2=nbCols+(built_emptytest ? 1 : 0);
     int* row2Index = new int[nbCols2];
     for(Index i=0;i<nbCols2;i++){row2Index[i]=i;}
     double* row2Vals = built_emptytest ? (new double[nbCols2]) : nullptr;
     this->addFacetToCoinBuild(buildObject,fct,0,row2Index,row2Vals);
     delete row2Index;
     if (built_emptytest) delete row2Vals;
     model->addRows(buildObject);
     this->status=1<<CHANGED;
   }
   return nbRows;
}

void LPclp::setObjective(const Row &objvect) {
    this->objvect=objvect;
    if (built && !built_emptytest) { /* change the values of the objective */
       for (Index i =0; i<nbCols; i++)
            model->setObjectiveCoefficient(i, objvect[i]);
       this->status=1<<CHANGED;
    }
}

void LPclp::setActive(Index cst, bool state) {
    assert_release(!cststat[cst][REMOVED]);
    if (state) {
       if (!cststat[cst][INACTIVE]) return;
       cststat[cst][INACTIVE]=false;
       if (built) {
           const auto &fct = (*Afacets)[cst];
           model->setRowUpper(cst,fct->second.get_rhs());
           if (fct->second.is_eqcst()) 
		model->setRowLower(cst,fct->second.get_rhs());
       }
       this->status=1<<CHANGED;
    } else {
       if (cststat[cst][INACTIVE]) return;
       cststat[cst][INACTIVE]=true;
       if (built) {
            model->setRowUpper(cst,DBL_MAX);
            model->setRowLower(cst,-DBL_MAX);
       }
       this->status=1<<CHANGED;
    }
}

void LPclp::update_model_bbox() {
   for (Index i=0;i<nbCols;i++) {
      model->setColumnBounds(i,bbox[i].lb(), bbox[i].ub());
   }
}


/******* coin-or building and configuration **********/

void LPclp::buildModel(bool emptytest) {
   assert_release(model==nullptr);
   model = new ClpSimplex();
   
   model->setLogLevel(0);
   model->setOptimizationDirection(-1); /* Maximize */

   model->setMaximumIterations(this->maxIteration);
   model->setMaximumSeconds(this->timeout);
   model->setPrimalTolerance(this->tolerance);
   model->setDualTolerance(this->tolerance);

   /* set dimension */
   Index nbCols2=nbCols+(emptytest ? 1 : 0);
   model->resize(0,nbCols2);
   this->update_model_bbox();
   if (emptytest) {
      model->setColumnBounds(nbCols,-COIN_DBL_MAX,COIN_DBL_MAX);
   }

   /* adding the rows using BuildObject ? */
   built_emptytest=emptytest;
   CoinBuild buildObject;
   int* row2Index = new int[nbCols2];
   for(Index i=0;i<nbCols2;i++){row2Index[i]=i;}
   double* row2Vals = emptytest ? (new double[nbCols2]) : nullptr;
   for(Index k=0; k<nbRows; k++) {
         this->addFacetToCoinBuild(buildObject,
		*((*Afacets)[k]),cststat[k],row2Index,row2Vals);
   }
   model->addRows(buildObject);
   delete row2Vals;
   delete row2Index;

   /* set the objective */
   for (Index i =0; i<nbCols; i++) {
      model->setObjectiveCoefficient(i, (emptytest ? 0.0 : objvect[i]));
   }
   if (emptytest) 
      model->setObjectiveCoefficient(nbCols, 1.0);
   built=true;
}

/* adding/removing emptiness testing on a built model */
void LPclp::setModel(bool emptytest) {
   if (!built) {
      buildModel(emptytest);
      return;
   }
   if (emptytest==built_emptytest) return;
   if (emptytest) {
      int* col2Index=new int[nbRows];
      for(Index i=0;i<nbRows;i++){col2Index[i]=i;}
      double *col2Vals = new double[nbRows];
      for(Index i=0;i<nbRows;i++)
  	 { col2Vals[i]=(*Afacets)[i]->second.is_eqcst() ? 0.0 : 1.0; }
      model->addColumn
	(nbRows,col2Index,col2Vals,-COIN_DBL_MAX,COIN_DBL_MAX,1.0);
      for (Index i =0; i<nbCols; i++)
          model->setObjectiveCoefficient(i, 0.0);
   } else {
      int nb=(int) nbCols;
      model->deleteColumns(1,&nb);
      for (Index i =0; i<nbCols; i++)
          model->setObjectiveCoefficient(i, objvect[i]);
   }
   built_emptytest=emptytest;
}


/* reinitialise the solution */
void LPclp::reset_solution() {
   status=0;
   rowBasis.assign(nbRows+nbCols,false);  
   Valobj=Interval();
   primalSol.setZero();
   dualSol.setZero();
   primalRay.setZero();
   dualRay.setZero();
}


/********** Solving functions *******************/

/***** launch one solving, returns the status *******/
LPclp::lp_result_stat LPclp::solve(bool checkempty, int option) {
   this->setModel(checkempty);
   reset_solution();
   model->primal(0,option | 1);
   int stat = model->status();
   if (stat==-1 || stat>2) { 
	status=1<<ERROR_LPCOIN; 
        return status;
   }
   return testGeneralApproach(model->status(),checkempty);
}

int LPclp::minimize_eqpolytope() {
    std::vector<IntervalRow> triMat;
    std::vector<Interval> triRhs;
    Index nbCsts=0;
    std::vector<Index> mainCol;
    for (Index i=0; i<nbCols; i++) {
       if (bbox[i].is_degenerated()) {
          Row r = Row::zero(nbCols);
          r[i]=1.0;
          triMat.push_back(std::move(r));
          triRhs.push_back(bbox[i].mid());
          mainCol.push_back(i);
          nbCsts++;
       }
    }
    for (Index i = 0; i<nbRows; i++) {
       if (!(*this->Afacets)[i]->second.is_eqcst()) continue; 
       if (cststat[i][REMOVED] || cststat[i][INACTIVE] || cststat[i][REDUNDANT])
           continue;
       triMat.emplace_back((*this->Afacets)[i]->first.get_row());
       triRhs.emplace_back((*this->Afacets)[i]->second.get_rhs());
       for (Index j=0;j<nbCsts;j++) {
          if (mainCol[j]<0) continue;
          Interval pivot = triMat[nbCsts][mainCol[j]]/triMat[j][mainCol[j]];
          triMat[nbCsts] -= triMat[j]*pivot;
          triMat[nbCsts][mainCol[j]]=0.0;
          triRhs[nbCsts] -= triRhs[j]*pivot;
       }
       Index bCol=-1;
       double bVal=0.0;
       for (Index j=0;j<nbCols;j++) {
          double a = triMat[nbCsts][j].mig();
          if (a>bVal) {
             bVal=a;
	     bCol=j;
          }        
       }
       if (bCol==-1) {
          if (triRhs[nbCsts].contains(0)) {
             cststat[i][REDUNDANT]=true;
             continue;
          }
          for (Index j=0;j<nbCols;j++) {
              if (triMat[nbCsts][j].mag()==0.0) continue;
              triRhs[nbCsts] -= bbox[j]*triMat[nbCsts][j];
              if (triRhs[nbCsts].contains(0)) break;
          }
          if (!triRhs[nbCsts].contains(0)) return -1;
       }
       mainCol.push_back(bCol);
       nbCsts++;
    }
    return nbCsts;
}

int LPclp::check_emptiness() {
   lp_result_stat stat=this->solve(true);
   if (stat[EMPTY]) return -1;
   return 0;
}

int LPclp::minimize_box() {
    /* minimize box */
    for (Index i=0;i<nbCols;i++) {
       Row a = Row::Zero(nbCols);
       a[i]=1.0;
       this->setObjective(a);
       lp_result_stat stat=this->solve(false,(i!=0 ? 4 : 0));
       if (stat[EMPTY]) return -1;
       if (stat[BOUNDED]) {
          bbox[i] = min(bbox[i],this->Valobj.ub());
       }
    }
    for (Index i=0;i<nbCols;i++) {
       Row a = Row::Zero(nbCols);
       a[i]=-1.0;
       this->setObjective(a);
       lp_result_stat stat=this->solve(false,4);
       if (stat[EMPTY]) return -1;
       if (stat[BOUNDED]) {
          bbox[i] = min(bbox[i],this->Valobj.ub());
       }
    }
    this->update_model_bbox();
    return 0;
}

int LPclp::minimize_polytope(const  Interval &tolerance, bool checkempty,
		bool checkbox) {
    int nb = this->minimize_eqpolytope();
    if (nb==-1) return -1;
    if (checkempty) 
       if (check_emptiness()==-1) return -1;
    if (checkbox)
        if (this->minimize_box()==-1) return -1;
    for (Index i=nbRows-1;i>=0;i--) {
       const auto &fct = (*this->Afacets)[i];
       if (cststat[i][REMOVED] || 
	   cststat[i][INACTIVE] || 
	   cststat[i][REDUNDANT])
           continue;
       if (fct->second.is_eqcst()) { nb++; continue; }
       this->setActive(i,false);
       this->setObjective(fct->first.get_row());
       lp_result_stat stat=this->solve(false,4);
       if (stat[EMPTY]) return -1;
       if (stat[BOUNDED]) {
           if (this->Valobj.ub()<=fct->second.get_rhs()+tolerance.ub() 
		&& this->Valobj.lb()<=fct->second.get_rhs()+tolerance.lb()) {
              cststat[i][REDUNDANT]=true;
              continue;
           }
       }
       nb++;
       this->setActive(i,true);
    }
    return nb;
}


/*************************************************************/
/****     Validation of the results of coin-or ***************/
/*************************************************************/


/************ correction of inverse matrix *******************/
/**** basis is a nb*N matrix       and basisInverse is  ******/
/**** an approximation of a nb*nb submatrix of basis    ******/
/**** the columns i part of the submatrix are defined   ******/
/**** by (*revbasics)[i]!=1 ( if nb!=N )                ******/
/**** modify basisInverse to contain the "true" inverse ******/
/**** of the submatrix of basis                         ******/
/**** if A is the submatrix and B is the approximate    ******/
/**** inverse, we compute Id-BA = Err                   ******/
/**** then basisInverse = (Id+Eps)*B with               ******/
/**** Eps = Err*(Id+Err+Err^2+...)                      ******/

/**** code for approximation of infinite sum enclosure ****/
/**** should be in codac2 ?                 	       ****/
IntervalMatrix LPclp::infinite_sum_enclosure (const IntervalMatrix& M) {
   size_t n = M.cols();
   Matrix TmpM = M.mig();
#if 0
   // computing the mag of the eigen matrix
   auto sit = M.reshaped().cbegin();
   auto dit = TmpM.reshaped().begin();
   for (;sit!=M.reshaped().cend(); ++sit,++dit) {
       (*dit) = (*sit).mag();
   }
#endif

   // Floyd algorithm below
   for (size_t k=0;k<n;k++) {
      double mkk = TmpM(k,k);
      if (mkk>0.0) {
        for (size_t r=0;r<n;r++) {
           if (r==k) continue;
           if (TmpM(r,k)>0.0 && TmpM(r,k)<mkk)
	      TmpM(r,k) = mkk;
	}
        for (size_t c=0;c<n;c++) {
           if (c==k) continue;
           if (TmpM(k,c)>0.0 && TmpM(k,c)<mkk)
	      TmpM(k,c) = mkk;
	}
      }
      for (size_t c=0;c<n;c++) {
         if (c==k) continue;
         if (TmpM(k,c)==0.0) continue;
         double mkc=TmpM(k,c);
         for (size_t r=0;r<n;r++) { 
            if (r==k) continue;
            if (TmpM(r,k)==0.0) continue;
            double vm = std::max(mkc,TmpM(r,k));
            if (TmpM(r,c)<vm) TmpM(r,c)=vm;
         }
      }
   }
   
   // and now the overall computation 
   IntervalMatrix res=IntervalMatrix::Zero(n,n);
   Interval v(-1.0,1.0);
   for (size_t c=0;c<n;c++) {
      for (size_t r=0;r<n;r++) {
         if (TmpM(r,c)==0.0) continue;
         Interval cnt=0.0;
         for (size_t k=0;k<n;k++) {
            cnt += std::min(TmpM(r,k),TmpM(k,c));
         }
         if (cnt.ub()==0.0) {
            res(r,c)=M(r,c); 
         } else {
            res(r,c)=(v*TmpM(r,c))/(1.0-cnt);
         }
      }
   }
   return res;
}


void LPclp::correctBasisInverse(IntervalMatrix &basisInverse,
           const IntervalMatrix &basis) const {
     IntervalMatrix Id=IntervalMatrix::Identity(basis.rows(),basis.cols());
     IntervalMatrix ErrorMat=Id;
     ErrorMat -= basisInverse*basis;

    IntervalMatrix E = LPclp::infinite_sum_enclosure(ErrorMat);
    IntervalMatrix Eps = ErrorMat*(Id+E);
			 /* or ErrorMat*(Id+Errormat*(Id+E))? */
    basisInverse = (Id+Eps)*basisInverse;
}   

Interval LPclp::dotprodrhs(const IntervalRow &d) const {
   Interval r(0.0);
   for (Index i=0;i<nbRows;i++) {
      r += d[i]*(*this->Afacets)[i]->second.get_rhs();
   }
   for (Index i=0;i<nbCols;i++) {
      r += d[i+nbRows]*bbox[i];
   }
   return r;
}


/* 
  try to check a product dualVect * facets(extended) + initSum = 0,
  correcting it if needed using basisInverse
  note : initially, dualVect is punctual !
 */
BoolInterval LPclp::checkDualForm(IntervalRow &dualVect,
		IntervalRow &initSum, const IntervalMatrix &basisInverse,
                const std::vector<Index> &wholeBasis, Index nbRowsInBasis,
                bool checkempty) const {
   for (Index i=0;i<nbRows;i++) {
       initSum.head(nbCols)+= dualVect[i]*(*Afacets)[i]->first.get_row();
       if (checkempty && !(*Afacets)[i]->second.is_eqcst())
		 initSum[nbCols]+=dualVect[i];
   }
   for (Index i=0;i<nbCols;i++) {
       initSum[i]+= dualVect[i+nbRows];
   }
   /* calcul de "l'erreur" à appliquer sur dualVect */
   IntervalRow error = initSum*basisInverse;
   BoolInterval ok=BoolInterval::TRUE;
   for (Index i=0;i<nbRowsInBasis;i++) {
       int r = wholeBasis[i];
       dualVect[r] -= error[i]; /* modify the ray with the error,
			       the value must stay positive... */
       if ((*Afacets)[r]->second.is_eqcst())
			 continue; /* ... except for equalities */
       if (dualVect[r].ub()<0.0) { return BoolInterval::FALSE; }
       if (dualVect[r].lb()<0.0) { ok = BoolInterval::UNKNOWN; }
   }
   /* vérifie les variables  colonnes, qui doivent
		  être libres si bbox est unbounded,
		 ainsi que les contraintes inactives */
   for (Index i=nbRowsInBasis;i<error.size();i++) {
     if (error[i].mig()==0.0) continue;
     int r = wholeBasis[i];
     dualVect[r] -= error[i];
     int c = r - nbRows;
     if (bbox[c].ub()==+oo) {
       if (dualVect[r].lb()>0.0) { return BoolInterval::FALSE; }
       if (dualVect[r].ub()>0.0) { ok = BoolInterval::UNKNOWN; }
     }
     if (bbox[c].lb()==-oo) {
       if (dualVect[r].ub()<0.0) { return BoolInterval::FALSE; }
       if (dualVect[r].lb()<0.0) { ok = BoolInterval::UNKNOWN; }
     }
   }
   for (Index i=0;i<nbRows;i++) {
     if ((cststat[i][INACTIVE] || cststat[i][REMOVED]) &&
        (dualVect[i].mig()>0.0)) { ok=BoolInterval::FALSE; break; }
   }
   return ok;
}



/* this function tests all part sequentially. It starts
   from the result of CLP solve */
LPclp::lp_result_stat LPclp::testGeneralApproach(int stat,bool checkempty) {
    /* 1. récupérer les données :
          * wholeBasis : la base utilisée pour inverser contraintes
            et RHS
            (à la fois lignes, et colonnes)
            contient nbCols(+1 si checkempty) éléments
          * basis : la matrice de la "base"
                  nbCols(+1) * nbCols(+1)
          * basisInverse : son inverse
     */
    Index nbCols2 = nbCols+(checkempty ? 1 : 0);
    std::vector<Index> wholeBasis; 
    wholeBasis.reserve(nbCols2);
    Index nbRowsInBasis=0;
    /* row basics */
    for (Index i=0;i<nbRows;i++) {
	if (model->getRowStatus(i)!=ClpSimplex::Status::basic) {
	   wholeBasis.push_back(i);
           nbRowsInBasis++;
           rowBasis[i]=true;
        } 
        dualSol[i]=model->dualRowSolution()[i];
    }
    /* col basics */
    double slack=0.0; /* only for emptiness checking */
    for (Index i=0;i<nbCols2;i++) {
	if (model->getColumnStatus(i)!=ClpSimplex::Status::basic) {
	 	wholeBasis.push_back(i+nbRows);
                rowBasis[i+nbRows]=true;
        }
        if (i<nbCols) {
             dualSol[i+nbRows]=model->dualColumnSolution()[i];
             primalSol[i]=model->primalColumnSolution()[i];
        }
        else slack = model->primalColumnSolution()[i];
    }
    assert((Index) wholeBasis.size()==nbCols2);
    
    std::vector<Index> revbasics(nbRows,-1);
    int *basics=new int[nbRows];
    model->getBasics(basics); 
    for (Index i=0;i<nbRows;i++) {
       if (cststat[i][REMOVED])
             continue;
       if (basics[i]<nbCols2) revbasics[basics[i]]=i;
		/* revbasis[c] = l'indice de la ligne dans 
		model->BInvRow qui "inverse" la c-ème colonne de Amat */
    }
    delete basics;

    IntervalMatrix basis=
		IntervalMatrix::Zero(nbCols2,nbCols2);
    IntervalMatrix basisInverse=
		IntervalMatrix::Zero(nbCols2,nbCols2);

    /* fill the basis with only the rows first */
    for (Index c=0;c<nbCols;c++) {
       for (Index i=0;i<nbRowsInBasis;i++) {
          int r = wholeBasis[i];
          basis(i,c)=(*this->Afacets)[r]->first.get_row()[c];
       }
    }
    if (checkempty) {
       for (Index i=0;i<nbRowsInBasis;i++) basis(i,nbCols)=1.0;
    }
    /* fill basisInverse with the inverse computed by CLP */
    double* invCol = new double[nbRows];
    for (Index c=0;c<nbCols2;c++) {
       if (revbasics[c]==-1) continue;
       model->getBInvRow(revbasics[c],invCol); 
       for (Index k=0;k<nbRowsInBasis;k++) {
           basisInverse(c,k)=invCol[wholeBasis[k]];
       }
       /* complete basisInverse for columns in the basis:
        * in this case, we want 
		basisInverse(c)*basis.col(wholeBasis[k]-nbRows)=0 
		with basis(k,wholeBasis[k]-nbRows)=1 (put later) */
       for (Index k=nbRowsInBasis;k<nbCols2;k++) {
           basisInverse(c,k)=
	    (-basisInverse.row(c).dot(basis.col(wholeBasis[k]-nbRows))).mid();
       }
    }
    /* now put the 1.0 for the remaining rows of basis 
		(columns of basisInverse) */
    for (Index i=nbRowsInBasis;i<nbCols2;i++) {
         int a = wholeBasis[i];
         basis(i,a-nbRows)=1.0;
         basisInverse(a-nbRows,i)=1.0;
    }

    /* correction of basisInverse to get the "true" inverse */
    correctBasisInverse(basisInverse,basis);

#ifdef CODAC2_EXTENSION_POLYTOPE_CLP_DEBUG
    std::cout << "stat: " << stat << "\n";
    std::cout << "basis: " << basis << "\n";
    std::cout << "basisInverse: " << basisInverse << "\n";
    std::cout << "primalSol: " <<  primalSol << "\n";
    std::cout << "dualSol: " <<  dualSol << "\n";
    std::cout << "slack: " <<  slack << "\n";
#endif
    /*
       Before checking all cases, we adapt the problem when
       checkempty = true. Normally stat!=1
       2 cases : stat=2, means that the polytope in
         unbounded ; if slack<0.0, to check the nonemptiness we must
         change primalSol. Furthermore, we fill primalRay already.
       if checkempty = true and we have unboundness, we may have
       the problem of having a negative slackvariable => modify it.
      */
    if (checkempty) {
        assert_release(stat!=1);
        if (stat==2) {
           double *ray = model->unboundedRay();
           assert_release(ray!=nullptr);
           assert_release(ray[nbCols]>0.0);
           for (Index i=0;i<nbCols;i++) primalRay[i]=ray[i];
#ifdef CODAC2_EXTENSION_POLYTOPE_CLP_DEBUG
	   std::cout << "Primal ray : " << primalRay << " " << ray[nbCols] << "\n";
#endif
           if (slack<0.1) primalSol += ((slack+1.0)/ray[nbCols])*primalRay;
           delete ray; 
        } else {  /* stat=0 */
          /* slack>=0 : normally, it's ok; slack<0.0 => empty ,
               dual value should prove the empty case */
           if (slack<0.0) { 
              dualRay=dualSol;
              stat=1;
           }
        }
    }

    /* basis rows : nb constraints, ordonnées par wholeBasis
                    des variables, = wholeBasis[i] - nbRows 
             cols : variables */
    /* basisinverse rows : variables
                    cols : nb constraints, ordonnées par wholeBasis
                           variables, wholeBasis[i]-nbRows.
              obj*basisInverse : combinaisons de contraintes pour 
				 obtenir obj
              basisInverse*rhs : point qui satisfait rhs */

    /* stat = 0 : CLP said optimal 
       stat = 1 : CLP said empty
       stat = 2 : CLP said unbounded */
    if (stat==1) { /* test only dual infeasibility */
         /* if checkempty=true, we already have the ray,
            otherwise we build it */
         double *ray=nullptr;
         if (!checkempty) {
           double *ray = model->infeasibilityRay(true);
           assert_release(ray!=nullptr);  /* false => bug of CLP */
             /* true for infeasibilityRay means that we have also
                the negation of the "error" term after the rows */
           for (Index i=0;i<nbRows+nbCols;i++) { dualRay[i]=ray[i]; }
             /* compute the error term , and reverse it with basis. */
             /* val : csts => val * Amat : variables (~0)  
	 	val*Amat*basisInverse => basisInverse.cols */
         }
#ifdef CODAC2_EXTENSION_POLYTOPE_CLP_DEBUG
         std::cout << "dualRay: " << dualRay << "\n";
#endif
         IntervalRow errorNeum=IntervalRow::Zero(nbCols2);
         if (checkempty) {
            errorNeum[nbCols]=-1.0;
         }
         BoolInterval ok = this->checkDualForm(dualRay, errorNeum,
		basisInverse, wholeBasis, nbRowsInBasis, checkempty);
#ifdef CODAC2_EXTENSION_POLYTOPE_CLP_DEBUG
         std::cout << "dualRay with error: " << dualRay << "\n";
         std::cout << "ok: " << ok << "\n";
#endif
        if ((ok & BoolInterval::TRUE)==BoolInterval::TRUE) {
          Interval product = this->dotprodrhs(dualRay);
          if (checkempty) Valobj=Interval(-oo,product.ub());         
          if (product.lb()>=0.0) ok=BoolInterval::FALSE;
          else if (product.ub()>=0) ok=BoolInterval::UNKNOWN;
        }
#ifdef CODAC2_EXTENSION_POLYTOPE_CLP_DEBUG
         std::cout << "ok after check1: " << ok << "\n";
#endif
        if (ok==BoolInterval::TRUE) {
           status[EMPTY]=true;
	   if (!checkempty) Valobj.set_empty();
           delete ray;
           return status;
        } 
        /* attempt with Neumaier */
        if (!checkempty && !bbox.is_unbounded()) {
            /* restart with the ray */
            for (Index i=0;i<nbRows+nbCols;i++) { 
                if (i<nbRows && 
			!(*Afacets)[i]->second.is_eqcst())
			 assert_release(ray[i]<0.0);
		dualRay[i]=ray[i]; 
            }
            Interval productNeum = this->dotprodrhs(dualRay);
            productNeum -= errorNeum.head(nbCols).dot(bbox);
            if (productNeum.ub()<0.0) {
              status[EMPTY]=true;
              status[EMPTY_BBOX]=true;
	      if (!checkempty) Valobj.set_empty();
              else Valobj=Interval(-oo,productNeum.ub());
              delete ray;
              return status;
            }
        } else {
	/* other approach : we try to compensate with other values */ 
          IntervalMatrix comp(nbRows,nbCols2);
          {
            IntervalVector v = IntervalVector::Ones(nbRows);
            for (Index j=0;j<nbRows;j++) {
              comp.row(j)=((IntervalRow)((*Afacets)[j]->first.get_row()))*
			basisInverse.topRows(nbCols);
              if ((*Afacets)[j]->second.is_eqcst()) v[j]=0.0;
            }
            if (checkempty) comp += v*basisInverse.row(nbCols);
          }
          /* for each epsilon value in dualSol, try to find a complement
             in comp */
          std::vector<Interval> okrow(nbRows,Interval(0.0,oo));
          int nbok=nbRows;
          for (Index j=0;j<nbRows;j++) {
             if (cststat[j][INACTIVE] || cststat[j][REMOVED]) {
			okrow[j].set_empty(); nbok--; continue; 
             }
             if (rowBasis[j]) { okrow[j].set_empty(); nbok--; continue; }
             if ((*Afacets)[j]->second.is_eqcst()) { okrow[j].init(); }
          }
          for (Index i=0;i<nbRowsInBasis;i++) {
             if (nbok==0) break;
             int r = wholeBasis[i];
             for (Index j=0;j<nbRows;j++) {
                if (okrow[j].is_empty()) continue;
                if (comp(j,i).ub()>0.0) { 
                   if (dualRay[r].lb()<=0.0) {
                        okrow[j].set_empty(); nbok--; continue;
                   } else {
			Interval u(dualRay[r].lb());
			u /= comp(j,i).ub();
                        okrow[j] &= Interval((*Afacets)[j]->second.is_eqcst() ?
				 -oo : 0.0,
				u.lb());
			if (okrow[j].is_empty()) { nbok--; continue; }
                   }
                } else if (comp(j,i).ub()==0.0) {
                   if (dualRay[r].lb()<0.0) {
                        okrow[j].set_empty(); nbok--; 
                   } 
	           continue;
                } else {
                   if (dualRay[r].lb()>=0.0) continue;
	   	   Interval u(dualRay[r].lb());
		   u /= comp(j,i).ub();
		   okrow[j] &= Interval(u.ub(),oo);
                   if (okrow[j].is_empty()) { nbok--; continue; }
                }
             }
             if (!checkempty) dualRay[r]=ray[r];
             else dualRay[r]=dualSol[r];
          }
          /* recompute errorNeum */
          if (nbok>0) {
             for (Index i=0;i<nbRows+nbCols;i++) {
                if (!checkempty) dualRay[i]=ray[i];
                    else dualRay[i]=dualSol[i];
             }
             for (Index j=0;j<nbRows;j++) {
                if (okrow[j].is_empty()) continue;
#ifdef CODAC2_EXTENSION_POLYTOPE_CLP_DEBUG
		std::cout << " corrected okrow " << j << " : " << okrow[j] << "\n";
#endif
                dualRay[j]=okrow[j].lb(); /* in theory, lb should
			be sufficient. In practice... that's a problem */
                if (dualRay[j].lb()==-oo) dualRay[j]=okrow[j].ub()-1.0;
                else
                if (dualRay[j].ub()>okrow[j].ub()) dualRay[j]=okrow[j].mid();
		break;
             }
             errorNeum=IntervalRow::Zero(nbCols2);
             for (Index i=0;i<nbRows;i++) {
                errorNeum.head(nbCols)+= dualRay[i]*(*Afacets)[i]->first.get_row();
                if (checkempty && !(*Afacets)[i]->second.is_eqcst()) 
				errorNeum[nbCols]+=dualRay[i];
             }
             if (checkempty) {
                errorNeum[nbCols]-=1.0;
             }
             IntervalRow error = errorNeum*basisInverse;
             ok=BoolInterval::TRUE;
             for (Index i=0;i<nbRowsInBasis;i++) {
                int r = wholeBasis[i];
                dualRay[r] -= error[i];
 		if (!(*Afacets)[r]->second.is_eqcst())
			 dualRay[r] &= Interval(0.0,oo);
//                if (dualSol[r].ub()<0.0) { ok = BoolInterval::FALSE; break; }
//                if (dualSol[r].lb()<0.0) { ok = BoolInterval::UNKNOWN; break; }
             }
             /* on refait les vérification sur les variables colonnes */
             if ((ok & BoolInterval::TRUE)==BoolInterval::TRUE) { 
                 /* vérifie les variables  colonnes, qui doivent
                  être libres, ainsi que les contraintes inactives */
                 for (Index i=nbRowsInBasis;i<nbCols2;i++) {
                    if (error[i].mig()==0.0) continue;
                    int r = wholeBasis[i];
                    dualRay[r] -= error[i];
                    int c = r - nbRows;
                    if (bbox[c].ub()==+oo) {
                      if (dualRay[r].lb()>0.0) 
			{ ok=BoolInterval::FALSE; break; }
                      if (dualRay[r].ub()>0.0) 
			{ ok = BoolInterval::UNKNOWN; }
                    }
                    if (bbox[c].lb()==-oo) {
                      if (dualRay[r].ub()<0.0) 
			{ ok= BoolInterval::FALSE; break; }
                      if (dualRay[r].lb()<0.0) 
		        { ok = BoolInterval::UNKNOWN; }
     	            }
                }
                if ((ok & BoolInterval::TRUE)==BoolInterval::TRUE) 
                for (Index i=0;i<nbRows;i++) {
	           if ((cststat[i][INACTIVE] || cststat[i][REMOVED]) &&
		       (dualRay[i].mig()>0.0))
				 { ok=BoolInterval::FALSE; break; }
                } 
             }
             if ((ok & BoolInterval::TRUE)==BoolInterval::TRUE) {
                Interval product = this->dotprodrhs(dualRay);
                if (checkempty) Valobj=Interval(-oo,product.ub());         
                if (product.lb()>=0.0) ok=BoolInterval::FALSE;
                else if (product.ub()>=0) ok=BoolInterval::UNKNOWN;
             }
#ifdef CODAC2_EXTENSION_POLYTOPE_CLP_DEBUG
             std::cout << "ok after check2: " << ok << "\n";
             std::cout << "dualRay: " << dualRay << "\n";
#endif
             if (ok==BoolInterval::TRUE) {
               status[EMPTY]=true;
	       if (!checkempty) Valobj.set_empty();
               delete ray;
               return status;
             } 
           }  
        } 
        delete ray;
        if (ok==BoolInterval::UNKNOWN) {
           status[EMPTY_APPROX]=true;
        } else {
           status[ERROR_PRIMAL_CHECK]=true;
        }
        /* we should try to bound the result, even if empty check failed */
        /* for that we use basisInverse directly */
        if (!checkempty)
        { 
           IntervalRow dSol= objvect*basisInverse;
           dualSol.setZero();
           bool ok=true;
           for (Index i=0;i<nbRowsInBasis;i++) {
               int r = wholeBasis[i];
               if (!(*Afacets)[r]->second.is_eqcst() && dSol[i].lb()<0.0)
			 { ok=false; break; }
               dualSol[r]=dSol[i];
           }
           if (ok) {
              for (Index i=nbRowsInBasis;i<nbCols;i++) {
                 if (dSol[i].mag()==0.0) continue;
                 int r = wholeBasis[i];
                 dualSol[r]=dSol[i];
                 int c = r-nbRows;
                 if (dSol[i].ub()>0.0 && bbox[c].ub()==+oo) { ok=false; break; }
                 if (dSol[i].lb()<0.0 && bbox[c].lb()==-oo) { ok=false; break; }
              }
           }
           if (ok) {
              Interval product = this->dotprodrhs(dualSol);
              status[BOUNDED]=true;
              Valobj = min(Valobj,product);
           }
        }
        return status;
    }
    /* following : stat=0 or stat=2, hence we look for nonemptiness 
       primalSol is supposed to give a good place */
    {
       /* 1) we modify primalSol to get a box around the "optimal" place */
       /* note: should be done only for checkempty=false
          (we do not _really_ want a box around the "center" of the polytope */
       IntervalVector errorPart(nbRows);
       for (Index i=0;i<nbRows;i++) {
           errorPart[i] = (*Afacets)[i]->first.get_row()*primalSol - 
			(*Afacets)[i]->second.get_rhs();
       } 
       if (!checkempty) {
          IntervalVector errorInBasis= IntervalVector::Zero(nbCols2);
       /* except for equalities,
	  error must be always negative, the real error is the positive part */
          for (Index i=0;i<nbRowsInBasis;i++) {
              int r = wholeBasis[i];
              if (cststat[r][INACTIVE] || cststat[r][REMOVED]) continue;
              if ((*Afacets)[r]->second.is_eqcst()) { 
		  errorInBasis[i] = errorPart[r]; continue; 
	      }
              errorInBasis[i]=max(errorPart[r],Interval::zero());
			 /* negative is ok */
          }
          for (Index i=nbRowsInBasis;i<nbCols;i++) {
             int c = wholeBasis[i]-nbRows;
             /* note : primalSol is punctual */
             if (primalSol[c].ub()>bbox[c].ub()) 
		errorInBasis[i] = primalSol[c]-bbox[c].ub();
             else if (primalSol[c].lb()<bbox[c].lb()) 
		errorInBasis[i] = primalSol[c]-bbox[c].lb();
          }
       /* compute the correction */
          IntervalVector correction=basisInverse*errorInBasis;
          primalSol -= correction;
          for (Index i=0;i<nbRows;i++) {
             errorPart[i] = (*Afacets)[i]->first.get_row()*primalSol 
			- (*Afacets)[i]->second.get_rhs();
          } 
       }

       /* check the different constraints 
	  if (!checkempty),
	  do not check for the basis, as it will be ~0 but is really 0
          by construction (which is why we do not use primalSol.mid()
          or equivalent : we know that there exists a point in primalSol
          where the basis is 0, but that's all */
#ifdef CODAC2_EXTENSION_POLYTOPE_CLP_DEBUG
       std::cout << "primalSol with correction : " << primalSol << "\n";
       std::cout << "errorPart : " << errorPart << "\n";
#endif
       Interval mxVal(-1.0);
       for (Index i=0;i<nbCols;i++) {
          if (!checkempty && rowBasis[i+nbRows]) continue;
          if (primalSol[i].is_subset(bbox[i])) continue;
          if (primalSol[i].is_disjoint(bbox[i])) { mxVal=1.0; break; }
          mxVal=Interval(-1.0,1.0);
       }
       if (mxVal.lb()<=0.0) {
         for (Index i=0;i<nbRows;i++) {
             if (cststat[i][INACTIVE] || cststat[i][REMOVED]) continue;
             if (!checkempty && rowBasis[i]) continue;
             if ((*Afacets)[i]->second.is_eqcst()) errorPart[i]=abs(errorPart[i]);
             mxVal=max(errorPart[i],mxVal);
             if (mxVal.lb()>0.0) break;
          }
       }
       if (mxVal.ub()<=0.0) {
           status[NOTEMPTY]=true;
           Valobj = max(Valobj,objvect.dot(primalSol));
       } else if (mxVal.lb()<=0.0) {
           if (!checkempty) {
              /* change of approach, we slightly move 
                 primalSol and hope to prove the non-emptiness
                 we use another vector to save primalSol and get
		 it back if needed */
              IntervalVector primalSol2=primalSol & bbox;
              bool ok=true;
              if (!primalSol2.is_empty()) {
                 primalSol2=primalSol2.mid();
                 for (Index i=0;i<nbRows;i++) {
                    if ((*Afacets)[i]->second.is_eqcst()) continue; 	
			/* we should not even try, in this case */
                    Interval err = (*Afacets)[i]->first.get_row()*primalSol2 
 			- (*Afacets)[i]->second.get_rhs();
                    if (err.ub()>=0.0) 
  		       primalSol2 -= (1.1*err.ub()/(*Afacets)[i]->first.get_row().squaredNorm())*(*Afacets)[i]->first.get_row();
                 }
                 primalSol2=primalSol2 & bbox;
                 if (!primalSol2.is_empty()) {
                   primalSol2=primalSol2.mid();
                   for (Index i=0;i<nbRows;i++) {
                     errorPart[i] = (*Afacets)[i]->first.get_row()*primalSol2 
		 	- (*Afacets)[i]->second.get_rhs();
                   }   
#ifdef CODAC2_EXTENSION_POLYTOPE_CLP_DEBUG
                   std::cout << "new primalSol : " << primalSol2 << "\n";
                   std::cout << "errorPart 2 : " << errorPart << "\n";
#endif
                   for (Index i=0;i<nbRows;i++) {
                      if (cststat[i][INACTIVE] || cststat[i][REMOVED]) continue;
                      if (errorPart[i].ub()>0.0) {
                        ok=false; /* failed */
                        break;
                      }
                      if ((*Afacets)[i]->second.is_eqcst() && errorPart[i].lb()<0.0) {
                        ok=false; /* failed */
                        break;
                      }
                   }
                 } else ok=false;
              } else ok=false;
              if (ok) {
                status[NOTEMPTY]=true;
                primalSol = primalSol2;
                Valobj = max(Valobj,objvect*primalSol);
               
              } else {
                status[NOTEMPTY_APPROX]=true;
              }
           } else {
              status[NOTEMPTY_APPROX]=true;
           }
       } else {
	   status[ERROR_PRIMAL_CHECK]=true;
       }
   } 
   if (stat==2) { 
	 /* check unboundedness. This is the last step 
            (alternative is stat==0) */
         if (!checkempty) { /* needs to build the ray */
            double *ray = model->unboundedRay();
            assert_release(ray!=nullptr);
            for (Index i=0;i<nbCols;i++) primalRay[i]=ray[i];
            delete(ray);
         }
         IntervalVector errorPart(nbRows);
         for (Index i=0;i<nbRows;i++) {
            errorPart[i] = (*Afacets)[i]->first.get_row()*primalRay;
         }  
         if (!checkempty) {
	    /* if !checkempty, the ray "follows" some facets,
               so we need to guarantee that it will not cross it */
            IntervalVector errorInBasis= IntervalVector::Zero(nbCols);
            for (Index i=0;i<nbRowsInBasis;i++) {
                int r = wholeBasis[i];
                if (cststat[r][INACTIVE] || cststat[r][REMOVED]) continue;
                if ((*Afacets)[r]->second.is_eqcst())
	             { errorInBasis[i] = errorPart[r]; continue; }
                errorInBasis[i]=max(errorPart[r],Interval::zero());
			 /* negative is ok */
            }
            for (Index i=nbRowsInBasis;i<nbCols;i++) {
                int r = wholeBasis[i];
                int c = r -nbRows;
                if (bbox[c].ub()<+oo) 
                    if (bbox[c].lb()>-oo) 
                      errorInBasis[i]=primalRay[r];
                    else
                      errorInBasis[i]=max(primalRay[r],Interval::zero());
                else 
		if (bbox[c].lb()>-oo) 
                   errorInBasis[i]=min(primalRay[r],Interval::zero());
            }
            /* compute the correction */
            IntervalVector correction=basisInverse*errorInBasis;
            primalRay -= correction;
            for (Index i=0;i<nbRows;i++) {
               errorPart[i] = (*Afacets)[i]->first.get_row()*primalRay;
            }  
         }
         Interval mxVal(-1.0);
         for (Index i=0;i<nbCols;i++) {
             if (!checkempty && rowBasis[i+nbRows]) continue;
             if (bbox[i].ub()<+oo) mxVal = max(primalSol[i],mxVal);
             if (bbox[i].lb()>-oo) mxVal = max(-primalSol[i],mxVal);
             if (mxVal.lb()>0.0) break;
         }
         if (mxVal.lb()<=0.0) {
            for (Index i=0;i<nbRows;i++) {
               if (cststat[i][INACTIVE] || cststat[i][REMOVED]) continue;
               if (!checkempty && rowBasis[i]) continue;
               if ((*Afacets)[i]->second.is_eqcst()) errorPart[i]=abs(errorPart[i]);
               mxVal=max(errorPart[i],mxVal);
               if (mxVal.lb()>0.0) break;
            }
         }
         if (mxVal.lb()<=0.0) {
              if (mxVal.ub()>0.0 && !checkempty) {
                 /* change of approach, we slightly move 
                    the value and hope to prove the validity
                    we use another vector to save primalSol and get
		    it back if needed */
		 /* for bbox, we create a IntervalVector for possible 
			primalRay */
                  IntervalVector admissibleRay(nbCols);
                  for (Index i=0;i<nbCols;i++) {
                     if (bbox[i].ub()<+oo) 
			admissibleRay[i] = Interval(0,+oo);
                     if (bbox[i].lb()>-oo) 
			admissibleRay[i] &= Interval(-oo,0);
                  }
                  IntervalVector primalRay2=primalRay & admissibleRay;
                  if (!primalRay2.is_empty()) {
                     primalRay2 = primalRay2.mid();
                     for (Index i=0;i<nbRows;i++) {
                        if ((*Afacets)[i]->second.is_eqcst()) continue; 	
			/* we should not even try, in this case */
                        Interval err = (*Afacets)[i]->first.get_row()*primalRay2;
                        if (err.ub()>=0.0) {
		            primalRay2 -= (1.1*err.ub()/(*Afacets)[i]->first.get_row().squaredNorm())*(*Afacets)[i]->first.get_row();
                        }
                    }
                    primalRay2=primalRay2 & admissibleRay;
                    if (!primalRay2.is_empty()) {
                      primalRay2.mid();
                      for (Index i=0;i<nbRows;i++) {
                        errorPart[i] = (*Afacets)[i]->first.get_row()*primalRay2;
                      }   
#ifdef CODAC2_EXTENSION_POLYTOPE_CLP_DEBUG
                      std::cout << "new primalRay : " << primalRay2 << "\n";
                      std::cout << "errorPart 2 : " << errorPart << "\n";
#endif
                      bool ok=true;
                      for (Index i=0;i<nbRows;i++) {
                         if (cststat[i][INACTIVE] || cststat[i][REMOVED]) continue;
                         if (errorPart[i].ub()>0.0) {
                           ok=false; /* failed */
                           break;
                         }
                      }
                      if (ok) {
                        Interval evolobj = 
			    checkempty ? Interval(1) : objvect.dot(primalRay2);
                        if (evolobj.lb()>0.0) {
                           status[UNBOUNDED]=true;
                           primalRay = primalRay2;
                           return status;
		        }
                     }
		  }
                }
             }
             Interval evolobj = checkempty ? Interval(1) :
			objvect.dot(primalRay);
             if (mxVal.ub()<=0.0 && evolobj.lb()>0.0)
	        status[UNBOUNDED]=true;
	     else if (evolobj.ub()>0.0)
	        status[UNBOUNDED_APPROX]=true;
             else status[ERROR_DUAL_CHECK]=true;
          } else {
             status[ERROR_DUAL_CHECK]=true;
          }  
          return status;
   }
   {
      /* finally, check boundedness */
             /* compute the error term , and reverse it with basis. */
             /* val : csts => val * Amat : variables => 
		val*Amat*basisInverse => basisInverse.cols */
      /* if checkempty, the goal may be to prove the global
         boundedness of the polytope */
      BoolInterval ok;
      IntervalRow errorNeum;
      if (checkempty) {
          errorNeum=IntervalRow::Zero(nbCols2);
          errorNeum[nbCols]= -1.0;
      } else {
          errorNeum=-objvect;
      }
      ok  = this->checkDualForm(dualSol, errorNeum,
		basisInverse, wholeBasis, nbRowsInBasis, checkempty);
      if (ok==BoolInterval::TRUE) {
         status[BOUNDED]=true;
	 Valobj = min(Valobj,this->dotprodrhs(dualSol));
         return status;
      } 
      /* Neumaier's approach , only with !checkempty */
      if (!checkempty && !bbox.is_unbounded()) {
          /* rebuild the solution */
          for (Index i=0;i<nbRows;i++) {
             dualSol[i]=model->dualRowSolution()[i];
             if (!(*Afacets)[i]->second.is_eqcst())
			 assert_release(model->dualRowSolution()[i]>=0.0);
          }
          for (Index i=0;i<nbCols;i++) {
             dualSol[i+nbRows]=model->dualColumnSolution()[i];
          }
	  Interval productNeum = this->dotprodrhs(dualSol);
          productNeum -= errorNeum.head(nbCols).dot(bbox);
          status[BOUNDED_BBOX]=true;
          status[BOUNDED]=true;
	  Valobj = min(Valobj,productNeum);
          return status;
      } else if (!checkempty) {
	/* other approach : we try to compensate with other values */ 
          /* for each epsilon value in dualSol, try to find a complement
             in comp */
          IntervalMatrix comp(nbRows,nbCols);
          {
            for (Index j=0;j<nbRows;j++) {
              comp.row(j)=((IntervalRow)((*Afacets)[j]->first.get_row()))*
			basisInverse;
            }
          }
          std::vector<Interval> okrow(nbRows,Interval(0.0,oo));
          int nbok=nbRows;
          for (Index j=0;j<nbRows;j++) {
             if (cststat[j][INACTIVE] || cststat[j][REMOVED]) {
			okrow[j].set_empty(); nbok--; continue; 
             }
             if (rowBasis[j]) { okrow[j].set_empty(); nbok--; continue; }
             if ((*Afacets)[j]->second.is_eqcst()) { okrow[j].init(); }
          }
          for (Index i=0;i<nbRowsInBasis;i++) {
             if (nbok==0) break;
             int r = wholeBasis[i];
             for (Index j=0;j<nbRows;j++) {
                if (okrow[j].is_empty()) continue;
                if (comp(j,i).ub()>0.0) { 
                   if (dualSol[r].lb()<=0.0) {
                        okrow[j].set_empty(); nbok--; continue;
                   } else {
			Interval u(dualSol[r].lb());
			u /= comp(j,i).ub();
                        okrow[j] &= Interval((*Afacets)[j]->second.is_eqcst() ? -oo : 0.0,
				u.lb());
			if (okrow[j].is_empty()) { nbok--; continue; }
                   }
                } else if (comp(j,i).ub()==0.0) {
                   if (dualSol[r].lb()<0.0) {
                        okrow[j].set_empty(); nbok--; 
                   } 
	           continue;
                } else {
                   if (dualSol[r].lb()>=0.0) continue;
	   	   Interval u(dualSol[r].lb());
		   u /= comp(j,i).ub();
		   okrow[j] &= Interval(u.ub(),oo);
                   if (okrow[j].is_empty()) { nbok--; continue; }
                }
             }
             dualSol[r]=model->dualRowSolution()[r];
          }
          /* recompute errorNeum */
          if (nbok>0) {
             for (Index i=0;i<nbRows;i++) {
                dualSol[i]=model->dualRowSolution()[i];
             }
             for (Index i=0;i<nbCols;i++) {
                dualSol[i+nbRows]=model->dualColumnSolution()[i];
             }
             for (Index j=0;j<nbRows;j++) {
                if (okrow[j].is_empty()) continue;
#ifdef CODAC2_EXTENSION_POLYTOPE_CLP_DEBUG
		std::cout << " correct dualsol : " << okrow[j] << "\n";
#endif
                dualSol[j]=okrow[j].lb(); /* in theory, lb should
			be sufficient. In practice... that's a problem */
                if (dualSol[j].ub()>okrow[j].ub()) dualSol[j]=okrow[j].mid();
		break;
             }
             errorNeum=-objvect;
             for (Index i=0;i<nbRows;i++) {
                errorNeum += dualSol[i]*(*Afacets)[i]->first.get_row();
             }
             IntervalRow error = errorNeum*basisInverse;
             ok=BoolInterval::TRUE;
             for (Index i=0;i<nbRowsInBasis;i++) {
                int r = wholeBasis[i];
                dualSol[r] -= error[i];
 		if (!(*Afacets)[r]->second.is_eqcst())
				 dualRay[r] &= Interval(0.0,oo);
//                if (dualSol[r].ub()<0.0) { ok = BoolInterval::FALSE; break; }
//                if (dualSol[r].lb()<0.0) { ok = BoolInterval::UNKNOWN; break; }
             }
             /* on refait les vérification sur les variables colonnes */
             if ((ok & BoolInterval::TRUE)==BoolInterval::TRUE) { 
                 /* vérifie les variables  colonnes, qui doivent
                  être libre, ainsi que les contraintes inactives */
                for (Index i=nbRowsInBasis;i<nbCols2;i++) {
                    if (error[i].mig()==0.0) continue;
                    int r = wholeBasis[i];
                    dualSol[r] -= error[i];
                    int c = r - nbRows;
                    if (bbox[c].ub()==+oo) {
                      if (dualSol[r].lb()>0.0) 
			{ ok=BoolInterval::FALSE; break; }
                      if (dualSol[r].ub()>0.0) 
			{ ok = BoolInterval::UNKNOWN; }
                    }
                    if (bbox[c].lb()==-oo) {
                      if (dualSol[r].ub()<0.0) 
			{ ok= BoolInterval::FALSE; break; }
                      if (dualSol[r].lb()<0.0) 
		        { ok = BoolInterval::UNKNOWN; }
     	            }
                }
              } 
              if ((ok & BoolInterval::TRUE)==BoolInterval::TRUE) {
                for (Index i=0;i<nbRows;i++) {
	           if ((cststat[i][INACTIVE] || cststat[i][REMOVED]) &&
		       (dualSol[i].mig()>0.0))
				 { ok=BoolInterval::FALSE; break; }
                } 
             }
          }
          if (ok==BoolInterval::TRUE) { 
             status[BOUNDED]=true;
	     Valobj = min(Valobj,this->dotprodrhs(dualSol));
             return status;
          } 
      }
      if (ok==BoolInterval::UNKNOWN) {
         status[BOUNDED_APPROX]=true;
      } else {
         status[ERROR_DUAL_CHECK]=true;
      }
   }  
   return status;
}   

/** simplified writing of lp_result_stat for debugging **/
std::ostream& print_lp_stat(std::ostream& os, const LPclp::lp_result_stat x) {
    if (x[LPclp::CHANGED]) {
        os << "Not_computed"; return os;
    }
    if (x[LPclp::ERROR_LPCOIN]) {
        os << "CLP_failed"; return os;
    }
    /* primal */
    if (x[LPclp::EMPTY]) {
        os << "Empty"; 
        if (x[LPclp::EMPTY_BBOX]) os << "(bbox)";
        return os; 
    } else if (x[LPclp::EMPTY_APPROX]) {
        os << "Empty(app)-";
    } else if (x[LPclp::NOTEMPTY]) {
        os << "NonEmpty-";
    } else if (x[LPclp::NOTEMPTY_APPROX]) {
        os << "NonEmpty(app)-";
    } else if (x[LPclp::ERROR_PRIMAL_CHECK]) {
        os << "Error(prim)";
    } else os << "?(prim)";
    if (x[LPclp::UNBOUNDED]) {
        os << "Unbounded"; 
    } else if (x[LPclp::UNBOUNDED_APPROX]) {
        os << "Unbounded(app)";
    } else if (x[LPclp::BOUNDED]) {
        os << "Bounded";
        if (x[LPclp::BOUNDED_BBOX]) os << "(bbox)";
    } else if (x[LPclp::BOUNDED_APPROX]) {
        os << "Bounded(app)";
    } else if (x[LPclp::ERROR_DUAL_CHECK]) {
        os << "Error(dual)";
    } else os << "?(dual)";
    return os;
}

}

