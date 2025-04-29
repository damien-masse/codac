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
#include "codac2_clp.h"

#undef CODAC2_EXTENSION_POLYTOPE_CLP_DEBUG

namespace codac2 {

/**********  Problem building funtions ***********/

LPclp::LPclp(const Matrix &mat, const Row &objvect,
                const Vector &rhsvect) : 
   nbRows(mat.rows()), nbCols(mat.cols()), 
   Amat(mat), 
   objvect(objvect),
   rhsvect(rhsvect), cststat(nbRows,0), 
   rowBasis(nbRows,false), built(false), model(nullptr),
   primalSol(nbCols), dualSol(nbRows), primalRay(nbCols), dualRay(nbRows),
   bbox(nbCols)
{
}

LPclp::LPclp(const Matrix &mat, 
                const Vector &rhsvect) : 
   nbRows(mat.rows()), nbCols(mat.cols()), 
   Amat(mat),
   objvect(Row::Zero(nbCols)),
   rhsvect(rhsvect), cststat(nbRows,0), 
   rowBasis(nbRows,false), built(false), model(nullptr),
   primalSol(nbCols), dualSol(nbRows), primalRay(nbCols), dualRay(nbRows),
   bbox(nbCols)
{
   
}

LPclp::LPclp(Index dim, Index nbcsts) :
   nbRows(nbcsts), nbCols(dim), 
   Amat(Matrix::Zero(nbRows,nbCols)), 
   objvect(Row::Zero(nbCols)),
   rhsvect(Vector::Zero(nbRows)),
   cststat(nbRows,1<<REMOVED),
   rowBasis(nbRows,false), built(false), model(nullptr),
   primalSol(nbCols), dualSol(nbRows), primalRay(nbCols), dualRay(nbRows),
   bbox(nbCols)
{
}

LPclp::~LPclp() {
   delete(model);
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
       if (built) model->setRowUpper(cst,rhsvect[cst]);
       this->status=1<<CHANGED;
    } else {
       if (cststat[cst][INACTIVE]) return;
       cststat[cst][INACTIVE]=true;
       if (built) model->setRowUpper(cst,DBL_MAX);
       this->status=1<<CHANGED;
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
   for (Index i=0;i<nbCols2;i++) {
      model->setColumnBounds(i,-COIN_DBL_MAX, COIN_DBL_MAX);
   }

   /* adding the rows using BuildObject ? */
   CoinBuild buildObject;
   int* row2Index = new int[nbCols2];
   for(Index i=0;i<nbCols2;i++){row2Index[i]=i;}
   double* row2Vals = new double[nbCols2];
   for(Index k=0; k<nbRows; k++) {
         double bnd = ((cststat[k][REMOVED] || cststat[k][INACTIVE]
			|| cststat[k][REDUNDANT]) ? 
			COIN_DBL_MAX : rhsvect[k]);
         
         for (Index j=0;j<nbCols;j++) { 
		/* FIXME : frustrating approach for the model,
                   but Amat is not in "correct" storage order.
		   (or we consider the transpose order ? */
           row2Vals[j]=Amat(k,j);
         }
         if (emptytest) { row2Vals[nbCols]=1.0; }
         buildObject.addRow(nbCols2, row2Index, row2Vals,
                     -COIN_DBL_MAX, bnd);
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
   built_emptytest=emptytest;
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
      for(Index i=0;i<nbRows;i++){col2Vals[i]=1.0;}
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
//   rowBasis.reset();  
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


int LPclp::minimize_polytope(const  Interval &tolerance, bool checkempty) {
    if (checkempty) {
       lp_result_stat stat=this->solve(true);
       if (stat[EMPTY]) return -1;
    }
    bool solvedOnce=checkempty;
    int nb=0;
    for (Index i=nbRows-1;i>=0;i--) {
       if (cststat[i][REMOVED] || 
	   cststat[i][INACTIVE] || 
	   cststat[i][REDUNDANT])
           continue;
       this->setActive(i,false);
       this->setObjective(this->getCst(i));
       lp_result_stat stat=this->solve(false,(solvedOnce ? 4 : 0));
       if (stat[EMPTY]) return -1;
       if (stat[BOUNDED]) {
           if (this->Valobj.ub()<=rhsvect[i]+tolerance.ub() 
		&& this->Valobj.lb()<=rhsvect[i]+tolerance.lb()) {
              cststat[i][REDUNDANT]=true;
              continue;
           }
       }
       nb++;
       this->setActive(i,true);
       solvedOnce=true;
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
        } else {
           rowBasis[i]=false;
        }
        dualSol[i]=model->dualRowSolution()[i];
    }
    /* col basics */
    double slack=0.0; /* only for emptiness checking */
    for (Index i=0;i<nbCols2;i++) {
	if (model->getColumnStatus(i)!=ClpSimplex::Status::basic)  
	 	wholeBasis.push_back(i+nbRows);
//      dualVals[i+nbRows]=model->dualColumnSolution()[i]; /* useless for us*/
        if (i<nbCols) primalSol[i]=model->primalColumnSolution()[i];
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
          basis(i,c)=Amat(r,c);
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
           double *ray = model->infeasibilityRay(false);
           assert_release(ray!=nullptr);  /* false => bug of CLP */
             /* true for infeasibilityRay means that we have also
                the negation of the "error" term after the rows,
                but we will compute it ourselves.
                unless the variables are bounded, this is not needed */
           for (Index i=0;i<nbRows;i++) { dualRay[i]=ray[i]; }
             /* compute the error term , and reverse it with basis. */
             /* val : csts => val * Amat : variables (~0)  
	 	val*Amat*basisInverse => basisInverse.cols */
         }
#ifdef CODAC2_EXTENSION_POLYTOPE_CLP_DEBUG
         std::cout << "dualRay: " << dualRay << "\n";
#endif
         IntervalRow errorNeum=IntervalRow::Zero(nbCols2);
         errorNeum.head(nbCols)=dualRay*Amat;
         if (checkempty) errorNeum[nbCols]=dualRay.sum()-1.0;
         IntervalRow error = errorNeum*basisInverse;
         BoolInterval ok=BoolInterval::TRUE;
         for (Index i=0;i<nbRowsInBasis;i++) {
             int r = wholeBasis[i];
             dualRay[r] -= error[i]; /* modify the ray with the error,
			    the value must stay positive */
             if (dualRay[r].ub()<0.0) { ok = BoolInterval::FALSE; break; }
             if (dualRay[r].lb()<0.0) { ok = BoolInterval::UNKNOWN; }
         }
#ifdef CODAC2_EXTENSION_POLYTOPE_CLP_DEBUG
         std::cout << "dualRay with error: " << dualRay << "\n";
         std::cout << "ok: " << ok << "\n";
#endif
         if ((ok & BoolInterval::TRUE)==BoolInterval::TRUE)
	 { /* vérifie les variables  colonnes, qui doivent
		  être libre, ainsi que les contraintes inactives */
            for (Index i=nbRowsInBasis;i<nbCols2;i++) {
              if (error[i].mig()>0.0) {
 	         ok=BoolInterval::FALSE;
	         break;
              }
            } 
            if ((ok & BoolInterval::TRUE)==BoolInterval::TRUE) 
               for (Index i=0;i<nbRows;i++) {
	          if ((cststat[i][INACTIVE] || cststat[i][REMOVED]) &&
		     (dualRay[i].mig()>0.0)) { ok=BoolInterval::FALSE; break; }
               }
	}
        if ((ok & BoolInterval::TRUE)==BoolInterval::TRUE) {
          Interval product = dualRay.dot(rhsvect);
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
        if (!bbox.is_unbounded()) {
            /* restart with the ray */
            for (Index i=0;i<nbRows;i++) { 
                assert_release(ray[i]<0.0);
		dualRay[i]=ray[i]; 
            }
	    Interval productNeum = dualRay.dot(rhsvect);
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
          IntervalMatrix comp = ((IntervalMatrix)(Amat))*basisInverse.topRows(nbCols);
          if (checkempty) {
             IntervalVector v = IntervalVector::Ones(nbRows);
             comp += v*basisInverse.row(nbCols);
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
                        okrow[j] &= Interval(0.0,
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
             for (Index i=0;i<nbRows;i++) {
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
                if (dualRay[j].ub()>okrow[j].ub()) dualRay[j]=okrow[j].mid();
		break;
             }
             errorNeum.head(nbCols)=dualRay*Amat;
             if (checkempty) errorNeum[nbCols]=dualRay.sum()-1.0;
             error = errorNeum*basisInverse;
             ok=BoolInterval::TRUE;
             for (Index i=0;i<nbRowsInBasis;i++) {
                int r = wholeBasis[i];
                dualRay[r] -= error[i];
 		dualRay[r] &= Interval(0.0,oo);
//                if (dualSol[r].ub()<0.0) { ok = BoolInterval::FALSE; break; }
//                if (dualSol[r].lb()<0.0) { ok = BoolInterval::UNKNOWN; break; }
             }
             /* on refait les vérification sur les variables colonnes */
             if ((ok & BoolInterval::TRUE)==BoolInterval::TRUE) { 
                 /* vérifie les variables  colonnes, qui doivent
                  être libre, ainsi que les contraintes inactives */
                for (Index i=nbRowsInBasis;i<nbCols2;i++) {
                  if (error[i].mig()>0.0) {
 	            ok=BoolInterval::FALSE;
	            break;
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
                Interval product = dualRay.dot(rhsvect);
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
           Interval product = dSol.dot(rhsvect);
           dualSol.setZero();
           bool ok=true;
           for (Index i=0;i<nbRowsInBasis;i++) {
               if (dSol[i].lb()<0.0) { ok=false; break; }
               dualSol[wholeBasis[i]]=dSol[i];
           }
           if (ok) {
              for (Index i=nbRowsInBasis;i<nbCols;i++) {
                 if (dSol[i].lb()<0.0 || dualSol[i].ub()>0.0) {
                    if (bbox.is_unbounded()) { ok=false; break; }
                    status[BOUNDED_BBOX]=true;
                    product -= dSol[i].lb() * 
				bbox[wholeBasis[i]-nbRows];
                 }
              }
           }
           if (ok) {
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
       IntervalVector errorPart = Amat*primalSol - rhsvect;
       if (!checkempty) {
          IntervalVector errorInBasis= IntervalVector::Zero(nbCols2);
       /* error must be always negative, the real error is the positive part */
          for (Index i=0;i<nbRowsInBasis;i++) {
              int r = wholeBasis[i];
              if (cststat[r][INACTIVE] || cststat[r][REMOVED]) continue;
              errorInBasis[i]=max(errorPart[r],Interval::zero());
			 /* negative is ok */
          }
       /* compute the correction */
          IntervalVector correction=basisInverse*errorInBasis;
          primalSol -= correction;
          errorPart = Amat*primalSol - rhsvect;
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
       for (Index i=0;i<nbRows;i++) {
           if (cststat[i][INACTIVE] || cststat[i][REMOVED]) continue;
           if (!checkempty && rowBasis[i]) continue;
#ifdef CODAC2_EXTENSION_POLYTOPE_CLP_DEBUG
           std::cout << "check " << i << "\n";
#endif
           mxVal=max(errorPart[i],mxVal);
           if (mxVal.lb()>0.0) break;
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
              IntervalVector primalSol2=primalSol.mid();
              for (Index i=0;i<nbRows;i++) {
                  Interval err = Amat.row(i)*primalSol2 - rhsvect[i];
                  if (err.ub()>=0.0) 
		     primalSol2 -= (1.1*err.ub()/Amat.row(i).squaredNorm())*Amat.row(i);
              }
              primalSol2=primalSol2.mid();
              errorPart=Amat*primalSol2 - rhsvect;
#ifdef CODAC2_EXTENSION_POLYTOPE_CLP_DEBUG
              std::cout << "new primalSol : " << primalSol2 << "\n";
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
         IntervalVector errorPart = Amat*primalRay;
         if (!checkempty) {
	    /* if !checkempty, the ray "follows" some facets,
               so we need to guarantee that it will not cross it */
            IntervalVector errorInBasis= IntervalVector::Zero(nbCols);
            for (Index i=0;i<nbRowsInBasis;i++) {
                int r = wholeBasis[i];
                if (cststat[r][INACTIVE] || cststat[r][REMOVED]) continue;
                errorInBasis[i]=max(errorPart[r],Interval::zero());
			 /* negative is ok */
            }
            /* compute the correction */
            IntervalVector correction=basisInverse*errorInBasis;
            primalRay -= correction;
            errorPart = Amat*primalRay;
         }
         Interval mxVal(-1.0);
         for (Index i=0;i<nbRows;i++) {
             if (cststat[i][INACTIVE] || cststat[i][REMOVED]) continue;
             if (!checkempty && rowBasis[i]) continue;
             mxVal=max(errorPart[i],mxVal);
             if (mxVal.lb()>0.0) break;
          }
          if (mxVal.lb()<=0.0) {
              if (mxVal.ub()>0.0 && !checkempty) {
                 /* change of approach, we slightly move 
                    the value and hope to prove the validity
                    we use another vector to save primalSol and get
		    it back if needed */
                  IntervalVector primalRay2=primalRay.mid();
                  for (Index i=0;i<nbRows;i++) {
                      Interval err = Amat.row(i)*primalRay2;
                      if (err.ub()>=0.0) {
		          primalRay2 -= (1.1*err.ub()/Amat.row(i).squaredNorm())*Amat.row(i);
                      }
                  }
                  primalRay2=primalRay2.mid();
                  errorPart=Amat*primalRay2;
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
                     Interval evolobj = objvect.dot(primalRay2);
                     if (evolobj.lb()>0.0) {
                        status[UNBOUNDED]=true;
                        primalRay = primalRay2;
                        return status;
		     }
                  }
             }
             Interval evolobj = objvect.dot(primalRay);
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
      IntervalRow errorNeum=IntervalRow::Zero(nbCols2);
      if (checkempty) {
         errorNeum.head(nbCols)=dualSol*Amat;
         errorNeum[nbCols]=dualSol.sum()-1.0;
      } else {
         errorNeum = dualSol*Amat-objvect;
      }
      IntervalRow error = errorNeum*basisInverse;
      BoolInterval ok=BoolInterval::TRUE;
      for (Index i=0;i<nbRowsInBasis;i++) {
          int r = wholeBasis[i];
          dualSol[r] -= error[i];
          if (dualSol[r].ub()<0.0) { ok = BoolInterval::FALSE; break; }
          if (dualSol[r].lb()<0.0) { ok = BoolInterval::UNKNOWN; }
      }
      if ((ok & BoolInterval::TRUE)==BoolInterval::TRUE) { 
            /* vérifie les variables  colonnes, qui doivent
               être libre, ainsi que les contraintes inactives */
            for (Index i=nbRowsInBasis;i<nbCols2;i++) {
               if (error[i].mig()>0.0) {
 	          ok=BoolInterval::FALSE;
	          break;
               }
            } 
            if ((ok & BoolInterval::TRUE)==BoolInterval::TRUE) 
              for (Index i=0;i<nbRows;i++) {
	         if ((cststat[i][INACTIVE] || cststat[i][REMOVED]) &&
		     (dualSol[i].mig()>0.0)) { ok=BoolInterval::FALSE; break; }
              }
      }
      if (ok==BoolInterval::TRUE) {
         status[BOUNDED]=true;
	 Valobj = min(Valobj,dualSol.dot(rhsvect));
         return status;
      } 
      /* Neumaier's approach , only with !checkempty */
      if (!bbox.is_unbounded()) {
          /* rebuild the solution */
          for (Index i=0;i<nbRows;i++) {
             dualSol[i]=model->dualRowSolution()[i];
          }
	  Interval productNeum = dualSol.dot(rhsvect);
          productNeum -= errorNeum.head(nbCols).dot(bbox);
          status[BOUNDED_BBOX]=true;
          status[BOUNDED]=true;
	  Valobj = min(Valobj,productNeum);
          return status;
      } else if (!checkempty) {
	/* other approach : we try to compensate with other values */ 
          /* for each epsilon value in dualSol, try to find a complement
             in comp */
          IntervalMatrix comp = ((IntervalMatrix)(Amat))*basisInverse;
          std::vector<Interval> okrow(nbRows,Interval(0.0,oo));
          int nbok=nbRows;
          for (Index j=0;j<nbRows;j++) {
             if (cststat[j][INACTIVE] || cststat[j][REMOVED]) {
			okrow[j].set_empty(); nbok--; continue; 
             }
             if (rowBasis[j]) { okrow[j].set_empty(); nbok--; continue; }
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
                        okrow[j] &= Interval(0.0,
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
             errorNeum = dualSol*Amat-objvect;
             error = errorNeum*basisInverse;
             ok=BoolInterval::TRUE;
             for (Index i=0;i<nbRowsInBasis;i++) {
                int r = wholeBasis[i];
                dualSol[r] -= error[i];
 		dualSol[r] &= Interval(0.0,oo);
//                if (dualSol[r].ub()<0.0) { ok = BoolInterval::FALSE; break; }
//                if (dualSol[r].lb()<0.0) { ok = BoolInterval::UNKNOWN; break; }
             }
             /* on refait les vérification sur les variables colonnes */
             if ((ok & BoolInterval::TRUE)==BoolInterval::TRUE) { 
                 /* vérifie les variables  colonnes, qui doivent
                  être libre, ainsi que les contraintes inactives */
                for (Index i=nbRowsInBasis;i<nbCols2;i++) {
                  if (error[i].mig()>0.0) {
 	            ok=BoolInterval::FALSE;
	            break;
                 }
              } 
              if ((ok & BoolInterval::TRUE)==BoolInterval::TRUE) 
                for (Index i=0;i<nbRows;i++) {
	           if ((cststat[i][INACTIVE] || cststat[i][REMOVED]) &&
		       (dualSol[i].mig()>0.0))
				 { ok=BoolInterval::FALSE; break; }
                } 
             }
             if (ok==BoolInterval::TRUE) { 
                status[BOUNDED]=true;
	        Valobj = min(Valobj,dualSol.dot(rhsvect));
                return status;
             } 
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

