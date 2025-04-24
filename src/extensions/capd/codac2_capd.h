/** 
 *  \file codac2_capd.h
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Maël Godard
 *  \copyright  Copyright 2024 CODAC Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <capd/intervals/lib.h>
#include <capd/vectalg/lib.h>
#include "codac2_IntervalVector.h"
#include "codac2_Matrix.h"
#include "codac2_IntervalMatrix.h"

namespace codac2
{
  /**
   * \brief Converts a CODAC Interval object into an CAPD Interval object
   * 
   * \param x the CODAC Interval to be converted
   * \return the converted CAPD Interval
   */
  capd::Interval to_capd(const codac2::Interval& x);

  /**
   * \brief Converts an CAPD Interval object into a CODAC Interval object
   * 
   * \param x the CAPD Interval to be converted
   * \return the converted CODAC Interval
   */
  codac2::Interval to_codac(const capd::Interval& x);

  /**
   * \brief Converts a CODAC IntervalVector object into an CAPD IVector object
   * 
   * \param x the CODAC IntervalVector to be converted
   * \return the converted CAPD IVector
   */
  capd::IVector to_capd(const codac2::IntervalVector& x);

  /**
   * \brief Converts an CAPD IVector object into a CODAC IntervalVector object
   * 
   * \param x the CAPD IVector to be converted
   * \return the converted CODAC IntervalVector
   */
  codac2::IntervalVector to_codac(const capd::IVector& x);

  /**
   * \brief Converts a CODAC IntervalMatrix object into an CAPD IMatrix object
   * 
   * \param x the CODAC IntervalMatrix to be converted
   * \return the converted CAPD IMatrix
   */
  capd::IMatrix to_capd(const codac2::IntervalMatrix& x);

  /**
   * \brief Converts an CAPD IMatrix object into a CODAC IntervalMatrix object
   * 
   * \param x the CAPD IMatrix to be converted
   * \return the converted CODAC IntervalMatrix
   */
  codac2::IntervalMatrix to_codac(const capd::IMatrix& x);

}