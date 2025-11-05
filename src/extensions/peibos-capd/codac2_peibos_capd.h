/** 
 *  \file codac2_peibos_capd.h
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Maël Godard
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include "codac2_peibos.h"
#include <capd/capdlib.h>
#include "codac2_capd.h"
// #include <capd/poincare/lib.h>
#include <vector>

namespace codac2
{
  std::vector<Parallelepiped> PEIBOS(const capd::IMap& i_map, double tf, const AnalyticFunction<VectorType>& psi_0, const std::vector<OctaSym>& Sigma, double epsilon, bool verbose = false);
  std::vector<Parallelepiped> PEIBOS(const capd::IMap& i_map, double tf, const AnalyticFunction<VectorType>& psi_0, const std::vector<OctaSym>& Sigma, double epsilon, const Vector& offset, bool verbose = false);
}