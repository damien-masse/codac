/** 
 *  \file codac2_directed_ctc.h
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <string>
#include "codac2_AnalyticType.h"
#include "codac2_template_tools.h"
#include "codac2_IntervalRow.h"

namespace codac2
{
  struct AddOp
  {
    // Unary operations
    static Interval fwd(const Interval& x1);
    static ScalarType fwd(const ScalarType& x1);
    static void bwd(const Interval& y, Interval& x1);

    static IntervalVector fwd(const IntervalVector& x1);
    static VectorType fwd(const VectorType& x1);
    static void bwd(const IntervalVector& y, IntervalVector& x1);

    static IntervalMatrix fwd(const IntervalMatrix& x1);
    static MatrixType fwd(const MatrixType& x1);
    static void bwd(const IntervalMatrix& y, IntervalMatrix& x1);

    // Binary operations
    static Interval fwd(const Interval& x1, const Interval& x2);
    static ScalarType fwd(const ScalarType& x1, const ScalarType& x2);
    static void bwd(const Interval& y, Interval& x1, Interval& x2);

    static IntervalVector fwd(const IntervalVector& x1, const IntervalVector& x2);
    static VectorType fwd(const VectorType& x1, const VectorType& x2);
    static void bwd(const IntervalVector& y, IntervalVector& x1, IntervalVector& x2);

    static IntervalMatrix fwd(const IntervalMatrix& x1, const IntervalMatrix& x2);
    static MatrixType fwd(const MatrixType& x1, const MatrixType& x2);
    static void bwd(const IntervalMatrix& y, IntervalMatrix& x1, IntervalMatrix& x2);
  };

  struct SubOp
  {
    // Unary operations
    static Interval fwd(const Interval& x1);
    static ScalarType fwd(const ScalarType& x1);
    static void bwd(const Interval& y, Interval& x1);

    static IntervalVector fwd(const IntervalVector& x1);
    static VectorType fwd(const VectorType& x1);
    static void bwd(const IntervalVector& y, IntervalVector& x1);

    static IntervalMatrix fwd(const IntervalMatrix& x1);
    static MatrixType fwd(const MatrixType& x1);
    static void bwd(const IntervalMatrix& y, IntervalMatrix& x1);

    // Binary operations
    static Interval fwd(const Interval& x1, const Interval& x2);
    static ScalarType fwd(const ScalarType& x1, const ScalarType& x2);
    static void bwd(const Interval& y, Interval& x1, Interval& x2);

    static IntervalVector fwd(const IntervalVector& x1, const IntervalVector& x2);
    static VectorType fwd(const VectorType& x1, const VectorType& x2);
    static void bwd(const IntervalVector& y, IntervalVector& x1, IntervalVector& x2);

    static IntervalMatrix fwd(const IntervalMatrix& x1, const IntervalMatrix& x2);
    static MatrixType fwd(const MatrixType& x1, const MatrixType& x2);
    static void bwd(const IntervalMatrix& y, IntervalMatrix& x1, IntervalMatrix& x2);
  };

  struct MulOp
  {
    static Interval fwd(const Interval& x1, const Interval& x2);
    static ScalarType fwd(const ScalarType& x1, const ScalarType& x2);
    static void bwd(const Interval& y, Interval& x1, Interval& x2);

    static IntervalVector fwd(const Interval& x1, const IntervalVector& x2);
    static VectorType fwd(const ScalarType& x1, const VectorType& x2);
    static void bwd(const IntervalVector& y, Interval& x1, IntervalVector& x2);

    static IntervalVector fwd(const IntervalVector& x1, const Interval& x2);
    static VectorType fwd(const VectorType& x1, const ScalarType& x2);
    static void bwd(const IntervalVector& y, IntervalVector& x1, Interval& x2);

    static Interval fwd(const IntervalRow& x1, const IntervalVector& x2);
    //static ScalarType fwd(const RowType& x1, const VectorType& x2); // RowType not yet defined
    static void bwd(const Interval& y, IntervalRow& x1, IntervalVector& x2);

    static IntervalVector fwd(const IntervalMatrix& x1, const IntervalVector& x2);
    static VectorType fwd(const MatrixType& x1, const VectorType& x2);
    static void bwd(const IntervalVector& y, IntervalMatrix& x1, IntervalVector& x2);
  };

  struct DivOp
  {
    static Interval fwd(const Interval& x1, const Interval& x2);
    static ScalarType fwd(const ScalarType& x1, const ScalarType& x2);
    static void bwd(const Interval& y, Interval& x1, Interval& x2);

    static IntervalVector fwd(const IntervalVector& x1, const Interval& x2);
    static VectorType fwd(const VectorType& x1, const ScalarType& x2);
    static void bwd(const IntervalVector& y, IntervalVector& x1, Interval& x2);
  };

  struct PowOp
  {
    static Interval fwd(const Interval& x1, const Interval& x2);
    static ScalarType fwd(const ScalarType& x1, const ScalarType& x2);
    static void bwd(const Interval& y, Interval& x1, Interval& x2);
  };

  struct SqrOp
  {
    static Interval fwd(const Interval& x1);
    static ScalarType fwd(const ScalarType& x1);
    static void bwd(const Interval& y, Interval& x1);
  };

  struct SqrtOp
  {
    static Interval fwd(const Interval& x1);
    static ScalarType fwd(const ScalarType& x1);
    static void bwd(const Interval& y, Interval& x1);
  };

  struct ExpOp
  {
    static Interval fwd(const Interval& x1);
    static ScalarType fwd(const ScalarType& x1);
    static void bwd(const Interval& y, Interval& x1);
  };

  struct LogOp
  {
    static Interval fwd(const Interval& x1);
    static ScalarType fwd(const ScalarType& x1);
    static void bwd(const Interval& y, Interval& x1);
  };

  struct CosOp
  {
    static Interval fwd(const Interval& x1);
    static ScalarType fwd(const ScalarType& x1);
    static void bwd(const Interval& y, Interval& x1);
  };

  struct SinOp
  {
    static Interval fwd(const Interval& x1);
    static ScalarType fwd(const ScalarType& x1);
    static void bwd(const Interval& y, Interval& x1);
  };

  struct TanOp
  {
    static Interval fwd(const Interval& x1);
    static ScalarType fwd(const ScalarType& x1);
    static void bwd(const Interval& y, Interval& x1);
  };

  struct AcosOp
  {
    static Interval fwd(const Interval& x1);
    static ScalarType fwd(const ScalarType& x1);
    static void bwd(const Interval& y, Interval& x1);
  };

  struct AsinOp
  {
    static Interval fwd(const Interval& x1);
    static ScalarType fwd(const ScalarType& x1);
    static void bwd(const Interval& y, Interval& x1);
  };

  struct AtanOp
  {
    static Interval fwd(const Interval& x1);
    static ScalarType fwd(const ScalarType& x1);
    static void bwd(const Interval& y, Interval& x1);
  };

  struct Atan2Op
  {
    static Interval fwd(const Interval& x1, const Interval& x2);
    static ScalarType fwd(const ScalarType& x1, const ScalarType& x2);
    static void bwd(const Interval& y, Interval& x1, Interval& x2);
  };

  struct CoshOp
  {
    static Interval fwd(const Interval& x1);
    static ScalarType fwd(const ScalarType& x1);
    static void bwd(const Interval& y, Interval& x1);
  };

  struct SinhOp
  {
    static Interval fwd(const Interval& x1);
    static ScalarType fwd(const ScalarType& x1);
    static void bwd(const Interval& y, Interval& x1);
  };

  struct TanhOp
  {
    static Interval fwd(const Interval& x1);
    static ScalarType fwd(const ScalarType& x1);
    static void bwd(const Interval& y, Interval& x1);
  };

  struct AbsOp
  {
    static Interval fwd(const Interval& x1);
    static ScalarType fwd(const ScalarType& x1);
    static void bwd(const Interval& y, Interval& x1);
  };

  struct ComponentOp
  {
    static Interval fwd(const IntervalVector& x1, Index i);
    static ScalarType fwd(const VectorType& x1, Index i);
    static void bwd(const Interval& y, IntervalVector& x1, Index i);
  };

  struct SubvectorOp
  {
    static IntervalVector fwd(const IntervalVector& x1, Index i, Index j);
    static VectorType fwd(const VectorType& x1, Index i, Index j);
    static void bwd(const IntervalVector& y, IntervalVector& x1, Index i, Index j);
  };

  struct VectorOp
  {
    template<typename... X>
      requires (std::is_base_of_v<Interval,X> && ...)
    static IntervalVector fwd(const X&... x)
    {
      return IntervalVector({x...});
    }

    template<typename... X>
      requires (std::is_base_of_v<ScalarType,X> && ...)
    static VectorType fwd(const X&... x)
    {
      IntervalMatrix d(sizeof...(X),std::get<0>(std::tie(x...)).da.cols());
      Index i = 0;
      ((d.row(i++) = x.da), ...);

      bool def_domain = true;
      ((def_domain &= x.def_domain), ...);
      
      return {
        fwd(x.m...),
        fwd(x.a...),
        d,
        def_domain
      };
    }

    template<typename... X>
      requires (std::is_base_of_v<Interval,X> && ...)
    static void bwd(const IntervalVector& y, X&... x)
    {
      Index i = 0;
      ((x &= y[i++]), ...);
    }
  };

  struct MatrixOp
  {
    static void fwd_i(IntervalMatrix& m, const IntervalVector& x, Index i);

    template<typename... X>
      requires (std::is_base_of_v<Domain,X> && ...)
    static IntervalMatrix fwd(const X&... x)
    {
      throw std::runtime_error("MatrixOp not fully implemented yet");
      IntervalMatrix m(1, sizeof...(X));
      Index i = 0;
      (MatrixOp::fwd_i(m, x, i++), ...);
      return m;
    }

    template<typename... X>
      requires (std::is_base_of_v<VectorType,X> && ...)
    static MatrixType fwd(const X&... x)
    {
      throw std::runtime_error("MatrixOp not fully implemented yet");
      return {
        IntervalMatrix({x.m...}),
        IntervalMatrix({x.a...}),
        IntervalMatrix(0,0), // not supported yet for matrices
        true // todo with variadic
      };
    }

    template<typename... X>
      requires (std::is_base_of_v<IntervalVector,X> && ...)
    static void bwd(const IntervalMatrix& y, X&... x)
    {
      throw std::runtime_error("MatrixOp not fully implemented yet");
      Index i = 0;
      ((x &= y.col(i++)), ...);
    }
  };

  struct DetOp
  {
    // For matrices
    static Interval fwd(const IntervalMatrix& x);
    static ScalarType fwd(const MatrixType& x);
    static void bwd(const Interval& y, IntervalMatrix& x);

    // For two vectors (merged into a 2×2 matrix)
    static Interval fwd(const IntervalVector& x1, const IntervalVector& x2);
    static ScalarType fwd(const VectorType& x1, const VectorType& x2);
    static void bwd(const Interval& y, IntervalVector& x1, IntervalVector& x2);

    // For three vectors (merged into a 3×3 matrix)
    static Interval fwd(const IntervalVector& x1, const IntervalVector& x2, const IntervalVector& x3);
    static ScalarType fwd(const VectorType& x1, const VectorType& x2, const VectorType& x3);
    static void bwd(const Interval& y, IntervalVector& x1, IntervalVector& x2, IntervalVector& x3);
  };
}