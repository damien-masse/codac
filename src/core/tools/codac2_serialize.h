/** 
 *  \file codac2_serialize.h
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <fstream>
#include "codac2_Vector.h"
#include "codac2_SampledTraj.h"

namespace codac2
{
  // Trivial types (int, float, double...)

    template <typename T>
      requires (std::is_trivially_copyable_v<T>)
    void serialize(std::ofstream& f, const T& x)
    {
      f.write(reinterpret_cast<const char*>(&x), sizeof(T));
    }

    template <typename T>
      requires (std::is_trivially_copyable_v<T>)
    void deserialize(std::ifstream& f, T& x)
    {
      f.read(reinterpret_cast<char*>(&x), sizeof(T));
    }

  // Vector

    void serialize(std::ofstream& f, const Vector& x)
    {
      Index size = x.size();
      f.write(reinterpret_cast<const char*>(&size), sizeof(size));
      for(Index i = 0 ; i < size ; i++)
        serialize(f,x[i]);
    }

    void deserialize(std::ifstream& f, Vector& x)
    {
      Index size;
      f.read(reinterpret_cast<char*>(&size), sizeof(size));
      x.resize(size);
      for(Index i = 0 ; i < size ; i++)
        deserialize(f,x[i]);
    }

  // SampledTraj<T>

    template <typename T>
    void serialize(std::ofstream& f, const SampledTraj<T>& x)
    {
      Index size = x.nb_samples();
      f.write(reinterpret_cast<const char*>(&size), sizeof(size));

      for(const auto& [ti, xi] : x)
      {
        serialize(f, ti);
        serialize(f, xi);
      }
    }

    template <typename T>
    void deserialize(std::ifstream& f, SampledTraj<T>& x)
    {
      Index size;
      f.read(reinterpret_cast<char*>(&size), sizeof(size));

      for(Index i = 0 ; i < size ; i++)
      {
        double ti; deserialize(f,ti);
        T xi; deserialize(f,xi);
        x.set(ti,xi);
      }
    }
}