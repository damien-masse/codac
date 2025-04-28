/** 
 *  \file codac2_FunctionArgsList.h
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <vector>
#include "codac2_VarBase.h"

namespace codac2
{
  /**
   * \class FunctionArgsList
   * \brief A container class to manage a collection of function arguments.
   */
  class FunctionArgsList : public std::vector<std::shared_ptr<VarBase>>
  {
    public:

      /**
       * \brief Default constructor. It creates an empty list of arguments.
       */
      FunctionArgsList()
      { }

      /**
       * \brief Copy constructor.
       * 
       * Each argument is deep-copied using its ``arg_copy()`` method.
       * 
       * \param args The ``FunctionArgsList`` to copy from.
       */
      FunctionArgsList(const FunctionArgsList& args)
        : std::vector<std::shared_ptr<VarBase>>(args.size())
      {
        size_t i = 0;
        for(const auto& arg : args)
          (*this)[i++] = arg->arg_copy();
      }

      /**
       * \brief Constructor to initialize a ``FunctionArgsList`` from an initializer list.
       * 
       * This constructor creates a ``FunctionArgsList`` from an initializer list of references 
       * to ``VarBase`` objects. Each argument is deep-copied by calling its ``arg_copy()`` on the referenced 
       * objects.
       * 
       * \param args An initializer list of references to ``VarBase`` objects.
       */
      FunctionArgsList(std::initializer_list<std::reference_wrapper<VarBase>> args)
      {
        for(const auto& arg : args)
          push_back(arg.get().arg_copy());
      }

      /**
       * \brief Calculates the total size of the function arguments,
       * as the sum of the sizes of each argument.
       * 
       * \return The sum of the sizes of all arguments.
       */
      Index total_size() const
      {
        Index n = 0;
        for(const auto& ai : *this)
          n += ai->size();
        return n;
      }
  };
}