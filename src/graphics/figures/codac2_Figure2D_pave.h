/** 
 *  \file codac2_drawwhilepaving.h
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include "codac2_Ctc.h"
#include "codac2_Sep.h"
#include "codac2_Figure2D.h"
#include "codac2_template_tools.h"

namespace codac2
{
  void init_fig(std::shared_ptr<Figure2D>& fig, const IntervalVector& x0);

  inline std::function<void(Figure2D&,const IntervalVector&,const StyleProperties&)> cartesian_drawing()
  {
    return [](Figure2D& fig, const IntervalVector& x, const StyleProperties& s) {
      fig.draw_box(x,s);
    };
  }

  inline std::function<void(Figure2D&,const IntervalVector&,const StyleProperties&)> polar_drawing()
  {
    return [](Figure2D& fig, const IntervalVector& x, const StyleProperties& s) {
      fig.draw_pie({0,0},x[0],x[1],s);
    };
  }
  
  /**
   * \brief Draws the paving while it is being computed, with custom drawing function
   * 
   * \param x0 Initial box
   * \param c Contractor to use
   * \param eps Accuracy of the paving algorithm (the undefined boxes will have their max_diam <= eps)
   * \param fig Figure to draw on (optionnal, if not provided, the default figure will be used)
   * \param draw_box Optional custom drawing function (for instance, if one want to draw in polar coordinates)
   */  
  template<typename C>
    requires IsCtcBaseOrPtr<C,IntervalVector>
  inline void draw_while_paving(const IntervalVector& x0, const C& c, double eps,
    const std::function<void(Figure2D&,const IntervalVector&,const StyleProperties&)>& draw_box,
    std::shared_ptr<Figure2D> fig = nullptr)
  {
    assert_release(eps > 0.);
    assert_release(c.size() >= 2 && "cannot reveal 1d contractors");
    
    if(!fig)
      fig = DefaultFigure::selected_fig();
    init_fig(fig, x0);

    clock_t t_start = clock();

    if(x0.size() > 2)
      draw_box(*fig, x0, StyleProperties::outside());

    std::list<IntervalVector> l { x0 };
    Index n = 0;

    while(!l.empty())
    {
      IntervalVector x = l.front(), prev_x = x;
      l.pop_front();

      ctc(c).contract(x);

      if(x0.size() == 2)
        for(const auto& bi : prev_x.diff(x))
          draw_box(*fig, bi, StyleProperties::outside());

      if(!x.is_empty())
      {
        if(x.max_diam() < eps)
        {
          n++;
          draw_box(*fig, x, StyleProperties::boundary());
        }

        else
        {
          auto p = x.bisect_largest();
          l.push_back(p.first); l.push_back(p.second);
        }
      }
    }

    printf("Computation time: %.4fs, %ld boxes\n", (double)(clock()-t_start)/CLOCKS_PER_SEC, n);
  }
  
  /**
   * \brief Draws the paving while it is being computed, in default Cartesian coordinates
   * 
   * \param x0 Initial box
   * \param c Contractor to use
   * \param eps Accuracy of the paving algorithm (the undefined boxes will have their max_diam <= eps)
   * \param fig Figure to draw on (optionnal, if not provided, the default figure will be used)
   */  
  template<typename C>
    requires IsCtcBaseOrPtr<C,IntervalVector>
  inline void draw_while_paving(const IntervalVector& x0, const C& c, double eps, std::shared_ptr<Figure2D> fig = nullptr)
  {
    return draw_while_paving(x0, c, eps, cartesian_drawing(), fig);
  }

  /**
   * \brief Draws the paving while it is being computed, with custom drawing function
   * 
   * \param x0 Initial box
   * \param s Separator to use
   * \param eps Accuracy of the paving algorithm (the undefined boxes will have their max_diam <= eps)
   * \param draw_box Custom drawing function (for instance, if one want to draw in polar coordinates)
   * \param fig Figure to draw on (optionnal, if not provided, the default figure will be used)
   */
  template<typename S>
    requires IsSepBaseOrPtr<S>
  inline void draw_while_paving(const IntervalVector& x0, const S& s, double eps,
    const std::function<void(Figure2D&,const IntervalVector&,const StyleProperties&)>& draw_box,
    std::shared_ptr<Figure2D> fig = nullptr)
  {
    assert_release(eps > 0.);
    assert_release(size_of(s) >= 2 && "cannot reveal 1d separators");

    if(!fig)
      fig = DefaultFigure::selected_fig();
    init_fig(fig, x0);
    
    clock_t t_start = clock();

    std::list<IntervalVector> l { x0 };

    while(!l.empty())
    {
      IntervalVector x = l.front();
      l.pop_front();

      auto x_sep = sep(s).separate(x);
      auto boundary = x_sep.inner & x_sep.outer;

      for(const auto& bi : x.diff(x_sep.inner))
        draw_box(*fig, bi, StyleProperties::inside());

      for(const auto& bi : x.diff(x_sep.outer))
        draw_box(*fig, bi, StyleProperties::outside());

      if(!boundary.is_empty())
      {
        if(boundary.max_diam() <= eps)
          draw_box(*fig, boundary, StyleProperties::boundary());

        else
        {
          auto p = boundary.bisect_largest();
          l.push_back(p.first); l.push_back(p.second);
        }
      }
    }
    
    printf("Computation time: %.4fs\n", (double)(clock()-t_start)/CLOCKS_PER_SEC);
  }

  /**
   * \brief Draws the paving while it is being computed, in default Cartesian coordinates
   * 
   * \param x0 Initial box
   * \param s Separator to use
   * \param eps Accuracy of the paving algorithm (the undefined boxes will have their max_diam <= eps)
   * \param fig Figure to draw on (optionnal, if not provided, the default figure will be used)
   */
  template<typename S>
    requires IsSepBaseOrPtr<S>
  inline void draw_while_paving(const IntervalVector& x0, const S& s, double eps, std::shared_ptr<Figure2D> fig = nullptr)
  {
    return draw_while_paving(x0, s, eps, cartesian_drawing(), fig);
  }
}