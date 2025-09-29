/** 
 *  codac2_drawwhilepaving.cpp
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include "codac2_drawwhilepaving.h"

using namespace std;
using namespace codac2;

namespace codac2
{
  void init_fig(std::shared_ptr<Figure2D>& fig, const IntervalVector& x0)
  {
    if(!fig)
    {
      fig = std::make_shared<Figure2D>(DEFAULT_FIG_NAME, GraphicOutput::VIBES);
      fig->set_window_properties({20,20}, {800,800});
      fig->set_as_default();
    }

    fig->set_axes(axis(0,x0[0],"x_1"), axis(1,x0[1],"x_2"));
    fig->auto_scale();
  }
}