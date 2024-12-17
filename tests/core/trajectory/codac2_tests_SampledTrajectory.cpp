/** 
 *  Codac tests
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <catch2/catch_test_macros.hpp>
#include <codac2_SampledTrajectory.h>
#include <codac2_Approx.h>
#include <codac2_Figure2D.h>

using namespace std;
using namespace codac2;

TEST_CASE("SampledTrajectory")
{
  SampledTrajectory<Vector> x({
    { 0.25, {-0.5,0.5} },
    { 1., {0,0} },
    { 2., {1,0} },
    { 3., {1,1} },
    { 4., {-1,1} },
    { 5., {-1,-1} },
    { 6., {2,-1} }
  });

  CHECK(x.tdomain() == Interval(0.25,6));
  CHECK(x.size() == 2);
  CHECK(x.nb_samples() == 7);
  CHECK(!x.is_empty());
  CHECK(x.codomain() == IntervalVector({{-1,2},{-1,1}}));
  CHECK(x(0.25) == Vector({-0.5,0.5}));
  CHECK(x(1.) == Vector({0,0}));
  CHECK(x(6.) == Vector({2,-1}));
  // Interpolations:
  CHECK(Approx(x(0.75)) == Vector({-1/6.,1/6.}));
  CHECK(x(1.5) == Vector({0.5,0}));
  CHECK(x(5.5) == Vector({0.5,-1}));
  CHECK(x(Interval(1,4)) == IntervalVector({{-1,1},{0,1}}));

  auto x_sampled = x.sampled(0.1);
  CHECK(x_sampled.tdomain() == Interval(0.25,6));
  CHECK(x_sampled.size() == 2);
  CHECK(x_sampled.nb_samples() > 5*x.nb_samples()); // approx
  CHECK(!x_sampled.is_empty());
  CHECK(x_sampled.codomain() == IntervalVector({{-1,2},{-1,1}}));
  CHECK(x_sampled(0.25) == Vector({-0.5,0.5}));
  CHECK(x_sampled(1.) == Vector({0,0}));
  CHECK(x_sampled(6.) == Vector({2,-1}));
  // Interpolations:
  CHECK(Approx(x_sampled(0.75)) == Vector({-1/6.,1/6.}));
  CHECK(Approx(x_sampled(1.5)) == Vector({0.5,0}));
  CHECK(Approx(x_sampled(5.5)) == Vector({0.5,-1}));
  CHECK(x_sampled(Interval(1,4)) == IntervalVector({{-1,1},{0,1}}));

  //DefaultView::set_window_properties({75,75},{700,700});
  //DefaultView::draw_trajectory(x, Color::blue());
  //DefaultView::draw_trajectory(x_sampled, Color::red());
}