/** 
 *  Codac tests
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <catch2/catch_test_macros.hpp>
#include <codac2_transformations.h>
#include <codac2_math.h>
#include <codac2_AnalyticTraj.h>
#include <codac2_Figure2D.h>
#include <codac2_Approx.h>

using namespace std;
using namespace codac2;

TEST_CASE("Affine transformation")
{
  ScalarVar t;
  AnalyticFunction f1({t}, { cos(2*t)+0.1*cos(10*t), 2*sin(t)+0.1*sin(10*t) });

  double a = PI + PI/3, b = 2.5;
  Vector T({ -3.2, 2.5 });
  //Matrix Rot {{cos(a),-sin(a)},{sin(a),cos(a)}};
  //Rot *= b;
  AnalyticFunction<VectorType> f2({t}, {
    b*cos(a)*f1(t)[0]-b*sin(a)*f1(t)[1] + T[0] + 0.05*cos(100*t),
    b*sin(a)*f1(t)[0]+b*cos(a)*f1(t)[1] + T[1] + 0.05*sin(100*t)
  });

  SampledTraj<Vector> src = AnalyticTraj(f1,{-1,3}).sampled(0.01);
  SampledTraj<Vector> dst = AnalyticTraj(f2,{-1,3}).sampled(0.01);

  //DefaultView::set_window_properties({75,75},{700,700});
  //DefaultView::draw_trajectory(src, Color::blue());
  //DefaultView::draw_trajectory(dst, Color::red());

  Eigen::Affine2d tr = affine_transformation(src, dst);
  SampledTraj<Vector> src2;
  for(const auto& [ti,src_i] : src)
    src2.set(ti, tr*Eigen::Vector2d(src_i));

  //DefaultView::draw_trajectory(src2, Color::purple());

  double scale = tr.linear().col(0).norm();
  Eigen::Matrix2d rotation_matrix = tr.linear() / scale;
  double rotation_angle = atan2(rotation_matrix(1, 0), rotation_matrix(0, 0));
  Eigen::Vector2d translation = tr.translation();

  CHECK(Approx(scale,1e-3) == b);
  CHECK(Approx(rotation_angle,1e-3) == (a-2*PI));
  CHECK(Approx(Vector(translation),1e-3) == T);

  //cout << "  Scale = " << scale << " (" << b << ")" << endl;
  //cout << "  Angle = " << rotation_angle << " (" << a-2*PI << ")" << endl;
  //cout << "  Translation = " << translation << " (" << T << ")" << endl;
  //cout << "End." << endl;
}