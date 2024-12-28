/** 
 *  Codac tests
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <catch2/catch_test_macros.hpp>
#include <codac2_CtcInverse.h>
#include <codac2_analytic_operations.h>
#include <codac2_Figure2D.h>
#include <codac2_Approx.h>
#include <codac2_pave.h>
#include <codac2_Subpaving.h>

using namespace std;
using namespace codac2;

TEST_CASE("CtcInverse")
{
  {
    ScalarVar x,y;
    AnalyticFunction f { {x,y}, x-y };
    CtcInverse<Interval> c(f, Interval(0.));
    CHECK(c.function().input_size() == 2);

    Interval a,b;

    c.contract(a,b);
    CHECK(a == Interval());
    CHECK(b == Interval());

    a = Interval(2,oo);
    b = Interval(2,oo);
    c.contract(a,b);
    CHECK(a == Interval(2,oo));
    CHECK(b == Interval(2,oo));

    a = Interval(1,2);
    b = Interval(1,2);
    c.contract(a,b);
    CHECK(a == Interval(1,2));
    CHECK(b == Interval(1,2));

    a = Interval(10);
    b = Interval(10);
    c.contract(a,b);
    CHECK(a == Interval(10));
    CHECK(b == Interval(10));

    a = Interval(1,5);
    b = Interval(8,9);
    c.contract(a,b);
    CHECK(a == Interval::empty());
    CHECK(b == Interval::empty());

    a = Interval(1,5);
    b = Interval(-oo,oo);
    c.contract(a,b);
    CHECK(a == Interval(1,5));
    CHECK(b == Interval(1,5));
  }

  {
    VectorVar x(2);
    AnalyticFunction f { {x}, x[0]-x[1] };
    CtcInverse_<Interval> c(f, {0});

    //pave(IntervalVector({{-10,10},{-10,10}}), c, 0.1);

    IntervalVector b(2);

    c.contract(b);
    CHECK(b == IntervalVector(2));

    b = IntervalVector({{1,2},{1,2}});
    c.contract(b);
    CHECK(b == IntervalVector({{1,2},{1,2}}));

    b = IntervalVector({{1,oo},{1,oo}});
    c.contract(b);
    CHECK(b == IntervalVector({{1,oo},{1,oo}}));

    b = IntervalVector({10,10});
    c.contract(b);
    CHECK(b == IntervalVector({10,10}));

    b = IntervalVector({{1,5},{8,9}});
    c.contract(b);
    CHECK(b == IntervalVector::empty(2));
  }

  {
    VectorVar x(2);
    AnalyticFunction f { {x}, vec(x[0]-x[1]) };
    CtcInverse_<IntervalVector> c(f, {0});

    //pave(IntervalVector({{-10,10},{-10,10}}), c, 0.1);

    IntervalVector b(2);

    c.contract(b);
    CHECK(b == IntervalVector(2));

    b = IntervalVector({{1,2},{1,2}});
    c.contract(b);
    CHECK(b == IntervalVector({{1,2},{1,2}}));

    b = IntervalVector({{1,oo},{1,oo}});
    c.contract(b);
    CHECK(b == IntervalVector({{1,oo},{1,oo}}));

    b = IntervalVector({10,10});
    c.contract(b);
    CHECK(b == IntervalVector({10,10}));

    b = IntervalVector({{1,5},{8,9}});
    c.contract(b);
    CHECK(b == IntervalVector::empty(2));
  }

  {
    VectorVar x(2);
    AnalyticFunction f { {x}, vec(x[0],sqr(x[0]/7.)+sqr(x[1]/3.)) };
    CtcInverse_<IntervalVector> c(f, {{0,oo},{1}});

    //pave(IntervalVector({{-10,10},{-10,10}}), c, 0.1);

    {
      IntervalVector b({{0,0.8},{-2.28,-1.56}});
      //DefaultView::draw_box(b);
      c.contract(b);
      CHECK(b == IntervalVector::empty(2));
      //DefaultView::draw_box(b,Color::blue());
    }

    {
      IntervalVector b({{4,5.4},{-0.05,2.45}});
      //DefaultView::draw_box(b);
      c.contract(b);
      CHECK(Approx(b,1e-2) == IntervalVector({{4.0397,5.40},{1.9089,2.45}}));
      //DefaultView::draw_box(b,Color::blue());
    }

    {
      IntervalVector b({{6.25,6.7},{0.9,1.85}});
      //DefaultView::draw_box(b);
      c.contract(b);
      CHECK(Approx(b,1e-2) == IntervalVector({{6.25,6.67},{0.9,1.351}}));
      //DefaultView::draw_box(b,Color::blue());
    }

    {
      IntervalVector b({{-6,-5},{0,2}});
      //DefaultView::draw_box(b);
      c.contract(b);
      CHECK(b == IntervalVector::empty(2));
      //DefaultView::draw_box(b,Color::blue());
    }

    {
      IntervalVector b({{2,3},{-1,1}});
      //DefaultView::draw_box(b);
      c.contract(b);
      CHECK(b == IntervalVector::empty(2));
      //DefaultView::draw_box(b,Color::blue());
    }
  }
}

TEST_CASE("ParabolasExample")
{
  ScalarVar u;
  ScalarVar v;
  VectorVar a(4);
  VectorVar a0(3); VectorVar a1(3); VectorVar a2(3);
  VectorVar b0(3); VectorVar b1(3); VectorVar b2(3);

  AnalyticFunction b({u,a0,a1,a2}, sqr(1-u)*a0+2*u*(1-u)*a1+sqr(u)*a2);
  AnalyticFunction f({u,v,a0,a1,a2,b0,b1,b2}, (1-v)*b(u,a0,a1,a2)+v*b(u,b0,b1,b2));

  AnalyticFunction h({a},
     f(a[0],a[1],
       Vector({0,0,0}),Vector({1,0,1}),Vector({2,0,0}),
       Vector({0,2,0}),Vector({1,2,1}),Vector({2,2,0}))
    -f(a[2],a[3],
      Vector({0,0,0.55}),Vector({0,1,-0.45}),Vector({0,2,0.55}),
      Vector({2,0,0.55}),Vector({2,1,-0.45}),Vector({2,2,0.55}))
  );

  CtcInverse_ ctc(h, IntervalVector::zero(3));
  IntervalVector x0 {{0,1},{0,1},{0,0.2},{0,0.2}};
  //draw_while_paving(x0, ctc, 0.001);
  auto p = pave(x0, ctc, 0.01);
  auto cs = p.connected_subsets();
  CHECK(cs.size() == 1);
  CHECK(Approx(cs.begin()->box(),1e-4) == IntervalVector({
    {0.13244,0.201099},
    {0.131459,0.202575},
    {0.132274,0.200001},
    {0.132283,0.200001}
  }));
}