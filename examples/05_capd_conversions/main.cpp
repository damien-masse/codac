#include "capd/capdlib.h"
#include <codac>
using namespace std;
using namespace codac2;
 
int main()
{
  cout<<"Interval"<<endl;
  cout<<"---------"<<endl;
  Interval codac_interval(0., 2.);
  cout << "codac=" << codac_interval << endl;
  capd::Interval capd_interval=to_capd(codac_interval);
  cout << "capd=" << capd_interval << endl;
  cout << "back to codac=" << to_codac(capd_interval) << endl <<endl;

  cout<<"IntervalVector"<<endl;
  cout<<"--------------"<<endl;
  IntervalVector codac_interval_vector({Interval(0., 2.), Interval(1., 3.)});
  cout << "codac=" << codac_interval_vector << endl;
  capd::IVector capd_interval_vector=to_capd(codac_interval_vector);
  cout << "capd=" << capd_interval_vector << endl;
  cout << "back to codac=" << to_codac(capd_interval_vector) << endl <<endl;

  cout<<"IntervalMatrix"<<endl;
  cout<<"-------------"<<endl;
  IntervalMatrix codac_interval_matrix({{Interval(0., 2.), Interval(1., 3.)}, {Interval(2., 4.), Interval(3., 5.)}});
  cout << "codac=" << codac_interval_matrix << endl;
  capd::IMatrix capd_interval_matrix=to_capd(codac_interval_matrix);
  cout << "capd=" << capd_interval_matrix << endl;
  cout << "back to codac=" << to_codac(capd_interval_matrix) << endl <<endl;

  
}