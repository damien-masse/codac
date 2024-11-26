/** 
 *  codac2_Ellipsoid.cpp
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Morgan Louédec
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include "codac2_Ellipsoid_utils.h"
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/KroneckerProduct>

using namespace std;
using codac2::Matrix;
using codac2::Vector;
using codac2::BoolInterval;

namespace codac2
{
    Matrix solve_discrete_lyapunov(const Matrix& a,const Matrix& q)
    {
        // implementation of the scipy solver for the discrete lyapunov equation (real matrix only)
        // works well under dimension 10
        // https://github.com/scipy/scipy/blob/v1.14.1/scipy/linalg/_solvers.py#L235-L323
        // Solves the discrete Lyapunov equation :math:`AXA^H - X + Q = 0`
        assert(a.rows() == a.cols());
        assert(a.rows() == q.rows());
        assert(a.cols() == q.cols());

        Eigen::MatrixXd lhs = Eigen::KroneckerProduct(a, a);
        lhs = Eigen::MatrixXd::Identity(lhs.rows(),lhs.cols()) - lhs;
        Eigen::MatrixXd x = lhs.colPivHouseholderQr().solve((Eigen::VectorXd)q.reshaped());
        return Matrix(x.reshaped(q.rows(),q.cols()));
    }

    BoolInterval stability_analysis(const AnalyticFunction<VectorOpValue> &f, int alpha_max, Ellipsoid &e, Ellipsoid &e_out)
    {
        assert_release(f.args().size() == 1 && "f must have only one arg");

        // get the Jacobian of f at the origin
        int n = f.input_size();
        Vector origin(Eigen::VectorXd::Zero(n));
        Matrix J = f.diff(IntervalVector(origin)).mid();

        // solve the axis aligned discrete lyapunov equation J.T * P * J − P = −J.T * J
        Matrix P = solve_discrete_lyapunov(J.transpose(),J.transpose()*J); // TODO solve the Lyapunov equation !!!
        Matrix G0((P.inverse()).sqrt());
        int alpha = 0;

        while(alpha <= alpha_max)
        {
            e = Ellipsoid(origin, std::pow(10,-alpha) * G0);
            e_out = nonlinear_mapping(e,f);
            cout << "\nwith alpha = " << alpha << endl;
            cout << "e is\n" << e << endl;
            cout << "e_out is\n" << e_out << endl;

            if(e_out.is_concentric_subset(e) == BoolInterval::TRUE)
            {
                cout << "The system is stable" << endl;
                cout << "Domain of attraction :\n" << e_out << endl;
                return BoolInterval::TRUE;
            }
            alpha++;
        }
        cout << "The method is not able to conclude on the stability" << endl;
        return BoolInterval::UNKNOWN;
    }
}