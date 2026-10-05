#ifndef _QUCOOPMATCH_H_
#define _QUCOOPMATCH_H_

/*#include <iostream> qucoop.h already included this stuff
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <set>
#include <algorithm>
#include <random>
#include <limits>
#include "Eigen/Dense"*/

using namespace std;
using namespace Eigen;

// ============================================
// GLOBAL TYPES (ROW-MAJOR EVERYWHERE)
// ============================================

using Mat = Matrix<double, Dynamic, Dynamic, RowMajor>;
using Vec = Matrix<double, Dynamic, 1>;

// ============================================
// Utilities
// ============================================

// Generate T matrices (row-major flattening)
Mat tt(int n) { //function name t() already in use by qucoop.h so compiler asks me to change it
    vector<Vec> Tlist;
    Tlist.reserve(n * (n - 1) / 2);

    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j) {
            Mat Ti = Mat::Identity(n, n);
            Ti(i,i)=0; Ti(j,j)=0;
            Ti(i,j)=1; Ti(j,i)=1;

            Vec flat(n*n);
            int k=0;
            for(int r=0;r<n;++r)
                for(int c=0;c<n;++c)
                    flat(k++) = Ti(r,c);

            Tlist.push_back(flat);
        }

    int n_trans = Tlist.size();
    Mat T(n*n, n_trans);
    for(int i=0;i<n_trans;++i)
        T.col(i) = Tlist[i];

    return T;
}

// ============================================
// Matrix Helpers
// ============================================

Vec flatten(const Mat& M) {
    Vec v(M.rows()*M.cols());
    int k=0;
    for(int r=0;r<M.rows();++r)
        for(int c=0;c<M.cols();++c)
            v(k++) = M(r,c);
    return v;
}

Mat reshape_vec_to_mat(const Vec& v, int rows, int cols) {
    if(v.size()!=rows*cols)
        throw runtime_error("reshape size mismatch");

    Mat M(rows,cols);
    int k=0;
    for(int r=0;r<rows;++r)
        for(int c=0;c<cols;++c)
            M(r,c)=v(k++);
    return M;
}

Mat kron(const Mat& A, const Mat& B) {
    int ar=A.rows(), ac=A.cols();
    int br=B.rows(), bc=B.cols();

    Mat C(ar*br, ac*bc);

    for(int i=0;i<ar;++i)
        for(int j=0;j<ac;++j) {
            double a = A(i,j);
            for(int r=0;r<br;++r)
                for(int c=0;c<bc;++c)
                    C(i*br+r, j*bc+c) = a * B(r,c);
        }

    return C;
}

double objective(const Vec &P, const Mat &Q) {
    return (P.transpose() * Q * P)(0,0);
}

// ============================================
// Simulated Annealing
// ============================================

Vec simulated_annealing(const Mat &Dij, const Vec &Dii,
                        int num_reads=50, int num_sweeps=100,
                        double beta_min=0.1, double beta_max=5.0)
{
    int n = Dii.size();
    mt19937 rng(42);
    uniform_real_distribution<double> U(0.0,1.0);

    auto energy = [&](const Vec &q){
        return (q.transpose()*Dij*q)(0,0) + Dii.dot(q);
    };

    double bestE = numeric_limits<double>::infinity();
    Vec bestQ = Vec::Zero(n);

    for(int read=0; read<num_reads; ++read){
        Vec q(n);
        for(int i=0;i<n;++i) q(i)=rng()%2;

        double E = energy(q);

        for(int sweep=0;sweep<num_sweeps;++sweep){
            double frac = double(sweep)/num_sweeps;
            double beta = beta_min * pow(beta_max/beta_min, frac);

            for(int i=0;i<n;++i){
                int idx=rng()%n;
                q(idx)=1-q(idx);

                double newE = energy(q);
                double dE = newE - E;

                if(dE<0 || U(rng)<exp(-beta*dE))
                    E=newE;
                else
                    q(idx)=1-q(idx);
            }
        }

        double finalE = energy(q);
        if(finalE<bestE){
            bestE=finalE;
            bestQ=q;
        }
    }

    return bestQ;
}

// ============================================
// Solver
// ============================================

class QuCOOPMatch {
public:
    Mat Q;
    int n;
    Mat T;
    int n_trans;
    Mat P_sol;
	bool print;
    QuCOOPMatch(const Mat &Q_in, bool p) {
        Q = Q_in;
        n = (int) sqrt(Q.rows());
        T = tt(n);
        n_trans = T.cols();
        P_sol = Mat::Identity(n,n);
		print = p;
    }

    void solve(int maxiter=15) {

        // Symmetrize (IMPORTANT)
        Mat Qsym = 0.5*(Q + Q.transpose());

        SelfAdjointEigenSolver<Mat> es(Qsym);
        double lambda_min = es.eigenvalues()(0);
//		double min_eig = Qsym.eigenvalues().real().minCoeff(); //qucoop
//cout << min_eig << " ==? " << lambda_min << "eigval0\n\n";
		lambda_min = (lambda_min < 0 ? lambda_min : -lambda_min);


        Mat Qshift = Q - 10.0 * lambda_min * Mat::Identity(n*n,n*n);
//		Mat Qshift = Q - 100000.0 * lambda_min * Mat::Identity(n*n,n*n);

        Mat I = Mat::Identity(n,n);
        Mat P = I;

        double sol = objective(flatten(P), Q);

        for(int iter=0; iter<maxiter; ++iter){




            Vec I_flat = flatten(I);

            Mat Pj = kron(I, P.transpose());
            Mat Qj = Pj.transpose()*Qshift*Pj;

            Mat Aj = T;
            for(int i=0;i<n_trans;++i)
                Aj.col(i) -= I_flat;

            Mat Dij = Aj.transpose()*Qj*Aj;

            RowVectorXd tmp = I_flat.transpose()*Qj*Aj;
            Vec Dii = (2.0*tmp).transpose();

            Vec q = simulated_annealing(Dij, Dii);

            Vec delta_flat = I_flat + Aj*q;
            Mat delta = reshape_vec_to_mat(delta_flat, n, n);

            P = delta * P;

            double new_sol = objective(flatten(P), Q);
			if (print)
				cout << "Iter " << iter << " | Objective = " << new_sol << " | Norm(P)^2 = " << P.squaredNorm() << endl;

            if(abs(new_sol-sol)<1e-12)
                break;
            sol=new_sol;
        }
        P_sol = P;
		if (print)
			cout << "\nFinal P:\n" << P_sol << endl;
    }
};

#endif
