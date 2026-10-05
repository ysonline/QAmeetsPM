#ifndef _QGM_H_
#define _QGM_H_

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

//using namespace std;
//using namespace Eigen;

double objective3(const VectorXd &P_flat, const MatrixXd &Q) {
    return P_flat.transpose() * Q * P_flat;
}

// -----------------------------
// Simulated annealing
// -----------------------------
VectorXd simulated_annealing(const MatrixXd &Dij, const VectorXd &Dii,
                             int num_reads=50, int num_sweeps=100,
                             double beta_min=0.1, double beta_max=5.0)
{
    int n = Dii.size();
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> realDist(0.0,1.0);

    auto energy = [&](const VectorXd &q){ return q.transpose()*Dij*q + Dii.dot(q); };

    double bestE = std::numeric_limits<double>::infinity();
    VectorXd bestQ = VectorXd::Zero(n);

    for(int read=0; read<num_reads; ++read){
        VectorXd q(n);
        for(int i=0;i<n;++i) q(i) = rng()%2;

        double E = energy(q);

        for(int sweep=0; sweep<num_sweeps; ++sweep){
            double frac = double(sweep)/num_sweeps;
            double beta = beta_min * pow(beta_max/beta_min, frac);

            for(int i=0;i<n;++i){
                int idx = rng()%n;
                q(idx) = 1.0 - q(idx);
                double newE = energy(q);
                double dE = newE - E;
                if(dE<0 || realDist(rng) < exp(-beta*dE)) E = newE;
                else q(idx) = 1.0 - q(idx); // revert
            }
        }
        double finalE = energy(q);
        if(finalE<bestE){ bestE = finalE; bestQ = q; }
    }
    return bestQ;
}

// -----------------------------
// QGM class
// -----------------------------
class QGM {
public:
    MatrixXd Q;           // full Q (n^2 x n^2)
    int n;                // problem dimension
    MatrixXd P_sol;       // n x n perm matrix solution
//    string method;        // "baseline" or "row_wise"
	bool print;
    QGM(const MatrixXd &Q_in, bool p){
        Q = Q_in;
        n = int(sqrt(Q.rows()));
        P_sol = MatrixXd::Identity(n,n);
//        method = method_in;
		print = p;
    }

    // Build A matrix (2n x n^2) using column-major flattening
    MatrixXd build_A() {
        int n2 = n*n;

        // Row constraints (sum across each row = 1)
        MatrixXd A1 = MatrixXd::Zero(n, n2);
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                int colIdx = j*n + i; // column-major flattening
                A1(i, colIdx) = 1.0;
            }
        }

        // Column constraints (sum across each column = 1)
        MatrixXd A2 = MatrixXd::Zero(n, n2);
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                int colIdx = i*n + j; // column-major flattening
                A2(i, colIdx) = 1.0;
            }
        }

        // Stack vertically
        MatrixXd A(2*n, n2);
        A.block(0,0,n,n2) = A1;
        A.block(n,0,n,n2) = A2;
        return A;
    }

    // Build BQM matrices
    pair<MatrixXd, VectorXd> get_bqm_matrices() {
        int n2 = n*n;
        VectorXd C = Q.diagonal();
        MatrixXd W = Q;
        for(int i=0;i<n2;i++) W(i,i) = 0.0;

        MatrixXd A = build_A();
        VectorXd b = VectorXd::Ones(2*n);

        MatrixXd Q_bqm;
        VectorXd qvec;

		int method = 1;//3;//2;//1;//2;3;
        if(method == 1){//"baseline"){
            
			
			
			
			//double Lambda_baseline = W.cwiseAbs().sum() / 2.0;
			double Lambda_baseline = W.cwiseAbs().sum() / 20.0; //less penalty


            Q_bqm = W + Lambda_baseline * (A.transpose() * A);
            VectorXd btA = (b.transpose() * A).transpose();
            qvec = C - 2.0 * Lambda_baseline * btA;
        } else if(method == 2){//"row_wise"){
            VectorXd optimizing = W.rowwise().sum() + W.colwise().sum().transpose() - W.diagonal();
            optimizing = optimizing.cwiseAbs();
            double MaxGrad = optimizing.maxCoeff();

            VectorXd Lambda_row_wise(2*n);
            for(int r=0;r<2*n;++r){
                double maxval = -std::numeric_limits<double>::infinity();
                for(int i=0;i<n2;++i){
                    double val = A(r,i) * optimizing(i);
                    if(val > maxval) maxval = val;
                }
                Lambda_row_wise(r) = maxval + 0.5 * MaxGrad;
            }

            MatrixXd LambdaA(2*n, n2);
            for(int r=0;r<2*n;++r)
                for(int i=0;i<n2;++i)
                    LambdaA(r,i) = Lambda_row_wise(r) * A(r,i);

            Q_bqm = W + (LambdaA.transpose() * A);

            VectorXd bt_LA = VectorXd::Zero(n2);
            for(int i=0;i<n2;++i){
                double acc = 0.0;
                for(int r=0;r<2*n;++r) acc += b(r) * (Lambda_row_wise(r) * A(r,i));
                bt_LA(i) = acc;
            }
            qvec = C - 2.0 * bt_LA;
        } else {
            double Lambda_baseline = W.cwiseAbs().sum() / 2.0;
            Q_bqm = W + Lambda_baseline * (A.transpose() * A);
            VectorXd btA = (b.transpose() * A).transpose();
            qvec = C - 2.0 * Lambda_baseline * btA;
        }

        return {Q_bqm, qvec};
    }

    vector<double> solve(bool simulated_anneal=true) {
        int n2 = n*n;
        // Identity start
        MatrixXd P0 = MatrixXd::Identity(n,n);
        VectorXd P0_flat = Map<VectorXd>(P0.data(), n2); // column-major
        vector<double> history;
        history.push_back(objective3(P0_flat, Q));

        // Get BQM
        MatrixXd Dij; VectorXd Dii;
        tie(Dij, Dii) = get_bqm_matrices();

        // Solve BQM
        VectorXd q_sol = simulated_annealing(Dij, Dii);

        // Reshape to n x n
        MatrixXd Pmat = Map<MatrixXd>(q_sol.data(), n, n); // column-major
        P_sol = Pmat;

        VectorXd P_flat = Map<VectorXd>(P_sol.data(), n2);
        double sol_val = objective3(P_flat, Q);
        history.push_back(sol_val);

		if (print)
			cout << "QGM | objective = " << sol_val << endl //QGM method = " << method << 
				 << "Approx permutation matrix (0/1 entries):\n" << P_sol << endl;
        return history;
    }
};

#endif
