#ifndef _QMATCH_H_
#define _QMATCH_H_

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

// =====================================================
// TYPES
// =====================================================

using Mat = Matrix<double, Dynamic, Dynamic, RowMajor>;
using Vec = Matrix<double, Dynamic, 1>;
using IntVec = vector<int>;

// =====================================================
// HELPERS
// =====================================================

/*Vec flatten(const Mat& M) { already exists in qucoopmatch.h
    Vec v(M.rows() * M.cols());

    int k = 0;
    for(int i=0;i<M.rows();++i)
        for(int j=0;j<M.cols();++j)
            v(k++) = M(i,j);

    return v;
}*/

double objective2(const Mat& P, const Mat& Q) {
    Vec p = flatten(P);
    return (p.transpose() * Q * p)(0,0);
}

Mat permutation_matrix(const IntVec& perm) {
    int n = perm.size();

    Mat P = Mat::Zero(n,n);

    for(int i=0;i<n;++i)
        P(i, perm[i]) = 1.0;

    return P;
}

Vec table_to_vector(const IntVec& table) {
    int n = table.size();

    Vec v = Vec::Zero(n*n);

    for(int i=0;i<n;++i)
        v(n*i + table[i]) = 1.0;

    return v;
}

// =====================================================
// SIMPLE BINARY SA
// =====================================================

Vec simulated_annealing2(
    const Mat& Q,
    const Vec& q,
    int num_reads = 500,
    int sweeps = 100,
    double beta_min = 0.1,
    double beta_max = 5.0)
{
    int n = q.size();

    mt19937 rng(42);
    uniform_real_distribution<double> U(0.0,1.0);

    auto energy = [&](const Vec& x){
        return (x.transpose()*Q*x)(0,0) + q.dot(x);
    };

    double bestE = numeric_limits<double>::infinity();
    Vec bestX = Vec::Zero(n);

    for(int r=0;r<num_reads;++r){

        Vec x(n);

        for(int i=0;i<n;++i)
            x(i)=rng()%2;

        double E = energy(x);

        for(int s=0;s<sweeps;++s){

            double frac = double(s)/sweeps;
            double beta = beta_min * pow(beta_max/beta_min, frac);

            for(int k=0;k<n;++k){

                int idx = rng()%n;

                x(idx)=1-x(idx);

                double newE = energy(x);
                double dE = newE - E;

                if(dE < 0 || U(rng) < exp(-beta*dE))
                    E = newE;
                else
                    x(idx)=1-x(idx);
            }
        }

        if(E < bestE){
            bestE = E;
            bestX = x;
        }
    }

    return bestX;
}

// =====================================================
// QMATCH
// =====================================================

class QMatch {

public:

    Mat Q;
    int n;
	bool print;
    Mat P_sol;

    QMatch(const Mat& Qin, bool p) {

        Q = Qin;

        n = (int) sqrt(Q.rows());

        P_sol = Mat::Identity(n,n);

		print = p;
    }

    // -------------------------------------------------

    pair<Mat,Vec> get_couplings(
        const IntVec& current_perm,
        const vector<pair<int,int>>& cycles)
    {
        int dim = cycles.size();

        Mat W = Mat::Zero(dim, dim);
        Vec c = Vec::Zero(dim);

        Vec constant = table_to_vector(current_perm);

        for(int i=0;i<dim;++i){

            auto cycle1 = cycles[i];

            Vec v1 = Vec::Zero(n*n);

            {
                int a = cycle1.first;
                int b = cycle1.second;

                int idxA = find(
                    current_perm.begin(),
                    current_perm.end(),
                    a) - current_perm.begin();

                int idxB = find(
                    current_perm.begin(),
                    current_perm.end(),
                    b) - current_perm.begin();

                v1(n*idxA + a) = -1;
                v1(n*idxA + b) =  1;

                v1(n*idxB + b) = -1;
                v1(n*idxB + a) =  1;
            }

            for(int j=0;j<dim;++j){

                auto cycle2 = cycles[j];

                Vec v2 = Vec::Zero(n*n);

                {
                    int a = cycle2.first;
                    int b = cycle2.second;

                    int idxA = find(
                        current_perm.begin(),
                        current_perm.end(),
                        a) - current_perm.begin();

                    int idxB = find(
                        current_perm.begin(),
                        current_perm.end(),
                        b) - current_perm.begin();

                    v2(n*idxA + a) = -1;
                    v2(n*idxA + b) =  1;

                    v2(n*idxB + b) = -1;
                    v2(n*idxB + a) =  1;
                }

                W(i,j) = (v1.transpose() * Q * v2)(0,0);

                if(i==0)
                    c(j) += (constant.transpose()*Q*v2)(0,0);
            }

            c(i) += (v1.transpose()*Q*constant)(0,0);
        }

        return {W,c};
    }

    // -------------------------------------------------

    IntVec update_table(
        const IntVec& table,
        const vector<pair<int,int>>& cycles,
        const Vec& decisions)
    {
        IntVec result = table;

        for(size_t i=0;i<cycles.size();++i){

            if(decisions(i) > 0.5){

                int a = cycles[i].first;
                int b = cycles[i].second;

                int idxA = find(
                    table.begin(),
                    table.end(),
                    a) - table.begin();

                int idxB = find(
                    table.begin(),
                    table.end(),
                    b) - table.begin();

                swap(result[idxA], result[idxB]);
            }
        }

        return result;
    }

    // -------------------------------------------------

    vector<pair<int,int>>
    next_pairs(const vector<pair<int,int>>& pairs)
    {
        vector<pair<int,int>> out;

        if(pairs.size() > 2){

            out.push_back({
                pairs[0].second,
                pairs[2].second
            });

            for(size_t k=2;k<pairs.size()-1;++k){

                out.push_back({
                    pairs[k-1].first,
                    pairs[k+1].second
                });
            }

            out.push_back({
                pairs[pairs.size()-2].first,
                pairs.back().first
            });
        }
        else {

            out.push_back({
                pairs[0].second,
                pairs[0].first
            });
        }

        return out;
    }

    // -------------------------------------------------

    pair<IntVec, double> optimize(
        IntVec initial_perm = {})
    {
        IntVec perm;
        double E;

        if(initial_perm.empty()){

            perm.resize(n);

            for(int i=0;i<n;++i)
                perm[i]=i;
        }
        else
            perm = initial_perm;

        // ---------------------------------------------
        // Build cycle lists
        // ---------------------------------------------

        vector<vector<pair<int,int>>> cycle_lists;

        vector<pair<int,int>> init_pairs;

        for(int i=0;i<n/2;++i)
            init_pairs.push_back({i, n-i-1});

        shuffle(
            init_pairs.begin(),
            init_pairs.end(),
            mt19937(random_device{}()));

        cycle_lists.push_back(init_pairs);

        auto current_pairs = next_pairs(init_pairs);

        for(int i=0;i<n-2;++i){

            cycle_lists.push_back(current_pairs);

            current_pairs = next_pairs(current_pairs);
        }

        // ---------------------------------------------
        // Optimization
        // ---------------------------------------------

        for(auto& cycles : cycle_lists){

            pair<Mat,Vec> result = get_couplings(perm, cycles);

            Mat W = result.first;
            Vec b = result.second;

            Vec decisions = simulated_annealing2(W,b);

            perm = update_table(
                perm,
                cycles,
                decisions);

            E =
                (table_to_vector(perm).transpose()
                 * Q
                 * table_to_vector(perm))(0,0);
            if(print)                
				cout << "Current energy: " << E << endl;
        }

        return make_pair(perm, E);
    }

    // -------------------------------------------------

    vector<double> solve(
		int iterations = 20)
{
    vector<double> history;

    double init_obj =
        objective2(Mat::Identity(n,n), Q);

    history.push_back(init_obj);

    IntVec current_perm;

    double E = init_obj;

    int stagnant_iterations = 0;

    const double tol = 1e-12;

    for(int i=0;i<iterations;++i){

        pair<IntVec,double> result =
            optimize(current_perm);

        current_perm = result.first;

        double newE = result.second;

        history.push_back(newE);

        if (print)
            cout << "Iteration "
                 << i
                 << " done - Energy: "
                 << newE << "\n";

        // -----------------------------------------
        // Check convergence
        // -----------------------------------------

        if(abs(newE - E) < tol)
            stagnant_iterations++;
        else
            stagnant_iterations = 0;

        if(stagnant_iterations >= 3){

            if(print)
                cout << "\nStopping early: "
                     << "energy unchanged for "
                     << 3
                     << " consecutive iterations.\n";

            break;
        }

        E = newE;
    }

    P_sol = permutation_matrix(current_perm);

    double final_obj =
        objective2(P_sol, Q);

    history.push_back(final_obj);

    return history;
}
};

#endif
