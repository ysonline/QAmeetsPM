#ifndef _QUCOOP_H_
#define _QUCOOP_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include "Eigen/Dense"
#include <vector>
#include <cmath>
#include <tuple>
#include <set>
#include <algorithm>
#include <random>

using namespace std;
using namespace Eigen;

// ============================================
// Utilities
// ============================================

// Generate T matrices for transpositions (row-major flattening)
MatrixXd t(int n) {
    vector<VectorXd> Tlist;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            MatrixXd Ti = MatrixXd::Identity(n, n);
            Ti(i, i) = 0;
            Ti(j, j) = 0;
            Ti(i, j) = 1;
            Ti(j, i) = 1;

            // Flatten row-major to match NumPy default
            VectorXd Ti_flat(n * n);
            for (int r = 0; r < n; ++r)
                for (int c = 0; c < n; ++c)
                    Ti_flat[r * n + c] = Ti(r, c);
            Tlist.push_back(Ti_flat);
        }
    }

    int n_trans = Tlist.size();
    MatrixXd T(n * n, n_trans);
    for (int k = 0; k < n_trans; ++k)
        T.col(k) = Tlist[k];
    return T;
}

// Evaluate permutation from x and T (row-major)
tuple<VectorXd, MatrixXd> eval_p(const MatrixXd &T, const VectorXd &x) {
    int n = (int) sqrt(T.rows()); //used n*n as # rows so sqrt here will be an int for sure (no loss of data)
    int n_trans = T.cols(); //= n * (n - 1) / 2 if no variable is set to fixed landmark matches, e.g. fixed_vars_1.empty(), o/w it is of shorter length
    MatrixXd I = MatrixXd::Identity(n, n);
    MatrixXd P = I;
    MatrixXd dP = MatrixXd::Zero(n * n, n_trans);

    for (int i = 0; i < n_trans; ++i) {
        MatrixXd Ti(n, n);
        for (int r = 0; r < n; ++r)
            for (int c = 0; c < n; ++c)
                Ti(r, c) = T(r * n + c, i);
        P = P * (x(i) * (Ti - I) + I);

        MatrixXd dp = I;
        for (int j = 0; j < i; ++j) {
            MatrixXd Tj(n, n);
            for (int r = 0; r < n; ++r)
                for (int c = 0; c < n; ++c)
                    Tj(r, c) = T(r * n + c, j);
            dp = dp * (x(j) * (Tj - I) + I);
        }

        MatrixXd Ti_i(n, n);
        for (int r = 0; r < n; ++r) //this nested loop can be replaced with the single MatrixXd Ti_i = Ti; statement (disable MatrixXd Ti_i(n, n); declaration as well)
            for (int c = 0; c < n; ++c)
                Ti_i(r, c) = T(r * n + c, i);
        dp = dp * (Ti_i - I);
//if (Ti != Ti_i) exit(0); cout << Ti << "\nvs.\n" << Ti_i << "\n\n"; never exits

        for (int j = i + 1; j < n_trans; ++j) {
            MatrixXd Tj(n, n);
            for (int r = 0; r < n; ++r)
                for (int c = 0; c < n; ++c)
                    Tj(r, c) = T(r * n + c, j);
            dp = dp * (x(j) * (Tj - I) + I);
        }

        VectorXd dp_flat(n * n);
        for (int r = 0; r < n; ++r)
            for (int c = 0; c < n; ++c)
                dp_flat[r * n + c] = dp(r, c);
        dP.col(i) = dp_flat;
    }

    VectorXd P_flat(n * n);
    for (int r = 0; r < n; ++r)
        for (int c = 0; c < n; ++c)
            P_flat[r * n + c] = P(r, c);

    return make_tuple(P_flat, dP);
}

// Objective
double objective(const VectorXd &P, const MatrixXd &Q) {
    return P.transpose() * Q * P;
}

// ============================================
// BQM Fixing Utility
// ============================================

struct BQMReduced {
    MatrixXd Dij;
    VectorXd Dii;
    double constant;
};

BQMReduced fix_variables(const MatrixXd &Dij, const VectorXd &Dii,
                         const vector<int> &fixed_to_1,
                         const vector<int> &fixed_to_0) {
    int n = Dij.rows();
    VectorXd fixed_vals = VectorXd::Zero(n);
    vector<int> fixed_indices;
    for (int i : fixed_to_1) { fixed_vals(i) = 1.0; fixed_indices.push_back(i); }
    for (int i : fixed_to_0) { fixed_vals(i) = 0.0; fixed_indices.push_back(i); }

    set<int> fixed_set(fixed_indices.begin(), fixed_indices.end());

    double constant = 0.0;
    VectorXd new_Dii = Dii;

    for (int i = 0; i < n; ++i) {
        if (fixed_set.count(i)) {
            constant += Dii(i) * fixed_vals(i);
            for (int j = 0; j < n; ++j)
                constant += 0.5 * Dij(i, j) * fixed_vals(i) * fixed_vals(j);
        } else {
            for (int j = 0; j < n; ++j)
                if (fixed_set.count(j))
                    new_Dii(i) += Dij(i, j) * fixed_vals(j);
        }
    }

    vector<int> free_idx;
    for (int i = 0; i < n; ++i)
        if (!fixed_set.count(i))
            free_idx.push_back(i);

    int m = free_idx.size();
    MatrixXd Dij_reduced(m, m);
    VectorXd Dii_reduced(m);

    for (int i = 0; i < m; ++i) {
        Dii_reduced(i) = new_Dii(free_idx[i]);
        for (int j = 0; j < m; ++j)
            Dij_reduced(i, j) = Dij(free_idx[i], free_idx[j]);
    }

    return {Dij_reduced, Dii_reduced, constant};
}

// ============================================
// QuCOOP Class
// ============================================

class QuCOOP {
public:
    MatrixXd Q;
    int n;
    MatrixXd T;
    int n_trans;
    MatrixXd P_sol;
	double sol;
    vector<int> fixed_vars_1;
    vector<int> fixed_vars_0;
	bool print;
    QuCOOP(const MatrixXd &Q_in, bool p, const vector<pair<int,int>> &fixed_assignments = {}) {
        Q = Q_in;
        n = (int) sqrt(Q.rows()); //used n*n as # rows so sqrt here will be an int for sure (no loss of data)
        T = t(n);
        n_trans = T.cols(); //= n * (n - 1) / 2 if no variable is set to fixed landmark matches, e.g. fixed_vars_1.empty(), o/w it is of shorter length
        P_sol = MatrixXd::Identity(n, n);
        set_fixed_variables(fixed_assignments);
		print = p;
    }

    void set_fixed_variables(const vector<pair<int,int>> &fixed_assignments) {
        vector<pair<int,int>> assignments;
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j)
                assignments.push_back({i, j});

        set<int> fixed_indices;
        for (auto pair : fixed_assignments)
            fixed_indices.insert(pair.first), fixed_indices.insert(pair.second);

        fixed_vars_1.clear();
        fixed_vars_0.clear();
        for (size_t idx = 0; idx < assignments.size(); ++idx) {
            auto pair = assignments[idx];
            if (find(fixed_assignments.begin(), fixed_assignments.end(), pair) != fixed_assignments.end())
                fixed_vars_1.push_back(idx);
            else if (fixed_indices.count(pair.first) || fixed_indices.count(pair.second))
                fixed_vars_0.push_back(idx);
        }
    }

	double startTemp = 2.0; //will be updated by the caller based on maxGeoDist (which makes it a litlle bit compatible w/ delta but not sure if this is a good solution)
	VectorXd qPrev;
	VectorXd SA(const VectorXd& Dii, const MatrixXd& Dij, int nMaxIters = 50, int nSweeps = 100, double coolingRate = 0.99)
	{
        //solves the problem bestQ.T * (Dij + diag Dii) * bestQ, bestQ in {0, 1} for the input Dij and Dii matrices via simulated annealing and returns the solution vector bestQ
		//note that bestQ.T * Dij * bestQ + Dii.dot(q) gives the same objective function (and avoid diagonalization so tiny faster)
	
//		std::mt19937 rng(42); //fixed seed (10 calls to realDist(rng) always generates 0.796543 0.183435 0.779691 0.59685 0.445833 0.0999749 0.459249 0.333709 0.142867 0.650888; similarly 7 calls to rng() % 2 always generates 0 1 0 0 0 1 0)
		std::mt19937 rng( (unsigned) time(NULL) ); //non-fixed seed that changes with every execution (recommended)
		std::uniform_real_distribution<double> realDist(0.0, 1.0);

		auto energy = [&](const VectorXd &q){ return q.transpose() * Dij * q + Dii.dot(q); };

		int m = Dii.size(), bitIdx;
		VectorXd bestQ = VectorXd::Zero(m);
		double E, newE, delta, temp = startTemp, bestE = INF, finalE;
		for (int i = 0; i < nMaxIters; i++) //multistart: each of nMaxIters iterations initialized by a random q vector tries to update bestQ below
		{
			VectorXd q(m);
			for (int j = 0; j < m; ++j)
				q(j) = rng() % 2; //random initialization
//				q(j) = qPrev(j); //alternatively initialize with the q_free of the previous solve() iteration (not a good idea 'cos multistarts are supposed to start randomly at each i iteration instead of a fixed qPrev start; bitIdx below still random but fixed qPrev is too much)
//				q(j) = (realDist(rng) < 0.75 ? qPrev(j) : rng() % 2); //inspire from qPrev instead of directly using it (75% of the time)

			E = energy(q); //MatrixXd Dii_Diag = Dii.asDiagonal(); if (fabs(q.transpose() * (Dij + Dii_Diag) * q - E)>TINY) exit(0); never exits

			temp = startTemp; //alternatively never reset it here, i.e., disable this line
			for (int sweep = 0; sweep < nSweeps; ++sweep)  //cooling-down: temp-based
			{
				//flip bits and check the new energy (revert in reject state)
				for (int b = 0; b < m; b++)
				{
					bitIdx = rng() % m, q(bitIdx) = 1 - q(bitIdx); //1 becomes 0 and 0 becomes 1 hence the flip (bitIdx is a random int in [0,m) interval)
					newE = energy(q), delta = newE - E;
					if (delta < 0.0 || realDist(rng) < exp(-delta / temp)) //newE is better/smaller or temparature-based random access is granted (accept state)
						E = newE;
					else
						q(bitIdx) = 1 - q(bitIdx); //revert
				}				
				temp *= coolingRate; //cool down
//				if (i % 10 == 0 && sweep % 25 == 0)
//					cout << "iterInner: " << i << "\tcost: " << E << "\ttemp: " << temp << endl;
			}
			finalE = energy(q); //energy based on the flipped q
			if (finalE < bestE)
				bestE = finalE, bestQ = q;
//			if (i % 10 == 0)
//				cout << "iterationOuter: " << i << "\tcost: " << finalE << "\tbestCost: " << bestE << "\n\n";
		}
//for (int i = 0; i < m; ++i) bestQ(i) = rng() % 2;//same as (int) (rand() % 2); //discard SA and use a random solution (for debugging only)
//for (int i = 0; i < m; ++i) cout << bestQ(i) << " ";cout <<"returnedSA\n\n";
		return bestQ;
	}

	VectorXd simulatedAnnealing(const MatrixXd &Dij, const VectorXd &Dii, int nMaxIters = 50, int nSweeps = 100, double beta_min = 0.1, double beta_max = 5.0)
	{
		//same as SA except this one does not have tempature vs. delta compatibility problems

//		std::mt19937 rng(42); //fixed seed (10 calls to realDist(rng) always generates 0.796543 0.183435 0.779691 0.59685 0.445833 0.0999749 0.459249 0.333709 0.142867 0.650888; similarly 7 calls to rng() % 2 always generates 0 1 0 0 0 1 0)
		std::mt19937 rng( (unsigned) time(NULL) ); //non-fixed seed that changes with every execution (recommended)
		std::uniform_real_distribution<double> realDist(0.0, 1.0);

		auto energy = [&](const VectorXd &q){ return q.transpose() * Dij * q + Dii.dot(q); };

		int m = Dii.size(), bitIdx;
		VectorXd bestQ = VectorXd::Zero(m);
		double E, newE, delta, bestE = INF, finalE;
		for (int i = 0; i < nMaxIters; i++) //multistart: each of nMaxIters iterations initialized by a random q vector tries to update bestQ below
		{
			VectorXd q(m);
			for (int j = 0; j < m; ++j)
				q(j) = rng() % 2; //random initialization
//				q(j) = qPrev(j); //alternatively initialize with the q_free of the previous solve() iteration (not a good idea 'cos multistarts are supposed to start randomly at each i iteration instead of a fixed qPrev start; bitIdx below still random but fixed qPrev is too much)
//				q(j) = (realDist(rng) < 0.75 ? qPrev(j) : rng() % 2); //inspire from qPrev instead of directly using it (75% of the time)

			E = energy(q); //MatrixXd Dii_Diag = Dii.asDiagonal(); if (fabs(q.transpose() * (Dij + Dii_Diag) * q - E)>TINY) exit(0); never exits

			for (int sweep = 0; sweep < nSweeps; ++sweep) //cooling-down: beta-based
			{
				double frac = double(sweep) / nSweeps;
				double beta = beta_min * pow(beta_max / beta_min, frac);

				//flip bits and check the new energy (revert in reject state)
				for (int b = 0; b < m; b++)
				{
					bitIdx = rng() % m, q(bitIdx) = 1 - q(bitIdx); //1 becomes 0 and 0 becomes 1 hence the flip (bitIdx is a random int in [0,m) interval)
					newE = energy(q), delta = newE - E;
					if (delta < 0 || realDist(rng) < exp(-beta*delta)) //newE is better/smaller or temparature-based random access is granted (accept state)
						E = newE;
					else
						q(bitIdx) = 1 - q(bitIdx); //revert
				}
//				if (i % 10 == 0 && sweep % 25 == 0)
//					cout << "iterInner: " << i << "\tcost: " << E << "\tbeta: " << beta << endl;
			}
			finalE = energy(q); //energy based on the flipped q
			if (finalE < bestE)
				bestE = finalE, bestQ = q;
//			if (i % 10 == 0)
//				cout << "iterationOuter: " << i << "\tcost: " << finalE << "\tbestCost: " << bestE << "\n\n";
		}
//for (int i = 0; i < m; ++i) bestQ(i) = rng() % 2;//same as (int) (rand() % 2); //discard SA and use a random solution (for debugging only)
//for (int i = 0; i < m; ++i) cout << bestQ(i) << " ";cout <<"returnedSimAnn\n\n";
		return bestQ;
	}

    // ================================
    // Main solver (with iteration loop)
    // ================================
    void solve(int maxiter = 15) {
        VectorXd x0 = VectorXd::Zero(n_trans); //means identity permutation matrix for the first j=0 iteration (may be a very bad initialization so consider multi-starts)
//		std::mt19937 rng( (unsigned) time(NULL) ); //non-fixed seed that changes with every execution (recommended)	
//		for (int i = 0; i < n_trans; i++)
//			x0(i) = rng() % 2;

        double min_eig = Q.eigenvalues().real().minCoeff(), sol_prev = 1e20;
		MatrixXd Q_mod = Q - 10 * min_eig * MatrixXd::Identity(n * n, n * n);
		//MatrixXd Q_mod = Q - 20 * min_eig * MatrixXd::Identity(n * n, n * n); //for bigger problems (big n), if results strange, increase 10 to, e.g., 20
////energy(VectorXd q){ return q.transpose() * Dij * q + Dii.dot(q); }; //when n=5, q.size is 10 and Dij is 10 x 10 and Dii.size is 10 (where 10 is n(n-1)/2)
////objective(VectorXd P, MatrixXd Q) { return P.transpose() * Q * P; } //when n=5, P.size is 25 and Q is 25 x 25; summary energy() and objective() return different values

        for (int j = 0; j < maxiter; ++j) {
            VectorXd Px0; MatrixXd dPx0;
            tie(Px0, dPx0) = eval_p(T, x0);
            MatrixXd Pc = Px0 - dPx0 * x0;

            MatrixXd Dij = dPx0.transpose() * Q_mod * dPx0;
            VectorXd Dii = (2 * Pc.transpose() * Q_mod * dPx0).transpose();
			//size of q_free = size of Dii = n * (n - 1) / 2 if no variable of q are set to fixed values, e.g. fixed_vars_1.empty(), o/w it is of shorter length

			if (j == 0) //for the initialization in SA()/simulatedAnnealing()
			{
				qPrev.resize( Dii.size() );
				for (int i = 0; i < Dii.size(); ++i)
					qPrev(i) = (int) (rand() % 2);
//cout << Dii << "\n\n" << Dij << "\n\n";
			}//*/

            auto reduced = fix_variables(Dij, Dii, fixed_vars_1, fixed_vars_0);

//			VectorXd q_free = SA(reduced.Dii, reduced.Dij); //simulated annealing gives the best q_free based on the current reduced.Dii and Dij
			VectorXd q_free = simulatedAnnealing(reduced.Dij, reduced.Dii); //better alternative w/o any temparature problem (-delta / temp in SA() was problematic as delta~7000 temp~2)
			qPrev = q_free; //for the initialization in SA()/simulatedAnnealing() (not a good idea 'cos multistarts are supposed to start randomly at each iteration)

            // --- Reconstruct full q vector
            VectorXd q_full = VectorXd::Zero(n_trans);
            set<int> fixed_set(fixed_vars_1.begin(), fixed_vars_1.end());
			fixed_set.insert(fixed_vars_0.begin(), fixed_vars_0.end());

			vector<int> free_idx;
            for (int i = 0; i < n_trans; ++i)
                if (! fixed_set.count(i))
                    free_idx.push_back(i);

			for (size_t k = 0; k < free_idx.size(); ++k)                
				q_full(free_idx[k]) = q_free(k);
            for (int idx : fixed_vars_1)
                q_full(idx) = 1;
            for (int idx : fixed_vars_0)
                q_full(idx) = 0;

            x0 = q_full;

            VectorXd P_vec; MatrixXd dP_dummy;
            tie(P_vec, dP_dummy) = eval_p(T, x0);
//			P_sol = Map<MatrixXd>(P_vec.data(), n, n); //old code
			P_sol = Map<MatrixXd>(P_vec.data(), n, n).transpose();
            sol = objective(P_vec, Q);
			if (print)
				cout << "\nIter " << j << " | Objective = " << sol << " | Norm(P)^2 = " << P_sol.squaredNorm() << endl;

			if (fabs(sol - sol_prev) < 1e-8)
				break;
            sol_prev = sol;
        }
		if (print)
			cout << "\nFinal Permutation Matrix:\n" << P_sol << endl;
    }

	void identityMap() {
        VectorXd x0 = VectorXd::Zero(n_trans); //means identity permutation matrix for the first j=0 iteration (may be a very bad initialization so consider multi-starts)
        
		VectorXd P_vec; MatrixXd dP_dummy;
        tie(P_vec, dP_dummy) = eval_p(T, x0);
//			P_sol = Map<MatrixXd>(P_vec.data(), n, n); //old code
		P_sol = Map<MatrixXd>(P_vec.data(), n, n).transpose();
        sol = objective(P_vec, Q);
        cout << "\nIdentityMap | Objective = " << sol << " | Norm(P)^2 = " << P_sol.squaredNorm() << endl;
        cout << "\nFinal Permutation Matrix:\n" << P_sol << endl; //identity for sure (since x0 above is all 0, i.e., 0 swap applied to the identity permutation matrix)
    }
};

#endif
