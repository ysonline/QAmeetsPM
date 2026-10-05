import numpy as np
import matplotlib.pyplot as plt
import dimod
from dwave.system import DWaveSampler, LazyFixedEmbeddingComposite
import neal
import utils
import logging

np.set_printoptions(linewidth=100)

class QuCOOP:
    def __init__(self, Q, fixed_assignments=None):
        self.Q = Q
        self.n = int(np.sqrt(Q.shape[0]))
        self.T = utils.t(self.n)
        self.n_trans = self.T.shape[1]
        self.P_sol = np.eye(self.n)
        self.fixed_variables = self.get_fixed_variables(fixed_assignments)

    def get_fixed_variables(self, fixed_assignments):
        fixed_assignments = [] if fixed_assignments is None else fixed_assignments
        fixed_assignments = [(np.minimum(*pair), np.maximum(*pair)) for pair in fixed_assignments]
        assignments = [(i, j) for i in range(self.n) for j in range(i + 1, self.n)]

        # Step 1: Get indices of fixed assignments
        fixed_variables = [assignments.index(pair) for pair in fixed_assignments if pair in assignments]

        # Step 2: Extract all indices involved in fixed assignments
        fixed_indices = set(i for pair in fixed_assignments for i in pair)

        # Step 3: Find all assignment pairs that involve any fixed index but are not already fixed
        related_variables = [
            idx for idx, pair in enumerate(assignments)
            if (pair[0] in fixed_indices or pair[1] in fixed_indices) and idx not in fixed_variables
        ]
        return fixed_variables, related_variables

    def solve(self, maxiter=15, simulated_anneal=True):
        x0 = np.zeros(self.n_trans).reshape(self.n_trans, 1)

        Q = self.Q - 10*np.min(np.linalg.eigvalsh(self.Q)) * np.eye(self.n*self.n)  # Q-matrix with penalty term

        P_sol = utils.eval_p(self.T, x0)[0].reshape(self.n, self.n)
        sol = utils.objective(P_sol, self.Q)
        history = [sol,]

        # Do nothing if Q is zero
        if np.array_equal(Q, np.zeros_like(Q)):
            return P_sol, history

        # Select solver and sampler
        if simulated_anneal:
            sampler = neal.SimulatedAnnealingSampler()           # Simulated annealing
        else:
            qpu_advantage = DWaveSampler(solver={'topology__type': 'pegasus'})
            sampler = LazyFixedEmbeddingComposite(qpu_advantage)  # Quantum annealing

        for j in range(maxiter):
            Px0, dPx0 = utils.eval_p(self.T, x0)
            Pc = Px0 - dPx0 @ x0

            Dij = dPx0.T @ Q @ dPx0
            Dii = (2 * Pc.T @ Q @ dPx0).flatten()

            W = {(i, j): Dij[i, j] for i in range(self.n_trans) for j in range(self.n_trans)}
            c = {i: Dii[i] for i in range(self.n_trans)}
            bqm = dimod.BinaryQuadraticModel(c, W, vartype=dimod.BINARY)
            for idx in self.fixed_variables[0]:
                bqm.fix_variable(idx, 1)
            for idx in self.fixed_variables[1]:
                bqm.fix_variable(idx, 0)

            if simulated_anneal:
                response = sampler.sample(bqm, num_reads=100, num_sweeps=2000)
            else:
                response = sampler.sample(bqm, num_reads=50)
            # dwave.inspector.show(response)  # If quantum annealing and if needed

            best = response.first.sample

            q_free = np.array(list(best.values()))
            q_full = np.zeros(self.n_trans, dtype=int)
            fixed_one = set(self.fixed_variables[0])
            fixed_zero = set(self.fixed_variables[1])
            all_fixed = fixed_one | fixed_zero
            free_indices = [i for i in range(self.n_trans) if i not in all_fixed]
            for idx, val in zip(free_indices, q_free):
                q_full[idx] = val
            for idx in fixed_one:
                q_full[idx] = 1
            for idx in fixed_zero:
                q_full[idx] = 0

            q = q_full.reshape(self.n_trans, 1)

            x0 = q
            # x0 = q + .5 * (q - x0)
            P_sol = utils.eval_p(self.T, x0)[0].reshape(self.n, self.n)
            sol = utils.objective(P_sol, self.Q)
            history.append(sol)

            logging.info(f"\nIter {j}"
                         f"\nPermutation matrix found: \nP = \n{P_sol}"
                         f"\nObjective = {sol}, norm P = {int(np.linalg.norm(P_sol) ** 2)}, ")

            if history[-1] == history[-2]:
                break

        self.P_sol = P_sol
        return history


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(message)s')

    # Create random problem
    seed = 2
    np.random.seed(seed)
    n = 10
    P = np.eye(n)
    np.random.shuffle(P)
    P_star = P

    ref_points = np.random.random((2, n))
    tmp_points = ref_points @ P_star

    A = np.array([[np.linalg.norm(ref_points[:, i] - ref_points[:, j]) for i in range(n)] for j in range(n)])
    B = np.array([[np.linalg.norm(tmp_points[:, i] - tmp_points[:, j]) for i in range(n)] for j in range(n)])
    Q = -np.kron(A, B)


    # Set seed back to None
    np.random.seed(None)

    # Start Solver
    qap = QuCOOP(Q, fixed_assignments=[])
    history = qap.solve()

    # ResultsFLFalse
    # ResultsFLFalse
    print(f'\nResults:            {(qap.P_sol @ np.arange(n)).astype(int)}, objective: {utils.objective(qap.P_sol, Q)}'
          f'\nGround through:     {(P_star @ np.arange(n)).astype(int)}, objective: {utils.objective(P_star, Q)}.')
    plt.plot(history, label='Objective')
    plt.hlines(utils.objective(P_star, Q), 0, len(history) - 1, linestyles='--', color='r', label='Objective target')
    plt.legend()
    plt.show()
