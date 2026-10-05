import numpy as np
import dimod
import logging
from dwave.system import DWaveSampler, LazyFixedEmbeddingComposite
import neal
import utils


class QMatch:
    def __init__(self, Q):
        self.Q = Q
        self.n = int(np.sqrt(Q.shape[0]))
        self.P_sol = np.eye(self.n)

    def table_to_vector(self, table):
        vector = np.zeros(self.n ** 2)
        for count, value in enumerate(table):
            vector[self.n * count + value] = 1
        return vector

    def get_couplings(self, current_perm, new_perm):
        dim = len(new_perm)
        Wneu = np.zeros((dim, dim))
        cneu = np.zeros((dim))
        constant = self.table_to_vector(current_perm)

        for i, cycle in enumerate(new_perm):
            for j, cycle2 in enumerate(new_perm):
                v1, v2 = np.zeros((self.n ** 2)), np.zeros((self.n ** 2))
                for k in range(len(cycle)):
                    idx = current_perm.index(cycle[k])
                    v1[self.n * idx + cycle[k]] = -1
                    v1[self.n * idx + cycle[(k + 1) % len(cycle)]] = 1
                for k in range(len(cycle2)):
                    idx = current_perm.index(cycle2[k])
                    v2[self.n * idx + cycle2[k]] = -1
                    v2[self.n * idx + cycle2[(k + 1) % len(cycle2)]] = 1
                Wneu[i, j] = v1.T @ self.Q @ v2
                if i == 0:
                    cneu[j] += constant.T @ self.Q @ v2
            cneu[i] += v1.T @ self.Q @ constant

        return Wneu, cneu

    def update_table(self, table, cycles, decisions):
        result = table.copy()
        for i, cycle in enumerate(cycles):
            if decisions[i] == 1:
                indices = [table.index(x) for x in cycle]
                for j, idx in enumerate(indices):
                    result[idx] = cycle[(j + 1) % len(cycle)]
        return result

    def anneal(self, Q, q, simulated_anneal=True):
        J = {(i, j): Q[i, j] for i in range(Q.shape[0]) for j in range(Q.shape[0])}
        bias = {i: qi for i, qi in enumerate(q.flatten())}
        bqm = dimod.BinaryQuadraticModel(bias, J, vartype=dimod.BINARY)

        if simulated_anneal:
            sampler = neal.SimulatedAnnealingSampler()  # Simulated annealing
        else:
            qpu_advantage = DWaveSampler(solver={'topology__type': 'pegasus'})
            sampler = LazyFixedEmbeddingComposite(qpu_advantage)  # Quantum annealing


        if simulated_anneal:
            response = sampler.sample(bqm, num_reads=100, num_sweeps=2000)
        else:
            response = sampler.sample(bqm, num_reads=50)

        best = response.first.sample

        return best

    def next_pairs(self, pairs):
        next_pairs = []

        next_pairs.append((pairs[0][0], pairs[1][1]))
        if len(pairs) > 2:
            next_pairs.append((pairs[0][1], pairs[2][1]))
            for k in range(2, len(pairs) - 1):
                next_pairs.append((pairs[k - 1][0], pairs[k + 1][1]))
            next_pairs.append((pairs[-2][0], pairs[-1][0]))
        else:
            next_pairs.append((pairs[0][1], pairs[1][0]))
        return next_pairs

    def optimize(self, initial_perm=None, cycles=None, verbose=True):
        if initial_perm is None:
            perm_table = list(range(self.n))
        else:
            perm_table = list(initial_perm)

        cycle_lists = cycles or []
        init_pairs = [[i, self.n - i - 1] for i in range(self.n // 2)]
        np.random.shuffle(init_pairs)
        init_pairs = np.array(init_pairs).reshape(-1, 2).tolist()
        cycle_lists.append(init_pairs)

        current_pairs = self.next_pairs(init_pairs)
        for _ in range(self.n - 2):
            cycle_lists.append(current_pairs)
            current_pairs = self.next_pairs(current_pairs)

        for cycles in cycle_lists:
            Wq, bq = self.get_couplings(perm_table, cycles)
            result = self.anneal(Wq, bq)
            perm_table = self.update_table(perm_table, cycles, result)
            energy = self.table_to_vector(perm_table).T @ self.Q @ self.table_to_vector(perm_table)
            if verbose:
                print(f"Updated energy: {energy}")

        return np.array(perm_table)

    def solve(self, number_of_iterations=20, verbose=False):

        P_sol = np.eye(self.n)
        sol = utils.objective(P_sol, self.Q)
        history = [sol,]

        permutations = []
        current_perm = None

        for iteration in range(number_of_iterations):
            current_perm = self.optimize(initial_perm=current_perm, verbose=verbose)
            permutations.append(current_perm)

            self.P_sol = utils.p(permutations[-1])
            sol = utils.objective(self.P_sol, self.Q)
            history.append(sol)

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
    qap = QMatch(Q)
    history = qap.solve()

    # ResultsFLFalse
    print(f'\nResults:            {(qap.P_sol @ np.arange(n)).astype(int)}, objective: {utils.objective(qap.P_sol, Q)}'
          f'\nGround through:     {(P_star @ np.arange(n)).astype(int)}, objective: {utils.objective(P_star, Q)}.')