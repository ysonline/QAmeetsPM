import numpy as np
import dimod
from dwave.system import DWaveSampler, LazyFixedEmbeddingComposite
import neal
import utils
import logging

np.set_printoptions(linewidth=100)

class QGM:
    def __init__(self, Q):
        self.Q = Q
        self.n = int(np.sqrt(Q.shape[0]))
        self.P_sol = np.eye(self.n)
        self.method = 'baseline'


    def get_bqm(self):

        C = np.diag(self.Q)
        W = self.Q - np.diag(C)

        # Prepare the permutation matrix-constraint
        A1 = np.kron(np.eye(self.n), np.ones((1, self.n)))
        A2 = np.kron(np.ones((1, self.n)), np.eye(self.n))
        A = np.concatenate((A1, A2), axis=0)

        b = np.ones((self.n * 2, 1))

        if self.method == 'baseline':

            # BASELINE weights
            Lambda_baseline = np.sum(np.abs(W)) / 2
            Q = W + Lambda_baseline * np.matmul(A.T, A)
            q = C.T - 2 * Lambda_baseline * np.matmul(b.T, A)

        if self.method == 'row_wise':
            # ROWWISE weights
            optimizing = np.abs(np.sum(W, axis=0)) + np.abs(np.sum(W, axis=1)) - np.abs(np.diag(W))
            MaxGrad = np.max(optimizing)
            Lambda_row_wise = np.max(A * optimizing, axis=1) + 1 / 2 * MaxGrad
            Q = W + np.matmul((Lambda_row_wise[:, np.newaxis] * A).T, A)
            q = C.T - 2 * np.matmul(b.T, Lambda_row_wise[:, np.newaxis] * A)


        J = {(i, j): Q[i, j] for i in range(self.n ** 2) for j in range(self.n ** 2)}
        bias = {i: qi for i, qi in enumerate(q.flatten())}
        bqm = dimod.BinaryQuadraticModel(bias, J, vartype=dimod.BINARY)
        return bqm


    def solve(self, simulated_anneal=True):

        P_sol = np.eye(self.n)
        sol = utils.objective(P_sol, self.Q)
        history = [sol,]

        # Select solver and sampler
        if simulated_anneal:
            sampler = neal.SimulatedAnnealingSampler()  # Simulated annealing
        else:
            qpu_advantage = DWaveSampler(solver={'topology__type': 'pegasus'})
            sampler = LazyFixedEmbeddingComposite(qpu_advantage)  # Quantum annealing

        bqm = self.get_bqm()
        num_reads = 500 if simulated_anneal else 50
        response = sampler.sample(bqm, num_reads=num_reads)

        # dwave.inspector.show(response)
        best = response.first.sample
        self.P_sol = np.array(list(best.values())).reshape(self.n, self.n)


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
    qap = QGM(Q)
    history = qap.solve()

    # ResultsFLFalse
    print(f'\nResults:            {(qap.P_sol @ np.arange(n)).astype(int)}, objective: {utils.objective(qap.P_sol, Q)}'
          f'\nGround through:     {(P_star @ np.arange(n)).astype(int)}, objective: {utils.objective(P_star, Q)}.')
