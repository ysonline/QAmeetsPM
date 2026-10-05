import numpy as np
import matplotlib.pyplot as plt
import dimod
from dwave.system import DWaveSampler, LazyFixedEmbeddingComposite
import neal
import utils
import logging

np.set_printoptions(linewidth=100)

class QuCM:
    def __init__(self, Q, fixed_assignments=None):
        self.Q = Q
        self.n = int(np.sqrt(Q.shape[0]))
        self.T = utils.t(self.n)
        self.n_trans = self.T.shape[1]
        self.P_sol = np.eye(self.n)


    def solve(self, maxiter=15, simulated_anneal=True):
        Q = self.Q - 10*np.min(np.linalg.eigvalsh(self.Q)) * np.eye(self.n*self.n)  # Q-matrix with penalty term

        I = np.eye(self.n)
        P_sol = I
        sol = utils.objective(P_sol, self.Q)
        history = [sol, ]

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
            I_flatt = I.reshape(self.n*self.n, 1)
            Pj = np.kron(I, P_sol.T)
            Qj = Pj.T @ Q @ Pj

            Aj = self.T - I_flatt

            Dij = Aj.T @ Qj @ Aj
            Dii = (2 * I_flatt.T @ Qj @ Aj).flatten()

            W = {(i, j): Dij[i, j] for i in range(self.n_trans) for j in range(self.n_trans)}
            c = {i: Dii[i] for i in range(self.n_trans)}
            bqm = dimod.BinaryQuadraticModel(c, W, vartype=dimod.BINARY)

            if simulated_anneal:
                response = sampler.sample(bqm, num_reads=100, num_sweeps=2000)
            else:
                response = sampler.sample(bqm, num_reads=50)
            # dwave.inspector.show(response)  # If quantum annealing and if needed

            best = response.first.sample

            q = np.array(list(best.values())).reshape(self.n_trans, 1)

            P_sol = (I_flatt + Aj @ q).reshape(self.n, self.n) @ P_sol
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
    qap = QuCM(Q, fixed_assignments=[])
    history = qap.solve()

    # ResultsFLFalse
    # ResultsFLFalse
    print(f'\nResults:            {(qap.P_sol @ np.arange(n)).astype(int)}, objective: {utils.objective(qap.P_sol, Q)}'
          f'\nGround through:     {(P_star @ np.arange(n)).astype(int)}, objective: {utils.objective(P_star, Q)}.')
    plt.plot(history, label='Objective')
    plt.hlines(utils.objective(P_star, Q), 0, len(history) - 1, linestyles='--', color='r', label='Objective target')
    plt.legend()
    plt.show()
