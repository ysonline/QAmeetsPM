import numpy as np
import scipy as sp

import polyscope as ps


# Functions for permutation matrices

def p(col_ind):
    col_ind = col_ind.flatten()
    n = len(col_ind)
    P = np.zeros((n, n))
    for i in range(n):
        P[i, col_ind[i]] = 1
    return P


def permutation_idx2mat(col_ind):
    col_ind = col_ind.flatten()
    n = len(col_ind)
    P = np.zeros((n, n))
    for i in range(n):
        P[i, col_ind[i]] = 1
    return P


def t(n):
    T = []
    for i in range(n):
        for j in range(i + 1, n):
            Ti = np.eye(n)
            Ti[i, i] = Ti[j, j] = 0
            Ti[i, j] = Ti[j, i] = 1
            T.append(Ti.reshape(n*n, 1))
    T = np.hstack(T)
    return T


def eval_p(T, x, P_rec=None):
    n = int(np.sqrt(T.shape[0]))
    n_trans = T.shape[1]
    I = np.eye(n)
    P_rec = I if P_rec is None else P_rec

    P = P_rec
    dP = np.zeros_like(T)
    for i in range(n_trans):
        P = P @ (x[i, 0] * (T[:, i].reshape(n, n) - I) + I)

        # Eval dP
        dp = P_rec
        for j in range(i):
            dp = dp @ (x[j, 0] * (T[:, j].reshape(n, n) - I) + I)
        dp = dp @ (T[:, i] - I.flatten()).reshape(n, n)
        for j in range(i+1, n_trans):
            dp = dp @ (x[j, 0] * (T[:, j].reshape(n, n) - I) + I)
        dP[:, i] = dp.flatten()

    P = P.reshape(n*n, 1)
    return P, dP


def objective(P, Q):
    P = P.reshape(-1)[:, None]
    return (P.T @ Q @ P).squeeze()


def transverse_field(n):
    sigma_x = np.array([[0, 1], [1, 0]])
    sigma_x_i = lambda i: sp.sparse.kron(sp.sparse.eye(2 ** i), sp.sparse.kron(sigma_x, sp.sparse.eye(2 ** (n - i - 1))))

    h0 = sp.sparse.csr_matrix((2**n, 2**n))
    for i in range(n):
        h0 += sigma_x_i(i)
    return h0


def problem_hamiltonian(couplers, biases, n):
    biases = biases.flatten()
    sigma_z = np.array([[1, 0], [0, -1]])
    sigma_z_i = lambda i: sp.sparse.kron(sp.sparse.eye(2 ** i),
                                         sp.sparse.kron(sigma_z, sp.sparse.eye(2 ** (n - i - 1))))
    sigma_z_ij = lambda i, j: sp.sparse.eye(2 ** n) if i == j \
        else sp.sparse.kron(sp.sparse.eye(2 ** i),
                            sp.sparse.kron(sigma_z,
                                           sp.sparse.kron(sp.sparse.eye(2 ** (j - i - 1)),
                                                          sp.sparse.kron(sigma_z, sp.sparse.eye(2 ** (n - j - 1))))))

    h1 = sp.sparse.csr_matrix((2 ** n, 2 ** n))
    for i in np.arange(n):
        h1 += biases[i] * sigma_z_i(i)
        for j in np.arange(n):
            h1 += couplers[i, j] * sigma_z_ij(np.minimum(i, j), np.maximum(i, j))
    return h1


def evaluateCorrespondences(C, Xgeodesics, Ygeodesics):
    def detect_outliers(scores):
        scores_sorted = sorted(scores)
        diffs = [scores_sorted[i + 1] - scores_sorted[i] for i in range(len(scores_sorted) - 1)]

        if not diffs:
            return []
        base_diff = diffs[0]
        threshold = base_diff

        for i, diff in enumerate(diffs):
            if diff > threshold:
                return scores_sorted[i + 1:]  # Outliers start here
        return []  # No outliers found

    scores = np.zeros(C.shape[0])

    # calc score for each correspondence
    for i in range(C.shape[0]):
        for j in range(i+1, C.shape[0]):
            scores[i] += np.abs(Xgeodesics[C[i, 0], C[j, 0]] - Ygeodesics[C[i, 1], C[j, 1]])
            scores[j] += np.abs(Xgeodesics[C[i, 0], C[j, 0]] - Ygeodesics[C[i, 1], C[j, 1]])

    # get worst matches
    bad_idx = detect_outliers(scores)

    return bad_idx, scores


def visualize(meshS, meshT, matching, bad_idx=[]):
    # VISUALIZATION VISUALIZATION
    ps.init()

    # Offset for Mesh T
    offset = np.array([-100., 0, 0])
    vertsT = np.array(meshT.vertices)
    vertsS = np.array(meshS.vertices) + offset
    facesS = np.array(meshS.faces)
    facesT = np.array(meshT.faces)

    # Register both meshes
    ps.register_surface_mesh("Mesh S", vertsS, facesS, smooth_shade=True)
    ps.register_surface_mesh("Mesh T", vertsT, facesT, smooth_shade=True)

    # Create line segments between matching vertices
    curve_points = []
    curve_edges = []
    sphere_points = []
    bad_sphere_points = []

    for idx, (i, j) in enumerate(matching):
        p1 = vertsS[meshS.samples[i]]
        p2 = vertsT[meshT.samples[j]]

        # Add to line network
        curve_points.extend([p1, p2])
        curve_edges.append([2 * idx, 2 * idx + 1])

        # Add sphere centers
        if i in bad_idx:
            bad_sphere_points.extend([p1, p2])
        else:
            sphere_points.extend([p1, p2])

    curve_points = np.array(curve_points)
    curve_edges = np.array(curve_edges)
    sphere_points = np.array(sphere_points)
    bad_sphere_points = np.array(bad_sphere_points)

    # Register the correspondence lines
    ps.register_curve_network("Correspondences", curve_points, curve_edges, color=(0.2, 0.9, 0.2), radius=0.003)

    # Register spheres at matched vertices
    # if len(sphere_points) > 0:
    #     ps.register_point_cloud("Match Points", sphere_points, radius=0.01, color=(0.3, 1.0, 0.3))
    # if len(bad_sphere_points) > 0:
    #     ps.register_point_cloud("Bad Match Points", bad_sphere_points, radius=0.01, color=(1., 0.3, 0.3))

    ps.register_point_cloud("Bad Match Points S", vertsS[meshS.samples], radius=0.01, color=(1., 0.3, 0.3))
    ps.register_point_cloud("Bad Match Points T", vertsT[meshT.samples], radius=0.01, color=(1., 0.3, 0.3))

    # Show the scene
    ps.show()

    return
