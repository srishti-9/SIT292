"""Core PageRank routines (row-stochastic convention: r G = r).

Notation follows the report:
    Adj : adjacency matrix, Adj[i, j] = 1 if page i links to page j
    H   : row-stochastic transition matrix (dangling rows -> uniform 1/n)
    G   : Google matrix, G = d H + (1 - d) (1/n) 1 1^T
"""

from __future__ import annotations

import numpy as np


def build_transition_matrix(adj: np.ndarray) -> np.ndarray:
    """Turn an adjacency matrix into a row-stochastic matrix H.

    Rows with no outgoing links (dangling nodes) are replaced by the
    uniform distribution 1/n.
    """
    adj = np.asarray(adj, dtype=float)
    if adj.ndim != 2 or adj.shape[0] != adj.shape[1]:
        raise ValueError("adjacency matrix must be square")
    n = adj.shape[0]
    outdeg = adj.sum(axis=1)
    H = np.empty_like(adj)
    for i in range(n):
        H[i] = 1.0 / n if outdeg[i] == 0 else adj[i] / outdeg[i]
    return H


def google_matrix(H: np.ndarray, d: float = 0.85) -> np.ndarray:
    """Build G = d H + (1 - d) S with S the uniform teleportation matrix."""
    if not 0.0 < d < 1.0:
        raise ValueError("damping factor d must satisfy 0 < d < 1")
    n = H.shape[0]
    return d * H + (1.0 - d) * np.ones((n, n)) / n


def pagerank_power(
    G: np.ndarray,
    tol: float = 1e-10,
    max_iter: int = 200,
    return_history: bool = False,
):
    """Power iteration r <- r G, stopping when ||r_new - r||_1 < tol.

    Returns (r, iterations), plus the array of iterates if return_history.
    """
    n = G.shape[0]
    r = np.full(n, 1.0 / n)
    history = [r.copy()]
    k = 0
    for k in range(1, max_iter + 1):
        r_next = r @ G
        history.append(r_next.copy())
        done = np.linalg.norm(r_next - r, ord=1) < tol
        r = r_next
        if done:
            break
    if return_history:
        return r, k, np.array(history)
    return r, k


def pagerank_eigen(G: np.ndarray) -> np.ndarray:
    """PageRank as the eigenvector of G^T for eigenvalue 1, normalised."""
    eigvals, eigvecs = np.linalg.eig(G.T)
    idx = np.argmin(np.abs(eigvals - 1.0))
    r = np.real(eigvecs[:, idx])
    return r / r.sum()


def pagerank_linear_solve(H: np.ndarray, d: float = 0.85) -> np.ndarray:
    """Solve r (I - d H) = (1 - d) u directly, u uniform."""
    n = H.shape[0]
    u = np.full(n, (1.0 - d) / n)
    # r (I - dH) = u  <=>  (I - dH)^T r^T = u^T
    return np.linalg.solve((np.eye(n) - d * H).T, u)