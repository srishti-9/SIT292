"""Reproduce every number in the report's four-page example."""

import numpy as np

from pagerank import (
    build_transition_matrix,
    google_matrix,
    pagerank_eigen,
    pagerank_linear_solve,
    pagerank_power,
)

np.set_printoptions(precision=4, suppress=True)

ADJ = np.array(
    [
        [0, 1, 1, 0],  # A -> B, C
        [0, 0, 1, 0],  # B -> C
        [1, 0, 0, 0],  # C -> A
        [1, 1, 1, 0],  # D -> A, B, C
    ],
    dtype=float,
)


def rank(adj, d=0.85):
    G = google_matrix(build_transition_matrix(adj), d)
    return pagerank_power(G)[0]


def main():
    H = build_transition_matrix(ADJ)
    G = google_matrix(H, 0.85)
    r, iters, hist = pagerank_power(G, return_history=True)

    print("H =\n", H)
    print("G =\n", G)
    print(f"\nPower iteration converged in {iters} iterations")
    print("Power iteration :", r)
    print("Eigenvector     :", pagerank_eigen(G))
    print("Linear solve    :", pagerank_linear_solve(H, 0.85))

    print("\nIteration table (k, r_A..r_D):")
    for k in (0, 1, 2, 3, 4, 5, 10, 20):
        print(f"{k:>3}", hist[k])

    base = rank(ADJ)

    add = ADJ.copy()
    add[0, 3] = 1  # A -> D
    r_add = rank(add)
    print("\nAdd A->D   :", r_add, " change:", r_add - base)

    rem = ADJ.copy()
    rem[2, 0] = 0  # remove C -> A (C becomes dangling)
    r_rem = rank(rem)
    print("Remove C->A:", r_rem, " change:", r_rem - base)

    print("\nDamping sweep:")
    for d in (0.50, 0.85, 0.95):
        print(f"d={d:.2f}", rank(ADJ, d))


if __name__ == "__main__":
    main()