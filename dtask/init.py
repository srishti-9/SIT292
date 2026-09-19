import numpy as np
import pytest

from pagerank import (
    build_transition_matrix,
    google_matrix,
    pagerank_eigen,
    pagerank_linear_solve,
    pagerank_power,
)

ADJ = np.array(
    [[0, 1, 1, 0], [0, 0, 1, 0], [1, 0, 0, 0], [1, 1, 1, 0]], dtype=float
)


def rank(adj, d=0.85):
    return pagerank_power(google_matrix(build_transition_matrix(adj), d))[0]


def test_H_is_row_stochastic():
    assert np.allclose(build_transition_matrix(ADJ).sum(axis=1), 1.0)


def test_dangling_row_becomes_uniform():
    adj = ADJ.copy()
    adj[2, 0] = 0
    assert np.allclose(build_transition_matrix(adj)[2], 0.25)


def test_G_positive_and_row_stochastic():
    G = google_matrix(build_transition_matrix(ADJ))
    assert (G >= 0.15 / 4 - 1e-12).all()
    assert np.allclose(G.sum(axis=1), 1.0)


def test_report_baseline_vector():
    assert np.allclose(rank(ADJ), [0.3732, 0.2068, 0.3825, 0.0375], atol=5e-5)


def test_three_methods_agree():
    H = build_transition_matrix(ADJ)
    G = google_matrix(H)
    r = pagerank_power(G)[0]
    assert np.allclose(r, pagerank_eigen(G), atol=1e-8)
    assert np.allclose(r, pagerank_linear_solve(H), atol=1e-8)


def test_probability_vector_and_fixed_point():
    G = google_matrix(build_transition_matrix(ADJ))
    r = pagerank_power(G)[0]
    assert r.min() > 0 and np.isclose(r.sum(), 1.0)
    assert np.allclose(r @ G, r, atol=1e-9)


def test_add_link_A_to_D():
    adj = ADJ.copy()
    adj[0, 3] = 1
    assert np.allclose(rank(adj), [0.3558, 0.1775, 0.3284, 0.1383], atol=5e-5)


def test_remove_link_C_to_A():
    adj = ADJ.copy()
    adj[2, 0] = 0
    assert np.allclose(rank(adj), [0.1712, 0.2440, 0.4514, 0.1334], atol=5e-5)


@pytest.mark.parametrize(
    "d, expected",
    [
        (0.50, [0.3141, 0.2244, 0.3365, 0.1250]),
        (0.85, [0.3732, 0.2068, 0.3825, 0.0375]),
        (0.95, [0.3910, 0.2022, 0.3943, 0.0125]),
    ],
)
def test_damping_table(d, expected):
    assert np.allclose(rank(ADJ, d), expected, atol=5e-5)


@pytest.mark.parametrize("d", [0.3, 0.5, 0.85, 0.95])
def test_D_equals_teleportation_floor(d):
    assert np.isclose(rank(ADJ, d)[3], (1 - d) / 4, atol=1e-8)


def test_second_eigenvalue_bounded_by_d():
    for d in (0.5, 0.85, 0.95):
        G = google_matrix(build_transition_matrix(ADJ), d)
        mags = np.sort(np.abs(np.linalg.eigvals(G)))[::-1]
        assert np.isclose(mags[0], 1.0)
        assert mags[1] <= d + 1e-9


def test_invalid_inputs():
    with pytest.raises(ValueError):
        google_matrix(np.eye(3), d=1.0)
    with pytest.raises(ValueError):
        build_transition_matrix(np.ones((2, 3)))