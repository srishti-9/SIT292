"""PageRank via power iteration and direct eigendecomposition."""

from .core import (
    build_transition_matrix,
    google_matrix,
    pagerank_eigen,
    pagerank_linear_solve,
    pagerank_power,
)

__all__ = [
    "build_transition_matrix",
    "google_matrix",
    "pagerank_power",
    "pagerank_eigen",
    "pagerank_linear_solve",
]
__version__ = "1.0.0"