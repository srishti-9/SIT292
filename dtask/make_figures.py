"""Generate the figures used in the report (saved to report/figures/)."""

from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from pagerank import build_transition_matrix, google_matrix, pagerank_power  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "report" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

ADJ = np.array(
    [[0, 1, 1, 0], [0, 0, 1, 0], [1, 0, 0, 0], [1, 1, 1, 0]], dtype=float
)
NODES = ["A", "B", "C", "D"]
NAVY = "#142a5a"

plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})


def G_for(adj, d):
    return google_matrix(build_transition_matrix(adj), d)


def convergence():
    """L1 error of power iteration vs k, against the d^k bound."""
    fig, ax = plt.subplots(figsize=(3.4, 2.5))
    for d, style in [(0.5, "-"), (0.85, "-"), (0.95, "-")]:
        G = G_for(ADJ, d)
        r_star = pagerank_power(G, tol=1e-15, max_iter=2000)[0]
        _, _, hist = pagerank_power(G, tol=1e-15, max_iter=200, return_history=True)
        err = np.abs(hist - r_star).sum(axis=1)
        ax.semilogy(err, style, label=f"$d={d}$", lw=1.2)
        ax.semilogy(d ** np.arange(len(err)), ":", color="gray", lw=0.7)
    ax.set_ylim(1e-12, 2)
    ax.set_xlim(0, 100)
    ax.set_xlabel("iteration $k$")
    ax.set_ylabel(r"$\|\mathbf{r}^{(k)}-\mathbf{r}\|_1$")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "convergence.pdf")
    plt.close(fig)


def damping_sweep():
    """PageRank of each page as d varies."""
    ds = np.linspace(0.01, 0.99, 99)
    ranks = np.array([pagerank_power(G_for(ADJ, d))[0] for d in ds])
    fig, ax = plt.subplots(figsize=(3.4, 2.5))
    for j, name in enumerate(NODES):
        ax.plot(ds, ranks[:, j], label=name, lw=1.3)
    ax.axvline(0.85, color="gray", ls=":", lw=0.8)
    ax.set_xlabel("damping factor $d$")
    ax.set_ylabel("PageRank")
    ax.legend(frameon=False, ncol=4, loc="upper center")
    ax.set_ylim(0, 0.5)
    fig.tight_layout()
    fig.savefig(OUT / "damping_sweep.pdf")
    plt.close(fig)


def link_changes():
    """Grouped bars: baseline vs added link vs removed link."""
    base = pagerank_power(G_for(ADJ, 0.85))[0]
    add = ADJ.copy()
    add[0, 3] = 1
    rem = ADJ.copy()
    rem[2, 0] = 0
    r_add = pagerank_power(G_for(add, 0.85))[0]
    r_rem = pagerank_power(G_for(rem, 0.85))[0]
    x = np.arange(4)
    w = 0.26
    fig, ax = plt.subplots(figsize=(3.4, 2.5))
    ax.bar(x - w, base, w, label="baseline", color=NAVY)
    ax.bar(x, r_add, w, label=r"add $A\to D$", color="#5b8dd6")
    ax.bar(x + w, r_rem, w, label=r"remove $C\to A$", color="#c0a060")
    ax.set_xticks(x, NODES)
    ax.set_ylabel("PageRank")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "link_changes.pdf")
    plt.close(fig)


if __name__ == "__main__":
    convergence()
    damping_sweep()
    link_changes()
    print("Figures written to", OUT)