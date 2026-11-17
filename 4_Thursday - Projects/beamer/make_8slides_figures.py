"""Figures for Project3_8slides.tex that differ from the Part 1 notebook:
- SIRS cartoon labelled with rates only (beta, gamma, omega), not propensities
- the three "8 SSA runs vs ODE" panels as separate, fully labelled images
Same code, parameters and random seed (2026) as Project3_Part1_SIRS_SSA_vs_ODE.ipynb.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.integrate import solve_ivp

plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})

beta, gamma, omega = 0.01, 0.15, 0.025
T_MAX = 365.0


# ---------------- Cartoon: rates only ----------------
fig, ax = plt.subplots(figsize=(8, 3.6))
ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 4.2)
pos = {"S": (1.2, 1.2), "I": (5, 3.0), "R": (8.8, 1.2)}
col = {"S": "#4C6A86", "I": "#BF4E2C", "R": "#228066"}
for k, (x, y) in pos.items():
    ax.add_patch(plt.Rectangle((x-0.8, y-0.45), 1.6, 0.9, fc="white", ec=col[k], lw=2.5))
    ax.text(x, y, k, ha="center", va="center", fontsize=18, fontweight="bold", color=col[k])
arrow = dict(arrowstyle="-|>", lw=2.2, mutation_scale=18)
ax.annotate("", xy=(4.15, 2.75), xytext=(1.9, 1.55), arrowprops=dict(color=col["I"], **arrow))
ax.text(2.3, 2.45, r"infection  $\beta$", color=col["I"], ha="center", fontsize=13)
ax.annotate("", xy=(8.1, 1.55), xytext=(5.85, 2.75), arrowprops=dict(color=col["R"], **arrow))
ax.text(7.75, 2.45, r"recovery  $\gamma$", color=col["R"], ha="center", fontsize=13)
ax.annotate("", xy=(2.0, 0.95), xytext=(8.0, 0.95),
            arrowprops=dict(color=col["S"], connectionstyle="arc3,rad=0.25", **arrow))
ax.text(5, 0.12, r"loss of immunity  $\omega$", color=col["S"], ha="center", fontsize=13)
ax.set_title("SIRS model: compartments and rates", fontsize=12)
fig.savefig("img/s2_cartoon_rates.png", dpi=200, bbox_inches="tight")
plt.close(fig)


# ---------------- ODE and SSA (as in the notebook) ----------------
def sirs_rhs(t, y, beta, gamma, omega):
    S, I, R = y
    return [-beta*S*I + omega*R, beta*S*I - gamma*I, gamma*I - omega*R]

def solve_sirs(N, I0=1, t_max=T_MAX, n=2000):
    t_eval = np.linspace(0, t_max, n)
    sol = solve_ivp(sirs_rhs, (0, t_max), [N-I0, I0, 0], args=(beta, gamma, omega),
                    t_eval=t_eval, rtol=1e-9, atol=1e-12)
    return sol.t, sol.y

def endemic_I(N):
    return max(0.0, omega/(omega+gamma) * (N - gamma/beta))

def ssa_sirs(N, I0=1, t_max=T_MAX, rng=None):
    S, I, R = N - I0, I0, 0
    t = 0.0
    times, states = [t], [(S, I, R)]
    while t < t_max and I > 0:
        a = (beta*S*I, gamma*I, omega*R)
        a0 = a[0] + a[1] + a[2]
        t += rng.exponential(1.0/a0)
        u = rng.random() * a0
        if u < a[0]:
            S -= 1; I += 1
        elif u < a[0] + a[1]:
            I -= 1; R += 1
        else:
            R -= 1; S += 1
        times.append(t); states.append((S, I, R))
    return np.array(times), np.array(states)


# ---------------- Three separate, labelled SSA panels ----------------
rng = np.random.default_rng(2026)          # same seed and call order as the notebook
rng_small = np.random.default_rng(10)      # separate stream for N = 10, so the other panels stay unchanged
for N in [10, 20, 50, 100]:
    if N == 10:
        rng_N = rng_small
    else:
        rng_N = rng
    fig, ax = plt.subplots(figsize=(4.6, 3.1))
    t, y = solve_sirs(N)
    for _ in range(8):
        ts, ss = ssa_sirs(N, rng=rng_N)
        ax.step(np.append(ts, T_MAX), np.append(ss[:, 1], ss[-1, 1]), where="post", lw=0.8, alpha=0.7)
    ax.plot(t, y[1], "k", lw=2.5)
    ax.axhline(endemic_I(N), color="k", ls="--", lw=1)
    ax.set_title(f"N = {N} people: 8 SSA runs vs ODE", fontsize=13, fontweight="bold")
    ax.set_xlabel("time (days)", fontsize=12)
    ax.set_ylabel("infected people, $I(t)$", fontsize=12)
    handles = [Line2D([], [], color="k", lw=2.5, label="ODE solution"),
               Line2D([], [], color="k", ls="--", lw=1, label=f"ODE equilibrium $I^*$ = {endemic_I(N):.2f}"),
               Line2D([], [], color="tab:orange", lw=1, label="SSA runs (one colour per run)")]
    ax.legend(handles=handles, fontsize=9.5, loc="upper right", framealpha=0.9)
    ax.tick_params(labelsize=10.5)
    fig.tight_layout()
    fig.savefig(f"img/s5_ssa_labelled_N{N}.png", dpi=200)
    plt.close(fig)
print("done")
