"""Figures for Project3_presentation.tex, drawn large enough to be read on a slide.

The models, parameters and random seeds are the same as in the notebooks
Project3_Part1_SIRS_SSA_vs_ODE.ipynb and Project3_Part2_TwoStrain_SSA_vs_ODE_backup_full_model.ipynb,
so the numbers on the slides are the numbers of the notebooks.
Run from the beamer folder:  python make_slide_figures.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

plt.rcParams.update({"font.size": 15, "axes.titlesize": 16, "axes.labelsize": 15, "legend.fontsize": 12.5,
                     "xtick.labelsize": 13, "ytick.labelsize": 13, "axes.grid": True, "grid.alpha": 0.3,
                     "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 200})
RED, BLUE, GREEN, SLATE, GOLD, PURPLE = "#BF4E2C", "#3B6FB6", "#228066", "#4C6A86", "#B3802F", "#7C4A9C"
OUT = "img/"

# ======================================================================
# Part 1: SIRS (same code and seeds as the Part 1 notebook)
# ======================================================================
beta, gamma, omega, T_MAX = 0.01, 0.15, 0.025, 365.0
N_crit = gamma / beta

def endemic_I(N):
    return max(0.0, omega/(omega+gamma) * (N - gamma/beta))

def solve_sirs(N, I0=1, t_max=T_MAX, n=2000):
    rhs = lambda t, y: [-beta*y[0]*y[1] + omega*y[2], beta*y[0]*y[1] - gamma*y[1], gamma*y[1] - omega*y[2]]
    sol = solve_ivp(rhs, (0, t_max), [N-I0, I0, 0], t_eval=np.linspace(0, t_max, n), rtol=1e-9, atol=1e-12)
    return sol.t, sol.y

def ssa_sirs(N, I0=1, t_max=T_MAX, rng=None):
    S, I, R = N - I0, I0, 0
    t = 0.0; times, states = [t], [(S, I, R)]
    while t < t_max and I > 0:
        a = (beta*S*I, gamma*I, omega*R); a0 = a[0] + a[1] + a[2]
        t += rng.exponential(1.0/a0); u = rng.random() * a0
        if u < a[0]: S -= 1; I += 1
        elif u < a[0] + a[1]: I -= 1; R += 1
        else: R -= 1; S += 1
        times.append(t); states.append((S, I, R))
    return np.array(times), np.array(states)

# --- Slide: SIRS cartoon ---
fig, ax = plt.subplots(figsize=(9, 3.6)); ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(-0.3, 4.2)
pos1 = {"S": (1.2, 1.2), "I": (5, 3.1), "R": (8.8, 1.2)}; col1 = {"S": SLATE, "I": RED, "R": GREEN}
for k_, (x, y_) in pos1.items():
    ax.add_patch(plt.Rectangle((x-0.8, y_-0.45), 1.6, 0.9, fc="white", ec=col1[k_], lw=3))
    ax.text(x, y_, k_, ha="center", va="center", fontsize=24, fontweight="bold", color=col1[k_])
st_ = dict(arrowstyle="-|>", lw=2.8, mutation_scale=22)
ax.annotate("", xy=(4.15, 2.85), xytext=(1.9, 1.6), arrowprops=dict(color=RED, **st_))
ax.text(1.9, 2.7, r"infection  $\beta SI$", color=RED, ha="center", fontsize=17)
ax.annotate("", xy=(8.1, 1.6), xytext=(5.85, 2.85), arrowprops=dict(color=GREEN, **st_))
ax.text(8.1, 2.7, r"recovery  $\gamma I$", color=GREEN, ha="center", fontsize=17)
ax.annotate("", xy=(2.0, 0.95), xytext=(8.0, 0.95), arrowprops=dict(color=SLATE, connectionstyle="arc3,rad=0.25", **st_))
ax.text(5, -0.15, r"loss of immunity  $\omega R$", color=SLATE, ha="center", fontsize=17)
plt.tight_layout(); plt.savefig(OUT + "slide_p1_cartoon.png", bbox_inches="tight"); plt.close()

# --- Slide: above the threshold, small populations still lose the infection ---
N_list = [10, 12, 15, 18, 20, 25, 30, 40, 50, 70, 100]
p_alive = []
for N in N_list:                                   # 500 runs per N, seed = N (as in the notebook)
    rng = np.random.default_rng(N)
    p_alive.append(np.mean([ssa_sirs(N, rng=rng)[1][-1, 1] > 0 for _ in range(500)]))
print("Part 1, P(I>0 at day 365):", dict(zip(N_list, np.round(p_alive, 3))))

fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw=dict(width_ratios=[1.15, 1]))
ax = axs[0]; N = 50
rng = np.random.default_rng(2026)
for N_skip in (20,):                                # the notebook draws N = 20 first with the same generator
    for _ in range(8): ssa_sirs(N_skip, rng=rng)
for j in range(8):
    ts, ss = ssa_sirs(N, rng=rng)
    ax.step(np.append(ts, T_MAX), np.append(ss[:, 1], ss[-1, 1]), where="post", lw=1.2, alpha=0.75)
t, y = solve_sirs(N)
ax.plot(t, y[1], "k", lw=3, label="ODE")
ax.axhline(endemic_I(N), color="k", ls="--", lw=1.2, label=f"ODE endemic level $I^*$ = {endemic_I(N):.0f}")
ax.set(xlabel="days", ylabel="infected people $I$", title=f"N = {N}: 8 stochastic runs (colours) vs ODE")
ax.legend(loc="upper right")
ax = axs[1]
ax.plot([10, N_crit, N_crit, 100], [0, 0, 1, 1], "k--", lw=2.5, label="ODE: endemic for $N>15$")
ax.plot(N_list, p_alive, "o-", color=RED, lw=2.5, ms=8, label="stochastic model")
ax.text(N_crit + 2, 0.86, "ODE threshold $N_c=\\gamma/\\beta=15$", fontsize=12.5, color="grey")
ax.set(xlabel="population size $N$", ylabel="fraction of runs still infected\nat day 365", title="Persistence over 500 runs", ylim=(-0.05, 1.08))
ax.legend(loc="center left", bbox_to_anchor=(0.12, 0.5))
plt.tight_layout(); plt.savefig(OUT + "slide_p1_persistence.png"); plt.close()

# ======================================================================
# Part 2: two strains with cross infection (same code and seeds as the Part 2 notebook)
# ======================================================================
mu = 0.02
beta1, gamma1, omega1 = 120.0, 52.0, 0.2
beta2, gamma2, omega2 = 100.0, 26.0, 0.5
R0_1, R0_2 = beta1/(gamma1+mu), beta2/(gamma2+mu)

def two_strain(t, y, c=1.0):
    """c = 1: the sheet's code (cross infection); c = 0: equations (2.1)-(2.5) as printed on the sheet."""
    S, I1, I2, R1, R2 = y
    return [mu*(1 - S) - beta1*S*I1 - beta2*S*I2 + omega1*R1 + omega2*R2,
            beta1*(S + c*R2)*I1 - (gamma1 + mu)*I1,
            beta2*(S + c*R1)*I2 - (gamma2 + mu)*I2,
            gamma1*I1 - c*beta2*R1*I2 - (omega1 + mu)*R1,
            gamma2*I2 - c*beta1*R2*I1 - (omega2 + mu)*R2]

def ssa_two_strain(N, y0, t_max, dt=0.01, rng=None, cross=1.0):
    """Gillespie direct method, 14 events, state stored on a time grid; also returns extinction times.
    cross = 0 switches the two cross infection events off (equations as printed on the sheet)."""
    S, I1, I2, R1, R2 = [int(v) for v in y0]
    grid = np.arange(0.0, t_max + dt/2, dt); out = np.zeros((len(grid), 5), dtype=np.int64)
    gi, t, k = 0, 0.0, 0; ext = [np.nan, np.nan]; U = rng.random((200_000, 2)); log = np.log
    while t < t_max:
        while gi < len(grid) and grid[gi] <= t:
            out[gi] = (S, I1, I2, R1, R2); gi += 1
        a = (mu*N, mu*S, mu*I1, mu*I2, mu*R1, mu*R2, beta1*S*I1/N, beta2*S*I2/N, cross*beta1*R2*I1/N, cross*beta2*R1*I2/N,
             gamma1*I1, gamma2*I2, omega1*R1, omega2*R2)
        a0 = sum(a)
        if k == len(U): U = rng.random((200_000, 2)); k = 0
        u1, u2 = U[k]; k += 1
        t += -log(u1) / a0; x = u2 * a0; c = 0.0
        for e in range(14):
            c += a[e]
            if x < c: break
        if   e == 0:  S += 1
        elif e == 1:  S -= 1
        elif e == 2:  I1 -= 1
        elif e == 3:  I2 -= 1
        elif e == 4:  R1 -= 1
        elif e == 5:  R2 -= 1
        elif e == 6:  S -= 1; I1 += 1
        elif e == 7:  S -= 1; I2 += 1
        elif e == 8:  I1 += 1; R2 -= 1
        elif e == 9:  I2 += 1; R1 -= 1
        elif e == 10: I1 -= 1; R1 += 1
        elif e == 11: I2 -= 1; R2 += 1
        elif e == 12: S += 1; R1 -= 1
        else:         S += 1; R2 -= 1
        if I1 == 0 and np.isnan(ext[0]): ext[0] = t
        if I2 == 0 and np.isnan(ext[1]): ext[1] = t
    while gi < len(grid):
        out[gi] = (S, I1, I2, R1, R2); gi += 1
    return grid, out, ext

N = 10_000; y0_ode = [0.99, 0.005, 0.005, 0.0, 0.0]; y0_stoch = [9900, 50, 50, 0, 0]
sol = solve_ivp(two_strain, (0, 50), y0_ode, t_eval=np.linspace(0, 50, 5001), method="LSODA", rtol=1e-10, atol=1e-16)

# --- Slide: two strain cartoon ---
fig, ax = plt.subplots(figsize=(11, 6.2)); ax.axis("off"); ax.set_xlim(-1.6, 12.2); ax.set_ylim(-0.6, 10.4)
pos = dict(S=(1.5, 5.0), I1=(5.4, 8.0), R1=(10.4, 8.0), I2=(5.4, 2.0), R2=(10.4, 2.0))
col = dict(S=SLATE, I1=RED, R1=GREEN, I2=PURPLE, R2=GOLD)
name = dict(S="S", I1="$I_1$", R1="$R_1$", I2="$I_2$", R2="$R_2$")
for k_, (x, y_) in pos.items():
    ax.add_patch(plt.Rectangle((x-0.85, y_-0.5), 1.7, 1.0, fc="white", ec=col[k_], lw=3, zorder=5))
    ax.text(x, y_, name[k_], ha="center", va="center", fontsize=22, fontweight="bold", color=col[k_], zorder=6)
def arrow(p, q, color, rad=0.0, ls="-", lw=2.6):
    ax.annotate("", xy=q, xytext=p, zorder=3, arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, ls=ls, mutation_scale=20,
                connectionstyle=f"arc3,rad={rad}", shrinkA=0, shrinkB=0))
arrow((2.35, 5.5), (4.55, 7.7), RED);    ax.text(1.4, 6.9, "infection", color=RED, fontsize=15)
arrow((6.25, 8.0), (9.55, 8.0), GREEN);  ax.text(7.9, 8.3, "recovery", color=GREEN, fontsize=15, ha="center")
arrow((2.35, 4.5), (4.55, 2.3), PURPLE); ax.text(1.4, 2.7, "infection", color=PURPLE, fontsize=15)
arrow((6.25, 2.0), (9.55, 2.0), GOLD);   ax.text(7.9, 1.3, "recovery", color=GOLD, fontsize=15, ha="center")
arrow((10.4, 8.5), (1.5, 5.5), SLATE, rad=0.30);  ax.text(6.0, 10.0, "immunity wanes", color=SLATE, fontsize=15, ha="center")
arrow((10.4, 1.5), (1.5, 4.5), SLATE, rad=-0.30); ax.text(6.0, -0.35, "immunity wanes", color=SLATE, fontsize=15, ha="center")
arrow((9.6, 2.5), (6.2, 7.5), RED, ls="--", lw=2.8)
arrow((9.6, 7.5), (6.2, 2.5), PURPLE, ls="--", lw=2.8)
ax.text(9.0, 5.0, "cross\ninfection", ha="center", va="center", fontsize=16, bbox=dict(fc="white", ec="none"), zorder=7)
arrow((-1.4, 5.0), (0.65, 5.0), "grey"); ax.text(-1.5, 5.35, "births", color="grey", fontsize=16)
ax.text(-1.5, 9.9, "+ deaths from every class", color="grey", fontsize=15)
plt.tight_layout(); plt.savefig(OUT + "slide_p2_cartoon.png", bbox_inches="tight"); plt.close()

# --- Slide: with 10,000 people both models agree (coexistence) ---
_, st, _ = ssa_two_strain(N, y0_stoch, 50, rng=np.random.default_rng(2026))
g50 = np.arange(0.0, 50 + 0.005, 0.01)
fig, ax = plt.subplots(figsize=(12, 4.4))
ax.plot(g50, st[:, 1]/N, color=RED, alpha=0.45, lw=0.8, label="strain 1, stochastic")
ax.plot(g50, st[:, 2]/N, color=BLUE, alpha=0.45, lw=0.8, label="strain 2, stochastic")
ax.plot(sol.t, sol.y[1], color=RED, lw=3, label=f"strain 1, ODE  ($I_1^*$ = {sol.y[1,-1]:.3f})")
ax.plot(sol.t, sol.y[2], color=BLUE, lw=3, label=f"strain 2, ODE  ($I_2^*$ = {sol.y[2,-1]:.3f})")
ax.set(xlabel="years", ylabel="fraction infected", ylim=(0, 0.36), xlim=(0, 50))
ax.legend(ncol=2, loc="upper right")
plt.tight_layout(); plt.savefig(OUT + "slide_p2_coexistence.png"); plt.close()
print("Part 2, N=10,000: SSA mean over years 10-50:", (st[g50 >= 10, 1]/N).mean().round(3), (st[g50 >= 10, 2]/N).mean().round(3),
      "| ODE:", sol.y[1, -1].round(3), sol.y[2, -1].round(3))

# --- Slide: zero is final, and few carriers means a high risk ---
Nv = 300
g_v, o_v, e_v = ssa_two_strain(Nv, [Nv-4, 2, 2, 0, 0], 50, rng=np.random.default_rng(3))
print(f"Part 2, N={Nv} run: strains lost at years {e_v[0]:.2f} and {e_v[1]:.2f}")

def lost_before_established(k, b, g, rng, established=100):
    I, S = k, N - k
    while 0 < I < established:
        a_inf, a_rem = b*S*I/N, (g + mu)*I
        if rng.random() < a_inf/(a_inf + a_rem): I += 1; S -= 1
        else: I -= 1
    return I == 0
rng = np.random.default_rng(0); ks = np.arange(1, 9)
p_sim = {s: np.array([np.mean([lost_before_established(k, b, g, rng) for _ in range(2000)]) for k in ks])
         for s, b, g in ((1, beta1, gamma1), (2, beta2, gamma2))}
print("Part 2, P(extinction | k):", {s: p_sim[s].round(3).tolist() for s in p_sim})

fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw=dict(width_ratios=[1.2, 1]))
ax = axs[0]
ax.step(g_v, o_v[:, 1], color=RED, where="post", lw=1.2, label="strain 1, stochastic")
ax.step(g_v, o_v[:, 2], color=BLUE, where="post", lw=1.2, label="strain 2, stochastic")
ax.plot(sol.t, sol.y[1]*Nv, color=RED, ls="--", lw=2.5, label="strain 1, ODE")
ax.plot(sol.t, sol.y[2]*Nv, color=BLUE, ls="--", lw=2.5, label="strain 2, ODE")
ax.annotate("both reach 0 here\nand never come back", xy=(e_v[1], 0), xytext=(e_v[1] + 4, 80),
            arrowprops=dict(arrowstyle="->", lw=1.5), fontsize=13.5)
ax.set(xlabel="years", ylabel="infected people", title=f"One run with N = {Nv} people", ylim=(-3, 200), xlim=(0, 30))
ax.legend(loc="upper right", fontsize=11.5)
ax = axs[1]
kk = np.linspace(1, 8, 100)
for s, c, R0 in ((1, RED, R0_1), (2, BLUE, R0_2)):
    ax.semilogy(kk, (1/R0)**kk, color=c, lw=2.5, label=f"theory $(1/\\mathcal{{R}}_0)^k$, strain {s}")
    ok = p_sim[s] >= 5/2000                         # at least 5 extinctions out of 2000 runs, otherwise only noise
    ax.semilogy(ks[ok], p_sim[s][ok], "o", color=c, ms=9, mfc="white", mew=2.2, label=f"simulation, strain {s}")
ax.set(xlabel="number $k$ of infected people", ylabel="probability the strain dies out", title="Few carriers = high risk", ylim=(1e-3, 1))
ax.legend(fontsize=11.5, loc="upper right")
plt.tight_layout(); plt.savefig(OUT + "slide_p2_zero_is_final.png"); plt.close()

# --- Slide: smaller populations lose strains (numbers of notebook section D4, 40 runs of 50 years per N) ---
N_vals = np.array([100, 200, 400, 700, 1000])
p_lost = np.array([1.00, 1.00, 0.60, 0.08, 0.02])         # identical for strain 1 and strain 2 in every row
fig, ax = plt.subplots(figsize=(12, 4.4))
ax.plot(N_vals, p_lost, "o-", color=PURPLE, lw=3, ms=10, label="stochastic model (40 runs per N)")
ax.axhline(0, color="k", lw=3, label="ODE: both strains persist, for every N")
for n_, p_ in zip(N_vals, p_lost):
    ax.annotate(f"{p_:.0%}", (n_, p_), textcoords="offset points", xytext=(8, 8), fontsize=14)
ax.set(xlabel="population size $N$", ylabel="fraction of runs losing\nthe strains within 50 years", ylim=(-0.08, 1.12), xlim=(50, 1050))
ax.legend(loc="center right")
plt.tight_layout(); plt.savefig(OUT + "slide_p2_population_size.png"); plt.close()

# --- Slide: equations as printed on the sheet (no cross infection), N = 10,000: the ODE keeps strain 2, the SSA loses it ---
sol_p = solve_ivp(two_strain, (0, 50), y0_ode, args=(0.0,), t_eval=np.linspace(0, 50, 50001), method="LSODA", rtol=1e-10, atol=1e-30)
i_pk = np.argmax(sol_p.y[2]); i_lo = i_pk + np.argmin(sol_p.y[2][i_pk:i_pk + 8000])
print(f"Part 2, printed equations, ODE: strain 2 lowest {sol_p.y[2][i_lo]*N:.2f} people at year {sol_p.t[i_lo]:.2f}, "
      f"{sol_p.y[2][5000]*N:.0f} people at year 5, {sol_p.y[2][-1]*N:.0f} at year 50; strain 1 at year 1: {sol_p.y[1][1000]*N:.1e} people")
T_P, n_runs = 5.0, 200
rng = np.random.default_rng(2026); runs_p = [ssa_two_strain(N, y0_stoch, T_P, rng=rng, cross=0.0) for _ in range(n_runs)]
g_p = runs_p[0][0]
alive2 = np.mean([o[:, 2] > 0 for _, o, _ in runs_p], axis=0)
alive1 = np.mean([o[:, 1] > 0 for _, o, _ in runs_p], axis=0)
print(f"Part 2, printed equations, {n_runs} SSA runs: strain 2 still present at year 1: {alive2[100]:.1%}, at year 5: {alive2[-1]:.1%}; "
      f"strain 1 at year 1: {alive1[100]:.1%}")

fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw=dict(width_ratios=[1.2, 1]))
ax = axs[0]; m = sol_p.t <= 3
for j, (g_, o_, e_) in enumerate(runs_p[:8]):
    keep = (o_[:, 2] > 0) & (g_ <= 3)
    ax.step(g_[keep], o_[keep, 2], where="post", color=BLUE, alpha=0.45, lw=1.0, label="strain 2, 8 stochastic runs" if j == 0 else None)
    if e_[1] < 3: ax.plot(e_[1], 1, "x", color="k", ms=9, mew=2.2, label="stochastic run hits 0" if j == 0 else None)
ax.semilogy(sol_p.t[m], sol_p.y[2][m]*N, color=BLUE, lw=3, label="strain 2, ODE")
ax.semilogy(sol_p.t[m], np.maximum(sol_p.y[1][m]*N, 1e-12), color=RED, lw=3, label="strain 1, ODE (outcompeted)")
ax.axhline(1, color="grey", ls=":", lw=2); ax.text(2.05, 1.4, "one person", color="grey", fontsize=13)
ax.annotate(f"ODE: {sol_p.y[2][i_lo]*N:.2f} of a person,\nthen back up", xy=(sol_p.t[i_lo], sol_p.y[2][i_lo]*N),
            xytext=(1.25, 0.02), arrowprops=dict(arrowstyle="->", lw=1.5), fontsize=13.5)
ax.set(xlabel="years", ylabel="infected people (log scale)", title=f"N = {N:,}: first 3 years", ylim=(1e-3, 1e6), xlim=(0, 3))
ax.legend(loc="upper right", fontsize=11)
ax = axs[1]
ax.plot([0, T_P], [1, 1], color=BLUE, ls="--", lw=3, label=f"ODE: strain 2 persists ({sol_p.y[2][-1]*N:.0f} people)")
ax.plot(g_p, alive2, color=BLUE, lw=3, label=f"stochastic model ({n_runs} runs)")
ax.annotate(f"{alive2[-1]:.0%}", (T_P, alive2[-1]), textcoords="offset points", xytext=(-38, 10), fontsize=15, color=BLUE, fontweight="bold")
ax.set(xlabel="years", ylabel="fraction of runs where\nstrain 2 is still present", title="Does strain 2 survive?", ylim=(-0.05, 1.12), xlim=(0, T_P))
ax.legend(loc="center right", fontsize=11.5)
plt.tight_layout(); plt.savefig(OUT + "slide_p2_printed_equations.png"); plt.close()
print("done")
