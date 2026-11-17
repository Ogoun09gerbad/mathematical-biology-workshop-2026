"""Part 2: ODE and Gillespie SSA for the system (2.1)-(2.5) of the sheet, and for the sheet's code (with cross infection)."""
import numpy as np, time
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

mu = 0.02; b1, g1, w1 = 120.0, 52.0, 0.2; b2, g2, w2 = 100.0, 26.0, 0.5
N = 10_000; T = 50
y0_ode = [0.99, 0.005, 0.005, 0.0, 0.0]; y0_cnt = [9900, 50, 50, 0, 0]

def rhs(t, y, cross):
    S, I1, I2, R1, R2 = y
    c = 1.0 if cross else 0.0
    return [mu*(1 - S) - b1*S*I1 - b2*S*I2 + w1*R1 + w2*R2,
            b1*S*I1 + c*b1*R2*I1 - (g1 + mu)*I1,
            b2*S*I2 + c*b2*R1*I2 - (g2 + mu)*I2,
            g1*I1 - c*b2*R1*I2 - (w1 + mu)*R1,
            g2*I2 - c*b1*R2*I1 - (w2 + mu)*R2]

def ode(cross, t_end=T, n=50001):
    return solve_ivp(rhs, (0, t_end), y0_ode, args=(cross,), t_eval=np.linspace(0, t_end, n),
                     method="LSODA", rtol=1e-10, atol=1e-30)

def ssa(cross, t_max=T, dt=0.01, rng=None, stop=False):
    S, I1, I2, R1, R2 = y0_cnt
    grid = np.arange(0.0, t_max + dt/2, dt); out = np.zeros((len(grid), 5), dtype=np.int64)
    gi, t, k = 0, 0.0, 0; ext = [np.nan, np.nan]; U = rng.random((200_000, 2)); log = np.log
    while t < t_max:
        while gi < len(grid) and grid[gi] <= t:
            out[gi] = (S, I1, I2, R1, R2); gi += 1
        if stop and I1 == 0 and I2 == 0: break
        a = (mu*N, mu*S, mu*I1, mu*I2, mu*R1, mu*R2, b1*S*I1/N, b2*S*I2/N,
             b1*R2*I1/N if cross else 0.0, b2*R1*I2/N if cross else 0.0, g1*I1, g2*I2, w1*R1, w2*R2)
        a0 = sum(a)
        if a0 == 0: break
        if k == len(U): U = rng.random((200_000, 2)); k = 0
        u1, u2 = U[k]; k += 1
        t += -log(u1)/a0; x = u2*a0; c = 0.0
        for e in range(14):
            c += a[e]
            if x < c: break
        if   e == 0: S += 1
        elif e == 1: S -= 1
        elif e == 2: I1 -= 1
        elif e == 3: I2 -= 1
        elif e == 4: R1 -= 1
        elif e == 5: R2 -= 1
        elif e == 6: S -= 1; I1 += 1
        elif e == 7: S -= 1; I2 += 1
        elif e == 8: I1 += 1; R2 -= 1
        elif e == 9: I2 += 1; R1 -= 1
        elif e == 10: I1 -= 1; R1 += 1
        elif e == 11: I2 -= 1; R2 += 1
        elif e == 12: S += 1; R1 -= 1
        else: S += 1; R2 -= 1
        if I1 == 0 and np.isnan(ext[0]): ext[0] = t
        if I2 == 0 and np.isnan(ext[1]): ext[1] = t
    while gi < len(grid):
        out[gi] = (S, I1, I2, R1, R2); gi += 1
    return grid, out, ext

plt.rcParams.update({"font.size": 12, "axes.grid": True, "grid.alpha": 0.3, "figure.dpi": 110})
RED, BLUE = "crimson", "royalblue"
labels = {False: "System (2.1)-(2.5) of the sheet (no cross infection)", True: "Sheet's code (with cross infection)"}
res = {}
for cross in (False, True):
    sol = ode(cross)
    g, st, ext = ssa(cross, rng=np.random.default_rng(2026))
    res[cross] = (sol, g, st, ext)
    print(f"\n=== {labels[cross]} ===")
    print("ODE at t=50  [S,I1,I2,R1,R2] =", np.array2string(sol.y[:, -1], precision=5))
    i2 = np.argmax(sol.y[2]); im = i2 + np.argmin(sol.y[2][i2:i2+8000])
    print(f"ODE strain 2: peak {sol.y[2][i2]:.3f} at t={sol.t[i2]:.2f}; minimum after peak {sol.y[2][im]:.2e} (= {sol.y[2][im]*N:.2f} people) at t={sol.t[im]:.2f}")
    print(f"ODE strain 1 at t=1 y: {sol.y[1][1000]:.1e} ({sol.y[1][1000]*N:.1e} people)")
    print("SSA seed 2026: extinction year strain 1, strain 2 =", np.round(ext, 3), "| final counts", st[-1], "| total pop min/max", st.sum(1).min(), st.sum(1).max())
    # 40 SSA runs, 50 years
    t0 = time.time(); rng = np.random.default_rng(7); E = []
    for _ in range(40): E.append(ssa(cross, rng=rng, dt=0.5, stop=True)[2])
    E = np.array(E)
    print(f"40 SSA runs of 50 years ({time.time()-t0:.0f} s): strain 1 lost in {np.mean(~np.isnan(E[:,0])):.0%} (median year {np.nanmedian(E[:,0]) if np.any(~np.isnan(E[:,0])) else float('nan'):.2f}), "
          f"strain 2 lost in {np.mean(~np.isnan(E[:,1])):.0%} (median year {np.nanmedian(E[:,1]) if np.any(~np.isnan(E[:,1])) else float('nan'):.2f}), both lost in {np.mean(~np.isnan(E).any(1) & ~np.isnan(E).any(1)) if False else np.mean(~np.isnan(E[:,0]) & ~np.isnan(E[:,1])):.0%}")
    res[cross] += (E,)

# Figure 1 (asked by the sheet): ODE and SSA side by side, 50 years, for each system
for cross, fname in ((False, "fig1_system_2.1-2.5_ODE_vs_SSA.png"), (True, "fig1b_with_cross_infection_ODE_vs_SSA.png")):
    sol, g, st, ext, E = res[cross]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
    a1.plot(sol.t, sol.y[1], RED, label="Strain 1 ($I_1$)"); a1.plot(sol.t, sol.y[2], BLUE, label="Strain 2 ($I_2$)")
    a1.set(title="Deterministic ODE", xlabel="Years", ylabel="Infected fraction"); a1.legend()
    a2.plot(g, st[:, 1]/N, RED, alpha=.7, label="Strain 1 ($I_1$)"); a2.plot(g, st[:, 2]/N, BLUE, alpha=.7, label="Strain 2 ($I_2$)")
    a2.set(title=f"Stochastic Gillespie SSA (N = {N:,}, one run)", xlabel="Years"); a2.legend()
    fig.suptitle(labels[cross]); plt.tight_layout(); plt.savefig(fname, dpi=130); plt.close()

# Figure 2: first 5 years, in people (log scale), ODE and 10 SSA runs
for cross, fname in ((False, "fig2_system_2.1-2.5_first_5_years.png"), (True, "fig2b_with_cross_infection_first_5_years.png")):
    sol = res[cross][0]; rng = np.random.default_rng(11)
    runs = [ssa(cross, t_max=5, dt=0.005, rng=rng) for _ in range(10)]
    fig, axs = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
    for ax, k, col, nm in ((axs[0], 1, RED, "strain 1"), (axs[1], 2, BLUE, "strain 2")):
        for j, (gg, oo, ee) in enumerate(runs):
            keep = oo[:, k] > 0
            ax.step(gg[keep], oo[keep, k], where="post", color=col, alpha=.4, lw=.9, label="SSA runs (people)" if j == 0 else None)
            if not np.isnan(ee[k-1]) and ee[k-1] < 5: ax.plot(ee[k-1], 1, "kx", ms=8, mew=2, label="SSA hits 0" if j == 0 else None)
        m = sol.t <= 5
        ax.semilogy(sol.t[m], np.maximum(sol.y[k][m]*N, 1e-12), "k", lw=2.4, label="ODE $\times N$")
        ax.axhline(1, color="grey", ls=":"); ax.set(title=f"{nm}: infected people (log scale)", xlabel="Years", ylim=(1e-8, 1e4)); ax.legend(loc="lower right")
    fig.suptitle(labels[cross] + ": first 5 years"); plt.tight_layout(); plt.savefig(fname, dpi=130); plt.close()
print("\nfigures saved")
