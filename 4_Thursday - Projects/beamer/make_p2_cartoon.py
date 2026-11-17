"""Two strain cartoon (img/p2fig1.png) for Project3_8slides.tex.
Same drawing as cell 3 of Project3_Part2_TwoStrain_SSA_vs_ODE_backup_full_model.ipynb,
with the death flow (rate mu) drawn out of every compartment.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C = dict(S="#4C6A86", I1="#BF4E2C", R1="#228066", I2="#7C4A9C", R2="#B3802F")
DEATH = "#555555"

fig, ax = plt.subplots(figsize=(11, 6.8)); ax.axis("off"); ax.set_xlim(-1.4, 12.6); ax.set_ylim(-1.2, 11.0)
pos = dict(S=(1.5, 5.0), I1=(5.2, 8.0), R1=(10.2, 8.0), I2=(5.2, 2.0), R2=(10.2, 2.0))
name = dict(S="S", I1="$I_1$", R1="$R_1$", I2="$I_2$", R2="$R_2$")
HW, HH = 0.8, 0.45
for k, (x, y) in pos.items():
    ax.add_patch(plt.Rectangle((x-HW, y-HH), 2*HW, 2*HH, fc="white", ec=C[k], lw=2.6, zorder=5))
    ax.text(x, y, name[k], ha="center", va="center", fontsize=17, fontweight="bold", color=C[k], zorder=6)
def arrow(p, q, color, rad=0.0, ls="-", lw=2.2):
    ax.annotate("", xy=q, xytext=p, zorder=3,
                arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, ls=ls, mutation_scale=16,
                                connectionstyle=f"arc3,rad={rad}", shrinkA=0, shrinkB=0))
# same strain loops (solid)
arrow((2.3, 5.45), (4.4, 7.7), C["I1"]);   ax.text(2.1, 7.0, r"$\beta_1 S I_1/N$", color=C["I1"], fontsize=11)
arrow((6.0, 8.0), (9.4, 8.0), C["R1"]);    ax.text(7.7, 8.25, r"$\gamma_1 I_1$", color=C["R1"], fontsize=11, ha="center")
arrow((2.3, 4.55), (4.4, 2.3), C["I2"]);   ax.text(2.1, 2.7, r"$\beta_2 S I_2/N$", color=C["I2"], fontsize=11)
arrow((6.0, 2.0), (9.4, 2.0), C["R2"]);    ax.text(7.7, 1.45, r"$\gamma_2 I_2$", color=C["R2"], fontsize=11, ha="center")
# waning back to S (solid grey-blue arches over and under)
arrow((10.2, 8.5), (1.5, 5.5), C["S"], rad=0.30);   ax.text(6.0, 10.15, r"$\omega_1 R_1$ (waning)", color=C["S"], fontsize=11, ha="center")
arrow((10.2, 1.5), (1.5, 4.5), C["S"], rad=-0.30);  ax.text(6.0, -0.55, r"$\omega_2 R_2$ (waning)", color=C["S"], fontsize=11, ha="center")
# cross infections (dashed): R2 -> I1 and R1 -> I2, straight so that the direction is obvious
arrow((9.55, 2.45), (5.75, 7.55), C["I1"], ls="--", lw=2.4)
arrow((9.55, 7.55), (5.75, 2.45), C["I2"], ls="--", lw=2.4)
ax.text(9.55, 4.2, r"$R_2\to I_1$:  $\beta_1 R_2 I_1/N$", color=C["I1"], fontsize=11, ha="center", bbox=dict(fc="white", ec="none", pad=1.5), zorder=7)
ax.text(9.55, 5.8, r"$R_1\to I_2$:  $\beta_2 R_1 I_2/N$", color=C["I2"], fontsize=11, ha="center", bbox=dict(fc="white", ec="none", pad=1.5), zorder=7)
# births
arrow((-1.2, 5.0), (0.7, 5.0), "gray");  ax.text(-1.3, 5.4, r"birth $\mu N$", color="gray", fontsize=10)
# deaths: one short grey arrow out of every compartment
arrow((0.9, 4.55), (0.9, 3.5), DEATH, lw=1.8);      ax.text(0.78, 3.75, r"$\mu S$", color=DEATH, fontsize=11, ha="right")
arrow((5.0, 7.55), (5.0, 6.6), DEATH, lw=1.8);      ax.text(4.88, 6.85, r"$\mu I_1$", color=DEATH, fontsize=11, ha="right")
arrow((5.0, 2.45), (5.0, 3.4), DEATH, lw=1.8);      ax.text(4.88, 3.0, r"$\mu I_2$", color=DEATH, fontsize=11, ha="right")
arrow((11.0, 8.0), (12.1, 8.0), DEATH, lw=1.8);     ax.text(11.55, 8.2, r"$\mu R_1$", color=DEATH, fontsize=11, ha="center")
arrow((11.0, 2.0), (12.1, 2.0), DEATH, lw=1.8);     ax.text(11.55, 2.2, r"$\mu R_2$", color=DEATH, fontsize=11, ha="center")
ax.text(-1.3, -0.9, r"dark grey arrows: death at rate $\mu$ from every class", color=DEATH, fontsize=9.5, ha="left")
ax.set_title("Two strain model with waning immunity (solid: same strain, dashed: cross infection)", fontsize=12)
fig.savefig("img/p2fig1.png", dpi=200, bbox_inches="tight")
plt.close(fig)
