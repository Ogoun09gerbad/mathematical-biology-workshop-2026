"""Build the README banner (assets/banner.png): title on the left, two-strain compartment
diagram on the right, faint Gillespie SSA trajectories vs ODE underneath.

Usage (from the repository root):
    python assets/make_banner.py "4_Thursday - Projects/beamer/img/slide_p2_cartoon.png" assets/banner.png
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.image import imread

cartoon_path, out_path = sys.argv[1], sys.argv[2]

C_S, C_I1, C_I2, INK = "#4C6A86", "#BF4E2C", "#7C4A9C", "#1F2A37"

# --- a quick SIRS Gillespie simulation (N=200) and its ODE, for the background ---
N, beta, gamma, omega, T = 200, 0.5, 0.1, 0.02, 400.0

def ssa(seed):
    r = np.random.default_rng(seed)
    S, I, R, t = N - 5, 5, 0, 0.0
    ts, Is = [0.0], [I]
    while t < T and I > 0:
        a = np.array([beta * S * I / N, gamma * I, omega * R])
        a0 = a.sum()
        t += r.exponential(1 / a0)
        k = r.choice(3, p=a / a0)
        if k == 0: S, I = S - 1, I + 1
        elif k == 1: I, R = I - 1, R + 1
        else: R, S = R - 1, S + 1
        ts.append(t); Is.append(I)
    ts.append(T); Is.append(Is[-1])
    return np.array(ts), np.array(Is)

tt = np.linspace(0, T, 4000); dt = tt[1] - tt[0]
s, i, r = (N - 5) / N, 5 / N, 0.0
ode = []
for _ in tt:
    ode.append(i * N)
    ds = -beta * s * i + omega * r; di = beta * s * i - gamma * i; dr = gamma * i - omega * r
    s, i, r = s + dt * ds, i + dt * di, r + dt * dr
ode = np.array(ode)

# --- figure ---
W, H = 16, 5
fig = plt.figure(figsize=(W, H), dpi=200)
fig.patch.set_facecolor("white")

# SSA trajectories vs ODE as a decorative band under the title, fading out on the right
band = fig.add_axes([0.045, 0.04, 0.43, 0.25]); band.axis("off")
band.set_xlim(0, T); band.set_ylim(-3, max(ode) * 1.08)
for seed, col in zip([1, 2, 3, 4, 5, 6], [C_I1, C_I2, C_S, C_I1, C_I2, C_S]):
    ts, Is = ssa(seed)
    band.step(ts, Is, where="post", color=col, lw=1.0, alpha=0.30)
band.plot(tt, ode, color=INK, lw=2.4, alpha=0.55, ls="--")
fade = np.ones((2, 256, 4)); fade[..., 3] = np.clip((np.linspace(0, 1, 256) - 0.55) / 0.45, 0, 1)
band.imshow(fade, extent=[0, T, -3, max(ode) * 1.08], aspect="auto", zorder=10)

# title block
fig.text(0.045, 0.79, "Mathematical Biology", fontsize=37, fontweight="bold", color=INK, va="center")
fig.text(0.045, 0.65, "Workshop 2026", fontsize=37, fontweight="bold", color=C_I1, va="center")
fig.text(0.047, 0.49, "Stochastic (Gillespie SSA) vs deterministic (ODE) models",
         fontsize=15, color="#3B4A5A", va="center")
fig.text(0.047, 0.40, "Epidemiology  ·  Immunology  ·  Python notebooks  ·  Projects",
         fontsize=12.5, color="#6B7785", va="center")

# compartment diagram on the right
img = imread(cartoon_path)
h, w = img.shape[:2]
card_h = 0.84
card_w = card_h * (w / h) * (H / W)
x0 = 0.975 - card_w
card = fig.add_axes([x0, (1 - card_h) / 2, card_w, card_h])
card.imshow(img, interpolation="lanczos")
card.set_xticks([]); card.set_yticks([])
for sp in card.spines.values():
    sp.set_visible(False)

fig.savefig(out_path, dpi=200, facecolor="white")
print(out_path)
