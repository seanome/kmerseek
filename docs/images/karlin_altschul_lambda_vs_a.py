import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C = 2.0

def karlin_altschul_lambda(a, penalty=C):
    """Positive root of a e^x + (1-a) e^(-penalty x) = 1, or 0 when none exists."""
    b = 1.0 - a
    if a <= 0 or a - penalty * b >= 0:
        return 0.0
    f = lambda x: a * math.exp(x) + b * math.exp(-penalty * x)
    lo, hi = 0.0, 1.0
    while f(hi) <= 1.0:
        hi *= 2
    for _ in range(200):
        mid = (lo + hi) / 2
        if f(mid) > 1.0:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2

xs = [i / 1000 for i in range(1, 1000)]
ys = [karlin_altschul_lambda(a) for a in xs]
boundary = C / (1 + C)

fig, ax = plt.subplots(figsize=(7.5, 4.6), dpi=150)
ax.plot(xs, ys, color="#1f77b4", lw=2.2, label="λ at mismatch penalty C = 2")
ax.axvline(boundary, color="#d62728", ls="--", lw=1.6,
           label="a = C / (1 + C) = 2/3: expected score of a random position is 0")
ex = [(0.5, "balanced pair\na = 0.5, λ = %.2f" % karlin_altschul_lambda(0.5), (0.33, 0.55)),
      (0.85, "two hydrophobic sequences\na = 0.85, λ = 0, never significant", (0.0, 0.25))]
for a, text, off in ex:
    lam = karlin_altschul_lambda(a)
    ax.plot(a, lam, "o", color="#111111", ms=7, zorder=5,
            label="example pair" if a == 0.5 else None)
    ax.annotate(text, (a, lam), xytext=(a + off[0], lam + off[1]),
                ha="center", fontsize=9,
                arrowprops=dict(arrowstyle="-", color="#555555", lw=0.8))

ax.set_xlim(0, 1)
ax.set_ylim(0, 2.6)
ax.set_xlabel("a: chance that a random query position and a random target position\n"
              "fall in the same class (from this pair's own class frequencies)")
ax.set_ylabel("λ (per unit of raw score S)")
ax.set_title("λ falls to 0 as the two compositions make matching likely by chance",
             fontsize=11, loc="left", pad=62)
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), frameon=False, fontsize=9, ncol=1)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
out = "docs/images/karlin_altschul_lambda_vs_a"
fig.savefig(out + ".png", bbox_inches="tight")
fig.savefig(out + ".svg", bbox_inches="tight")
print(karlin_altschul_lambda(0.5), math.log((1 + 5 ** 0.5) / 2))
