import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C = 2.0

def karlin_altschul_lambda(match_prob, penalty=C):
    """Positive root of u e^x + (1-u) e^(-penalty x) = 1, u = match_prob, or 0 when none exists."""
    mismatch_prob = 1.0 - match_prob
    if match_prob <= 0 or match_prob - penalty * mismatch_prob >= 0:
        return 0.0
    f = lambda x: match_prob * math.exp(x) + mismatch_prob * math.exp(-penalty * x)
    lower, upper = 0.0, 1.0
    while f(upper) <= 1.0:
        upper *= 2
    for _ in range(200):
        mid = (lower + upper) / 2
        if f(mid) > 1.0:
            upper = mid
        else:
            lower = mid
    return (lower + upper) / 2

xs = [i / 1000 for i in range(1, 1000)]
ys = [karlin_altschul_lambda(u) for u in xs]
boundary = C / (1 + C)

fig, ax = plt.subplots(figsize=(7.5, 4.6), dpi=150)
ax.plot(xs, ys, color="#1f77b4", lw=2.2, label="λ at mismatch penalty C = 2")
ax.axvline(boundary, color="#d62728", ls="--", lw=1.6,
           label="u = C / (1 + C) = 2/3: expected score of a random position is 0")
ex = [(0.5, "balanced pair\nu = 0.5, λ = %.2f" % karlin_altschul_lambda(0.5), (0.33, 0.55)),
      (0.85, "two hydrophobic sequences\nu = 0.85, λ = 0, never significant", (0.0, 0.25))]
for u, text, off in ex:
    lam = karlin_altschul_lambda(u)
    ax.plot(u, lam, "o", color="#111111", ms=7, zorder=5,
            label="example pair" if u == 0.5 else None)
    ax.annotate(text, (u, lam), xytext=(u + off[0], lam + off[1]),
                ha="center", fontsize=9,
                arrowprops=dict(arrowstyle="-", color="#555555", lw=0.8))

ax.set_xlim(0, 1)
ax.set_ylim(0, 2.6)
ax.set_xlabel("u: chance that a random query position and a random target position\n"
              "fall in the same class (from this pair's own class frequencies)")
ax.set_ylabel("λ (per unit of raw score S)")
ax.set_title("λ falls to 0 as the two compositions make matching likely by chance",
             fontsize=11, loc="left", pad=62)
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), frameon=False, fontsize=9, ncol=1)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
out = "docs/images/karlin_altschul_lambda_vs_u"
fig.savefig(out + ".png", bbox_inches="tight")
fig.savefig(out + ".svg", bbox_inches="tight")
print(karlin_altschul_lambda(0.5), math.log((1 + 5 ** 0.5) / 2))
