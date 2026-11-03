"""How many chance runs reach E <= c between random sequences, under four ways of counting
Pr(same class). A calibrated rule gives c regions per pair at E <= c, so the printed ratio
observed / (c x pairs) should be near 1; below 1 the E-values are too large (too few chance
regions pass), above 1 too small.

    python scripts/region_pr_same_null.py --classes 2 --pairs 4000 --min-run 10

Each pair is a query and a target of `--length` positions, each drawn from its own class
composition (Dirichlet with all parameters 5, so 2 classes span about 0.2-0.8). Every
maximal run of agreeing positions on every diagonal, at least `--min-run` long, is scored
with kmerseek's run E-value, pr_same^L ((1 - pr_same) m n + (m - 1) + n) for one target, at:

  protein        sum_c p_c q_c of the two whole sequences
  region         sum_c p_c q_c of the run's own positions (one class on both sides, so
                 sum_c f_c^2)
  region a=...   1 / (sum_c sqrt(f_c))^2, f the run's class shares with a pseudo-counts per
                 class spread by p_c q_c / sum p q: kmerseek's `region_pr_same`

Results quoted in kmerseek's doc comments (2026-10-20, numpy 1.26, seed = number of classes),
ratios at E <= 0.01, 0.1, 1:
  2 classes, 4000 pairs, runs >= 10: protein 0.59-0.78, region 0.19-0.25,
      a=0/class 0.43-0.57, a=2/class 0.43-0.55
  6 classes, 6000 pairs, runs >= 4: protein 0.36-0.65, region 0.01-0.03,
      a=0/class 0.02-0.12, a=2/class 0.22-0.23
"""

import argparse

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--classes", type=int, default=2)
    parser.add_argument("--pairs", type=int, default=4000)
    parser.add_argument("--length", type=int, default=300)
    parser.add_argument("--min-run", type=int, default=10)
    parser.add_argument("--per-class", type=float, nargs="+", default=[0, 1, 2, 5])
    args = parser.parse_args()

    C, m = args.classes, args.length
    rng = np.random.default_rng(C)
    rules = ["protein", "region"] + [f"region a={a:g}/class" for a in args.per_class]
    cuts = [0.01, 0.1, 1.0]
    observed = {rule: np.zeros(len(cuts)) for rule in rules}
    for _ in range(args.pairs):
        p_q, p_t = rng.dirichlet(np.full(C, 5.0), 2)
        q = rng.choice(C, m, p=p_q)
        t = rng.choice(C, m, p=p_t)
        f_q = np.bincount(q, minlength=C) / m
        f_t = np.bincount(t, minlength=C) / m
        pr_protein = float(f_q @ f_t)
        shares = f_q * f_t / pr_protein
        n_classes = int((shares > 0).sum())
        for d in range(-(m - 1), m):
            qs, ts = (0, d) if d >= 0 else (-d, 0)
            span = min(m - qs, m - ts)
            agree = q[qs : qs + span] == t[ts : ts + span]
            edges = np.diff(np.concatenate([[0], agree.astype(int), [0]]))
            for s, e in zip(np.where(edges == 1)[0], np.where(edges == -1)[0]):
                run = e - s
                if run < args.min_run:
                    continue
                counts = np.bincount(q[qs + s : qs + e], minlength=C)
                f = counts / run
                prs = [pr_protein, float(f @ f)]
                for a in args.per_class:
                    alpha = a * n_classes
                    prs.append(1.0 / np.sqrt((counts + alpha * shares) / (run + alpha)).sum() ** 2)
                for rule, p in zip(rules, prs):
                    e_value = ((1 - p) * m * m + (m - 1) + m) * p**run
                    observed[rule] += np.array([e_value <= c for c in cuts])
    print(f"{C} classes, {args.pairs} pairs of {m} positions, runs >= {args.min_run}")
    for rule in rules:
        ratios = " ".join(f"E<={c:g}: {o / (c * args.pairs):.2f}" for c, o in zip(cuts, observed[rule]))
        print(f"  {rule:20s} {ratios}")


if __name__ == "__main__":
    main()
