#!/usr/bin/env python3
"""
Plot the scaling CSV emitted by scaling.f90.

TODO 5b/5c: this script is deliberately incomplete. Fill in the metric
computation so the plot shows what you actually want to argue.

Usage:
    ./plot_scaling.py scaling_strong.csv scaling_weak.csv -o scaling.png

Input rows look like:
    strong,ranks,4,nlocal,1000000,time,1.234560E-01,sweeps,200
"""

import argparse
import sys


def read_csv(path):
    """Parse the harness's key,value CSV into a list of dicts."""
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(",")
            mode = parts[0]
            kv = dict(zip(parts[1::2], parts[2::2]))
            try:
                rows.append({
                    "mode": mode,
                    "ranks": int(kv["ranks"]),
                    "nlocal": int(kv["nlocal"]),
                    "time": float(kv["time"]),
                    "sweeps": int(kv["sweeps"]),
                })
            except (KeyError, ValueError):
                print(f"skipping malformed line: {line}", file=sys.stderr)
    return sorted(rows, key=lambda r: r["ranks"])


def strong_metrics(rows):
    """
    TODO 5b: return (ranks, speedup, efficiency, karpflatt) lists.

        speedup    S = t[0] / t[i]
        efficiency E = S / (P / P[0])
        Karp-Flatt e = (1/S - 1/P) / (1 - 1/P)     for P > 1

    Use the SAME formulas as parallel_metrics in scaling.f90 -- if the
    Fortran and the Python disagree, you will not know which to trust.
    """
    ranks = [r["ranks"] for r in rows]
    speedup, eff, karp = [], [], []
    # TODO 5b
    return ranks, speedup, eff, karp


def weak_metrics(rows):
    """
    TODO 5c: weak-scaling efficiency is t[0] / t[i] -- ideal is a FLAT line
    at 1.0, not a rising one. Do not reuse the strong-scaling formula here;
    that mistake is what produces the classic misleading scaling plot.
    """
    ranks = [r["ranks"] for r in rows]
    eff = []
    # TODO 5c
    return ranks, eff


def _fmt(seq, i):
    return f"{seq[i]:>12.4f}" if i < len(seq) else f"{'--':>12}"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("csv", nargs="+", help="CSV files from the harness")
    ap.add_argument("-o", "--out", default="scaling.png")
    args = ap.parse_args()

    strong, weak = [], []
    for path in args.csv:
        for row in read_csv(path):
            (weak if row["mode"] == "weak" else strong).append(row)

    if not strong and not weak:
        sys.exit("no usable rows found -- did the harness run?")

    # Text report first: it needs no dependencies, and you should be able to
    # read your own numbers before you look at a picture.
    if strong:
        print("=== strong scaling ===")
        print(f"{'ranks':>8} {'time (s)':>12} {'speedup':>12} "
              f"{'efficiency':>12} {'karp-flatt':>12}")
        _, sp, ef, kf = strong_metrics(strong)
        for i, r in enumerate(strong):
            print(f"{r['ranks']:>8} {r['time']:>12.6f} "
                  f"{_fmt(sp, i)} {_fmt(ef, i)} {_fmt(kf, i)}")
        if not sp:
            print("  (metrics empty -- TODO 5b)")

    if weak:
        print("\n=== weak scaling ===")
        print(f"{'ranks':>8} {'time (s)':>12} {'efficiency':>12}")
        _, ef = weak_metrics(weak)
        for i, r in enumerate(weak):
            print(f"{r['ranks']:>8} {r['time']:>12.6f} {_fmt(ef, i)}")
        if not ef:
            print("  (metrics empty -- TODO 5c)")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("\n(matplotlib not installed -- text report only)")
        print("  pip install matplotlib   to get the figure")
        return

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    if strong:
        ranks, sp, _, _ = strong_metrics(strong)
        if sp:
            axes[0].plot(ranks, sp, "o-", label="measured")
            axes[0].plot(ranks, [r / ranks[0] for r in ranks], "k--",
                         label="ideal")
    axes[0].set_xscale("log", base=2)
    axes[0].set_yscale("log", base=2)
    axes[0].set_xlabel("MPI ranks")
    axes[0].set_ylabel("speedup")
    axes[0].set_title("Strong scaling")
    axes[0].legend()
    axes[0].grid(alpha=0.3)

    if weak:
        ranks, ef = weak_metrics(weak)
        if ef:
            axes[1].plot(ranks, ef, "s-", label="measured")
    axes[1].axhline(1.0, color="k", ls="--", label="ideal")
    axes[1].set_xscale("log", base=2)
    axes[1].set_ylim(0, 1.15)
    axes[1].set_xlabel("MPI ranks")
    axes[1].set_ylabel("parallel efficiency")
    axes[1].set_title("Weak scaling")
    axes[1].legend()
    axes[1].grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(args.out, dpi=140)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
