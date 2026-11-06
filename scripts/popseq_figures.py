#!/usr/bin/env python3
"""Figures for the fungal population-sequencing catalog (issue #10).

Reads only catalog/summary/*.csv.gz and catalog/kingdom/fungi_species_summary.csv.gz.
Writes, for each figure, three files into figures/:

    <name>.light.png   for a light page
    <name>.dark.png    for a dark page
    <name>.pdf         vector, light, for print and slides

Light and dark are two renders of the same figure with the palette stepped for
each surface, not an automatic inversion. Colors are the validated default
categorical palette; the six-slot subset used here was checked with the dataviz
validator on 2026-11-06 and passes every adjacent-pair gate in both modes.
Three light slots sit below 3:1 contrast on the light surface, so every series
carries a direct label as well as a legend. Identity is never colour alone.

Every figure states its n and the date of the ENA snapshot it came from.

    uv run --extra viz python scripts/popseq_figures.py
"""
from __future__ import annotations

import csv
import gzip
import os
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import FuncFormatter, LogLocator, MultipleLocator  # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SUMMARY = os.path.join(HERE, "catalog", "summary")
KINGDOM = os.path.join(HERE, "catalog", "kingdom")
OUT = os.path.join(HERE, "figures")

# Validated default categorical palette, slots 1-6, stepped per surface.
THEMES = {
    "light": dict(
        surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e", ink3="#8a8884",
        grid="#e4e3df",
        series=["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"],
        seq="#2a78d6", seq_light="#9ec5f4",
    ),
    "dark": dict(
        surface="#1a1a19", ink="#f0efec", ink2="#c3c2b7", ink3="#8e8d85",
        grid="#33332f",
        series=["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#008300"],
        seq="#3987e5", seq_light="#184f95",
    ),
}

SNAPSHOT = "ENA snapshot 2026-11-06"


def read(path):
    with gzip.open(path, "rt", newline="") as fh:
        first = fh.readline()
        if not first.startswith("#"):
            fh.seek(0)
        return list(csv.DictReader(fh))


def orig_by_text(pts, items, ax):
    """The unmoved axes-fraction position of each label, in the order of pts."""
    pos = {}
    for x, y, text in items:
        fx, fy = ax.transLimits.transform(ax.transScale.transform((x, y)))
        pos[text] = (fy, fx, text)
    return [pos[p[2]] for p in pts]


def place_end_labels(ax, items, t, min_gap=0.055, dx=6, push="down"):
    """Write one label at each series end, pushed apart so none overlaps.

    Positions are resolved in axes fractions, so the de-collision is correct on a
    log axis too, where equal data gaps are not equal screen gaps.
    """
    pts = []
    for x, y, text in items:
        fx, fy = ax.transLimits.transform(ax.transScale.transform((x, y)))
        pts.append([fy, fx, text])
    pts.sort(reverse=(push == "down"))
    sign = -1 if push == "down" else 1
    for i in range(1, len(pts)):
        if (pts[i - 1][0] - pts[i][0]) * sign < min_gap:
            pts[i][0] = pts[i - 1][0] + sign * min_gap
    moved = {id(p): p for p in pts}
    for (fy, fx, text), (oy, ox, _) in zip(pts, orig_by_text(pts, items, ax)):
        ax.annotate(text, (fx, fy), xycoords="axes fraction",
                    textcoords="offset points", xytext=(dx, 0),
                    fontsize=7.5, color=t["ink2"], va="center", ha="left",
                    annotation_clip=False)
        if abs(fy - oy) > 1e-6:
            # The label was pushed off its mark to avoid an overlap. A thin
            # leader says which mark it belongs to.
            ax.plot([ox, fx], [oy, fy], transform=ax.transAxes,
                    color=t["ink3"], linewidth=0.6, zorder=1,
                    clip_on=False)


def thousands(x, _pos=None):
    if x >= 1_000_000:
        return f"{x / 1_000_000:g}M"
    if x >= 1000:
        return f"{x / 1000:g}k"
    return f"{x:g}"


def style(fig, axes, t):
    fig.patch.set_facecolor(t["surface"])
    for ax in axes:
        ax.set_facecolor(t["surface"])
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(t["grid"])
            ax.spines[side].set_linewidth(0.8)
        ax.tick_params(colors=t["ink2"], labelsize=8, length=3, width=0.8)
        ax.xaxis.label.set_color(t["ink2"])
        ax.yaxis.label.set_color(t["ink2"])
        ax.title.set_color(t["ink"])
        ax.grid(True, color=t["grid"], linewidth=0.7, alpha=0.9)
        ax.set_axisbelow(True)


def save(fig, name, mode, t):
    os.makedirs(OUT, exist_ok=True)
    for ax in fig.axes:
        for art in (ax.title, ax._left_title, ax._right_title):
            art.set_color(t["ink"])
        ax.xaxis.label.set_color(t["ink2"])
        ax.yaxis.label.set_color(t["ink2"])
    fig.savefig(os.path.join(OUT, f"{name}.{mode}.png"), dpi=200,
                facecolor=t["surface"], bbox_inches="tight")
    if mode == "light":
        fig.savefig(os.path.join(OUT, f"{name}.pdf"),
                    facecolor=t["surface"], bbox_inches="tight")
    plt.close(fig)


def footer(fig, text, t):
    fig.text(0.0, -0.015, text, fontsize=7, color=t["ink3"], ha="left", va="top")


# ---------------------------------------------------------------- figure 1
def fig_cumulative(mode, t):
    rows = read(os.path.join(SUMMARY, "by_year.csv.gz"))
    by_year = defaultdict(int)
    per_sp = defaultdict(lambda: defaultdict(int))
    totals = defaultdict(int)
    for r in rows:
        y, n = int(r["year"]), int(r["runs"])
        if not 2008 <= y <= 2026:
            continue
        by_year[y] += n
        per_sp[r["species_name"]][y] += n
        totals[r["species_name"]] += n

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    style(fig, (ax1, ax2), t)

    years = sorted(by_year)
    cum, run = [], 0
    for y in years:
        run += by_year[y]
        cum.append(run)
    ax1.plot(years, cum, color=t["seq"], linewidth=2.0, solid_capstyle="round")
    ax1.fill_between(years, cum, color=t["seq"], alpha=0.12, linewidth=0)
    ax1.annotate(f"{run:,} runs", (years[-1], cum[-1]), textcoords="offset points",
                 xytext=(-6, 8), ha="right", fontsize=9, color=t["ink"], weight="bold")
    ax1.set_title("All fungal WGS runs, cumulative", fontsize=10, loc="left", pad=10)
    ax1.set_xlabel("first public (year)")
    ax1.yaxis.set_major_formatter(FuncFormatter(thousands))

    top = sorted(totals, key=totals.get, reverse=True)[:6]
    ends = []
    for i, sp in enumerate(top):
        c = t["series"][i]
        ys, run = sorted(per_sp[sp]), 0
        cv = []
        for y in ys:
            run += per_sp[sp][y]
            cv.append(run)
        ax2.plot(ys, cv, color=c, linewidth=2.0, solid_capstyle="round", label=sp)
        ends.append((ys[-1], cv[-1], f"{sp.split()[0][0]}. {' '.join(sp.split()[1:])}"))
    ax2.set_yscale("log")
    ax2.set_xlim(2008, 2027)
    place_end_labels(ax2, ends, t)
    ax2.set_title("Six largest species, cumulative (log)", fontsize=10, loc="left", pad=10)
    ax2.set_xlabel("first public (year)")
    ax2.yaxis.set_major_formatter(FuncFormatter(thousands))
    for ax in (ax1, ax2):
        ax.xaxis.set_major_locator(MultipleLocator(4))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, p: f"{int(v)}"))
    leg = ax2.legend(frameon=False, fontsize=7, loc="upper left", labelcolor=t["ink2"])
    for txt in leg.get_texts():
        txt.set_color(t["ink2"])

    footer(fig, f"n = {sum(by_year.values()):,} runs with a public date, 2008-2026. "
                f"{SNAPSHOT}. Direct labels repeat the legend; identity is never colour alone.", t)
    save(fig, "01_cumulative_runs", mode, t)
    return sum(by_year.values())


# ---------------------------------------------------------------- figure 2
def fig_rank(mode, t):
    rows = read(os.path.join(KINGDOM, "fungi_species_summary.csv.gz"))
    iso = sorted((int(r["isolates"]) for r in rows), reverse=True)
    rank = range(1, len(iso) + 1)

    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    style(fig, (ax,), t)
    ax.plot(rank, iso, color=t["seq"], linewidth=2.0, solid_capstyle="round")
    ax.set_xscale("log")
    ax.set_yscale("log")

    n100 = sum(1 for v in iso if v >= 100)
    n10k = sum(1 for v in iso if v >= 10000)
    ax.axhline(100, color=t["ink3"], linewidth=1.0, linestyle=(0, (4, 3)))
    ax.annotate(f"100 isolates — {n100} species at or above",
                (len(iso), 100), textcoords="offset points", xytext=(-4, 6),
                ha="right", fontsize=8, color=t["ink2"])
    ax.annotate(f"{n10k} species over 10,000", (1, iso[0]),
                textcoords="offset points", xytext=(6, -4),
                fontsize=8, color=t["ink2"], va="top")

    ax.set_xlabel("species rank (log)")
    ax.set_ylabel("distinct isolates (log)")
    ax.set_title("Isolates per species: a few huge, a long tail of ones",
                 fontsize=10, loc="left", pad=10)
    ax.yaxis.set_major_formatter(FuncFormatter(thousands))
    ax.xaxis.set_major_formatter(FuncFormatter(thousands))
    ax.xaxis.set_major_locator(LogLocator(base=10))

    footer(fig, f"n = {len(iso):,} species taxids, {sum(iso):,} isolates. "
                f"{SNAPSHOT}. One series, so no legend.", t)
    save(fig, "02_isolates_rank", mode, t)
    return len(iso), n100, n10k


# ---------------------------------------------------------------- figure 3
def fig_coverage(mode, t, n=25):
    rows = read(os.path.join(SUMMARY, "coverage.csv.gz"))
    rows = [r for r in rows if r["median_cov"] and r["runs_with_estimate"]]
    rows.sort(key=lambda r: int(r["runs_with_estimate"]), reverse=True)
    rows = rows[:n][::-1]

    fig, ax = plt.subplots(figsize=(7.8, 0.26 * len(rows) + 1.8))
    style(fig, (ax,), t)
    y = range(len(rows))
    for i, r in enumerate(rows):
        lo, md, hi = float(r["p10_cov"]), float(r["median_cov"]), float(r["p90_cov"])
        ax.plot([lo, hi], [i, i], color=t["seq_light"], linewidth=3.2,
                solid_capstyle="round", zorder=2)
        ax.plot([md], [i], marker="o", markersize=5.5, color=t["seq"],
                markeredgecolor=t["surface"], markeredgewidth=1.4, zorder=3)
    ax.axvline(15, color=t["ink3"], linewidth=1.0, linestyle=(0, (4, 3)), zorder=1)
    ax.set_xscale("log")
    ax.set_yticks(list(y))
    ax.set_yticklabels([f"{r['species_name']}  ({int(r['runs_with_estimate']):,})"
                        for r in rows], fontsize=7.5, color=t["ink2"])
    ax.set_xlabel("estimated depth, x (log) — bar is p10 to p90, dot is the median")
    ax.set_title("Sequencing depth, 25 species with the most runs",
                 fontsize=10, loc="left", pad=10)
    ax.annotate(" 15x", (15, len(rows) - 0.45), fontsize=8, color=t["ink2"],
                ha="left", va="bottom")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, p: f"{v:g}"))
    ax.grid(axis="y", visible=False)

    tot = sum(int(r["runs_with_estimate"]) for r in rows)
    footer(fig, f"n = {tot:,} runs with a base count and a reference genome size, "
                f"across {len(rows)} species. Count per species in brackets. "
                f"{SNAPSHOT}.", t)
    save(fig, "03_coverage", mode, t)
    return tot


# ---------------------------------------------------------------- figure 4
def fig_clonality(mode, t):
    rows = read(os.path.join(SUMMARY, "clonality.csv.gz"))
    rows = [r for r in rows if int(r["isolates"]) >= 100]
    x = [int(r["isolates"]) for r in rows]
    y = [int(r["largest_study"]) / int(r["isolates"]) * 100 for r in rows]

    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    style(fig, (ax,), t)
    ax.scatter(x, y, s=26, color=t["seq"], alpha=0.55,
               edgecolors=t["surface"], linewidths=0.8, zorder=3)
    ax.set_xscale("log")
    ax.set_xlabel("distinct isolates in the species (log)")
    ax.set_ylabel("largest single study, % of the species")
    ax.set_ylim(0, 103)
    ax.set_title("How much of each species comes from one study",
                 fontsize=10, loc="left", pad=10)
    ax.xaxis.set_major_formatter(FuncFormatter(thousands))

    label = sorted(rows, key=lambda r: int(r["isolates"]), reverse=True)[:5]
    place_end_labels(ax, [(int(r["isolates"]),
                           int(r["largest_study"]) / int(r["isolates"]) * 100,
                           r["species_name"]) for r in label], t,
                     min_gap=0.055, push="up")

    half = sum(1 for v in y if v >= 50)
    footer(fig, f"n = {len(rows)} species with 100 or more isolates. "
                f"{half} of them take half or more of their isolates from a single study, "
                f"so a run count is not a strain count. {SNAPSHOT}.", t)
    save(fig, "04_clonality", mode, t)
    return len(rows), half


# ---------------------------------------------------------------- figure 5
def fig_geography(mode, t, n=25):
    rows = read(os.path.join(SUMMARY, "geography.csv.gz"))
    # Placeholder strings are not places. They are counted and reported in the
    # footer, not drawn as if they were a country.
    NOT_A_PLACE = {"missing", "not provided", "not collected", "not applicable",
                   "unknown", "na", "none", "null", "n/a", "uncalculated",
                   "not determined", "restricted access"}
    tot = defaultdict(int)
    placeholder = 0
    for r in rows:
        c = r["country"].strip()
        if c.lower() in NOT_A_PLACE:
            placeholder += int(r["isolates"])
            continue
        tot[c] += int(r["isolates"])
    top = sorted(tot.items(), key=lambda kv: kv[1], reverse=True)[:n][::-1]

    fig, ax = plt.subplots(figsize=(7.2, 0.25 * len(top) + 1.8))
    style(fig, (ax,), t)
    ax.barh([k for k, _ in top], [v for _, v in top], height=0.62,
            color=t["seq"], zorder=3)
    for i, (_, v) in enumerate(top):
        ax.annotate(f" {v:,}", (v, i), fontsize=7.5, color=t["ink2"], va="center")
    ax.set_xlabel("isolates with a recorded country")
    ax.set_title(f"Where the isolates come from, top {n} countries",
                 fontsize=10, loc="left", pad=10)
    ax.tick_params(axis="y", labelsize=7.5)
    ax.xaxis.set_major_formatter(FuncFormatter(thousands))
    ax.grid(axis="y", visible=False)
    ax.set_xlim(0, max(v for _, v in top) * 1.12)

    placed = sum(tot.values())
    footer(fig, f"n = {placed:,} isolates with a named country, across {len(tot):,} "
                f"countries; the top {n} hold {sum(v for _, v in top) / placed * 100:.0f}% "
                f"of them. A further {placeholder:,} isolates carry a placeholder "
                f'("missing", "not provided") and are not drawn. {SNAPSHOT}.', t)
    save(fig, "05_geography", mode, t)
    return placed, len(tot)


def main():
    stats = {}
    for mode, t in THEMES.items():
        stats["runs"] = fig_cumulative(mode, t)
        stats["rank"] = fig_rank(mode, t)
        stats["cov"] = fig_coverage(mode, t)
        stats["clon"] = fig_clonality(mode, t)
        stats["geo"] = fig_geography(mode, t)
    for k, v in stats.items():
        print(f"{k}: {v}", file=sys.stderr)
    print(f"wrote {len(os.listdir(OUT))} files to figures/", file=sys.stderr)


if __name__ == "__main__":
    main()
