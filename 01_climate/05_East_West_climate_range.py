#!/usr/bin/env python3
"""
Author: Iris Sammarco
Date: 16/10/2026

Compare climatic variation between eastern and western populations of
Anastatica using Climate_PC1, Climate_PC2, and Climate_PC3.

PURPOSE
-------
This script assesses whether eastern and western populations cover
comparable ranges of climatic variation, helping evaluate whether
differences in climate-associated gene expression might coincide with
differences in climatic breadth.

INPUT
-----
/Anastatica/metadata/population_metadata.csv

Grouping variable:
    Site_of_distribution_range_2-level  (East / West)

Climate variables:
    Climate_PC1, Climate_PC2, Climate_PC3

OUTPUT
------
Saved under:
/Anastatica/climate/East_West

Tables and report:
    1. climate_summary_by_region.csv
    2. east_west_climate_comparisons.csv
    3. climate_range_overlap.csv
    4. climate_range_analysis_report.txt

Figures (each saved as both PNG and PDF):
    5. climate_PC_distributions.(png|pdf)   boxplots
    6. climate_PC_ranges.(png|pdf)          observed ranges
    7. climate_PC_merged.(png|pdf)          boxplots (a) on top,
                                            observed ranges (b) below

ANALYSIS
--------
- Uses one climate record per population.
- Reports descriptive statistics for each region and climate PC.
- Performs Welch's t-test and Mann-Whitney U test.
- Calculates the observed range intersection and range overlap ratio.
- Produces plots showing distributions and observed ranges.

The tests are exploratory and do not establish causality.
The range overlap ratio describes overlap of observed intervals,
not the overlap of the full climate distributions.

REQUIREMENTS
------------
Python 3, pandas, numpy, scipy, matplotlib.

Install missing dependencies with:
    pip install pandas numpy scipy matplotlib

Run with:
mamba activate python
python 05_East_West_climate_range.py
"""

from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import ttest_ind, mannwhitneyu


# ============================================================
# 1. PATHS AND SETTINGS
# ============================================================

INPUT_FILE = Path(
    "/Anastatica/metadata/population_metadata.csv"
)

OUTPUT_DIR = Path(
    "/Anastatica/climate/East_West"
)

GROUP_COLUMN = "Site_of_distribution_range_2-level"
POPULATION_COLUMN = "Population"

CLIMATE_COLUMNS = [
    "Climate_PC1",
    "Climate_PC2",
    "Climate_PC3",
]

REGIONS = ["East", "West"]

# Colours used for each region in ALL figures.
# REGION_COLORS: boxes, points, range lines and range end caps.
# REGION_DOT_COLORS: darker shade of the same hue, used for the median
# dot on the range plot so it stays visible on top of the line.
REGION_COLORS = {
    "East": "#1f77b4",   # blue
    "West": "#ff7f0e",   # orange
}
REGION_DOT_COLORS = {
    "East": "#0b3c5d",   # dark blue
    "West": "#a04a00",   # dark orange
}

# Figure formats written for every figure.
FIGURE_FORMATS = ["png", "pdf"]

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# 2. LOAD AND VALIDATE DATA
# ============================================================

if not INPUT_FILE.exists():
    raise FileNotFoundError(f"Input file not found: {INPUT_FILE}")

metadata = pd.read_csv(INPUT_FILE)

required_columns = [
    POPULATION_COLUMN,
    GROUP_COLUMN,
    *CLIMATE_COLUMNS,
]

missing_columns = [
    column for column in required_columns
    if column not in metadata.columns
]

if missing_columns:
    raise ValueError(
        f"Missing required columns: {missing_columns}"
    )

metadata[GROUP_COLUMN] = (
    metadata[GROUP_COLUMN].astype("string").str.strip()
)

for column in CLIMATE_COLUMNS:
    metadata[column] = pd.to_numeric(
        metadata[column], errors="coerce"
    )

# Retain only the requested East and West categories.
data = metadata[
    metadata[GROUP_COLUMN].isin(REGIONS)
].copy()

if data.empty:
    raise ValueError(
        f"No East/West records found in {GROUP_COLUMN}."
    )

if data[POPULATION_COLUMN].isna().any():
    raise ValueError(
        "Missing population identifiers found in East/West records. "
        "Resolve these before analysis."
    )

# Check whether duplicated population records have conflicting
# region assignments or climate values.
for population, group in data.groupby(POPULATION_COLUMN):
    if group[GROUP_COLUMN].nunique() > 1:
        raise ValueError(
            f"Population {population} is assigned to multiple regions."
        )

    for column in CLIMATE_COLUMNS:
        n_unique = group[column].dropna().nunique()
        if n_unique > 1:
            raise ValueError(
                f"Population {population} has conflicting values for "
                f"{column}. Inspect the input metadata before proceeding."
            )

# Climate is a population-level variable. Collapse repeated rows to
# one record per population, retaining the unique climate value.
population_data = (
    data.groupby(POPULATION_COLUMN, as_index=False)
    .agg({
        GROUP_COLUMN: "first",
        **{column: "first" for column in CLIMATE_COLUMNS},
    })
)

print("\nPopulation counts by region:")
print(
    population_data[GROUP_COLUMN]
    .value_counts()
    .reindex(REGIONS, fill_value=0)
)

# Warn if climate PCs contain missing values.
for column in CLIMATE_COLUMNS:
    n_missing = population_data[column].isna().sum()
    if n_missing:
        warnings.warn(
            f"{column}: {n_missing} populations have missing values. "
            "They will be excluded from analyses of this PC."
        )


# ============================================================
# 3. DESCRIPTIVE STATISTICS
# ============================================================

summary_rows = []

for climate in CLIMATE_COLUMNS:
    for region in REGIONS:
        values = population_data.loc[
            population_data[GROUP_COLUMN] == region, climate
        ].dropna()

        if values.empty:
            continue

        summary_rows.append({
            "Climate_PC": climate,
            "Region": region,
            "n_populations": len(values),
            "mean": values.mean(),
            "sd": values.std(ddof=1),
            "median": values.median(),
            "Q1": values.quantile(0.25),
            "Q3": values.quantile(0.75),
            "IQR": values.quantile(0.75) - values.quantile(0.25),
            "minimum": values.min(),
            "maximum": values.max(),
            "range": values.max() - values.min(),
        })

summary = pd.DataFrame(summary_rows)

summary.to_csv(
    OUTPUT_DIR / "climate_summary_by_region.csv",
    index=False,
    float_format="%.6f",
)


# ============================================================
# 4. EAST-WEST COMPARISONS
# ============================================================

comparison_rows = []

for climate in CLIMATE_COLUMNS:
    east = population_data.loc[
        population_data[GROUP_COLUMN] == "East", climate
    ].dropna().to_numpy()

    west = population_data.loc[
        population_data[GROUP_COLUMN] == "West", climate
    ].dropna().to_numpy()

    if len(east) < 2 or len(west) < 2:
        warnings.warn(
            f"Insufficient data for statistical tests of {climate}."
        )
        continue

    # Welch's t-test does not assume equal group variances.
    welch = ttest_ind(east, west, equal_var=False)

    # Two-sided Mann-Whitney U test.
    mw = mannwhitneyu(
        east, west, alternative="two-sided"
    )

    east_mean = np.mean(east)
    west_mean = np.mean(west)
    east_sd = np.std(east, ddof=1)
    west_sd = np.std(west, ddof=1)

    # Ratio of observed ranges; descriptive only.
    east_range = np.ptp(east)
    west_range = np.ptp(west)

    if east_range > 0:
        range_ratio_west_east = west_range / east_range
    else:
        range_ratio_west_east = np.nan

    comparison_rows.append({
        "Climate_PC": climate,
        "n_East": len(east),
        "n_West": len(west),
        "mean_East": east_mean,
        "mean_West": west_mean,
        "mean_difference_East_minus_West": east_mean - west_mean,
        "sd_East": east_sd,
        "sd_West": west_sd,
        "range_East": east_range,
        "range_West": west_range,
        "range_ratio_West_over_East": range_ratio_west_east,
        "Welch_t": welch.statistic,
        "Welch_p": welch.pvalue,
        "Mann_Whitney_U": mw.statistic,
        "Mann_Whitney_p": mw.pvalue,
    })

comparisons = pd.DataFrame(comparison_rows)

# Adjust the three p-values within each test family using
# the Benjamini-Hochberg procedure.
def benjamini_hochberg(p_values):
    """Return Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(p_values, dtype=float)
    result = np.full(p.shape, np.nan)

    valid = np.isfinite(p)
    pv = p[valid]

    if len(pv) == 0:
        return result

    order = np.argsort(pv)
    sorted_p = pv[order]
    m = len(sorted_p)

    adjusted = sorted_p * m / np.arange(1, m + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)

    restored = np.empty(m)
    restored[order] = adjusted
    result[valid] = restored

    return result


if not comparisons.empty:
    comparisons["Welch_FDR"] = benjamini_hochberg(
        comparisons["Welch_p"]
    )
    comparisons["Mann_Whitney_FDR"] = benjamini_hochberg(
        comparisons["Mann_Whitney_p"]
    )

comparisons.to_csv(
    OUTPUT_DIR / "east_west_climate_comparisons.csv",
    index=False,
    float_format="%.6g",
)


# ============================================================
# 5. OBSERVED CLIMATE-RANGE OVERLAP
# ============================================================

overlap_rows = []

for climate in CLIMATE_COLUMNS:
    east = population_data.loc[
        population_data[GROUP_COLUMN] == "East", climate
    ].dropna().to_numpy()

    west = population_data.loc[
        population_data[GROUP_COLUMN] == "West", climate
    ].dropna().to_numpy()

    if len(east) == 0 or len(west) == 0:
        continue

    east_min, east_max = np.min(east), np.max(east)
    west_min, west_max = np.min(west), np.max(west)

    intersection_min = max(east_min, west_min)
    intersection_max = min(east_max, west_max)

    intersection_width = max(
        0.0, intersection_max - intersection_min
    )

    union_min = min(east_min, west_min)
    union_max = max(east_max, west_max)
    union_width = union_max - union_min

    # Intersection length divided by union length.
    # This is an interval-overlap metric, not a distribution-overlap
    # coefficient. It is undefined if both groups have zero range.
    overlap_ratio = (
        intersection_width / union_width
        if union_width > 0 else np.nan
    )

    overlap_rows.append({
        "Climate_PC": climate,
        "East_min": east_min,
        "East_max": east_max,
        "East_range": east_max - east_min,
        "West_min": west_min,
        "West_max": west_max,
        "West_range": west_max - west_min,
        "intersection_min": intersection_min,
        "intersection_max": intersection_max,
        "intersection_width": intersection_width,
        "union_width": union_width,
        "range_overlap_ratio": overlap_ratio,
        "range_overlap_percent": (
            100 * overlap_ratio
            if np.isfinite(overlap_ratio) else np.nan
        ),
    })

overlap = pd.DataFrame(overlap_rows)

overlap.to_csv(
    OUTPUT_DIR / "climate_range_overlap.csv",
    index=False,
    float_format="%.6f",
)


# ============================================================
# 6. PLOTTING HELPERS
# ============================================================

plt.rcParams.update({
    "font.size": 11,
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,   # keep text editable in PDF
    "ps.fonttype": 42,
})


def save_figure(fig, basename):
    """Save a figure in every format listed in FIGURE_FORMATS."""
    for fmt in FIGURE_FORMATS:
        fig.savefig(
            OUTPUT_DIR / f"{basename}.{fmt}",
            bbox_inches="tight",
        )
    plt.close(fig)


def add_panel_letter(ax, letter):
    """Add a bold panel letter at the top-left of an axis."""
    ax.text(
        -0.14, 1.08,
        letter,
        transform=ax.transAxes,
        fontsize=18,
        fontweight="bold",
        va="bottom",
        ha="left",
    )


def draw_distribution(ax, climate, show_title=True):
    """Boxplot with individual populations for one climate PC."""
    values_by_region = [
        population_data.loc[
            population_data[GROUP_COLUMN] == region, climate
        ].dropna().to_numpy()
        for region in REGIONS
    ]

    bp = ax.boxplot(
        values_by_region,
        tick_labels=REGIONS,
        showfliers=False,
        widths=0.5,
    )

    # Colour every element of each box (box, median, whiskers, caps)
    # with the region colour. Each box has 2 whiskers and 2 caps.
    for i, region in enumerate(REGIONS):
        colour = REGION_COLORS[region]
        bp["boxes"][i].set_color(colour)
        bp["medians"][i].set_color(colour)
        for line in bp["whiskers"][2 * i:2 * i + 2]:
            line.set_color(colour)
        for line in bp["caps"][2 * i:2 * i + 2]:
            line.set_color(colour)

    # Show each population as a point (fixed seed -> same jitter in
    # every figure).
    rng = np.random.default_rng(42)

    for position, (region, values) in enumerate(
        zip(REGIONS, values_by_region), start=1
    ):
        jitter = rng.uniform(-0.07, 0.07, size=len(values))
        ax.scatter(
            position + jitter,
            values,
            color=REGION_COLORS[region],
            alpha=0.75,
            s=22,
            zorder=3,
        )

    if show_title:
        ax.set_title(climate.replace("_", " "))
    ax.set_ylabel("Climate PC score")
    ax.set_xlabel("Distribution range")
    ax.grid(axis="y", alpha=0.25)


def draw_ranges(ax, climate, show_title=True):
    """Observed min-max range (line) and median (dot) per region."""
    for position, region in enumerate(REGIONS):
        values = population_data.loc[
            population_data[GROUP_COLUMN] == region, climate
        ].dropna().to_numpy()

        if len(values) == 0:
            continue

        low, high = values.min(), values.max()
        median = np.median(values)

        colour = REGION_COLORS[region]
        dot_colour = REGION_DOT_COLORS[region]

        # Range line (minimum to maximum).
        ax.hlines(
            y=position,
            xmin=low,
            xmax=high,
            color=colour,
            linewidth=4,
        )
        # End caps at minimum and maximum.
        ax.plot(
            [low, high],
            [position, position],
            "|",
            color=colour,
            markersize=12,
            markeredgewidth=2,
        )
        # Median dot (darker shade of the same colour).
        ax.plot(
            median,
            position,
            "o",
            color=dot_colour,
            markersize=7,
            zorder=3,
        )

    ax.set_yticks([0, 1])
    ax.set_yticklabels(REGIONS)
    ax.set_ylim(-0.5, 1.5)
    if show_title:
        ax.set_title(climate.replace("_", " "))
    ax.set_xlabel("Climate PC score")
    ax.set_ylabel("Distribution range")
    ax.grid(axis="x", alpha=0.25)


# ============================================================
# 7. FIGURE: CLIMATE DISTRIBUTIONS (BOXPLOTS)
# ============================================================

fig, axes = plt.subplots(
    1, 3, figsize=(13, 4.5), constrained_layout=True
)

for ax, climate in zip(axes, CLIMATE_COLUMNS):
    draw_distribution(ax, climate)

save_figure(fig, "climate_PC_distributions")


# ============================================================
# 8. FIGURE: OBSERVED CLIMATE RANGES
# ============================================================

fig, axes = plt.subplots(
    1, 3, figsize=(13, 4.5), constrained_layout=True
)

for ax, climate in zip(axes, CLIMATE_COLUMNS):
    draw_ranges(ax, climate)

save_figure(fig, "climate_PC_ranges")


# ============================================================
# 9. MERGED FIGURE: (a) BOXPLOTS ON TOP, (b) RANGES BELOW
# ============================================================
# Lines in panel b = minimum to maximum observed range; dots = medians.

fig, axes = plt.subplots(
    2, 3, figsize=(13, 9), constrained_layout=True
)

for col, climate in enumerate(CLIMATE_COLUMNS):
    # Top row: boxplots (PC names shown as column titles).
    draw_distribution(axes[0, col], climate, show_title=True)
    # Bottom row: observed ranges (column titles already on top row).
    draw_ranges(axes[1, col], climate, show_title=False)

add_panel_letter(axes[0, 0], "a")
add_panel_letter(axes[1, 0], "b")

save_figure(fig, "climate_PC_merged")


# ============================================================
# 10. WRITE A TEXT REPORT
# ============================================================

report_path = OUTPUT_DIR / "climate_range_analysis_report.txt"

with open(report_path, "w", encoding="utf-8") as report:
    report.write(
        "EAST-WEST CLIMATE RANGE ANALYSIS\n"
        "================================\n\n"
    )

    report.write(
        f"Input file: {INPUT_FILE}\n"
        f"Output directory: {OUTPUT_DIR}\n\n"
    )

    report.write("Population counts:\n")

    for region in REGIONS:
        n = (
            population_data[GROUP_COLUMN] == region
        ).sum()
        report.write(f"  {region}: {n} populations\n")

    report.write(
        "\nDESCRIPTIVE STATISTICS\n"
        "----------------------\n"
    )

    report.write(summary.to_string(index=False))
    report.write("\n\n")

    report.write(
        "EAST-WEST STATISTICAL COMPARISONS\n"
        "---------------------------------\n"
    )

    if comparisons.empty:
        report.write("No valid statistical comparisons available.\n")
    else:
        report.write(comparisons.to_string(index=False))

    report.write("\n\nOBSERVED RANGE OVERLAP\n")
    report.write("----------------------\n")

    if overlap.empty:
        report.write("No valid range comparisons available.\n")
    else:
        report.write(overlap.to_string(index=False))

    report.write(
        "\n\nINTERPRETATION NOTES\n"
        "--------------------\n"
        "- Range is maximum minus minimum. Range overlap describes "
        "observed intervals, not full probability distributions.\n"
        "- Non-significant tests do not establish equivalence between "
        "regions.\n"
        "- These analyses describe climate distributions and do not "
        "establish the cause of differences in gene expression.\n"
        "- Climate PC scores are interpreted along their supplied axes; "
        "PC1, PC2, and PC3 represent distinct climatic gradients.\n"
    )

print("\nAnalysis completed.")
print(f"Results saved to: {OUTPUT_DIR}")
print("\nGenerated files:")

for path in sorted(OUTPUT_DIR.iterdir()):
    if path.is_file():
        print(f"  {path.name}")