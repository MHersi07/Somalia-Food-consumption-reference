"""
sihbs_consumption.py
====================

Core analysis pipeline for:

    Hersi, M. A. & Fiidow, O. A.
    "A National Food Consumption Reference for Food Safety Risk Assessment
     in Somalia"

Implements, in the order described in the manuscript Methods:

    1. Unit-code resolution and non-standard-unit (NSU) conversion
    2. Data-quality screening for implausible quantities
    3. Per-capita consumption calculation (among consumers)
    4. Person-weighted population-level per-capita quantities
    5. Survey weighting and stratified variance estimation
    6. Design effect (DEFF) sensitivity analysis and hypothesis testing

The module is deliberately dependency-light (pandas, numpy, scipy) and has
no side effects on import.

Data source
-----------
Somalia Integrated Household Budget Survey (SIHBS) 2022, Somali National
Bureau of Statistics, available from the SNBS microdata portal
(https://microdata.nbs.gov.so) subject to registration and approval by SNBS.
The microdata are NOT redistributed in this repository; see README.md for
access and preparation instructions.

Licence: MIT (code). See LICENSE.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

# ----------------------------------------------------------------------
# Constants documented in the manuscript Methods
# ----------------------------------------------------------------------

#: Recall period of the SIHBS food consumption module, in days.
RECALL_DAYS = 7

#: Implausibility screen: maximum per-capita daily quantity for a single
#: SOLID item, in grams. Total daily food intake across all items is of the
#: order of 1-2 kg, so a single item exceeding 1 kg/person/day cannot be
#: correct. Sits well above the 99.9th percentile of the observed
#: distribution (686 g/cap/day).
SOLID_IMPLAUSIBLE_G = 1000.0

#: Implausibility screen: maximum per-capita daily quantity for a single
#: LIQUID item, in millilitres, reflecting the upper limit of physiological
#: fluid intake. Deliberately conservative: every liquid item other than
#: mineral water has an observed maximum below 1,500 mL/cap/day.
LIQUID_IMPLAUSIBLE_ML = 5000.0

#: Average number of sampled households per Enumeration Area (601 EAs,
#: ~12 households each), used for the Kish design effect.
CLUSTER_SIZE_M = 12

#: Primary assumed intra-cluster correlation for the sensitivity analysis.
PRIMARY_ICC = 0.10

#: The eleven purposively selected headline items reported in the main text.
#: Names follow the official SIHBS item-code labels verbatim (including the
#: source codebook's own capitalisation inconsistencies) so that they remain
#: valid join keys against the SNBS microdata.
HEADLINE_ITEMS = [
    "Sugar - white",
    "Cooking oil",
    "Onion",
    "Wheat flour - white",
    "Black tea (in bulk)",
    "Rice - White [Local]",
    "Goat Meat - Fresh",
    "Camel meat - fresh",
    "Maize grain",
    "Camel Milk - Fresh",
    "Sorghum grain",
]

#: All meat item codes, used for the "any meat" comparator in the livestock
#: linkage analysis.
MEAT_ITEMS = [
    "Goat Meat - Fresh",
    "Camel meat - fresh",
    "Chicken Meat - frozen",
    "Sheep Meat - Fresh",
    "Chicken Meat - fresh",
    "Cow Meat [Beef] - Fresh",
    "Shark meat, fresh",
    "Other Meat products (specify)",
]


# ----------------------------------------------------------------------
# Loading
# ----------------------------------------------------------------------

def load_sihbs(path: str) -> dict[str, pd.DataFrame]:
    """Load the prepared SIHBS analysis file.

    Parameters
    ----------
    path
        Path to the prepared SIHBS 2022 analysis file (see README, section 3,
        for how to prepare this from the SIHBS microdata release).

    Returns
    -------
    dict
        Keys: ``households``, ``food``, ``items``, ``livestock``.
    """
    return {
        "households": pd.read_excel(path, sheet_name="Household_Master"),
        "food": pd.read_excel(path, sheet_name="Food_Consumption_Long"),
        "items": pd.read_excel(path, sheet_name="Food_Item_Dictionary"),
        "livestock": pd.read_excel(path, sheet_name="Livestock_Own"),
    }


# ----------------------------------------------------------------------
# Step 1-2: resolution and implausibility screening
# ----------------------------------------------------------------------

def build_household_item_table(
    food: pd.DataFrame,
    households: pd.DataFrame,
    solid_max: float = SOLID_IMPLAUSIBLE_G,
    liquid_max: float = LIQUID_IMPLAUSIBLE_ML,
) -> pd.DataFrame:
    """Aggregate the long food file to one row per household-item.

    Resolved gram and millilitre quantities are summed within household-item,
    per-capita daily quantities are derived, and the pre-specified
    implausibility screen is applied as a boolean ``implausible`` column.

    The screen flags QUANTITIES only. Flagged households are deliberately
    retained for coverage estimation, because the household did consume the
    item and only the reported quantity is in error; dropping them from the
    coverage denominator would introduce a different bias.
    """
    pos = food[food["cr15_04quantity"] > 0].copy()
    pos["g"] = pos["resolved_grams"].fillna(0)
    pos["mL"] = pos["resolved_mL"].fillna(0)

    agg = (
        pos.groupby(["hid", "item", "item_name"], as_index=False)
        .agg(g=("g", "sum"), mL=("mL", "sum"))
    )

    hhsize = households.set_index("hid")["hhsize"]
    agg["hhsize"] = agg["hid"].map(hhsize)
    agg["pc_g"] = agg["g"] / (agg["hhsize"] * RECALL_DAYS)
    agg["pc_mL"] = agg["mL"] / (agg["hhsize"] * RECALL_DAYS)
    agg["implausible"] = (agg["pc_g"] > solid_max) | (agg["pc_mL"] > liquid_max)
    return agg


def screening_summary(agg: pd.DataFrame) -> pd.DataFrame:
    """Per-item count of records removed by the implausibility screen."""
    flagged = agg[agg["implausible"]]
    out = (
        flagged.groupby("item_name")
        .size()
        .sort_values(ascending=False)
        .rename("n_records_excluded")
        .reset_index()
    )
    return out


# ----------------------------------------------------------------------
# Weighted estimators
# ----------------------------------------------------------------------

def weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    return float((values * weights).sum() / weights.sum())


def weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    """Weighted quantile via the inverse of the cumulative-weight function.

    Note on definition: this is the standard survey-weighted quantile (inverse
    empirical CDF with linear interpolation on cumulative weight). It differs
    marginally from ``numpy.median``'s convention, which averages the two
    central order statistics for even n: for [1..10] with equal weights this
    returns 5.0 where numpy returns 5.5. The difference is negligible at the
    sample sizes used here (thousands of households) and this definition is the
    one used for all published medians and P95 values.
    """
    order = np.argsort(values)
    v, w = np.asarray(values)[order], np.asarray(weights)[order]
    cum = np.cumsum(w) / w.sum()
    return float(np.interp(q, cum, v))


def weighted_proportion(
    group: pd.DataFrame,
    member_ids: set,
    weight_col: str = "wgt",
    strata: tuple = ("region_n", "ea_type_n"),
) -> tuple[float, float]:
    """Weighted proportion of households in ``group`` that are in ``member_ids``.

    Returns ``(proportion, linearised_standard_error)``.

    The variance is the Taylor linearisation of the Hajek ratio estimator under
    stratified sampling::

        u_i      = w_i (y_i - p)
        Var(p)   = (1 / W^2) * SUM_h [ n_h/(n_h - 1) * SUM_{i in h} (u_i - ubar_h)^2 ]

    where h indexes the design strata (region x residence type; 49 non-empty
    strata, minimum size 48, so no singleton-stratum problem arises), and
    W = sum of weights.

    Stratification is incorporated because SIHBS is an explicitly stratified
    design and discarding that information inflates the variance without
    justification. **Clustering is not** incorporated: the public-use files omit
    the EA/PSU identifier, so it cannot be estimated. That omission is
    anti-conservative and is instead parameterised through the design effect;
    see :func:`deff_adjusted_test` and :func:`breakdown_deff`. Do not treat the
    unadjusted interval as a full design-based confidence interval.

    Passing ``strata=None`` falls back to the unstratified estimator, which is
    provided only for sensitivity comparison.
    """
    y = group["hid"].isin(member_ids).astype(float).to_numpy()
    w = group[weight_col].to_numpy()
    W = w.sum()
    p = (w * y).sum() / W
    u = w * (y - p)

    if strata is None or not all(s in group.columns for s in strata):
        return float(p), float(np.sqrt((u ** 2).sum()) / W)

    keys = list(zip(*[group[s] for s in strata]))
    frame = pd.DataFrame({"_k": keys, "_u": u})
    var = 0.0
    for _, g in frame.groupby("_k", sort=False):
        uh = g["_u"].to_numpy()
        nh = len(uh)
        if nh < 2:
            continue
        var += nh / (nh - 1) * ((uh - uh.mean()) ** 2).sum()
    return float(p), float(np.sqrt(var) / W)


# ----------------------------------------------------------------------
# Step 3-4: coverage and per-capita quantities
# ----------------------------------------------------------------------

def coverage(
    agg: pd.DataFrame, households: pd.DataFrame, item_name: str,
    weight_col: str = "wgt",
) -> dict:
    """Weighted household coverage (%) for one item, with unweighted n."""
    ids = set(agg.loc[agg["item_name"] == item_name, "hid"])
    p, se = weighted_proportion(households, ids, weight_col)
    return {
        "item_name": item_name,
        "n_consuming": len(ids),
        "coverage_pct": 100 * p,
        "se_pct": 100 * se,
    }


def among_consumers_quantity(
    agg: pd.DataFrame, households: pd.DataFrame, item_name: str, liquid: bool = False
) -> dict:
    """Per-capita daily quantity among consuming households (screened)."""
    col = "pc_mL" if liquid else "pc_g"
    d = agg[(agg["item_name"] == item_name) & (~agg["implausible"]) & (agg[col] > 0)]
    if d.empty:
        return {"item_name": item_name, "n": 0}
    w = households.set_index("hid")["wgt"].reindex(d["hid"]).to_numpy()
    v = d[col].to_numpy()
    return {
        "item_name": item_name,
        "n": int(len(d)),
        "mean": weighted_mean(v, w),
        "sd": float(np.sqrt(weighted_mean((v - weighted_mean(v, w)) ** 2, w))),
        "median": weighted_quantile(v, w, 0.50),
        "p95": weighted_quantile(v, w, 0.95),
    }


def population_quantity(
    agg: pd.DataFrame, households: pd.DataFrame, item_name: str, liquid: bool = False
) -> float:
    """Person-weighted population-level per-capita daily quantity.

    Implements Equation 2 of the manuscript::

        sum_h (w_h * Q_h) / (RECALL_DAYS * sum_h (w_h * n_h))

    where the denominator sums over ALL households, so that non-consuming
    households correctly contribute zero.
    """
    col = "mL" if liquid else "g"
    denom = RECALL_DAYS * (households["wgt"] * households["hhsize"]).sum()
    d = agg[(agg["item_name"] == item_name) & (~agg["implausible"]) & (agg[col] > 0)]
    if d.empty:
        return float("nan")
    w = households.set_index("hid")["wgt"].reindex(d["hid"]).to_numpy()
    return float((w * d[col].to_numpy()).sum() / denom)


# ----------------------------------------------------------------------
# Step 6: design effect sensitivity and hypothesis testing
# ----------------------------------------------------------------------

def deff_from_icc(icc: float, m: int = CLUSTER_SIZE_M) -> float:
    """Kish design effect: DEFF = 1 + (m - 1) * ICC."""
    return 1.0 + (m - 1) * icc


def deff_adjusted_test(
    group1: pd.DataFrame,
    group2: pd.DataFrame,
    member_ids: set,
    deff: float = None,
    weight_col: str = "wgt",
) -> dict:
    """Two-sided z test on the difference between two weighted proportions.

    Linearised variances are multiplied by ``deff`` to approximate the effect
    of the unobserved clustering. Returns the difference in percentage points,
    its confidence interval, the z statistic and the exact P value.
    """
    if deff is None:
        deff = deff_from_icc(PRIMARY_ICC)
    p1, s1 = weighted_proportion(group1, member_ids, weight_col)
    p2, s2 = weighted_proportion(group2, member_ids, weight_col)
    k = np.sqrt(deff)
    se_diff = np.sqrt((s1 * k) ** 2 + (s2 * k) ** 2)
    z = (p1 - p2) / se_diff
    p_value = float(2 * stats.norm.sf(abs(z)))
    return {
        "group1_pct": 100 * p1,
        "group2_pct": 100 * p2,
        "diff_pp": 100 * (p1 - p2),
        "ci_low_pp": 100 * ((p1 - p2) - 1.96 * se_diff),
        "ci_high_pp": 100 * ((p1 - p2) + 1.96 * se_diff),
        "z": float(z),
        "p_value": p_value,
        "deff": float(deff),
    }


def breakdown_deff(
    group1: pd.DataFrame,
    group2: pd.DataFrame,
    member_ids: set,
    weight_col: str = "wgt",
    alpha: float = 0.05,
    upper: float = 2000.0,
) -> float:
    """Design effect at which a comparison would cease to reach ``alpha``.

    Solved by bisection. A breakdown DEFF far above the plausible range
    (ICC 0.05-0.20, i.e. DEFF 1.55-3.20 at m = 12) indicates the comparison
    is robust to the absent clustering term. Returns ``upper`` if the
    comparison never loses significance below that bound.
    """
    p1, s1 = weighted_proportion(group1, member_ids, weight_col)
    p2, s2 = weighted_proportion(group2, member_ids, weight_col)

    def p_at(deff: float) -> float:
        k = np.sqrt(deff)
        z = (p1 - p2) / np.sqrt((s1 * k) ** 2 + (s2 * k) ** 2)
        return 2 * stats.norm.sf(abs(z))

    if p_at(1.0) >= alpha:
        return float("nan")  # already non-significant
    if p_at(upper) < alpha:
        return upper
    lo, hi = 1.0, upper
    for _ in range(200):
        mid = (lo + hi) / 2
        if p_at(mid) < alpha:
            lo = mid
        else:
            hi = mid
    return float(lo)


def weighted_trend_test(
    households: pd.DataFrame,
    member_ids: set,
    group_col: str = "wealth_quintile",
    weight_col: str = "wgt_adj_consagg",
    deff: float = None,
) -> dict:
    """Variance-weighted linear trend across ordered groups (e.g. quintiles)."""
    if deff is None:
        deff = deff_from_icc(PRIMARY_ICC)
    sub = households[households[group_col].notna()]
    xs, ys, ses = [], [], []
    for g, grp in sub.groupby(group_col):
        p, se = weighted_proportion(grp, member_ids, weight_col)
        xs.append(float(g)); ys.append(p); ses.append(se * np.sqrt(deff))
    x, y, se = np.array(xs), np.array(ys), np.array(ses)
    w = 1.0 / se ** 2
    xb = np.average(x, weights=w)
    yb = np.average(y, weights=w)
    slope = (w * (x - xb) * (y - yb)).sum() / (w * (x - xb) ** 2).sum()
    se_slope = np.sqrt(1.0 / (w * (x - xb) ** 2).sum())
    z = slope / se_slope
    return {
        "slope_pp_per_group": 100 * float(slope),
        "z": float(z),
        "p_value": float(2 * stats.norm.sf(abs(z))),
        "deff": float(deff),
    }

# ======================================================================
#  RUN THE ANALYSIS
# ======================================================================
if __name__ == "__main__":
    import sys

    # ---- Colab: prompt for the prepared data file -------------------
    try:
        from google.colab import files          # type: ignore
        print("Select your prepared SIHBS data file...")
        uploaded = files.upload()
        DATA_PATH = list(uploaded.keys())[0]
    except Exception:
        # Running locally: pass the path as the first argument
        DATA_PATH = sys.argv[1] if len(sys.argv) > 1 else "prepared_sihbs_data.xlsx"
    print("Using:", DATA_PATH, "\n")

    d = load_sihbs(DATA_PATH)
    hh, food, items, livestock = d["households"], d["food"], d["items"], d["livestock"]
    agg = build_household_item_table(food, hh)

    print("SAMPLE")
    print(f"  Households          {hh['hid'].nunique():>10,}   [paper: 7,212]")
    print(f"  Consumption records {len(food):>10,}   [paper: 152,999]")
    print(f"  Population          {hh['wgt_pop'].sum():>10,.0f}   [paper: 15,361,943]")
    print(f"  Welfare aggregate   {hh['wealth_quintile'].notna().sum():>10,}   [paper: 7,152]")

    print("\nIMPLAUSIBILITY SCREEN")
    scr = screening_summary(agg)
    print(f"  Flagged {int(agg['implausible'].sum())} records over {len(scr)} items "
          f"in {agg.loc[agg['implausible'],'hid'].nunique()} households   [paper: 188 / 12 / 184]")

    print("\nHOUSEHOLD COVERAGE (Table 1)")
    for it in HEADLINE_ITEMS:
        r = coverage(agg, hh, it)
        print(f"  {it:24s} {r['coverage_pct']:5.1f}%  "
              f"({r['coverage_pct']-1.96*r['se_pct']:.1f}-{r['coverage_pct']+1.96*r['se_pct']:.1f})"
              f"   n={r['n_consuming']:,}")

    print("\nPER-CAPITA QUANTITY AMONG CONSUMERS (Table 5)")
    for it in ["Sugar - white", "Wheat flour - white", "Rice - White [Local]",
               "Maize grain", "Sorghum grain", "Onion"]:
        r = among_consumers_quantity(agg, hh, it)
        print(f"  {it:24s} mean {r['mean']:6.1f}  median {r['median']:6.1f}  P95 {r['p95']:6.1f}")

    print("\nPOPULATION-LEVEL QUANTITY (Table 9, Equation 2)")
    cereal = ["Wheat flour - white", "Rice - White [Local]", "Maize grain", "Sorghum grain"]
    for it in ["Sugar - white"] + cereal + ["Onion"]:
        print(f"  {it:24s} {population_quantity(agg, hh, it):6.2f} g/cap/day")
    print(f"  {'Cereal total (4 items)':24s} "
          f"{sum(population_quantity(agg, hh, i) for i in cereal):6.1f} g/person/day   [paper: 140.3]")

    print("\nLIVESTOCK LINKAGE (Table 3)")
    owners  = set(livestock.loc[livestock["animalid"] == "CAMELS", "hid"])
    cm      = set(agg.loc[agg["item_name"] == "Camel meat - fresh", "hid"])
    ck      = set(agg.loc[agg["item_name"] == "Camel Milk - Fresh", "hid"])
    anymeat = set(agg.loc[agg["item_name"].isin(MEAT_ITEMS), "hid"])
    for stratum, g in hh.groupby("ea_type_n"):
        print(f"  {stratum:8s} own {100*weighted_proportion(g, owners)[0]:5.1f}%  "
              f"meat {100*weighted_proportion(g, cm)[0]:5.1f}%  "
              f"milk {100*weighted_proportion(g, ck)[0]:5.1f}%  "
              f"any meat {100*weighted_proportion(g, anymeat)[0]:5.1f}%")
    nom = hh[hh["ea_type_n"] == "Nomadic"]
    o, n = nom[nom["hid"].isin(owners)], nom[~nom["hid"].isin(owners)]
    print(f"  Nomadic camel owners     milk {100*weighted_proportion(o, ck)[0]:5.1f}%  "
          f"meat {100*weighted_proportion(o, cm)[0]:5.1f}%")
    print(f"  Nomadic non-owners       milk {100*weighted_proportion(n, ck)[0]:5.1f}%  "
          f"meat {100*weighted_proportion(n, cm)[0]:5.1f}%")

    print(f"\nFORMAL TESTS, DEFF = {deff_from_icc(PRIMARY_ICC):.2f} (Table 4)")
    urban = hh[hh["ea_type_n"] == "Urban"]
    hq = hh[hh["wealth_quintile"].notna()]
    q1, q5 = hq[hq["wealth_quintile"] == 1], hq[hq["wealth_quintile"] == 5]
    mz = set(agg.loc[agg["item_name"] == "Maize grain", "hid"])
    for label, g1, g2, ids, wc in [
        ("Camel meat: urban vs nomadic",        urban, nom, cm, "wgt"),
        ("Camel milk: nomadic owners vs non",   o,     n,   ck, "wgt"),
        ("Maize: Q1 vs Q5",                     q1,    q5,  mz, "wgt_adj_consagg"),
        ("Camel milk: urban vs nomadic",        urban, nom, ck, "wgt"),
        ("Camel meat: nomadic owners vs non",   o,     n,   cm, "wgt"),
    ]:
        r  = deff_adjusted_test(g1, g2, ids, weight_col=wc)
        bd = breakdown_deff(g1, g2, ids, weight_col=wc)
        bd_s = ">2000" if bd >= 2000 else ("n/a" if np.isnan(bd) else f"{bd:.1f}")
        print(f"  {label:36s} {r['diff_pp']:+6.1f} pp "
              f"({r['ci_low_pp']:+.1f},{r['ci_high_pp']:+.1f})  "
              f"z={r['z']:6.2f}  P={r['p_value']:.3g}  breakdownDEFF={bd_s}")
    tr = weighted_trend_test(hh, mz)
    print(f"  {'Maize wealth trend':36s} {tr['slope_pp_per_group']:+6.2f} pp/quintile"
          f"                z={tr['z']:6.2f}  P={tr['p_value']:.3g}")

    print("\nFULL 164-ITEM REFERENCE (Supplementary Table S1)")
    s1 = []
    for _, row in items.iterrows():
        name = row["item_name"]
        sub  = agg[agg["item"] == row["item_code"]]
        if len(sub):
            c = coverage(agg, hh, name)
            liquid = sub["mL"].sum() > sub["g"].sum()
            q = population_quantity(agg, hh, name, liquid=liquid)
            unit = "mL/cap/day" if liquid else "g/cap/day"
        else:
            c, q, unit = {"n_consuming": 0, "coverage_pct": 0.0}, np.nan, "n/a"
        s1.append({"item_code": int(row["item_code"]), "item_name": name,
                   "n_consuming": c["n_consuming"],
                   "coverage_pct": round(c["coverage_pct"], 1),
                   "population_quantity": None if pd.isna(q) else round(q, 2),
                   "unit": unit,
                   "n_excluded": int(sub["implausible"].sum()) if len(sub) else 0,
                   "headline": "Yes" if name in HEADLINE_ITEMS else ""})
    s1df = (pd.DataFrame(s1)
            .sort_values(["coverage_pct", "item_name"], ascending=[False, True])
            .reset_index(drop=True))
    s1df.insert(0, "rank", range(1, len(s1df) + 1))
    s1df.to_csv("table_S1_all_164_items.csv", index=False)
    print(f"  {len(s1df)} items written to table_S1_all_164_items.csv   [paper: 164]")

    print("\nDone.")
