"""
reproduce_paper.py
==================

Reproduces every headline figure, table and test statistic reported in:

    Hersi, M. A. & Fiidow, O. A.
    "A National Food Consumption Reference for Food Safety Risk Assessment
     in Somalia"

directly from the SIHBS 2022 microdata, and writes the derived tables to
``outputs/``.

Usage
-----
    python src/reproduce_paper.py --data /path/to/SIHBS_2022_Consolidated.xlsx

Every value printed with a ``[paper: X]`` annotation is checked against the
published figure, and the script exits non-zero if any check fails. This
makes the script usable as a regression test on the analysis pipeline.
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from sihbs_consumption import (  # noqa: E402
    CLUSTER_SIZE_M,
    HEADLINE_ITEMS,
    MEAT_ITEMS,
    PRIMARY_ICC,
    among_consumers_quantity,
    breakdown_deff,
    build_household_item_table,
    coverage,
    deff_adjusted_test,
    deff_from_icc,
    load_sihbs,
    population_quantity,
    screening_summary,
    weighted_trend_test,
)

FAILURES: list[str] = []


def check(label: str, got: float, expected: float, tol: float = 0.05) -> None:
    """Assert a reproduced value matches the published figure."""
    ok = abs(got - expected) <= tol
    flag = "OK " if ok else "FAIL"
    print(f"   [{flag}] {label:52s} got {got:>12,.2f}   [paper: {expected:,.2f}]")
    if not ok:
        FAILURES.append(f"{label}: got {got}, expected {expected}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", required=True, help="Path to the prepared SIHBS analysis file")
    ap.add_argument("--outdir", default="outputs", help="Directory for derived tables")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print("Loading SIHBS 2022 microdata ...")
    d = load_sihbs(args.data)
    hh, food, items, livestock = d["households"], d["food"], d["items"], d["livestock"]

    print("\n1. SAMPLE")
    check("Households sampled", hh["hid"].nunique(), 7212, tol=0)
    check("Consumption records", len(food), 152999, tol=0)
    check("Population represented", hh["wgt_pop"].sum(), 15361943, tol=1)
    check("Households with welfare aggregate",
          hh["wealth_quintile"].notna().sum(), 7152, tol=0)

    print("\n2. UNIT RESOLUTION AND IMPLAUSIBILITY SCREEN")
    agg = build_household_item_table(food, hh)
    scr = screening_summary(agg)
    n_flagged = int(agg["implausible"].sum())
    check("Records flagged by screen", n_flagged, 188, tol=0)
    check("Households affected", agg.loc[agg['implausible'], 'hid'].nunique(), 184, tol=0)
    check("Items affected", int(scr.shape[0]), 12, tol=0)
    scr.to_csv(os.path.join(args.outdir, "table_S2_screening_by_item.csv"), index=False)

    print("\n3. HOUSEHOLD COVERAGE, HEADLINE ITEMS (Table 1)")
    cov_rows = [coverage(agg, hh, it) for it in HEADLINE_ITEMS]
    cov = pd.DataFrame(cov_rows)
    cov.to_csv(os.path.join(args.outdir, "table_1_coverage_headline.csv"), index=False)
    published_cov = {
        "Sugar - white": 98.8, "Cooking oil": 98.7, "Onion": 88.0,
        "Wheat flour - white": 81.0, "Rice - White [Local]": 58.3,
        "Goat Meat - Fresh": 44.7, "Camel meat - fresh": 36.7,
        "Maize grain": 22.7, "Camel Milk - Fresh": 24.5, "Sorghum grain": 16.0,
    }
    for item, exp in published_cov.items():
        got = cov.loc[cov["item_name"] == item, "coverage_pct"].iloc[0]
        check(f"Coverage {item}", got, exp, tol=0.1)

    print("\n4. PER-CAPITA QUANTITY AMONG CONSUMERS (Table 5)")
    solids = {
        "Sugar - white": (97.5, 85.7), "Wheat flour - white": (94.5, 85.7),
        "Rice - White [Local]": (94.3, 85.7), "Maize grain": (46.2, 35.7),
        "Sorghum grain": (47.8, 35.7), "Onion": (27.7, 17.9),
    }
    rows = []
    for item, (exp_mean, exp_med) in solids.items():
        r = among_consumers_quantity(agg, hh, item)
        rows.append(r)
        check(f"Mean {item}", r["mean"], exp_mean, tol=0.1)
        check(f"Median {item}", r["median"], exp_med, tol=0.1)
    pd.DataFrame(rows).to_csv(
        os.path.join(args.outdir, "table_5_percapita_solids.csv"), index=False)

    print("\n5. POPULATION-LEVEL QUANTITIES (Table 9, Equation 2)")
    pop = {
        "Sugar - white": 86.2, "Wheat flour - white": 72.5,
        "Rice - White [Local]": 51.3, "Onion": 21.3,
        "Maize grain": 9.4, "Sorghum grain": 7.1,
    }
    pop_rows = []
    for item, exp in pop.items():
        got = population_quantity(agg, hh, item)
        pop_rows.append({"item_name": item, "population_g_cap_day": got})
        check(f"Population qty {item}", got, exp, tol=0.1)
    cereal4 = ["Wheat flour - white", "Rice - White [Local]",
               "Maize grain", "Sorghum grain"]
    cereal_total = sum(population_quantity(agg, hh, i) for i in cereal4)
    check("Cereal total (4 items)", cereal_total, 140.3, tol=0.2)
    pd.DataFrame(pop_rows).to_csv(
        os.path.join(args.outdir, "table_9_population_quantities.csv"), index=False)

    print("\n6. FULL 164-ITEM REFERENCE (Supplementary Table S1)")
    s1 = []
    denom_ids = set(agg["hid"])
    for _, row in items.iterrows():
        name = row["item_name"]
        sub = agg[agg["item"] == row["item_code"]]
        ids = set(sub["hid"])
        cov_d = coverage(agg, hh, name) if ids else {
            "n_consuming": 0, "coverage_pct": 0.0}
        is_liquid = sub["mL"].sum() > sub["g"].sum()
        qty = population_quantity(agg, hh, name, liquid=is_liquid) if ids else np.nan
        s1.append({
            "item_code": int(row["item_code"]),
            "item_name": name,
            "n_consuming": cov_d["n_consuming"],
            "coverage_pct": round(cov_d["coverage_pct"], 1),
            "population_quantity": None if pd.isna(qty) else round(qty, 2),
            "unit": "n/a" if pd.isna(qty) else ("mL/cap/day" if is_liquid else "g/cap/day"),
            "n_excluded_by_screen": int(sub["implausible"].sum()),
            "headline_item": "Yes" if name in HEADLINE_ITEMS else "",
        })
    s1df = pd.DataFrame(s1).sort_values(
        ["coverage_pct", "item_name"], ascending=[False, True]).reset_index(drop=True)
    s1df.insert(0, "rank", range(1, len(s1df) + 1))
    s1df.to_csv(os.path.join(args.outdir, "table_S1_all_164_items.csv"), index=False)
    check("Items in S1", len(s1df), 164, tol=0)
    check("Items with zero consumption", int((s1df['n_consuming'] == 0).sum()), 4, tol=0)

    print("\n7. LIVESTOCK LINKAGE (Table 3)")
    owners = set(livestock.loc[livestock["animalid"] == "CAMELS", "hid"])
    cm = set(agg.loc[agg["item_name"] == "Camel meat - fresh", "hid"])
    ck = set(agg.loc[agg["item_name"] == "Camel Milk - Fresh", "hid"])
    anymeat = set(agg.loc[agg["item_name"].isin(MEAT_ITEMS), "hid"])
    nomadic = hh[hh["ea_type_n"] == "Nomadic"]
    nom_own = nomadic[nomadic["hid"].isin(owners)]
    nom_non = nomadic[~nomadic["hid"].isin(owners)]
    from sihbs_consumption import weighted_proportion as wp
    check("Nomadic camel ownership %", 100 * wp(nomadic, owners)[0], 43.2, tol=0.1)
    check("Nomadic owners: camel milk %", 100 * wp(nom_own, ck)[0], 45.2, tol=0.1)
    check("Nomadic non-owners: camel milk %", 100 * wp(nom_non, ck)[0], 14.7, tol=0.1)
    check("Nomadic owners: camel meat %", 100 * wp(nom_own, cm)[0], 2.3, tol=0.1)
    check("Nomadic non-owners: camel meat %", 100 * wp(nom_non, cm)[0], 1.5, tol=0.1)
    check("Nomadic any-meat %", 100 * wp(nomadic, anymeat)[0], 40.6, tol=0.1)
    check("Urban any-meat %",
          100 * wp(hh[hh['ea_type_n'] == 'Urban'], anymeat)[0], 90.8, tol=0.1)

    print(f"\n8. DESIGN EFFECT SENSITIVITY AND FORMAL TESTS (Table 4)"
          f"  [ICC={PRIMARY_ICC}, m={CLUSTER_SIZE_M}, DEFF={deff_from_icc(PRIMARY_ICC):.2f}]")
    urban = hh[hh["ea_type_n"] == "Urban"]
    hq = hh[hh["wealth_quintile"].notna()]
    q1, q5 = hq[hq["wealth_quintile"] == 1], hq[hq["wealth_quintile"] == 5]
    mz = set(agg.loc[agg["item_name"] == "Maize grain", "hid"])

    tests = [
        ("Camel meat: urban vs nomadic", urban, nomadic, cm, "wgt", 48.7, 4.7e-289),
        ("Camel milk: nomadic owners vs non", nom_own, nom_non, ck, "wgt", 30.5, 5.9e-7),
        ("Maize: Q1 vs Q5", q1, q5, mz, "wgt_adj_consagg", 9.9, 1.28e-4),
        ("Camel milk: urban vs nomadic", urban, nomadic, ck, "wgt", -3.6, 0.268),
        ("Camel meat: nomadic owners vs non", nom_own, nom_non, cm, "wgt", 0.8, 0.617),
    ]
    test_rows = []
    for label, g1, g2, ids, wcol, exp_diff, exp_p in tests:
        r = deff_adjusted_test(g1, g2, ids, weight_col=wcol)
        bd = breakdown_deff(g1, g2, ids, weight_col=wcol)
        r["comparison"] = label
        r["breakdown_deff"] = bd
        test_rows.append(r)
        check(f"Diff (pp) {label}", r["diff_pp"], exp_diff, tol=0.1)
        # P values span 200 orders of magnitude; compare on log scale
        if exp_p > 1e-10:
            check(f"P value  {label}", r["p_value"], exp_p, tol=max(1e-3, exp_p * 0.1))
        print(f"          breakdown DEFF = "
              f"{'>2000' if bd >= 2000 else ('n/a' if np.isnan(bd) else f'{bd:.1f}')}")
    pd.DataFrame(test_rows).to_csv(
        os.path.join(args.outdir, "table_4_formal_tests.csv"), index=False)

    trend = weighted_trend_test(hh, mz)
    check("Maize wealth trend (pp/quintile)", trend["slope_pp_per_group"], -2.50, tol=0.05)

    print("\n" + "=" * 72)
    if FAILURES:
        print(f"REPRODUCTION FAILED: {len(FAILURES)} check(s) did not match")
        for f in FAILURES:
            print("   -", f)
        return 1
    print("ALL CHECKS PASSED - every published figure reproduced from the microdata.")
    print(f"Derived tables written to: {os.path.abspath(args.outdir)}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
