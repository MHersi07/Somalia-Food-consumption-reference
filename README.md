# A National Food Consumption Reference for Food Safety Risk Assessment in Somalia

Analysis code accompanying:
> Hersi, M. A. & Fiidow, O. A. *A National Food Consumption Reference for Food
> Safety Risk Assessment in Somalia.*

Running `src/reproduce_paper.py` regenerates **every figure, table and test
statistic reported in the manuscript** and verifies each against the published
value, exiting non-zero if any check fails.

---

## Contents

```
src/sihbs_consumption.py    Analysis pipeline
src/reproduce_paper.py      End-to-end reproduction and verification
METHODS_NOTES.md            Detail on the estimators
CITATION.cff                How to cite this code
outputs/                    Empty; populated by running the script
requirements.txt            Dependencies
LICENSE                     MIT (code only)
```

---

## 1. Obtaining the microdata

The SIHBS 2022 microdata are the property of the **Somali National Bureau of
Statistics (SNBS)** and are available from the SNBS microdata portal
(<https://microdata.nbs.gov.so>), **subject to registration and approval by
SNBS**. They are not included here and we are not permitted to redistribute
them; access must be requested directly from SNBS.

Once approved, you will receive the SIHBS 2022 files in Stata `.dta` format.

## 2. Files you need from the release

| SIHBS module                                                | Provides                                                                          |
| ------------------------------------------------------------ | ----------------------------------------------------------------------------------|
| Food consumption (7-day recall, long format; ~153,000 rows) | Household identifier, food item code and label, quantity consumed, unit code      |
| Household roster / cover                                    | Household size, region, residence type (urban / rural / nomadic), sampling weight |
| Consumption aggregate                                       | Real per-capita expenditure decile and its adjusted weight                        |
| Livestock ownership                                         | Household identifier and animal type                                              |

Variable names as released:

- **Food module:** `cr15_04quantity` (quantity consumed in the last 7 days),
  `cr15_04units` (unit code 1–14), plus the item code and label.
- **Household:** `hhsize`, `wgt`, `wgt_pop` (= `wgt` × `hhsize`), `region_n`,
  `ea_type_n`.
- **Consumption aggregate:** `pcer_dec` (expenditure decile),
  `wgt_adj_consagg`.
- **Livestock:** animal type identifier (camels are needed for the ownership
  linkage).

## 3. Preparing the data

**a. Join the modules.** Merge the food consumption file, household roster and
consumption aggregate on the household identifier. Collapse the expenditure
decile `pcer_dec` into quintiles to create `wealth_quintile`. Keep livestock as
a separate table.

**b. Resolve the metric unit codes.** Quantities are recorded against 14 unit
codes (`cr15_04units`). Codes 1–4 are metric and convert directly:

| Code | Unit        | Conversion           |
| ---- | ----------- | --------------------- |
| 1    | Grams       | as reported           |
| 2    | Kilograms   | × 1000 → grams        |
| 3    | Litres      | × 1000 → millilitres  |
| 4    | Millilitres | as reported           |

These cover 75.7% of all consumption records. Among the headline items: sugar,
rice, goat meat and camel meat are 100% metric; wheat flour and cooking oil
99.9%; camel milk 97%; maize and sorghum 75–78%; onion 59% kg + 21% piece +
19% gram.

**c. Convert the non-standard units.** Codes 5–14 (cans, tins, bottles,
bunches, plates, pieces) are converted using the official national conversion
factors published by SNBS in the *Somalia Poverty Report 2023*, Appendix 2.1,
Table B1, derived from a dedicated Unit of Measurement survey conducted in
October 2021. Obtain that table from the source document; it is not reproduced
here.

Two properties of the source table determine how it should be applied:

- **Solids only.** Table B1's factors for liquid items are denominated in
  grams rather than millilitres, implying a density the source does not
  state. Apply the table to solid items only — which achieves complete
  resolution for solids — and report liquid items from metric-unit records
  alone (97–99.9% of records), leaving the small non-metric remainder
  unconverted.
- **No meat entries.** Goat meat and camel meat have no entry, consistent
  with the empirical finding that they are reported exclusively in metric
  units.

Write the results into two columns, `resolved_grams` and `resolved_mL`,
leaving the other blank for each record.

**d. Assemble.** Supply four tables to `load_sihbs()`. The default reader
expects them as four sheets in a single file (`Household_Master`,
`Food_Consumption_Long`, `Food_Item_Dictionary`, `Livestock_Own`). If you hold
them as separate `.dta` or `.csv` files, replace `load_sihbs()` with your own
loader returning the same four DataFrames — nothing downstream depends on the
file format.

Required columns after preparation:

- **Households:** `hid`, `hhsize`, `wgt`, `wgt_pop`, `wgt_adj_consagg`,
  `region_n`, `ea_type_n`, `wealth_quintile`
- **Food (long):** `hid`, `item`, `item_name`, `cr15_04quantity`,
  `resolved_grams`, `resolved_mL`
- **Item dictionary:** `item_code`, `item_name`
- **Livestock:** `hid`, `animalid`

## 4. Running

Requires Python 3.9 or later.

```
pip install -r requirements.txt
python src/reproduce_paper.py --data /path/to/prepared_sihbs_data.xlsx
```

Expected output ends with:

```
ALL CHECKS PASSED - every published figure reproduced from the microdata.
```

The script verifies 47 published quantities: sample counts, household
coverage for all 11 headline items, per-capita quantities among consumers,
person-weighted population-level quantities, the full 164-item reference, the
livestock-linkage results, and the design effect sensitivity analysis.

---

## Method summary

| Step                                        | Function                                                       |
| -------------------------------------------- | ---------------------------------------------------------------|
| Unit resolution, implausibility screening   | `build_household_item_table`, `screening_summary`              |
| Household coverage                          | `coverage`                                                     |
| Per-capita quantity among consumers         | `among_consumers_quantity`                                     |
| Person-weighted population quantity         | `population_quantity`                                          |
| Survey weighting, stratified variance       | `weighted_proportion`                                          |
| Design effect sensitivity, hypothesis tests | `deff_adjusted_test`, `breakdown_deff`, `weighted_trend_test` |

**Implausibility screen.** Records exceeding 1,000 g/cap/day for a solid item
or 5,000 mL/cap/day for a liquid are excluded from **quantity** estimates but
retained in **coverage**, because the household did consume the item and only
the reported quantity is erroneous.

**Variance.** Confidence intervals use Taylor linearisation of the Hájek
ratio estimator under stratified sampling, with region × residence as the
design strata (49 non-empty strata, minimum size 48). Clustering is not
incorporated because the released files omit the EA/PSU identifier; that
omission is parameterised through the design effect. See `METHODS_NOTES.md`.

---

## Citation

If you use this code, please cite it via `CITATION.cff` (GitHub's "Cite this
repository" button generates the formatted citation).

---

## Licence

The **code** is released under the MIT Licence (see `LICENSE`). Nothing else
is covered by it. This repository does not redistribute the SIHBS microdata,
the SIHBS unit codebook, the non-standard-unit conversion factors, or derived
result tables; all must be obtained from their original sources.
