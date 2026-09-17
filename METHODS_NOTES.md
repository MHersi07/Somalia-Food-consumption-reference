# Methods notes

Supplementary documentation for the analysis pipeline. This file records the
source-of-truth values that the code depends on, so that the derivation can be
reconstructed independently.

## 1. Unit codebook

The SIHBS consumption module records quantities against 14 unit codes, captured
in questionnaire question `cr15_04units` ("Q15.04A2 UNITS"). Codes 1-4 are
metric (grams, kilograms, litres, millilitres) and are directly usable; codes
5-14 are non-standard units (tins, bottles, cans, bunches, plates, pieces).

**The code list itself is not reproduced here.** It is the property of the
Somali National Bureau of Statistics and is published in the SIHBS
questionnaire, available from the SNBS microdata portal:

- SNBS. *Questionnaire: Somali Integrated Household Budget Survey (SIHBS)*
  (2022). Somali National Bureau of Statistics.
  <https://microdata.nbs.gov.so> (registration and approval required)

Coverage of metric codes, computed from the microdata: codes 1-4 account for
75.7% of all 152,999 consumption records. For the headline items: sugar, rice,
goat meat and camel meat are 100% metric; wheat flour and cooking oil 99.9%;
camel milk 97%; maize and sorghum 75-78%; onion 59% kg + 21% piece + 19% gram.

## 2. Non-standard-unit (NSU) conversion

Non-metric records are resolved using the official national conversion factors
prepared by SNBS from a dedicated Unit of Measurement survey conducted in
October 2021, published as:

- SNBS & World Bank. *Somalia Poverty Report 2023*, Appendix 2.1, Table B1
  ("Non-standard unit conversion factors (NSU to grams) for consumed items").

**These conversion factors are not reproduced here.** They are the property of
SNBS and the World Bank; users should obtain them from the source document
above and supply them to the pipeline via the `resolved_grams` / `resolved_mL`
columns (see README, Data availability).

Two properties of the source table that affect the analysis, and which are
documented in the manuscript:

- **Liquid density caveat.** Table B1's factors for liquid items are
  denominated in grams, not millilitres, implying an assumed density the source
  does not state. The analysis therefore applies the NSU table only to solid
  items (100% resolution) and reports liquid items using metric-unit records
  only (97-99.9% of records), leaving the small non-metric remainder
  unconverted.
- **No meat entries.** Goat meat and camel meat have no NSU entry in either
  Poverty Report table, confirming the empirical finding that these items are
  reported exclusively in metric units by design.

## 3. Weight variables

| Variable | Used for |
|---|---|
| `wgt` | Household-level statistics: coverage, per-capita quantity distributions |
| `wgt_pop` (= `wgt` × `hhsize`) | National population total (15,361,943) |
| `wgt_adj_consagg` | Wealth-quintile-stratified estimates (7,152 households with a valid welfare aggregate) |

**Documented pitfall.** Merging `wgt_pop` onto the person-level roster and
summing per person double-applies household size, yielding 119.9 million
instead of the correct 15,361,943.

## 4. Implausibility screen

Pre-specified, absolute, item-independent bounds chosen to remove only values
that are physiologically impossible rather than merely high:

- **Solids:** > 1,000 g/cap/day for a single item. Total daily food intake
  across all items is of the order of 1–2 kg, so a single item exceeding 1 kg
  per person per day cannot be correct. The bound lies well above the 99.9th
  percentile of the observed distribution (686 g/cap/day).
- **Liquids:** > 5,000 mL/cap/day for a single item, reflecting the upper limit
  of physiological fluid intake. Deliberately conservative: every liquid item
  other than mineral water has an observed maximum below 1,500 mL/cap/day.

Result: 188 of 152,999 household-item records (0.123%), across 184 households
and 12 items — mineral water (127), pumpkin (38), spinach (13), green pepper
(2), and one record each for tomato, cow peas, mango, onion, white potatoes,
cabbages, spaghetti and oranges (8 records across eight items).
127 + 38 + 13 + 2 + 8 = 188.

Flagged records are excluded from **quantity** estimates only and retained in
all **coverage** estimates.

## 5. Design effect sensitivity

The public-use SIHBS files omit the EA/PSU identifier, so clustering variance
cannot be estimated directly. The Kish design effect is used instead:

    DEFF = 1 + (m − 1) · ICC,    m = 12 households per EA (601 EAs)

with ICC 0.05–0.20 adopted as a plausible, deliberately conservative range for
household food consumption outcomes, giving DEFF 1.55–3.20 and confidence
intervals 1.24–1.79 times wider than the stratified-only intervals. ICC = 0.10
(DEFF = 2.10) is the primary specification.

Empirical support for this range:

- Geyer, J., Davis, M. & Narayan, T. Intracluster correlation coefficients of
  household economic and agricultural outcomes in Mozambique.
  *Evaluation Review* **40**, 526–545 (2016).
- Seidenfeld, D., Handa, S., de Hoop, T. & Morey, M. Intraclass correlations
  values in international development: evidence across commonly studied
  domains in sub-Saharan Africa. *Evaluation Review* **47**, 786–819 (2023).

For each substantively interpreted comparison the pipeline also reports the
**breakdown DEFF**: the design effect at which the comparison would cease to
reach p < 0.05. A breakdown DEFF far above the plausible range indicates the
comparison is robust to the absent clustering term.

Note that with m = 12 and ICC bounded above by 1.0, the maximum attainable
DEFF under equal cluster sizes is 12.0. Breakdown values above that are
therefore unattainable in principle, not merely implausible. (With unequal
cluster sizes the effective m exceeds 12, so this ceiling is not absolute.)
