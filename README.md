# Derived outputs

This folder is intentionally empty in the published repository.

Running the reproduction script populates it with the derived tables:

```bash
python src/reproduce_paper.py --data /path/to/prepared_sihbs_data.xlsx
```

produces:

| File | Corresponds to |
|---|---|
| `table_1_coverage_headline.csv` | Table 1, household coverage of the headline items |
| `table_5_percapita_solids.csv` | Table 5, per-capita quantities among consumers |
| `table_9_population_quantities.csv` | Table 9, person-weighted population quantities |
| `table_4_formal_tests.csv` | Table 4, formal tests and design effect sensitivity |
| `table_S1_all_164_items.csv` | Supplementary Table S1, full 164-item reference |
| `table_S2_screening_by_item.csv` | Records removed per item by the implausibility screen |

The tables are regenerated from the microdata rather than shipped, so that the
repository contains code only and nothing derived from data that is not the
authors' to redistribute. The published versions of these tables are available
with the paper and its Supplementary Information.
