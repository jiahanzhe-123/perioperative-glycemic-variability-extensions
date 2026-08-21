# Phase 3B public builder

`01_public_aggregate_landscape.R` rebuilds the public aggregate version of the
analytic-context landscape from the 18-row SD-GV-only table in
`results/phase3b/figure_3_source_data.csv`.

The public builder reads aggregate rows only. It does not read patient/stay
records, same-patient glucose values, real-data bootstrap replicates, or local
configuration, and it does not refit any model. Figure 2's patient-level
scatter and Bland–Altman source data remain in the controlled local production
package and are intentionally not part of the public release tier.

Run from the repository root:

```bash
Rscript analysis_extensions/scripts/phase3b/01_public_aggregate_landscape.R
```

The builder writes the public aggregate figure to
`analysis_extensions/figures/phase3b_public/`.
