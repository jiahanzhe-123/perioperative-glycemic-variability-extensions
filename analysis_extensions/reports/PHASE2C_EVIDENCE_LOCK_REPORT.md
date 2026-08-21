# PHASE2C Evidence Lock Report

Phase 2C assembled a provenance-bearing evidence package from the final Phase 1.6 lineage and locked Phase 2A/2B outputs. No new inferential analysis was performed in this phase.

## Material Passport

- Workspace: `JAHA_v5_analysis_extensions`.
- Upstream gates: `GO_PHASE2_MEASUREMENT_CONTEXT_ANALYSES`, `GO_PHASE2B_ANALYTIC_CONTEXT_LANDSCAPE`, and the Phase 2B `GO_PHASE3_MANUSCRIPT_REWRITE` gate were required before assembly.
- Primary registry: `outputs/evidence_lock/master_evidence_registry.csv`.
- Locked landscape: `outputs/evidence_lock/analytic_context_landscape_LOCKED.csv`.
- Figure source lock: `outputs/evidence_lock/figure_source_map_LOCKED.csv`.
- Reproducibility manifest: `outputs/evidence_lock/evidence_lock_manifest.json`.
- The submitted manuscript and controlled inputs were read-only provenance sources; they were not modified.

## 1. Can every proposed main numerical claim be traced to a locked result or an extension result?

Yes, subject to the Phase 2C QC gate. The master registry assigns a unique evidence ID, population, estimand/effect type, numeric fields, source file, source row/key, script origin, lineage status, and manuscript location to every proposed main numerical claim. The primary traceability chain is: MIMIC Model A/B/C and final context summary; same-patient Phase 2A agreement and source-model reproduction; paired bootstrap; Phase 2B log-scale and count-imbalance outputs; and the 18-row Phase 2B analytic-context landscape. The original mean-glucose/absolute-risk result is also registered but explicitly marked supplementary.

## 2. Are the proposed main claims reproducible from frozen scripts and inputs?

Yes for the locked evidence package. The manifest records SHA-256 hashes for the Phase 2A/2B machine outputs, candidate figures, scripts, final-lineage inputs, submitted provenance files, registries, maps, and rewrite brief. Reproduction here means re-reading and verifying the locked artifacts; Phase 2C does not rerun models. Candidate figure bitmap regeneration remains a Phase 3 production task because the submitted figure builders are not present in the public-backup repository.

## 3. What are the core 5–10 results that deserve main-manuscript emphasis?

1. Final context identity and estimands: MIMIC primary N=10,561/296; INSPIRE I2 N=1,353/27 and corrected 48-hour N=1,511/31; eICU M3/M4 N=7,115/130. Source IDs: `MIMIC_PRIMARY_COHORT_30D`, `INSPIRE_I2_24H_30D`, `INSPIRE_CORRECTED_48H_30D`, `EICU_M3`, `EICU_M4`.
2. Same-patient GV agreement on the original scale: POCT mean 28.6203, laboratory mean 34.6317, Pearson r 0.4893, Spearman rho 0.3269; mean difference -6.0114 mg/dL with limits -69.9302 to 57.9075 mg/dL. Evidence IDs: `P2A_AGREEMENT_PRIMARY_GV_MEAN_DIFFERENCE`, `P2A_AGREEMENT_PRIMARY_GV_PEARSON`, `P2A_AGREEMENT_PRIMARY_GV_SPEARMAN`.
3. Scale dependence of disagreement: primary positive-pair effective N=447, mean log ratio 0.0303, median -0.0012, geometric mean ratio 1.0308; exponentiated limits 0.1060 to 10.0207; log-ratio versus log-pair-mean Pearson/Spearman -0.3198/-0.3324. Evidence IDs: `P2B_LOG_PRIMARY_*`; this sensitivity does not replace the original scale.
4. Identical-patient mortality coefficients: POCT HR 0.7589 (N=409/49) versus laboratory HR 0.9709 (N=409/49). Evidence IDs: `MIMIC_SAMEPATIENT_POCT_HR` and `MIMIC_SAMEPATIENT_LAB_HR`; the identical N/events comparison is marked in the candidate figure source lock.
5. Paired coefficient-difference sensitivity: delta-beta -0.2464 (95% percentile interval -0.6323 to -0.0551), exponentiated ratio 0.7816, bootstrap median delta-beta -0.2687; B=2,000, successful=2,000, failed=0. Evidence IDs: `P2A_BOOTSTRAP_OBSERVED_DELTA_BETA` and `P2A_BOOTSTRAP_COEFFICIENT_RATIO`; source-definition descriptive evidence, not a causal interaction.
6. MIMIC adjustment context: `MIMIC_MODEL_A_30D`, `MIMIC_MODEL_B_30D`, and `MIMIC_MODEL_C_30D` in conceptual order; the final primary Model B is HR 0.9794 (0.9225–1.0398), P=0.4952, with Model C as the locked context shift.
7. Source-specific MIMIC context: priority series, POCT-only, central-laboratory-only, blood-gas-only, common-source, and the same-patient pair in the `MIMIC_SOURCE_*` rows. The count-imbalance closure shows delta-count mean 2.2279 (SD 2.7576), median 2, IQR 4, q5–q95 -1–7; source-defined GV disagreement Pearson/Spearman -0.0204/0.0556. Evidence ID: `P2B_COUNT_PRIMARY_GV_CORRELATION`; source-specific span is `P2B_SPAN_PRIMARY_NOT_ESTIMABLE`.
8. Analytic-context landscape and context shifts: all `LANDSCAPE_*` rows are separated into MIMIC adjustment, MIMIC source, INSPIRE timing/adjustment, and eICU adjustment families with HR/RR preserved. Locked shifts are MIMIC B→C delta log effect 0.0490, eICU M3→M4 -0.0218, and INSPIRE I2→I3 0.0261; no pooling or heterogeneity test.

## 4. Which analyses belong in the supplement?

The supplement should contain the N=453 paired sensitivity; full original- and log-scale agreement details; mean-glucose paired context; all count correlations and the source-span NOT_ESTIMABLE rows; the three specification-shift diagnostics; the full machine-readable/tabular 18-row landscape that underlies Candidate Figure 3; source-restricted cohort details; alternative GV metrics; MIMIC 365-day PH and time-varying material; MICE/model diagnostics; the original mean-glucose/absolute-risk contrast; and the full INSPIRE coverage/timing boundary audit. The evidence map records this assignment row by row.

## 5. Which prior main-manuscript results should be demoted?

The original mean-glucose conditioning/absolute-risk Figure 3 result (`ORIG_MIMIC_ABSOLUTE_RISK_30D`) should move to the supplement: it is valid secondary evidence but is not the revised story's central estimand. The MIMIC 365-day average-HR sequence should move to the supplement and carry the documented PH violation (`ORIG_MIMIC_365D_PH_DIAGNOSTIC`; `MIMIC_MODEL_B_365D`). The old cross-database non-pooled display should be refocused into the conceptual landscape rather than treated as a common validation effect. Withdrawn INSPIRE 365-day and historical coverage-gate results remain excluded.

## 6. Should source dependence be the centerpiece?

Yes. The same-patient comparison provides the cleanest evidence that changing the observed measurement source can change derived GV and its mortality coefficient while holding the patient/event cohort fixed. The interpretation remains measurement-context dependence, not source accuracy and not a causal source-by-exposure interaction.

## 7. Is ‘sampling opportunity’ too strong for the title?

Yes. It is too strong as a title-level causal-sounding explanation because source-specific span was not estimable and measurement-count imbalance alone did not account for source-defined GV disagreement. Sampling opportunity should remain a bounded context dimension in the Results and Supplement, not the title's headline mechanism.

## Findings that could weaken the framing

The broad original-scale limits of agreement and broad multiplicative ratio limits (`P2A_AGREEMENT_PRIMARY_GV_MEAN_DIFFERENCE`; `P2B_LOG_PRIMARY_RATIO_LOA_LOWER`; `P2B_LOG_PRIMARY_RATIO_LOA_UPPER`) limit claims of precise interchangeability. Source-specific span is not estimable (`P2B_SPAN_PRIMARY_NOT_ESTIMABLE`), the count closure is descriptive rather than explanatory (`P2B_COUNT_PRIMARY_GV_CORRELATION`), INSPIRE estimates are event-limited, and eICU preserves aggregate rather than event-level glucose inputs. These findings weaken generalization and mechanism claims, but they do not materially overturn the narrower measurement-context framing.

## 8. What exact title is recommended?

**Measurement-Context Dependence of Routine-Care Glycemic Variability After Cardiac Surgery: A Multidatabase Cohort Study**

## 9. What exact one-sentence central claim is recommended?

In routine-care cardiac-surgery data, SD-based glycemic variability was measurement-context dependent: measurement source materially altered both derived GV values and the estimated mortality association in identical patients, whereas observed measurement-count differences alone did not account for the source-defined disagreement and source-specific span was not estimable.

## 10. What contradictions exist between the old manuscript and the revised evidence package?

- The old primary null/mean-glucose framing is not numerically invalid for final MIMIC Model B, but it is no longer the central interpretation; it is subordinate to the source-dependence evidence.
- The old absolute-risk Figure 3 remains a valid complete-case secondary contrast, but its manuscript role changes from main figure to supplement.
- The old 365-day average-HR material remains an audit result only because final provenance documents PH violation; it must not be presented as an uncomplicated long-horizon primary effect.
- The submitted Figure 5 provenance included a harmonized MIMIC comparison frame of N=8,117/128, whereas the final Phase 2B landscape uses the final MIMIC primary/source contexts and the locked eICU M3/M4 N=7,115/130 frame. These are different non-pooled analytic contexts, not a single cross-database estimand and not a biological contradiction.
- The INSPIRE final I1/I2/I3 and corrected 48-hour results are retained with the final v5 administrative-censoring and coverage-attestation boundaries; withdrawn 365-day results and historical wide-mask/discharge-censoring outputs are not part of the revised evidence.
- The candidate figures are provenance-locked but remain candidate production artifacts; missing public-backup figure builders are a Phase 3 production task, not a license to introduce manually typed numbers.

## Lock interpretation

The total evidence supports ‘measurement-context dependence’ as the revised manuscript's central framing. It does not support the stronger claim that sampling intensity explains the GV association, a reference-standard hierarchy between POCT and laboratory measurements, causal source interactions, database pooling, or biological heterogeneity claims.

GO_PHASE3_MANUSCRIPT_REWRITE
