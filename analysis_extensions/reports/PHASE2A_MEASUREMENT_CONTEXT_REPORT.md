# Phase 2A — Measurement-Context Analyses

Final decision: **GO_PHASE2B_ANALYTIC_CONTEXT_LANDSCAPE**

## Scope and lineage

All analyses use the final JAHA v5 lineage declared in `config.yaml`. No historical extension input was used, and no manuscript or authoritative source repository was modified by this task.

Phase 1.6 hard gate: `GO_PHASE2_MEASUREMENT_CONTEXT_ANALYSES`; new modeling remains outside the scope of this report.

## A. Same-patient measurement-source agreement

The primary comparison is the 452-patient final-target intersection; the 453-patient all-paired cohort is a sensitivity denominator. Paired differences are POCT minus central-laboratory values.

- Primary GV SD: N=452; POCT mean=28.6203; laboratory mean=34.6317; Pearson r=0.4893; Spearman rho=0.3269; paired mean difference=-6.0114; 95% limits of agreement [-69.9302, 57.9075].
- Primary GV difference versus pair mean: Pearson r=-0.5278; Spearman rho=-0.3380.
- The negative difference-versus-pair-mean correlations indicate scale-dependent/proportional disagreement: POCT minus laboratory GV tends to become more negative at larger paired GV values; this is not consistent with a constant offset alone.
- All-paired GV sensitivity: N=453; paired mean difference=-5.9386; Pearson r=0.4881; Spearman rho=0.3263.
- Primary mean-glucose comparison: POCT mean=151.0934; laboratory mean=155.1507; Pearson r=0.8051; paired mean difference=-4.0573.

These are diagnostic agreement and source-comparability results. They do not establish a gold-standard device, analytical validity, or pure analytical measurement error.

## B. Paired source-definition coefficient difference

The authoritative final source-model contract was reproduced for both source definitions at N=409 and 49 events. Observed HR per 10-unit GV SD: POCT 0.7589; laboratory 0.9709; observed delta-beta=-0.2464; HR ratio=0.7816.
Paired bootstrap: B=2000; successful=2000; failed=0 (0.0000%). Delta-beta median=-0.2687; SD=0.1477; percentile 95% CI [-0.6323, -0.0551].
The percentile interval excludes zero, so the source-specific coefficient difference is not negligible descriptively in this sensitivity analysis; it is not evidence of a causal interaction.
This paired difference is a descriptive sensitivity analysis of source-definition coefficients, not a causal interaction or effect-modification test.

## C. Sampling-process dependence

MIMIC-IV uses the final 10,561-patient / 296-event target. eICU uses the final aggregate M3/M4 input of 7,115 patients / 130 events; event-level eICU sampling-process analysis was not available. INSPIRE was conservatively not assigned a new sampling-process analysis because an unambiguous validated sampling-count/span field was not established.

MIMIC GV versus glucose count: Pearson r=0.1508; Spearman rho=0.1756; GV versus span: Pearson r=0.1248; Spearman rho=0.1376.
eICU aggregate GV versus glucose count: Pearson r=0.1761; Spearman rho=0.2176; GV versus span: Pearson r=0.1388; Spearman rho=0.1507.

## D. Locked specification shifts

These are sampling-process specification shifts reported as before/after effects and delta log-effect. No attenuation percentage is calculated.

- MIMIC-IV: MIMIC Model B HR 0.9794 [0.9225, 1.0398] to MIMIC Model C HR 1.0286 [0.9671, 1.0940]; delta log-effect=0.0490; null crossing=TRUE.
- eICU-CRD: eICU M3 RR 1.1577 [1.0246, 1.3082] to eICU M4 RR 1.1328 [0.9902, 1.2960]; delta log-effect=-0.0218; null crossing=FALSE.
- INSPIRE: INSPIRE I2 HR 0.9046 [0.6526, 1.2539] to INSPIRE I3 HR 0.9285 [0.6724, 1.2822]; delta log-effect=0.0261; null crossing=FALSE.

## E. Interpretation boundary and Phase 2B decision

The combined evidence is an analytic-context landscape: observed GV associations can depend on glucose-source definition, the opportunity to observe glucose, and locked specification choices. These analyses do not identify a causal measurement-error mechanism, establish interchangeability, or support pooled cross-database inference/heterogeneity claims.

Results that weaken a stronger interpretation are the only-moderate same-patient GV agreement with wide limits of agreement, the source-specific coefficient difference, modest rather than dominant count/span correlations, and the preserved aggregate-only eICU structure. These findings support contextual caution rather than a claim that one measurement process fully explains the outcome association.

Decision: **GO_PHASE2B_ANALYTIC_CONTEXT_LANDSCAPE**.
