# Phase 2B — Analytic-Context Landscape and Source-Dependence Closure

## Material Passport

- Phase: 2B analytic-context landscape
- Lineage: JAHA v5 final lineage only
- Inputs: Phase 1.6 validated context, locked final result files, final same-patient source file
- Status: reproducible extension package; no manuscript or authoritative code modification
- Decision: Phase 3 manuscript rewrite is permitted only under the bounded framing below

## Final decision: GO_PHASE3_MANUSCRIPT_REWRITE

## 1. Multiplicative-scale source agreement

Primary log-scale comparison used N=447 positive paired observations. Mean log ratio=0.0303; median log ratio=-0.0012; geometric mean ratio=1.0308.
95% limits on the log scale [-2.2440, 2.3047]; exponentiated ratio limits [0.1060, 10.0207].
Correlation of log ratio with log pair mean: Pearson r=-0.3198; Spearman rho=-0.3324. Original-scale difference/pair-mean correlations were Pearson r=-0.5278 and Spearman rho=-0.3380. The log-scale sensitivity retains the direction of the original-scale scale-dependence signal and does not overturn the Phase 2A disagreement conclusion.

The log-scale analysis is a sensitivity representation, not a replacement for the original-scale Bland–Altman analysis. It does not turn either routine-care source into a reference standard.

## 2. Sampling-count/span imbalance versus source GV disagreement

For the primary 452-patient cohort, POCT minus laboratory source-count difference had mean=2.2279; SD=2.7576; median=2.0000; IQR=4.0000; Q5/Q95=[-1.0000, 7.0000].
Correlation of delta GV with delta count: primary Pearson r=-0.0204; Spearman rho=0.0556; sensitivity Pearson r=-0.0222; Spearman rho=0.0520.
This near-zero source-difference/count-difference relationship is distinct from the positive-but-modest patient-level GV-versus-count correlations reported in Phase 2A.
Source-specific span is `NOT_ESTIMABLE`: the validated final same-patient source file contains source counts but no source-specific span fields. No undocumented span was reconstructed. Therefore the available count component does not establish that sampling imbalance explains source-GV disagreement, and the full count/span question remains incomplete rather than solved.

These are descriptive source-process comparisons, not causal adjustment models and not evidence that count or span is a confounder.

## 3. Analytic-context landscape

The master table contains 18 ordered, SD-GV-only rows grouped by context family: MIMIC adjustment (A/B/C), MIMIC measurement source (priority, POCT-only, laboratory-only, blood-gas-only, common-source, and the two same-patient rows), INSPIRE timing/adjustment (I1/I2/I3 and corrected 48-hour landmark), and eICU adjustment (M1-M4). Rows are presented in conceptual order, never sorted by effect size.

Same-patient outcome-model estimates were POCT HR=0.7589 and laboratory HR=0.9709; their paired log-effect difference was -0.2464 (ratio descriptor 0.7816).
For comparison, MIMIC B→C delta log-effect=0.0490; eICU M3→M4=-0.0218; INSPIRE I2→I3=0.0261. The largest absolute movement among the prespecified comparisons is the same-patient source definition: delta log-effect=-0.2464.
No common summary effect, pooled estimate, or formal cross-database heterogeneity statistic was calculated. eICU estimates remain RR and MIMIC/INSPIRE estimates remain HR.

## 4. Does the evidence support measurement-context dependence as the central framing?

Yes, with a bounded formulation. The evidence supports describing routine-care SD-GV as measurement-context dependent: source-defined GV shows only moderate within-patient agreement, disagreement is scale-dependent, source definition changes the estimated mortality coefficient in an identical-patient comparison, and simple count relationships are positive but modest.

The evidence does not support the stronger claim that sampling intensity explains the GV association. No single observed process fully explains the heterogeneity, and source-specific differences are not causal interactions.

## 5. Main manuscript versus supplement

Main manuscript candidate content:
- A compact same-patient source-dependence panel: original-scale scatter/Bland–Altman, with the log-scale result identified as sensitivity rather than replacement.
- The identical-patient POCT versus laboratory mortality coefficients and paired delta-beta result, explicitly labelled descriptive source-definition sensitivity.
- A compact four-panel analytic-context landscape or selected context-family summary, retaining HR/RR labels and conceptual ordering.

Supplement candidate content:
- The full 18-row machine-readable landscape and provenance map.
- Full log-scale agreement statistics, source-count imbalance distributions/correlations, the NOT_ESTIMABLE span boundary, all source-specific rows, and QC/hash evidence.
- The complete candidate figures and exact locked-result selectors.

## 6. Findings that materially weaken the framing

The moderate Pearson/Spearman source agreement and broad original-scale limits of agreement weaken any claim of interchangeable GV sources. The log-scale sensitivity does not erase the disagreement pattern. The validated file lacks source-specific span, so the sampling-imbalance closure is incomplete. The count relationships are modest, eICU remains aggregate-only, and the source-specific coefficient difference cannot be interpreted causally. These are reasons for a cautious measurement-context framing, not reasons to call POCT wrong or laboratory true.

## 7. Final candidate title and central claim

Candidate title: **Routine-Care Glycemic Variability After Cardiac Surgery: Dependence on Measurement Source, Sampling Opportunity, and Analytic Context**

Central claim: **In routine-care cardiac-surgery data, SD-based glycemic variability is a measurement-context-dependent biomarker whose source, observation opportunity, and analytic specification can change its apparent mortality association, while no single observed process fully explains the differences.**

GO_PHASE3_MANUSCRIPT_REWRITE
