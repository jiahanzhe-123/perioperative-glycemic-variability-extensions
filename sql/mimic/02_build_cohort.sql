-- 02_build_cohort.sql
-- Index time is the earliest primary cardiac procedure chartdate. MIMIC-IV
-- procedures_icd contains chartdate rather than a procedure charttime, so the
-- time precision is explicitly retained as date_only.

BEGIN;
SET LOCAL statement_timeout = 0;
SET LOCAL search_path = public;

DROP TABLE IF EXISTS mimic_custom.cardiac_surgery_cohort_v2;
DROP TABLE IF EXISTS mimic_custom.cardiac_surgery_cohort_all_candidates_v2;

CREATE TABLE mimic_custom.cardiac_surgery_cohort_all_candidates_v2 AS
WITH procedure_events AS (
    SELECT p.subject_id,
           p.hadm_id,
           p.seq_num,
           p.chartdate::timestamp AS surgery_time,
           'date_only'::text AS surgery_time_precision,
           p.icd_code::text AS procedure_icd_code,
           p.icd_version AS procedure_icd_version,
           cb.procedure_title,
           cb.procedure_group,
           cb.open_surgery_candidate,
           cb.catheter_based_candidate,
           cb.include_primary,
           CASE cb.procedure_group
             WHEN 'cabg' THEN 1
             WHEN 'open_valve' THEN 2
             WHEN 'aortic_surgery' THEN 3
             WHEN 'transplant_vad' THEN 4
             WHEN 'congenital_cardiac' THEN 5
             WHEN 'transcatheter_valve' THEN 6
             WHEN 'pci' THEN 7
             WHEN 'other_cardiac' THEN 8
             ELSE 99
           END AS procedure_priority
    FROM mimiciv_hosp.procedures_icd p
    JOIN mimic_custom.cardiac_surgery_codebook_v2 cb
      ON cb.icd_version = p.icd_version
     AND cb.icd_code = btrim(p.icd_code::text)
     AND cb.include_primary
), procedure_with_icu AS (
    SELECT pe.*,
           a.admittime,
           a.dischtime,
           a.deathtime,
           a.hospital_expire_flag,
           pat.gender,
           pat.anchor_age,
           pat.anchor_year,
           pat.dod,
           (pat.anchor_age + EXTRACT(YEAR FROM a.admittime)::integer - pat.anchor_year)::integer AS age_at_admission,
           icu.stay_id,
           icu.first_careunit,
           icu.last_careunit,
           icu.intime AS matched_icu_intime,
           icu.outtime AS matched_icu_outtime,
           icu.los AS matched_icu_los,
           CASE WHEN icu.stay_id IS NULL THEN false ELSE true END AS has_icu_match,
           CASE WHEN icu.stay_id IS NULL THEN 'no ICU stay overlaps or follows procedure date within 48 hours'
                WHEN icu.intime <= pe.surgery_time AND icu.outtime >= pe.surgery_time THEN 'ICU stay covers date-only procedure index time'
                WHEN icu.intime >= pe.surgery_time THEN 'nearest ICU stay begins after date-only procedure index time'
                ELSE 'nearest ICU stay overlaps preceding date-only procedure window' END AS selected_stay_reason
    FROM procedure_events pe
    LEFT JOIN mimiciv_hosp.admissions a
      ON a.subject_id = pe.subject_id AND a.hadm_id = pe.hadm_id
    LEFT JOIN mimiciv_hosp.patients pat
      ON pat.subject_id = pe.subject_id
    LEFT JOIN LATERAL (
        SELECT i.stay_id,
               i.first_careunit,
               i.last_careunit,
               i.intime,
               i.outtime,
               i.los
        FROM mimiciv_icu.icustays i
        WHERE i.subject_id = pe.subject_id
          AND i.hadm_id = pe.hadm_id
          AND i.intime <= pe.surgery_time + interval '48 hours'
          AND i.outtime >= pe.surgery_time - interval '24 hours'
        ORDER BY CASE WHEN i.intime <= pe.surgery_time AND i.outtime >= pe.surgery_time THEN 0
                      WHEN i.intime >= pe.surgery_time THEN 1 ELSE 2 END,
                 abs(EXTRACT(EPOCH FROM (i.intime - pe.surgery_time))),
                 i.stay_id
        LIMIT 1
    ) icu ON true
), procedure_summary AS (
    SELECT subject_id,
           hadm_id,
           count(*)::integer AS n_cardiac_procedures,
           string_agg(DISTINCT procedure_group, '|' ORDER BY procedure_group) AS procedure_group_all,
           bool_or(procedure_group = 'cabg') AS cabg_flag,
           bool_or(procedure_group = 'pci') AS pci_flag,
           bool_or(procedure_group = 'open_valve') AS open_valve_flag,
           bool_or(procedure_group = 'transcatheter_valve') AS transcatheter_valve_flag,
           bool_or(procedure_group = 'aortic_surgery') AS aortic_surgery_flag,
           bool_or(procedure_group = 'transplant_vad') AS transplant_vad_flag,
           bool_or(procedure_group = 'congenital_cardiac') AS congenital_cardiac_flag,
           bool_or(open_surgery_candidate) AS open_surgery_flag,
           bool_or(catheter_based_candidate) AS catheter_based_flag,
           bool_or(open_surgery_candidate) AND bool_or(catheter_based_candidate) AS mixed_open_and_catheter_flag,
           bool_and(procedure_group = 'pci') AS pure_pci_flag
    FROM procedure_with_icu
    GROUP BY subject_id, hadm_id
), ranked AS (
    SELECT pw.*,
           ps.n_cardiac_procedures,
           ps.procedure_group_all,
           ps.cabg_flag,
           ps.pci_flag,
           ps.open_valve_flag,
           ps.transcatheter_valve_flag,
           ps.aortic_surgery_flag,
           ps.transplant_vad_flag,
           ps.congenital_cardiac_flag,
           ps.open_surgery_flag,
           ps.catheter_based_flag,
           ps.mixed_open_and_catheter_flag,
           ps.pure_pci_flag,
           row_number() OVER (
             PARTITION BY pw.subject_id, pw.hadm_id
             ORDER BY pw.has_icu_match DESC, pw.surgery_time, pw.procedure_priority,
                      pw.seq_num, pw.procedure_icd_code
           ) AS candidate_rank_within_hadm,
           row_number() OVER (
             PARTITION BY pw.subject_id
             ORDER BY pw.has_icu_match DESC, pw.surgery_time, pw.procedure_priority,
                      pw.seq_num, pw.hadm_id, pw.procedure_icd_code
           ) AS candidate_rank_within_subject
    FROM procedure_with_icu pw
    JOIN procedure_summary ps USING (subject_id, hadm_id)
)
SELECT subject_id,
       hadm_id,
       stay_id,
       age_at_admission,
       (age_at_admission >= 18) AS adult_flag,
       (age_at_admission >= 75) AS age_ge_75,
       gender,
       admittime AS admittime_ts,
       dischtime AS dischtime_ts,
       deathtime,
       dod,
       hospital_expire_flag,
       surgery_time,
       surgery_time_precision,
       matched_icu_intime AS first_icu_intime,
       matched_icu_outtime AS matched_icu_outtime,
       procedure_group AS procedure_group_main,
       procedure_group_all,
       procedure_icd_code AS procedure_icd_code_main,
       procedure_icd_version AS procedure_icd_version_main,
       procedure_title AS procedure_title_main,
       n_cardiac_procedures,
       cabg_flag,
       pci_flag,
       open_valve_flag,
       transcatheter_valve_flag,
       aortic_surgery_flag,
       transplant_vad_flag,
       congenital_cardiac_flag,
       false AS support_only_flag,
       open_surgery_flag,
       catheter_based_flag,
       mixed_open_and_catheter_flag,
       pure_pci_flag AS exclude_pci_sensitivity_flag,
       (EXTRACT(EPOCH FROM (matched_icu_intime - surgery_time)) / 3600.0)::double precision AS surgery_to_icu_hours,
       (EXTRACT(EPOCH FROM (surgery_time - matched_icu_intime)) / 3600.0)::double precision AS icu_intime_to_surgery_hours,
       selected_stay_reason,
       candidate_rank_within_hadm,
       candidate_rank_within_subject,
       has_icu_match,
       include_primary AS primary_cardiac_procedure_flag,
       adm_to_surgery_days_placeholder
FROM (
    SELECT r.*,
           (EXTRACT(EPOCH FROM (r.surgery_time - r.admittime)) / 86400.0)::double precision AS adm_to_surgery_days_placeholder
    FROM ranked r
) q;

CREATE INDEX cardiac_surgery_all_v2_subject_idx
    ON mimic_custom.cardiac_surgery_cohort_all_candidates_v2 (subject_id);
CREATE INDEX cardiac_surgery_all_v2_hadm_idx
    ON mimic_custom.cardiac_surgery_cohort_all_candidates_v2 (hadm_id);
CREATE INDEX cardiac_surgery_all_v2_stay_idx
    ON mimic_custom.cardiac_surgery_cohort_all_candidates_v2 (stay_id);

CREATE TABLE mimic_custom.cardiac_surgery_cohort_v2 AS
WITH eligible AS (
    SELECT a.*,
           row_number() OVER (
             PARTITION BY subject_id
             ORDER BY surgery_time, candidate_rank_within_hadm, hadm_id, stay_id
           ) AS subject_first_rank
    FROM mimic_custom.cardiac_surgery_cohort_all_candidates_v2 a
    WHERE a.adult_flag
      AND a.has_icu_match
      AND a.primary_cardiac_procedure_flag
)
SELECT e.*,
       e.adm_to_surgery_days_placeholder AS adm_to_surgery_days,
       true AS default_first_cardiac_surgery_icu_stay,
       'procedure_icd_chartdate; nearest eligible ICU stay; first eligible stay per subject'::text AS cohort_source
FROM eligible e
WHERE subject_first_rank = 1;

ALTER TABLE mimic_custom.cardiac_surgery_cohort_v2
    DROP COLUMN adm_to_surgery_days_placeholder;
ALTER TABLE mimic_custom.cardiac_surgery_cohort_v2
    RENAME COLUMN matched_icu_outtime TO last_icu_outtime;

CREATE INDEX cardiac_surgery_cohort_v2_subject_idx
    ON mimic_custom.cardiac_surgery_cohort_v2 (subject_id);
CREATE INDEX cardiac_surgery_cohort_v2_hadm_idx
    ON mimic_custom.cardiac_surgery_cohort_v2 (hadm_id);
CREATE INDEX cardiac_surgery_cohort_v2_stay_idx
    ON mimic_custom.cardiac_surgery_cohort_v2 (stay_id);

COMMIT;

\copy (SELECT * FROM mimic_custom.cardiac_surgery_cohort_all_candidates_v2 ORDER BY subject_id, hadm_id, surgery_time, procedure_icd_code_main) TO 'outputs/qc/cardiac_surgery_cohort_all_candidates_v2.csv' WITH (FORMAT csv, HEADER true, ENCODING 'UTF8');
\copy (SELECT * FROM mimic_custom.cardiac_surgery_cohort_v2 ORDER BY subject_id) TO 'outputs/qc/cardiac_surgery_cohort_v2.csv' WITH (FORMAT csv, HEADER true, ENCODING 'UTF8');
