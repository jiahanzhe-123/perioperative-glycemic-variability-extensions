-- 01_build_cardiac_surgery_codebook.sql
-- Keyword/code-pattern screen only. The resulting table is a review artifact;
-- broad matches are not silently treated as a gold-standard phenotype.

BEGIN;
SET LOCAL statement_timeout = 0;
SET LOCAL search_path = public;

DROP TABLE IF EXISTS mimic_custom.cardiac_surgery_codebook_v2;
CREATE TABLE mimic_custom.cardiac_surgery_codebook_v2 (
    icd_version smallint NOT NULL,
    icd_code text NOT NULL,
    procedure_title text NOT NULL,
    procedure_group text,
    include_primary boolean NOT NULL,
    exclude_reason text,
    match_keyword text,
    open_surgery_candidate boolean NOT NULL,
    catheter_based_candidate boolean NOT NULL,
    support_only boolean NOT NULL,
    diagnostic_only boolean NOT NULL,
    confidence text NOT NULL,
    review_status text NOT NULL,
    codebook_rule_version text NOT NULL DEFAULT 'v2_keyword_screen_2026-07-18'
);

WITH src AS (
    SELECT icd_version,
           btrim(icd_code::text) AS icd_code,
           btrim(long_title) AS procedure_title,
           lower(regexp_replace(btrim(long_title), '\s+', ' ', 'g')) AS t
    FROM mimiciv_hosp.d_icd_procedures
), flags AS (
    SELECT s.*,
           (t ~ '(intra[- ]aortic balloon|intraaortic balloon|extracorporeal membrane|\becmo\b|mechanical ventilation|ventilator|dialysis|hemodialysis|\bpacemaker\b|cardioverter|defibrillator|temporary transvenous pacemaker)'
             AND t !~ '(heart transplant|ventricular assist|heart assist|total artificial heart)') AS is_support,
           (t ~ '(catheterization|angiograph|angiogram|diagnostic|cardiac mapping|biopsy of heart|echocardiography)'
             AND t !~ '(angioplasty|stent|atherectomy|bypass|valve replacement|valve repair|transcatheter|therapeutic)') AS is_diagnostic,
           (t ~ '(coronary artery bypass|aortocoronary bypass|bypass.*coronary|internal mammary-coronary|heart revascularization)') AS is_cabg,
           (t ~ '(percutaneous transluminal coronary|coronary.*stent|stent.*coronary|coronary angioplasty|coronary atherectomy|dilation of coronary artery|extirpation of coronary artery obstruction|removal of coronary artery obstruction)') AS is_pci,
           (t ~ '(transcatheter|transapical|endovascular.*replacement.*valve|percutaneous.*valve|\btavr\b|\btavi\b)') AS is_trans_valve,
           (t ~ '(open heart.*valv|open and other replacement.*valve|replacement.*valve|repair.*valve|valvuloplasty)') AS is_open_valve,
           (t ~ '(aortic|aorta|thoracic aorta|ascending aorta|aortic root|aortic arch)'
             AND t ~ '(repair|replacement|resection|excision|aneurysm|dissection|bypass|arch|root|ascending)'
             AND t !~ '(angiograph|angiography|catheterization|diagnostic)') AS is_aortic,
           (t ~ '(heart transplant|ventricular assist|heart assist|total artificial heart)') AS is_transplant_vad,
           (t ~ '(congenital|septal defect|tetralogy|fontan|coarctation|patent ductus|atrial septal|ventricular septal|anomalous pulmonary)') AS is_congenital,
           (t ~ '(\bheart\b|cardiac|pericard|myocard|coronary|valve)'
             AND t ~ '(repair|replacement|resection|excision|operation|revascularization|implant|graft|valvotomy|valvuloplasty)') AS is_other_cardiac
    FROM src s
), classified AS (
    SELECT f.*,
           CASE
             WHEN is_support THEN 'support_only'
             WHEN is_diagnostic THEN 'diagnostic_only'
             WHEN is_cabg THEN 'cabg'
             WHEN is_pci THEN 'pci'
             WHEN is_trans_valve THEN 'transcatheter_valve'
             WHEN is_open_valve THEN 'open_valve'
             WHEN is_aortic THEN 'aortic_surgery'
             WHEN is_transplant_vad THEN 'transplant_vad'
             WHEN is_congenital THEN 'congenital_cardiac'
             WHEN is_other_cardiac THEN 'other_cardiac'
             ELSE NULL
           END AS grp,
           concat_ws('; ',
             CASE WHEN is_cabg THEN 'cabg_keyword' END,
             CASE WHEN is_pci THEN 'pci_keyword' END,
             CASE WHEN is_open_valve THEN 'open_valve_keyword' END,
             CASE WHEN is_trans_valve THEN 'transcatheter_valve_keyword' END,
             CASE WHEN is_aortic THEN 'aortic_keyword' END,
             CASE WHEN is_transplant_vad THEN 'transplant_vad_keyword' END,
             CASE WHEN is_congenital THEN 'congenital_keyword' END,
             CASE WHEN is_support THEN 'support_only_keyword' END,
             CASE WHEN is_diagnostic THEN 'diagnostic_only_keyword' END,
             CASE WHEN is_other_cardiac THEN 'other_cardiac_keyword' END
           ) AS matched_by
    FROM flags f
)
INSERT INTO mimic_custom.cardiac_surgery_codebook_v2
    (icd_version, icd_code, procedure_title, procedure_group, include_primary,
     exclude_reason, match_keyword, open_surgery_candidate,
     catheter_based_candidate, support_only, diagnostic_only, confidence,
     review_status)
SELECT icd_version,
       icd_code,
       procedure_title,
       grp,
       (grp IS NOT NULL AND grp NOT IN ('support_only','diagnostic_only')) AS include_primary,
       CASE WHEN grp = 'support_only' THEN 'support-only procedure; not a primary cardiac surgery criterion'
            WHEN grp = 'diagnostic_only' THEN 'diagnostic catheterization/imaging/biopsy only'
            WHEN grp IS NULL THEN 'no cardiac procedure keyword match'
            WHEN grp = 'other_cardiac' THEN 'broad cardiac keyword; manual review required'
            ELSE NULL END AS exclude_reason,
       COALESCE(NULLIF(matched_by,''), 'no_match') AS match_keyword,
       (grp IN ('cabg','open_valve','aortic_surgery','transplant_vad','congenital_cardiac','other_cardiac')) AS open_surgery_candidate,
       (grp IN ('pci','transcatheter_valve')) AS catheter_based_candidate,
       (grp = 'support_only') AS support_only,
       (grp = 'diagnostic_only') AS diagnostic_only,
       CASE WHEN grp IN ('cabg','pci','open_valve','transcatheter_valve','transplant_vad') THEN 'high'
            WHEN grp IN ('aortic_surgery','congenital_cardiac') THEN 'moderate'
            WHEN grp = 'other_cardiac' THEN 'low'
            ELSE 'not_applicable' END AS confidence,
       CASE WHEN grp IN ('cabg','pci','open_valve','transcatheter_valve','transplant_vad') THEN 'auto_candidate'
            WHEN grp IS NULL THEN 'not_candidate'
            ELSE 'manual_review' END AS review_status
FROM classified
WHERE grp IS NOT NULL;

CREATE INDEX cardiac_surgery_codebook_v2_code_idx
    ON mimic_custom.cardiac_surgery_codebook_v2 (icd_version, icd_code);

COMMIT;

\copy (SELECT * FROM mimic_custom.cardiac_surgery_codebook_v2 ORDER BY icd_version, procedure_group, icd_code) TO 'outputs/qc/cardiac_surgery_codebook_review.csv' WITH (FORMAT csv, HEADER true, ENCODING 'UTF8');
