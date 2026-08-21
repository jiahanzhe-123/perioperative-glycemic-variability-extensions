-- Step 2b: 载入候选字符串并按规则分类
DROP TABLE IF EXISTS eicu_cardio_validation.cand_terms_raw;
CREATE TABLE eicu_cardio_validation.cand_terms_raw(
  source_table text, source_column text, original_value text, frequency int);
\copy eicu_cardio_validation.cand_terms_raw FROM '/tmp/cand_terms_raw.tsv' DELIMITER E'\t' CSV HEADER

DROP TABLE IF EXISTS eicu_cardio_validation.cand_terms_classified;
CREATE TABLE eicu_cardio_validation.cand_terms_classified AS
WITH n AS (
  SELECT source_table, source_column, original_value, frequency,
         lower(regexp_replace(trim(original_value), '\s+', ' ', 'g')) AS normalized_value
  FROM eicu_cardio_validation.cand_terms_raw
),
c AS (
  SELECT *,
    -- 非开放操作(优先判断)
    (normalized_value ~ 'tavr|tavi|transcatheter|percutaneous|ptca|\bpci\b|angioplasty|stent|mitraclip|teer') AS is_transcath,
    (normalized_value ~ 'cardiac cath|heart cath|coronary angiograph|angiogram|cardiac catheterization|catheterization' ) AND normalized_value !~ 'cabg|bypass|valve replacement|valve repair|transplant' AS is_cath_only,
    (normalized_value ~ 'ablation|electrophysiolog|\beps\b|cardioversion|pacemaker|aicd|\bicd\b|defibrillator implant|device implant') AS is_ep_device,
    (normalized_value ~ 'iabp|balloon pump|intra-aortic balloon|ecmo|extracorporeal') AND normalized_value !~ 'cabg|valve|transplant' AS is_support_only,
    (normalized_value ~ 'abdominal aort|femoral|popliteal|aorto-?femoral|aorto-?iliac|aorto-?bifem|peripheral|carotid|renal bypass|mesenteric|periphery' ) AS is_noncardiac_vasc,
    (normalized_value ~ 'liver transplant|kidney|renal transplant|lung transplant|kidney-pancreas|pancreas transplant') AS is_noncardiac_tx,
    (normalized_value ~ 'intracranial|subarachnoid|cerebral|hysterectomy|lymph node|oophorectomy|myxoma, non|noncardiac' ) AS is_noncardiac_other,
    (normalized_value ~ 'minimally invasive|mid-cabg|midcabg|off-pump|robotic') AS is_minimally_invasive,
    -- 开放手术类别
    (normalized_value ~ 'heart transplant|cardiac transplant|\blvad\b|\brvad\b|ventricular assist|\bvad\b|assist device|artificial heart') AS is_tx_vad,
    (normalized_value ~ 'cabg|coronary artery bypass|coronary bypass|aortocoronary bypass') AS is_cabg,
    (normalized_value ~ 'aortic valve|mitral valve|tricuspid valve|pulmonary valve|pulmonic valve|valve replacement|valve repair|valve surgery|valvuloplasty, open|open valv|bentall|ross procedure|valve sparing|homograft|prosthetic valve') AS is_valve,
    (normalized_value ~ 'thoracic aort|ascending aort|aortic arch|descending aort|aortic aneurysm|aortic dissection|aneurysm, dissecting aortic|aortic repair|aortic root|aortic replacement|elephant trunk|aortic graft') AND normalized_value !~ 'abdominal' AS is_aortic,
    (normalized_value ~ 'atrial septal|ventricular septal|\basd\b|\bvsd\b|congenital|tetralogy|fontan|coarctation|patent ductus') AS is_congenital,
    (normalized_value ~ 'open-heart|open heart|open cardiac|sternotomy|ventricular aneurysm|left ventricular aneurysm|cardiac aneurysm|maze procedure|atrial myxoma|cardiac tumor|myomectomy|pericardial window|pericardiectomy|cardiectomy' ) AS is_other_open
  FROM n
)
SELECT *,
  CASE
    WHEN is_tx_vad THEN 'TRANSPLANT_VAD'
    WHEN is_minimally_invasive AND is_cabg THEN 'EXCLUDE_NONOPEN'
    WHEN is_cabg AND is_valve THEN 'CABG'          -- 联合手术,主类别 CABG,stay 级再标 multi
    WHEN is_cabg THEN 'CABG'
    WHEN is_valve THEN 'OPEN_VALVE'
    WHEN is_aortic THEN 'OPEN_AORTIC'
    WHEN is_congenital THEN 'CONGENITAL'
    WHEN is_other_open THEN 'OTHER_OPEN'
    WHEN is_transcath OR is_cath_only OR is_ep_device OR is_support_only
      OR is_noncardiac_vasc OR is_noncardiac_tx OR is_noncardiac_other THEN 'EXCLUDE_NONOPEN'
    ELSE 'UNCERTAIN'
  END AS proposed_category,
  (is_cabg::int + is_valve::int + is_aortic::int + is_tx_vad::int + is_congenital::int) > 1 AS multi_category_flag
FROM c;
\o /tmp/cand_terms_classified.tsv
SELECT source_table, source_column, original_value, normalized_value, frequency, proposed_category, multi_category_flag
FROM eicu_cardio_validation.cand_terms_classified
ORDER BY proposed_category, frequency DESC;
\o
