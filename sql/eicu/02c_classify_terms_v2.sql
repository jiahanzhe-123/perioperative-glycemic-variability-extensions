-- Step 2c: 修正版分类(修正 \b 为 [[:<:]],并优先使用 eICU 层级路径结构)
DROP TABLE IF EXISTS eicu_cardio_validation.cand_terms_classified;
CREATE TABLE eicu_cardio_validation.cand_terms_classified AS
WITH n AS (
  SELECT source_table, source_column, original_value, frequency,
         lower(regexp_replace(trim(original_value), '\s+', ' ', 'g')) AS normalized_value
  FROM eicu_cardio_validation.cand_terms_raw
),
c AS (
  SELECT *,
    -- ===== 路径结构标记 =====
    (normalized_value LIKE '%non-operative%') AS path_nonoperative,
    (normalized_value ~ 'cardiovascular\|cardiac surgery\|' OR normalized_value LIKE 'cardiovascular|cardiac surgery|%') AS path_cardiac_surgery,
    (normalized_value ~ 'cardiovascular\|post vascular surgery\|' OR normalized_value LIKE 'cardiovascular|post vascular surgery|%'
     OR normalized_value ~ 'cardiovascular\|vascular surgery\|' OR normalized_value LIKE 'cardiovascular|vascular surgery|%') AS path_vascular_surgery,
    (normalized_value LIKE 'pulmonary|post thoracic surgery%' OR normalized_value LIKE 'pulmonary|thoracic surgery%') AS path_thoracic_surgery,
    (normalized_value LIKE 'admission diagnosis|all diagnosis|operative%') AS path_adm_operative,
    (normalized_value LIKE 'admission diagnosis|all diagnosis|non-operative%') AS path_adm_nonoperative,
    -- ===== 非开放操作 =====
    (normalized_value ~ 'tavr|tavi|transcatheter|percutaneous|ptca|[[:<:]]pci[[:>:]]|angioplasty|stent|mitraclip|[[:<:]]teer[[:>:]]') AS is_transcath,
    (normalized_value ~ 'cardiac cath|heart cath|coronary angiograph|angiogram|cardiac catheterization|catheterization|cardiac angiograph') AS is_cath,
    (normalized_value ~ 'ablation or mapping|catheter ablation|electrophysiolog|[[:<:]]eps[[:>:]]|cardioversion|pacemaker|aicd|[[:<:]]icd[[:>:]]|defibrillator|cardiac hardware') AS is_ep_device,
    (normalized_value ~ 'iabp|balloon pump|intraaortic balloon|intra-aortic balloon|[[:<:]]ecmo[[:>:]]|extracorporeal') AS is_support_only,
    (normalized_value ~ 'abdominal aort|femoral|popliteal|aorto-?femoral|aorto-?iliac|aortobifem|peripheral|carotid|renal bypass|mesenteric|endarterectomy|graft for dialysis|infected vascular') AS is_noncardiac_vasc,
    (normalized_value ~ 'liver transplant|kidney|renal transplant|lung transplant|kidney-pancreas|pancreas transplant') AS is_noncardiac_tx,
    (normalized_value ~ 'intracranial|subarachnoid|cerebral|thoracotomy|thoracoscopic|lobectomy|wedge resection|decortication|pleurodesis') AS is_noncardiac_other,
    (normalized_value ~ 'minimally invasive|mid-cabg|midcabg|off-pump|robotic|keyhole') AS is_minimally_invasive,
    (normalized_value ~ 'consultation|diagnostic ultrasound|echocardiograph|foley|nasogastric|chest tube|surgical drains|central venous|vascular catheter|arterial line|swan|pulmonary artery catheter') AS is_nonsurg_noise,
    (normalized_value ~ 'thrombectomy|embolectomy|dilatation') AS is_vascular_proc,
    -- ===== 开放手术类别 =====
    (normalized_value ~ 'heart transplant|cardiac transplant|[[:<:]]lvad[[:>:]]|[[:<:]]rvad[[:>:]]|[[:<:]]bivad[[:>:]]|ventricular assist|[[:<:]]vad[[:>:]]|assist device|artificial heart|heart transplant complication') AS is_tx_vad,
    (normalized_value ~ 'cabg|coronary artery bypass|coronary bypass|aortocoronary bypass|internal mammary') AS is_cabg,
    (normalized_value ~ 'aortic valve|mitral valve|tricuspid valve|pulmonary valve|pulmonic valve|valve replacement|valve repair|valve surgery|open valv|bentall|ross procedure|valve sparing|prosthetic valve|valve replacement or repair') AS is_valve,
    ((normalized_value ~ 'thoracic aort|ascending aort|aortic arch|descending thoracic|aortic dissection|dissecting aortic|aneurysm resection/repair-thoracic|aortic root|aortic replacement|elephant trunk')
      AND normalized_value !~ 'abdominal') AS is_aortic,
    (normalized_value ~ 'atrial septal|ventricular septal|[[:<:]]asd[[:>:]]|[[:<:]]vsd[[:>:]]|congenital|tetralogy|fontan|coarctation|patent ductus|septal defect') AS is_congenital,
    (normalized_value ~ 'open-heart|open heart|sternotomy|ventricular aneurysm|maze|pericardial window|pericardiectomy|pericardial surgery|intracardiac|redo cardiac surgery|complications of previous open|mediastinitis|unstable sternum|post-pericardiotomy|phrenic nerve|tumor removal, intracardiac|surgery for congenital heart disease') AS is_other_open,
    -- ===== 模糊 =====
    (normalized_value ~ 'cardiovascular surgery, other|vascular surgery, other|aneurysms, repair of other|aneurysm/pseudoaneurysm, other|aneurysm resection / repair|vascular bypass|pericardial effusion|tamponade') AS is_ambiguous
  FROM n
)
SELECT *,
  CASE
    -- 优先级:结构性排除 -> 移植/VAD -> 微创排除 -> 各开放类别 -> 模糊 -> 非开放 -> UNCERTAIN
    WHEN path_nonoperative OR path_adm_nonoperative THEN 'EXCLUDE_NONOPEN'
    WHEN is_transcath AND NOT (is_cabg OR is_valve) THEN 'EXCLUDE_NONOPEN'
    WHEN is_tx_vad THEN 'TRANSPLANT_VAD'
    WHEN is_minimally_invasive AND is_cabg THEN 'EXCLUDE_NONOPEN'
    WHEN is_cabg THEN 'CABG'
    WHEN is_valve THEN 'OPEN_VALVE'
    WHEN is_aortic THEN 'OPEN_AORTIC'
    WHEN is_congenital THEN 'CONGENITAL'
    WHEN is_other_open AND NOT (is_ep_device OR is_cath) THEN 'OTHER_OPEN'
    WHEN path_vascular_surgery OR path_thoracic_surgery THEN 'EXCLUDE_NONOPEN'
    WHEN is_noncardiac_vasc OR is_noncardiac_tx OR is_noncardiac_other THEN 'EXCLUDE_NONOPEN'
    WHEN is_ep_device OR is_cath OR is_support_only OR is_nonsurg_noise OR is_vascular_proc THEN 'EXCLUDE_NONOPEN'
    WHEN path_cardiac_surgery AND is_ep_device THEN 'EXCLUDE_NONOPEN'
    WHEN is_ambiguous THEN 'UNCERTAIN'
    WHEN path_cardiac_surgery THEN 'OTHER_OPEN'
    ELSE 'UNCERTAIN'
  END AS proposed_category,
  ((is_cabg::int + is_valve::int + is_aortic::int + is_tx_vad::int + is_congenital::int) > 1
   OR normalized_value ~ 'cabg with|cabg and valve|cabg redo|aortic / mitral|aortic and mitral|double valve') AS multi_category_flag
FROM c;
