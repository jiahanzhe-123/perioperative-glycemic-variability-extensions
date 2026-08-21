-- Step 2d: 导出最终候选术语表 outputs/candidate_cardiac_surgery_terms.csv
\o /tmp/candidate_cardiac_surgery_terms.csv
SELECT source_table, source_column, original_value, normalized_value, frequency,
       proposed_category,
       CASE proposed_category
         WHEN 'CABG' THEN 'INCLUDE'
         WHEN 'OPEN_VALVE' THEN 'INCLUDE'
         WHEN 'OPEN_AORTIC' THEN 'INCLUDE'
         WHEN 'TRANSPLANT_VAD' THEN 'INCLUDE'
         WHEN 'CONGENITAL' THEN 'INCLUDE'
         WHEN 'OTHER_OPEN' THEN 'INCLUDE'
         WHEN 'UNCERTAIN' THEN 'REVIEW_DEFAULT_EXCLUDE'
         WHEN 'EXCLUDE_NONOPEN' THEN 'EXCLUDE'
       END AS proposed_action,
       CASE
         WHEN proposed_category='CABG' AND multi_category_flag THEN 'CABG 相关字符串;联合手术标记 multi'
         WHEN proposed_category='CABG' THEN 'CABG 相关字符串'
         WHEN proposed_category='OPEN_VALVE' AND multi_category_flag THEN '开放瓣膜手术字符串;联合手术标记 multi'
         WHEN proposed_category='OPEN_VALVE' THEN '开放瓣膜手术字符串'
         WHEN proposed_category='OPEN_AORTIC' THEN '开放胸主动脉手术字符串(已排除腹主动脉)'
         WHEN proposed_category='TRANSPLANT_VAD' THEN '心脏移植/VAD 字符串'
         WHEN proposed_category='CONGENITAL' THEN '先心病矫治手术字符串'
         WHEN proposed_category='OTHER_OPEN' THEN '其他明确开放心脏手术或开放心脏术后并发症再手术'
         WHEN proposed_category='EXCLUDE_NONOPEN' AND normalized_value ~ 'minimally invasive|mid-cabg' THEN '微创 CABG,非正中开胸,按方案排除'
         WHEN proposed_category='EXCLUDE_NONOPEN' THEN '经导管/非心脏/血管/器械/支持性操作,按方案排除'
         ELSE '无明确开放手术证据(疾病诊断、药物、检查、非手术操作或部位不明的移植/动脉瘤),默认不纳入,仅人工核查用'
       END AS reason,
       multi_category_flag
FROM eicu_cardio_validation.cand_terms_classified
ORDER BY proposed_category, frequency DESC;
\o
