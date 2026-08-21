-- Step 1b: 关键分类字段取值与缺失审计
\pset footer off
\o /tmp/audit_patient_cat.txt
SELECT 'gender' AS field, gender AS val, count(*) FROM patient GROUP BY 2 ORDER BY 3 DESC;
SELECT 'ethnicity', ethnicity, count(*) FROM patient GROUP BY 2 ORDER BY 3 DESC;
SELECT 'unittype', unittype, count(*) FROM patient GROUP BY 2 ORDER BY 3 DESC;
SELECT 'unitadmitsource', unitadmitsource, count(*) FROM patient GROUP BY 2 ORDER BY 3 DESC;
SELECT 'hospitaladmitsource', hospitaladmitsource, count(*) FROM patient GROUP BY 2 ORDER BY 3 DESC;
SELECT 'unitdischargestatus', unitdischargestatus, count(*) FROM patient GROUP BY 2 ORDER BY 3 DESC;
SELECT 'hospitaldischargestatus', hospitaldischargestatus, count(*) FROM patient GROUP BY 2 ORDER BY 3 DESC;
SELECT 'unitstaytype', unitstaytype, count(*) FROM patient GROUP BY 2 ORDER BY 3 DESC;
SELECT 'unitvisitnumber', unitvisitnumber::text, count(*) FROM patient GROUP BY 2 ORDER BY 2;
SELECT 'age_special', CASE WHEN age ~ '^[0-9]+$' THEN 'numeric' ELSE age END, count(*) FROM patient GROUP BY 2 ORDER BY 3 DESC;
\o /tmp/audit_apache.txt
SELECT 'apacheversion', apacheversion, count(*) FROM apachepatientresult GROUP BY 2 ORDER BY 3 DESC;
SELECT 'electivesurgery', electivesurgery::text, count(*) FROM apachepredvar GROUP BY 2 ORDER BY 2;
SELECT 'admitsource_code', admitsource::text, count(*) FROM apachepredvar GROUP BY 2 ORDER BY 2;
\o
