-- 00_schema_inventory.sql
-- Catalog-only inventory. It intentionally uses pg_class.reltuples estimates for
-- very large raw tables and does not force a full COUNT(*) scan.

BEGIN;
SET LOCAL statement_timeout = 0;
SET LOCAL lock_timeout = '5s';
SET LOCAL search_path = public;

DROP TABLE IF EXISTS mimic_custom.schema_inventory_v2;
CREATE TABLE mimic_custom.schema_inventory_v2 (
    inventory_ts timestamptz NOT NULL DEFAULT clock_timestamp(),
    table_schema text NOT NULL,
    table_name text NOT NULL,
    table_exists boolean NOT NULL,
    approx_rows bigint,
    n_columns integer,
    column_list text,
    row_count_method text NOT NULL,
    notes text
);

WITH required(table_schema, table_name, notes) AS (
    VALUES
      ('mimiciv_hosp','patients','demographics and date of death'),
      ('mimiciv_hosp','admissions','admission-level dates and in-hospital mortality'),
      ('mimiciv_hosp','transfers','hospital/ICU movement timeline'),
      ('mimiciv_hosp','procedures_icd','procedure codes and chartdate'),
      ('mimiciv_hosp','d_icd_procedures','procedure code dictionary'),
      ('mimiciv_hosp','diagnoses_icd','diagnoses by admission'),
      ('mimiciv_hosp','d_icd_diagnoses','diagnosis code dictionary'),
      ('mimiciv_hosp','labevents','hospital laboratory events'),
      ('mimiciv_hosp','d_labitems','laboratory item dictionary'),
      ('mimiciv_hosp','prescriptions','hospital medication orders'),
      ('mimiciv_hosp','emar','medication administration records'),
      ('mimiciv_hosp','emar_detail','administration details'),
      ('mimiciv_hosp','pharmacy','pharmacy orders'),
      ('mimiciv_hosp','services','service timeline'),
      ('mimiciv_hosp','omr','outpatient measurements'),
      ('mimiciv_icu','icustays','ICU stay grain'),
      ('mimiciv_icu','chartevents','ICU charted numeric/text events'),
      ('mimiciv_icu','d_items','ICU item dictionary'),
      ('mimiciv_icu','inputevents','ICU inputs and infusions'),
      ('mimiciv_icu','procedureevents','ICU procedures and support'),
      ('mimiciv_icu','outputevents','ICU outputs'),
      ('mimiciv_icu','datetimeevents','ICU date/time events')
), catalog AS (
    SELECT n.nspname AS table_schema,
           c.relname AS table_name,
           c.reltuples::bigint AS approx_rows
    FROM pg_class c
    JOIN pg_namespace n ON n.oid = c.relnamespace
    WHERE c.relkind IN ('r','p','v','m')
), cols AS (
    SELECT table_schema, table_name,
           count(*)::integer AS n_columns,
           string_agg(column_name || ':' || data_type, ', ' ORDER BY ordinal_position) AS column_list
    FROM information_schema.columns
    GROUP BY table_schema, table_name
)
INSERT INTO mimic_custom.schema_inventory_v2
    (table_schema, table_name, table_exists, approx_rows, n_columns, column_list,
     row_count_method, notes)
SELECT r.table_schema,
       r.table_name,
       (c.table_name IS NOT NULL),
       c.approx_rows,
       co.n_columns,
       co.column_list,
       CASE WHEN c.table_name IS NULL THEN 'not_applicable'
            ELSE 'pg_class.reltuples_estimate' END,
       r.notes
FROM required r
LEFT JOIN catalog c USING (table_schema, table_name)
LEFT JOIN cols co USING (table_schema, table_name)
ORDER BY r.table_schema, r.table_name;

COMMIT;

\copy (SELECT * FROM mimic_custom.schema_inventory_v2 ORDER BY table_schema, table_name) TO 'outputs/qc/schema_inventory.csv' WITH (FORMAT csv, HEADER true, ENCODING 'UTF8');
