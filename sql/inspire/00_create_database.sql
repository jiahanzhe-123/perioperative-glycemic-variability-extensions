-- 00_create_database.sql
-- INSPIRE v1.4.2 导入:schemas + DDL
-- 字段类型唯一依据:schema.csv(Number→double precision, String/M-F/1-5/Binary→text,
-- Relative Time→bigint,依据为发布方类型名与实际整型分钟格式的核对,见 schema_inventory.csv)。
-- 原始表只读;派生结果一律写 derived/analysis/audit。
BEGIN;

CREATE SCHEMA IF NOT EXISTS raw;
CREATE SCHEMA IF NOT EXISTS meta;
CREATE SCHEMA IF NOT EXISTS derived;
CREATE SCHEMA IF NOT EXISTS analysis;
CREATE SCHEMA IF NOT EXISTS audit;

-- ---------- raw:忠实镜像(含导入元数据) ----------
DROP TABLE IF EXISTS raw.operations;
CREATE TABLE raw.operations (
    op_id              double precision,
    subject_id         double precision,
    hadm_id            double precision,
    case_id            double precision,
    opdate             bigint,
    age                double precision,
    sex                text,
    weight             double precision,
    height             double precision,
    race               text,
    asa                text,
    emop               text,
    department         text,
    antype             text,
    icd10_pcs          text,
    orin_time          bigint,
    orout_time         bigint,
    opstart_time       bigint,
    opend_time         bigint,
    admission_time     bigint,
    discharge_time     bigint,
    anstart_time       bigint,
    anend_time         bigint,
    cpbon_time         bigint,
    cpboff_time        bigint,
    icuin_time         bigint,
    icuout_time        bigint,
    inhosp_death_time  bigint,
    allcause_death_time bigint
);

DROP TABLE IF EXISTS raw.diagnosis;
CREATE TABLE raw.diagnosis (
    subject_id         double precision,
    chart_time         bigint,
    icd10_cm           text
);

DROP TABLE IF EXISTS raw.vitals;
CREATE TABLE raw.vitals (
    op_id              double precision,
    subject_id         double precision,
    chart_time         bigint,
    item_name          text,
    value              double precision
);

DROP TABLE IF EXISTS raw.ward_vitals;
CREATE TABLE raw.ward_vitals (
    subject_id         double precision,
    chart_time         bigint,
    item_name          text,
    value              double precision
);

DROP TABLE IF EXISTS raw.labs;
CREATE TABLE raw.labs (
    subject_id         double precision,
    chart_time         bigint,
    item_name          text,
    value              double precision
);

DROP TABLE IF EXISTS raw.medications;
-- 注意:schema.csv 列出 9 列,但实际文件头为 5 列(subject_id,chart_time,route,drug_name,atc_code)。
-- 以实际文件头为准;不一致已记入 audit.schema_file_mismatch。
CREATE TABLE raw.medications (
    subject_id         double precision,
    chart_time         bigint,
    route              text,
    drug_name          text,
    atc_code           text
);

-- schema.csv 与实际文件头差异台账
DROP TABLE IF EXISTS audit.schema_file_mismatch;
CREATE TABLE audit.schema_file_mismatch (
    object_name text,
    issue       text,
    detail      text,
    recorded_at timestamptz DEFAULT now()
);
INSERT INTO audit.schema_file_mismatch (object_name, issue, detail) VALUES
 ('raw.medications', 'schema.csv vs file header',
  'schema.csv lists subject_id,chart_time,drug_name,drug_name2,drug_name3,route,atc_code,atc_code2,atc_code3 (9 cols); actual file header is subject_id,chart_time,route,drug_name,atc_code (5 cols). Imported per actual header.');

-- ---------- meta:导入与字典 ----------
DROP TABLE IF EXISTS meta.load_log;
CREATE TABLE meta.load_log (
    table_name   text PRIMARY KEY,
    source_file  text,
    source_sha256 text,
    csv_rows     bigint,
    db_rows      bigint,
    match_ok     boolean,
    loaded_at    timestamptz DEFAULT now()
);

-- 字段类型台账(schema.csv 为唯一依据)
DROP TABLE IF EXISTS meta.schema_inventory;
CREATE TABLE meta.schema_inventory (
    table_name  text,
    variable    text,
    schema_type text,
    description text,
    PRIMARY KEY (table_name, variable)
);

DROP TABLE IF EXISTS meta.parameter_dictionary;
CREATE TABLE meta.parameter_dictionary (
    table_name  text,
    label       text,
    unit        text,
    description text,
    PRIMARY KEY (table_name, label)
);

DROP TABLE IF EXISTS meta.department_dictionary;
CREATE TABLE meta.department_dictionary (
    abbrev      text PRIMARY KEY,
    full_name   text
);

COMMIT;
