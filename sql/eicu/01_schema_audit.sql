-- Step 1: schema audit
\pset footer off
\o /tmp/audit_tables.csv
SELECT table_name FROM information_schema.tables WHERE table_schema='public' ORDER BY table_name;
\o /tmp/audit_columns.csv
SELECT table_name, ordinal_position, column_name, data_type, character_maximum_length
FROM information_schema.columns
WHERE table_schema='public'
ORDER BY table_name, ordinal_position;
\o
