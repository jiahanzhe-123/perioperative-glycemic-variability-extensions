#!/usr/bin/env bash
# 01_import_tables.sh — INSPIRE v1.4.2 导入(gzip | psql COPY FROM STDIN)
# 用法:bash 01_import_tables.sh(需要 docker 容器 inspire-pg 运行中)
set -euo pipefail
D="PGV("inspire_raw")"
PSQL="docker exec -i inspire-pg psql -U postgres -d inspire_v142 -v ON_ERROR_STOP=1 -q"

echo "[$(date '+%F %T')] 创建 schemas 与 DDL"
$PSQL < "$(dirname "$0")/../sql/00_create_database.sql"

import_gz () {
  local file="$1" table="$2"
  local csv_rows
  csv_rows=$(( $(gunzip -c "$D/$file" | wc -l | tr -d ' ') - 1 ))
  echo "[$(date '+%F %T')] 导入 $table($file,CSV 数据行 $csv_rows)"
  gunzip -c "$D/$file" | $PSQL -c "COPY raw.$table FROM STDIN WITH (FORMAT csv, HEADER true, NULL '')"
}

import_gz operations.csv.gz   operations
import_gz diagnosis.csv.gz    diagnosis
import_gz labs.csv.gz         labs
import_gz medications.csv.gz  medications
import_gz vitals.csv.gz       vitals
import_gz ward_vitals.csv.gz  ward_vitals

echo "[$(date '+%F %T')] 导入 meta 字典"
$PSQL -c "\COPY meta.schema_inventory (table_name, variable, schema_type, description) FROM STDIN WITH (FORMAT csv, HEADER true, NULL '')" < \
  <(python3 -c "
import csv
rows=list(csv.reader(open('$D/schema.csv', newline='')))
rows=[r for r in rows if r and any(x.strip() for x in r)]
cur=None
for r in rows[1:]:
    if len(r)<4: continue
    t,var,typ,desc=r[0].strip(),r[1].strip(),r[2].strip(),r[3].strip()
    if t: cur=t
    print(f'{cur},{var},{typ},\"{desc.replace(chr(34), chr(34)*2)}\"')
")
$PSQL -c "\COPY meta.parameter_dictionary (table_name, label, unit, description) FROM STDIN WITH (FORMAT csv, HEADER true, NULL '')" < "$D/parameters.csv"
$PSQL -c "\COPY meta.department_dictionary (abbrev, full_name) FROM STDIN WITH (FORMAT csv, HEADER true, NULL '')" < "$D/department.csv"

echo "[$(date '+%F %T')] 计数核对"
for t in operations diagnosis labs medications vitals ward_vitals; do
  csv_rows=$(( $(gunzip -c "$D/$t.csv.gz" | wc -l | tr -d ' ') - 1 ))
  db_rows=$($PSQL -Atc "SELECT count(*) FROM raw.$t")
  ok=$([ "$csv_rows" = "$db_rows" ] && echo true || echo false)
  echo "$t: csv=$csv_rows db=$db_rows match=$ok"
  sha=$(shasum -a 256 "$D/$t.csv.gz" | awk '{print $1}')
  $PSQL -c "INSERT INTO meta.load_log (table_name, source_file, source_sha256, csv_rows, db_rows, match_ok)
            VALUES ('$t', '$t.csv.gz', '$sha', $csv_rows, $db_rows, $ok)
            ON CONFLICT (table_name) DO UPDATE SET csv_rows=$csv_rows, db_rows=$db_rows, match_ok=$ok, source_sha256='$sha', loaded_at=now()"
done
echo "[$(date '+%F %T')] 导入完成"
