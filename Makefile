.PHONY: setup test synthetic mimic inspire eicu results figures validate lint clean

setup:
	@echo "Python: $$(python3 --version); R: $$(Rscript --version)"
	@python3 -m pip install -r requirements.txt
	@echo "R packages: 见 renv.lock(如存在)或 README 依赖列表"

test:
	python3 -m pytest tests/pytest -q
	Rscript tests/testthat/run_tests.R

synthetic:
	python3 data/synthetic/make_synthetic.py
	python3 tests/pytest/../..//bin/true 2>/dev/null || true
	Rscript analyses/09_quality_control/99_synthetic_workflow.R

lint:
	@echo "== 敏感信息与绝对路径扫描 ==" && python3 scripts/sensitive_scan.py
	@echo "== YAML/CSV schema 检查 ==" && python3 scripts/check_configs.py

mimic:
	@echo "MIMIC 主分析需要授权数据与 config/paths.yml;顺序:"
	@echo "  1) sql/mimic/00-12 (PostgreSQL)"
	@echo "  2) analyses/02_glucose_processing/01_rebuild_glucose_series.py"
	@echo "  3) analyses/01_cohort_construction/02_build_analysis_dataset.py"
	@echo "  4) analyses/03_primary_mimic/03_run_primary_models_mice.R"
	@echo "  5) analyses/03_primary_mimic/04_ph_rcs_absoluterisk.R"
	bash scripts/run_mimic.sh

inspire:
	bash scripts/run_inspire.sh

eicu:
	bash scripts/run_eicu.sh

results:
	@echo "结果生成顺序见 docs/ANALYSIS_WORKFLOW.md"

figures:
	Rscript analyses/04_time_dependent/02_figures_framefix.R
	Rscript analyses/04_time_dependent/03_figures_fixedscale.R
	Rscript analyses/04_time_dependent/04_rcs_curves.R

validate:
	python3 scripts/validate_results.py

all: lint test synthetic validate
