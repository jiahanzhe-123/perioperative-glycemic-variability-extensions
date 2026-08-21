#!/usr/bin/env Rscript
# Shared Phase 2A helpers. Authoritative data paths are resolved from the
# Phase 1.6 config.yaml; no historical extension path is used here.

phase2a_script_arg <- commandArgs(trailingOnly = FALSE)
phase2a_file_arg <- grep("^--file=", phase2a_script_arg, value = TRUE)
if (length(phase2a_file_arg) == 0L) stop("Phase 2A scripts must be run with Rscript")
phase2a_script_path <- normalizePath(sub("^--file=", "", phase2a_file_arg[1]), mustWork = TRUE)
phase2a_workspace <- normalizePath(file.path(dirname(phase2a_script_path), "../.."), mustWork = TRUE)

phase2a_parse_value <- function(x) {
  x <- trimws(x)
  if ((startsWith(x, "\"") && endsWith(x, "\"")) ||
      (startsWith(x, "'") && endsWith(x, "'"))) return(substr(x, 2L, nchar(x) - 1L))
  low <- tolower(x)
  if (low %in% c("true", "yes")) return(TRUE)
  if (low %in% c("false", "no")) return(FALSE)
  if (low %in% c("null", "none", "na", "~")) return(NA)
  if (grepl("^[-+]?[0-9]+$", x)) return(as.integer(x))
  if (grepl("^[-+]?(?:[0-9]+\\.[0-9]*|[0-9]*\\.[0-9]+)(?:[eE][-+]?[0-9]+)?$", x, perl = TRUE)) return(as.numeric(x))
  x
}

phase2a_load_config <- function() {
  cfg_path <- file.path(phase2a_workspace, "config.yaml")
  lines <- readLines(cfg_path, warn = FALSE, encoding = "UTF-8")
  cfg <- list()
  for (raw in lines) {
    line <- trimws(sub("#.*$", "", raw))
    if (!nzchar(line)) next
    parts <- strsplit(line, ":", fixed = TRUE)[[1]]
    if (length(parts) < 2L) stop("Malformed flat config line: ", raw)
    key <- trimws(parts[1])
    value <- paste(parts[-1], collapse = ":")
    cfg[[key]] <- phase2a_parse_value(value)
  }
  required <- c("public_code_root", "analysis_record_root", "controlled_package_root",
                "manuscript_provenance_root", "output_root", "random_seed")
  missing <- required[!vapply(required, function(k) !is.null(cfg[[k]]), logical(1))]
  if (length(missing)) stop("Config missing keys: ", paste(missing, collapse = ", "))
  cfg$workspace_root <- phase2a_workspace
  cfg
}

phase2a_cfg_path <- function(cfg, key) {
  value <- path.expand(as.character(cfg[[key]]))
  if (grepl("^/", value)) return(normalizePath(value, mustWork = FALSE))
  normalizePath(file.path(cfg$workspace_root, value), mustWork = FALSE)
}

phase2a_output_path <- function(cfg, ...) {
  dir.create(phase2a_cfg_path(cfg, "output_root"), recursive = TRUE, showWarnings = FALSE)
  normalizePath(file.path(phase2a_cfg_path(cfg, "output_root"), ...), mustWork = FALSE)
}

phase2a_read <- function(path, ...) {
  if (!file.exists(path)) stop("Missing declared input: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, ...)
}

phase2a_require_columns <- function(data, columns, label) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) stop(label, " missing columns: ", paste(missing, collapse = ", "))
}

phase2a_truth <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "1", "yes", "y", "t")
}

phase2a_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

phase2a_write <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(data, path, row.names = FALSE, na = "NA", quote = TRUE)
}

phase2a_open_log <- function(cfg, name) {
  path <- file.path(cfg$workspace_root, "logs", name)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  sink(path, split = TRUE)
  path
}

phase2a_close_log <- function() {
  while (sink.number() > 0L) sink()
}

phase2a_phase1_gate <- function(cfg) {
  gate <- jsonlite::fromJSON(file.path(cfg$workspace_root, "outputs", "qc", "final_gate.json"))
  identical(as.character(gate$decision), "GO_PHASE2_MEASUREMENT_CONTEXT_ANALYSES") && isTRUE(gate$modeling_performed == FALSE)
}

phase2a_target_ids <- function(cfg) {
  context <- phase2a_read(phase2a_output_path(cfg, "machine_readable", "final_measurement_context.csv"))
  context$patient_or_episode_id[context$database == "MIMIC-IV" & context$analysis_context == "MIMIC_primary_24h_day1_landmark"]
}

phase2a_source_data <- function(cfg) {
  root <- phase2a_cfg_path(cfg, "analysis_record_root")
  base <- phase2a_read(phase2a_cfg_path(cfg, "mimic_analysis_base"))
  feat <- phase2a_read(phase2a_cfg_path(cfg, "mimic_features_priority"))
  sp <- phase2a_read(phase2a_cfg_path(cfg, "mimic_samepatient_source"))
  base$stay_id <- as.character(base$stay_id)
  feat$stay_id <- as.character(feat$stay_id)
  sp$stay_id <- as.character(sp$stay_id)
  base$landmark_eligible <- phase2a_truth(base$landmark_eligible)
  list(root = root, base = base, feat = feat, sp = sp)
}
