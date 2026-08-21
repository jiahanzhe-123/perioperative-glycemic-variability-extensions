#!/usr/bin/env Rscript
# Shared Phase 2B helpers. All source paths come from the existing final-lineage
# config.yaml; this layer never resolves the historical extension branch.

phase2b_arg <- commandArgs(trailingOnly = FALSE)
phase2b_file_arg <- grep("^--file=", phase2b_arg, value = TRUE)
if (length(phase2b_file_arg) == 0L) stop("Phase 2B scripts must be run with Rscript")
phase2b_script_path <- normalizePath(sub("^--file=", "", phase2b_file_arg[1]), mustWork = TRUE)
phase2b_workspace <- normalizePath(file.path(dirname(phase2b_script_path), "../.."), mustWork = TRUE)

source(file.path(phase2b_workspace, "scripts", "phase2a", "phase2a_common.R"))
suppressMessages(library(jsonlite))

phase2b_load_config <- function() phase2a_load_config()
phase2b_output_path <- function(cfg, ...) phase2a_output_path(cfg, ...)
phase2b_read <- function(path, ...) phase2a_read(path, ...)
phase2b_write <- function(data, path) phase2a_write(data, path)
phase2b_num <- function(x) phase2a_num(x)

phase2b_require_gates <- function(cfg) {
  phase1_ok <- phase2a_phase1_gate(cfg)
  gate_path <- phase2b_output_path(cfg, "qc", "phase2a_final_gate.json")
  if (!file.exists(gate_path)) stop("PHASE2A_GATE_MISSING")
  gate <- jsonlite::fromJSON(gate_path)
  phase2_ok <- identical(as.character(gate$decision), "GO_PHASE2B_ANALYTIC_CONTEXT_LANDSCAPE") &&
    isTRUE(gate$hard_failure == FALSE) && isTRUE(gate$modeling_performed == FALSE)
  if (!phase1_ok) stop("PHASE1_6_GATE_FAIL")
  if (!phase2_ok) stop("PHASE2A_GATE_FAIL")
  invisible(TRUE)
}

phase2b_sha256 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  out <- system2("shasum", c("-a", "256", path), stdout = TRUE, stderr = FALSE)
  if (!length(out)) return(NA_character_)
  strsplit(out[1], "[[:space:]]+")[[1]][1]
}

phase2b_origin <- function(path, selector, origin_type = "LOCKED_RESULT") {
  data.frame(
    origin_type = origin_type,
    origin_path = path,
    origin_sha256 = phase2b_sha256(path),
    row_selector = selector,
    stringsAsFactors = FALSE
  )
}

phase2b_scalar <- function(x, default = NA_real_) {
  z <- suppressWarnings(as.numeric(x)[1])
  if (!length(z) || !is.finite(z)) default else z
}

phase2b_fmt <- function(x, digits = 4L) {
  z <- phase2b_scalar(x)
  if (!is.finite(z)) return("NA")
  formatC(z, format = "f", digits = digits)
}
