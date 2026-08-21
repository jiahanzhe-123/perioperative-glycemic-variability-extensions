# PGV 配置访问(R)。
# 规则:
# - 外部数据键(患者级/授权目录)必须来自 config/paths.yml(本机,gitignored);
#   模板见 config/paths.example.yml。未配置时明确报错,绝不回退到仓库内目录。
# - 仓库内部键(synthetic/home)有仓库相对默认值,无需配置。
pgv_repo_root <- function(start = getwd()) {
  d <- normalizePath(start, mustWork = FALSE)
  repeat {
    if (file.exists(file.path(d, "config", "paths.example.yml"))) return(d)
    nd <- dirname(d)
    if (nd == d) stop("repository root (config/paths.example.yml) not found above ", start)
    d <- nd
  }
}
.pgv_read_yml <- function(path) {
  out <- list()
  if (!file.exists(path)) return(out)
  for (ln in readLines(path, warn = FALSE)) {
    ln <- sub("#.*$", "", ln)
    if (!grepl("^[A-Za-z0-9_]+\\s*:", ln)) next
    key <- sub("\\s*:.*$", "", ln)
    val <- sub("^[A-Za-z0-9_]+\\s*:\\s*", "", ln)
    val <- gsub('^["\']|["\']$', "", trimws(val))
    out[[key]] <- val
  }
  out
}
.pgv_cfg <- local({
  cfg <- .pgv_read_yml(file.path(pgv_repo_root(), "config", "paths.yml"))
  cfg <- Filter(function(v) nzchar(v) && !startsWith(v, "/path/to/"), cfg)
  lapply(cfg, path.expand)
})
.pgv_internal <- list(synthetic = file.path("data", "synthetic"), home = ".")
PGV <- function(key) {
  v <- .pgv_cfg[[key]]
  if (!is.null(v) && nzchar(v)) return(normalizePath(v, mustWork = FALSE))
  if (key %in% names(.pgv_internal))
    return(file.path(pgv_repo_root(), .pgv_internal[[key]]))
  stop("PGV: key '", key, "' is not configured. External data paths must be set ",
       "in config/paths.yml (copy config/paths.example.yml and edit). ",
       "Formal code never falls back to in-repo data directories.")
}
PGV_OUT <- function(sub = "results") file.path(pgv_repo_root(), sub)
PGV_SEED <- function() 20260726L
