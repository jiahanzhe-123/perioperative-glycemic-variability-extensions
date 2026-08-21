#!/usr/bin/env Rscript

# Public Phase 3B builder.
#
# This script intentionally consumes aggregate rows only. It does not access
# controlled data, patient-level source values, bootstrap replicates, or local
# configuration, and it does not refit a model.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (!length(file_arg)) stop("Run with Rscript so the repository can be resolved")
script_path <- normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)
repo <- normalizePath(file.path(dirname(script_path), "../../.."), mustWork = TRUE)
ext <- file.path(repo, "analysis_extensions")
input <- file.path(ext, "results", "phase3b", "figure_3_source_data.csv")
out_dir <- file.path(ext, "figures", "phase3b_public")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

landscape <- read.csv(input, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("context_id", "database", "context_family", "effect_type", "N", "events",
              "estimate", "lower95", "upper95", "figure_panel", "figure_order",
              "same_patient_pair")
missing <- setdiff(required, names(landscape))
if (length(missing)) stop("Missing aggregate columns: ", paste(missing, collapse = ", "))
if (nrow(landscape) != 18L) stop("Expected 18 locked aggregate rows; found ", nrow(landscape))
if (any(grepl("stay_id|subject_id|hadm_id|patient_id|bootstrap", names(landscape), ignore.case = TRUE))) {
  stop("Patient-level or replicate-like column detected in public landscape input")
}

short_labels <- c(
  MIMIC_ADJUSTMENT_A_30D = "Model A",
  MIMIC_ADJUSTMENT_B_30D = "Model B",
  MIMIC_ADJUSTMENT_C_30D = "Model C",
  MIMIC_SOURCE_PRIORITY_30D = "Priority series",
  MIMIC_SOURCE_POCT_ONLY_30D = "POCT-only",
  MIMIC_SOURCE_CENTRAL_LAB_ONLY_30D = "Central-lab-only",
  MIMIC_SOURCE_BLOOD_GAS_ONLY_30D = "Blood-gas-only",
  MIMIC_SOURCE_COMMON_30D = "Common-source",
  MIMIC_SOURCE_SAME_PATIENT_POCT_30D = "Same-patient POCT",
  MIMIC_SOURCE_SAME_PATIENT_LAB_30D = "Same-patient laboratory",
  INSPIRE_I1_24H_30D = "I1",
  INSPIRE_I2_24H_30D = "I2",
  INSPIRE_I3_24H_30D = "I3",
  INSPIRE_I2_48H_CORRECTED_30D = "Corrected 48-hour",
  EICU_M1_24H_AGGREGATE = "M1",
  EICU_M2_24H_AGGREGATE = "M2",
  EICU_M3_24H_AGGREGATE = "M3",
  EICU_M4_24H_AGGREGATE = "M4"
)
if (any(!landscape$context_id %in% names(short_labels))) {
  stop("Unmapped context_id in public landscape input")
}
landscape$label <- paste0(unname(short_labels[landscape$context_id]), "\nN=",
                          format(landscape$N, big.mark = ","), "/", landscape$events)
landscape$label <- factor(landscape$label, levels = rev(landscape$label[order(landscape$figure_order)]))
landscape$colour <- ifelse(landscape$same_patient_pair, "firebrick", "steelblue4")
landscape$shape <- ifelse(landscape$same_patient_pair, 17, 16)

theme_public <- theme_classic(base_size = 8) +
  theme(
    plot.title = element_text(face = "bold", size = 9),
    axis.text.y = element_text(size = 5.7, colour = "black"),
    axis.text.x = element_text(size = 6.5, colour = "black"),
    axis.title = element_text(size = 7),
    panel.grid = element_blank(),
    plot.margin = margin(4, 20, 4, 4, unit = "pt")
  )

make_panel <- function(data, title, effect_label) {
  data <- data[order(data$figure_order), , drop = FALSE]
  data$label <- factor(data$label, levels = rev(data$label))
  ggplot(data, aes(x = estimate, y = label)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey25", linewidth = 0.4) +
    geom_segment(aes(x = lower95, xend = upper95, y = label, yend = label, colour = colour), linewidth = 0.65) +
    geom_point(aes(colour = colour, shape = shape), size = 2) +
    scale_colour_identity() +
    scale_shape_identity() +
    scale_x_log10(limits = c(0.45, 1.65), breaks = c(0.5, 0.75, 1, 1.25, 1.5)) +
    labs(title = title, x = paste0(effect_label, " (log scale; reference = 1)"), y = NULL) +
    theme_public
}

panels <- split(landscape, landscape$figure_panel)
p1 <- make_panel(panels[["MIMIC_adjustment"]], "MIMIC adjustment context", "HR")
p2 <- make_panel(panels[["MIMIC_source"]], "MIMIC measurement-source context", "HR")
p3 <- make_panel(panels[["INSPIRE_timing_adjustment"]], "INSPIRE timing/adjustment context", "HR")
p4 <- make_panel(panels[["eICU_adjustment"]], "eICU adjustment context", "RR")
figure <- (p1 | p2) / (p3 | p4) +
  plot_annotation(title = "Public aggregate analytic-context landscape",
                  subtitle = "Conceptual order; no pooling or common summary effect") &
  theme(plot.title = element_text(face = "bold", size = 11, hjust = 0))

ggsave(file.path(out_dir, "Figure_3_public.png"), figure, width = 183, height = 142, units = "mm", dpi = 300)
ggsave(file.path(out_dir, "Figure_3_public.pdf"), figure, width = 183, height = 142, units = "mm", device = grDevices::pdf)
svglite::svglite(file.path(out_dir, "Figure_3_public.svg"), width = 183 / 25.4, height = 142 / 25.4)
print(figure)
dev.off()

writeLines(c(
  "Public builder completed from 18 aggregate SD-GV rows.",
  "No model was refit; no patient-level source values were read.",
  "MIMIC/INSPIRE rows are HRs and eICU rows are RRs; results are not pooled."
), file.path(out_dir, "BUILD_NOTE.txt"))
message("PHASE3B_PUBLIC_AGGREGATE_FIGURE_DONE")
