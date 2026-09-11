# ==============================================================================
# run_all.R
# Runs the full clinical ML pipeline end to end, knitting every notebook in
# strict numerical order (01 -> 04). Each notebook depends on objects
# serialized by the previous one (data/split, data/processed, data/models),
# so they must be executed in sequence, never skipped or reordered.
#
# Usage:
#   - RStudio: open this file and click "Source" (or Ctrl/Cmd + Shift + S)
#   - Terminal: Rscript run_all.R
# ==============================================================================

# Anchor the working directory to the project root regardless of where the
# script is launched from, using the .Rproj file as the anchor point.
if (!requireNamespace("rprojroot", quietly = TRUE)) {
  install.packages("rprojroot")
}
root <- rprojroot::find_root(rprojroot::has_file("clinical-ml-pipeline.Rproj"))
setwd(root)

# Restore the exact package versions locked in renv.lock before running
# anything, so the pipeline behaves identically to the original development
# environment.
if (requireNamespace("renv", quietly = TRUE)) {
  renv::restore(prompt = FALSE)
}

# Notebooks to knit, in strict execution order
notebooks <- c(
  "notebooks/01_EDA_and_preprocessing.Rmd",
  "notebooks/02_unsupervised_modeling.Rmd",
  "notebooks/03_supervised_modeling.Rmd",
  "notebooks/04_clinical_synthesis.Rmd"
)

# Knit each notebook in a fresh environment, stopping immediately if one
# fails so a broken step never silently propagates into the next notebook.
for (nb in notebooks) {
  message("== Rendering ", nb, " ==")
  rmarkdown::render(nb, envir = new.env())
  message("== Done: ", nb, " ==\n")
}

message("Pipeline finished successfully. All four notebooks were knitted. See /notebooks for the HTML outputs.")
