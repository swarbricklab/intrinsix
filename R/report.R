# Opt-in self-contained HTML report for a pam50_run(). Renders the bundled
# rmarkdown template (inst/rmd/intrinsix_report.Rmd) with the run's pam50_result +
# config passed in via `params`, re-drawing the package's own plot functions so the
# figures match the saved PDFs. Off by default; threaded through pam50_run() and
# `intrinsix run --report`.

# Render the intrinsix HTML report for a pam50_result.
#
# @param res A pam50_result (from pam50_subtype()/pam50_run()).
# @param output_file Path for the HTML file to write.
# @param cfg Optional parsed run config list (for context in the report).
# @param title Report title.
# @return The output path, invisibly; NULL (with a warning) if rmarkdown/knitr or
#   the template are unavailable.
# @noRd
.pam50_report <- function(res, output_file, cfg = NULL, title = "intrinsix PAM50 report") {
  if (!requireNamespace("rmarkdown", quietly = TRUE) ||
      !requireNamespace("knitr", quietly = TRUE)) {
    warning("report = TRUE needs the 'rmarkdown' and 'knitr' packages; skipping report.",
            call. = FALSE)
    return(invisible(NULL))
  }
  rmd_src <- system.file("rmd", "intrinsix_report.Rmd", package = "intrinsix")
  if (!nzchar(rmd_src) || !file.exists(rmd_src)) {
    warning("report template not found (inst/rmd/intrinsix_report.Rmd); skipping report.",
            call. = FALSE)
    return(invisible(NULL))
  }

  output_file <- normalizePath(output_file, winslash = "/", mustWork = FALSE)
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

  # Render from a copy in a temp dir so knitr's intermediate artefacts never land
  # in the installed inst/rmd (which may be read-only) or the output dir.
  work <- tempfile("intrinsix_report_")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  rmd_tmp <- file.path(work, "intrinsix_report.Rmd")
  file.copy(rmd_src, rmd_tmp, overwrite = TRUE)

  rmarkdown::render(
    rmd_tmp,
    output_format = rmarkdown::html_document(self_contained = TRUE, toc = TRUE,
                                             toc_float = TRUE, theme = "readable"),
    output_file = basename(output_file),
    output_dir = dirname(output_file),
    params = list(result = res, cfg = cfg, title = title),
    quiet = TRUE,
    envir = new.env(parent = globalenv())
  )
  invisible(output_file)
}
