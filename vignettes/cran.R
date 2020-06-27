# Don't run everything on CRAN (it takes too long)
vignette.eval.full = !as.logical(Sys.getenv("_R_CHECK_LIMIT_CORES_", "FALSE"))
if (!vignette.eval.full) {
	message("Partial vignette evaluation")
}