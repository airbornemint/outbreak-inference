# Don't run everything on CRAN (it takes too long)
cran = as.logical(Sys.getenv("_R_CHECK_LIMIT_CORES_", "FALSE"))
vignette.eval = list(full=!cran, limitCores=ifelse(cran, 2, NA))
if (!vignette.eval$full) {
	message("Partial vignette evaluation")
}
