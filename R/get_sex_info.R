#' get probe data for sex inference
#'
#' @param sdf SigDF object.
#' @platform character. HM450, EPIC, or EPICv2.
#' @param cutoff a cutoff value of log2(medianY/medianX) for sex inference (default: -3.5).
#' @returns a named list with pred.sex (predicted sex), medianX (median intensity of X probes), and
#' medianY (median intensity of Y probes).
#'
#' @export
#'
get_sex_info <- function(sdf, platform, cutoff=-3.5) {
  stopifnot(is(sdf, "SigDF"))
  if (platform == "EPICv2") {
    pids <- lapply(str_split(sdf$Probe_ID, "\\_"), function(x) x[1]) |> unlist() |> unique()
    platform <- "EPIC"
  } else {
    pids <- sdf$Probe_ID
  }

  cleanY <- sesameDataGet(paste0(
    platform,'.probeInfo'))$chrY.clean
  xLinked <- sesameDataGet(paste0(
    platform,'.probeInfo'))$chrX.xlinked

  medY <- median(totalIntensities(sdf[pids %in% cleanY, ]))
  medX <- median(totalIntensities(sdf[pids %in% xLinked, ]))
  sex <-  if (log2(medY / medX) < cutoff) "FEMALE" else "MALE"
  list(pred.sex=sex, medianX=medX, medianY=medY)
}
