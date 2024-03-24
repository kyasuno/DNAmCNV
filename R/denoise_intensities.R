#' Denoise intensities using median filters recursively
#'
#' @param x numeric vector.
#' @param chrArm a factor of chromosome arms of the same size as x.
#' @param k window size to be used for the sliding window (actually half-window size).
#' @param verbose logical.
#' @returns numeric vector.
#' @export
#'
denoise_intensities <- function(x, chrArm, k, n.cores, verbose=FALSE) {
  arms <- chrArm %>% unique()
  resi <- x
  for (i in 1:k) {
    if (verbose) {
      message(Sys.time(), ": step ", i)
    }

    res <- parallel::mclapply(arms, function(ax) {
      keep <- chrArm == ax
      runmed(resi[keep], k=2*i + 1, endrule="median")
    }, mc.cores=n.cores) |> unlist()
    resi <- res
  }
  res
}

