#' Denoise intensities using median filters recursively or non-recursively
#'
#' @param x numeric vector.
#' @param chrArm a factor of chromosome arms of the same size as x.
#' @param k window size to be used for the sliding window (actually half-window size).
#' @param recursive logical. Perform median filters recursively or non-recursively (default: FALSE).
#' @param verbose logical.
#' @returns numeric vector.
#' @export
#'
denoise_intensities <- function(x, chrArm, k, n.cores, recursive=FALSE, verbose=FALSE) {
  arms <- chrArm %>% unique()

  if (recursive) {
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
  } else {
    res <- parallel::mclapply(arms, function(ax) {
      keep <- chrArm == ax
      runmed(x[keep], k=2*k + 1, endrule="median")
    }, mc.cores=n.cores) |> unlist()
  }

  res
}

