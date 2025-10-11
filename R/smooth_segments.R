#' Smooth segments
#'
#' @description Merge consecutive segments in the same copy number state (loss, copy-neutral or gain).
#' @param segs list. Output of run_segmentation or annotate_segments.
#' @param loss.lrr numeric. Cutoff value of lrr for losses (default: -0.1). If NULL, segments are not classified.
#' @param gain.lrr numeric. Cutoff value of lrr for gains (default: 0.1). If NULL, segments are not classified.
#' @returns list. The same list with additional tbl named smoothed.
#' @export
#'
smooth_segments <- function(segs, loss.lrr=-0.1, gain.lrr=0.1) {
  newsegs <- segs$segments |>
    dplyr::mutate(
      CN=dplyr::if_else(mean.lrr < loss.lrr, "LOSS", dplyr::if_else(mean.lrr > gain.lrr, "GAIN", "CN")),
      chrArm=paste0(seqnames, arm)
    )
  res <- lapply(split(newsegs, newsegs$chrArm), function(df) {
    if (nrow(df) == 0) { return(NULL) }
    if (nrow(df) == 1) {
      return(
        df |>
          dplyr::mutate(nMerged=1L) |>
          dplyr::select(seqnames, arm, start, end, sizeMb, n.markers, mean.lrr, CN, nMerged)
      )
    }
    out <- tibble::tibble()
    ## the first segment
    CN <- df$CN[1]
    chrom <- df$seqnames[1]
    arm <- df$arm[1]
    start <- df$start[1]
    end <- df$end[1]
    #sizeMb <- df$sizeMb[1]
    n.mk <- df$n.markers[1]
    lrr <- df$mean.lrr[1]
    for (i in 2:nrow(df)) {
      CN.i <- df$CN[i]
      if (CN.i == CN) {
        # merge segments
        arm <- c(arm, df$arm[i])
        end <- df$end[i]
        #sizeMb <- c(sizeMb, df$sizeMb[i])
        n.mk <- c(n.mk, df$n.markers[i])
        lrr <- c(lrr, df$mean.lrr[i])
        if (i == nrow(df)) {
          # we must stop here
          tmp <- tibble::tibble(
            seqnames=chrom,
            arm=paste(unique(arm), collapse=""),
            start=start, end=end, sizeMb=(end - start + 1)/1e6,
            n.markers=sum(n.mk),
            mean.lrr=weighted.mean(lrr, n.mk),
            CN=CN,
            nMerged=length(n.mk)
          )
          out <- dplyr::bind_rows(out, tmp)
        } # otherwise, continue
      } else {
        # record the previous segment
        tmp <- tibble::tibble(
          seqnames=chrom,
          arm=paste(unique(arm), collapse=""),
          start=start, end=end, sizeMb=(end - start + 1)/1e6,
          n.markers=sum(n.mk),
          mean.lrr=weighted.mean(lrr, n.mk),
          CN=CN,
          nMerged=length(n.mk)
        )
        out <- dplyr::bind_rows(out, tmp)
        if (i < nrow(df)) {
          # reset values
          CN <- df$CN[i]
          arm <- df$arm[i]
          start <- df$start[i]
          end <- df$end[i]
          n.mk <- df$n.markers[i]
          lrr <- df$mean.lrr[i]
        } else {
          # record the last segment
          tmp <- df[i, ] |>
            dplyr::mutate(nMerged=1L) |>
            dplyr::select(seqnames, arm, start, end, sizeMb, n.markers, mean.lrr, CN, nMerged)
          out <- dplyr::bind_rows(out, tmp)
        }
      }
    }
    out
  }) |>
    dplyr::bind_rows() |> dplyr::arrange(seqnames, arm)
  segs <- c(segs, list(smoothed=res))
  segs
}
