#' Merge nearby probes
#'
#' @param crs list. Output of `tangent_normalization`.
#' @param min.gapwidth integer. For `GenomicRanges::reduce` function (default: 500L).
#' @param max.probes integer. Maximum number of probes to be merged within a set of probes. Probes are split into nearly equally
#' `ceiling(k/max.probes)` groups, where k is the number of probes in the group.
#' @param n.cores integer. Number of cores to be used to run this function (default: 1).
#' @returns tbl. seqnames, start, end, width, nProbes, chrArm and median log R ratio (lrr).
#' @export
#'
merge_nearby_probes <- function(crs, min.gapwidth, max.probes, n.cores=1L) {

  df <- crs$probeCoords |>
    dplyr::mutate(
      width=end - start + 1,
      nProbes=1L,
      lrr=crs$lrr
    )
  gr <- GenomicRanges::GRanges(
    seqnames=df$seqnames,
    IRanges::IRanges(start=df$start, end=df$end)
  )
  gr.red <- GenomicRanges::reduce(gr, min.gapwidth=min.gapwidth, with.revmap=TRUE)

  binsize <- lapply(gr.red$revmap, length) |> unlist()

  # For binsize = 1, no need to merge
  idx1 <- gr.red$revmap[binsize == 1] |> unlist()
  df1 <- df[idx1, ] |>
    dplyr::select(seqnames, chrArm, start, end, width, nProbes, lrr)

  # analyze others
  revmap <- gr.red$revmap[binsize > 1]

  merged <- parallel::mclapply(revmap, function(x) {
        n <- length(x)
        if (n <= max.probes) {
          return(
            tibble::tibble(
              seqnames=crs$probeCoords$seqnames[x[1]],
              chrArm=crs$probeCoords$chrArm[x[1]],
              start=crs$probeCoords$start[x[1]] ,
              end=crs$probeCoords$end[x[n]],
              width=end - start + 1,
              nProbes=n,
              lrr=median(crs$lrr[x])
            )
          )
        }
        out <- purrr::map_dfr(split(x, ceiling(seq(n)/max.probes)), function(idx) {
          k <- length(idx)
          tibble::tibble(
            seqnames=crs$probeCoords$seqnames[idx[1]],
            chrArm=crs$probeCoords$chrArm[idx[1]],
            start=crs$probeCoords$start[idx[1]],
            end=crs$probeCoords$end[idx[k]],
            width=end - start + 1,
            nProbes=length(idx),
            lrr=median(crs$lrr[idx])
          )
        })
        out
      }, mc.cores=n.cores) |> dplyr::bind_rows()

  out <- dplyr::bind_rows(df1, merged) |>
    dplyr::arrange(seqnames, start) |>
    dplyr::mutate(
      start=as.integer(start),
      end=as.integer(end),
      width=as.integer(width),
      lrr=as.numeric(lrr)
    )

  return(
    #list(merged=out, shift=shift)
    out
  )
}
