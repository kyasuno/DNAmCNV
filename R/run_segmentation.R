#' Run segmentation using copynumber
#'
#' @param probeData tbl. seqnames (chromosomes chr1, ..., chrX), start, end, chrArm (chr1p, ..., chrXq), lrr (log R ratio).
#' @param denoise.method character. none, winsorize or rmf (recursive median filter).
#' @param kmin integer. minimum number of probes (or bins) in each segment used by copynumber::pcf (default: 5).
#' @param gamma numeric. A penalty parameter used by copynumber::pcf for each discontinuity in the curve (default: 40).
#' @param normalize logical. Whether the copy number measurements should be scaled by the sample residual standard error (default: TRUE).
#' When denoise.method = "winsorize", normalize = TRUE is recommended while when denoise.method = "rmf", normalize = FALSE is recommended.
#' @param k integer. window size to be used for the sliding window (actually half-window size) for copynumber::winsorize or denoise_intensities.
#' @param tau numeric. Winsorization threshold, default is 1.5.
#' @param adjust.baseline logical. Shift factor is calculated using limma::weighted.median of lrr among potentially
#' copy neutral segments (abs(lrr) < 0.1) (default: FALSE).
#' @param use.n.probes logical. If TRUE, the number of probes (bins) will be used as weights in the baseline adjusment.
#' If FALSE, the physical length of segments will be used as weights. (Default: TRUE)
#' @param n.cores integer. number of cores to be used (default: 1), n.cores > 1 is effective only for \code{denoise.method = "rmv"}.
#' @param verbose logical.
#' @returns list.
#' @export
#'
run_segmentation <- function(probeData, denoise.method=c("none", "winsorize", "rmf"),
                             kmin=5, gamma=40, normalize=TRUE, k=25, tau=1.5, adjust.baseline=FALSE,
                             use.n.probes=TRUE,
                             n.cores=1L, verbose=FALSE) {
  denoise.method <- match.arg(denoise.method)
  probeData <- probeData |>
    dplyr::mutate(
      chrom=sub("chr", "", as.character(seqnames)),
      arms=dplyr::if_else(grepl("p$", chrArm), "p", "q"),
      position=end # use end position
    )
  # denoise signals
  if (denoise.method == "none") {
    probeData$lrr.denoised <- round(probeData$lrr, digits=6)
  } else if (denoise.method == "winsorize") {
    probeData$lrr.denoised <- copynumber::winsorize(
      data=probeData |> dplyr::select(chrom, position, lrr) |> as.data.frame(),
      arms=probeData$arms, k=k, tau=tau, method="mad", assembly="hg38", digits=6,
      verbose=FALSE
    )$lrr
  } else if (denoise.method == "rmf") {
    probeData$lrr.denoised <- denoise_intensities(
      probeData$lrr, probeData$chrArm, k=k, n.cores=n.cores, verbose=FALSE
    )
  }

  # pcf
  segs <- copynumber::pcf(data=probeData |>
                            dplyr::select(chrom, position, lrr.denoised) |>
                            as.data.frame(),
                          arms=probeData$arms, assembly="hg38", digits=6,
                          normalize=normalize, kmin=kmin, gamma=gamma, verbose=FALSE) |>
    tibble::as_tibble() |>
    dplyr::rename(seqnames=chrom) |>
    dplyr::select(-sampleID) |>
    dplyr::mutate(
      seqnames=paste0("chr", seqnames) |> factor(levels=paste0("chr", c(1:22, "X")))
    )
  # for binned data, adjust start coordinates
  if (any(probeData$end != probeData$start)) {
    tmp <- probeData |> dplyr::select(seqnames, position, start) |> dplyr::rename(start.pos=position)
    segs <- left_join(segs, tmp, by=c("seqnames", "start.pos")) |>
      dplyr::rename(end=end.pos)
  } else {
    segs <- segs |> dplyr::rename(start=start.pos, end=end.pos)
  }
  segs <- segs |>
    dplyr::mutate(sizeMb=as.numeric(end - start + 1) / 1e6) |>
    dplyr::select(seqnames, arm, start, end, sizeMb, n.probes, mean) |>
    dplyr::rename(n.markers=n.probes)

  # check any bias in lrr for copy-neutral segment
  cns <- segs |> dplyr::filter(abs(mean) < 0.1)
  if (use.n.probes) {
    cns.median <- with(cns, limma::weighted.median(mean, n.markers))
  } else {
    cns.median <- with(cns, limma::weighted.median(mean, sizeMb))
  }


  # calculate Z score for lrr
  segs.gr <- GenomicRanges::GRanges(
    seqnames=segs$seqnames,
    IRanges::IRanges(start=segs$start, end=segs$end)
  )
  probe.gr <- GenomicRanges::GRanges(
    seqnames=probeData$seqnames,
    IRanges::IRanges(start=probeData$start, end=probeData$end)
  )
  X <- GenomicRanges::findOverlaps(segs.gr, probe.gr)
  qx <- S4Vectors::queryHits(X)
  sx <- S4Vectors::subjectHits(X)
  sx.list <- split(x=sx, f=qx)
  lrr.stats <- map_dfr(sx.list, function(idx) {
    x <- probeData$lrr.denoised[idx]
    tibble::tibble(
      mean.lrr=mean(x - cns.median, na.rm=TRUE),
      sd.lrr=sd(x - cns.median, na.rm=TRUE) #mad(x, na.rm=TRUE) / qnorm(3/4)
    ) |>
      dplyr::mutate(
        Z.lrr=dplyr::if_else(sd.lrr > 0, mean.lrr / sd.lrr, NA_real_)
      )
  })
  segs <- dplyr::bind_cols(segs, lrr.stats) |>
    dplyr::mutate(seg.id=1L:n())

  probeData$seg.id <- NA_integer_
  for (i in seq_along(sx.list)) {
    probeData$seg.id[sx.list[[i]]] <- i
  }
  probeData <- probeData |>
    dplyr::select(-chrom, -arms, -position)

  return(
    list(probeData=probeData, segments=segs, cns.shift=cns.median)
  )
}
