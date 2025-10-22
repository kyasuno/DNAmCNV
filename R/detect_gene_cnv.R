#' Detect gene-level CNVs
#'
#' @description Detect gene-level (or any local region) CNVs. Adopted from conumee2.0
#' (https://github.com/hovestadtlab/conumee2.git)
#' @param segdata list. Output of \code{run_segmentation}.
#' @param crs list. Output of \code{tangent_normalization}. Necessary only when probe-binning was used (default: NULL).
#' @param genes GenomicRanges object (predefined gene list is available: \code{data(genes.hg38)})
#' @param extend.bp integer. Extend gene interval for \code{+/- extend.bp} (default: \code{500000L}).
#' @param conf numeric. Confidence level to calculate to determine the log2-threshold for
#' high-level alterations. Choose between \code{0.95} and \code{0.99}. Default to \code{0.99}.
#' @param R numeric. Parameter for the \code{bootRanges} function.
#' The number of bootstrap samples to generate. Default to \code{100}.
#' @param blockLength numeric. Parameter for the \code{bootRanges} function. The length
#' (in basepairs) of the blocks for segmented block bootstrapping. Default to \code{500000}.
#' @param proportionLength logical. From the \code{nullranges} package: for the segmented block
#' bootstrap, whether to use scaled block lengths, (scaling by the proportion of the segmentation
#' state out of the total genome length)
#' @returns list. genes, tbl of gene-level result. amp/del_threshold, thresholds for amplification
#' and deletions determined by bootRanges.
#' @export
#'
detect_gene_cnv <- function(segs, crs=NULL, genes, extend.bp=500000, conf = 0.99, R = 100, blockLength = 500000, proportionLength = TRUE) {
  data(hg38.seqinfo)

  seq <- rep(segs$segments$mean.lrr, segs$segments$n.markers)
  km <- kmeans(seq, centers=3, nstart=25)

  chroms <- paste0("chr", c(1:22, "X"))
  bins.df <- segs$probeData |>
    dplyr::mutate(
      seqnames=as.character(seqnames) |> factor(levels=chroms),
      lrr.adjusted=lrr.denoised - segs$cns.shift,
      state=km$cluster
    )

  # This part assumes that all the bins are connected.
  # In our case, we need to merge consecutive bins if the states are the same
  # seg <- lapply(seq_len(3), function(s){ # Combine nearby regions within same states
  #   bins.s <- bins |> dplyr::filter(state == s)
  #   gr.s <- GenomicRanges::GRanges(
  #     seqnames=bins.s$seqnames,
  #     ranges=IRanges::IRanges(start=bins.s$bin.start, bins.s$bin.end)
  #   )
  #   x <- IRanges::reduce(gr.s)
  #   mcols(x)$state <- s
  #   x
  # })
  #
  # seg <- c(seg[[1]], seg[[2]], seg[[3]])


  seg.df <- lapply(chroms, function(chrom) {
    bins.chr <- bins.df |>
      dplyr::filter(seqnames == chrom) |>
      dplyr::mutate(stateChange=c(0, diff(state)))
    idx <- which(bins.chr$stateChange != 0)
    start.idx <- c(1, idx)
    end.idx <- c(idx - 1, nrow(bins.chr))
    bins.chr$group.id <- NA_integer_
    for (i in seq_along(start.idx)) {
      bins.chr$group.id[start.idx[i]:end.idx[i]] <- i
    }
    bins.chr
  }) |> dplyr::bind_rows() |>
    dplyr::group_by(seqnames, group.id, state) |>
    dplyr::summarise(start=min(start), end=max(end), .groups="drop") |>
    arrange(state) # as in CNV.focal

  seg <- GenomicRanges::GRanges(
    seqnames=seg.df$seqnames,
    ranges=IRanges::IRanges(start=seg.df$start, end=seg.df$end),
    state=seg.df$state
  )
  seqinfo(seg) <- hg38.seqinfo[paste0("chr", c(1:22,"X"))]

  bins <- GenomicRanges::GRanges(
    seqnames=bins.df$seqnames,
    ranges=IRanges::IRanges(start=bins.df$start, end=bins.df$end),
    lrr.adjusted=bins.df$lrr.adjusted,
    state=bins.df$state
  )
  seqinfo(bins) <- hg38.seqinfo[paste0("chr", c(1:22,"X"))]

  boots <- nullranges::bootRanges(bins, blockLength = blockLength, R = R,
                                  seg = seg, proportionLength = proportionLength)

  # determine cutoff values for deletions and amplifications
  t <- round(length(boots) * (1-conf), digits = 0)/2
  del_t <- round(sort(boots$lrr.adjusted)[t], digits = 3)
  amp_t <- round(sort(boots$lrr.adjusted, decreasing =  TRUE)[t], digits = 3)

  if (is.null(crs)) {
    probes.gr <- GenomicRanges::GRanges(
      seqnames=segs$probeData$seqnames,
      ranges=IRanges::IRanges(start=segs$probeData$start, end=segs$probeData$end)
    )
    probe_lrr <- segs$probeData$lrr
  } else {
    probes.gr <- GenomicRanges::GRanges(
      seqnames=crs$probeCoords$seqnames,
      ranges=IRanges::IRanges(start=crs$probeCoords$start, end=crs$probeCoords$end)
    )
    probe_lrr <- crs$lrr
  }

  if (extend.bp > 0L) {
    BiocGenerics::start(genes) <- BiocGenerics::start(genes) - extend.bp
    BiocGenerics::end(genes) <- BiocGenerics::end(genes) + extend.bp
  }

  ov <- GenomicRanges::findOverlaps(genes, probes.gr)
  qx <- queryHits(ov)
  sx <- subjectHits(ov)

  genes.lrr <- vapply(split(probe_lrr[sx], qx), median, FUN.VALUE = numeric(1)) - segs$cns.shift
  genes.nprobes <- vapply(split(probe_lrr[sx], qx), length, FUN.VALUE = numeric(1))
  genes.w.data <- genes[unique(qx)]
  mcols(genes.w.data)$n.probes <- genes.nprobes
  mcols(genes.w.data)$lrr <- genes.lrr
  mcols(genes.w.data)$amplified <- genes.lrr >= amp_t
  mcols(genes.w.data)$deleted <- genes.lrr <= del_t

  out <- genes.w.data |> tibble::as_tibble() |>
    dplyr::select(-strand)

  list(genes=out, amp_threshold=amp_t, del_threshold=del_t)

}
