#' Annotate segments
#'
#' @param segs list. Output of run_segmentation.
#' @returns list. The same list with additional columns in segs$segments.
#' @export
#'
annotate_segments <- function(segs) {
  cytoband <- sesameData::sesameData_getGenomeInfo("hg38")$cytoBand |>
    tibble::as_tibble() |>
    dplyr::filter(chrom %in% paste0("chr", c(1:22, "X"))) |>
    dplyr::mutate(
      chrom=as.character(chrom) |> factor(levels=paste0("chr", c(1:22, "X"))),
      chromStart=chromStart+1L
    )
  cyto.gr <- GenomicRanges::GRanges(
    seqnames=cytoband$chrom,
    IRanges::IRanges(start=cytoband$chromStart, end=cytoband$chromEnd),
    gieStain=cytoband$gieStain,
  )
  acen.gr <- reduce(cyto.gr[cytoband$gieStain == "acen"])

  segs.gr <- GenomicRanges::GRanges(
    seqnames=segs$segments$seqnames,
    IRanges::IRanges(start=segs$segments$start, end=segs$segments$end),
  )

  ov <- findOverlaps(segs.gr, acen.gr, type="within")
  if (length(ov) > 0) {
    qx <- S4Vectors::queryHits(ov)
    segs$segments <- segs$segments |>
      dplyr::mutate(withinCentromere=1:n() %in% qx)
  } else {
    segs$segments$withinCentromere <- FALSE
  }
  return(segs)
}
