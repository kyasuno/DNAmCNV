#' Get probe bins
#'
#' @param crs list. Output of `tangent_normalization`.
#' @param bin.size integer. In base pairs (default: 1000L).
#' @param hg38.seqinfo Use `data("hg38.seqinfo")`.
#' @param centromere.hg38 Use `data("centromere.hg38")`.
#' @returns tbl. seqnames, start, end, width, nProbes, chrArm and median log R ratio (lrr).
#' @export
#'
get_probe_bins <- function(crs, bin.size=1000L, hg38.seqinfo, centromere.hg38) {

  tiled <- GenomicRanges::tileGenome(hg38.seqinfo,
                                     tilewidth = bin.size, cut.last.tile.in.chrom = TRUE)
  keep <- as.character(seqnames(tiled)) %in% paste0("chr", c(1:22, "X"))
  tiled <- tiled[keep]

  # exclude centromeres
  centro.gr <- GenomicRanges::GRanges(
    seqnames=centromere.hg38$chrom,
    IRanges::IRanges(start=centromere.hg38$chromStart, end=centromere.hg38$chromEnd)
  )

  qx <- queryHits(findOverlaps(tiled, centro.gr))
  tiled <- tiled[-qx]

  gr <- GenomicRanges::GRanges(
    seqnames=crs$probeCoords$seqnames,
    IRanges::IRanges(start=crs$probeCoords$start, end=crs$probeCoords$end)
  )

  mcols(tiled)$nProbes <- countOverlaps(tiled, gr)
  keep <- tiled$nProbes > 0
  tiled <- tiled[keep]

  o1 <- findOverlaps(query = tiled, subject = gr)
  bins <- queryHits(o1)
  prbs <- subjectHits(o1)

  chrarm <- as.character(crs$probeCoords$chrArm)

  mcols(tiled)$lrr <- vapply(split(crs$lrr[prbs], bins), median, FUN.VALUE=numeric(1))
  mcols(tiled)$chrArm <- vapply(split(chrarm[prbs], bins), FUN=function(x) x[1], FUN.VALUE=character(1))
  mcols(tiled)$chrArm <- factor(mcols(tiled)$chrArm, levels=levels(crs$probeCoords$chrArm))
  mcols(tiled)$true.start <- vapply(split(crs$probeCoords$start[prbs], bins), dplyr::first, FUN.VALUE=double(1))
  mcols(tiled)$true.end <- vapply(split(crs$probeCoords$end[prbs], bins), dplyr::last, FUN.VALUE=double(1))


  df <- tiled |> as_tibble() |> dplyr::select(-strand) |>
    dplyr::mutate(width=true.end - true.start + 1L) |>
    dplyr::rename(bin.start=start, bin.end=end) |>
    dplyr::rename(start=true.start, end=true.end) |>
    dplyr::select(seqnames, start, end, width, nProbes, chrArm, lrr, bin.start, bin.end) |>
    dplyr::mutate(start=as.integer(start), end=as.integer(end))
  return(df)
}
