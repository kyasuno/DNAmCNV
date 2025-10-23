#' Annotate segments
#'
#' @param segs list. Output of run_segmentation.
#' @param target_genes GRanges object containing genes of interest (default: NULL).
#' @returns list. The same list with additional columns in segs$segments.
#' @export
#'
annotate_segments <- function(segs, target_genes=NULL) {
  cytoband <- sesameData::sesameData_getGenomeInfo("hg38")$cytoBand |>
    tibble::as_tibble() |>
    dplyr::filter(chrom %in% paste0("chr", c(1:22, "X"))) |>
    dplyr::mutate(
      name=paste0(sub("chr", "", chrom), name)
    ) |>
    dplyr::mutate(
      chrom=as.character(chrom) |> factor(levels=paste0("chr", c(1:22, "X"))),
      chromStart=chromStart+1L
    )
  cyto.gr <- GenomicRanges::GRanges(
    seqnames=cytoband$chrom,
    IRanges::IRanges(start=cytoband$chromStart, end=cytoband$chromEnd),
    name=cytoband$name
  )

  segs.gr <- GenomicRanges::GRanges(
    seqnames=segs$segments$seqnames,
    IRanges::IRanges(start=segs$segments$start, end=segs$segments$end),
    seg.id=segs$segments$seg.id,
    mean.lrr=segs$segments$mean.lrr,
    sizeMb=segs$segments$sizeMb,
    n.markers=segs$segments$n.markers,
    Location=paste0(segs$segments$seqnames, ":", segs$segments$start, "-", segs$segments$end)
  )

  ## annotate cytoband
  message(Sys.time(), ": annotate cytoband")
  ov <- GenomicRanges::findOverlaps(segs.gr, cyto.gr)
  qH <- S4Vectors::queryHits(ov)
  sH <- S4Vectors::subjectHits(ov)
  band.list <- split(x=mcols(cyto.gr)[["name"]][sH], f=qH)
  segs$segments$cytoband <- NA_character_
  segs$segments$cytoband[unique(qH)] <- vapply(band.list, function(blst) {
    if (length(blst) == 1) {
      return(blst)
    } else {
      bs <- blst[1]
      be <- blst[length(blst)]
      paste(bs, be, sep="-")
    }
  }, FUN.VALUE = character(1))

  ## annotate segments with genes of interest
  message(Sys.time(), ": annotate segments with genes of interest")
  if (!is.null(target_genes)) {
    if (!(class(target_genes) == "GRanges")) {
      stop("target_genes must be provided as a GRanges object.")
    }
    if (!any(names(mcols(target_genes)) == "gene_name")) {
      stop("target_genes should contain a column 'gene_name' for gene identifier.")
    }
    ov <- GenomicRanges::findOverlaps(segs.gr, target_genes)
    qH <- S4Vectors::queryHits(ov)
    sH <- S4Vectors::subjectHits(ov)
    gene.list <- split(x=sH, f=qH)
    segs$segments$genesOfInterest <- NA_character_
    segs$segments$genesOfInterest[unique(qH)] <- vapply(gene.list, function(x) {
      if (length(x) == 1) {
        return(mcols(target_genes)$gene_name[x])
      } else {
        paste(mcols(target_genes)$gene_name[x], collapse=",")
      }
    }, FUN.VALUE = character(1))
  }

  return(segs)
}
