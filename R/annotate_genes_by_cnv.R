#' Annotate genes by CNV status
#'
#' @description This function can take long time to complete for a large number of genes.
#' @param segs list. Output of run_segmentation.
#' @param genes GRanges object containing genes or transcripts. If `NULL`,
#' the program calls `sesameData::sesameData_getTxnGRanges("hg38", merge2gene=TRUE)`.
#' @param n.cores integer. Specify how many cores to be used for computation (default: 1).
#' @returns list. An expanded list that includes 'genes' annotated with overlapping segments.
#' @export
#'
annotate_genes_by_cnv <- function(segs, genes=NULL, n.cores=1) {

  segs.gr <- GenomicRanges::GRanges(
    seqnames=segs$segments$seqnames,
    IRanges::IRanges(start=segs$segments$start, end=segs$segments$end),
    seg.id=segs$segments$seg.id,
    mean.lrr=segs$segments$mean.lrr,
    sizeMb=segs$segments$sizeMb,
    n.markers=segs$segments$n.markers,
    Location=paste0(segs$segments$seqnames, ":", segs$segments$start, "-", segs$segments$end)
  )

  ## annotate genes with segment id and mean.lrr
  message(Sys.time(), ": annotate genes by segments")

  if (is.null(genes)) {
    genes <- sesameData::sesameData_getTxnGRanges("hg38", merge2gene=TRUE)
    keep <- as.character(GenomeInfoDb::seqnames(genes)) %in% paste0("chr", c(1:22, "X"))
    genes <- genes[keep]
  }


  ov <- GenomicRanges::findOverlaps(genes, segs.gr)
  qH <- S4Vectors::queryHits(ov)
  sH <- S4Vectors::subjectHits(ov)
  seg.list.per.gene <- split(x=sH, f=qH)

  message(length(seg.list.per.gene), " genes to be annotated")

  if (n.cores == 1) {
    ann <- lapply(seg.list.per.gene, function(lst) {
      tibble::tibble(
        seg.id=paste(mcols(segs.gr)$seg.id[lst], collapse=","),
        mean.lrr=if (length(lst) == 1) mcols(segs.gr)$mean.lrr[lst] else NA_real_,
        sizeMb=if (length(lst) == 1) mcols(segs.gr)$sizeMb[lst] else NA_real_,
        n.markers=if (length(lst) == 1) mcols(segs.gr)$n.markers[lst] else NA_real_,
        Location=if (length(lst) == 1) mcols(segs.gr)$Location[lst] else NA_character_,
      )
    }) |> dplyr::bind_rows()

  } else {
    ann <- parallel::mclapply(seg.list.per.gene, function(lst) {
      tibble::tibble(
        seg.id=paste(mcols(segs.gr)$seg.id[lst], collapse=","),
        mean.lrr=if (length(lst) == 1) mcols(segs.gr)$mean.lrr[lst] else NA_real_,
        sizeMb=if (length(lst) == 1) mcols(segs.gr)$sizeMb[lst] else NA_real_,
        n.markers=if (length(lst) == 1) mcols(segs.gr)$n.markers[lst] else NA_real_,
        Location=if (length(lst) == 1) mcols(segs.gr)$Location[lst] else NA_character_,
      )
    }, mc.cores=n.cores) |> dplyr::bind_rows()
  }

  mcols(genes)$seg.id <- NA_character_
  mcols(genes)$seg.id[unique(qH)] <- ann$seg.id
  mcols(genes)$mean.lrr <- NA_real_
  mcols(genes)$mean.lrr[unique(qH)] <- ann$mean.lrr
  mcols(genes)$sizeMb <- NA_real_
  mcols(genes)$sizeMb[unique(qH)] <- ann$sizeMb
  mcols(genes)$n.markers <- NA_real_
  mcols(genes)$n.markers[unique(qH)] <- ann$n.markers
  mcols(genes)$Location <- NA_character_
  mcols(genes)$Location[unique(qH)] <- ann$Location

  genes <- genes |> tibble::as_tibble()

  segs <- c(segs, list(genes=genes))
  message(Sys.time(), ": done.")

  return(segs)
}
