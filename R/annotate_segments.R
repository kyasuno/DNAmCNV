#' Annotate segments
#'
#' @param segs list. Output of run_segmentation.
#' @param target_genes GRanges object containing genes of interest (default: NULL).
#' @param annotate_genes logical. If TRUE, annotate genes with segment id and mean.lrr.
#' @param genes_for_annotation GRanges object containing genes/transcripts. If `NULL`,
#' the program calls `sesameData::sesameData_getTxnGRanges("hg38", merge2gene=TRUE)`.
#' @returns list. The same list with additional columns in segs$segments or an expanded list
#' that includes genes annotated with overlapping segments.
#' @export
#'
annotate_segments <- function(segs, target_genes=NULL, annotate_genes=TRUE, genes_for_annotation=NULL) {
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
  ov <- GenomicRanges::findOverlaps(segs.gr, cyto.gr)
  qH <- S4Vectors::queryHits(ov)
  sH <- S4Vectors::subjectHits(ov)
  band.list <- split(x=sH, f=qH)
  segs$segments$cytoband <- NA_character_
  segs$segments$cytoband[unique(qH)] <- lapply(band.list, function(blst) {
    if (length(blst) == 1) {
      return(blst)
    } else {
      bs <- blst[1]
      be <- blst[length(blst)]
      paste(bs, be, sep="-")
    }
  })

  ## annotate segments with genes of interest
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
    segs$segments$genesOfInterest[unique(qH)] <- lapply(gene.list, function(x) {
      if (length(x) == 1) {
        return(x)
      } else {
        paste(mcols(target_genes)$gene_name[x], collapse=",")
      }
    })
  }

  ## annotate genes with segment id and mean.lrr
  if (!annotate_genes) {
    return(segs)
  }

  if (!is.null(genes_for_annotation)) {
    genes <- genes_for_annotation
  } else {
    genes <- sesameData::sesameData_getTxnGRanges("hg38", merge2gene=TRUE)
  }

  keep <- seqnames(genes) %in% paste0("chr", c(1:22, "X"))
  genes <- genes[keep]

  ov <- GenomicRanges::findOverlaps(genes, segs.gr)
  qH <- S4Vectors::queryHits(ov)
  sH <- S4Vectors::subjectHits(ov)
  seg.list.per.gene <- split(x=sH, f=qH)
  mcols(genes)$seg.id <- NA_character_ # when boundaries of segments overlap the gene, they must be reported
  mcols(genes)$seg.id[unique(qH)] <- vapply(seg.list.per.gene, function(lst) {
    paste(mcols(segs.gr)$seg.id[lst], collapse=",")
  }, FUN.VALUE = character(1))
  mcols(genes)$mean.lrr <- NA_real_
  mcols(genes)$mean.lrr[unique(qH)] <- vapply(seg.list.per.gene, function(lst) {
    if (length(lst) > 1) {
      return(NA_real_)
    }
    mcols(segs.gr)$mean.lrr[lst]
  }, FUN.VALUE = numeric(1))
  mcols(genes)$sizeMb <- NA_real_
  mcols(genes)$sizeMb[unique(qH)] <- vapply(seg.list.per.gene, function(lst) {
    if (length(lst) > 1) {
      return(NA_real_)
    }
    mcols(segs.gr)$sizeMb[lst]
  }, FUN.VALUE = numeric(1))
  mcols(genes)$n.markers <- NA_real_
  mcols(genes)$n.markers[unique(qH)] <- vapply(seg.list.per.gene, function(lst) {
    if (length(lst) > 1) {
      return(NA_real_)
    }
    mcols(segs.gr)$n.markers[lst]
  }, FUN.VALUE = numeric(1))
  mcols(genes)$Location <- NA_character_
  mcols(genes)$Location[unique(qH)] <- vapply(seg.list.per.gene, function(lst) {
    if (length(lst) > 1) {
      return(NA_character_)
    }
    mcols(segs.gr)$Location[lst]
  }, FUN.VALUE = character(1))

  genes <- genes |> tibble::as_tibble()

  segs <- c(segs, list(genes=genes))

  return(segs)
}
