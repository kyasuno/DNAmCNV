#' Make a panel of normals from SigDF data
#'
#' @param sdfs list. A named list of SigDF data.
#' @param platform character. HM450, EPIC, or EPICv2.
#' @param use.mask logical. Use SeSAMe mask to filter signals (default: TRUE).
#' @param max.p.missing.sample numeric. Maximum proportion of missing data per sample (default: 0.2).
#' @param percentile.var numeric. Remove top `percentile.var` percentile of highest variance (default: 1 percentile).
#' @param n.cores integer. Number of cores to be used to process data.
#' @param median.normalize logical. If TRUE, divide the intensity by their median (default: TRUE).
#' @param exclude.acen logical. If TRUE, probes within gieStain = "acen" will be removed (default: FALSE).
#' @param exclude.regions GenomicRanges object. Regions to be removed from the analysis (default: NULL).
#' A predefined set is available: data(exclude_regions.hg38).
#' @returns list. probeCoords, GenomicRanges object for probe coordinates. pred.sex, predicted sex.
#' data, a matrix of total intensities that were divided by sample median and log2 transformed.
#' @export
#'
make_pon <- function(sdfs, platform=c("HM450", "EPIC", "EPICv2"),
                     use.mask=TRUE,
                     max.p.missing.sample=0.2, percentile.var=1, n.cores=1L, median.normalize=TRUE,
                     exclude.acen=FALSE, exclude.regions=NULL) {
  platform <- match.arg(platform)
  # Note: 0-based coordinates
  genome <- sesameData::sesameData_check_genome(NULL, platform)
  genomeInfo <- sesameData::sesameData_getGenomeInfo(genome)

  # process cytoband information and remove p arms of acrocentric chromosomes
  data(chrom_arm.hg38)

  probeCoords <- sesameData::sesameData_getManifestGRanges(platform, genome = genome)
  x <- IRanges::mergeByOverlaps(probeCoords, chrom_arm.hg38)
  probeCoords <- x[["probeCoords"]]
  probeCoords$chrArm <- x[["chrArm"]]
  rm(x)

  message("Building a panel of normals for ", platform, " and ", genome, " genome.")

  # remove probes in PAR region
  if (genome == "hg38") {
    message("Removing probes in PAR region")
    parx <- GenomicRanges::GRanges(
      seqnames=c("chrX", "chrX"),
      IRanges::IRanges(start=c(10001, 155701383)-1, end=c(2781479, 156030895))
    )
    probeCoords <- probeCoords[!IRanges::overlapsAny(probeCoords, parx)]
  }

  # remove probes on chrY and chrM
  message("Removing probes in chrY or chrM")
  remove <- as.character(GenomeInfoDb::seqnames(probeCoords)) %in% c("chrY", "chrM")
  probeCoords <- probeCoords[!remove]

  # remove probes in gieStain = "acen"
  if (exclude.acen) {
    data(acen.hg38)
    message("Removing probes in gieStain = acen")
    probeCoords <- probeCoords[!IRanges::overlapsAny(probeCoords, acen.hg38)]

  }
  # remove probes in exclude.region
  if (!is.null(exclude.regions)) {
    message("Removing probes in exclue.regions")
    probeCoords <- probeCoords[!IRanges::overlapsAny(probeCoords, exclude.regions)]
  }

  # keep only CpG probes (as in conumee)
  message("Keep only CpG probes")
  keep <- grepl("^cg", names(probeCoords))
  probeCoords <- probeCoords[keep]

  # infer sex
  message("Predicting sex of samples")
  pred.sex <- parallel::mclapply(sdfs, function(sdf) sesame::inferSex(sdf), mc.cores=n.cores) |> unlist()

  # calculate M+U
  MU <- parallel::mclapply(sdfs, function(sdf) {
    mu <- sesame::totalIntensities(sdf, mask=FALSE)
    if (use.mask) {
      mu[sdf$mask] <- NA_real_
    }
    mu
  }, mc.cores=n.cores)
  MU <- do.call(cbind, MU)
  MU <- MU[names(probeCoords), , drop=FALSE]

  # remove samples with high missingness
  rmv <- apply(MU, 2, function(x) sum(is.na(x))/length(x)) >= max.p.missing.sample
  if (sum(rmv) > 0) {
    MU <- MU[, !rmv, drop=FALSE]
    message(sum(rmv), " samples are excluded for high missingness (>= ", max.p.missing.sample, ")")
    pred.sex <- pred.sex[!rmv]
  }

  nF <- sum(pred.sex == "FEMALE")
  nM <- sum(pred.sex == "MALE")
  message("Building a panel of normals consisting of ", nF, " females and ", nM, " males.")

  # remove probes with one or more missing data
  remove <- apply(MU, 1, function(x) sum(is.na(x))/length(x)) > 0
  message("Removing ", sum(remove), " probes with missing data")
  probeCoords <- probeCoords[!remove]
  MU <- MU[!remove, , drop=FALSE]

  # remove probes in the top k percentile of highest variance
  W <- matrixStats::rowVars(log2(MU), na.rm=TRUE)
  keep <- W <= quantile(W, probs=(100 - percentile.var)/100)

  message("Removing ", sum(!keep), " highly variable probes with top ", percentile.var, " percentile of its variance")
  probeCoords <- probeCoords[keep]
  MU <- MU[keep, , drop=FALSE]

  if (median.normalize) {
    # Divide MU by sample median
    sample.med <- matrixStats::colMedians(MU, na.rm=TRUE)
    MU <- t(t(MU)/sample.med)
  }
  # Then log2 transform using safeLog2
  epsilon <- 1e-9
  log2epsilon <- log2(epsilon)
  MU <- apply(MU, 2, function(x) dplyr::if_else(x < epsilon, log2epsilon, log2(x), missing=NA_real_))
  rownames(MU) <- names(probeCoords)

  probenames <- names(probeCoords)

  probeCoords <- probeCoords |>
    as.data.frame() |>
    tibble::as_tibble() |>
    dplyr::mutate(Probe_ID=probenames, start=start + 1) |>
    dplyr::select(Probe_ID, seqnames, start, end, strand, chrArm)

  message("Building a panel of normals is completed with ", nrow(MU), " probes and ", ncol(MU), " samples.")

  return(
    list(probeCoords=probeCoords, pred.sex=pred.sex, data=MU, median.normalized=median.normalize)
  )
}
