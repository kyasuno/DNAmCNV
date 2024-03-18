#' Make a panel of normals from SigDF data
#'
#' @param sdfs list. A named list of SigDF data.
#' @param max.p.missing.sample numeric. Maximum proportion of missing data per sample (default: 0.2).
#' @param max.p.missing.probe numeric. Maximum proportion of missing data per probe (default: 0.01).
#' @param percentile.var numeric. Remove top `percentile.var` percentile of highest variance (default: 1 percentile).
#' @param n.cores integer. Number of cores to be used to process data.
#' @returns list. probeCoords, GenomicRanges object for probe coordinates. pred.sex, predicted sex.
#' data, a matrix of total intensities that were divided by sample median and log2 transformed.
#' @export
#'
make_pon <- function(sdfs, max.p.missing.sample=0.2, max.p.missing.probe=0.01, percentile.var=1, n.cores=1L) {
  # Note: 0-based coordinates
  genome <- sesameData::sesameData_check_genome(NULL, platform)
  probeCoords <- sesameData::sesameData_getManifestGRanges(platform, genome = genome)
  # remove probes in PAR region
  if (genome == "hg38") {
    parx <- GenomicRanges::GRanges(
      seqnames=c("chrX", "chrX"),
      IRanges::IRanges(start=c(10001, 155701383)-1, end=c(2781479, 156030895))
    )
    probeCoords <- probeCoords[!IRanges::overlapsAny(probeCoords, parx)]
  }

  # remove probes on chrY and chrM
  remove <- seqnames(probeCoords) %in% c("chrY", "chrM")
  probeCoords <- probeCoords[!remove]

  # infer sex
  pred.sex <- parallel::mclapply(sdfs, function(sdf) sesame::inferSex(sdf), mc.cores=n.cores) |> unlist()

  # calculate M+U
  MU <- parallel::mclapply(sdfs, function(sdf) {
    mu <- sesame::totalIntensities(sdf, mask=FALSE)
    mu[sdf$mask] <- NA_real_
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

  # remove probes with high missingness
  keep <- apply(MU, 1, function(x) sum(is.na(x))/length(x)) < max.p.missing.probe
  message("Removing ", sum(!keep), " probes with high missingness (>= ", max.p.missing.probe, ")")
  probeCoords <- probeCoords[keep]
  MU <- MU[keep, , drop=FALSE]

  # remove probes in the top k percentile of highest variance
  W <- rowVars(log2(MU), na.rm=TRUE)
  keep <- W < quantile(W, probs=(100 - percentile.var)/100)
  message("Removing ", sum(!keep), " highly variable probes with top ", percentile.var, " percentile of its variance")
  probeCoords <- probeCoords[keep]
  MU <- MU[keep, , drop=FALSE]

  # Divide MU by sample median and then log2 transform using safeLog2
  sample.med <- colMedians(MU)
  MU <- t(t(MU)/sample.med)
  log2epsilon <- log2(1e-9)
  MU <- apply(MU, 2, function(x) dplyr::if_else(x < epsilon, log2epsilon, log2(x), missing=NA_real_))

  message("Building a panel of normals is completed with ", nrow(MU), " probes and ", ncol(MU), " samples.")

  return(
    list(probeCoords=probeCoords, pred.sex=pred.sex, data=MU)
  )
}
