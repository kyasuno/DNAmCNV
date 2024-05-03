#' Tangent normalization of probe intensities
#'
#' @param sdf SigDF object.
#' @param sex character. Sex of the sample. If NULL, sex is estimated using `sesame::inferSex(sdf)`.
#' @param pon list. panel of normal data.
#' @returns list. probeCoords, probe coordinates. obs, observed signal (`sesame::totalIntensities(sdf)`). pred, predicted (fitted) signal.
#' lrr, log R ratio (`obs - pred`). shift, shift factor of the baseline signal.
#' @export
#'
tangent_normalization <- function(sdf, sex=NULL, pon) {
  # infer sex if not given
  if (is.null(sex)) {
    sex <- sesame::inferSex(sdf)
    message("This sample is predicted ", sex)
  } else {
    sex <- toupper(sex)
  }
  # calculate total intensities
  mu <- sesame::totalIntensities(sdf, mask=FALSE)

  # extract probes both in sdf and pon
  keep <- !sdf$mask & !is.na(mu) & (names(mu) %in% pon$probeCoords$Probe_ID)
  mu <- mu[keep]

  remove <- !(pon$probeCoords$Probe_ID %in% names(mu))
  pon$probeCoords <- pon$probeCoords[!remove, ]
  pon$data <- pon$data[!remove, ]

  # sort by coordinates
  mu <- mu[pon$probeCoords$Probe_ID]

  # normalize tumor signal
  if (pon$median.normalized) {
    med.mu <- median(mu, na.rm=TRUE)
    mu <- mu / med.mu
  }
  epsilon <- 1e-9
  log2epsilon <- log2(epsilon)
  mu <- dplyr::if_else(mu < epsilon, log2epsilon, log2(mu))

  # check if target sample is included in PoN
  rr <- cor(mu, pon$data, use="pairwise.complete.obs")[1, ] < 0.99
  if (any(!rr)) {
    message("target sample seems to be in the PoN. The corresponding sample is removed from PoN, ",
            "resulting in ", sum(rr), " PoN samples.")
    pon$data <- pon$data[, rr, drop=FALSE]
    pon$pred.sex <- pon$pred.sex[rr]
  }

  ## split data by autosomes vs chrX
  chrX <- pon$probeCoords$seqnames == "chrX"
  fita <- lm(y ~ ., data=data.frame(y=mu[!chrX], pon$data[!chrX, ]), na.action=na.omit)
  fitx <- lm(y ~ ., data=data.frame(y=mu[chrX], pon$data[chrX, pon$pred.sex == sex]))

  pred <- c(predict(fita), predict(fitx))

  if (!identical(names(pred), names(mu))) {
    stop("There could be missing data in PoN.")
  }
  lrr <- mu - pred
  names(lrr) <- names(mu)

  # shift <- optim(0, function(s) median(abs(lrr - s), na.rm = TRUE),
  #                method = "Brent", lower = -100, upper = 100)$par

  return(
    list(probeCoords=pon$probeCoords, observed=mu, predicted=pred, lrr=lrr) #, shift=shift)
  )
}
