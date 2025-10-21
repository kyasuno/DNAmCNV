#' Smooth segments
#'
#' @description Merge consecutive segments in the same copy number state (loss, copy-neutral or gain).
#' @param segs list. Output of run_segmentation or annotate_segments.
#' @param loss.lrr numeric. Cutoff value of lrr for losses (default: -0.1). If NULL, segments are not classified.
#' @param gain.lrr numeric. Cutoff value of lrr for gains (default: 0.1). If NULL, segments are not classified.
#' @param use.pcf logical. Use PCF method within each class (default: FALSE).
#' @param gamma numeric. penalty for each discontinuity in the curve (default: 40).
#' @returns list. The same list with additional tbl named smoothed.
#' @export
#'
smooth_segments <- function(segs, loss.lrr=-0.1, gain.lrr=0.1, use.pcf=FALSE, gamma=40) {
  newsegs <- segs$segments |>
    dplyr::mutate(
      CN=dplyr::if_else(mean.lrr < loss.lrr, "LOSS", dplyr::if_else(mean.lrr > gain.lrr, "GAIN", "CN")),
      chrArm=paste0(seqnames, arm)
    )
  res <- lapply(split(newsegs, newsegs$chrArm), function(df) {
    if (nrow(df) == 0) { return(NULL) }
    if (nrow(df) == 1) {
      return(
        df |>
          dplyr::mutate(nMerged=1L) |>
          dplyr::select(seqnames, arm, start, end, sizeMb, n.markers, mean.lrr, CN, nMerged)
      )
    }
    out <- tibble::tibble()
    ## the first segment
    CN <- df$CN[1]
    chrom <- df$seqnames[1]
    arm <- df$arm[1]
    start <- df$start[1]
    end <- df$end[1]
    #sizeMb <- df$sizeMb[1]
    n.mk <- df$n.markers[1]
    lrr <- df$mean.lrr[1]
    for (i in 2:nrow(df)) {
      CN.i <- df$CN[i]
      if (CN.i == CN) {
        # merge segments
        # arm <- c(arm, df$arm[i]) # no need to update arm anymore as we analyze each arm
        # end <- df$end[i]
        start <- c(start, df$start[i])
        end <- c(end, df$end[i])
        #sizeMb <- c(sizeMb, df$sizeMb[i])
        n.mk <- c(n.mk, df$n.markers[i])
        lrr <- c(lrr, df$mean.lrr[i])
        if (i == nrow(df)) {
          # we must stop here
          # find optimal smoothing
          if (use.pcf) {
            tmp <- find_optimal_partition(chrom=chrom, arm=arm, start=start, end=end,
                                          lrr=lrr, n.mk=n.mk, CN=CN, gamma=gamma)
          } else {
            tmp <- tibble::tibble(
              seqnames=chrom,
              arm=arm,
              start=min(start),
              end=max(end),
              sizeMb=NA_real_,
              n.markers=sum(n.mk),
              mean.lrr=weighted.mean(lrr, n.mk),
              CN=CN,
              nMerged=length(start)
            ) |>
              dplyr::mutate(
                sizeMb=(end - start + 1)/1e6
              )
          }

          out <- dplyr::bind_rows(out, tmp)
        } # otherwise, continue
      } else {
        # record the previous segment
        if (use.pcf) {
          tmp <- find_optimal_partition(chrom=chrom, arm=arm, start=start, end=end,
                                        lrr=lrr, n.mk=n.mk, CN=CN, gamma=gamma)
        } else {
          tmp <- tibble::tibble(
            seqnames=chrom,
            arm=arm,
            start=min(start),
            end=max(end),
            sizeMb=NA_real_,
            n.markers=sum(n.mk),
            mean.lrr=weighted.mean(lrr, n.mk),
            CN=CN,
            nMerged=length(start)
          ) |>
            dplyr::mutate(
              sizeMb=(end - start + 1)/1e6
            )
        }

        out <- dplyr::bind_rows(out, tmp)
        if (i < nrow(df)) {
          # reset values
          CN <- df$CN[i]
          arm <- df$arm[i]
          start <- df$start[i]
          end <- df$end[i]
          n.mk <- df$n.markers[i]
          lrr <- df$mean.lrr[i]
        } else {
          # record the last segment
          tmp <- df[i, ] |>
            dplyr::mutate(nMerged=1L) |>
            dplyr::select(seqnames, arm, start, end, sizeMb, n.markers, mean.lrr, CN, nMerged)
          out <- dplyr::bind_rows(out, tmp)
        }
      }
    }
    out
  }) |>
    dplyr::bind_rows() |> dplyr::arrange(seqnames, arm)

  if (any(names(segs) == "smoothed")) {
    message("Previous smoothing result will be removed")
    segs$smoothed <- NULL
  }

  segs <- c(segs, list(smoothed=res))
  segs
}

find_optimal_partition <- function(chrom, arm, start, end, lrr, n.mk, CN, gamma) {
  n <- length(lrr)
  if (n == 1) {
    return(
      tibble::tibble(
        seqnames=chrom,
        arm=arm,
        start=start, #min(start),
        end=end, #max(end),
        sizeMb=NA_real_,
        n.markers=n.mk, #sum(n.mk),
        mean.lrr=lrr, #weighted.mean(lrr, n.mk),
        CN=CN,
        nMerged=1L #length(x0:x1)
      ) |>
        dplyr::mutate(
          sizeMb=(end - start + 1)/1e6
        )
    )
  }
  m <- partitions::compositions(n)
  # each column represents the number of elements of each partition
  mx <- apply(m, 2, function(x) x[x > 0], simplify = FALSE)
  mx <- lapply(mx, function(x) {
    # get start-end indices from the set sizes
    starts <- c(1, 1+cumsum(x)[1:(length(x)-1)])
    ends <- cumsum(x)
    size <- ends - starts + 1
    keep <- size > 0
    cbind(start=starts[keep], end=ends[keep])
  })
  loss <- lapply(mx, function(mm) {
    S <- nrow(mm)
    nI <- sum(n.mk)
    z <- apply(mm, 1, function(idx) {
      (sum(lrr[idx]*n.mk[idx]))^2
    })
    - sum(z)/nI + gamma*S
  }) |> unlist()
  p <- which.min(loss)
  mopt <- mx[[p]]
  nsegments <- nrow(mopt)
  if (nsegments == 1) {
    out <- tibble::tibble(
      seqnames=chrom,
      arm=arm,
      start=min(start), end=max(end), sizeMb=NA_real_,
      n.markers=sum(n.mk),
      mean.lrr=weighted.mean(lrr, n.mk),
      CN=CN,
      nMerged=length(n.mk)
    ) |>
      dplyr::mutate(
        sizeMb=(end - start + 1)/1e6
      )
  } else {
    out <- purrr::map_dfr(seq(nsegments), function(i) {
      x0 <- mopt[i, 1]
      x1 <- mopt[i, 2]
      tibble::tibble(
        seqnames=chrom,
        arm=arm,
        start=min(start[x0:x1]), end=max(end[x0:x1]), sizeMb=NA_real_,
        n.markers=sum(n.mk[x0:x1]),
        mean.lrr=weighted.mean(lrr[x0:x1], n.mk[x0:x1]),
        CN=CN,
        nMerged=length(x0:x1)
      ) |>
        dplyr::mutate(
          sizeMb=(end - start + 1)/1e6
        )
    })
  }
  out
}
