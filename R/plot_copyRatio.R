#' @title Plot copy ratio data
#'
#' @param segs segment data (seqnames, start, end, mean.lrr, Z.lrr).
#' @param crs denoised copy ratio data (seqnames, start, end, lrr).
#' @param smoothed.segs smoothed segments.
#' @param hg38.seqinfo Use `data("hg38.seqinfo")`.
#' @param centromere.hg38 Use `data("centromere.hg38")`.
#' @param ymax Maximum copy ratio (default: 3).
#' @param space.skip Space between chromosomes (default: 0.01)
#' @param loss.lrr numeric. Cutoff value of lrr for losses (default: -0.1). If NULL, segments are not classified.
#' @param gain.lrr numeric. Cutoff value of lrr for gains (default: 0.1). If NULL, segments are not classified.
#' @returns A ggplot object
#' @import ggplot2
#' @export
#'
plot_copyRatio <- function(segs, crs=NULL, smoothed.segs=NULL,
                           hg38.seqinfo, centromere.hg38, ymax=1,
                           space.skip=0.01, loss.lrr=NULL, gain.lrr=NULL) {

  cols <- RColorBrewer::brewer.pal(8,"Dark2")
  col.p <- RColorBrewer::brewer.pal(8,"Paired")
  thickness <- 0.005

  hg38.seqinfo <- hg38.seqinfo[paste0("chr", c(1:22,"X"))]

  # centromere data
  centromere.hg38 <- centromere.hg38 |>
    dplyr::filter(chrom %in% GenomeInfoDb::seqlevels(hg38.seqinfo)) |>
    dplyr::rename(seqnames=chrom, start=chromStart, end=chromEnd) |>
    dplyr::mutate(start=start + 1)
  c.genome <- get_genomic_coord(centromere.hg38, hg38.seqinfo, space.skip=space.skip)
  centromere.hg38 <- c.genome$data

  p <- ggplot()

  ## probe intensities

  if (!is.null(crs)) {
    crs <- crs |> dplyr::filter(seqnames %in% GenomeInfoDb::seqlevels(hg38.seqinfo))
    crs.genome <- get_genomic_coord(crs, hg38.seqinfo, space.skip=space.skip)
    crs <- crs.genome$data
    crs <- crs |>
      dplyr::mutate(
        lrr=dplyr::if_else(lrr > ymax, ymax - thickness,
                           dplyr::if_else(lrr < - ymax, -ymax + thickness, lrr))
      )
    p <- p +
      geom_point(
        data=crs, aes(x=.start, y=lrr),
        colour=col.p[3], size=1/ggplot2::.pt, alpha=0.2,
        inherit.aes=FALSE
      )
  }

  ## segments

  segs <- segs |> dplyr::filter(seqnames %in% GenomeInfoDb::seqlevels(hg38.seqinfo))
  segs.genome <- get_genomic_coord(segs, hg38.seqinfo, space.skip=space.skip)
  segs <- segs.genome$data |>
    dplyr::mutate(
      lrr=dplyr::if_else(mean.lrr > ymax, ymax - thickness,
                         dplyr::if_else(mean.lrr < - ymax, -ymax + thickness, mean.lrr))
    )

  anyLoss <- if (!is.null(loss.lrr)) any(segs$mean.lrr < loss.lrr) else FALSE
  anyGain <- if (!is.null(gain.lrr)) any(segs$mean.lrr > gain.lrr) else FALSE

  if (anyLoss | anyGain) {
    p <- p +
      geom_hline(yintercept=0, linewidth=0.5, colour="gray") +
      geom_rect(
        data=segs |> dplyr::filter(mean.lrr <= gain.lrr & mean.lrr >= loss.lrr),
        aes(xmin=.start, xmax=.end, ymin=lrr - thickness, ymax=lrr + thickness),
        colour="#E6AB02", fill="#E6AB02", na.rm=TRUE, inherit.aes=FALSE
      ) +
      geom_rect(
        data=segs |> dplyr::filter(mean.lrr > gain.lrr | mean.lrr < loss.lrr),
        aes(xmin=.start, xmax=.end, ymin=lrr - thickness, ymax=lrr + thickness),
        colour="#A6761D", fill="#A6761D", na.rm=TRUE, inherit.aes=FALSE
      )
  } else {
    p <- p +
      geom_hline(yintercept=0, linewidth=0.5, colour="gray") +
      geom_rect(
        data=segs,
        aes(xmin=.start, xmax=.end, ymin=lrr - thickness, ymax=lrr + thickness),
        colour="#E6AB02", fill="#E6AB02", na.rm=TRUE, inherit.aes=FALSE
      )
  }

  if (!is.null(smoothed.segs)) {
    smoothed.segs <- smoothed.segs |> dplyr::filter(seqnames %in% GenomeInfoDb::seqlevels(hg38.seqinfo))
    smoothed.segs.genome <- get_genomic_coord(smoothed.segs, hg38.seqinfo, space.skip=space.skip)
    smoothed.segs <- smoothed.segs.genome$data |>
      dplyr::mutate(
        lrr=dplyr::if_else(mean.lrr > ymax, ymax - thickness,
                           dplyr::if_else(mean.lrr < - ymax, -ymax + thickness, mean.lrr))
      )
    p <- p + geom_rect(
      data=smoothed.segs |> dplyr::filter(CN == "CN"),
      aes(xmin=.start, xmax=.end, ymin=lrr - thickness, ymax=lrr + thickness),
      colour="#CCCCCC", fill="#CCCCCC", na.rm=TRUE, inherit.aes=FALSE
    ) +
      geom_rect(
        data=smoothed.segs |> dplyr::filter(CN != "CN"),
        aes(xmin=.start, xmax=.end, ymin=lrr - thickness, ymax=lrr + thickness),
        colour="#000000", fill="#000000", na.rm=TRUE, inherit.aes=FALSE
      )
  }

  p <- p +
    scale_y_continuous(limits=c(-ymax,ymax), breaks=scales::breaks_pretty(n=7), expand=c(0,0)) +
    labs(y="R", x="")

  # centromere
  p <- p +
    scale_x_continuous(breaks = (centromere.hg38$.start + centromere.hg38$.end) / 2,
                       labels=as.character(centromere.hg38$seqnames),
                       expand=c(0,0)) +
    theme(
      panel.background = element_blank(),
      panel.grid.major.y = element_line(colour="gray90", linewidth=0.3),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10)
    )

  return(p)
}



