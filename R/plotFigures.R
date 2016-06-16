#' Plot Kaplan-Meier Curve
#'
#' This function plots Kaplan-Meier Curve for survival analysis.
#'
#' @param clinical a survival object created by Surv function in survival package
#' @param labels a vector containing subtyping labels of patients
#' @param annot a character indicating the annotation showed up in the plot
#' @param color a vector containing colors used for different subtypes
#' @param font a numeric indicating the font size used in the plot
#' @param xlab a character indicating the label of x-axis (default: "Follow up (weeks)")
#' @param ylab a character indicating the label of y-axis (default: "DFS (prob.)")
#' @return a ggplot2 object of the plot
#' @import ggplot2
#' @export
#' @examples
#' clinical <- Surv(t.rfs, e.rfs)
#' color.cms <- c("#E69E00","#0070B0","#CA78A6", "#009C73")
#' plotKMCurve(clinical, labels, "GSE39582", color.cms)
plotKMCurve <- function(clinical, labels, annot = NULL, color = NULL, font = 10, xlab = "Follow up (weeks)", ylab = "DFS (prob.)"){
  obj <-  clinical ~ labels
  surv <- survival::survfit(obj)
  survstats <- survival::survdiff(obj)
  survstats$p.value <- 1 - pchisq(survstats$chisq, length(survstats$n) - 1)

  if(is.null(color)) {
    if(length(unique(labels)) < 3) {
      color <- RColorBrewer::brewer.pal(3, "Set1")
    } else {
      color <- RColorBrewer::brewer.pal(length(unique(labels)), "Set1")
    }
    color <- color[1:length(unique(labels))]
  }

  p <- GGally::ggsurv(surv, surv.col = color, cens.col = color, xlab = xlab, ylab = ylab) +
    annotate("text", x = Inf, y = Inf, label = paste0("P==", fancy_scientific(survstats$p.value, 3)), hjust = 1.2, vjust = 2, size=font/3, parse = TRUE) +
    theme_bw() +
    theme(text = element_text(size=font),
          legend.title=element_blank(),
          axis.text=element_text(size=font),
          axis.title=element_text(size=font*1.2))

  # annotate data set name
  if(!is.null(annot)) p <- p + annotate("text", x = 0, y = 0, label = annot, hjust = 0, vjust = 0, size=font/3)

  p
}
