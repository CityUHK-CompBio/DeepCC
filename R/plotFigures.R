#' Plot Kaplan-Meier Curve
#'
#' This function plots Kaplan-Meier Curve for survival analysis.
#'
#' @param clinical a survival object created by Surv function in survival package
#' @param labels a vector containing subtyping labels of patients
#' @param data_set a character indicating the name of the data set
#' @param color a vector containing colors used for different subtypes
#' @param font a numeric indicating the font size used in the plot
#' @return a ggplot2 object of the plot
#' @import ggplot2
#' @export
#' @examples
#' clinical <- Surv(t.rfs, e.rfs)
#' color.cms <- c("#E69E00","#0070B0","#CA78A6", "#009C73")
#' plotKMCurve(clinical, labels, "GSE39582", color.cms)
plotKMCurve <- function(clinical, labels, data_set = NULL, color = NULL, font = 10){
  obj <-  clinical ~ labels
  surv <- survival::survfit(obj)
  survstats <- survival::survdiff(obj)
  survstats$p.value <- 1 - pchisq(survstats$chisq, length(survstats$n) - 1)

  if(is.null(color)) color <- RColorBrewer::brewer.pal(length(unique(labels)), "Set1")[1:length(unique(labels))]

  p <- GGally::ggsurv(surv, surv.col = color, cens.col = color) +
    annotate("text", x = Inf, y = Inf, label = paste0("P==", fancy_scientific(survstats$p.value, 3)), hjust = 1.2, vjust = 2, size=font/2, parse = TRUE) +
    xlab("Follow up (weeks)") + ylab("DFS (prob.)") +
    theme_bw() + theme(text = element_text(size=font), legend.title=element_blank())

  # annotate data set name
  if(!is.null(data_set)) p <- p + annotate("text", x = 0, y = 0, label = data_set, hjust = 0, vjust = 0, size=font/2)

  p
}
