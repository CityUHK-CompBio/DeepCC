#' Cross Validation of DeepCC Model
#'
#' This function performs cross validation of DeepCC Mode on the training data.
#'
#' @param fs a data.frame containing functional spectra of patients (each row presents one sample)
#' @param labels a character vector containing training lables
#' @param fold a integer indicating the fold number of cross validation (default: 5)
#' @return a numeric indicating error rate in a single run
#' @export
#' @examples
#' cross_validataion(tcga.fs, tcga.labels)
cross_validataion <- function(fs, labels, fold = 5) {
  fs <- fs[!is.na(labels), ]
  labels <- na.omit(labels)

  n <- length(labels)
  testidx <- sample(n, n*(1/fold))

  trainData.fs <- fs[-testidx, ]
  trainData.labels <- labels[-testidx]

  testData.fs <- fs[testidx, ]
  testData.labels <- labels[testidx]

  deepcc.model <- trainDeepCCModel(trainData.fs, trainData.labels)
  pred.lables <- getDeepCCLabels(deepcc.model, testData.fs, 0)

  error_rate <- 1-mean(as.character(testData.labels) == as.character(pred.lables))
  error_rate
}

#' Get Gene Sets
#'
#' This function extract a list of gene sets from gmt file.
#'
#' @param file filename of the gmt file
#' @return a list containing gene sets by EntrezID
#' @export
#' @examples
#' msigdbv51 <- get_gene_sets("msigdb.v5.1.entrez.gmt")

get_gene_sets <- function(file) {
  msig <- GSEABase::getGmt(file, geneIdType=EntrezIdentifier())
  GSEABase::geneIds(msig)
}


#' Visualization of samples
#'
#' This function visualize samples.
#'
#' @param data
#' @param labels
#' @param color
#' @param guide_fill
#' @return a list containing gene sets by EntrezID
#' @import ggplot2 cowplot
#' @export
vis_samples <- function(data, labels, color, guide_fill="legend") {
  normalise <- function(x, na.rm = TRUE) {
    ranx <- range(x, na.rm = na.rm)
    (x - ranx[1]) / diff(ranx)
  }

  calcAvgSilhouetteWidth <- function(data, labels) {
    data <- data[!is.na(labels), ]
    labels <- na.omit(labels)
    labels <- as.numeric(as.factor(labels))
    summary(cluster:::silhouette.default.R(labels, dist(data)))$avg.width
  }

  asw <- calcAvgSilhouetteWidth(data, labels)
  asw_df <- data.frame(x=0.8, y=-0.1, label=paste("ASW =", round(asw, 3)))

  pc <- prcomp(data)
  r2 <- data.frame(PC1=normalise(pc$x[,1]), PC2=normalise(pc$x[,2]), Class=factor(labels))
  r2 <- r2[!is.na(r2$Class), ]

  ggplot(r2, aes(fill=Class, x=PC1, y=PC2)) +
    stat_density2d(geom="tile", aes(x=PC1, y=PC2, alpha=..density..), contour=F, h=c(0.3, 0.3)) +
    scale_alpha_continuous(range = c(0, 1), guide=F) +
    scale_color_manual(values=color, guide=F) +
    scale_fill_manual(values=color, guide=guide_fill) +
    geom_point(aes(colour=Class, x=PC1, y=PC2), size=2, alpha=0.3, stroke = 0) +
    scale_x_continuous(expand=c(0,0), limits = c(-0.2, 1.2)) + scale_y_continuous(expand=c(0,0), limits = c(-0.2, 1.2)) +
    theme(text = element_text(family = "Arial"),
          legend.title = element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          strip.background = element_rect(fill = "white",  linetype = "blank")) +
    geom_text(data=asw_df, aes(x,y,label=label), inherit.aes=FALSE)
}
