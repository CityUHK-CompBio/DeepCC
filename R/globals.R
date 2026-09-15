# Column names referenced inside ggplot2 aesthetics and foreach expressions.
# Registering them keeps R CMD check from reporting false-positive notes about
# undefined global variables.
utils::globalVariables(c(
  "Class", "PC1", "PC2", "..density..", "x", "y", "label", "idx"
))
