# Global variables declaration for R CMD check
# This file declares global variables used in NSE (Non-Standard Evaluation)
# contexts such as ggplot2 aes() and dplyr operations

utils::globalVariables(c(
  # ggplot2 aesthetics
  "x", "y", "x1", "y1", "x2", "y2",
  "celltype", "Celltype",
  "src_x", "src_y", "dest_x", "dest_y",
  "tf", "Expressed_genes",
  "cluster", "num", "type", "group",
  "..level..", "..density..",
  
  # foreach loop variables
  "i", "j",
  
  # C++ function references
  "cpp_fast_permutation"
))
