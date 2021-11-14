library(tidyverse)
library(ggvenn)

set.seed(20190708)
genes <- paste("gene",1:1000,sep="")

x<-list(
  Dmm = dmm_path,
  Sham= sham_path)

ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 2
)

selected_path<-dmm_path[!(dmm_path %in% sham_path)]

selected_path
