library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
EnhancedVolcano(tT_sham,
                lab = rownames(tT_sham),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                FCcutoff = 0.5)
