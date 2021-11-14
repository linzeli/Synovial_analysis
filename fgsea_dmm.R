library(fgsea)
library('org.Mm.eg.db')
set.seed(42)
symbols<-row.names(tT_dmm)
ent_id<-mapIds(org.Mm.eg.db, symbols, 'ENTREZID', 'SYMBOL')
ent_id.df <- as.data.frame(ent_id)
tT_dmm<-cbind(tT_dmm, ent_id.df)
tT_dmm<-na.omit(tT_dmm)
ranks_dmm <- tT_dmm$logFC
names(ranks_dmm)<-tT_dmm$ent_id

head(ranks_dmm)
fgseaRes_dmm <- fgseaMultilevel(pathwaysH, ranks_dmm, minSize=15, maxSize = 500)
valid_dmm<-subset(fgseaRes_dmm, padj<=0.05)
valid_dmm
dmm_path<-valid_dmm$pathway
dmm_edge<-valid_dmm$leadingEdge

