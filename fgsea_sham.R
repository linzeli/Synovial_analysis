library(fgsea)
library('org.Mm.eg.db')
set.seed(42)
symbols<-row.names(tT_sham)
ent_id<-mapIds(org.Mm.eg.db, symbols, 'ENTREZID', 'SYMBOL')
ent_id.df <- as.data.frame(ent_id)
tT_sham<-cbind(tT_sham, ent_id.df)
tT_sham<-na.omit(tT_sham)
ranks_sham <- tT_sham$logFC
names(ranks_sham)<-tT_sham$ent_id
head(ranks_sham)
fgseaRes_sham <- fgseaMultilevel(pathwaysH, ranks_sham, minSize=15, maxSize = 500)
min(fgseaRes_sham$padj)
valid_sham<-subset(fgseaRes_sham, padj<=0.05)
valid_sham
sham_path<-valid_sham$pathway
sham_edge<-valid_sham$leadingEdge
sham_edge[1]
