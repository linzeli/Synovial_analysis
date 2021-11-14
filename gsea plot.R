library('org.Mm.eg.db')

heat_immune<-read.csv("/Users/linze/PycharmProjects/untitled/datacleaning/immune_heat.csv",check.names=FALSE,header = TRUE,row.names="pathway")
heat_immune_col_name<-colnames(heat_immune)
heat_immune_col_symbol<-mapIds(org.Mm.eg.db, heat_immune_col_name,'SYMBOL', 'ENTREZID')
colnames(heat_immune)<-heat_immune_col_symbol
gene_list_immune<-read.csv("/Users/linze/PycharmProjects/untitled/datacleaning/immune_occ.csv",check.names=FALSE,header = TRUE)
path_list_immnue<-<-read.csv("/Users/linze/PycharmProjects/untitled/datacleaning/immune_stat.csv",check.names=FALSE,row.names="pathway")

heat_h<-read.csv("/Users/linze/PycharmProjects/untitled/datacleaning/Hallmark_heat.csv",check.names=FALSE,header = TRUE,row.names="pathway")
heat_h_col_name<-colnames(heat_h)
heat_h_col_symbol<-mapIds(org.Mm.eg.db, heat_h_col_name,'SYMBOL', 'ENTREZID')
colnames(heat_h)<-heat_h_col_symbol

gene_list_h<-read.csv("/Users/linze/PycharmProjects/untitled/datacleaning/Hallmark_occ.csv",check.names=FALSE,header = TRUE)
h_gene_col_symbol_list<-mapIds(org.Mm.eg.db, colnames(gene_list_h),'SYMBOL', 'ENTREZID')
colnames(gene_list_h)<-heat_h_col_symbol_list

path_list_h<-<-read.csv("/Users/linze/PycharmProjects/untitled/datacleaning/Hallmark_stat.csv",check.names=FALSE,row.names="pathway")



heat_go<-read.csv("/Users/linze/PycharmProjects/untitled/datacleaning/Hallmark_heat.csv",check.names=FALSE,header = TRUE,row.names="pathway")
heat_go_col_name<-colnames(heat_go)
heat_go_col_symbol<-mapIds(org.Mm.eg.db, heat_go_col_name,'SYMBOL', 'ENTREZID')
colnames(heat_go)<-heat_go_col_symbol

gene_list_go<-read.csv("/Users/linze/PycharmProjects/untitled/datacleaning/go_occ.csv",check.names=FALSE,header = TRUE)
go_gene_col_symbol_list<-mapIds(org.Mm.eg.db, colnames(gene_list_go),'SYMBOL', 'ENTREZID')
colnames(gene_list_go)<-go_gene_col_symbol_list

path_list_go<-read.csv("/Users/linze/PycharmProjects/untitled/datacleaning/go_stat.csv",check.names=FALSE,row.names="pathway")

library(ComplexUpset)

upset(data = heat_h, intersect = heat_h_col_symbol, 
      name="upsetplot for dmm sepcific immune signiture", 
      min_size = 0,
      width_ratio = 0.125) +
  labs(title = "",
       caption = "")
heatmap(as.matrix(heat_h))
heatmap(as.matrix(heat_immune), scale = "none", Rowv = NA, Colv = NA, main = "HeatMap Example") 
