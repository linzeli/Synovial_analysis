#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO

gset <- getGEO("GSE99738", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL21163", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "0123321023011032"
sml <- strsplit(gsms, split="")[[1]]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

df<-read.csv("/Users/linze/PycharmProjects/untitled/datacleaning/GSE99738_mRNA.csv",header = TRUE,row.names="GENE_SYMBOL")
df = subset(df, select = -c(MAD,variance,Entropy))

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("dmm1","sham1","dmm6","sham6"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(df, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts_dmm1v6 <- c(paste(groups[3],"-",groups[1],sep=""))
cts_sham1v6 <- c(paste(groups[4],"-",groups[2],sep=""))

cont.matrix_1 <- makeContrasts(contrasts=cts_dmm1v6, levels=design)
fit_dmm <- contrasts.fit(fit, cont.matrix_1)

cont.matrix_2 <- makeContrasts(contrasts=cts_sham1v6, levels=design)
fit_sham <- contrasts.fit(fit, cont.matrix_2)

# compute statistics and table of top significant genes
fit_dmm <- eBayes(fit_dmm, 0.01)
tT_dmm <- topTable(fit_dmm, adjust="fdr", sort.by="B",number=5000)
tT_dmm<-subset(tT_dmm, adj.P.Val<=0.05)
tT_dmm<-subset(tT_dmm,abs(logFC)>=0.5)
tT_dmm_up<-subset(tT_dmm,logFC>=0.5)
tT_dmm_down<-subset(tT_dmm,logFC<=-0.5)
fit_sham <- eBayes(fit_sham, 0.01)
tT_sham <- topTable(fit_sham, adjust="fdr", sort.by="B", number=5000)
tT_sham<-subset(tT_sham, adj.P.Val<=0.05)
tT_sham<-subset(tT_sham,abs(logFC)>=0.5)
tT_sham_up<-subset(tT_sham,logFC>=0.5)
tT_sham_down<-subset(tT_sham,logFC<=-0.5)
dmm_name<-row.names(tT_dmm)
sham_name<-row.names(tT_sham)
dif_dmm<-setdiff(as.vector(dmm_name),as.vector(sham_name))
dif_sham<-setdiff(as.vector(sham_name),as.vector(dmm_name))
dmm_spec<-tT_dmm[dif_dmm,]
sham_spec<-tT_sham[dif_sham,]
selected_name<-c(dmm_name,dif_sham)
select_gene<-df[selected_name,]
write.csv(select_gene,'/Users/linze/PycharmProjects/untitled/datacleaning/selected.csv')

