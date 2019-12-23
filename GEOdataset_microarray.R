library(edgeR)
library(limma)
library(Biobase)
library(GEOquery)
library(stringr)
library(sva)
library(oligo)

# load series and platform data from GEO
gset <- getGEO("GSE11375", GSEMatrix =TRUE,getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

ex <- exprs(gset)
pheno <- gset@phenoData@data
batch<-gset@phenoData@data[["title"]]

b1_pt <- grepl("group-1-pt", batch)
b2_pt <- grepl("group-2-pt", batch)
b1_ctl <- grepl("group-1-ctl", batch)
b2_ctl <- grepl("group-2-ctl", batch)
batch1<- grepl("group-1-", batch)
batch2<- grepl("group-2-", batch)
group1_pt <- batch[b1_pt]
group2_pt<- batch[b2_pt]
group1_ctl <- batch[b1_ctl]
group2_ctl <- batch[b2_ctl]
group1 <- batch[batch1]
group2 <- batch[batch2]

############################
#####Preprocessing#########
############################
ex = backgroundCorrect.matrix(ex, Eb=NULL, method="auto", offset=0, printer=NULL, normexp.method="saddle", verbose=TRUE)

f1 <- kOverA(5, 200)
ffun <- filterfun(f1)
wh1 <- genefilter(ex, ffun) 
sum(wh1) #this shouldn't change the number of samples, just the one of genes
ex_filt <-exprs(gset)[wh1,]

######################################################################
#differential expression between controls and patients in group 1
######################################################################c=['Array','Batch','MOI','Covariate2','Covariate3','ISS_grp_h','Shock3Cat','ShockYN_h','BD_grp_h']

type <- factor(c(rep("Controls", times = 10),
                 rep("Patients", times = 79)))
ex_group1 <- ex_filt[,1:89]

#Group 1 only
final_ex <- removeBatchEffect(ex_group1)
control_patients <- type
design <- model.matrix(~ -1 + control_patients)
colnames(design) <- levels(control_patients)
fit <- lmFit(final_ex, design)
contrasts = makeContrasts( "Controls-Patients", levels= colnames(coef(fit)))# design )
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2, trend=TRUE)
summary(decideTests(fit2, method="separate"))
boxplot(final_ex[1,]~type,col = "orange", main = "difference in a gene patients vs control")

#####TOP GENES############
top_table = topTable(fit2, coef=1, number=Inf)
df = data.frame("exp" = row.names(top_table[1:10,]), "group" = control_patients)
boxplot(df$exp ~ df$group)

# order samples by group
ex <- ex_group1[ , order(type)]
sml <- type[order(type)]
fl <- as.factor(type)
labels <- c("controls_1","patients_1")

# set parameters and draw the plot
palette(c("#dfeaf4","#f4dfdf", "#AABBCC"))
dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE11375", '/', annotation(gset), " selected samples", sep ='')
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, fill=palette(), bty="n")


#how to check which are the genes that are upregulated in here
library("annotate")
library("hgu133plus2.db")
top_10 <- top_table[1:10,]
#df_top10 <- data.frame("exp" = ex_group1[1,], "group" = control_patients)
rownames(top_10)
top_ex_g1 <- select(hgu133plus2.db, rownames(top_10), c("SYMBOL","ENTREZID", "GENENAME"))

DE_ge_g1 <-top_ex_g1$PROBEID[1:10]
DE_ge_g1 <- c("208858_s_at", "212346_s_at","212346_s_at","37652_at", "227002_at", "217991_x_at", "217729_s_at", "202624_s_at","212056_at", "205718_at")
exp_top_g1 <- final_ex[DE_ge_g1,]
heatmap.2(exp_top_g1,trace = 'none', scale='row', labCol = type, labRow = top_ex_g1$SYMBOL[1:10])

library("EMA")
ID <- featureNames(gset)
symbol <- getSYMBOL(ID,"hgu133plus2.db")
symbol

#################################
#Group2 only btw ctr and patients
#################################
ex_group2 <-ex_filt[,90:184]
type2 <- factor(c(rep("Controls", times = 16),
                 rep("Patients", times = 79)))

final_ex_g2 <- removeBatchEffect(ex_group2)
control_patients_g2 <- type2
design_2 <- model.matrix(~ -1 + control_patients_g2)
colnames(design_2) <- levels(control_patients_g2)
fit_2b <- lmFit(final_ex_g2, design_2)
contrasts_g2 = makeContrasts( "Controls-Patients", levels= colnames(coef(fit_2b)))# design )
fit3 <- contrasts.fit(fit_2b, contrasts_g2)
fit3 <- eBayes(fit3, trend=TRUE)
results <- summary(decideTests(fit3, method="separate"))

top_table_2 <- topTable(fit3, coef=1, number=Inf)
top_10_2 <- top_table_2[1:10,]
rownames(top_10_2)
top_ex_g2 <- select(hgu133plus2.db, rownames(top_10_2), c("SYMBOL","ENTREZID", "GENENAME"))

#heatmap
DE_ge_g2 <- top_ex_g2$PROBEID[1:10]
DE_ge_g2 <-c("202917_s_at", "209604_s_at", "203535_at", "226969_at", "205006_s_at", "219024_at", "228282_at", "244511_at", "209711_at", "227934_at")
top_ge_g2 <- final_ex_g2[DE_ge_g2,]
heatmap.2(top_ge_g2,trace = 'none', scale='row', labCol = type2, labRow = top_ex_g2$SYMBOL[1:10])


######################
##diff btw gr1 and g2
######################


#ex_filt is what we will use

type_total <- factor(c(rep("Controls", times = 10),
                 rep("Patients", times = 79), 
                 rep("Controls", times = 16), 
                 rep("Patients_2", times = 79)))

total <- removeBatchEffect(ex_filt)
control_patients_t <- type_total
design_t <- model.matrix(~control_patients_t)
colnames(design_t) <- levels(control_patients_t)
fit_t <- lmFit(total, design_t)
contrasts_t = makeContrasts(c("Patients-Patients_2"), levels= colnames(coef(fit_t)))# design )
fit_t2 <- contrasts.fit(fit_t, contrasts_t)
fit_t2 <- eBayes(fit_t2, trend=TRUE)
summary(decideTests(fit_t2, method="separate"))
top_table_t <- topTable(fit_t2, coef=1, lfc=2, number=200)
top_table_t2 <-top_table_t[1:10, ]
top_all <- select(hgu133plus2.db, rownames(top_table_t2), c("SYMBOL","ENTREZID", "GENENAME"))

top_genes_v <-top_all$PROBEID[1:30]
###We can filter the ex matrix to have a small matrix that only contains 
###the top 10 differemtially expressed genes
top_DE_genes <- c("235461_at" ,  "243295_at", "211947_s_at", "229010_at", "213998_s_at", "211944_at" ,"242352_at" , "210218_s_at",
                   "219878_s_at", "219679_s_at")
top_expr <-total[top_DE_genes,]
heatmap.2(top_expr, trace = 'none', scale='row', labCol = all_types, labRow = top_all$SYMBOL[1:10])

#Check if my two control groups are the same
all_types <- factor(c(rep("Controls", times = 10),
                       rep("Patients", times = 79), 
                       rep("Controls_2", times = 16), 
                       rep("Patients_2", times = 79)))
plotMDS(total,top = 100, label= all_types)

############################
#####PCA on group1 ctl and patients
library(factoextra)
pca <- prcomp(t(final_ex), scale = TRUE)
fviz_eig(pca)
summary(pca)
plot(pca$x[,1],pca$x[,2], xlab="PC1", ylab = "PC2", main = "PC1 / PC2 - plot")

fviz_pca_ind(pca, geom.ind = "point", pointshape = 21, 
             pointsize = 2, 
             col.ind = "black", 
             palette = "jco", 
             addEllipses = TRUE,
             label = "var",
             col.var = "red",
             repel = TRUE,) +
  ggtitle("2D PCA-plot") +
  theme(plot.title = element_text(hjust = 0.5))

#Just to try 
library(kernlab)
kpc <- kpca(t(final_ex))

plot(rotated(kpc),xlab="1st Principal Component",ylab="2nd Principal Component")




