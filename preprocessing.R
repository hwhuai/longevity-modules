library(org.Ce.eg.db)
library(DESeq2) # version 1.24.0
library(mixOmics)
library(RUVSeq)
library(RColorBrewer)
library(pheatmap)
library(PoiClaClu)
library(EDASeq)
library(mixOmics)
library(dplyr)

#This is a test preprocessing samples
setwd('/Users/weihanhuai1/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/clock analysis/preprocessing')
source("../Functions/networkFunctions-extras-18.R")
source("../Functions/labelPoints2-01.R");
source("../Functions/outlierRemovalFunctions.R")
source("../Functions/preprocessing-General-013.R")
source("../Functions/GNVFunctions-016.R")
source("../Functions/RUVseq covariance removal.R")

dir.create("RData", recursive = TRUE);
dir.create("Results", recursive = TRUE);
dir.create("Plots", recursive = TRUE);

#read data
data = read.csv('celegans_rawdata.csv', row.names = 1)
GeneAnnot = AnnotationDbi::select(org.Ce.eg.db, 
                                         keys = keys(org.Ce.eg.db,keytype="ENTREZID"),
                                         columns = c('ENTREZID','SYMBOL','GENENAME',"ENSEMBL"),
                                         keytypes = "ENSEMBL"
)


GeneAnnot<-GeneAnnot %>% filter(!is.na(GeneAnnot$ENSEMBL))
GeneAnnot<-GeneAnnot %>% group_by(ENSEMBL) %>% mutate(x=n_distinct(ENSEMBL))%>% filter(x==1)

#Dataset 1: N2 data with different ages
N2 = data[,c(1:80)]
N2 = sapply(N2, as.integer)
rownames(N2) <- rownames(data)

index<-match(GeneAnnot$ENSEMBL,rownames(N2))
any(is.na(index))

N2 = N2[index,]
rownames(N2) = GeneAnnot$SYMBOL
N2 = N2[which(rowSums(N2) != 0),]

#Dataset 2: drug data
drug = data[,c(145:180,220:236)]
drug = sapply(drug, as.integer)
rownames(drug) <- rownames(data)

index2 = match(GeneAnnot$ENSEMBL,rownames(drug))
any(is.na(index2))

drug = drug[index2,]
rownames(drug) <- GeneAnnot$SYMBOL
drug = drug[which(rowSums(drug) != 0),]
colnames(drug) <- colnames(drug) %>% gsub('\\.','-',.)



#Load metadata
metadata = read.csv('metadata.csv', row.names = 1)
metadata_N2 = metadata[which(metadata$Sample_id %in% colnames(N2)),]
metadata_drug = metadata[c(145:180,220:236),]
metadata_drug[c(19:20,35:36),2] = 'N2_D6'

#PCA
pca = prcomp(t(N2), center = T)

clock_pca = data.frame(metadata_N2, pca$x[,1:9])

#ggplot visualization
ggplot(clock_pca, aes(x = PC1, y = PC2, color = Batch, shape = condition))+
  geom_point()+
  scale_shape_manual(values=0:50) +
  scale_colour_continuous()+
  scale_fill_discrete() +
  theme_bw()     

#Put into DESeq2
ddsTxi_N2 = DESeqDataSetFromMatrix(N2,
                                   colData = metadata_N2,
                                   design = ~1)
keep_N2 <- rowSums(cpm(ddsTxi_N2) >= 0.5) >= 2
ddsTxi_N2 = ddsTxi_N2[keep_N2,]
vsd = vst(ddsTxi_N2)
N2_vst = assay(vsd)

conditions = metadata_N2$condition
diff = makeGroups(conditions)



plotHeatmapRUV_k4DEG_optimal <- function(countMatrix, plot_PCA = F, names = conditions, differences = diff,
                                         controlGenes = housekeeping, optimal_performance = 1){
  set.seed(123)
  optimal_k = 1
  x = countMatrix
  x = x[which(rowVars(x) != 0), ]
  x = x[which(rowSums(x) != 0), ]
  # cv.out = Classify.cv(x, as.numeric(as.factor(names)), rhos = 0, nfolds = 2, type="deseq", transform=TRUE)
  out = Classify(t(x),as.numeric(as.factor(names)),rho=0, xte = NULL) # cv.out$bestrh
  tabela = table(out$ytehat,as.numeric(as.factor(names)))
  print(tabela)
  performance = sum(diag(tabela))/length(names)
  if (performance == 1){
    print("It does not have any unwanted source of variation.")
    return(x)
  } else {
    print("It does have unwanted sources of variation.")
    while (performance < optimal_performance) {
      setRUVs <- RUVs(round(countMatrix), k = optimal_k, scIdx = differences, cIdx = controlGenes)
      x = setRUVs$normalizedCounts
      x = x[which(rowVars(x) != 0), ]
      x = x[which(rowSums(x) != 0), ]
      # cv.out = Classify.cv(x, as.numeric(as.factor(names)), rhos = NULL, nfolds = 3, type="deseq", transform=TRUE)
      out = Classify(t(x),as.numeric(as.factor(names)),rho=0, xte = NULL) # cv.out$bestrho
      tabela = table(out$ytehat,as.numeric(as.factor(names)))
      print(cbind(out$ytehat,as.numeric(as.factor(names))))
      print(tabela)
      performance = sum(diag(tabela))/length(names)
      print(optimal_k)
      print(performance)
      optimal_k = optimal_k + 1
    }
  }
  if (plot_PCA == T){
    plotIndiv(pca(t(x), ncomp = 2, center = T, scale = F), group = names, ind.names = T, legend = TRUE, title = 'PCA',
              ellipse = T, centroid = F)
  }
  poisd <- PoissonDistance(t(x), type = "deseq")
  samplePoisDistMatrix <- as.matrix( poisd$dd )
  rownames(samplePoisDistMatrix) <- names # paste( names, batch, sep=" - " )
  colnames(samplePoisDistMatrix) <- NULL
  pheatmap(samplePoisDistMatrix,
           clustering_distance_rows = poisd$dd,
           clustering_distance_cols = poisd$dd) #,
  # col = colors)
  return(setRUVs)
}
#Removing sources of variations
countsRUVsDownstream_N2 = plotHeatmapRUV_k4DEG_optimal(as.matrix(N2_vst), T, controlGenes = rownames(N2_vst))

library(sva)
N2_EBLM = empiricalBayesLM(t(N2_vst), removedCovariates = metadata_N2$Batch, retainedCovariates = metadata_N2$condition)$adjustedData

mod_N2 = model.matrix(~as.factor(condition), data=metadata_N2)
mod0_N2 = model.matrix(~1, data = metadata_N2)
N2.sva.n <- num.sv(dat = N2_vst, mod = mod_N2)

sva_N2 = sva(N2_vst, mod = mod_N2, mod0 = mod0_N2, n.sv = N2.sva.n)

N2_EBLM_sva = empiricalBayesLM(t(N2_vst), removedCovariates = sva_N2$sv[,1:2], retainedCovariates = metadata_N2$condition)$adjustedData


N2_filter = iterativeOutlierRemoval.withDetails(t(N2_vst), Zcut = 6, nTopProbes = 9000, robustScale = T, scaleData = T)
N2_EBLM_sva_filtered = N2_EBLM_sva[which(rownames(N2_EBLM_sva) %in% N2_filter$details$Step.1$sample),]
N2_EBLM_sva = N2_EBLM_sva_filtered
metadata_N2 =  metadata_N2[which(metadata_N2$Sample_id %in% rownames(N2_EBLM_sva_filtered)),]


plotIndiv(pca(t(N2_vst), ncomp = 2, center = T, scale = F), group = metadata_N2$condition, ind.names = F, legend = TRUE, title = 'PCA',
          ellipse = T, centroid = F)

#N2_EBLM is the final form of hte first dataset.
pdf('Plots/PCA_N2.pdf')
par(mfrow = c(1,3))
par(mar = c(3.5, 3.5, 1.5, 1))
par(mgp = c(2.1, 0.8, 0))
par(cex = 1);
plotIndiv(pca(N2_EBLM_sva, ncomp = 3, center = T, scale = F), group = metadata_N2$condition, ind.names = F, legend = TRUE, title = 'PCA',
          ellipse = T, centroid = F)
dev.off();



#Drug data
ddsTxi_drug = DESeqDataSetFromMatrix(drug,
                                     colData = metadata_drug,
                                     design = ~1)
keep_drug <- rowSums(cpm(ddsTxi_drug) >= 0.5) >= 2
ddsTxi_drug = ddsTxi_drug[keep_drug,]
drug_count_W = as.matrix(as.data.frame(counts(ddsTxi_drug)))
vsd_drug = vst(ddsTxi_drug)
drug_vst = assay(vsd_drug)

drug_filter = iterativeOutlierRemoval.withDetails(t(drug_vst), Zcut = 6, nTopProbes = 9000, robustScale = T, scaleData = T)
drug_filtered = drug_vst[,which(colnames(drug_vst) %in% drug_filter$details$Step.2$sample)]
drug_count_W = drug_count_W[,which(colnames(drug_count_W) %in% colnames(drug_filtered))]
drug_vst = drug_vst[,which(colnames(drug_vst) %in% colnames(drug_filtered))]

metadata_drug_filtered = metadata_drug[which(metadata_drug$Sample_id %in% colnames(drug_count_W)),]

conditions = metadata_drug_filtered$condition
diff = makeGroups(conditions)
countsRUVsDownstream_drug = plotHeatmapRUV_k4DEG_optimal(drug_count_W, T, controlGenes = rownames(drug_count_W), optimal_performance = 0.96)

#drug_EBLM = empiricalBayesLM(t(drug_filtered), removedCovariates = metadata_drug_filtered$Batch, fitToSamples = metadata_drug_filtered$condition=='N2_D6')$adjustedData
drug_EBLM = t(drug_vst)
drug_EBLM_filtered = t(na.omit(t(drug_EBLM)))

mod_drug = model.matrix(~as.factor(condition), data=metadata_drug_filtered)
mod0_drug = model.matrix(~1, data = metadata_drug_filtered)
drug.sva.n <- num.sv(dat = drug_vst, mod = mod_drug)

sva_drug = sva(drug_vst, mod = mod_drug, mod0 = mod0_drug, n.sv = drug.sva.n)

drug_EBLM_sva = empiricalBayesLM(drug_EBLM_filtered, removedCovariates = sva_drug$sv[,1], retainedCovariates = metadata_drug_filtered$condition)$adjustedData
drug_EBLM_sva = t(na.omit(t(drug_EBLM_sva)))

plotIndiv(pca(drug_EBLM_sva, ncomp = 2, center = T, scale = F), group = metadata_drug_filtered$condition, ind.names = F, legend = TRUE, title = 'PCA',
          ellipse = T, centroid = F)

# Check the correspondence between original outlier Z and total sample expression
pdf('Plots/Z1_N2.pdf')
Z1 = N2_filter$details$Step.1$Z
te = rowSums(t(N2_vst))
plot(te, Z1, log = 'x', ylim = c(-8,8), main = 'Z outlier removal for N2')
abline(h=c(-6,6), col = "red")
dev.off()

pdf('Plots/Z1_drug.pdf')
Z1 = drug_filter$details$Step.1$Z
te = rowSums(t(assay(vsd_drug)));
plot(te, Z1, log = 'x', ylim = c(-8,18), main = 'Z outlier removal for Drug iteration 1')
abline(h=c(-6,6), col = "red")
dev.off()

pdf('Plots/Z2_drug.pdf')
Z2 = drug_filter$details$Step.2$Z
te2 = rowSums(drug_EBLM_filtered);
plot(te2, Z2, log = 'x', ylim = c(-8,8), main = 'Z outlier removal for Drug iteration 2')
abline(h=c(-6,6), col = "red")
dev.off()

pdf('Plots/PCA_drugs.pdf')
par(mfrow = c(1,3))
par(mar = c(3.5, 3.5, 1.5, 1))
par(mgp = c(2.1, 0.8, 0))
par(cex = 1);
plotIndiv(pca(drug_EBLM_sva, ncomp = 3, center = T, scale = F), group = c(metadata_drug_filtered$condition), ind.names = F, legend = TRUE, title = 'PCA',
          ellipse = T, centroid = F)
dev.off();


#Get gene annotations for each dataset
geneAnnot_N2 = GeneAnnot[GeneAnnot$SYMBOL %in% colnames(N2_EBLM_sva),]
geneAnnot_drug = GeneAnnot[GeneAnnot$SYMBOL %in% colnames(drug_EBLM_sva),]


ageLevelTranslation.paper = cbind(c("N2_D6", "CF_D6", "N2_D8", "N2_D11", "N2_D14", "N2_D18"),
                                   c("N2 Day6", "Daf-16 mutant Day6","N2 Day8", "N2 Day11", "N2 Day14", "N2 Day18"));
ageOrder.paper = c("N2 Day6", "Daf-16 mutant Day6", "N2 Day8", "N2 Day11", "N2 Day14", "N2 Day18");
ageOrder = translateUsingTable(ageOrder.paper, ageLevelTranslation.paper[, c(2,1)]);
ageLevels = ageOrder;


bin_N2 = binarizeCategoricalColumns.pairwise(metadata_N2['condition'], convertColumns = 'condition', levelOrder = list(ageLevels),
                                          includePrefix = FALSE, prefixSep = "");
bin_N2 = bin_N2[,c(1:5,10:15)]
metadata_N2 = cbind(metadata_N2, bin_N2, sva_N2$sv[,1:3], countsRUVsDownstream_N2$W)
names(metadata_N2)[22:24] = c('SVA.1','SVA.2','SVA.3')

drugLevelTranslation.paper = cbind(c('N2_D6','N2_Allan_D6','N2_Met_D6','N2_Psora_D6','N2_Rap_D6','N2_Rif_D6','eat2_D6','N2_RapAllan_D6',
                                     'N2_RapPsora_D6','N2_PsoraAllan_D6','N2_RifAllan_D6','N2_RifPsora_D6','N2_RapMet_D6','N2_RapRif_D6','N2_RapRifPsora_D6','N2_RifPsoraAllan_D6','N2_RifRapAllan_D6'),
                                   c('N2 (20.5)', 'Allantonin (24.12)','Metformin (24)','Psora-4 (23.5)','Rapamycin (24.4)','Rifampicin (27)','Eat-2 mutant (25)', 'RapAllan (19.9)', 'RapPsora (23.5)',
                                     'PsoraAllan (21.2)','RifAllan (26.5)','RifPsora (30)','RapMet (27.25)','RapRif (29.2)','RapRifPsora (25.8)','RifPsoraAllan (34.85)','RifRapAllan (34.95)'))
drugOrder.paper = c('N2 (20.5)', 'Allantonin (24.12)','Metformin (24)','Psora-4 (23.5)','Rapamycin (24.4)','Rifampicin (27)','Eat-2 mutant (25)', 'RapAllan (19.9)', 'RapPsora (23.5)',
                    'PsoraAllan (21.2)','RifAllan (26.5)','RifPsora (30)','RapMet (27.25)','RapRif (29.2)','RapRifPsora (25.8)','RifPsoraAllan (34.85)','RifRapAllan (34.95)')
drugOrder = translateUsingTable(drugOrder.paper, drugLevelTranslation.paper[, c(2,1)]);
drugLevels = drugOrder;

bin_drug = binarizeCategoricalColumns.pairwise(metadata_drug_filtered['condition'], convertColumns = 'condition', levelOrder = list(drugLevels),
                                               includePrefix = FALSE, prefixSep = "");
metadata_drug_filtered = cbind(metadata_drug_filtered, bin_drug, sva_drug$sv[,1:4], countsRUVsDownstream_drug$W)
names(metadata_drug_filtered)[147:150] = c('SVA.1','SVA.2','SVA.3','SVA.4')
metadata_drug_filtered = metadata_drug_filtered[,c(1:26,32, 34:35,40:41,51,60:61,63,66:67,71:72,76:78,80,85:86,88:90,91,110,116,124,130:131,134:135,141,143,147:153)]

metadata_N2$condition_paper = translateUsingTable(metadata_N2$condition, ageLevelTranslation.paper)
metadata_drug_filtered$condition_paper = translateUsingTable(metadata_drug_filtered$condition, drugLevelTranslation.paper)

N2_count = N2[which(rownames(N2) %in% colnames(N2_EBLM_sva)),]
drug_count = drug_count_W[which(rownames(drug_count_W) %in% colnames(drug_EBLM_sva)),]

save(N2_EBLM_sva, metadata_N2, drug_EBLM_sva, metadata_drug_filtered, geneAnnot_N2, geneAnnot_drug, N2_count, drug_count, file = 'RData/preprocessing.RData')
