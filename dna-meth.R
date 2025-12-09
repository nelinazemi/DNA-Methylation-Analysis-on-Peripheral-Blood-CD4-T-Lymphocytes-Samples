## ===========================
## ---- Install BiocManager if needed ----
## ===========================
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

## ===========================
## ---- List of required packages ----
## ===========================
bioc_pkgs <- c(
  "limma",
  "minfi",
  "missMethyl",
  "DMRcate",
  "mCSEA",
  "mCSEAdata",
  "Gviz",
  "GenomicRanges",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "TxDb.Hsapiens.UCSC.hg19.knownGene",
  "org.Hs.eg.db"
)

cran_pkgs <- c(
  "RColorBrewer",
  "stringr"
)

## ===========================
## ---- Install missing CRAN packages ----
## ===========================
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

## ===========================
## ---- Install missing Bioconductor packages ----
## ===========================
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)
  }
}

## ===========================
## ---- Load packages ----
## ===========================
library(limma)
library(minfi)
library(missMethyl)
library(DMRcate)
library(mCSEA)
library(mCSEAdata)
library(Gviz)
library(GenomicRanges)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(RColorBrewer)
library(stringr)

## ===========================
## ---- Load IDAT files ----
## ===========================
dataDirectory <- "path/GSE205495"
list.files(dataDirectory, recursive = TRUE)

files <- list.files(dataDirectory, pattern="idat$", full.names=TRUE)
head(files)
gsm_ids <- unique(gsub("_(.*)", "", basename(files)))
gsm_ids

sample_group <- ifelse(grepl("368[0-3]$", gsm_ids), "Case", "Control")
targets <- data.frame(
  Sample_Name = gsm_ids,
  Sample_Group = sample_group,
  stringsAsFactors = FALSE
)
targets$Basename <- file.path(dataDirectory, targets$Sample_Name)
targets

RGset <- read.metharray.exp(dataDirectory)
RGset
pData(RGset)
getManifest(RGset)

## ===========================
## ---- Raw preprocessing ----
## ===========================
MSet <- preprocessRaw(RGset)
MSet

head(getMeth(MSet)[,1:3])
head(getUnmeth(MSet)[,1:3])

ratioSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
ratioSet

Gset <- mapToGenome(ratioSet)
Gset

beta <- getBeta(Gset)
head(beta)
m <- getM(Gset)
head(m)
cn <- getCN(Gset)
head(cn)

## ===========================
## ---- Quality control ----
## ===========================
qc <- getQC(MSet)
plotQC(qc)

detP <- detectionP(RGset)
head(detP)
barplot(colMeans(detP), las=2, cex.names=0.8, ylab="Mean detection p-values")
abline(h=0.05,col="red")

phenoData <- pData(MSet)
densityPlot(MSet, sampGroups = phenoData$Sample_Group)

## ===========================
## ---- Functional normalization ----
## ===========================
mSetSq <- preprocessFunnorm(RGset)

par(mfrow=c(1,2))
plotBetasByType(MSet[,1],main="Raw")
typeI <- getProbeInfo(MSet, type = "I")[, c("Name","nCpG")]
typeII <- getProbeInfo(MSet, type = "II")[, c("Name","nCpG")]
probeTypes <- rbind(typeI, typeII)
probeTypes$Type <- rep(x = c("I", "II"), times = c(nrow(typeI), nrow(typeII)))
plotBetasByType(getBeta(mSetSq)[,1], probeTypes = probeTypes, main="Normalized",)

## ===========================
## ---- Filtering failed probes ----
## ===========================
detP <- detP[match(featureNames(mSetSq),rownames(detP)),]
keep <- rowSums(detP < 0.01) == ncol(mSetSq)
table(keep)
mSetSqFlt <- mSetSq[keep,]
mSetSqFlt

mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
mSetSqFlt

## ===========================
## ---- Differential methylation (DMPs) ----
## ===========================
group <- factor(targets$Sample_Group) 
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
design

mVals <- getM(mSetSqFlt)
fit <- lmFit(mVals, design)
contMatrix <- makeContrasts(Case - Control, levels=design)
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

DMPs <- topTable(fit2, num=Inf, coef=1, sort.by="p")
head(DMPs)

annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annSub <- annEPIC[match(rownames(mVals), annEPIC$Name), ]
annSub
DMPs <- topTable(fit2, num=Inf, coef=1, genelist=annSub)
head(DMPs)

## ===========================
## ---- Plot example CpG ----
## ===========================
bVals <- getBeta(mSetSqFlt)
plotCpg(bVals, cpg = "cg11144986", pheno = targets$Sample_Group)

## ===========================
## ---- Differentially methylated regions (DMRs) ----
## ===========================
myAnnotation <- cpg.annotate(
  object = mVals,
  datatype = "array",
  what = "M",
  analysis.type = "differential",
  design = design,
  contrasts = TRUE,
  cont.matrix = contMatrix,
  coef = "Case - Control",
  arraytype = "EPICv1"
)
myAnnotation

DMRs <- dmrcate(myAnnotation, lambda = 1000, C = 2, pcutoff = 0.05)
DMRs

results.ranges <- extractRanges(DMRs)
results.ranges

cols <- ifelse(targets$Sample_Group == "Case", "red", "blue")
DMR.plot(
  ranges = results.ranges,
  dmr = 1,
  CpGs = mSetSqFlt,
  phen.col = cols,
  genome = "hg19"
)

## ===========================
## ---- mCSEA analysis ----
## ===========================
myRank <- DMPs$logFC
names(myRank) <- rownames(DMPs)

pheno <- as.data.frame(pData(mSetSqFlt))
pheno$Sample_Group <- targets$Sample_Group
pheno <- pheno[,"Sample_Group", drop=FALSE]
pheno

myResults <- mCSEATest(
  myRank,
  bVals,
  pheno,
  regionsTypes = "promoters",
  platform = "EPIC"
)

head(myResults$promoters)

## ===========================
## ---- Extract CpGs inside DMRs ----
## ===========================
anno <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

anno_sub <- anno[rownames(anno) %in% rownames(mSetSqFlt), ]

gr_probes <- GRanges(
  seqnames = anno_sub$chr,
  ranges   = IRanges(start = anno_sub$pos, end = anno_sub$pos),
  strand   = "*",
  probeID  = rownames(anno_sub)
)

gr_dmrs <- results.ranges 

hits <- findOverlaps(gr_probes, gr_dmrs)
sigCpGs <- unique(mcols(gr_probes)$probeID[queryHits(hits)])
length(sigCpGs)

allCpGs <- rownames(mSetSqFlt)

## ===========================
## ---- GO enrichment (missMethyl) ----
## ===========================
gst_GO <- gometh(
  sig.cpg = sigCpGs,
  all.cpg = allCpGs,
  collection = "GO"
)

topGSA(gst_GO, number = 10)

## ===========================
## ---- Save workspace ----
## ===========================
save.image()
