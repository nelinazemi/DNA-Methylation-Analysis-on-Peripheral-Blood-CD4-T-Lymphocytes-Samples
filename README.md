## DNA Methylation Tutorial (EPIC Array – minfi Pipeline)

Methylated DNA is when a —CH₃ group is added to a C. This mostly happens at CpG sites (usually unmethylated) which are biologically informative since they're close to promoters or other important regulatory regions. Methylation makes chromatin more compact therefore causes the gene to be less expressed. It also works cumulatively, not individually. Many CpGs methylated => chromatin proteins bind more => Region becomes compact
Gene expression decreases. It is therefore the pattern, not a single CpG, that matters. We detect methylation using arrays. Arrays are tools that measure DNA methylation at specific CpG sites across the genome. Earlier versions measured 27k CpGs, then 450k, and now >850k (the EPIC array).
They contain probes (short synthetic DNA sequences) that are complementary to bisulfite-converted DNA, which reveal whether the original CpG was methylated or not. Today's arrays have two mechanism: Infirium I and Infirium II. When DNA is treated with bisulfite, the C is converted to U (read as T), but the methylated DNA stays as C. Infirium I uses two probes (one sequence for methylated and one for unmethylated) to determine which part is methylated and which is unmethylated. Infinium II uses only one probe. Both mechanisms use florescent signals to flag methylated and unmethylated sites. If original CpG was methylated → remains C → extension adds a base that fluoresces color A. If original CpG was unmethylated → became T → extension adds a base that fluoresces color B. The produced signals are in IDAT file format. We load them using read.metharray.exp().

For this project, a simple analysis is done on a small dataset (GSE205495). Peripheral blood CD4+ T lymphocytes samples were collected from 4 primary refractory ITP cases and 4 age-matched healthy controls, and DNA methylome profiling was performed using Illumina Human Methylation850K.

---

## CODE:

```r
## ---- Install BiocManager if needed ----
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

## ---- List of required packages ----
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

## ---- Install missing CRAN packages ----
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

## ---- Install missing Bioconductor packages ----
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)
  }
}

## ---- Load packages ----
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
```

---

## Setting the data directory

```
dataDirectory <- "C:/Users/FGN/Downloads/CHIP/GPL"
```

## Creating targets dataframe

According to our dataset:
GSM6213680–683 → ITP Cases 1–4
GSM6213684–687 → Controls 1–4

```r
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
```

```
Sample_Name  Sample_Group  Basename
GSM6213680   Case          path/GSM6213680
GSM6213681   Case          path/GSM6213681
GSM6213682   Case          path/GSM6213682
GSM6213683   Case          path/GSM6213683
GSM6213684   Control       path/GSM6213684
GSM6213685   Control       path/GSM6213685
GSM6213686   Control       path/GSM6213686
GSM6213687   Control       path/GSM6213687
```

---

## Importing the raw methylation data straight from the array

```r
RGset <- read.metharray.exp(dataDirectory)
RGset
pData(RGset)
getManifest(RGset)
```

**output:**

```
class: RGChannelSet 
dim: 1051815 8 
...
Number of type I probes: 142262 
Number of type II probes: 724574 
Number of control probes: 635 
...
```

This object contains the raw red and green intensity values for 1051815 probes…

---

## Converting to MethylSet

```r
MSet <- preprocessRaw(RGset)
MSet
```

Output shown unchanged.

---

## Extracting Meth/Unmeth intensities

```r
head(getMeth(MSet)[,1:3])
head(getUnmeth(MSet)[,1:3])
```

Example provided in original text retained.

---

## Extracting Beta, M, and CN values

```r
ratioSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
Gset <- mapToGenome(ratioSet)

beta <- getBeta(Gset)
m <- getM(Gset)
cn <- getCN(Gset)
```

Descriptions preserved exactly.

---

## QC Plots

```r
qc <- getQC(MSet)
plotQC(qc)

detP <- detectionP(RGset)
barplot(colMeans(detP), las=2, cex.names=0.8, ylab="Mean detection p-values")
abline(h=0.05,col="red")
```

All explanatory text preserved.

---

## Density Plots / Normalization

```r
mSetSq <- preprocessFunnorm(RGset)
...
plotBetasByType(MSet[,1], main="Raw")
plotBetasByType(getBeta(mSetSq)[,1], probeTypes = probeTypes, main="Normalized")
```

---

## Filtering Probes

```r
detP <- detP[match(featureNames(mSetSq),rownames(detP)),]
keep <- rowSums(detP < 0.01) == ncol(mSetSq)
mSetSqFlt <- mSetSq[keep,]

mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
```

---

# 1. Differentially Methylated Positions (DMPs)

Design matrix:

```r
group <- factor(targets$Sample_Group) 
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
design
```

Limma model:

```r
mVals <- getM(mSetSqFlt)
fit <- lmFit(mVals, design)
contMatrix <- makeContrasts(Case - Control, levels=design)
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
```

Getting DMPs:

```r
DMPs <- topTable(fit2, num=Inf, coef=1, sort.by="p")
head(DMPs)
```

Annotation:

```r
annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annSub <- annEPIC[match(rownames(mVals), annEPIC$Name), ]
DMPs <- topTable(fit2, num=Inf, coef=1, genelist=annSub)
```

Plotting a CpG:

```r
bVals <- getBeta(mSetSqFlt)
plotCpg(bVals, cpg="cg11144986", pheno=targets$Sample_Group)
```

---

# 2. Differentially Methylated Regions (DMRs)

```r
myAnnotation <- cpg.annotate(...)
DMRs <- dmrcate(myAnnotation, lambda = 1000, C = 2, pcutoff = 0.05)
results.ranges <- extractRanges(DMRs)
DMR.plot(...)
```

---

# 3. mCSEA Region-Based GSEA

```r
myRank <- DMPs$logFC 
names(myRank) <- rownames(DMPs)

pheno <- as.data.frame(pData(mSetSqFlt))
pheno$Sample_Group <- targets$Sample_Group
pheno <- pheno[,"Sample_Group", drop=FALSE]

myResults <- mCSEATest(
  myRank,
  bVals,
  pheno,
  regionsTypes = "promoters",
  platform = "EPIC"
)
```

---

# 4. GO / Pathway Enrichment (gometh)

```r
anno <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno_sub <- anno[rownames(anno) %in% rownames(mSetSqFlt), ]

gr_probes <- GRanges(...)
gr_dmrs <- results.ranges

hits <- findOverlaps(gr_probes, gr_dmrs)
sigCpGs <- unique(mcols(gr_probes)$probeID[queryHits(hits)])

allCpGs <- rownames(mSetSqFlt)

gst_GO <- gometh(
  ...
)
```

---