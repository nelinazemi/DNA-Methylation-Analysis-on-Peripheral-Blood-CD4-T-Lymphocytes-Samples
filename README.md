# Methylated DNA Analysis Tutorial

## Overview

Methylated DNA is when a —CH₃ group is added to a C. This mostly happens at CpG sites (usually unmethylated), which are biologically informative since they're close to promoters or other important regulatory regions. Methylation makes chromatin more compact, therefore causing the gene to be less expressed. It also works cumulatively, not individually. Many CpGs methylated => chromatin proteins bind more => Region becomes compact => Gene expression decreases. It is therefore the pattern, not a single CpG, that matters.

We detect methylation using arrays. Arrays are tools that measure DNA methylation at specific CpG sites across the genome. Earlier versions measured 27k CpGs, then 450k, and now >850k (the EPIC array). They contain probes (short synthetic DNA sequences) that are complementary to bisulfite-converted DNA, which reveal whether the original CpG was methylated or not. Today's arrays have two mechanisms: Infinium I and Infinium II.

When DNA is treated with bisulfite, the C is converted to U (read as T), but the methylated DNA stays as C. Infinium I uses two probes (one sequence for methylated and one for unmethylated) to determine which part is methylated and which is unmethylated. Infinium II uses only one probe. Both mechanisms use fluorescent signals to flag methylated and unmethylated sites. If the original CpG was methylated → remains C → extension adds a base that fluoresces color A. If the original CpG was unmethylated → became T → extension adds a base that fluoresces color B. The produced signals are in IDAT file format. We load them using `read.metharray.exp()`.

For this project, a simple analysis is done on a small dataset (GSE205495). Peripheral blood CD4+ T lymphocytes samples were collected from 4 primary refractory ITP cases and 4 age-matched healthy controls, and DNA methylome profiling was performed using Illumina Human Methylation850K.

## Pipeline

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

```r
# Setting the data directory
dataDirectory <- "C:/Users/FGN/Downloads/CHIP/GPL"
```
```r
# Creating targets dataframe
According to our dataset:
GSM6213680–683 → ITP Cases 1–4
GSM6213684–687 → Controls 1–4
```
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
```
```r
targets
Sample_Name   Sample_Group   Basename
GSM6213680    Case          path/GSM6213680
GSM6213681    Case          path/GSM6213681
GSM6213682    Case          path/GSM6213682
GSM6213683    Case          path/GSM6213683
GSM6213684    Control       path/GSM6213684
GSM6213685    Control       path/GSM6213685
GSM6213686    Control       path/GSM6213686
GSM6213687    Control       path/GSM6213687
```
```r
# Importing the raw methylation data straight from the array
RGset <- read.metharray.exp(dataDirectory)
RGset
pData(RGset)
getManifest(RGset)
output:
class: RGChannelSet 
dim: 1051815 8 
metadata(0):
assays(2): Green Red
rownames(1051815): 1600101 1600111 ... 99810990 99810992
rowData names(0):
colnames(8): GSM6213680_204027220018_R05C01 GSM6213681_204027220018_R06C01 ...
  GSM6213686_204027220018_R03C01 GSM6213687_204027220018_R04C01
colData names(0):
Annotation
  array: IlluminaHumanMethylationEPIC
  annotation: ilm10b4.hg19
```
```r
Number of type I probes: 142262 
Number of type II probes: 724574 
Number of control probes: 635 
Number of SNP type I probes: 21 
Number of SNP type II probes: 38 
```
This object contains the raw red and green intensity values for 1051815 probes, directly from the Illumina EPIC methylation array before any processing.
Rows = CpG probes on the array which corresponds to: a specific CpG site, or a control probe or a SNP probe
Columns = samples

```r
# The object is next turned into MethylSet and converts the red/green signals to methylated and unmethylated labels and intensities (but still raw and unnormalized).
MSet <- preprocessRaw(RGset)
MSet
output:
class: MethylSet 
dim: 866091 8 
metadata(0):
assays(2): Meth Unmeth
rownames(866091): cg18478105 cg09835024 ... cg10633746 cg12623625
rowData names(0):
colnames(8): GSM6213680_204027220018_R05C01 GSM6213681_204027220018_R06C01 ...
  GSM6213686_204027220018_R03C01 GSM6213687_204027220018_R04C01
colData names(0):
Annotation
  array: IlluminaHumanMethylationEPIC
  annotation: ilm10b4.hg19
Preprocessing
  Method: Raw (no normalization or bg correction)
  minfi version: 1.50.0
  Manifest version: 0.3.0

head(getMeth(MSet)[,1:3])
head(getUnmeth(MSet)[,1:3])
```
```r
output:
 head(getMeth(MSet)[,1:3])
           GSM6213680_204027220018_R05C01 GSM6213681_204027220018_R06C01 GSM6213682_204027220018_R07C01
cg18478105                            541                            549                            721

head(getUnmeth(MSet)[,1:3])
           GSM6213680_204027220018_R05C01 GSM6213681_204027220018_R06C01 GSM6213682_204027220018_R07C01
cg18478105                          14753                          14201                          17334
```
for example in one sample:
Meth = 721 
Unmeth = 17334
these numbers mean:
How brightly the “methylated probe” glowed
vs
How brightly the “unmethylated probe” glowed
It's obvious that the unmethylated probe glowed more. So this specific CpG in this sample is definitely less methylated.

```r
# Extracting the Beta and M-values
ratioSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
ratioSet
output:
> ratioSet
class: RatioSet 
dim: 485512 95 
metadata(0):
assays(3): Beta M CN

# Getting genomic coordinates
Gset <- mapToGenome(ratioSet)
Gset

beta <- getBeta(Gset)
head(beta)
m <- getM(Gset)
head(m)
cn <- getCN(Gset)
head(cn)
```
Beta:
This extracts the beta values = methylation level for each CpG. Beta is calculated as:

\[
\text{Beta} = \frac{\text{Meth}}{\text{Meth} + \text{Unmeth}} \times 100
\]

Ranges from 0 (unmethylated) to 1 (fully methylated). Example:

cg13869341   0.877   0.774   0.798
Meaning: this CpG is very highly methylated in those samples.

M:
This extracts M-values, which represent:

\[
M = \log_2\left(\frac{\text{Meth}}{\text{Unmeth}}\right)
\]

Used for statistical tests (DESeq-like). Positive M → more methylated, Negative M → less methylated.

C:
This extracts copy number estimates, basically:
\[
CN = \text{Meth} + \text{Unmeth}
\]

it represents how much signal the probe produced. It’s used to detect low-quality probes, hybridization issues, and sometimes copy-number changes because abnormal intensity can indicate extra or missing DNA.

```r
# Plotting the two medians against each other, good samples tend to cluster together
qc <- getQC(MSet)
plotQC(qc)

# Plotting the mean detection p-value for each sample
detP <- detectionP(RGset)
barplot(colMeans(detP), las=2, cex.names=0.8, ylab="Mean detection p-values")
abline(h=0.05,col="red")
```
P-value tests whether each probe’s signal is distinguishable from background noise. For each probe in each sample, minfi compares the observed fluorescence intensity vs a background distribution (estimated from negative control probes -> technical noise -> unrelated to biological condition.) It performs a hypothesis test:

- H₀ (null): “This probe’s intensity is just background noise.”
- H₁: “This probe shows real signal.”

A low detection p-value → probe signal ≫ background → reliable measurement. A high detection p-value → probe ≈ background → unreliable, probably failed probe.

```r
# Plotting the beta-value density distribution
phenoData <- pData(MSet)
densityPlot(MSet, sampGroups = phenoData$Sample_Group)
```
The overall density distribution of Beta values for each sample is another useful metric to determine sample quality. Usually, one would expect to see most Beta values to be either close to 0 or 1, indicating most of the CpG sites in the sample are unmethylated or methylated.

```r
# Preprocessing raw data and visualizing it
mSetSq <- preprocessFunnorm(RGset)
par(mfrow=c(1,2))
plotBetasByType(MSet[,1],main="Raw")
typeI <- getProbeInfo(MSet, type = "I")[, c("Name","nCpG")]
typeII <- getProbeInfo(MSet, type = "II")[, c("Name","nCpG")]
probeTypes <- rbind(typeI, typeII)
probeTypes$Type <- rep(x = c("I", "II"), times = c(nrow(typeI), nrow(typeII)))
plotBetasByType(getBeta(mSetSq)[,1], probeTypes = probeTypes, main="Normalized",)
```
Normalization removes technical differences between Type I and Type II probes, so their beta-value distributions become comparable. In the normalized plot, it should look like the distributions are closer to each other. 

NOTE: There are various normalization methods according to the dataset type. Here we chose `preprocessFunnorm()` since our samples come from different individuals and contain real biological differences.

```r
# Removing samples with a high p-value
detP <- detP[match(featureNames(mSetSq),rownames(detP)),]
keep <- rowSums(detP < 0.01) == ncol(mSetSq)
table(keep)
mSetSqFlt <- mSetSq[keep,]
mSetSqFlt
```
The probes must be in the same order in the mSetSq and detP objects. We also get rid of values with a high p-value.

```r
# Removing SNPs
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
```
Because the presence of short nucleotide polymorphisms (or SNPs) inside the probe body can have important consequences on the downstream analysis, minfi offers the possibility to remove such probes with `dropLociWithSnps()` function.

Once normalization & filtering are done, we have clean Beta and M-values. From here, there are four different analysis goals, each using different aspects of the data.

## Analysis Goals

### 1. Differentially Methylated Positions (DMPs)
**Purpose:** to find individual CpG sites whose methylation differs.

```r
# Creating design matrix by encoding group membership
group <- factor(targets$Sample_Group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
design
```

```r
# Calculating M-vals and fitting the limma model to estimate the difference in methylation between Case and Control samples.
mVals <- getM(mSetSqFlt)
fit <- lmFit(mVals, design)
contMatrix <- makeContrasts(Case - Control, levels=design)
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
```

We only need one comparison (case vs control) to set up the contrast matrix. In methylation studies, we often have ≈800,000 CpGs while having maybe 6–10 samples; so for each CpG we estimate a variance. The main issue is, with so few samples, these variances are very noisy and unreliable. Some CpGs will look significant just because their variance was underestimated by chance. `eBayes()` pulls each CpG’s variance estimate toward a common average variance across all probes.

```r
# Getting differentially methylated positions
DMPs <- topTable(fit2, num=Inf, coef=1, sort.by="p")
head(DMPs)
```

```r
> head(DMPs)
               logFC    AveExpr         t      P.Value    adj.P.Val         B
cg11144986  3.885836  2.9584572  31.81617 3.744560e-10 0.0001700595 -3.163037
cg27586797  3.015580  2.0192165  31.48622 4.089477e-10 0.0001700595 -3.163371
cg08477332 -1.934019  1.2863126 -14.04383 3.480625e-07 0.0964937224 -3.226106
cg18697803  1.910400  1.9557880  13.38563 5.155165e-07 0.1071878702 -3.233648
cg12162195 -1.835956  0.5116746 -10.28428 4.317390e-06 0.7181486773 -3.288495
cg03609398 -1.738445 -1.5687543  -9.55069 7.747355e-06 0.8048764146 -3.308794
```

```r
# Getting genome coordinates + gene context for CpGs
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Adding gene annotation to the significant DMPs
annSub <- annEPIC[match(rownames(mVals), annEPIC$Name), ]
annSub
DMPs <- topTable(fit2, num=Inf, coef=1, genelist=annSub)
```

```r
> DMPs
             chr       pos strand       Name AddressA AddressB                                          ProbeSeqA
cg11144986 chr12 132935765      - cg11144986 77661550 21622264 ATAAATCTCATCTACAAATCCCATAAATATAAACAACAACCACCTCTACA
cg27586797  chr5  13664584      - cg27586797 79775386          AAAACTTTCTTAACATAACCCTTAAATTAAACACACACRAATACATATAC
```

```r
# Drawing per-sample beta-value plot for a single probe to visualize group differences
bVals <- getBeta(mSetSqFlt)
plotCpg(bVals, cpg="cg11144986", pheno=targets$Sample_Group)
```

### 2. Differentially Methylated Regions (DMRs)
**Purpose:** to detect regions of consistent methylation change instead of just one site (much more biologically meaningful).

```r
# Preparing limma results for the structure DMRcate needs
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
```

```
CpGannotated object describing 831693 CpG sites, with independent
CpG threshold indexed at fdr=0 and 2 significant CpG sites.
Since we have very few samples, only a few probes turned out to be significant.
```

```r
# Getting differentially methylated regions (DMR)
DMRs <- dmrcate(myAnnotation, lambda = 1000, C = 2, pcutoff = 0.05)
> DMRs
DMResults object with 2025 DMRs.
We have far more DMRs (2025) than DMPs (2) because regional smoothing boosts power.
```

```r
# Plotting the range and extracting genome coordinates for visualization
results.ranges <- extractRanges(DMRs)
results.ranges

cols <- ifelse(targets$Sample_Group == "Case", "red", "blue")

DMR.plot(
  ranges = results.ranges,
  dmr = 1,
  CpGs = mSetSqFlt,
  phen.col = cols,
  genome = "hg38"
)
```

### 3. mCSEA (region-based GSEA)
**Purpose:** Instead of the previous approach where the regions are defined according to heuristic distance rules (genomic proximity), we can define regions based on a shared function (promoters, cpg islands, etc). For this, we will use the package mCSEA. mCSEA is based on Gene Set Enrichment analysis (GSEA), a popular methodology for functional analysis.

```r
# Creating a named vector containing the rank metric (here: logFC) 
myRank <- DMPs$logFC 
names(myRank) <- rownames(DMPs)

# Reshaping the phenotype data to a format suitable for mCSEA
pheno <- as.data.frame(pData(mSetSqFlt))
pheno$Sample_Group <- targets$Sample_Group
pheno <- pheno[,"Sample_Group", drop=FALSE]
pheno

myResults <- mCSEATest(
  myRank,
  bVals,                   # Beta values
  pheno,
  regionsTypes = "promoters",   # or "cpgislands" or "geneBody"
  platform = "EPIC"             
)
```

```
767 DMRs found (padj < 0.05)

> head(myResults$promoters)
                 pval         padj   log2err         ES       NES size
BCOR     7.357473e-17 4.440347e-14 1.0574636 -0.6481874 -2.665805  143
S100A1   1.350620e-10 1.855099e-08 0.8266573 -0.8762404 -2.600072   26
NLGN4X   1.649483e-10 2.234768e-08 0.8266573 -0.9080207 -2.483241   19
RGN      8.991000e-12 1.593208e-09 0.8753251 -0.9420370 -2.476953   16
KCNQ1OT1 3.623672e-08 3.137785e-06 0.7195128 -0.7869523 -2.458419   31
RFPL2    8.723725e-09 8.601074e-07 0.7477397 -0.8778018 -2.450701   20
```

### 4. Gene Ontology / Pathway Enrichment (gometh)
**Purpose:** After obtaining a list of significantly differentially methylated CpG sites, one might wonder whether there is a (or multiple) specific biological pathways over-represented in this list. Gene-set analysis (GSA) is frequently used to discover meaningful biological patterns from lists of genes. It answers: "What biological processes are altered due to methylation change?"

```r
# Recovering the CpG IDs
anno <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Keeping only the probes present in the filtered dataset
anno_sub <- anno[rownames(anno) %in% rownames(mSetSqFlt), ]

# Getting probe coordinates
gr_probes <- GRanges(
  seqnames = anno_sub$chr,
  ranges   = IRanges(start = anno_sub$pos, end = anno_sub$pos),
  strand   = "*",
  probeID  = rownames(anno_sub)
)

# Converting to GRanges
gr_dmrs <- results.ranges 

# Matching CpGs to DMRs with genomic overlap
hits <- findOverlaps(gr_probes, gr_dmrs)
sigCpGs <- unique(mcols(gr_probes)$probeID[queryHits(hits)])
length(sigCpGs)

# Gene ontology (gometh)
allCpGs <- rownames(mSetSqFlt)

gst_GO <- gometh(
  sig.cpg = sigCpGs,
  all.cpg = allCpGs,
  collection = "GO"
)

> topGSA(gst_GO, number = 10)
           ONTOLOGY                                                             TERM  N DE         P.DE        FDR
GO:0045915       BP           positive regulation of catecholamine metabolic process  6  5 6.271948e-06 0.07060331
GO:0045964       BP                positive regulation of dopamine metabolic process  6  5 6.271948e-06 0.07060331
GO:0062237       BP                              protein localization to postsynapse 45 12 2.717284e-05 0.20392313
GO:0033240       BP                   positive regulation of amine metabolic process  8  5 6.041174e-05 0.34002746
GO:1903539       BP                    protein localization to postsynaptic membrane 44 11 1.114537e-04 0.50185391
GO:2000378       BP negative regulation of reactive oxygen species metabolic process 48 11 1.548601e-04 0.57482131
GO:0097120       BP                                 receptor localization to synapse 62 13 1.787221e-04 0.57482131
GO:1904020       BP         regulation of G protein-coupled receptor internalization  3  3 2.735938e-04 0.73259709
GO:0005587       CC                                          collagen type IV trimer  6  4 3.101745e-04 0.73259709
```
