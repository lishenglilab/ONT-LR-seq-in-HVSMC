# R version 4.0.2 (2020-06-22)
# IsoformSwitchAnalyzeR version=1.17.04

setwd("/Data/Flair/")
library(IsoformSwitchAnalyzeR)
library(dplyr)
library(data.table)

isoformCountMatrix <- fread("counts_matrix.tsv",sep = "\t",header = T,check.names = F,data.table = F)
colnames(isoformCountMatrix) <- c("isoform_id","m221-1","221-2","221-3","221-4","NC-1","NC-2","NC-3","NC-4","PDGF-1","PDGF-2","PDGF-3","PDGF-4","TGF-1","TGF-2","TGF-3","TGF-4")
isoformRepExpression <- fread("counts_matrix.tsv.tpm.tsv",sep = "\t",header = T,check.names = F,data.table = F)
colnames(isoformRepExpression) <- c("isoform_id","m221-1","221-2","221-3","221-4","NC-1","NC-2","NC-3","NC-4","PDGF-1","PDGF-2","PDGF-3","PDGF-4","TGF-1","TGF-2","TGF-3","TGF-4")
isoformCountMatrix$isoform_id <- gsub('[_][^_]+$',"",isoformCountMatrix$isoform_id)
isoformRepExpression$isoform_id <- gsub('[_][^_]+$',"",isoformRepExpression$isoform_id)
designMatrix <- data.frame(
  sampleID = colnames(isoformCountMatrix[-1]),
  condition = rep(c("a221","NC","PDGF","TGF"),each = 4)
)
isoformExonAnnoation <- "flair.collapse.isoforms.gtf"
isoformNtFasta <- "flair.collapse.isoforms.rename.fa"
comparisonsToMake <- data.frame(
  condition_1 = rep("NC",3),
  condition_2 = c("a221","PDGF","TGF")
)
### Create SwitchAnalyzeRlist
aSwitchList <- importRdata(
  isoformCountMatrix = isoformCountMatrix,
  isoformRepExpression = isoformRepExpression,
  designMatrix = designMatrix,
  isoformExonAnnoation = isoformExonAnnoation,
  isoformNtFasta = isoformNtFasta,
  comparisonsToMake = comparisonsToMake
)
summary(aSwitchList)
### Filter
aSwitchList <- preFilter( aSwitchList )
### Test for isoform switches
aSwitchList <- isoformSwitchTestDEXSeq( aSwitchList )
### If analysing (some) novel isoforms (else use CDS from ORF as explained in importRdata() )
aSwitchList <- addORFfromGTF( aSwitchList,pathToGTF = "~/project/pinghuaji/full_lengh_pinghuaji/reference/gencode.v38.primary_assembly.annotation.gtf")
aSwitchList <- analyzeNovelIsoformORF( aSwitchList,analysisAllIsoformsWithoutORF = T)
### Extract Sequences
aSwitchList <- extractSequence( aSwitchList,
                                pathToOutput = "isoform_switch/",
                                removeLongAAseq = T,
                                alsoSplitFastaFile = T)
### Summary
extractSwitchSummary(aSwitchList)
#Part2
### Add annotation
aSwitchList <- analyzeCPAT( aSwitchList,
                            pathToCPATresultFile = "isoform_switch/extdata/capt_result.txt",
                            codingCutoff = 0.725,removeNoncodinORFs = T)

aSwitchList <- analyzePFAM( aSwitchList,
                            pathToPFAMresultFile = "isoform_switch/extdata/pfam_results.txt")
aSwitchList <- analyzeSignalP( aSwitchList,
                               pathToSignalPresultFile = "isoform_switch/extdata/signalP_results.txt" )

aSwitchList <- analyzeIUPred2A(aSwitchList,
                pathToIUPred2AresultFile = "isoform_switch/extdata/iupred2a_result.txt") 

aSwitchList <- analyzeAlternativeSplicing( aSwitchList )

### Analyse consequences
aSwitchList <- analyzeSwitchConsequences( aSwitchList )

gene_id <- data.frame(gene_id = aSwitchList$isoformFeatures$gene_id)
id_map <- fread("~/project/pinghuaji/full_lengh_pinghuaji/reference/gencode_v38_genetype.txt",sep = "\t",header = T,check.names = F)
gene_id <- left_join(gene_id,id_map[,1:2],by = c("gene_id"="Gene_id"))
aSwitchList$isoformFeatures$gene_name <- gene_id$Gene_name
#### gene_q_value
gene_exp <- fread("gene.count.txt",sep = "\t",header = T,check.names = F,data.table = F)
rownames(gene_exp) <- gene_exp$Geneid
gene_exp$Geneid <- NULL
group_list <- c(rep(c("a221","NC","PDGF","TGF"),each = 4))
condition <- relevel(factor(group_list),ref = "NC")
colData <- data.frame(row.names = colnames(gene_exp), 
                      condition=condition)
dds <- DESeqDataSetFromMatrix(countData = gene_exp, colData = colData, design = ~ condition)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, name = "condition_a221_vs_NC")
summary(res)
deg <- as.data.frame(res)
deg <- na.omit(deg)
deg$gene_id <- rownames(deg)
deg$tmp <- paste(deg$gene_id,"a221",sep = "_")
aSwitchList$isoformFeatures$tmp <- paste(aSwitchList$isoformFeatures$gene_id,aSwitchList$isoformFeatures$condition_2,sep = "_")
gene <- data.frame(tmp = aSwitchList$isoformFeatures$tmp)
gene <- left_join(gene,deg[,c("tmp","padj")],by = "tmp")
res <- results(dds, name = "condition_PDGF_vs_NC")
summary(res)
deg <- as.data.frame(res)
deg <- na.omit(deg)
deg$gene_id <- rownames(deg)
deg$tmp <- paste(deg$gene_id,"PDGF",sep = "_")
gene <- left_join(gene,deg[,c("tmp","padj")],by = "tmp")
res <- results(dds, name = "condition_TGF_vs_NC")
summary(res)
deg <- as.data.frame(res)
deg <- na.omit(deg)
deg$gene_id <- rownames(deg)
deg$tmp <- paste(deg$gene_id,"TGF",sep = "_")
gene <- left_join(gene,deg[,c("tmp","padj")],by = "tmp")
gene[is.na(gene)] <- 0
for (i in 1:nrow(gene)) {
  gene$qvalue[i] <- sum(gene$padj.x[i],gene$padj.y[i],gene$padj[i])
}
identical(aSwitchList$isoformFeatures$tmp,gene$tmp)
aSwitchList$isoformFeatures$gene_q_value <- gene$qvalue
#### iso_q_value
iso_exp <- fread("counts_matrix.tsv",sep = "\t",header = T,check.names = F,data.table = F)
rownames(iso_exp) <- iso_exp$ids
iso_exp$ids <- NULL
group_list <- c(rep(c("a221","NC","PDGF","TGF"),each = 4))
condition <- relevel(factor(group_list),ref = "NC")
colData <- data.frame(row.names = colnames(iso_exp), 
                      condition=condition)
dds <- DESeqDataSetFromMatrix(countData = iso_exp, colData = colData, design = ~ condition)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, name = "condition_a221_vs_NC")
summary(res)
deg <- as.data.frame(res)
deg <- na.omit(deg)
rownames(deg) <- gsub('[_][^_]+$',"",rownames(deg))
deg$gene_id <- rownames(deg)
deg$tmp <- paste(deg$gene_id,"a221",sep = "_")
aSwitchList$isoformFeatures$tmp <- paste(aSwitchList$isoformFeatures$isoform_id,aSwitchList$isoformFeatures$condition_2,sep = "_")
gene <- data.frame(tmp = aSwitchList$isoformFeatures$tmp)
gene <- left_join(gene,deg[,c("tmp","padj")],by = "tmp")
res <- results(dds, name = "condition_PDGF_vs_NC")
summary(res)
deg <- as.data.frame(res)
deg <- na.omit(deg)
rownames(deg) <- gsub('[_][^_]+$',"",rownames(deg))
deg$gene_id <- rownames(deg)
deg$tmp <- paste(deg$gene_id,"PDGF",sep = "_")
gene <- left_join(gene,deg[,c("tmp","padj")],by = "tmp")
res <- results(dds, name = "condition_TGF_vs_NC")
summary(res)
deg <- as.data.frame(res)
deg <- na.omit(deg)
rownames(deg) <- gsub('[_][^_]+$',"",rownames(deg))
deg$gene_id <- rownames(deg)
deg$tmp <- paste(deg$gene_id,"TGF",sep = "_")
gene <- left_join(gene,deg[,c("tmp","padj")],by = "tmp")
gene[is.na(gene)] <- 0
for (i in 1:nrow(gene)) {
  gene$qvalue[i] <- sum(gene$padj.x[i],gene$padj.y[i],gene$padj[i])
}
identical(aSwitchList$isoformFeatures$tmp,gene$tmp)
aSwitchList$isoformFeatures$iso_q_value <- gene$qvalue
aSwitchList$isoformFeatures$tmp <- NULL
saveRDS(aSwitchList,"isoform_switch/Result/iso_switch.RDS")
iso_switch <- aSwitchList$isoformSwitchAnalysis
gene_type <- read.table("~/project/pinghuaji/full_lengh_pinghuaji/reference/gencode_v38_genetype.txt",sep = "\t",header = T)
gene_type <- gene_type[,c("gene_id","gene_name")]
gene_type <- unique(gene_type)
iso_switch <- left_join(iso_switch,gene_type)
iso_switch$gene_name[is.na(iso_switch$gene_name)] <- "Novel_gene"
write.table(iso_switch,"isoform_switch/Result/iso_switch_v3.txt",quote = F,row.names = F,sep = "\t")
### Visual analysis
# Indiviudal switches
switchPlotTopSwitches( aSwitchList ,pathToOutput = "isoform_switch/Result/test")
switchPlot(aSwitchList,gene = 'PAIP1',condition1 = 'a221',condition2 = 'NC')
### Summary
extractSwitchSummary(aSwitchList, filterForConsequences = TRUE)

extractSwitchOverlap(
  aSwitchList,
  filterForConsequences=TRUE,
  plotIsoforms = FALSE
)
extractConsequenceSummary(
  aSwitchList,
  consequencesToAnalyze='all',
  plotGenes = FALSE,           # enables analysis of genes (instead of isoforms)
  asFractionTotal = FALSE      # enables analysis of fraction of significant features
)
# Consequence Enrichment Analysis
extractConsequenceEnrichment(
  aSwitchList,
  consequencesToAnalyze='all',
  analysisOppositeConsequence = TRUE,
  returnResult = FALSE # if TRUE returns a data.frame with the summary statistics
)

extractConsequenceEnrichmentComparison(
  aSwitchList,
  consequencesToAnalyze=c('domains_identified','intron_retention','coding_potential'),
  analysisOppositeConsequence = TRUE,
  returnResult = FALSE # if TRUE returns a data.frame with the summary statistics
)

# Splicing Enrichment Analysis
extractSplicingEnrichment(
  aSwitchList,
  returnResult = FALSE # if TRUE returns a data.frame with the summary statistics
)

extractSplicingEnrichmentComparison(
  aSwitchList,
  splicingToAnalyze = c('A3','A5','ATSS','ATTS'), # the splicing highlighted above
  returnResult = FALSE # if TRUE returns a data.frame with the results
)

#Overview Plots
### Volcano like plot:
pdf("isoform_switch/Result/overview_volvano.pdf",width = 10,height = 5)
ggplot(data=aSwitchList$isoformFeatures, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
  geom_point(
    aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
    size=1
  ) +
  geom_hline(yintercept = -log10(0.05), linetype='dashed') + # default cutoff
  geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff
  facet_wrap( ~ condition_2) +
  #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
  scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
  labs(x='dIF', y='-Log10 ( Isoform Switch Q Value )') +
  theme_bw()
dev.off()
