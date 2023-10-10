# R version 4.0.2 (2020-06-22)
# tidyverse version=2.0.0
# DESeq2 version=1.30.1

library(DESeq2)
library(tidyverse)
library(reshape2)
library(data.table)
library(ggplot2)
exp_count = fread("Data/Flair/counts_matrix.tsv",sep = "\t",header = T,check.names = F,data.table = F)
rownames(exp_count) <- exp_count$ids
exp_count$ids = NULL
group_list <- c(rep("221",4),rep("NC",4),rep("PDGF",4),rep("TGF",4))
condition <- relevel(factor(group_list),ref = "NC") #set `NC` as reference level
colData <- data.frame(row.names = colnames(exp_count),
                      condition = condition)
dds <- DESeqDataSetFromMatrix(countData = exp_count, colData = colData, design = ~ condition)
dds <- DESeq(dds)
resultsNames(dds)
compare <- resultsNames(dds)[2:4]
for (i in compare) {
  res <- results(dds, name = i)
  res$padj <- ifelse(is.na(res$padj), 1, res$padj)
  resordered <- res[order(res$padj),]
  deg <- as.data.frame(resordered)
  FC_cutoff <- 1.5
  deg$Change <- ifelse(abs(deg$log2FoldChange)>log2(FC_cutoff) & deg$padj<0.05,ifelse(deg$log2FoldChange<0,"Down","Up"),"Stable")
  print(i)
  table(deg$Change)
  deg$transcript_id <- rownames(deg)
  deg <- dplyr::select(deg,transcript_id,everything())
  write.table(deg,paste0("Result/DESeq/",i),sep = "\t",row.names = F,quote = F)
  # volcano plot
  deg_sig=deg %>%
    arrange(padj) %>% 
    filter(Change != "Stable")

  deg_stable=deg %>%
    filter(Change == "Stable")
  

  deg_top_up=deg_sig %>%
    filter(Change == "Up")
  
  deg_top_down=deg_sig %>%
    filter(Change == "Down")
  deg_top=rbind(deg_top_up[1:5,],deg_top_down[1:5,])

  index=deg_sig$transcript_id %in% deg_top$transcript_id
  deg_sig=deg_sig[!index,]
  library(ggrepel)
  p <- ggplot(deg_stable, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(alpha = 0.6, size = 1.5,color="gray") + 
    geom_point(data=deg_sig,aes(color=Change),shape=19,alpha = 0.8,size = 4)+ 
    geom_point(data=deg_top,aes(fill=Change),shape=21,alpha = 0.8,size = 6)+ 
    scale_color_manual(values=c("#377EB8","#E41A1C"))+ 
    scale_fill_manual(values=c("Down" = "#377EB8","Up" = "#E41A1C"))+ 
    geom_text_repel(data=deg_top,(aes(label=transcript_id)),force=5,color="black",hjust = 0.5,vjust=1)+ 
    guides(color="none",fill="none")+ 
    theme_test()+ 
    geom_vline(xintercept=c(-log2(FC_cutoff),log2(FC_cutoff)),lty=4,col="black",lwd=0.8) + 
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) + 
    labs(x="log2(FoldChange)",y="-log10 adjusted P value")+ 
    theme(axis.text = element_text(size = 12, colour = "black"), 
          axis.title= element_text(size = 13, colour = "black"))
  
  ggsave(plot = p,filename = paste0(paste0("Figure/volcano/",i),".pdf"), width = 15, height = 10)
  gene_type <- read.table("reference/gencode_v38_genetype.txt",sep = "\t",header = T)
  deg_sig=deg %>%
    arrange(padj) %>% 
    filter(Change != "Stable")
  tmp <- deg_sig[str_detect(pattern = "ENSG",string = deg_sig$transcript_id) & !str_detect(pattern = "PAR",string = deg_sig$transcript_id),]
  tmp <- as.data.frame(str_split(tmp$transcript_id,pattern = "_",simplify = T))
  gene_select <- tmp$V2
  gene_select <- data.frame(Gene_id=unique(gene_select))
  gene_select <- left_join(gene_select,gene_type[,1:2],by = "Gene_id")
  gene_select <- na.omit(gene_select)
  gg <- enrichGO(gene = gene_select$Gene_name,
                 keyType = "SYMBOL",
                 OrgDb = org.Hs.eg.db,
                 ont = "ALL",
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH")
  p <- dotplot(gg)
  ggsave(plot = p,filename = paste0(paste0("Figure/GO/",i),".pdf"),width = 10,height = 8)
  gg <- as.data.frame(gg)
  write.table(gg,paste0("Result/Enrichment/GO/",i),sep = "\t",quote = F,row.names = F)
  gene_select <- bitr(gene_select$Gene_name,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
  kk <- enrichKEGG(gene = gene_select$ENTREZID,
                   organism = "hsa",
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH")
  p <- dotplot(kk)
  ggsave(plot = p,filename = paste0(paste0("Figure/KEGG/",i),".pdf"),width = 10,height = 8)
  kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  kk <- as.data.frame(kk)
  write.table(kk,paste0("Result/Enrichment/KEGG/",i),sep = "\t",quote = F,row.names = F)
}
