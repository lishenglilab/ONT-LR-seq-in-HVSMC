## Prepare blast data

library(dplyr)
library(tidyverse)
library(data.table)
library(ggplot2)
library(scales)
library(ggsci)
orf <- fread("Result/ORFfind/orffind.out.fa",header = F,data.table = F,sep = "\t")
df <- data.frame(transcript_orf = orf$V1[grepl(pattern = ">",orf$V1)],
                 orf_seq = orf$V1[!grepl(pattern = ">",orf$V1)])
df$transcript_id <- str_split(str_split(df$transcript_orf,pattern = "_",simplify = T,n = 2)[,2],pattern = ":",n = 2,simplify = T)[,1]
df$orf_len <- str_count(df$orf_seq)
df <- df %>% group_by(transcript_id) %>%
  top_n(n = 1 , wt = orf_len)
df <- df[!duplicated(df$transcript_id),]
id_map <- read.table("/Data/Flair/Flair_stringent/id_map.txt",sep = "\t",header = T,check.names = F)
df <- left_join(df,id_map,by = c("transcript_id"))
write.table(df,"Result/ORFfind/transcript_orf_longest.txt",sep = "\t",row.names = F,quote = F)
conditions <- list.files(paste0(work_dir,"Result/DESeq/"),pattern = "NC$")
for (i in conditions) {
  res <- fread("Result/DESeq/",i,sep = "\t",data.table = F,check.names = F)
  tmp <- str_split(res$raw_name,pattern = "_",simplify = T)
  tmp <- as.data.frame(tmp)
  for (j in 1:nrow(tmp)) {
    tmp$transcript[j] <- ifelse(tmp$V3[j] == "",tmp$V1[j],paste(tmp$V1[j],tmp$V2[j],sep = "_"))
  }
  res$transcript_short_name <- tmp$transcript
  colnames(res)
  res <- res[,c("transcript_short_name","raw_name","Change")]
  res <- left_join(res,id_map,by = c("transcript_short_name" = "transcript_id"))
  diff_transcript <- res[res$Change != "Stable",]
  df_sub <- df[df$gene_id %in% diff_transcript$gene_id,]
  df_sub <- df_sub %>% group_by(gene_id)
  df_sub <- df_sub[str_detect(df_sub$gene_id,pattern = "ENSG"),]
  df_sub$Change <- "Stable"
  for (tr in 1:nrow(df_sub)) {
    if (sum(res$transcript_short_name == df_sub$transcript_id[tr])>0) {
      df_sub$Change[tr] <- res[res$transcript_short_name == df_sub$transcript_id[tr],]$Change
    }
  }
  filepath <- paste0("Result/","ORFfind/",i)
  dir.create(filepath)
  write.table(df_sub,file = paste0(filepath,"/transcript_orf.txt"),sep = "\t",row.names = F,quote = F)
  for (tr in intersect(diff_transcript$transcript_short_name,df_sub$transcript_id)) {
    gene =  df_sub[df_sub$transcript_id == tr,]$gene_id
    df_sub2 <- df_sub[df_sub$gene_id == gene,]
    if (nrow(df_sub2[df_sub2$Change == "Stable",])>0) {
      ref.df <- df_sub2[df_sub2$Change == "Stable",]
      ref.df <- ref.df[grep("ENST",ref.df$transcript_id),]
      if (nrow(ref.df)>0) {
        ref.df <- ref.df[ref.df$orf_len == max(ref.df$orf_len),]
        ref.fa <- data.frame(matrix(ncol = 0,nrow = 0))
        for (j in 1:nrow(ref.df)) {
          tmp <- rbind(ref.df$transcript_orf[j],ref.df$orf_seq[j])
          ref.fa <- rbind(ref.fa,tmp)
        }
        query.df <- df_sub2[df_sub2$transcript_id == tr,]
        query.fa <- data.frame(matrix(nrow = 0,ncol = 0))
        for (j in 1:nrow(query.df)) {
          tmp <- rbind(query.df$transcript_orf[j],query.df$orf_seq[j])
          query.fa <- rbind(query.fa,tmp)
        }
        dir = paste0("Data/ORFfind/",i,"/",tr)
        if (!dir.exists(dir)) {
          dir.create(dir)
        }
        setwd(dir)
        write.table(ref.fa,"ref.fa",sep = "\t",col.names = F,row.names = F,quote = F)
        write.table(query.fa,"query.fa",sep = "\t",col.names = F,row.names = F,quote = F)
      }
    }
  }
}     