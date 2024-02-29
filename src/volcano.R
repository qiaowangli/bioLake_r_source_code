library(dplyr)
library(ggplot2)
library(arrow)
library(dplyr)



args <- commandArgs(trailingOnly = FALSE)

targetGene = args[7]
CutOffPvalue = as.numeric(args[8])
CutOfflogFC = as.numeric(args[9])
DATA_FLAG = args[10]
customTitle = args[11]
dataset = args[12]
geneTrackId = args[13]

print(dataset)
print(DATA_FLAG)




file_path_out = paste0("/var/www/bioLake/static/mRNAStage/tmp/", geneTrackId, "_volcanoResult.png")
file_path_raw = paste0(dirname(dataset), "/tTest/",targetGene,"_tTest.csv")
file_path_rank_p = paste0(dirname(dataset), "/tTest/",targetGene,"_pRank.csv")
file_path_rank_f = paste0(dirname(dataset), "/tTest/",targetGene,"_fRank.csv")
Raw_data = dataset

if(!file.exists(file_path_out)){

if(DATA_FLAG == 'False'){
  print(Raw_data)

  DS <- arrow::open_dataset(sources = Raw_data)              
  SO <- Scanner$create(DS)
  AT <- SO$ToTable()
  STEP_ONE_REL_unsort <- as.data.frame(AT)

  row_index <- which(STEP_ONE_REL_unsort$sample == targetGene)[1]
  targetRow = as.list(STEP_ONE_REL_unsort[row_index,])
  orderIndexList = order(unlist(targetRow)) 
  last_element = orderIndexList[length(orderIndexList)]
  orderIndexList = orderIndexList[-length(orderIndexList)]
  orderIndexList = c(last_element, orderIndexList)
  sorted_data <- STEP_ONE_REL_unsort[, orderIndexList]


  total_columns <- ncol(sorted_data)
  mid_column <- (total_columns + 1) / 2
  group1 <- 2:floor(mid_column) 
  group2 <- (floor(mid_column) + 1):(total_columns)

columns <- c("gene","pValue","foldChange") 
df <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(df) <- columns



for(i in 1:nrow(sorted_data)){
  print(i)
  
  if(sorted_data[i,'mean'] <= 0.2){
    next
  }

  geneName <- sorted_data[i,]$sample
  
  groupA <- as.vector(sorted_data[i, group1], 'numeric')
  groupB <- as.vector(sorted_data[i, group2], 'numeric')
  
  T_Test_Ans <- t.test(groupA, groupB, alternative = "two.sided", var.equal = FALSE)
  fold_change <- mean(groupA) / mean(groupB)
  df[nrow(df) + 1,] <- c(geneName, T_Test_Ans$p.value, fold_change)
}


#############################################################################################################################
STEP_ONE_REL = df

STEP_ONE_REL = subset(STEP_ONE_REL,STEP_ONE_REL$foldChange != "Inf")
STEP_ONE_REL$pValue = as.double(STEP_ONE_REL$pValue)
STEP_ONE_REL$foldChange = as.double(STEP_ONE_REL$foldChange)

write.csv(STEP_ONE_REL, file_path_raw)

STEP_ONE_REL$pValue_log <- log10(STEP_ONE_REL$pValue)

sorted_by_pValue <- STEP_ONE_REL %>% arrange(desc(pValue_log))

output_data <- data.frame(gene = sorted_by_pValue$gene, value = sorted_by_pValue$pValue_log)
write.csv(output_data, file_path_rank_p, row.names = FALSE, quote = FALSE)

sorted_by_foldChange_desc <- STEP_ONE_REL %>% arrange(desc(foldChange))
output_data <- data.frame(gene = sorted_by_pValue$gene, value = sorted_by_pValue$foldChange)

write.csv(output_data, file_path_rank_f, row.names = FALSE, quote = FALSE)


}else{

  STEP_ONE_REL <- read.csv(file_path_raw, row.names = NULL)
}

STEP_ONE_REL$pValueAdjusted <- p.adjust(STEP_ONE_REL$pValue, method = "BH")

STEP_ONE_REL$pValueAdjusted = STEP_ONE_REL$pValue

cut_off_pvalue = CutOffPvalue
cut_off_logFC = CutOfflogFC
STEP_ONE_REL$change = ifelse(-log10(STEP_ONE_REL$pValue) > -log10(cut_off_pvalue), ifelse(log2(STEP_ONE_REL$foldChange)> cut_off_logFC ,'Up',ifelse(log2(STEP_ONE_REL$foldChange) < - cut_off_logFC ,'Down','Stable')), 'Stable')

p <- ggplot(
  STEP_ONE_REL, aes(x = log2(foldChange), y = -log10(pValue), colour=change)) +
  geom_point(alpha=0.4, size=1.5) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  
  geom_vline(xintercept=c(-cut_off_logFC,cut_off_logFC),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.2), 
        legend.position="right", 
        legend.title = element_blank())


ggsave(file_path_out, plot=p, dpi=300, width=10, height=8, units="in")

}