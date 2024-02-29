library(ggplot2)
library(arrow)
library(pheatmap)
library(reshape2)
library(dplyr)


args <- commandArgs(trailingOnly = FALSE)

targetGene = args[7]
numSample = as.numeric(args[8])
numGene = as.numeric(args[9])
correlationMethod = args[10]
dataset = args[11]
geneTrackId = args[12]
taskID = args[13]
correlationMethod <- ifelse(args[10] %in% c('pearson', 'kendall', 'spearman'),'spearman')


heatmapUtility <- function(geneID) {

    DS <- arrow::open_dataset(sources = dataset)
                                  
    SO <- Scanner$create(DS)
    AT <- SO$ToTable()
    raw_data <- as.data.frame(AT)

   
    raw_data <- raw_data[, 1:(ncol(raw_data) - 2)]

    base_list <- raw_data[raw_data$sample == as.character(geneID), ][1, 2:ncol(raw_data)]
    base_list <- t(raw_data[raw_data$sample == as.character(geneID), 2:ncol(raw_data)])
    statList <- cor(t(raw_data[, 2:ncol(raw_data)]), base_list, method = correlationMethod)
    raw_data$correlation <- statList

    print('ok')
   


    return(raw_data)
    
}

heatmapMaker <- function(targetGene, numSample, numGene){

  numSample = as.numeric(numSample)
  numGene = as.numeric(numGene)

  full_data = heatmapUtility(targetGene)

  file_path_out_positive = paste0("/var/www/bioLake/static/mRNAStage/tmp/",geneTrackId, "_", numGene, "_rawDataHeatMapresult_positive",".png")
  file_path_out_negative = paste0("/var/www/bioLake/static/mRNAStage/tmp/",geneTrackId, "_", numGene, "_rawDataHeatMapresult_negative",".png")
  file_path_out_positive_data = paste0(dirname(dataset), "/usrData/", taskID, "/",targetGene, "_", numGene, "_heatmap_p",".csv")
  file_path_out_negative_data = paste0(dirname(dataset), "/usrData/", taskID, "/",targetGene, "_", numGene, "_heatmap_n",".csv")
  file_path_out_ranked_data = paste0(dirname(dataset), "/usrData/", taskID, "/",targetGene, "_heatmap_rank",".csv")



  numGene = ifelse(nrow(full_data) > numGene, numGene, nrow(full_data))
  numSample = ifelse(ncol(full_data) > numSample, numSample, ncol(full_data)-5)
  

  if(!file.exists(file_path_out_negative)){
    target_gene = full_data[full_data$sample == targetGene,]
    full_data = full_data[order(full_data$correlation,decreasing=FALSE),]
    Top_Data = full_data[1:numGene,]
    Top_Data <- subset(Top_Data, select = -c(mean, correlation, sd))
    target_gene <- subset(target_gene, select = -c(mean, correlation, sd))
    target_gene = target_gene[,order(as.numeric(target_gene[1,]),decreasing = F)]
    Sort_Top_Data = rbind(target_gene,Top_Data)
    if (length(unique(Sort_Top_Data[, ncol(Sort_Top_Data)])) < nrow(Sort_Top_Data)) {
      duplicated_rows = duplicated(Sort_Top_Data[, ncol(Sort_Top_Data)])
      second_occurrence = duplicated_rows & !duplicated(duplicated_rows)
      Sort_Top_Data = Sort_Top_Data[!second_occurrence, ]
    }
    rownames(Sort_Top_Data) <- Sort_Top_Data$sample
    colnames(Sort_Top_Data) <- colnames(Sort_Top_Data)
    geneExp_matrix <- as.matrix(Sort_Top_Data[1:numSample])
    rownames(geneExp_matrix) = Sort_Top_Data$sample
    min_bar = geneExp_matrix[1,1]
    geneExp_matrix[geneExp_matrix <= min_bar] <- min_bar
    write.csv(geneExp_matrix, file = file_path_out_positive_data)
    result = pheatmap(geneExp_matrix, fontsize=5, color=colorRampPalette(rev(c("red","white","blue")))(100), cluster_rows = F, cluster_cols = F, show_colnames = F,
                      main = paste(targetGene,"Heat Map of Negatively Correlated Genes"), legend = TRUE)
    ggsave(file_path_out_negative, plot=result, dpi=300, width=10, height=8, units="in")
  }


  if(!file.exists(file_path_out_positive)){
    target_gene = full_data[full_data$sample == targetGene,]
    full_data = full_data[order(full_data$correlation,decreasing=TRUE),]

    raw_data_rank <- data.frame(gene = full_data$sample, rank = full_data$correlation)
    output_data <- raw_data_rank %>% arrange(desc(full_data$correlation))
    write.csv(output_data, file = file_path_out_ranked_data, row.names = FALSE, quote = FALSE)

    Top_Data = full_data[1:numGene,]
    Top_Data = subset(Top_Data, select = -c(mean))
    Top_Data = subset(Top_Data, select = -c(correlation) )
    Top_Data = subset(Top_Data, select = -c(sd) )
    Sort_Top_Data = Top_Data[,order(as.numeric(Top_Data[1,]),decreasing = F)]
    rownames(Sort_Top_Data) <- Sort_Top_Data$sample
    colnames(Sort_Top_Data) <- colnames(Sort_Top_Data)
    geneExp_matrix <- as.matrix(Sort_Top_Data[1:numSample])
    rownames(geneExp_matrix) = Sort_Top_Data$sample
    min_bar = geneExp_matrix[1,1]
    geneExp_matrix[geneExp_matrix <= min_bar] <- min_bar
    write.csv(geneExp_matrix, file = file_path_out_negative_data)
    result = pheatmap(geneExp_matrix, fontsize=5, color=colorRampPalette(rev(c("red","white","blue")))(100), cluster_rows = F, cluster_cols = F, show_colnames = F,
          main = paste(targetGene,"Heat Map of Positively Correlated Genes"), legend = TRUE)

    ggsave(file_path_out_positive, plot=result, dpi=300, width=10, height=8, units="in")

  }
}

heatmapMaker(targetGene, numSample, numGene)


