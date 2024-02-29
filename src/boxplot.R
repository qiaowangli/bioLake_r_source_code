library(arrow)
library(ggplot2)
library(cowplot)
library(dplyr)

GroupConvertor <- function(groups, group_names, clinicalRawData) {
  clinicalRawData[2, ] = "other"
  
  for (i in seq_along(groups)) {
    group <- groups[[i]]
    group_name <- group_names[i]
    
    for (j in seq_along(group)) {
      value <- group[[j]]
      
      matching_col <- clinicalRawData[1, ] %in% value
      clinicalRawData[2, matching_col] = group_name
    }
  }
  
  return(clinicalRawData)
}


boxplotMaker <- function(targetGene, targetClinicalGroup, title, plotType, topPre, lowPre, list_str, task_position, targetOutput, taskID, raw_position){

  print(task_position)
  print(raw_position)

  DS <- arrow::open_dataset(sources = task_position)                        
  SO <- Scanner$create(DS)
  AT <- SO$ToTable()
  data <- as.data.frame(AT)
  data <- data[, 1:(ncol(data) - 2)]
  data <- as.data.frame(data)

  DS <- arrow::open_dataset(sources = raw_position)                            
  SO <- Scanner$create(DS)
  AT <- SO$ToTable()
  checkData <- as.data.frame(AT)
  checkData <- checkData[, 1:(ncol(checkData) - 2)]
  checkData <- as.data.frame(checkData)


rnaData = checkData[checkData['sample'] == targetGene,]
clinicalRawData = data[data['sampleID'] == targetClinicalGroup,]

targetRow <- unlist(clinicalRawData[1,], use.names = FALSE)[-1]
targetRow <- Filter(nzchar, targetRow)
numericChecker = any(is.na(as.numeric(targetRow)))
isNonCustom = any(list_str == "")


if(isNonCustom){
if (!numericChecker && plotType == "NonNumeric"){

  targetRow = clinicalRawData[1 , colSums(is.na(clinicalRawData))==0]
  targetRow = as.numeric(targetRow)


  if(is.na(as.numeric(lowPre)) || as.numeric(lowPre) >= 1){lowPre = 0.25}
  if(is.na(as.numeric(topPre)) || as.numeric(lowPre) >= 1){topPre = 0.25}

  lowPre = as.numeric(lowPre)
  topPre = as.numeric(topPre)

  q1 <- quantile(targetRow, probs = lowPre, na.rm=TRUE)
  q3 <- quantile(targetRow, probs = 1 - topPre, na.rm=TRUE)

  top_threshold <- q3
  bottom_threshold <- q1

  clinicalRawData[2,] <- ifelse(targetRow > top_threshold, "HIGH",
                    ifelse(targetRow < bottom_threshold, "LOW", "STABLE"))

  common <- intersect(colnames(clinicalRawData[2,]), colnames(rnaData))
  df3 <- rbind(clinicalRawData[2,common], rnaData[common])


  res=as.data.frame(t(df3))
  colnames(res) <- c('Clinical_group','C2')
  res$C2 = as.numeric(res$C2)

  tryCatch(
    expr = {
      res <- subset(res, res$Clinical_group != "")
      res$Clinical_group <- factor(res$Clinical_group, levels = unique(res$Clinical_group))
       

    },
    error = function(e){ 
        res$Clinical_group <- ifelse(res$Clinical_group == "", "N/A", res$Clinical_group)
    },
    warning = function(w){
        res$Clinical_group <- ifelse(res$Clinical_group == "", "N/A", res$Clinical_group)
    }
)



}else{
  common <- intersect(colnames(clinicalRawData[1,]), colnames(rnaData))
  df3 <- rbind(clinicalRawData[1,common], rnaData[common])

  res=as.data.frame(t(df3))
  colnames(res) <- c('Clinical_group','C2')
  res$C2 = as.numeric(res$C2)
  res <- res %>% arrange(as.numeric(Clinical_group))
  res$Clinical_group <- factor(res$Clinical_group, levels = unique(res$Clinical_group))

}
}else{

  GroupConvertResult = GroupConvertor(list_str, names(list_str), clinicalRawData)
  common <- intersect(colnames(GroupConvertResult[2,]), colnames(rnaData))
  df3 <- rbind(GroupConvertResult[2,common], rnaData[common])

  res=as.data.frame(t(df3))
  colnames(res) <- c('Clinical_group','C2')
  res$C2 = as.numeric(res$C2)

}



give.n <- function(x){
  return(c(y = median(x)*1.007, label = length(x))) 
}

mean.n <- function(x){
  return(c(y = median(x)*0.993, label = round(mean(x),2))) 
}

Group <- targetClinicalGroup
colnames(res) <- c(targetClinicalGroup,'C2')

p = ggplot(res, aes(x = .data[[Group]], y = C2, fill = .data[[Group]])) +
  geom_boxplot() +
  theme_classic() +
  stat_summary(fun.data = give.n, geom = "text") +
  stat_summary(fun.data = mean.n, geom = "text", colour = "red")+
  xlab(targetClinicalGroup) +
  ylab("Gene expression") +
  background_grid(major = "xy", minor = "none")


file_path_out = paste0("/var/www/bioLake/static/",targetOutput)

file_path_out_data = paste0(dirname(task_position), "/usrData/", taskID, "/",targetGene, "_", targetClinicalGroup, "_expression.csv")


write.csv(res, file = file_path_out_data)
ggsave(file_path_out, plot=p)
}