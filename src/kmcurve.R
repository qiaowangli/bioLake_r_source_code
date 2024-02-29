library(ggplot2)
library(cowplot)
library(survminer)
library(survival)
library(jsonlite)
library(arrow)


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




KMVurve <- function(targetGene, targetClinicalGroup, title, plotType, topPre, lowPre, targetEnd, list_str, task_position, targetOutput, taskID, raw_position, sur_position){


  DS <- arrow::open_dataset(sources = task_position)                        
  SO <- Scanner$create(DS)
  AT <- SO$ToTable()
  clinicalData <- as.data.frame(AT)
  clinicalData <- clinicalData[, 1:(ncol(clinicalData) - 2)]
  clinicalData <- as.data.frame(clinicalData)


  DS <- arrow::open_dataset(sources = raw_position)                            
  SO <- Scanner$create(DS)
  AT <- SO$ToTable()
  checkData <- as.data.frame(AT)
  checkData <- checkData[, 1:(ncol(checkData) - 2)]
  checkData <- as.data.frame(checkData)



  DS <- arrow::open_dataset(sources = sur_position)              
  SO <- Scanner$create(DS)
  AT <- SO$ToTable()
  survivalData <- as.data.frame(AT)
  survivalData <- survivalData[, 1:(ncol(survivalData) - 2)]
  survivalData <- as.data.frame(survivalData)

  isNonCustom = any(list_str == "")

  rnaData = checkData[checkData['sample'] == targetGene,]
  clinicalRawData = clinicalData[clinicalData['sampleID'] == targetClinicalGroup,]

  checkRow <- unlist(clinicalRawData[1,], use.names = FALSE)[-1]
  checkRow <- Filter(nzchar, checkRow)
  numericChecker = any(is.na(as.numeric(checkRow)))

  df3= NULL

  if(isNonCustom){
    if(substr(targetClinicalGroup, 1, 16) == "Gene_expression_"){
      if(targetClinicalGroup == "Gene_expression_quantile"){

        coptRnaData = rnaData
        tmp =coptRnaData[1,]

        # Calculate the quartiles
        q1 <- quantile(tmp[,2:ncol(tmp)], probs = 0.25, na.rm=TRUE)
        q2 <- quantile(tmp[,2:ncol(tmp)], probs = 0.5, na.rm=TRUE)
        q3 <- quantile(tmp[,2:ncol(tmp)], probs = 0.75, na.rm=TRUE)
        max_value <- max(tmp[,2:ncol(tmp)], na.rm=TRUE)


        coptRnaData[2,] = ifelse(tmp > q3, "Q4",
                          ifelse(tmp > q2, "Q3", 
                          ifelse(tmp > q1, "Q2", "Q1")))

        df3 = coptRnaData[,2:ncol(coptRnaData)]
        df3[c(1, 2), ] = df3[c(2, 1), ]



      }else if(targetClinicalGroup == "Gene_expression_middle"){

        coptRnaData = rnaData
        tmp =coptRnaData[1,]

        mean <- quantile(tmp[,2:ncol(tmp)], probs = 0.5, na.rm=TRUE)


        coptRnaData[2,] = ifelse(tmp > mean, "High Expression", "Low Expression")

        df3 = coptRnaData[,2:ncol(coptRnaData)]
        df3[c(1, 2), ] = df3[c(2, 1), ]
      }
    }
    else if (!numericChecker && plotType == "NonNumeric"){

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



      clinicalRawData[2,] = ifelse(targetRow > top_threshold, "HIGH",
                        ifelse(targetRow < bottom_threshold, "LOW", "STABLE"))
    
      common <- intersect(colnames(clinicalRawData[2,]), colnames(rnaData))
      df3 <- rbind(clinicalRawData[2,common], rnaData[common])


    }else{

      common <- intersect(colnames(clinicalRawData[1,]), colnames(rnaData))
      df3 <- rbind(clinicalRawData[1,common], rnaData[common])
    }
  }else{

    GroupConvertResult = GroupConvertor(list_str, names(list_str), clinicalRawData)
    common <- intersect(colnames(GroupConvertResult[2,]), colnames(rnaData))
    df3 <- rbind(GroupConvertResult[2,common], rnaData[common])


  }


  res=as.data.frame(t(df3))
  colnames(res) <- c('Clinical_group','value')
  res$value = as.numeric(res$value)
  rownames(res)= gsub("[.-]", "*", rownames(res))
  survivalData$sample= gsub("[.-]", "*", survivalData$sample)

  mergeData <- merge(res, survivalData, by.x = "row.names", by.y = "sample")

  time = paste0(targetEnd, ".time")

  selected_column <- subset(mergeData, select = get(time))


  mergeData['duration'] = selected_column
  mergeData['targetEnd'] = mergeData[targetEnd]

  colnames(mergeData)[colnames(mergeData) == 'Clinical_group'] <- targetClinicalGroup

  formulaStr <- paste("Surv(duration, targetEnd) ~", targetClinicalGroup)

  formulaObj <- as.formula(formulaStr)

  fit <- surv_fit(formulaObj, data = mergeData)

  p = ggsurvplot(fit, data = mergeData, censor.shape="|", censor.size = 3,size = 1,
    conf.int = TRUE,       
    pval = TRUE,          
    ggtheme = theme_bw(),
    legend = "right")

  p <- p + guides(color = guide_legend(direction = "vertical"))

  file_path_out = paste0("/var/www/bioLake/static/", targetOutput)

  file_path_out_data = paste0(dirname(task_position), "/usrData/", taskID, "/",targetGene, "_", targetClinicalGroup, "_survival.csv")

  write.csv(mergeData, file = file_path_out_data)

  ggsave_workaround <- function(g){survminer:::.build_ggsurvplot(x = g)}

  g_to_save <- ggsave_workaround(p)
  ggsave(file_path_out, plot=g_to_save, dpi=300, width=10, height=6, units="in")
}
