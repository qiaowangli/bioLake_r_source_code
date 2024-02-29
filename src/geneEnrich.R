library(R.utils)
library(clusterProfiler)
library(arrow)

options(download.file.method="auto", download.file.extra="-k -L")
R.utils::setOption("clusterProfiler.download.method","auto")

GO_database <- 'org.Hs.eg.db' #http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
KEGG_database <- 'hsa' ##http://www.genome.jp/kegg/catalog/org_list.html

args <- commandArgs(trailingOnly = FALSE)

dataPath = args[7]
geneID = args[8]
TopDown = args[9] # high/low
NumSample = as.numeric(args[10])
taskID = args[11]
task_position = args[12]


col_names <- c("gene","value")
data <- read.csv(dataPath, col.names = col_names)


if (!is.na(NumSample) && NumSample < nrow(data)) {
    if (TopDown == "low") {
        data = tail(data, NumSample)
    } else {
        data = head(data, NumSample)
    }
}


gene <- bitr(data$gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)

KEGG<-enrichKEGG(gene$ENTREZID,
                 organism = KEGG_database,
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)


KEGG_pathways_id <-KEGG@result$ID[KEGG@result$`pvalue`<0.05&KEGG@result$`qvalue`<0.05]
KEGG_pathways_des <-KEGG@result$Description[KEGG@result$`pvalue`<0.05&KEGG@result$`qvalue`<0.05]
first_column <- c(KEGG_pathways_id)
second_column <- c(KEGG_pathways_des)
df=data.frame(first_column, second_column)
names(df) =c('PathwayID', 'PathwayName')


mainPart = paste(task_position,'/', taskID, '/', 'KEGG_pathway_', geneID, '.csv', sep = "")
write.csv(df, mainPart)

mainPartstandby = paste("/var/www/bioLake/static/mRNAStage/tmp/KEGG_pathway_", geneID, '.csv', sep = "")
write.csv(df, mainPartstandby)


mainPartImage = paste("/var/www/bioLake/static/mRNAStage/tmp/KEGG_pathway_result_", geneID,sep = "")
FullPartImage = paste(mainPartImage, ".jpg",sep = "")

jpeg(FullPartImage, quality = 100, res = 130, width = 800, height = 600)
barplot(KEGG,showCategory = 10,title = 'KEGG Pathway')
dev.off()