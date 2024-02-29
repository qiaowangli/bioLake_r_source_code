library(enrichplot)
library(clusterProfiler)
library(R.utils)

options(download.file.method="auto", download.file.extra="-k -L")
R.utils::setOption("clusterProfiler.download.method","curl")

GO_database <- 'org.Hs.eg.db' #http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
KEGG_database <- 'hsa' ##http://www.genome.jp/kegg/catalog/org_list.html


args <- commandArgs(trailingOnly = FALSE)
dataPath = args[7]
pathwayID = args[9]

data <- read.csv(dataPath)
colnames(data) = c("SYMBOL","value")

gene <- bitr(data$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)

info_merge <- merge(data, gene, by = 'SYMBOL') 
GSEA_input <- info_merge$value
names(GSEA_input) = info_merge$ENTREZID
GSEA_input = sort(GSEA_input, decreasing = TRUE)

GSEA_KEGG <- gseKEGG(GSEA_input, organism = KEGG_database, pvalueCutoff = -log10(0.05))

mainPart = paste("/var/www/bioLake/static/mRNAStage/tmp/KEGG_pathway_result_", pathwayID,sep = "")

FullPart = paste(mainPart, ".jpg",sep = "")
jpeg(FullPart, quality = 100, res = 130, width = 800, height = 600)
gseaplot2(GSEA_KEGG,pathwayID,pvalue_table = FALSE, rel_heights=c(1, .2, .6))
ridgeplot(GSEA_KEGG) 
dev.off()
