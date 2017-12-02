library(mygene)
library(plyr)

setwd("/Users/joy/macpro_mount/data_alpha/home/jpai/lillie/go_analysis/")
genes <- scan("deg_pca_up_genes.txt", character())

hiv_genes <- read.table("hiv_interaction_genes.txt", stringsAsFactors = FALSE, col.names = "gene_symbol")

selected_go_table <- read.csv("interested_categories.csv", header=TRUE, stringsAsFactors = FALSE)
selected_go_ids <- selected_go_table[, 'GO.id']

gene_annots <- queryMany(genes, scopes='symbol', fields=c('go.BP', 'go.MF', 'type_of_gene', 'alias'), species='human', size=1, returnall=TRUE)

gene_annots$go.BP <- lapply(seq_along(gene_annots$response$go.BP), clean_mygene_output, "BP")
gene_annots$go.MF <- lapply(seq_along(gene_annots$response$go.MF), clean_mygene_output, "MF")



go_df <- do.call(rbind, c(gene_annots$go.BP, gene_annots$go.MF))
selected_bp <- subset(go_df, id %in% selected_go_ids)
  


clean_mygene_output <- function(idx, category) {
  clean_df <- data.frame()
  required_fields <- c("id","term","evidence","pubmed")
  mygene_df <- as.data.frame(gene_annots$response[idx, paste("go.",category,sep="")][[1]])
  
  if (nrow(mygene_df) == 0) {
    mygene_df <- as.data.frame(setNames(replicate(4,NA, simplify = F), c("evidence","id","term","pubmed")))
  } else if (!all(required_fields %in% colnames(mygene_df))) {
    missing_fields <- setdiff(required_fields, colnames(mygene_df))
    mygene_df[missing_fields] <- NA
  }
  
  clean_df <- cbind('query_gene'=genes[[idx]], mygene_df[,required_fields], 'go_category'=category)
  
  return(clean_df)
}

