gene_list = read.table(file = "application/DLPFCdata/151672/151672_doamin2_gene_list.txt")
library(clusterProfiler)
library(org.Hs.eg.db)
gene_list = unlist(gene_list)
# Convert gene symbols to Entrez IDs
entrez_ids <- bitr(gene_list, fromType = "SYMBOL",
                   toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Perform GO enrichment analysis
go_results <- enrichGO(gene         = entrez_ids$ENTREZID,
                       OrgDb        = org.Hs.eg.db,
                       keyType      = "ENTREZID",
                       ont          = "BP",  # Change to "MF" or "CC" for other ontologies
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable     = TRUE)

BP_result = go_results@result
p4 = dotplot(go_results, title = "Gene set enrichment analysis for biological process", 
             showCategory = 20)

library(ggplot2)
ggsave(p4, filename = "reproduce/img/FigureS28.pdf", width = 8, height = 8, units = "in",  bg = "white")
