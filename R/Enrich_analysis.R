
library(clusterProfiler)
library(org.Hs.eg.db)

# Example: GO enrichment
enrich_result <- enrichGO(gene = subset_htr, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.05)

enrich_result_all <- enrichGO(gene = subset_htr_gp_tf, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.05)
# View results
head(enrich_result_all)

# Plot results
dotplot(enrich_result_all)
