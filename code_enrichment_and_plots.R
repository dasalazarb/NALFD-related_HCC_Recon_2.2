library(clusterProfiler)
library(AnnotationHub)
library(pathview)
library(ggpubr)
library(ggplot2)
library(ggrepel)
library(enrichplot)
library(dplyr)
library(tidyr)
library(DOSE)
library(hgnc) #https://github.com/maialab/hgnc
(url <- latest_archive_url())
hgnc_dataset <- import_hgnc_dataset(url); hgnc_dataset %>% names(); head(hgnc_dataset); dim(hgnc_dataset)

## frist part: read and clean results
results <- readxl::read_xlsx(path = "G:/.shortcut-targets-by-id/1xYVRN5UwMFXmj3OdNbFJmArfObc9kWSH/Tutorial_2018-10/07_ResearchArticles/Research_article_VIII/results/resultados_background_GSE10142_MoA.xlsx") %>% 
  mutate(on_off_rev = case_when(reversible == "IN HCC RIGHT" & REVERSIBLE_ON_HCC == "ON IN HCC" ~ "ON in right", 
                                reversible == "IN HCC LEFT" & REVERSIBLE_ON_HCC == "ON IN HCC" ~ "ON in left",
                                reversible == "IN HCC RIGHT" & REVERSIBLE_ON_HCC == "OFF IN HCC" ~ "OFF in right", 
                                reversible == "IN HCC LEFT" & REVERSIBLE_ON_HCC == "OFF IN HCC" ~ "OFF in left", 
                                reversible == "IN HCC LEFT" & is.na(REVERSIBLE_ON_HCC) ~ "HCC > NORMAL in left", 
                                reversible == "IN HCC RIGHT" & is.na(REVERSIBLE_ON_HCC) ~ "HCC > NORMAL in right", 
                                TRUE ~ "NA_NA"), 
         irreversible = case_when(irreversible == "ON" ~ "ON in NAFLD-related HCC", 
                                  irreversible == "OFF" ~ "OFF in NAFLD-related HCC", 
                                  TRUE ~ irreversible), 
         reversible = gsub("IN ", "", reversible), ); head(results); dim(results)

irrev_results <- unique(results$irreversible); irrev_results <- irrev_results[!is.na(irrev_results) & irrev_results != "FALSE"]; irrev_results <- irrev_results[!grepl(">|<", irrev_results)]
# rev_results <- unique(results$reversible); rev_results <- rev_results[!is.na(rev_results)];
on_off_rev_results <- unique(results$on_off_rev); on_off_rev_results <- on_off_rev_results[which(on_off_rev_results != "NA_NA")]; on_off_rev_results <- on_off_rev_results[!grepl(">|<", on_off_rev_results)]

giveMeGenes <- function(irrev_results, irreversible) {
  list_irrev_results <- list()
  for(i in irrev_results) {
    a <- results[results[,irreversible] == i, "genes"]
    a <- a[!is.na(a),]
    list_irrev_results[i] <- strsplit(gsub("//s+", " ", gsub("\\[||\\]||\\'||,", "", paste(a$genes, collapse = " "))), split = " ")
    list_irrev_results[[i]] <- unique(list_irrev_results[[i]])
    list_irrev_results[[i]] <- data.table::data.table(hgnc_id = sort(list_irrev_results[[i]][list_irrev_results[[i]] != ""]))
  }
  return(list_irrev_results)
}

list_irrev_results <- giveMeGenes(irrev_results, "irreversible"); list_irrev_results
# list_rev_results <- giveMeGenes(rev_results, "reversible"); list_rev_results
list_on_off_rev_results <- giveMeGenes(on_off_rev_results, "on_off_rev"); list_on_off_rev_results

## second part: hgnc to gene symbol
list_irrev_results <- lapply(list_irrev_results, function(x) x %>% 
         inner_join(hgnc_dataset, by="hgnc_id"))

# list_rev_results <- lapply(list_rev_results, function(x) x %>% 
#                                inner_join(hgnc_dataset, by="hgnc_id"))

list_on_off_rev_results <- lapply(list_on_off_rev_results, function(x) x %>% 
                             inner_join(hgnc_dataset, by="hgnc_id"))

## third part: enrichment analysis
list_irrev_results <- lapply(list_irrev_results, function(x) bitr(x$symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db"))
list_irrev_results <- lapply(list_irrev_results, function(x) x$ENTREZID)

# list_rev_results <- lapply(list_rev_results, function(x) bitr(x$symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db"))
# list_rev_results <- lapply(list_rev_results, function(x) x$ENTREZID)

list_on_off_rev_results <- lapply(list_on_off_rev_results, function(x) bitr(x$symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db"))
list_on_off_rev_results <- lapply(list_on_off_rev_results, function(x) x$ENTREZID)

## fourth part: enrichment analysis
ck1 <- compareCluster(geneCluster = list_irrev_results, fun = "enrichKEGG", pvalueCutoff=1e-10)
# ck2 <- as.data.frame(ck1)
ck1 <- dotplot(ck1)
# ck2 <- compareCluster(geneCluster = list_rev_results, fun = "enrichKEGG", pvalueCutoff=0.05)
# # ck2 <- as.data.frame(ck2)
# ck2 <- dotplot(ck2)
ck3 <- compareCluster(geneCluster = list_on_off_rev_results, fun = "enrichKEGG", pvalueCutoff=1e-10)
# ck2 <- as.data.frame(ck2)
ck3 <- dotplot(ck3)

jpeg(filename = paste0("G:/.shortcut-targets-by-id/1xYVRN5UwMFXmj3OdNbFJmArfObc9kWSH/Tutorial_2018-10/07_ResearchArticles/Research_article_VIII/results/Figure 5.jpg"), 
     units = 'in', res = 300, width = 16, height = 5)
ggarrange(ck1, ck3,
          labels = c("a", "b"), font.label = list(size=20),
          ncol = 2, nrow = 1)
dev.off()

pdf(file = paste0("G:/.shortcut-targets-by-id/1xYVRN5UwMFXmj3OdNbFJmArfObc9kWSH/Tutorial_2018-10/07_ResearchArticles/Research_article_VIII/results/Figure 5.pdf"), 
     width = 16, height = 5)
ggarrange(ck1, ck3,
          labels = c("a", "b"), font.label = list(size=20),
          ncol = 2, nrow = 1)
dev.off()


for (item in names(list_irrev_results)) {
  edo <- enrichDGN(list_irrev_results[[item]])
  edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
  plotFinal <- heatplot(edox, showCategory=5) + 
    labs(title = item, 
         y = "Enriched terms", x = "Gene symbols") + 
    theme(text = element_text(size = 20))
  ggsave(plotFinal, width = 15, height = 8,filename = paste0("G:/.shortcut-targets-by-id/1xYVRN5UwMFXmj3OdNbFJmArfObc9kWSH/Tutorial_2018-10/07_ResearchArticles/Research_article_VIII/results/plot irrev ", gsub(">", "_great_", gsub("<", "_less_",item)), ".pdf"))
}

for (item in names(list_on_off_rev_results)) {
  edo <- enrichDGN(list_on_off_rev_results[[item]])
  edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
  plotFinal <- heatplot(edox, showCategory=5) + 
    labs(title = item, 
         y = "Enriched terms", x = "Gene symbols") + 
    theme(text = element_text(size = 20))
  ggsave(plotFinal, width = 15, height = 8, 
         filename = paste0("G:/.shortcut-targets-by-id/1xYVRN5UwMFXmj3OdNbFJmArfObc9kWSH/Tutorial_2018-10/07_ResearchArticles/Research_article_VIII/results/plot rev ", gsub(">", "_great_", gsub("<", "_less_",item)), ".pdf"))
}

### fifth part: KEGG enrichment
results_ <- results %>% 
  dplyr::select(genes, SB_HCC, normal, irreversible, reversible, on_off_rev) %>% 
  mutate(genes = gsub("\\[||\\]||\\s+||\\'", "", genes)) %>% 
  tidyr::separate_rows(genes,sep = ",") %>% 
  dplyr::rename(hgnc_id = genes) %>% 
  inner_join(hgnc_dataset %>% dplyr::select(hgnc_id, symbol), by="hgnc_id")

giveMeEnzymes <- function(results_, irrev_results, irreversible, hcc_or_normal) {
  list_irrev_results <- list()
  for(i in irrev_results) {
    a <- results_[results_[,irreversible] == i, c(hcc_or_normal,"symbol")]
    a <- a[!is.na(a$symbol),]; c <- bitr(a$symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    c <- c %>% dplyr::rename(symbol = SYMBOL); a <- a %>% inner_join(c, by="symbol") #%>% filter(!is.na(symbol) | symbol != "")
    b <- dplyr::pull(a, hcc_or_normal); c <- a$ENTREZID
    names(b) <- c
    list_irrev_results[[i]] <- b
  }
  return(list_irrev_results)
}

normal <- function(x) {
  x[x > 10] <- 10
  x[x < -10] <- -10
  return(x)
  # if (is.na((x - min(x)) / (max(x) - min(x)))) {
  #   return(x)
  # } else {
  #   return((x - min(x)) / (max(x) - min(x)))
  # }
}

list_irrev_results_SB_HCC <- giveMeEnzymes(results_, irrev_results, "irreversible", "SB_HCC")
list_irrev_results_SB_HCC <- lapply(list_irrev_results_SB_HCC, function(x) normal(x))

# list_irrev_results_normal <- giveMeEnzymes(results_, irrev_results, "irreversible", "normal")
# list_irrev_results_normal <- lapply(list_irrev_results_normal, function(x) normal(x))

list_on_off_rev_results_SB_CC <- giveMeEnzymes(results_, on_off_rev_results, "on_off_rev", "SB_HCC")
list_on_off_rev_results_SB_CC <- lapply(list_on_off_rev_results_SB_CC, function(x) normal(x))

# list_on_off_rev_results_normal <- giveMeEnzymes(results_, on_off_rev_results, "on_off_rev", "normal")
# list_on_off_rev_results_normal <- lapply(list_on_off_rev_results_normal, function(x) normal(x))
a <- data.frame(matrix("", nrow=1,ncol = 2))
colnames(a) <- c("pathway", "ec")
b <- data.frame(matrix("", nrow=1,ncol = 2))
colnames(b) <- c("pathway", "ec")

x <- org.Hs.egENZYME
# Get the entrez gene identifiers that are mapped to an EC number 
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
xx <- xx[!grepl("\\-", xx)]

for (i in names(list_irrev_results_SB_HCC)) {
  setwd("G:/.shortcut-targets-by-id/1xYVRN5UwMFXmj3OdNbFJmArfObc9kWSH/Tutorial_2018-10/07_ResearchArticles/Research_article_VIII/results/irrev_results_SB_HCC/")
  kk2 <- gseKEGG(geneList     = sort(list_irrev_results_SB_HCC[[i]], decreasing = T),
                 organism     = 'hsa',
                 # minGSSize    = 120,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)
  kk2 <- kk2[kk2$ID != "hsa01100",] # sacar mapa de todo el metabolismo humano
  for (j in kk2$ID) {
    hsa04110 <- pathview(gene.data  = sort(list_irrev_results_SB_HCC[[i]], decreasing = T),
                         pathway.id = j,
                         species    = "hsa",
                         limit      = list(gene=max(abs(sort(list_irrev_results_SB_HCC[[i]], decreasing = T))), cpd=1),
                         kegg.dir = "G:/.shortcut-targets-by-id/1xYVRN5UwMFXmj3OdNbFJmArfObc9kWSH/Tutorial_2018-10/07_ResearchArticles/Research_article_VIII/results/irrev_results_SB_HCC/")
    
    if (!is.null(hsa04110[[1]])) {
      w <- hsa04110[[1]] %>% 
        filter(!is.na(mol.data)) %>% 
        pull(kegg.names)
      
      a <- a %>% 
        bind_rows(data.frame(pathway = kk2$Description[kk2$ID == j], 
                             ec = unlist(lapply(w, function(v) xx[grepl(v, names(xx))]))) %>% 
                    distinct())
    }
  }

}; setwd("C:/Users/da.salazarb/Documents")

for (i in names(list_on_off_rev_results_SB_CC)) {
  setwd("G:/.shortcut-targets-by-id/1xYVRN5UwMFXmj3OdNbFJmArfObc9kWSH/Tutorial_2018-10/07_ResearchArticles/Research_article_VIII/results/on_off_rev_results_SB_CC/")
  kk2 <- gseKEGG(geneList     = sort(list_on_off_rev_results_SB_CC[[i]], decreasing = T),
                 organism     = 'hsa',
                 # minGSSize    = 120,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)
  kk2 <- kk2[kk2$ID != "hsa01100",] # sacar mapa de todo el metabolismo humano
  for (j in kk2$ID) {
    hsa04110 <- pathview(gene.data  = sort(list_on_off_rev_results_SB_CC[[i]], decreasing = T),
                         pathway.id = j,
                         species    = "hsa",
                         limit      = list(gene=max(abs(sort(list_on_off_rev_results_SB_CC[[i]], decreasing = T))), cpd=1),
                         kegg.dir = "G:/.shortcut-targets-by-id/1xYVRN5UwMFXmj3OdNbFJmArfObc9kWSH/Tutorial_2018-10/07_ResearchArticles/Research_article_VIII/results/on_off_rev_results_SB_CC/")
    if (!is.null(hsa04110[[1]])) {
      w <- hsa04110[[1]] %>% 
        filter(!is.na(mol.data)) %>% 
        pull(kegg.names)
      
      b <- b %>% 
        bind_rows(data.frame(pathway = kk2$Description[kk2$ID == j], 
                             ec = unlist(lapply(w, function(v) xx[grepl(v, names(xx))]))) %>% 
                    distinct())
    }
  }
}; setwd("C:/Users/da.salazarb/Documents")

write.csv(a, "G:/.shortcut-targets-by-id/1xYVRN5UwMFXmj3OdNbFJmArfObc9kWSH/Tutorial_2018-10/07_ResearchArticles/Research_article_VIII/results/irrev_enzymes_in_kegg.csv")
write.csv(b, "G:/.shortcut-targets-by-id/1xYVRN5UwMFXmj3OdNbFJmArfObc9kWSH/Tutorial_2018-10/07_ResearchArticles/Research_article_VIII/results/rev_enzymes_in_kegg.csv")

