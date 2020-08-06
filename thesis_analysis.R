# Libraries
library(viridis)
library(dplyr)
library(purrr)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(purrr)
library(reshape2)
library(vegan)
library(cluster)
library(maSigPro)
library(gridExtra)
library(gplots)
library(ggalt)
set.seed(42)

#### Functions ####
createPairedData <- function(df_meta, pair_list){
  
  # Create groups e.g. stool vs dental
  df_meta_pairs <- data.frame()
  for(i in 1:length(pair_list)){
    
    samples_one <- unique(df_meta$Sample.name[df_meta$sample_type == pair_list[[i]][1]])
    samples_two <- unique(df_meta$Sample.name[df_meta$sample_type == pair_list[[i]][2]])
    df_meta_pair <- df_meta[(df_meta$Sample.name %in% Reduce(intersect,list(samples_one, samples_two))) & (df_meta$sample_type %in% pair_list[[i]]),]
    
    df_meta_pair$group <- paste(pair_list[[i]][1], "vs.", pair_list[[i]][2])
    if(nrow(df_meta_pairs) == 0) {
      df_meta_pairs <- df_meta_pair
    } else {
      df_meta_pairs <- rbind(df_meta_pairs, df_meta_pair)
    }
  }
  
  # Order by Location and change characters in group
  df_meta_pairs <- df_meta_pairs[order(df_meta_pairs$Location),]
  df_meta_pairs$group <- as.character(df_meta_pairs$group)
  
  return(df_meta_pairs)
}

runTtest <- function(df_paired){
  
  df_paired_richness <- df_paired[!duplicated(paste0(df_paired$ID, df_paired$group)),]
  
  # T-test
  p_values <- c()
  ttest_groups <- df_paired_richness[!duplicated(paste0(df_paired_richness$Location,
                                                        df_paired_richness$group)),]
  for(i in 1:nrow(ttest_groups)){
    y <- df_paired_richness[df_paired_richness$Location == ttest_groups$Location[i] & 
                              df_paired_richness$group == ttest_groups$group[i],]
    y <- y[order(y$Sample.name),]
    type1 <- strsplit(ttest_groups$group[i], " vs. ")[[1]][1]
    type2 <- strsplit(ttest_groups$group[i], " vs. ")[[1]][2]
    y1 <- y$richness[y$sample_type == type1]
    y2 <- y$richness[y$sample_type == type2]
    p_values <- c(p_values, wilcox.test(y1, y2, paired = TRUE, alternative = "two.sided")$p.value)
  }
  ttest_groups$pvalue <- p_values
  asterisk <- rep(NA, length(p_values))
  asterisk[p_values < 0.05] <- "*"
  asterisk[p_values < 0.01] <- "**"
  asterisk[p_values < 0.001] <- "***"
  ttest_groups$asterisk <- asterisk
  
  # Order ttest_groups
  ttest_groups <- ttest_groups[order(ttest_groups$Location, ttest_groups$group),]
  ttest_groups <- ttest_groups[,c("Location", "group", "pvalue", "asterisk")]
  return(ttest_groups)
}

plotRichnessGraph <- function(df_paired_richness_group, ttest_group, cols) {
  set.seed(1) # for jitter
  g <- ggplot(df_paired_richness_group, aes(sample_type, richness)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(size = 0.8) +
    theme_classic() +
    ylab("Plasmid Richness") +
    xlab("") +
    ggtitle(ttest_group$Location) +
    geom_text(data = ttest_group, aes(label=asterisk),
              x = 1.5, y = max(df_paired_richness_group$richness)+10, size = 7,
              inherit.aes = FALSE) +
    theme(legend.text = element_text(size = 14), legend.title = element_text(size = 16))
  return(g)
}

plotMultipleRichnessGraphs <- function(ttest_groups, df_paired, cols){
  df_paired_richness <- df_paired[!duplicated(paste0(df_paired$ID, df_paired$group)),]
  g <- list()
  for(i in 1:nrow(ttest_groups)){
    df_paired_richness_group <- df_paired_richness[df_paired_richness$Location == ttest_groups$Location[i] & 
                                                     df_paired_richness$group == ttest_groups$group[i],]
    g[[i]] <- plotRichnessGraph(df_paired_richness_group, ttest_groups[i,], cols) +
      theme(legend.position = "none") + ylim(c(0, 50))
  }
  return(g)
}

getBestBlastHit <- function(blast_result, qseqid_id){
  blast_result_tmp <- blast_result[blast_result$qseqid == qseqid_id,]
  remove_ids <- c()
  for(i in 1:(nrow(blast_result_tmp)-1)){
    
    for(j in (i+1):nrow(blast_result_tmp)){
      
      if(blast_result_tmp$qend[i] >= blast_result_tmp$qstart[j] & blast_result_tmp$qend[j] >= blast_result_tmp$qstart[i]){
        overlap_prop <- (blast_result_tmp$qend[i] - blast_result_tmp$qstart[j]) / min(blast_result_tmp$length[i], blast_result_tmp$length[j])
      } else if(blast_result_tmp$qend[j] >= blast_result_tmp$qstart[i] & blast_result_tmp$qstart[j] <= blast_result_tmp$qend[i]){
        overlap_prop <- (blast_result_tmp$qend[j] - blast_result_tmp$qstart[i]) / min(blast_result_tmp$length[i], blast_result_tmp$length[j])
      } else if(blast_result_tmp$qend[i] >= blast_result_tmp$qend[j] & blast_result_tmp$qstart[i] <= blast_result_tmp$qstart[j]){
        overlap_prop <- (blast_result_tmp$qend[i] - blast_result_tmp$qstart[i]) / min(blast_result_tmp$length[i], blast_result_tmp$length[j])
      } else if(blast_result_tmp$qend[j] >= blast_result_tmp$qend[i] & blast_result_tmp$qstart[j] <= blast_result_tmp$qstart[i]){
        overlap_prop <- (blast_result_tmp$qend[j] - blast_result_tmp$qstart[j]) / min(blast_result_tmp$length[i], blast_result_tmp$length[j])
      } else {
        overlap_prop <- 0
      }
      
      if(overlap_prop > 0.2 & blast_result_tmp[i, "evalue"] != blast_result_tmp[j, "evalue"]){
        ind <- which.max(blast_result_tmp$evalue[c(i,j)])
        remove_id <- which(blast_result$qseqid == blast_result_tmp$qseqid[c(i,j)[ind]] & 
                             blast_result$sseqid == blast_result_tmp$sseqid[c(i,j)[ind]] &
                             blast_result$qstart == blast_result_tmp$qstart[c(i,j)[ind]] & 
                             blast_result$qend == blast_result_tmp$qend[c(i,j)[ind]])
        remove_ids <- unique(c(remove_ids, remove_id))
        
      } else if (overlap_prop > 0.2 & blast_result_tmp[i, "evalue"] == blast_result_tmp[j, "evalue"] & blast_result_tmp[i, "pident"] != blast_result_tmp[j, "pident"]){
        ind <- which.min(blast_result_tmp$pident[c(i,j)])
        remove_id <- which(blast_result$qseqid == blast_result_tmp$qseqid[c(i,j)[ind]] & 
                             blast_result$sseqid == blast_result_tmp$sseqid[c(i,j)[ind]] &
                             blast_result$qstart == blast_result_tmp$qstart[c(i,j)[ind]] & 
                             blast_result$qend == blast_result_tmp$qend[c(i,j)[ind]])
        remove_ids <- unique(c(remove_ids, remove_id))
      }
    }
  }
  return(remove_ids)
}


#### Read data ####
metadata <- read.csv("data/supplementary_metadata.csv", stringsAsFactors = FALSE)
samples <- read.delim("data/samples_all.txt", stringsAsFactors = FALSE, header = FALSE)
metadata <- metadata[metadata$ID %in% samples$V1 & metadata$Location %in% c("China", "USA"),]
counts_total <- readRDS("data/counts_total.RDS")
counts_total <- counts_total[,names(counts_total) %in% metadata$ID]
metaphlan <- read.csv("data/all_metaphlan.csv", stringsAsFactors = FALSE)
metadata$Location_sampletype <- paste(metadata$sample_type, "-", metadata$Location)
rownames(metadata) <- metadata$ID

# Summarise
metadata_catalog_summary <- metadata %>% group_by(Location, sample_type, Sample.name) %>%
  mutate(timepoint = rank(as.numeric(Visit_Number))) %>%
  group_by(Location, sample_type, timepoint) %>%
  summarise(n = n_distinct(ID))

# GIT site colours
cols <- plasma(length(unique(metadata$sample_type)), end = 0.8)
names(cols) <- sort(unique(metadata$sample_type))

# Cohort colours
cohort_cols <- c("grey", brewer.pal(9, "Blues")[c(5,7)], "gold", brewer.pal(9, "YlOrRd")[c(5,7)], brewer.pal(9, "RdPu")[c(3,5)])
names(cohort_cols) <- sort(unique(metadata$Location_sampletype)) 

# Combine metaplasmidspades
plasmid_res <- data.frame()
for (i in 1:length(metadata$ID)) {
  scaffold_res <- read.csv(paste0("data/scaffold_results/", metadata$ID[i], "_scaffolds_results_table.csv"), stringsAsFactors = FALSE, header = FALSE)
  scaffold_res$ID <- metadata$ID[i]
  plasmid_res <- rbind(plasmid_res, scaffold_res)
}
plasmid_res <- plasmid_res[plasmid_res$V2 == "Plasmid",]
plasmid_res$size <- as.numeric(sapply(plasmid_res$V1, function(x) strsplit(x, "_")[[1]][4]))
plasmid_res <- plasmid_res[order(plasmid_res$size, decreasing = FALSE),]

# # Network of plasmids
# count = 0
# for (i in 1:(nrow(plasmid_res)-1)) {
#   for (j in (i+1):nrow(plasmid_res)) {
#     count = count + 1
#   }
# }
# 
# plas_net <- data.frame(matrix(NA, nrow = count, ncol = 3))
# names(plas_net) <- c("plasmid1", "plasmid2", "weight")
# count = 0
# for (i in 1:(nrow(plasmid_res)-1)) {
#   plas1 <- plasmid_res[i,]
#   for (j in (i+1):nrow(plasmid_res)) {
#     plas2 <- plasmid_res[j,]
#     size_prop <- plas1$size / plas2$size
#     if (size_prop > 0.99) {
#       plas1_prot <- unique(strsplit(plas1$V4, " ")[[1]])
#       plas2_prot <- unique(strsplit(plas2$V4, " ")[[1]])
#       shared <- length(intersect(plas1_prot, plas2_prot)) / length(unique(c(plas1_prot, plas2_prot)))
#       count = count + 1
#       plas_net[count,] <- c(plas1$V1, plas2$V1, weight = size_prop * shared)
#     } else {
#       break
#     }
#   }
# }
# plas_net <- plas_net[!is.na(plas_net$weight),]
# saveRDS(plas_net, "data/plasmid_network.RDS")

plas_net <- readRDS("data/plasmid_network.RDS")

# # Cluster plasmids where weight > 0.8
# plas_cluster_list <- list()
# cluster_num <- 0
# for (i in 1:nrow(plas_net)) {
#   if (plas_net$weight[i] > 0.8) {
#     plas1 <- plas_net$plasmid1[i]
#     plas2 <- plas_net$plasmid2[i]
#     if (any(grepl(plas1, plas_cluster_list))) {
#       if (any(grepl(plas2, plas_cluster_list))) {
#         # Do nothing
#       } else {
#         plas_cluster_list[[grep(plas1, plas_cluster_list)]] <- c(plas_cluster_list[[grep(plas1, plas_cluster_list)]], plas2)
#       }
#     } else {
#       if (any(grepl(plas2, plas_cluster_list))) {
#         plas_cluster_list[[grep(plas2, plas_cluster_list)]] <- c(plas_cluster_list[[grep(plas2, plas_cluster_list)]], plas1)
#       } else {
#         plas_cluster_list[[paste0("PC_", as.character(cluster_num))]] <- c(plas1, plas2)
#         cluster_num <- cluster_num + 1
#       }
#     }
#   }
# }
# saveRDS(plas_cluster_list, "data/plasmid_cluster_list.RDS")

# Include plasmid clusters
plas_cluster_list <- readRDS("data/plasmid_cluster_list.RDS")
plasmid_cluster <- rep(NA, nrow(plasmid_res))
for (i in 1:length(plas_cluster_list)) {
  plasmid_cluster[plasmid_res$V1 %in% plas_cluster_list[[i]]] <- names(plas_cluster_list)[i]
}
plasmid_res$plasmid_cluster <- plasmid_cluster
plasmid_res$plasmid_name_cluster <- plasmid_res$plasmid_cluster
plasmid_res$plasmid_name_cluster[is.na(plasmid_res$plasmid_name_cluster)] <- plasmid_res$V1[is.na(plasmid_res$plasmid_name_cluster)]

plasmid_res <- plasmid_res %>%
  select(-c(V2,ID)) %>%
  group_by(plasmid_name_cluster) %>%
  mutate(av_size = mean(size))

# Plasmid breadth coverage
plasmid_bed <- data.frame()
for (i in 1:length(metadata$ID)) {
  if (file.info(paste0("data/coverage/", metadata$ID[i], ".bam.bedtools.coverage.txt"))$size > 0) {
    plasmid_cov <- read.delim(paste0("data/coverage/", metadata$ID[i], ".bam.bedtools.coverage.txt"), stringsAsFactors = FALSE, header = FALSE)
    plasmid_cov$ID <- metadata$ID[i]
    plasmid_bed <- rbind(plasmid_bed, plasmid_cov)
  }
}

# Filter by 75% breadth coverage 
plasmid_bed <- plasmid_bed[plasmid_bed$V7 >= 0.75,]

# Combine info
plasmid_prof <- plasmid_bed %>%
  select(V1, V4, V7, ID) %>%
  rename("cov_depth"="V4", "cov_breadth"="V7") %>%
  inner_join(plasmid_res, by = "V1") %>%
  rename("perc_plasmid_prot" = "V3", "plasmid_prot" = "V4")

# Add taxa
plasmid_taxa <- data.frame()
taxa_files <- list.files('data/plasmidfinder', pattern = "*.out", full.names = TRUE)
for (i in 1:length(taxa_files)) {
  if (file.info(taxa_files[i])$size > 0) {
    plasmid_rep <- read.delim(taxa_files[i], stringsAsFactors = FALSE, header = FALSE)
    plasmid_rep$type <- gsub(".out", "", gsub(".*plasmid_nonredundant_catalog_", "", taxa_files[i]))
    plasmid_taxa <- rbind(plasmid_taxa, plasmid_rep)
  }
}

plasmid_taxa <- plasmid_taxa %>%
  select(V1, V2, type) %>%
  rename("genome"="V2")

metadata <- metadata[metadata$ID %in% plasmid_prof$ID,]

# Summarise new metadata
metadata_summary <- metadata %>% group_by(Location, sample_type, Sample.name) %>%
  mutate(timepoint = rank(as.numeric(Visit_Number))) %>%
  group_by(Location, sample_type, timepoint) %>%
  summarise(n = n_distinct(ID))
metadata <- metadata %>% group_by(Location, sample_type, Sample.name) %>%
  mutate(timepoint = rank(as.numeric(Visit_Number))) %>%
  data.frame()
rownames(metadata) <- metadata$ID

# Plasmid profile
plasmid_prof <- left_join(plasmid_prof, plasmid_taxa, by = "V1") %>%
  rename("name"="V1") %>%
  left_join(metadata, by = "ID")

# Number of unique plasmid sequences
length(unique(plasmid_prof$name))

# Number of oral and stool samples
length(unique(metadata$ID[metadata$sample_type == "stool"]))
length(unique(metadata$ID[metadata$sample_type != "stool"]))

# Number of singletons
length(unique(plasmid_prof$name[is.na(plasmid_prof$plasmid_cluster)]))

# Number of plasmid sequences in clusters
length(unique(plasmid_prof$name)) - length(unique(plasmid_prof$name[is.na(plasmid_prof$plasmid_cluster)]))

# Number of clusters
length(unique(plasmid_prof$plasmid_name_cluster[!is.na(plasmid_prof$plasmid_cluster)]))

# Without USA longidutinal samples
plasmid_prof_sing <- plasmid_prof %>% 
  filter(timepoint == 1)

# Plasmid profile summary
plasmid_prof_summary <- plasmid_prof_sing %>%
  group_by(Location_sampletype, Location, sample_type) %>%
  mutate(total = n_distinct(ID)) %>%
  group_by(Location_sampletype, Location, sample_type, plasmid_name_cluster, total) %>%
  summarise(n = n_distinct(ID)) %>%
  mutate(perc = n/total*100)

# Graph of incidence frequency
tiff("figures/plasmid_count_summary.tiff", height = 500, width = 2000, res = 150)
ggplot(plasmid_prof_summary, aes(perc, fill = Location_sampletype)) +
  geom_histogram(binwidth = 10) +
  facet_grid(~Location_sampletype) +
  theme_classic() +
  scale_fill_manual("GIT site - Country", values = cohort_cols) +
  xlab("% samples") + ylab("Frequency")
dev.off()

#### BETA-DIVERSITY ####
plasmid_prof_top <- plasmid_prof_sing[plasmid_prof_sing$plasmid_name_cluster %in% plasmid_prof_summary$plasmid_name_cluster[plasmid_prof_summary$perc > 50],]
plasmid_prof_top$labels <- plasmid_prof_top$plasmid_name_cluster
plasmid_prof_top$labels[!is.na(plasmid_prof_top$type)] <- paste0(plasmid_prof_top$plasmid_name_cluster[!is.na(plasmid_prof_top$type)], " (", plasmid_prof_top$type[!is.na(plasmid_prof_top$type)], ")")
plasmid_read_counts <- dcast(plasmid_prof_top, labels ~ ID, value.var = "cov_depth", fun.aggregate = sum)
row.names(plasmid_read_counts) <- plasmid_read_counts$labels
plasmid_read_counts <- plasmid_read_counts[,-1]
plasmid_counts <- plasmid_read_counts
plasmid_counts[plasmid_counts > 0] <- 1

# Heatmap
tiff("figures/plasmid_incidence_heatmap.tiff", width = 1500, height = 2000, res = 150)
heatmap.2(as.matrix(plasmid_counts),
          margins = c(10,20),
          trace = "none",
          scale = "none",
          hclustfun = function(x) {hclust(x, method = "ward.D2")},
          dendrogram = "column",
          col =  colorRampPalette(brewer.pal(9, "PuRd")[c(1, 9)])(50),
          breaks = seq(min(plasmid_counts), max(plasmid_counts), length.out = 51),
          symbreaks = FALSE,
          key = FALSE, 
          lhei = c(1,8),
          ColSideColors = sapply(metadata[colnames(plasmid_counts), "Location_sampletype"], function(x) cohort_cols[x]),
          #labRow = as.expression(lapply(family_labels, function(a) bquote(italic(.(a))))),
          labCol = NA,
          cexCol = 0.5,
          ylab = "Plasmid singleton/cluster",
          main = NA
)
legend(x = 0.85, y = 1.05, xpd=TRUE, legend = levels(as.factor(metadata[colnames(plasmid_counts), "Location_sampletype"])),
       col = cohort_cols[names(cohort_cols) %in% sort(unique(metadata[colnames(plasmid_counts), "Location_sampletype"]))], bg = "white", box.col = "black",
       lty = 1, lwd = 5, cex = 0.5, title = "GIT Sites")
dev.off()

# Compute distance matrix and run NMDS on plasmid clusters
plasmid_read_prop <- dcast(plasmid_prof_sing, plasmid_name_cluster ~ ID, value.var = "cov_depth", fun.aggregate = length)
rownames(plasmid_read_prop) <- plasmid_read_prop$plasmid_name_cluster
plasmid_read_prop <- plasmid_read_prop[,-1]
cluster_samples_nmds <- metaMDS(t(plasmid_read_prop), distance = "bray", k = 2, trymax = 20)
df_cluster_samples_nmds <- as.data.frame(cluster_samples_nmds$points)

# Cluster raw data for check
cluster_samples_raw <- vegdist(t(plasmid_read_prop))

# Silhoette analysis of PAM (k-medoids)
avg_sil <- numeric(20)
for(k in 3:(length(avg_sil)+1)) {
  tmp <- silhouette(pam(df_cluster_samples_nmds[,c("MDS1", "MDS2")], k = k), df_cluster_samples_nmds[,c("MDS1", "MDS2")])
  avg_sil[k-1] <- mean(tmp[,3])
}

# Silhoette analysis of raw data for check
avg_sil_check <- numeric(20)
for(k in 3:(length(avg_sil_check)+1)) {
  tmp <- silhouette(pam(as.matrix(cluster_samples_raw), k = k), as.matrix(cluster_samples_raw))
  avg_sil_check[k-1] <- mean(tmp[,3])
}

# Group by silhouette width
samples_clust <- pam(df_cluster_samples_nmds[,c("MDS1", "MDS2")], which.max(avg_sil)+1)
samples_clust_check <- pam(as.matrix(cluster_samples_raw), which.max(avg_sil_check)+1)
df_cluster_samples_nmds$cluster <- as.factor(samples_clust$cluster[row.names(df_cluster_samples_nmds)])
df_cluster_samples_nmds$cluster_raw <- as.factor(samples_clust_check$cluster[row.names(df_cluster_samples_nmds)])
df_cluster_samples_nmds$ID <- row.names(df_cluster_samples_nmds)
df_cluster_samples_nmds$Sample.name <- as.character(sapply(df_cluster_samples_nmds$ID, function(x) metadata$Sample.name[metadata$ID == x]))
df_cluster_samples_nmds$Location <- sapply(df_cluster_samples_nmds$ID, function(x) metadata$Location[metadata$ID == x])
df_cluster_samples_nmds$sample_type <- sapply(df_cluster_samples_nmds$ID, function(x) metadata$sample_type[metadata$ID == x])
df_cluster_samples_nmds$Location_sampletype <- paste(df_cluster_samples_nmds$sample_type, "-", df_cluster_samples_nmds$Location)
df_cluster_samples_nmds$Age <- sapply(df_cluster_samples_nmds$ID, function(x) metadata$Age[metadata$ID == x])
df_cluster_samples_nmds$Sex <- sapply(df_cluster_samples_nmds$ID, function(x) metadata$Sex[metadata$ID == x])

# Plot NMDS by sample type and cluster
tiff("figures/nmds_clusters.tiff", width = 2000, height = 1000, res = 300)
ggplot(df_cluster_samples_nmds, aes(MDS1, MDS2, colour = Location_sampletype, shape = cluster)) +
  theme_classic() +
  geom_text(aes(label=cluster)) +
  scale_colour_manual("GIT Site - Country", values = cohort_cols, guide = guide_legend(override.aes = list(shape = 21, size = 4))) +
  xlab("NMDS 1") + ylab("NMDS 2")
dev.off()

# Check clusters with raw data
ggplot(df_cluster_samples_nmds, aes(MDS1, MDS2, colour = Location_sampletype, shape = cluster_raw)) +
  theme_classic() +
  geom_text(aes(label=cluster_raw)) +
  scale_colour_manual("Body Site", values = cohort_cols, guide = guide_legend(override.aes = list(shape = 21, size = 3))) +
  xlab("NMDS 1") + ylab("NMDS 2")

# Calculate proportion of samples
cluster_res <- df_cluster_samples_nmds %>% 
  group_by(cluster, Location, sample_type) %>% 
  summarise(n = n()) %>%
  group_by(cluster) %>%
  mutate(total_n = sum(n)) %>%
  mutate(prop_cluster = n/total_n*100, Location_sampletype = paste(sample_type, "-", Location))

# Calculate proportion of body sites
cluster_site <- cluster_res %>% 
  group_by(cluster, sample_type) %>%
  summarise(prop_cluster = signif(sum(prop_cluster), 3)) %>%
  mutate(summary = paste(sample_type, ": ", prop_cluster, "%")) %>%
  group_by(cluster) %>%
  summarise(summary = paste(summary, collapse = "; ")) %>%
  mutate(summary = paste0("Group ", cluster, " (", summary, ")"))

# Breakdown of groups
tiff("figures/nmds_group_breakdown.tiff", width = 1200, height = 1000, res = 300)
ggplot(cluster_res, aes(cluster, prop_cluster, fill = Location_sampletype)) +
  geom_bar(stat = "identity") +
  theme_classic() + xlab("Group") + ylab("% Samples") +
  scale_fill_manual("GIT Site - Country", values = cohort_cols)
dev.off()

# Venn diagram
unique_clusters <- unique(df_cluster_samples_nmds$cluster)

plas_group_list <- sapply(unique_clusters, function(x) {
  plas_tmp <- plasmid_read_prop[,colnames(plasmid_read_prop) %in% df_cluster_samples_nmds$ID[df_cluster_samples_nmds$cluster == x]]
  plas_group <- rownames(plas_tmp)[rowSums(plas_tmp) != 0]
  return(plas_group)
})

group_names <- paste("Group", unique_clusters, c("\n(mostly\ndental and\nsaliva)", "\n(mostly buccal mucosa)", "\n(stool)", "\n(mostly\ndorsum of tongue)"))
names(plas_group_list) <- group_names
max_group_length <- max(sapply(plas_group_list, function(x) length(x)))

plas_group_list <- lapply(plas_group_list, function(x) c(x, rep(NA, max_group_length - length(x))))

plas_group_all <- do.call(cbind, plas_group_list)

tiff("figures/venn_diagram.tiff", width = 600, height = 600)
suma2Venn(plas_group_all, cexil = 2, cexsn = 1.4)
dev.off()

#### Alpha-diversity ####
# Remove samples with lower than three plasmid clusters
plasmid_counts <- dcast(plasmid_prof_sing, ID ~ plasmid_name_cluster, value.var = "name", fun.aggregate = length)
row.names(plasmid_counts) <- plasmid_counts[,1]
plasmid_counts <- plasmid_counts[,-1]
remove_ids <- rownames(plasmid_counts)[rowSums(plasmid_counts > 0) <= 3]
plasmid_counts <- plasmid_counts[!rownames(plasmid_counts) %in% remove_ids,]
metadata_richness <- metadata[!metadata$ID %in% unique(c(remove_ids, metadata$ID[!metadata$ID %in% row.names(plasmid_counts)])),]

pair_list <- list(c("stool", "dental"), c("stool", "saliva"), c("dental", "saliva"),
                  c("stool", "dorsum of tongue"), c("stool", "buccal mucosa"),
                  c("dorsum of tongue", "buccal mucosa"), c("dorsum of tongue", "dental"), c("buccal mucosa", "dental"))
paired_metadata <- createPairedData(metadata[metadata$ID %in% row.names(plasmid_counts),], pair_list)

# Get richness from matrix
richness_paired <- data.frame(ID = row.names(plasmid_counts), richness = rowSums(plasmid_counts > 0), no_plasmids = rowSums(plasmid_counts)) %>%
  right_join(paired_metadata, by = "ID")

# Alpha-diversity vs. number of contigs
richness <- data.frame(ID = row.names(plasmid_counts), richness = rowSums(plasmid_counts > 0), no_plasmids = rowSums(plasmid_counts)) %>%
  right_join(metadata[metadata$ID %in% unique(paired_metadata$ID),], by = "ID")

linear_mod <- lm(richness ~ no_plasmids, richness)
summary(linear_mod)
richness$predict <- predict(linear_mod, richness)
tiff("figures/richness_no_plasmids.tiff", width = 1500, height = 750, res = 150)
ggplot(richness) +
  geom_point(aes(no_plasmids, richness, color = Location_sampletype)) +
  xlab("No. plasmids") + ylab("Plasmid Cluster/Singleton Richness") +
  theme_classic() +
  scale_color_manual("GIT Site - Country", values = cohort_cols) +
  geom_abline(color = "red", slope = linear_mod$coefficients[2], intercept = linear_mod$coefficients[1])
dev.off()

# For each group, remove samples with less than or equal to 100 plasmid contigs
unique_groups <- unique(richness_paired$group)
for (i in 1:length(unique_groups)) {
  remove_samples <- richness_paired$Sample.name[(richness_paired$no_plasmids <= 20 & richness_paired$group %in% unique_groups[i])]
  richness_paired <- richness_paired[!(richness_paired$Sample.name %in% remove_samples & richness_paired$group %in% unique_groups[i]),]
}
paired_metadata_summary <- paired_metadata %>% filter(ID %in% richness_paired$ID) %>%
  group_by(Location, group, sample_type) %>% summarise(n())

# Subsample matrix and calculate richness
richness_paired_ss <- data.frame()
unique_groups <- unique(richness_paired$group)
for (i in 1:length(unique_groups)) {

  group_ids <- unique(richness_paired$ID[richness_paired$group %in% unique_groups[i]])
  plasmid_counts_tmp <- plasmid_counts[rownames(plasmid_counts) %in% group_ids,]
  min_plasmids <- min(rowSums(plasmid_counts_tmp))

  plasmid_counts_tmp <- t(apply(plasmid_counts_tmp, 1, function(x) {
    while (sum(x) > min_plasmids) {
      ss_index <- sample(1:length(x), 1, prob = ifelse(x > 0, x/sum(x), 0))
      x[ss_index] <- x[ss_index] - 1
    }
    return(x)
  }))

  richness_paired_tmp <- data.frame(ID = rownames(plasmid_counts_tmp), richness = rowSums(plasmid_counts_tmp > 0)) %>%
    left_join(metadata_richness, by = "ID") %>%
    mutate(group = unique_groups[i])

  richness_paired_ss <- rbind(richness_paired_ss, richness_paired_tmp)
}

# T-test and graphs of subsampled data
richness_ttest_ss <- runTtest(richness_paired_ss)
richness_graphs_ss <- plotMultipleRichnessGraphs(richness_ttest_ss, richness_paired_ss, cols)
#richness_graphs_ss[[length(richness_graphs_ss)+1]] <- g_legend(plotRichnessGraph(richness_paired_ss, richness_ttest_ss, cols))

# Plot graph
lay <- rbind(c(1,2,3), c(4,5,6), c(7,8,9))
tiff("figures/alpha_diversity_subsampled.tiff", width = 2000, height = 2500, res = 250)
grid.arrange(grobs = richness_graphs_ss, layout_matrix = lay)
dev.off()

#### LONGITUDINAL SAMPLES ####
# Get samples with 3 or more timepoints
rm_samples <- c("SRS019025")
metadata_longus <- metadata[metadata$Location == "USA",] %>% 
  group_by(Sample.name, sample_type) %>%
  filter(!ID %in% rm_samples) %>%
  filter(any(timepoint == 3))

metadata_longus_summary <- metadata_longus %>%
  group_by(sample_type) %>%
  summarise(n = n_distinct(Sample.name))

plasmid_prof_long <- left_join(metadata_longus[,c("ID", "Sample.name", "timepoint")], plasmid_prof)

# Count average time points of plasmids
cluster_no_tp <- plasmid_prof_long %>% group_by(sample_type, Sample.name, plasmid_name_cluster) %>%
  summarise(no_tp = n_distinct(timepoint)) 

# Count all distinct plasmid clusters and singletons
cluster_all_summary <- plasmid_prof_long %>% group_by(sample_type, Sample.name) %>%
  summarise(total_plasmids = n_distinct(plasmid_name_cluster))

# Count clusters in at least x number of time points
cluster_count_tp <- map_df(.x = unique(cluster_no_tp$no_tp), .f = function(.x) {
  tmp <- data.frame(sample_type = rep(cluster_no_tp$sample_type[cluster_no_tp$no_tp == .x], .x),
                    Sample.name = rep(cluster_no_tp$Sample.name[cluster_no_tp$no_tp == .x], .x),
                    plasmid_name_cluster = rep(cluster_no_tp$plasmid_name_cluster[cluster_no_tp$no_tp == .x]),
                    no_tp = rep(cluster_no_tp$no_tp[cluster_no_tp$no_tp == .x], .x))
  tmp$tp <- c(1:.x)
  return(tmp)
})

cluster_count_tp_summary <- cluster_count_tp %>% group_by(sample_type, Sample.name, tp) %>%
  summarise(total_cluster_least = n_distinct(plasmid_name_cluster))

cluster_one <- cluster_count_tp_summary %>%
  filter(tp == 1) %>%
  rename(total_cluster = total_cluster_least) %>%
  select(-tp)

cluster_count_tp_summary <- left_join(cluster_count_tp_summary, cluster_one, by = c("sample_type", "Sample.name")) %>%
  mutate(cluster_frac = total_cluster_least/total_cluster)

# Plot
tiff("figures/longitudinal_cluster_count.tiff", height = 600, width = 1500, res = 100)
set.seed(1)
ggplot(cluster_count_tp_summary, aes(as.factor(tp), cluster_frac, fill = as.factor(sample_type))) +
  geom_boxplot() +
  geom_jitter() +
  facet_grid(~ sample_type, scale = "free", space = "free") +
  theme_classic() + xlab("No. timepoints") + ylab("Proportion of plasmid clusters/singletons") +
  theme(strip.background = element_blank(), strip.text.x = element_blank(), text = element_text(size = 16)) +
  scale_fill_manual("GIT Site", values = cols)
dev.off()

# Persistent and transient plasmid clusters and singletons
persist_graphs <- list()
transient_graphs <- list()
n_samples <- rep(NA, length(unique(plasmid_prof_long$sample_type)))
n_persist_clusters <- rep(NA, length(unique(plasmid_prof_long$sample_type)))
n_transient_clusters <- rep(NA, length(unique(plasmid_prof_long$sample_type)))
permanova_persist <- list()
permanova_transient <- list()
cluster_persist_all <- data.frame()
cluster_transient_all <- data.frame()
cluster_sharing_stats <- data.frame()

for (i in 1:length(unique(plasmid_prof_long$sample_type))) {
  cluster_persist <- left_join(plasmid_prof_long, cluster_no_tp, by = c("sample_type", "Sample.name", "plasmid_name_cluster")) %>%
    filter(no_tp >= 3 & sample_type == unique(plasmid_prof_long$sample_type)[i]) %>%
    group_by(ID) %>%
    filter(n_distinct(plasmid_name_cluster) > 1)
  
  # Select transient clusters
  cluster_transient <- left_join(plasmid_prof_long, cluster_no_tp, by = c("sample_type", "Sample.name", "plasmid_name_cluster")) %>%
    filter(no_tp < 3 & sample_type == unique(plasmid_prof_long$sample_type)[i]) %>%
    group_by(ID) %>%
    filter(n_distinct(plasmid_name_cluster) > 1)
  
  # Same individuals in transient and persistent
  cluster_persist <- cluster_persist %>%
    filter(Sample.name %in% intersect(unique(cluster_transient$Sample.name), unique(cluster_persist$Sample.name)))
  
  # Wilcoxon rank test of shared clusters
  sample_names <- unique(cluster_persist$Sample.name)
  cluster_persist_sharing <- rep(NA, length(sample_names))
  cluster_transient_sharing <- rep(NA, length(sample_names))
  for (j in 1:length(sample_names)) {
    cluster_persist_sample <- unique(cluster_persist$plasmid_name_cluster[cluster_persist$Sample.name == sample_names[j]])
    cluster_persist_others <- unique(cluster_persist$plasmid_name_cluster[cluster_persist$Sample.name != sample_names[j]])
    cluster_persist_sharing[j] <- sum(cluster_persist_sample %in% cluster_persist_others)/length(cluster_persist_sample)*100 
    
    cluster_transient_sample <- unique(cluster_transient$plasmid_name_cluster[cluster_transient$Sample.name == sample_names[j]])
    cluster_transient_others <- unique(cluster_transient$plasmid_name_cluster[cluster_transient$Sample.name != sample_names[j]])
    cluster_transient_sharing[j] <- sum(cluster_transient_sample %in% cluster_transient_others)/length(cluster_transient_sample)*100 
  }
  cluster_sharing_stats <- rbind(cluster_sharing_stats, data.frame(persist_med = median(cluster_persist_sharing), transient_med = median(cluster_transient_sharing),
                                                                   persist_iqr = IQR(cluster_persist_sharing), transient_iqr = IQR(cluster_transient_sharing),
                                                                   pval = wilcox.test(cluster_persist_sharing, cluster_transient_sharing)$p.value, sample_type = unique(plasmid_prof_long$sample_type)[i]))
  
  cluster_persist_all <- rbind(cluster_persist_all, as.data.frame(cluster_persist))
  
  cluster_transient <- cluster_transient %>%
    filter(Sample.name %in% intersect(unique(cluster_transient$Sample.name), unique(cluster_persist$Sample.name)))
  
  cluster_transient_all <- rbind(cluster_transient_all, as.data.frame(cluster_transient))
  
  # Cast for persisent NMDS 
  cluster_persist_cast <- dcast(cluster_persist, ID ~ plasmid_name_cluster, sum, value.var = "cov_depth") 
  rownames(cluster_persist_cast) <- cluster_persist_cast$ID
  cluster_persist_cast <- cluster_persist_cast[,names(cluster_persist_cast) != "ID"]
  
  # Run NMDS
  set.seed(1)
  nmds_persist <- metaMDS(cluster_persist_cast, distance = "bray", k = 2, trymax = 20)
  df_nmds_persist <- as.data.frame(nmds_persist$points)
  df_nmds_persist$ID <- row.names(df_nmds_persist)
  df_nmds_persist$Sample.name <- as.character(sapply(df_nmds_persist$ID, function(x) metadata$Sample.name[metadata$ID == x]))
  df_nmds_persist$Location <- sapply(df_nmds_persist$ID, function(x) metadata$Location[metadata$ID == x])
  
  persist_sample <- data.frame(Sample.name = sapply(rownames(cluster_persist_cast), function(x) unique(cluster_persist$Sample.name[cluster_persist$ID == x]))) 
  permanova_persist[[i]] <- adonis(cluster_persist_cast ~ Sample.name, data = persist_sample)
  
  # Summarise number of samples and clusters
  n_samples[i] <- length(unique(df_nmds_persist$Sample.name))
  n_persist_clusters[i] <- length(unique(cluster_persist$plasmid_name_cluster))
  
  # Plot NMDS
  persist_graphs[[i]] <- ggplot(df_nmds_persist, aes(MDS1, MDS2)) +
    geom_point(alpha = 0.5, colour = cols[unique(plasmid_prof_long$sample_type)[i]]) +
    geom_encircle(aes(fill = Sample.name), alpha=0.3, expand = 0) +
    theme_classic() +
    theme(text = element_text(size = 16)) +
    guides(fill = FALSE) +
    scale_fill_manual(values = rep(cols[unique(plasmid_prof_long$sample_type)[i]], length(unique(df_nmds_persist$Sample.name)))) +
    ggtitle(unique(plasmid_prof_long$sample_type)[i])
  
  # Cast for transient NMDS
  cluster_transient_cast <- dcast(cluster_transient, ID ~ plasmid_name_cluster, sum, value.var = "cov_depth") 
  rownames(cluster_transient_cast) <- cluster_transient_cast$ID
  cluster_transient_cast <- cluster_transient_cast[,names(cluster_transient_cast) != "ID"]
  
  # Run NMDS
  set.seed(1)
  nmds_transient <- metaMDS(cluster_transient_cast, distance = "bray", k = 2, trymax = 20)
  df_nmds_transient <- as.data.frame(nmds_transient$points)
  df_nmds_transient$ID <- row.names(df_nmds_transient)
  df_nmds_transient$Sample.name <- as.character(sapply(df_nmds_transient$ID, function(x) metadata$Sample.name[metadata$ID == x]))
  df_nmds_transient$Location <- sapply(df_nmds_transient$ID, function(x) metadata$Location[metadata$ID == x])
  
  # Plot NMDS
  transient_graphs[[i]] <- ggplot(df_nmds_transient, aes(MDS1, MDS2)) +
    geom_point(alpha = 0.5, colour = cols[unique(plasmid_prof_long$sample_type)[i]]) +
    geom_encircle(aes(fill = Sample.name), alpha=0.3, expand = 0) +
    theme_classic() +
    theme(text = element_text(size = 16)) +
    guides(fill = FALSE) +
    scale_fill_manual(values = rep(cols[unique(plasmid_prof_long$sample_type)[i]], length(unique(df_nmds_transient$Sample.name)))) +
    ggtitle(unique(plasmid_prof_long$sample_type)[i])
  
  # Summarise number of clusters
  n_transient_clusters[i] <- length(unique(cluster_transient$plasmid_name_cluster))
  
  transient_sample <- data.frame(Sample.name = sapply(rownames(cluster_transient_cast), function(x) unique(cluster_transient$Sample.name[cluster_transient$ID == x]))) 
  permanova_transient[[i]] <- adonis(cluster_transient_cast ~ Sample.name, data = transient_sample)
}

# Number of individuals
n_individuals <- cluster_persist_all %>% 
  group_by(sample_type) %>%
  summarise(n = n_distinct(Sample.name))

# Number of persistent clusters
n_persist_clusters <- cluster_persist_all %>% 
  group_by(sample_type) %>%
  summarise(n = n_distinct(plasmid_name_cluster))

# Number of transient clusters
n_transient_clusters <- cluster_transient_all %>% 
  group_by(sample_type) %>%
  summarise(n = n_distinct(plasmid_name_cluster))

tiff("figures/persistent_plasmid_clusters.tiff", width = 2000, height = 1500, res = 200)
grid.arrange(grobs = persist_graphs, layout_matrix = rbind(c(1,3), c(2,4)))
dev.off()

tiff("figures/transient_plasmid_clusters.tiff", width = 2000, height = 1500, res = 200)
grid.arrange(grobs = transient_graphs, layout_matrix = rbind(c(1,3), c(2,4)))
dev.off()

# Percentage of individuals where plasmid cluster is persistent, transient or doesn't exist
long_n <- metadata_longus %>%
  group_by(sample_type) %>%
  summarise(n = n_distinct(Sample.name))

cluster_all <- rbind(cluster_persist_all, cluster_transient_all)

cluster_summary <- cluster_no_tp %>%
  filter(paste0(Sample.name, sample_type) %in% unique(paste0(cluster_all$Sample.name, cluster_all$sample_type))) %>%
  group_by(sample_type, plasmid_name_cluster) %>%
  summarise(persistent = n_distinct(Sample.name[no_tp >= 2]), transient = n_distinct(Sample.name[no_tp < 2])) %>%
  left_join(long_n, by = "sample_type") %>%
  mutate(perc_pers = persistent/n*100, perc_trans = transient/n*100) %>%
  mutate(perc_na = (n - (persistent + transient))/n*100)

cluster_summary_melt <- rbind(data.frame(sample_type = cluster_summary$sample_type, plasmid_name_cluster = cluster_summary$plasmid_name_cluster, 
                                         perc = cluster_summary$perc_pers, group = "Persistent"),
                              data.frame(sample_type = cluster_summary$sample_type, plasmid_name_cluster = cluster_summary$plasmid_name_cluster, 
                                         perc = cluster_summary$perc_trans, group = "Transient"),
                              data.frame(sample_type = cluster_summary$sample_type, plasmid_name_cluster = cluster_summary$plasmid_name_cluster, 
                                         perc = cluster_summary$perc_na, group = "Absent"))
cluster_summary_melt$group <- factor(cluster_summary_melt$group, levels = c("Absent", "Transient", "Persistent"))

# Plot persistency of plasmid clusters
cluster_summary_graphs <- list()
cluster_ranks <- data.frame()
for (i in 1:length(unique(cluster_summary$sample_type))) {
  cluster_sample_type <- cluster_summary_melt[cluster_summary_melt$sample_type == unique(cluster_summary$sample_type)[i],]
  cluster_sample_type$plasmid_name_cluster <- factor(cluster_sample_type$plasmid_name_cluster,
                                                 levels = cluster_sample_type$plasmid_name_cluster[cluster_sample_type$group == "Absent"][order(cluster_sample_type$perc[cluster_sample_type$group == "Absent"],
                                                                                                                                            cluster_sample_type$perc[cluster_sample_type$group == "Transient"])])
  
  cluster_summary_graphs[[i]] <- ggplot(cluster_sample_type, aes(plasmid_name_cluster, perc, fill = group)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    theme(text = element_text(size = 16)) +
    scale_fill_manual("Status", values = c("white", "lightblue", "blue"), labels = c("Absent", "Transient", "Persistent")) +
    ylab("% Individuals") + xlab("Plasmid clusters/singletons") +
    theme(axis.text.x = element_blank()) + 
    ggtitle(unique(cluster_summary$sample_type)[i])
  
}
tiff("figures/plasmid_pers_trans_perc.tiff", width = 3000, height = 1000, res = 200)
grid.arrange(grobs = cluster_summary_graphs, layout_matrix = rbind(c(1,2),c(3,4)))
dev.off()

#### RESISTANCE PLASMIDS ####
# Combine ARG results
arg_res <- data.frame()
for (i in 1:length(metadata$ID)) {
  if (file.size(paste0("data/card/", metadata$ID[i], "_plasmid_card.out"))) {
    arg_tmp <- read.delim(paste0("data/card/", metadata$ID[i], "_plasmid_card.out"), stringsAsFactors = FALSE, header = FALSE)
    arg_tmp$ID <- metadata$ID[i]
    arg_res <- rbind(arg_res, arg_tmp)
  }
}

blast_colnames <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
names(arg_res) <- c(blast_colnames, "ID")

# Select hits
catalog_args <- read.delim("data/card/plasmid_nonredundant_catalog_card.out", stringsAsFactors = FALSE, header = FALSE)
names(catalog_args) <- blast_colnames
catalog_args_dup <- catalog_args[catalog_args$qseqid %in% catalog_args$qseqid[duplicated(catalog_args$qseqid)],]
catalog_args_sig <- catalog_args[!catalog_args$qseqid %in% catalog_args_dup$qseqid,]
remove_ids <- unlist(lapply(unique(catalog_args_dup$qseqid), function(x) getBestBlastHit(catalog_args_dup, x)))
catalog_args <- catalog_args_dup[-remove_ids,]
catalog_args <- rbind(catalog_args, catalog_args_sig)
arg_res <- left_join(catalog_args, arg_res)
  
# Add ARG information
card_aro_categories_index <- read.csv("db/CARD_DB/card-data/aro_categories_index.csv", stringsAsFactors = FALSE, sep = "\t")
card_index <- read.csv("db/CARD_DB/card-data/aro_index.csv", stringsAsFactors = FALSE, sep = "\t")
card_metadata <- full_join(card_index, card_aro_categories_index, by = "Protein.Accession")
card_metadata <- card_metadata[!duplicated(card_metadata$ARO.Accession),]
arg_res$ARO.Accession <- paste0("ARO:", gsub("_.*", "", gsub(".*_ARO:", "", arg_res$sseqid)))
arg_res <- left_join(arg_res, card_metadata, by = "ARO.Accession")
arg_res <- arg_res %>% select(ID, qseqid, pident, evalue, qstart, qend, bitscore, ARO.Name, AMR.Gene.Family, Drug.Class, Resistance.Mechanism) 
arg_res$Drug.Class <- sapply(arg_res$Drug.Class, function(x) paste(unique(strsplit(x, ",")[[1]]), collapse = ","))
arg_res$AMR.Gene.Family = sapply(arg_res$AMR.Gene.Family, function(x) paste(unique(strsplit(x, ",")[[1]]), collapse = ","))
arg_res$ARO.Name = sapply(arg_res$ARO.Name, function(x) paste(unique(strsplit(x, ",")[[1]]), collapse = ","))
arg_res$Resistance.Mechanism = sapply(arg_res$Resistance.Mechanism, function(x) paste(unique(strsplit(x, ",")[[1]]), collapse = ","))
arg_res$Drug.Class.alt <- arg_res$Drug.Class
arg_res$Drug.Class.alt[sapply(arg_res$Drug.Class.alt, function(x) str_count(x, ";")) > 1] <- "multidrug"
arg_res$Drug.Class.alt[arg_res$Resistance.Mechanism == "antibiotic efflux"] <- paste(arg_res$Drug.Class.alt[arg_res$Resistance.Mechanism == "antibiotic efflux"], "efflux")
arg_res$Drug.Class.alt <- gsub(";", " and ", arg_res$Drug.Class.alt)
arg_res <- arg_res %>% 
  group_by(ID, qseqid) %>%
  summarise_all(paste, collapse = ",") %>%
  select(qseqid, ARO.Name, AMR.Gene.Family, Drug.Class.alt, Resistance.Mechanism)

# Combine plasmid and ARG results
plasmid_args <- inner_join(plasmid_prof_sing, arg_res, by = c("ID", "name"="qseqid"), suffix = c(".plasmid", ".arg"))

# Number of unique plasmid clusters/singletons
n_uniq_plasmid_red <- length(unique(plasmid_prof$name))

# Number of unique plasmid clusters
n_uniq_clusters <- length(unique(plasmid_prof$plasmid_cluster)) - 1

# Number of unique plasmids in plasmid clusters
n_uniq_plasmid_cluster <- length(unique(plasmid_prof$name[!is.na(plasmid_prof$plasmid_cluster)]))

# Number of unique plasmid clusters/singletons
n_uniq_plasmids <- length(unique(plasmid_prof$plasmid_name_cluster))

# Proportion of unique plasmids with ARGs
n_uniq_plasmids_args <- length(unique(plasmid_args$plasmid_name_cluster))/length(unique(plasmid_prof$plasmid_name_cluster))

# Participant overview
plasmid_args %>% group_by(Location) %>%
  summarise(n_distinct(Sample.name))

# Prevalence of resistance plasmids
plasmid_args %>% 
  group_by(Location, sample_type, name) %>%
  summarise(n_args = n_distinct(ARO.Name))

# Prevalence of each ARG (for each plasmid)
plasmid_args %>% 
  group_by(Location, sample_type, ARO.Name, Drug.Class.alt, Resistance.Mechanism) %>%
  summarise(n = n_distinct(name))

plasmid_args_groups <- plasmid_args %>%
  group_by(plasmid_name_cluster, av_size, genome, type, Location, Sample.name, ARO.Name, Drug.Class.alt, Resistance.Mechanism) %>%
  summarise(n_sample_types = n_distinct(sample_type), sample_type_group = paste(sort(unique(sample_type)), collapse = "\n"))

# Create data for plasmid figure 1
plasmid_args <- left_join(plasmid_args, plasmid_args_groups, by = c("plasmid_name_cluster", "av_size", "genome", "type", "Location", "ARO.Name", "Drug.Class.alt", "Resistance.Mechanism", "Sample.name"))

args_summary_pnts <- plasmid_args %>%
  group_by(plasmid_name_cluster, av_size, genome, type, Location_sampletype, Location, sample_type, n_sample_types, sample_type_group, ARO.Name, Drug.Class.alt, Resistance.Mechanism) %>%
  summarise(labels = as.character(n_distinct(Sample.name))) %>%
  ungroup() %>% group_by(sample_type, n_sample_types, sample_type_group, ARO.Name, Drug.Class.alt, Resistance.Mechanism)

point_loc <- data.frame(Location_sample_type_group = sort(unique(paste0(args_summary_pnts$sample_type_group, "-", args_summary_pnts$Location))))
point_loc$loc <- c(1:nrow(point_loc))
point_loc$Location <- sapply(as.character(point_loc$Location_sample_type_group), function(x) strsplit(x, "-")[[1]][2])
point_loc$sample_type_group <- sapply(as.character(point_loc$Location_sample_type_group), function(x) strsplit(x, "-")[[1]][1])

point_loc$num_per_group <- sapply(as.character(point_loc$Location_sample_type_group), function(x) {
  location <- strsplit(x, "-")[[1]][2]
  sample_types <- strsplit(strsplit(x, "-")[[1]][1],"\n")[[1]]
  num_per_group <- metadata[metadata$sample_type %in% sample_types & metadata$Location %in% location,] %>%
    group_by(Sample.name) %>%
    summarise(n = n()) %>%
    filter(n >= length(sample_types)) %>%
    nrow()
  return(num_per_group)
})

point_loc$Location_sample_type_group <- gsub("-", "\n", point_loc$Location_sample_type_group)
point_loc$label <- paste0(point_loc$Location_sample_type_group, "\nn = ", as.character(point_loc$num_per_group))

args_summary_pnts <- left_join(args_summary_pnts, point_loc, by = c("Location", "sample_type_group"))

args_summary_pnts$x_labels <- paste0(args_summary_pnts$ARO.Name, "\n(", args_summary_pnts$Drug.Class.alt, ")")

# Colours
cohort_cols <- c("grey", brewer.pal(9, "Blues")[c(5,7)], "gold", brewer.pal(9, "YlOrRd")[c(5,7)], brewer.pal(9, "RdPu")[c(3,5)])
names(cohort_cols) <- sort(unique(metadata$Location_sampletype)) 

# Summary graph
tiff("figures/resistant_plasmids.tiff", height = 2500, width = 4000, res = 150)
ggplot(args_summary_pnts, aes(sample_type, loc)) +
  geom_label(aes(label = type), size = 2, nudge_x = -0.4, alpha = 0) +
  geom_point(aes(colour = Location_sampletype, size = log10(av_size))) +
  coord_flip() +
  geom_text(aes(label = labels), colour = "black", size = 3) +
  facet_grid(plasmid_name_cluster + x_labels ~ ., scale = "free", space = "free", switch = "both") +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 180, size = 6),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_colour_manual("Body Site - Country", 
                      values = cohort_cols[names(cohort_cols) %in% args_summary_pnts$Location_sampletype], 
                      labels = names(cohort_cols)[names(cohort_cols) %in% args_summary_pnts$Location_sampletype]) +
  scale_size(name = "log10(mean size [bp])") +
  xlab("") + ylab("Shared GIT sites") + 
  scale_y_continuous(breaks = c(1:nrow(point_loc)), labels = as.character(point_loc$label)) +
  guides(colour = guide_legend(override.aes = list(size=5)))
dev.off()


#### LONGITUDINAL RESISTANCE PLASMIDS ####
# Combine plasmid and ARG results
plasmid_args_long <- inner_join(plasmid_prof_long, arg_res, by = c("ID", "name"="qseqid"), suffix = c(".plasmid", ".arg"))

# Number of unique plasmid clusters/singletons
n_uniq_plasmid_red_long <- length(unique(plasmid_prof_long$name))

# Number of unique plasmid clusters
n_uniq_clusters_long <- length(unique(plasmid_prof_long$plasmid_cluster)) - 1

# Number of unique plasmids in plasmid clusters
n_uniq_plasmid_cluster_long <- length(unique(plasmid_prof_long$name[!is.na(plasmid_prof_long$plasmid_cluster)]))

# Number of unique plasmid clusters/singletons
n_uniq_plasmids_long <- length(unique(plasmid_prof_long$plasmid_name_cluster))

# Proportion of unique plasmids with ARGs
n_uniq_plasmids_args_long <- length(unique(plasmid_args_long$plasmid_name_cluster))/length(unique(plasmid_prof_long$plasmid_name_cluster))

# Participant overview
plasmid_args_long %>% group_by(Location) %>%
  summarise(n_distinct(Sample.name))

# Prevalence of resistance plasmids
plasmid_args_long %>% 
  group_by(Location, sample_type, name) %>%
  summarise(n_args = n_distinct(ARO.Name))

# Prevalence of each ARG (for each plasmid)
plasmid_args_long %>% 
  group_by(Location, sample_type, ARO.Name, Drug.Class.alt, Resistance.Mechanism) %>%
  summarise(n = n_distinct(name))

# Count average time points of plasmids
cluster_args_no_tp <- plasmid_args_long %>% 
  mutate(plasmid_arg = paste(plasmid_name_cluster, "-", ARO.Name)) %>%
  group_by(Location, sample_type, Sample.name, plasmid_arg, Drug.Class.alt, Resistance.Mechanism) %>%
  summarise(no_tp = n_distinct(timepoint)) 

# Count all distinct plasmid - args
cluster_args_summary <- plasmid_args_long %>% 
  mutate(plasmid_arg = paste(plasmid_name_cluster, "-", ARO.Name)) %>%
  group_by(sample_type, Sample.name) %>%
  summarise(total_plasmids = n_distinct(plasmid_arg))

# Count plasmid - args in at least x number of time points
cluster_args_count_tp <- map_df(.x = unique(cluster_args_no_tp$no_tp), .f = function(.x) {
  tmp <- data.frame(sample_type = rep(cluster_args_no_tp$sample_type[cluster_args_no_tp$no_tp == .x], .x),
                    Sample.name = rep(cluster_args_no_tp$Sample.name[cluster_args_no_tp$no_tp == .x], .x),
                    plasmid_arg = rep(cluster_args_no_tp$plasmid_arg[cluster_args_no_tp$no_tp == .x]),
                    no_tp = rep(cluster_args_no_tp$no_tp[cluster_args_no_tp$no_tp == .x], .x))
  tmp$tp <- c(1:.x)
  return(tmp)
})
  
# Proportion of plasmid singletons/clusters in timepoints
cluster_count_tp_summary <- cluster_args_count_tp %>% 
  group_by(sample_type, Sample.name, tp) %>%
  summarise(total_cluster_least = n_distinct(plasmid_arg))

cluster_one <- cluster_count_tp_summary %>%
  filter(tp == 1) %>%
  rename(total_cluster = total_cluster_least) %>%
  select(-tp)

cluster_count_tp_summary <- left_join(cluster_count_tp_summary, cluster_one, by = c("sample_type", "Sample.name")) %>%
  mutate(cluster_frac = total_cluster_least/total_cluster)

# Plot
tiff("figures/longitudinal_plasmid_args.tiff", height = 600, width = 1500, res = 100)
set.seed(1)
ggplot(cluster_args_no_tp, aes(as.factor(no_tp), cluster_frac, fill = as.factor(sample_type))) +
  geom_boxplot() +
  geom_jitter() +
  facet_grid(~ sample_type, scale = "free", space = "free") +
  theme_bw() + xlab("No. timepoints") + ylab("Proportion of plasmid singletons/clusters") +
  theme(strip.background = element_blank(), strip.text.x = element_blank(), text = element_text(size = 16)) +
  scale_fill_manual("GIT Site", values = cols)
dev.off()  

# Samples with 3 timepoints only
id_3no_tp <- metadata_longus %>%
  group_by(sample_type, Sample.name) %>%
  summarise(no_tp = n_distinct(timepoint)) %>%
  filter(no_tp >= 3) %>%
  inner_join(metadata) %>% ungroup()

# Summary of no. samples
id_3no_tp %>% group_by(sample_type, timepoint) %>%
  summarise(n_distinct(Sample.name))
  
# Persistent and transient plasmid singletons/clusters
plasmid_args_long_summary <- cluster_args_no_tp %>%
  inner_join(metadata_longus) %>%
  filter(ID %in% id_3no_tp$ID) %>%
  filter(no_tp %in% c(1,2,3)) %>%
  group_by(sample_type, plasmid_arg, no_tp) %>%
  summarise(n = n_distinct(Sample.name))

# Plot
tiff("figures/longitudinal_three_timepoints.tiff", height = 800, width = 1500, res = 120)
ggplot(plasmid_args_long_summary, aes(as.character(no_tp), n, fill = sample_type)) +
  geom_bar(stat = "identity") +
  facet_grid(~plasmid_arg, scale = "free", space = "free") +
  theme_bw() +
  scale_fill_manual("GIT Site", values = cols) +
  xlab("No. timepoints") + ylab("No. individuals") +
  scale_y_continuous(breaks = seq(0, 80, 5)) +
  theme(strip.text.x = element_text(angle = 90, size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
dev.off()

#### Taxonomic compositions ####
metaphlan_fam <- metaphlan[,colnames(metaphlan) %in% c("X", metadata$ID[metadata$timepoint == 1])]
rownames(metaphlan_fam) <- metaphlan_fam$X
metaphlan_fam <- metaphlan_fam[,-1]
metaphlan_fam <- metaphlan_fam[grepl("f__", rownames(metaphlan_fam)) & !grepl("g__", rownames(metaphlan_fam)),]
metaphlan_fam <- metaphlan_fam[rowSums(metaphlan_fam) > 1,]
family_labels <- gsub("_", " ", gsub(".*f__", "", rownames(metaphlan_fam)))

tiff("figures/metaphlan_family.tiff", width = 1500, height = 1500, res = 150)
heatmap.2(as.matrix(metaphlan_fam),
          margins = c(10,20),
          trace = "none",
          scale = "none",
          hclustfun = function(x) {hclust(x, method = "ward.D2")},
          dendrogram = "column",
          col =  colorRampPalette(brewer.pal(9, "PuRd")[c(1, 9)])(50),
          breaks = seq(min(metaphlan_fam), max(metaphlan_fam), length.out = 51),
          symbreaks = FALSE,
          keysize = 1, 
          #lhei = c(1,8),
          key.title = NA, key.xlab = "rel. abundance", key.ylab = NA,
          density.info = "none",
          ColSideColors = sapply(metadata[colnames(metaphlan_fam), "Location_sampletype"], function(x) cohort_cols[x]),
          labRow = as.expression(lapply(family_labels, function(a) bquote(italic(.(a))))),
          labCol = NA,
          cexCol = 0.5,
          ylab = "Family",
          main = NA
)
legend(x = 0.85, y = 1.05, xpd=TRUE, legend = levels(as.factor(metadata[colnames(metaphlan_fam), "Location_sampletype"])),
       col = cohort_cols[names(cohort_cols) %in% sort(unique(metadata[colnames(metaphlan_fam), "Location_sampletype"]))], bg = "white", box.col = "black",
       lty = 1, lwd = 5, cex = 0.5, title = "GIT Sites")
dev.off()

#### NMDS of metaphlan data
metaphlan_gen <- metaphlan[,colnames(metaphlan) %in% c("X", metadata$ID[metadata$timepoint == 1])]
rownames(metaphlan_gen) <- metaphlan_gen$X
metaphlan_gen <- metaphlan_gen[,-1]
metaphlan_gen <- metaphlan_gen[grepl("g__", rownames(metaphlan_gen)) & !grepl("s__", rownames(metaphlan_gen)),]
metaphlan_gen <- metaphlan_gen[grepl("k__Bacteria", row.names(metaphlan_gen)) | grepl("k__Archaea", row.names(metaphlan_gen)),]
metaphlan_gen <- metaphlan_gen[!grepl("unclassified", row.names(metaphlan_gen)),]
metaphlan_gen <- metaphlan_gen[!grepl("noname", row.names(metaphlan_gen)),]

metaphlan_nmds <- metaMDS(t(metaphlan_gen), distance = "bray", k = 2, trymax = 20)
df_metaphlan_nmds <- as.data.frame(metaphlan_nmds$points)
df_metaphlan_nmds$ID <- names(metaphlan_gen)
names(df_metaphlan_nmds) = c("MDS1", "MDS2", "ID")
df_metaphlan_nmds <- left_join(df_metaphlan_nmds, metadata, by = "ID")

# Procrustes analysis and plot
protest_res <- protest(df_cluster_samples_nmds[,c(1:2)], df_metaphlan_nmds[,c(1:2)], scale = TRUE)
protest_res$t0
tiff("figures/procrustes.tiff", height = 1000, width = 1200, res = 200)
plot(protest_res, cex = 0.5, ar.col = "grey", )
points(protest_res, display = "target", col = "red", cex = 0.5)
dev.off()

#### Enterobacteriaceae composition for individuals with ColE1 ####
entrb <- metaphlan[,colnames(metaphlan) %in% c("X", metadata$ID[metadata$Sample.name %in% unique(plasmid_args$Sample.name[plasmid_args$genome %in% c("ColRNAI_1__DQ298019")])])]
rownames(entrb) <- entrb$X
entrb <- entrb[,-1]
entrb <- entrb[grep("f__Enterobacteriaceae", rownames(entrb)),]
entrb_total <- c(as.matrix(entrb[!grepl("g__", rownames(entrb)),]))
#entrb <- entrb[,entrb_total != 0]
#entrb_total <- entrb_total[entrb_total != 0]
entrb <- entrb[grepl("s__", rownames(entrb)) & !grepl("t__", rownames(entrb)),]
entrb <- entrb[rowSums(entrb) != 0,]
entrb <- sweep(as.matrix(entrb), 2, entrb_total, "/")
entrb[is.nan(entrb)] <- 0
species_labels <- gsub("_", " ", gsub(".*s__", "", rownames(entrb)))

tiff("figures/enterobacteriaceae_species.tiff", width = 1000, height = 750, res = 150)
heatmap.2(as.matrix(entrb),
          margins = c(10,20),
          trace = "none",
          scale = "none",
          hclustfun = function(x) {hclust(x, method = "ward.D2")},
          dendrogram = "column",
          col =  colorRampPalette(brewer.pal(9, "PuRd")[c(1, 9)])(50),
          breaks = seq(min(entrb), max(entrb), length.out = 51),
          symbreaks = FALSE,
          keysize = 1, 
          #lhei = c(1,8),
          key.title = NA, key.xlab = "rel. abundance", key.ylab = NA,
          density.info = "none",
          ColSideColors = sapply(metadata[colnames(entrb), "sample_type"], function(x) cols[x]),
          labRow = as.expression(lapply(species_labels, function(a) bquote(italic(.(a))))),
          labCol = NA,
          cexCol = 0.5,
          ylab = "Species from the Enterobacteriaceae family",
          main = NA
)
legend(x = 0.85, y = 1.05, xpd=TRUE, legend = levels(as.factor(metadata[colnames(entrb), "sample_type"])),
       col = cols[names(cols) %in% sort(unique(metadata[colnames(entrb), "sample_type"]))], bg = "white", box.col = "black",
       lty = 1, lwd = 5, cex = 0.5, title = "GIT Sites")
dev.off()
