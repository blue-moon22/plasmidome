# Libraries
library(dplyr)
library(purrr)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(purrr)

metadata <- read.csv("data/supplementary_metadata.csv", stringsAsFactors = FALSE)
samples <- read.delim("data/samples_analysed.txt", stringsAsFactors = FALSE, header = FALSE)
metadata <- metadata[metadata$ID %in% samples$V1 & metadata$Location %in% c("China", "USA"),]

# Summarise
metadata %>% group_by(Location, sample_type) %>%
  summarise(n = n_distinct(ID))
metadata$Location_sampletype <- paste(metadata$sample_type, "-", metadata$Location)

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

# Phage breadth coverage
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

plasmid_prof <- left_join(plasmid_prof, plasmid_taxa, by = "V1") %>%
  rename("name"="V1") %>%
  left_join(metadata, by = "ID")

# Combine ARG results
arg_res <- data.frame()
for (i in 1:length(metadata$ID)) {
  if (file.size(paste0("data/card/", metadata$ID[i], "_card.out"))) {
    arg_tmp <- read.delim(paste0("data/card/", metadata$ID[i], "_card.out"), stringsAsFactors = FALSE, header = FALSE)
    arg_tmp$ID <- metadata$ID[i]
    arg_res <- rbind(arg_res, arg_tmp)
  }
}

# Add ARG information
card_aro_categories_index <- read.csv("db/CARD_DB/card-data/aro_categories_index.csv", stringsAsFactors = FALSE, sep = "\t")
card_index <- read.csv("db/CARD_DB/card-data/aro_index.csv", stringsAsFactors = FALSE, sep = "\t")
card_metadata <- full_join(card_index, card_aro_categories_index, by = "Protein.Accession")
card_metadata <- card_metadata[!duplicated(card_metadata$ARO.Accession),]

arg_res$ARO.Accession <- paste0("ARO:", gsub("_.*", "", gsub(".*_ARO:", "", arg_res$V2)))
arg_res <- left_join(arg_res, card_metadata, by = "ARO.Accession")
arg_res <- arg_res %>% select(V1, V3, V11, V12, ARO.Name, AMR.Gene.Family, Drug.Class, Resistance.Mechanism) %>%
  rename("identity"="V3", "evalue"="V11", "bitscore"="V12")

# Select hits
arg_res <- arg_res %>% group_by(V1) %>%
  filter(evalue == min(evalue)) %>%
  filter(bitscore == max(bitscore)) %>% summarise_all(paste, collapse = ",") %>%
  select(V1, ARO.Name, AMR.Gene.Family, Drug.Class, Resistance.Mechanism)

arg_res$Drug.Class <- sapply(arg_res$Drug.Class, function(x) paste(unique(strsplit(x, ",")[[1]]), collapse = ","))
arg_res$AMR.Gene.Family = sapply(arg_res$AMR.Gene.Family, function(x) paste(unique(strsplit(x, ",")[[1]]), collapse = ","))
arg_res$ARO.Name = sapply(arg_res$ARO.Name, function(x) paste(unique(strsplit(x, ",")[[1]]), collapse = ","))
arg_res$Resistance.Mechanism = sapply(arg_res$Resistance.Mechanism, function(x) paste(unique(strsplit(x, ",")[[1]]), collapse = ","))

# Combine plasmid and ARG results
plasmid_args <- inner_join(plasmid_prof, arg_res, by = c("name"="V1"), suffix = c(".plasmid", ".arg"))

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
  group_by(Location, sample_type, ARO.Name, Drug.Class, Resistance.Mechanism) %>%
  summarise(n = n_distinct(name))

plasmid_args_groups <- plasmid_args %>%
  group_by(plasmid_name_cluster, av_size, genome, type, Location, Sample.name, ARO.Name, Drug.Class, Resistance.Mechanism) %>%
  summarise(n_sample_types = n_distinct(sample_type), sample_type_group = paste(sort(unique(sample_type)), collapse = "\n"))

# Create data for plasmid figure 1
plasmid_args <- left_join(plasmid_args, plasmid_args_groups, by = c("plasmid_name_cluster", "av_size", "genome", "type", "Location", "ARO.Name", "Drug.Class", "Resistance.Mechanism", "Sample.name"))

args_summary_pnts <- plasmid_args %>%
  group_by(plasmid_name_cluster, av_size, genome, type, Location_sampletype, Location, sample_type, n_sample_types, sample_type_group, ARO.Name, Drug.Class, Resistance.Mechanism) %>%
  summarise(labels = as.character(n_distinct(Sample.name))) %>%
  ungroup() %>% group_by(sample_type, n_sample_types, sample_type_group, ARO.Name, Drug.Class, Resistance.Mechanism)

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

args_summary_pnts$Drug.Class.alt <- args_summary_pnts$Drug.Class
args_summary_pnts$Drug.Class.alt[sapply(args_summary_pnts$Drug.Class.alt, function(x) str_count(x, ";")) > 1] <- "multidrug"
args_summary_pnts$Drug.Class.alt[args_summary_pnts$Resistance.Mechanism == "antibiotic efflux"] <- paste(args_summary_pnts$Drug.Class.alt[args_summary_pnts$Resistance.Mechanism == "antibiotic efflux"], "efflux")
args_summary_pnts$Drug.Class.alt <- gsub(";", " and\n", args_summary_pnts$Drug.Class.alt)

args_summary_pnts$x_labels <- paste0(args_summary_pnts$ARO.Name, "\n(", args_summary_pnts$Drug.Class.alt, ")")

# Colours
cohort_cols <- c("grey", brewer.pal(9, "Blues")[c(5,7)], "gold", brewer.pal(9, "YlOrRd")[c(5,7)], brewer.pal(9, "RdPu")[c(3,5)])
names(cohort_cols) <- sort(unique(metadata$Location_sampletype)) 

# Summary graph
tiff("figures/resistant_plasmids.tiff", height = 2000, width = 3500, res = 150)
ggplot(args_summary_pnts, aes(sample_type, loc)) +
  geom_label(aes(label = type), size = 2, nudge_x = -0.4, alpha = 0) +
  geom_point(aes(colour = Location_sampletype, size = log10(av_size))) +
  coord_flip() +
  geom_text(aes(label = labels), colour = "black", size = 3) +
  facet_grid(x_labels + plasmid_name_cluster ~ ., scale = "free", space = "free", switch = "both") +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 180, size = 6),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_colour_manual("Body Site - Geographical Location", 
                      values = cohort_cols[names(cohort_cols) %in% args_summary_pnts$Location_sampletype], 
                      labels = names(cohort_cols)[names(cohort_cols) %in% args_summary_pnts$Location_sampletype]) +
  scale_size(name = "log10(size [bp])") +
  xlab("") + ylab("Shared GIT sites") + 
  scale_y_continuous(breaks = c(1:nrow(point_loc)), labels = as.character(point_loc$label)) +
  guides(colour = guide_legend(override.aes = list(size=5)))
dev.off()
