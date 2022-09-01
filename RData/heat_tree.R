library(metacoder)

tmp3<-subset_samples(physeq = tmp2,Host=="Infant") 
tmp3<-prune_taxa(taxa_sums(tmp3) > 1, tmp3) 

tmp3
tmp3<-aggregate_taxa(x =tmp3,level = "Genus",verbose = T)
tmp3
tmp3<-subset_taxa(physeq = tmp3,Genus!="Unknown")
tmp3
sam<-data.frame(sample_data(tmp3))
sam<-sam%>%mutate(controlled=HbA1C_1trym_2<7.2)


obj<-parse_phyloseq(tmp3)
obj
obj$data$otu_table <- zero_low_counts(obj, data = "otu_table", min_count = 5000)
no_reads <- rowSums(obj$data$otu_table[, sam$SampleID]) == 0
sum(no_reads)
obj <- filter_obs(obj, data = "otu_table", ! no_reads, drop_taxa = TRUE)
print(obj)
obj$data$tax_data <- calc_obs_props(obj, "otu_table")
obj$data$tax_abund <- calc_taxon_abund(obj, "otu_table",
                                       cols = sam$SampleID)
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = sam$SampleType, cols = sam$SampleID)
set.seed(314) # This makes the plot appear the same each time it is run 


obj$data$diff_table <- compare_groups(obj,
                                      data = "tax_abund",
                                      cols = sam$SampleID, # What columns of sample data to use
                                      groups = sam$controlled) # What category each sample is assigned to
print(obj$data$diff_table)


set.seed(999)
# heat_tree(obj, 
#           node_label = taxon_names,
#           node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
#           node_color = log2_median_ratio, # A column from `obj$data$diff_table`
#           node_color_interval = c(-2, 2), # The range of `log2_median_ratio` to display
#           node_color_range = c("cyan", "gray", "tan"), # The color palette used
#           node_size_axis_label = "OTU count",
#           node_color_axis_label = "Log 2 ratio of median proportions",
#           layout = "davidson-harel", # The primary layout algorithm
#           initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations
# 
# 

obj$data$diff_table <- compare_groups(obj, dataset = "tax_abund",
                                      cols = sam$SampleID, # What columns of sample data to use
                                      groups = paste0(sam$disease,"_",sam$SampleType)) # What category each sample is assigned to
print(obj$data$diff_table)

desc(unique(obj$data$diff_table$log2_median_ratio))

set.seed(1)
heat_tree_matrix(obj,
                 data = "diff_table",
                 node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
                 node_label = taxon_names,
                 node_color = log2_median_ratio, # A column from `obj$data$diff_table`
                 node_color_range = diverging_palette(), # The built-in palette for diverging data
                 node_color_trans = "linear", # The default is scaled by circle area
                 node_color_interval = c(-0.3, 0.3), # The range of `log2_median_ratio` to display
                 edge_color_interval = c(-0.3, 0.3), # The range of `log2_median_ratio` to display
                 node_size_axis_label = "Number of OTUs",
                 key_size = 0.65,
                 row_label_size = 24,
                 col_label_size = 24,
                 label_small_trees = F,
                 seed = 314,
                 node_color_axis_label = "Log2 ratio median proportions",
                 layout = "davidson-harel", # The primary layout algorithm
                 initial_layout = "reingold-tilford", # The layout algorithm that initializes node locations
                 output_file = "differential_heat_tree.pdf") # Saves the plot as a pdf file

