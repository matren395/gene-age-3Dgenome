##Katherine Fleck

library(ggplot2)
library(writexl)
library(rstatix)
library(dplyr)

######### Analysis of enhancers assigned to genes or Igen controls through domain or E-P looping definitions
#pipeline shown for GM12878, also completed with K562

setwd("/Users/katherinefleck/Library/Erceg_Lab/Data_files_for_3D_genome_eras_paper/new_CTCF_enhancer")

GM12878_counts <- read.table("GM12878_Disease_gene_dataset_enhancer_counts.txt", header = TRUE, sep = "\t")

#collapse mammal and primate genes into chordate category for disease analysis
GM12878_counts$disease_era <- ifelse(GM12878_counts$Era == "Ancient", "Ancient", 
                                     ifelse(GM12878_counts$Era == "Metazoan", "Metazoan", 
                                            ifelse(GM12878_counts$Era == "Chordate" | GM12878_counts$Era == "Mammal" | GM12878_counts$Era == "Primate", "Chordate", NA)))


# add a column for simple disease assignment, collapsing sex-linked with autosomal into just "non-disease", "recessive", and "dominant"
GM12878_counts$simple_disease <- ifelse(GM12878_counts$disease_gene_inheritance == "NonDisease" | GM12878_counts$disease_gene_inheritance == "NonDisease_SexLinked", "NonDisease", 
                                        ifelse(GM12878_counts$disease_gene_inheritance == "AR" | GM12878_counts$disease_gene_inheritance == "XLR", "Recessive", 
                                               ifelse(GM12878_counts$disease_gene_inheritance == "AD" | GM12878_counts$disease_gene_inheritance == "XLD", "Dominant", NA)))


# gene promoters that are not in domains from each cell line
GM12878_notTAD_IDs <- read.table("GM12878_gene_promoters_notin_TADs.bed")

# filter out genes from the counts files that do not have promoters in domains
GM12878_counts_filtered <- GM12878_counts[!(GM12878_counts$Protein_ID %in% GM12878_notTAD_IDs$V4),]

# bin enhancer counts in domains for plotting
GM12878_counts_filtered$TADs_binned <- cut(GM12878_counts_filtered$TAD, 
                                           breaks = c(-Inf, 0, 1, 2, 4, Inf),
                                           labels = c("0", "1", "2", "3-4", "5+"), 
                                           right = TRUE)

# bin E-P loop counts for plotting
GM12878_counts$loop_binned <- cut(GM12878_counts$loop, 
                                            breaks = c(-Inf, 1, 2, Inf),
                                            labels = c("1", "2", "3+"), 
                                            right = TRUE)


# joining Igens (with gene defined-enhancers) with genes in the same table
# load in the Igen information
GM12878_new_igens <- read.table("GM12878_igen_controls_gene_enhancer_counts.txt", header = TRUE, sep = "\t")

# add a column for control type
GM12878_new_igens$control_status <- ifelse(grepl("Non", GM12878_new_igens[[5]]), "igen_NonORF", "igen_ORF") 

# Igen promoters that are not in domains from each cell line
GM12878_notTAD_igens <- read.table("GM12878_igen_promoters_notin_TADs.bed")

# filter out Igens from the counts files that do not have promoter in domains
GM12878_new_igens_filtered <- GM12878_new_igens[!(GM12878_new_igens$ID %in% GM12878_notTAD_igens$V4),]

# do binning of counts for plotting
GM12878_new_igens_filtered$TADs_binned <- cut(GM12878_new_igens_filtered$TAD, 
                                              breaks = c(-Inf, 0, 1, 2, 4, Inf),
                                              labels = c("0", "1", "2", "3-4", "5+"), 
                                              right = TRUE)

GM12878_new_igens$loop_binned <- cut(GM12878_new_igens$loop, 
                                               breaks = c(-Inf, 1, 2, Inf),
                                               labels = c("1", "2", "3+"), 
                                               right = TRUE)


#combining genes and Igens together, done separately for loop data and domain data
GM12878_gene_columns <- GM12878_counts[ ,c("Protein_ID", "Era", "loop", "loop_binned")]
GM12878_new_igen_columns <- GM12878_new_igens[ , c("ID", "control_status", "loop", "loop_binned")]
colnames(GM12878_new_igen_columns) <- c("Protein_ID", "Era", "loop", "loop_binned")
GM12878_gene_new_controls <- rbind(GM12878_gene_columns, GM12878_new_igen_columns)

GM12878_gene_columns_filtered <- GM12878_counts_filtered[ ,c("Protein_ID", "Era", "TADs_binned")]
GM12878_new_igen_columns_filtered <- GM12878_new_igens_filtered[ , c("ID", "control_status", "TADs_binned")]
colnames(GM12878_new_igen_columns_filtered) <- c("Protein_ID", "Era", "TADs_binned")
GM12878_gene_new_controls_filtered <- rbind(GM12878_gene_columns_filtered, GM12878_new_igen_columns_filtered)


# graph the counts for each gene/control per era as box plots
categories <- c("Ancient", "Metazoan", "Chordate", "Mammal", "Primate", "igen_ORF", "igen_NonORF")

ggplot(GM12878_gene_new_controls_filtered, aes(x = factor(Era, categories), fill = TADs_binned)) +
  geom_bar(position = position_fill(reverse = TRUE)) +
  labs(title = "GM12878 Distribution of Eras with Enhancer Counts within Domains", 
       x = "Eras", y = "Proportion of genes with enhancer count", fill = "Enhancer number") +
  theme_minimal()

ggplot(GM12878_gene_new_controls_filtered, aes(x = factor(Era, categories), fill = TADs_binned)) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  labs(title = "GM12878 Distribution of Eras with Enhancer Counts within Domains", 
       x = "Eras", y = "Number of genes with enhancer count", fill = "Enhancer number") +
  theme_minimal()

ggplot(subset(GM12878_gene_new_controls, !loop %in% 0), aes(x = factor(Era, categories), fill = loop_binned)) +
  geom_bar(position = position_fill(reverse = TRUE)) +
  labs(title = "GM12878 Distribution of Eras with E-P Loop Counts", 
       x = "Era", y = "Proportion of genes with E-P loop count", fill = "E-P loop number") +
  theme_minimal()

ggplot(subset(GM12878_gene_new_controls, !loop %in% 0), aes(x = factor(Era, categories), fill = loop_binned)) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  labs(title = "GM12878 Distribution of Eras with E-P Loop Counts", 
       x = "Era", y = "Number of genes with E-P loop count", fill = "E-P loop number") +
  theme_minimal()


#####statistical analysis

chisq.test(table(GM12878_gene_new_controls_filtered$Era, GM12878_gene_new_controls_filtered$TADs_binned))
# results are significant, perform all pairwise comparisons to determine which specific comparisons are significant

GM12878_TAD_gene_new_control <- list()

for (i in 1:(length(categories) - 1)) {
  for (j in (i + 1):length(categories)) {
    # Filter the dataframe for the two categories
    subset_data <- subset(GM12878_gene_new_controls_filtered, Era %in% c(categories[i], categories[j]))
    
    # Perform the chi-squared test
    test_result <- chisq.test(table(subset_data$Era, subset_data$TADs_binned))
    
    # Store the result
    comparison <- paste(categories[i], "vs", categories[j])
    GM12878_TAD_gene_new_control[[comparison]] <- test_result
  }
}
results_GM12878_TAD_gene_new_control <- data.frame(
  Comparison = names(GM12878_TAD_gene_new_control),
  p_value = unlist(lapply(GM12878_TAD_gene_new_control, function(x) x$p.value))
)
results_GM12878_TAD_gene_new_control$BH_p_val_adjust <- p.adjust(results_GM12878_TAD_gene_new_control$p_value, method = "BH")

chisq.test(table(subset(GM12878_gene_new_controls, !loop %in% 0)$Era, subset(GM12878_gene_new_controls, !loop %in% 0)$loop_binned))


###### graphing by disease association

# by simple mode of inheritance
disease_eras <- c("Ancient", "Metazoan", "Chordate")
simple_class <- c("NonDisease", "Recessive", "Dominant")

ggplot(subset(GM12878_counts_filtered, Disease_included == TRUE), aes(x = factor(simple_disease, simple_class), fill = TADs_binned)) +
  geom_bar(position = position_fill(reverse = TRUE)) + 
  labs(title = "GM12878 Distribution of Genes with Enhancer Counts within Domains", 
       x = "Disease Association", y = "Proportion of genes with enhancer count", fill = "Enhancer number") +
  theme_minimal()
ggplot(subset(GM12878_counts_filtered, Disease_included == TRUE), aes(x = factor(simple_disease, simple_class), fill = TADs_binned)) +
  geom_bar(position = position_stack(reverse = TRUE)) + 
  labs(title = "GM12878 Distribution of Genes with Enhancer Counts within Domains", 
       x = "Disease Association", y = "Number of genes with enhancer count", fill = "Enhancer number") +
  theme_minimal()


#statistics
chisq.test(table(GM12878_counts_filtered$simple_disease[GM12878_counts_filtered$Disease_included == TRUE], GM12878_counts_filtered$TADs_binned[GM12878_counts_filtered$Disease_included == TRUE]))
# results are significant, perform all pairwise comparisons to determine which specific comparisons are significant

GM12878_disease_simpleTAD <- list()

for (i in 1:(length(simple_class) - 1)) {
  for (j in (i + 1):length(simple_class)) {
    # Filter the dataframe for the two categories
    subset_data <- subset(GM12878_counts_filtered, simple_disease %in% c(simple_class[i], simple_class[j]) & GM12878_counts_filtered$Disease_included == TRUE)
    
    # Perform the chi-squared test
    test_result <- chisq.test(table(subset_data$simple_disease, subset_data$TADs_binned))
    
    # Store the result
    comparison <- paste(simple_class[i], "vs", simple_class[j])
    GM12878_disease_simpleTAD[[comparison]] <- test_result
  }
}
results_GM12878_disease_simpleTAD <- data.frame(
  Comparison = names(GM12878_disease_simpleTAD),
  p_value = unlist(lapply(GM12878_disease_simpleTAD, function(x) x$p.value))
)
results_GM12878_disease_simpleTAD$BH_p_val_adjust <- p.adjust(results_GM12878_disease_simpleTAD$p_value, method = "BH")


# graphing by autosomal and sex-linked split mode of inheritance
inheritance <- c("NonDisease", "AR", "AD", "NonDisease_SexLinked", "XLR", "XLD")

ggplot(subset(GM12878_counts_filtered, Disease_included == TRUE), aes(x = factor(disease_gene_inheritance, inheritance), fill = TADs_binned)) +
  geom_bar(position = position_fill(reverse = TRUE)) + 
  labs(title = "GM12878 Distribution of Genes with Enhancer Counts within Domains", 
       x = "Disease Association", y = "Proportion of genes with enhancer count", fill = "Enhancer number") +
  theme_minimal()
ggplot(subset(GM12878_counts_filtered, Disease_included == TRUE), aes(x = factor(disease_gene_inheritance, inheritance), fill = TADs_binned)) +
  geom_bar(position = position_stack(reverse = TRUE)) + 
  labs(title = "GM12878 Distribution of Genes with Enhancer Counts within Domains", 
       x = "Disease Association", y = "Number of genes with enhancer count", fill = "Enhancer number") +
  theme_minimal()

#statistics
chisq.test(table(GM12878_counts_filtered$disease_gene_inheritance[GM12878_counts_filtered$Disease_included == TRUE], GM12878_counts_filtered$TADs_binned[GM12878_counts_filtered$Disease_included == TRUE]))
# results are significant, perform all pairwise comparisons to determine which specific comparisons are significant

GM12878_diseaseTAD <- list()

for (i in 1:(length(inheritance) - 1)) {
  for (j in (i + 1):length(inheritance)) {
    # Filter the dataframe for the two categories
    subset_data <- subset(GM12878_counts_filtered, disease_gene_inheritance %in% c(inheritance[i], inheritance[j]) & GM12878_counts_filtered$Disease_included == TRUE)
    
    # Perform the chi-squared test
    test_result <- chisq.test(table(subset_data$disease_gene_inheritance, subset_data$TADs_binned))
    
    # Store the result
    comparison <- paste(inheritance[i], "vs", inheritance[j])
    GM12878_diseaseTAD[[comparison]] <- test_result
  }
}
results_GM12878_diseaseTAD <- data.frame(
  Comparison = names(GM12878_diseaseTAD),
  p_value = unlist(lapply(GM12878_diseaseTAD, function(x) x$p.value))
)
results_GM12878_diseaseTAD$BH_p_val_adjust <- p.adjust(results_GM12878_diseaseTAD$p_value, method = "BH")
