##Katherine Fleck

library(tidyr)

######statistical comparisons for loop class distributions in each gene era or Igen control
# pipeline shown for GM12878, also performed for K562

# import contingency table for loop type counts for each era/Igen control
GM12878_gene_igen_loops_counts <- read.table("/Users/katherinefleck/Library/Erceg_Lab/Data_files_for_3D_genome_eras_paper/loop_distribution_byage/GM12878_gene_igen_loops_counts.txt", header = TRUE, sep = "\t")
colnames(GM12878_gene_igen_loops_counts) <- c("CTCF-CTCF only", "cis-regulatory only", "mixed E/P/CTCF", "other-E/P/CTCF", "other-other")
chisq.test(GM12878_gene_igen_loops_counts)
# results are significant, perform all pairwise comparisons to determine which specific comparisons are significant

files <- list("GM12878_ancient" = GM12878_ancient, "GM12878_metazoan" = GM12878_metazoan, "GM12878_chordate"=GM12878_chordate, "GM12878_mammal"=GM12878_mammal, "GM12878_primate"=GM12878_primate, "GM12878_igenORF"=GM12878_igenORF, "GM12878_igenNonORF"=GM12878_igenNonORF)

GM12878_loop_types <- list()

for (i in 1:(length(files) - 1)) {
  for (j in (i + 1):length(files)) {
    # Filter the dataframe for the two categories
    subset_data <- rbind(files[[i]], files[[j]])
    
    # Perform the chi-squared test
    test_result <- chisq.test(subset_data)
    
    # Store the result
    comparison <- paste(names(files)[i], "vs", names(files)[j])
    GM12878_loop_types[[comparison]] <- test_result
  }
}
results_GM12878_loop_types <- data.frame(
  Comparison = names(GM12878_loop_types),
  p_value = unlist(lapply(GM12878_loop_types, function(x) x$p.value))
)
results_GM12878_loop_types$BH_p_val_adjust <- p.adjust(results_GM12878_loop_types$p_value, method = "BH")
