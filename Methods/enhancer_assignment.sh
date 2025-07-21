##Katherine Fleck
#########Enhancer count per gene analysis

###pipeline was also completed for genes in K562 and Igen controls in both GM12878 and K562
###assigning enhancers to gene/Igen promoters within the same minimal TAD

#intersect promoters with domains to obtain the coordinates of domains that contain gene promoters
bedtools intersect -wa -wb -a GM12878_TADs_3column_sorted_lengths.bed -b Disease_gene_dataset_promoters_sorted.txt > GM12878_TADs_w_gene_promoter.bed

#filter to only keep the smallest domain with promoters
for ID in $(cat disease_protein_IDs.txt)
do
	grep -w "$ID" GM12878_TADs_w_gene_promoter.bed | sort -k4,4n | head -n 1 >> GM12878_TADs_w_gene_promoter_filtered.bed
done

#intersect enhancers with domain to obtain the coordinates of domains that contain enhancers
bedtools intersect -wa -wb -a GM12878_TADs_3column_sorted_lengths.bed -b GM12878_enhancers_bygenes.bed > GM12878_TADs_w_gene_enhancer.bed

#make files of just the enhancer ID #'s that are in domains
awk '{OFS="\t"} {print $8}' GM12878_TADs_w_gene_enhancer.bed | sort | uniq > GM12878_TADs_w_gene_enhancer_IDs.txt

#filter to only keep the smallest domain with enhancers
for ID in $(cat GM12878_TADs_w_gene_enhancer_IDs.txt)
do
	grep -w "$ID" GM12878_TADs_w_gene_enhancer.bed | sort -k4,4n | head -n 1 >> GM12878_TADs_w_gene_enhancer_filtered.bed
done

#reorganize files to run bedtools closest between enhancers and promoters that are in these domains
awk '{OFS="\t"} {print $5, $6, $7, $8, $9, $10, $11, $1, $2, $3, $4}' GM12878_TADs_w_gene_promoter_filtered.bed | sort -k1,1 -k2,2n > GM12878_TADs_w_gene_promoter_filtered_reordered.bed
awk '{OFS="\t"} {print $5, $6, $7, $8, $1, $2, $3, $4}' GM12878_TADs_w_gene_enhancer_filtered.bed | sort -k1,1 -k2,2n > GM12878_TADs_w_gene_enhancer_filtered_reordered.bed

#find closest promoter for each enhancer
bedtools closest -a GM12878_TADs_w_gene_enhancer_filtered_reordered.bed -b GM12878_TADs_w_gene_promoter_filtered_reordered.bed > GM12878_enhancer_closest_gene_promoter.bed

#filter to keep only rows where enhancer and promoter are in the same domain, matching the chr, start, and stop of domains
awk '{ if ($5 == $16 && $6 == $17 && $7 == $18) print $0}' GM12878_enhancer_closest_gene_promoter.bed > GM12878_enhancer_closest_gene_promoter_filter_byTAD.bed

#get count of number of times each protein ID (or Igen ID) appears as the count of enhancers assigned to that gene
for ID in $(cat disease_protein_IDs.txt)
do
	echo "$ID" >> temp.txt
	grep -cw "$ID" GM12878_enhancer_closest_gene_promoter_filter_byTAD.bed >> counts_temp.txt
done
paste temp.txt counts_temp.txt > GM12878_gene_enhancer_counts_byTAD.txt

#make list of promoters that are not in domains to exclude from graphing and statistical analysis
bedtools intersect -v -a Disease_gene_dataset_promoters_sorted.txt -b GM12878_TADs_3column_sorted_lengths.bed > GM12878_gene_promoters_notin_TADs.bed



###assigning E-P loops to gene promoters by looping interactions

#filter paired loop anchors with counts of promoters, enhancers, and CTCF/cohesin binding (see "loop_classification.sh") for loops with enhancers at one anchor and promoters at the other anchor
awk '{OFS = "\t"}{if ($4!=0 && $11!=0) {print $0}}' GM12878_paired_loop_anchors_genes_P_E_CTCF.txt > GM12878_E_P_paired_loop_anchors_bygenes.txt
awk '{OFS = "\t"}{if ($5!=0 && $10!=0) {print $0}}' GM12878_paired_loop_anchors_genes_P_E_CTCF.txt >> GM12878_E_P_paired_loop_anchors_bygenes.txt

#concatenate the anchors that contain promoters from the E-P filtered loop files together
awk '{OFS="\t"}{if ($5!=0 && $10!=0) {print $7, $8, $9} else if ($11!=0 && $4!=0) {print $1, $2, $3}}' GM12878_E_P_paired_loop_anchors_bygenes.txt > GM12878_promoter_anchors_bygenes.txt

#intersect these anchor files with promoters to get gene/Igen control IDs
bedtools intersect -wa -a Disease_gene_dataset_promoters_sorted.txt -b GM12878_promoter_anchors_bygenes.txt > Disease_gene_dataset_promoters_sorted.txt_GM12878_promoter_anchors_bygenes_inclusive.txt_Overlap_Sorted.bed

#count number of appearances of each protein/Igen control ID in these output files as the number of E-P loops assigned to each gene/control
for ID in $(cat disease_protein_IDs.txt)
do
	grep -cw "$ID" Disease_gene_dataset_promoters_sorted.txt_GM12878_promoter_anchors_bygenes.txt_Overlap_Sorted.bed >> GM12878_gene_promoters_looped_to_enhancers_proteinID_counts.bed_temp.txt
done
paste disease_protein_IDs.txt GM12878_gene_promoters_looped_to_enhancers_proteinID_counts.bed_temp.txt > GM12878_enhancer_looped_proteinID_counts.bed
