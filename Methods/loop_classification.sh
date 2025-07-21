##Katherine Fleck
#########Loop class intersections
#pipeline shown for GM12878, also completed with K562


# intersect each anchor of loops with gene promoters, Igen promoters, enhancers, or CTCF/cohesin bound sites
for feature in Disease_gene_dataset_promoters_sorted.txt igen_controls_IDs_promoters_sorted.txt GM12878_enhancers_bygenes.bed GM12878_CTCF_Rad21_bound.bed
do
	bedtools intersect -c -a GM12878_HiCCUPS_looplist_6columns_anchor1_expanded.bed -b ${feature} > GM12878_HiCCUPS_looplist_6columns_anchor1_expanded.bed_${feature}_count_overlap.bed
	bedtools intersect -c -a GM12878_HiCCUPS_looplist_6columns_anchor2_expanded.bed -b ${feature} > GM12878_HiCCUPS_looplist_6columns_anchor2_expanded.bed_${feature}_count_overlap.bed
done

# make files with counts of promoter, enhancer, and CTCF/cohesin intersects for each anchor
paste GM12878_HiCCUPS_looplist_6columns_anchor1_expanded.bed_Disease_gene_dataset_promoters_sorted.txt_count_overlap.bed GM12878_HiCCUPS_looplist_6columns_anchor1_expanded.bed_GM12878_enhancers_bygenes.bed_count_overlap.bed GM12878_HiCCUPS_looplist_6columns_anchor1_expanded.bed_GM12878_CTCF_Rad21_bound.bed_count_overlap.bed | awk '{OFS = "\t"}{print $1, $2, $3, $4, $8, $12}' > GM12878_loop_anchor1_genes_P_E_CTCF.txt
paste GM12878_HiCCUPS_looplist_6columns_anchor2_expanded.bed_Disease_gene_dataset_promoters_sorted.txt_count_overlap.bed GM12878_HiCCUPS_looplist_6columns_anchor2_expanded.bed_GM12878_enhancers_bygenes.bed_count_overlap.bed GM12878_HiCCUPS_looplist_6columns_anchor2_expanded.bed_GM12878_CTCF_Rad21_bound.bed_count_overlap.bed | awk '{OFS = "\t"}{print $1, $2, $3, $4, $8, $12}' > GM12878_loop_anchor2_genes_P_E_CTCF.txt
paste GM12878_HiCCUPS_looplist_6columns_anchor1_expanded.bed_igen_controls_IDs_promoters_sorted.txt_count_overlap.bed GM12878_HiCCUPS_looplist_6columns_anchor1_expanded.bed_GM12878_enhancers_bygenes.bed_count_overlap.bed GM12878_HiCCUPS_looplist_6columns_anchor1_expanded.bed_GM12878_CTCF_Rad21_bound.bed_count_overlap.bed | awk '{OFS = "\t"}{print $1, $2, $3, $4, $8, $12}' > GM12878_loop_anchor1_igenP_geneE_CTCF.txt
paste GM12878_HiCCUPS_looplist_6columns_anchor2_expanded.bed_igen_controls_IDs_promoters_sorted.txt_count_overlap.bed GM12878_HiCCUPS_looplist_6columns_anchor2_expanded.bed_GM12878_enhancers_bygenes.bed_count_overlap.bed GM12878_HiCCUPS_looplist_6columns_anchor2_expanded.bed_GM12878_CTCF_Rad21_bound.bed_count_overlap.bed | awk '{OFS = "\t"}{print $1, $2, $3, $4, $8, $12}' > GM12878_loop_anchor2_igenP_geneE_CTCF.txt

# pair anchors together
paste GM12878_loop_anchor1_genes_P_E_CTCF.txt GM12878_loop_anchor2_genes_P_E_CTCF.txt > GM12878_paired_loop_anchors_genes_P_E_CTCF.txt
paste GM12878_loop_anchor1_igenP_geneE_CTCF.txt GM12878_loop_anchor2_igenP_geneE_CTCF.txt > GM12878_paired_loop_anchors_igenP_geneE_CTCF.txt

#make copies of loop anchor pair files with reversed anchor order for intersecting with anchor 2
awk '{OFS="\t"}{print $7, $8, $9, $10, $11, $12, $1, $2, $3, $4, $5, $6}' GM12878_paired_loop_anchors_genes_P_E_CTCF.txt > GM12878_paired_LA_reversed_forgenes.txt
awk '{OFS="\t"}{print $7, $8, $9, $10, $11, $12, $1, $2, $3, $4, $5, $6}' GM12878_paired_loop_anchors_igenP_geneE_CTCF.txt > GM12878_paired_LA_reversed_forigens.txt

#intersect genes or Igen controls with loop anchors keeping full information of gene/control and loop anchors
bedtools intersect -wa -wb -a Gene_CDS_promoter.bed -b GM12878_paired_loop_anchors_genes_P_E_CTCF.txt > GM12878_genes_in_anchor1.txt
bedtools intersect -wa -wb -a Gene_CDS_promoter.bed -b GM12878_paired_LA_reversed_forgenes.txt > GM12878_genes_in_anchor2.txt

bedtools intersect -wa -wb -a Igen_CDS_promoter.bed -b GM12878_paired_loop_anchors_igenP_geneE_CTCF.txt > GM12878_igens_in_anchor1.txt
bedtools intersect -wa -wb -a Igen_CDS_promoter.bed -b GM12878_paired_LA_reversed_forigens.txt > GM12878_igens_in_anchor2.txt
#these files were then merged back together to have the coordinates of anchor 1, then anchor 2, then gene info in anchor 1, then gene info in anchor 2 for each loop

#filter tables by gene era or Igen control type
for era in Ancient Metazoan Chordate Mammal Primate igenORF igenNonORF
do
	awk '{OFS="\t"}{if ($0 ~ "$era") print $0}' GM12878_loop_gene_merge.txt > GM12878_loop_gene_merge_${era}.txt
	awk '{OFS="\t"}{if ($0 ~ "$era") print $0}' GM12878_loop_igen_merge.txt > GM12878_loop_igen_merge_${era}.txt
done

#filter by loop classification subtype for loop lists for each era/Igen control: CTCF-CTCF only, cis-regulatory elements only, mixed E/P/CTCF, other-E/P/CTCF, other-other
for file in GM12878_loop_gene_merge_ancient GM12878_loop_gene_merge_metazoan GM12878_loop_gene_merge_chordate GM12878_loop_gene_merge_mammal GM12878_loop_gene_merge_primate GM12878_loop_igen_merge_igenORF GM12878_loop_igen_merge_igenNonORF
do
	awk '{OFS = "\t"}{if ($6!=0 && $12!=0 && $4==0 && $5==0 && $10==0 && $11==0) {print $0}}' ${file}.txt > ${file}_CTCF_only.txt
	awk '{OFS = "\t"}{if ($6 == 0 && $12 == 0 && ($4 != 0 || $5 != 0) && ($10 != 0 || $11 != 0)) {print $0}}' ${file}.txt > ${file}_cis_reg_only.txt
	awk '{OFS = "\t"}{if (($4 != 0 || $5 != 0 || $6 !=0) && ($10 != 0 || $11 != 0 || $12 !=0)) {print $0}}' ${file}.txt | awk 'NR==FNR {lines[$0]; next} !($0 in lines)' ${file}_CTCF_only.txt - | awk 'NR==FNR {lines[$0]; next} !($0 in lines)' ${file}_cis_reg_only.txt - > ${file}_mixed_e_p_ctcf.txt
		# this is filtering out the lines that are CTCF-CTCF only and cis-regulatory elements only to leave just those loops with mixed E/P/CTCF
	awk '{OFS = "\t"}{if ($4==0 && $5==0 && $6==0 && ($10!=0 || $11!=0 || $12!=0)) {print $0}}' ${file}.txt > ${file}_e_p_ctcf_other.txt
	awk '{OFS = "\t"}{if ($10==0 && $11==0 && $12==0 && ($4!=0 || $5!=0 || $6!=0)) {print $0}}' ${file}.txt >> ${file}_e_p_ctcf_other.txt
	awk '{OFS = "\t"}{if ($4==0 && $5==0 && $6==0 && $10==0 && $11==0 && $12==0) {print $0}}' ${file}.txt > ${file}_other_other.txt
done
