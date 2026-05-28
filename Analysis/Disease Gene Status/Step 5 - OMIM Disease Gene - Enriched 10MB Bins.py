#!/usr/bin/env python
# coding: utf-8


## Daniel Marten
## Statistics for 10MB - Correlation & Significance

import numpy as np
import pandas as pd
import scipy
from scipy.stats.mstats import zscore
import statsmodels.stats.multitest
from scipy.stats import binomtest


# Set visualizations
pd.set_option("display.max_columns", None)


# Table containing mean read count data for 6 tissue categories, with disease gene status
df5 = pd.read_csv(
    "../Marten OMIM Disease Gene - 202503134/Tables/marten_s3a_omimdisease_20250314.tsv",
    sep="\t",
    index_col="Identifier",
)

# Map values for different grupings

# Remove Sex-Linked Resolution, merge together
inheritance_map = {
    "Non-Disease": "Non-Disease",
    "Non-Disease Sex-Linked": "Non-Disease",
    "AD": "Dominant",
    "XLD": "Dominant",
    "AR": "Recessive",
    "XLR": "Recessive",
    "Mixed Autosomal": "Mixed Dominant-Recessive",
    "Mixed Autosomal Digenic": "Mixed Dominant-Recessive Digenic",
    "Mixed Sex-Linked": "Mixed Dominant-Recessive",
    "AD/AR MOI": "AD/AR MOI",
    "Digenic": "Digenic",
    "EXCLUDED-Unannotated": "EXCLUDED-Unannotated",
    "EXCLUDED-Intergenic Control": "EXCLUDED-Intergenic Control",
}

# Keep Sex-Linked split out, added resolution
publication_map = {
    "Non-Disease": "Non-Disease Autosomal",
    "Non-Disease Sex-Linked": "Non-Disease Sex-Linked",
    "AD": "AD",
    "XLD": "XLD",
    "AR": "AR",
    "XLR": "XLR",
    "Mixed Autosomal": "EXCLUDED_Mixed Autosomal",
    "Mixed Autosomal Digenic": "EXCLUDED_Mixed Autosomal Digenic",
    "Mixed Sex-Linked": "EXCLUDED_Mixed Sex-Linked",
    "AD/AR MOI": "EXCLUDED_AD/AR MOI",
    "Digenic": "EXCLUDED_Digenic",
    "EXCLUDED-Unannotated": "EXCLUDED-Unannotated",
    "EXCLUDED-Intergenic Control": "EXCLUDED-Intergenic Control",
}

# Perform this mapping
df5["disease_gene_inheritance_merged"] = [
    inheritance_map[dgi] for dgi in df5["disease_gene_inheritance"]
]

df5["disease_gene_inheritance_publication"] = [
    publication_map[dgi] for dgi in df5["disease_gene_inheritance"]
]

df5["disease_gene_inheritance"] = [
    "Non-Disease Autosomal" if xi == "Non-Disease" else xi
    for xi in df5["disease_gene_inheritance"]
]

df5["evolutionary_category_3era"] = [
    xi.replace("5-Primate", "3-Chordate").replace("4-Mammal", "3-Chordate")
    for xi in df5["evolutionary_category"]
]


# Only work with relevant genes
# Not working with unannotated genes or intergenic controls
# Further, not working with: genes with digenic MOIs, mixed AD/AR MOIs, or mixed XLD/XLR MOIs
df5_nonexcl = df5[df5.GRCh38_Expression_Plotting == True]
df5_nonexcl["Autosomal_Status"] = [
    "Sex-Linked" if xi in ["chrX", "chrY"] else "Autosomal" for xi in df5_nonexcl["chr"]
]

# Setting up some count for Chromosomes, since Integer Sorting =/= String Sorting
df5_nonexcl["smart_chr"] = [
    int(xi.replace("chr", "").replace("X", "23").replace("Y", "24"))
    for xi in df5_nonexcl["chr"]
]

df5_nonexcl["10MB_bin"] = None

# De-Index and Bin by 10MB Bins

df6_hist_deindex = df5_nonexcl.reset_index(drop=False).sort_values(
    by=["chr", "cdsstarthg38", "cdsendhg38"], ascending=True
)

for xi, yi in df6_hist_deindex.iterrows():
    floor_cdsstarthg38_10MB = yi["cdsstarthg38"] // 10_000_000
    chr_bin = yi["chr"]
    df6_hist_deindex.loc[xi, "10MB_bin"] = f"{chr_bin}:{floor_cdsstarthg38_10MB}"


# We do include End Bins (with length less than 10MB) in our analyses for consistency

# Construct and format bin data by dominant, recessive, and non-disease status
bin_data_disease = (
    df6_hist_deindex.value_counts(["10MB_bin", "disease_gene_inheritance_merged"])
    .to_frame()
    .reset_index()
)

# Creat pivot table, keyed by 10MB bin, containing disease gene data
pivot_disease = pd.pivot(
    bin_data_disease.reset_index(drop=True),
    index="10MB_bin",
    columns=["disease_gene_inheritance_merged"],
    values="count",
)
pivot_disease["Total_Count"] = pivot_disease.sum(axis=1)


# Value counts for each bin by evolutionary category
bin_data_era = (
    df6_hist_deindex.value_counts(["10MB_bin", "evolutionary_category_3era"])
    .to_frame()
    .reset_index()
)


# Create pivot table, keyed by 10MB bin
pivot_era = pd.pivot(
    bin_data_era.reset_index(drop=True),
    index="10MB_bin",
    columns=["evolutionary_category_3era"],
    values="count",
)


# Join both tables. Now, each bin has data for evolutinoary and disease information
pivot_bin = pivot_disease.join(pivot_era)

# Fill None values with 0s
pivot_bin = pivot_bin.fillna(0)


# Calculate proportions for our metrics we are interest in
pivot_bin["dominant_proportion"] = pivot_bin["Dominant"] / pivot_bin["Total_Count"]
pivot_bin["recessive_proportion"] = pivot_bin["Recessive"] / pivot_bin["Total_Count"]
pivot_bin["nondisease_proportion"] = pivot_bin["Non-Disease"] / pivot_bin["Total_Count"]
pivot_bin["disease_proportion"] = (
    1 - pivot_bin["Non-Disease"] / pivot_bin["Total_Count"]
)

pivot_bin["ancient_proportion"] = pivot_bin["1-Ancient"] / pivot_bin["Total_Count"]
pivot_bin["non_ancient_proportion"] = (
    1 - pivot_bin["1-Ancient"] / pivot_bin["Total_Count"]
)

pivot_bin["disease_sum"] = pivot_bin.Dominant + pivot_bin.Recessive
pivot_bin["non_ancient_sum"] = pivot_bin["2-Metazoan"] + pivot_bin["3-Chordate"]


# Calculate: spearman correlations between different proportions (and Total Count) by 10MB bin
# Question: are bins being older correlation with bins being more disease-rich?

genome_wide_spearman = pivot_bin[
    [
        "dominant_proportion",
        "recessive_proportion",
        "nondisease_proportion",
        "disease_proportion",
        "non_ancient_proportion",
        "ancient_proportion",
        "Total_Count",
    ]
].corr(method="spearman")

genome_wide_spearman.to_csv("corr_stats_genome_wide_spearman_correlation.tsv", sep="\t")


# Calculate the same above but by individual chromosome, and only within chromosomes

pivot_bin["contig"] = [xi.split(":")[0] for xi in pivot_bin.index]

pivot_bin["contig_int"] = [
    int(xi.split(":")[0].replace("X", "23").replace("Y", "24").replace("chr", ""))
    for xi in pivot_bin.index
]

pivot_bin["bin_int"] = [float(xi.split(":")[1]) for xi in pivot_bin.index]

pivot_bin = pivot_bin.sort_values(by=["contig_int", "bin_int"])

for contig_i in pivot_bin.contig_int.unique():
    contig_spearman = pivot_bin[pivot_bin.contig_int == contig_i][
        [
            "nondisease_proportion",
            "disease_proportion",
            "non_ancient_proportion",
            "ancient_proportion",
            "Total_Count",
        ]
    ].corr(method="spearman")

    contig_str = str(contig_i).replace("23", "X").replace("24", "Y")
    contig_spearman.to_csv(
        f"corr_stats_{contig_str}_spearman_correlation.tsv", sep="\t"
    )


def do_binomial_significance(
    binom_input_df: pd.DataFrame,
    query_col: str,
    suffix: str,
    total_str: str = "Total_Count",
) -> pd.DataFrame:
    """
    Method to calculate if a given 10MB bin has more of some query kind of gene
    than expected, given its total count of genes
    and the prevalence of that kind of query gene across the whole genome.
    Note that this method does concern proportion.

    Null Hypothesis: given the rate of the query kind of gene, and the total count of genes,
    then we would expect to see this amount of the query kind of gene by random chance
    """

    ratio = sum(binom_input_df[query_col]) / sum(binom_input_df[total_str])

    for xi, yi in binom_input_df.iterrows():
        retvals = binomtest(
            int(yi[query_col]), int(yi[total_str]), p=ratio, alternative="greater"
        ).pvalue
        binom_input_df.loc[xi, f"binomial%_{query_col}_binom_p"] = retvals

    fixed_pvals = statsmodels.stats.multitest.multipletests(
        binom_input_df[f"binomial%_{query_col}_binom_p"], method="fdr_bh"
    )[1]
    binom_input_df[f"binomial%_{query_col}_binom_p_adj"] = fixed_pvals

    return binom_input_df


def do_count_significance(
    input_df: pd.DataFrame,
    query_col: str,
    suffix: str,
) -> pd.DataFrame:
    """
    Method to calculate the z-score for the distribution of some query type of gene
    in each 10MB bin. A p-value is returned from that, and then corrected
    via Benjamini-Hochberg correction. This checks to see if a given bin has more
    of some query kind of gene than would be expected, given the mean and standard
    deviation of the distribution of it across 10MB bins. This has no bearing
    on the proportion of this given kind of gene in each bin.

    Null Hypothesis: a given bin does not have a significantly greater count of some
    kind of query gene than is expected, given the distribution of that kind of query gene.


    This method also calls do_binomial_significance() above, prints known significance values,
    and returns a DataFrame with all results.

    """

    input_df[f"{query_col}_zscore_{suffix}"] = zscore(input_df[query_col])
    input_df[f"{query_col}_pvalue_{suffix}"] = 1 - scipy.stats.norm.cdf(
        abs(input_df[f"{query_col}_zscore_{suffix}"])
    )
    input_df[f"{query_col}_pvalue_corrected_{suffix}"] = (
        statsmodels.stats.multitest.multipletests(
            input_df[f"{query_col}_pvalue_{suffix}"], method="fdr_bh"
        )[1]
    )

    input_df = do_binomial_significance(
        binom_input_df=input_df,
        query_col=query_col,
        suffix=suffix,
        total_str="Total_Count",
    )

    sig_df_tc = input_df[
        input_df[f"{query_col}_pvalue_corrected_{suffix}"] < 0.05
    ].sort_values(by=f"{query_col}_pvalue_corrected_{suffix}")

    sig_df_binom = input_df[
        input_df[f"binomial%_{query_col}_binom_p_adj"] < 0.05
    ].sort_values(by=f"binomial%_{query_col}_binom_p_adj")

    def _print_sig(input_sig_df, named_test, sig_column_title) -> None:

        print(f"Significance shape for {named_test} count: ", input_sig_df.shape)

        with pd.option_context("display.max_rows", None):

            display(
                input_sig_df[input_sig_df[sig_column_title] < 0.05].sort_values(
                    by=sig_column_title
                )
            )

    _print_sig(
        sig_df_tc, "total disease count", f"{query_col}_pvalue_corrected_{suffix}"
    )

    _print_sig(
        sig_df_binom,
        "% disease genes of total genes",
        f"binomial%_{query_col}_binom_p_adj",
    )

    return input_df


# Examine significance for total disease genes per bin
# Note: 9 are found to have significantly more disease genes than expected for the total distribution
# and 0 are found to be significantly rich for disease genes as a proportioin of the total genes in a bin
pivot_with_sig_disease = do_count_significance(
    pivot_bin, query_col="disease_sum", suffix="overall"
)


# Display and write out the findings above
pivot_with_sig_disease[
    pivot_with_sig_disease.disease_sum_pvalue_corrected_overall < 0.05
][
    [
        "Dominant",
        "Non-Disease",
        "Recessive",
        "Total_Count",
        "dominant_proportion",
        "recessive_proportion",
        "nondisease_proportion",
        "disease_sum",
        "disease_sum_zscore_overall",
        "disease_sum_pvalue_overall",
        "disease_sum_pvalue_corrected_overall",
    ]
].to_csv(
    f"corr_stats_total_disease_gene_count_significance.tsv", sep="\t"
)


# Examine significance for dominant genes per bin
# Note: 7 are found to have significantly more disease genes than expected for the total distribution
# 1 bin is found to have more dominant genes than expected, given their prevalence across the whole genome
# and the total count of genes in the given bin (ChrX:10MB-20MB)
# This may be biased by general sex-linked genes being considered sex-linked dominant
pivot_with_sig_dominant = do_count_significance(
    pivot_bin, query_col="Dominant", suffix="overall"
)


# Display and write out the findings above

pivot_with_sig_dominant[
    pivot_with_sig_dominant["binomial%_Dominant_binom_p_adj"] < 0.05
][
    [
        "Dominant",
        "Non-Disease",
        "Recessive",
        "Total_Count",
        "dominant_proportion",
        "recessive_proportion",
        "nondisease_proportion",
        "disease_sum",
        "binomial%_Dominant_binom_p",
        "binomial%_Dominant_binom_p_adj",
    ]
].to_csv(
    f"corr_stats_dominant_binomial_gene_makeup_significance.tsv", sep="\t"
)
