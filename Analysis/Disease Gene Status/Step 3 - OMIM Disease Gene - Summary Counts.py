#!/usr/bin/env python
# coding: utf-8


# Daniel Marten
# Summary Stats for Disease Genes - Relevant Counts & Means

import matplotlib.pyplot as plt
import pandas as pd
from cmapPy.pandasGEXpress.parse_gct import parse as tpm_parser


# Display options
pd.set_option("display.max_columns", None)

DATE_NAME = "0905"


# Table containing mean read count data for 6 tissue categories, with disease gene status
df5 = pd.read_csv(
    "Marten OMIM Disease Gene - 202503134/Tables/marten_s3a_omimdisease_20250314.tsv",
    sep="\t",
    index_col="Identifier",
)

# Map values for different groupings

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

# Keep sex-linked split out, added resolution
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


# Create DataFrame prepared for melting
df5_nonexcl_melt = df5_nonexcl[
    [
        "evolutionary_category_3era",
        "disease_gene_inheritance_publication",
        "BRAIN",
        "ECTO",
        "MESO",
        "ENDO",
        "TESTIS",
        "OVARY",
    ]
]


# Create a generalized method, for two input IDs, to summarize expression across tissue categories


def return_stats(input_df, input_ids):
    """Summarize expression across tissue categories for given group IDs.

    Parameters
    ----------
    input_df : pandas.DataFrame
        Wide-form DataFrame containing germ-layer columns (BRAIN, ECTO, MESO, ENDO,
        TESTIS, OVARY) and grouping columns listed in `input_ids`.
    input_ids : list[str]
        Column names to group by (e.g., ["evolutionary_category_3era",
        "disease_gene_inheritance_publication"]).

    Returns
    -------
    pandas.DataFrame
        Multi-indexed table with per-group Gene_Count, Mean, Median, and STDev
        across germ layers.
    """

    df_melted = pd.melt(
        input_df,
        id_vars=input_ids,
        value_vars=["BRAIN", "ECTO", "MESO", "ENDO", "TESTIS", "OVARY"],
        var_name="Germ_Layer",
        value_name="Count",
    )

    group_ids = input_ids + ["Germ_Layer"]

    # Count of genes fitting each category
    gene_table = (
        df_melted.groupby(group_ids).count().rename({"Count": "Gene_Count"}, axis=1)
    )

    # Mean counts
    mean_table = df_melted.groupby(group_ids).mean().rename({"Count": "Mean"}, axis=1)
    # Median counts
    median_table = (
        df_melted.groupby(group_ids).median().rename({"Count": "Median"}, axis=1)
    )

    # Standard deviation
    stdev_table = df_melted.groupby(group_ids).std().rename({"Count": "STDev"}, axis=1)

    # Table to return
    ret_table = (
        gene_table.join(mean_table)
        .join(median_table, lsuffix="_mean", rsuffix="_median")
        .join(stdev_table, rsuffix="_stdev")
    )
    ret_table = ret_table[sorted(ret_table.columns)]

    return ret_table


# Summarize results for 3 Eras with Disease Information (including sex linked details)

ret_germ_resolution = return_stats(
    df5_nonexcl_melt,
    input_ids=["evolutionary_category_3era", "disease_gene_inheritance_publication"],
)

with pd.option_context(
    "display.max_rows",
    None,
):
    display(ret_germ_resolution)

ret_germ_resolution.to_csv(f"Table_RS_Fig1c_{DATE_NAME}.tsv", sep="\t")


# Summarize results for 3 Eras with Disease Information (without sex linked details)

df5_nonexcl_nongerm = df5_nonexcl[
    [
        "evolutionary_category_3era",
        "disease_gene_inheritance_merged",
        "BRAIN",
        "ECTO",
        "MESO",
        "ENDO",
        "TESTIS",
        "OVARY",
    ]
]

ret_simple_resolution = return_stats(
    df5_nonexcl_nongerm,
    input_ids=["evolutionary_category_3era", "disease_gene_inheritance_merged"],
)
with pd.option_context(
    "display.max_rows",
    None,
):
    display(ret_simple_resolution)

ret_simple_resolution.to_csv(f"Table_RS_Fig1b_{DATE_NAME}.tsv", sep="\t")


# Summarize results and counts without any era information

df5_nonexcl_nongerm_no_eras = df5_nonexcl[
    [
        "disease_gene_inheritance_merged",
        "BRAIN",
        "ECTO",
        "MESO",
        "ENDO",
        "TESTIS",
        "OVARY",
    ]
]

df5_nonexcl_no_eras = df5_nonexcl[
    [
        "disease_gene_inheritance_publication",
        "BRAIN",
        "ECTO",
        "MESO",
        "ENDO",
        "TESTIS",
        "OVARY",
    ]
]

# Without any sex-linked information
ret_table_disease_germ_simple = return_stats(
    df5_nonexcl_nongerm_no_eras, input_ids=["disease_gene_inheritance_merged"]
)
with pd.option_context(
    "display.max_rows",
    None,
):
    display(ret_table_disease_germ_simple)
ret_table_disease_germ_simple.to_csv(f"Table_RS_Fig2_3b_{DATE_NAME}.tsv", sep="\t")

# With sex-linked information
ret_table_disease_germ = return_stats(
    df5_nonexcl_no_eras, input_ids=["disease_gene_inheritance_publication"]
)
with pd.option_context(
    "display.max_rows",
    None,
):
    display(ret_table_disease_germ)
ret_table_disease_germ.to_csv(f"Table_RS_Fig2_3c_{DATE_NAME}.tsv", sep="\t")


# Build a generalized method to tabulate the proportion of each era by disease status
# Example: the Ancient era is X.Y% dominant disease genes


def value_count_builder(input_df, era_str, disease_str):
    """Compute per-era proportions for a disease-status stratification.

    Parameters
    ----------
    input_df : pandas.DataFrame
        Long- or wide-form table containing `era_str` and `disease_str` columns.
    era_str : str
        Column name representing evolutionary era (e.g., "evolutionary_category_3era").
    disease_str : str
        Column name representing disease grouping (e.g.,
        "disease_gene_inheritance_publication").

    Returns
    -------
    pandas.DataFrame
        Indexed by (era, disease) with columns:
        - Gene_Count: counts per stratum
        - Proportion: within-era proportion for each disease status.
    """

    era_count_dict = input_df.value_counts(era_str).to_dict()

    value_count_frame = (
        input_df.value_counts([era_str, disease_str])
        .to_frame()
        .rename(columns={"count": "Gene_Count"})
    )
    value_count_frame["Proportion"] = None

    for xi, yi in value_count_frame.iterrows():
        era_ci = era_count_dict[xi[0]]
        prop_i = yi["Gene_Count"] / era_ci
        value_count_frame.loc[xi, "Proportion"] = prop_i

    return value_count_frame.sort_index()


table_4c = value_count_builder(
    df5_nonexcl, "evolutionary_category_3era", "disease_gene_inheritance_publication"
)
table_4c.to_csv(f"Table_RS_4c_{DATE_NAME}.tsv", sep="\t")
display(table_4c)


## Additional table - above, but excluding non-disease

table_4c_extra = value_count_builder(
    df5_nonexcl[
        ~df5_nonexcl.disease_gene_inheritance_publication.str.contains("Non-Disease")
    ],
    "evolutionary_category_3era",
    "disease_gene_inheritance_publication",
)
table_4c_extra.to_csv(f"Table_RS_4c_extra_{DATE_NAME}.tsv", sep="\t")
display(table_4c_extra)


table_4b = value_count_builder(
    df5_nonexcl, "evolutionary_category_3era", "disease_gene_inheritance_merged"
)
table_4b.to_csv(f"Table_RS_4b_{DATE_NAME}.tsv", sep="\t")
display(table_4b)


table_5b = value_count_builder(
    df5_nonexcl, "evolutionary_category_3era", "disease_gene_binary"
)
table_5b.to_csv(f"Table_RS_5_3eras_{DATE_NAME}.tsv", sep="\t")
display(table_5b)


table_5_5era = value_count_builder(
    df5_nonexcl, "evolutionary_category", "disease_gene_binary"
)
table_5_5era.to_csv(f"Table_RS_5_5eras_{DATE_NAME}.tsv", sep="\t")
display(table_5_5era)


table_6 = value_count_builder(
    df5_nonexcl, "Autosomal_Status", "disease_gene_inheritance_merged"
)
table_6.to_csv(f"Table_RS_6_5eras_{DATE_NAME}.tsv", sep="\t")
display(table_6)

# Reading in table - produced in step 2 - used to create the returned histograms
# Adding relevant information and annotations

df6_hist = pd.read_csv("histogram_working_table.tsv", sep="\t").sort_values(
    by="smart_chr"
)
display(df6_hist)


# From NCBI GRCh38
# https://www.ncbi.nlm.nih.gov/grc/human/data

chromosome_lengths = {
    "1": 248956422,
    "2": 242193529,
    "3": 198295559,
    "4": 190214555,
    "5": 181538259,
    "6": 170805979,
    "7": 159345973,
    "8": 145138636,
    "9": 138394717,
    "10": 133797422,
    "11": 135086622,
    "12": 133275309,
    "13": 114364328,
    "14": 107043718,
    "15": 101991189,
    "16": 90338345,
    "17": 83257441,
    "18": 80373285,
    "19": 58617616,
    "20": 64444167,
    "21": 46709983,
    "22": 50818468,
    "X": 156040895,
    "Y": 57227415,
}


# Disease Counts and Makeup per Contig - by Dominant, Recessive, and Non-Disease
chr_sum_nongerm = df6_hist[["chr"]].value_counts().sort_index().to_dict()

chr_sum_nongerm_counts = (
    df6_hist[["chr", "disease_gene_inheritance_merged"]]
    .value_counts()
    .sort_index()
    .to_frame()
)

chr_sum_nongerm_counts["chr_proportion"] = None

for xi, yi in chr_sum_nongerm_counts.iterrows():
    chr_sum = chr_sum_nongerm[(xi[0],)]

    chr_sum_nongerm_counts.loc[xi, "chr_proportion"] = yi["count"] / chr_sum
    chr_sum_nongerm_counts.loc[xi, "per_1MB"] = yi["count"] / (
        chromosome_lengths[xi[0].split(":")[0].replace("chr", "")] / 1_000_000
    )
    chr_sum_nongerm_counts.loc[xi, "contig"] = int(
        xi[0].split(":")[0].replace("chr", "").replace("X", "23").replace("Y", "24")
    )

with pd.option_context(
    "display.max_rows",
    None,
):
    display(chr_sum_nongerm_counts.sort_values(by="contig"))

chr_sum_nongerm_counts.sort_values(by="contig").to_csv(
    f"Table_RS_Fig7_contig_disease_{DATE_NAME}.tsv", sep="\t"
)


# Disease Counts and Makeup per Contig - by Disease and Non-Disease
chr_sum_nongerm_binary = df6_hist[["chr"]].value_counts().sort_index().to_dict()

chr_sum_nongerm_counts_binary = (
    df6_hist[["chr", "disease_gene_binary"]].value_counts().sort_index().to_frame()
)

chr_sum_nongerm_counts_binary["chr_proportion"] = None

for xi, yi in chr_sum_nongerm_counts_binary.iterrows():
    chr_sum_binary = chr_sum_nongerm_binary[(xi[0],)]

    chr_sum_nongerm_counts_binary.loc[xi, "chr_proportion"] = (
        yi["count"] / chr_sum_binary
    )
    chr_sum_nongerm_counts_binary.loc[xi, "per_1MB"] = yi["count"] / (
        chromosome_lengths[xi[0].split(":")[0].replace("chr", "")] / 1_000_000
    )
    chr_sum_nongerm_counts_binary.loc[xi, "contig"] = int(
        xi[0].split(":")[0].replace("chr", "").replace("X", "23").replace("Y", "24")
    )

with pd.option_context(
    "display.max_rows",
    None,
):
    display(chr_sum_nongerm_counts_binary.sort_values(by="contig"))

chr_sum_nongerm_counts_binary.sort_values(by="contig").to_csv(
    f"Table_RS_Fig7_contig_disease_binary_{DATE_NAME}.tsv", sep="\t"
)


# Era Counts and Makeup per Contig
chr_sum_nongerm_era = df6_hist[["chr"]].value_counts().sort_index().to_dict()

chr_sum_nongerm_counts_era = (
    df6_hist[["chr", "evolutionary_category_3era"]]
    .value_counts()
    .sort_index()
    .to_frame()
)

chr_sum_nongerm_counts_era["chr_proportion"] = None

for xi, yi in chr_sum_nongerm_counts_era.iterrows():
    chr_sum = chr_sum_nongerm_era[(xi[0],)]

    chr_sum_nongerm_counts_era.loc[xi, "chr_proportion"] = yi["count"] / chr_sum
    chr_sum_nongerm_counts_era.loc[xi, "contig"] = int(
        xi[0].split(":")[0].replace("chr", "").replace("X", "23").replace("Y", "24")
    )

with pd.option_context(
    "display.max_rows",
    None,
):
    display(
        chr_sum_nongerm_counts_era.sort_values(
            by=["contig", "evolutionary_category_3era"]
        )
    )

chr_sum_nongerm_counts_era.sort_values(
    by=["contig", "evolutionary_category_3era"]
).to_csv(f"Table_RS_Fig7_contig_era_{DATE_NAME}.tsv", sep="\t")


df6_hist_deindex = df6_hist.reset_index(drop=False).sort_values(
    by=["chr", "cdsstarthg38", "cdsendhg38"], ascending=True
)

for xi, yi in df6_hist_deindex.iterrows():
    floor_cdsstarthg38_1MB = yi["cdsstarthg38"] // 1_000_000
    floor_cdsstarthg38_10MB = yi["cdsstarthg38"] // 10_000_000
    floor_cdsstarthg38_50MB = yi["cdsstarthg38"] // 50_000_000
    floor_cdsstarthg38_100MB = yi["cdsstarthg38"] // 100_000_000

    chr_bin = yi["chr"]

    df6_hist_deindex.loc[xi, "1MB_bin"] = f"{chr_bin}:{floor_cdsstarthg38_1MB}"
    df6_hist_deindex.loc[xi, "10MB_bin"] = f"{chr_bin}:{floor_cdsstarthg38_10MB}"
    df6_hist_deindex.loc[xi, "50MB_bin"] = f"{chr_bin}:{floor_cdsstarthg38_50MB}"
    df6_hist_deindex.loc[xi, "100MB_bin"] = f"{chr_bin}:{floor_cdsstarthg38_100MB}"


# Disease Counts and Makeup per 10MB Bin

bin_sum_nongerm = df6_hist_deindex[["10MB_bin"]].value_counts().sort_index().to_dict()

bin_sum_nongerm_counts = (
    df6_hist_deindex[["10MB_bin", "disease_gene_inheritance_merged"]]
    .value_counts()
    .sort_index()
    .to_frame()
)

bin_sum_nongerm_counts["bin_proportion"] = None

for xi, yi in bin_sum_nongerm_counts.iterrows():
    bin_sum = bin_sum_nongerm[(xi[0],)]

    bin_sum_nongerm_counts.loc[xi, "bin_proportion"] = yi["count"] / bin_sum
    bin_sum_nongerm_counts.loc[xi, "contig"] = int(
        xi[0].split(":")[0].replace("chr", "").replace("X", "23").replace("Y", "24")
    )
    bin_sum_nongerm_counts.loc[xi, "bin_num"] = float(xi[0].split(":")[-1])

with pd.option_context(
    "display.max_rows",
    None,
):
    display(
        bin_sum_nongerm_counts.sort_values(by=["contig", "bin_num"]).drop(
            ["contig", "bin_num"], axis=1
        )
    )

bin_sum_nongerm_counts.sort_values(by=["contig", "bin_num"]).to_csv(
    f"Table_RS_Fig8_bin_disease_{DATE_NAME}.tsv", sep="\t"
)


# Era Counts and Makeup per 10MB Bin
bin_sum_nongerm_era = (
    df6_hist_deindex[["10MB_bin"]].value_counts().sort_index().to_dict()
)

bin_sum_nongerm_counts_era = (
    df6_hist_deindex[["10MB_bin", "evolutionary_category_3era"]]
    .value_counts()
    .sort_index()
    .to_frame()
)

bin_sum_nongerm_counts_era["bin_proportion"] = None


for xi, yi in bin_sum_nongerm_counts_era.iterrows():
    bin_sum = bin_sum_nongerm_era[(xi[0],)]

    bin_sum_nongerm_counts_era.loc[xi, "bin_proportion"] = yi["count"] / bin_sum
    # bin_sum_nongerm_counts.loc[xi,'bin_order'] = list_order_dict[xi[0]]
    bin_sum_nongerm_counts_era.loc[xi, "contig"] = int(
        xi[0].split(":")[0].replace("chr", "").replace("X", "23").replace("Y", "24")
    )
    bin_sum_nongerm_counts_era.loc[xi, "bin_num"] = float(xi[0].split(":")[-1])

with pd.option_context(
    "display.max_rows",
    None,
):
    display(
        bin_sum_nongerm_counts_era.sort_values(by=["contig", "bin_num"]).drop(
            ["contig", "bin_num"], axis=1
        )
    )

bin_sum_nongerm_counts_era.to_csv(f"Table_RS_Fig8_bin_era_{DATE_NAME}.tsv", sep="\t")


# Do not include End Bins (with length less than 10MB) in 10MB localization plots for consistency
# Proofed and validated elsewhere and later on.

# This contains all of the 10MB bins and coordinates we will be working with, generated from the same data
# Read in dataframe of the final bin of each chromosome, as produced in Step #2
ends = pd.read_csv("bin_endpoint_data.tsv", sep="\t")

ends_idx = ends.set_index("bin_name")

end_list = ends_idx.index.to_list()
end_list_dict = {end_i.split(":")[0]: end_i.split(":")[1] for end_i in end_list}


ends_idx = ends.set_index("bin_name")


# Calculate counts per 10MB bin - including end bins

bin_sum_dict = bin_sum_nongerm_counts.to_dict()

per_bin_rd_counts = {}

for bin_i in bin_sum_nongerm_counts.reset_index(drop=False)["10MB_bin"].unique():

    per_bin_rd_counts[bin_i] = {}

    per_bin_rd_counts[bin_i]["contig"] = bin_sum_dict["contig"][bin_i, "Non-Disease"]
    per_bin_rd_counts[bin_i]["bin_num"] = bin_sum_dict["bin_num"][bin_i, "Non-Disease"]
    per_bin_rd_counts[bin_i]["is_end"] = bin_i in list(ends.bin_name)

    if not per_bin_rd_counts[bin_i]["is_end"]:
        per_bin_rd_counts[bin_i]["length"] = 10_000_000
    else:
        per_bin_rd_counts[bin_i]["length"] = ends_idx.loc[bin_i, "length"]

    # Dominant
    try:
        per_bin_rd_counts[bin_i]["dominant"] = bin_sum_dict["count"][bin_i, "Dominant"]
    except:
        print("No dominant in: ", bin_i)
        per_bin_rd_counts[bin_i]["dominant"] = 0

    # Recessive
    try:
        per_bin_rd_counts[bin_i]["recessive"] = bin_sum_dict["count"][
            bin_i, "Recessive"
        ]
    except:
        per_bin_rd_counts[bin_i]["recessive"] = 0

    # Non-Disease
    try:
        per_bin_rd_counts[bin_i]["non-disease"] = bin_sum_dict["count"][
            bin_i, "Non-Disease"
        ]
    except:
        #         print('No dominant in: ',bin_i)
        per_bin_rd_counts[bin_i]["non-disease"] = 0

    r_i = per_bin_rd_counts[bin_i]["recessive"]
    d_i = per_bin_rd_counts[bin_i]["dominant"]
    total_disease = r_i + d_i

    nd_i = per_bin_rd_counts[bin_i]["non-disease"]


# Starting logic for writing tables for our 10MB bin data
table_13 = pd.DataFrame.from_dict(per_bin_rd_counts, orient="index").sort_values(
    by=["contig", "bin_num"]
)

# Table 13 is sorted by contig and bin_num
# Populate Table 13 with bin lengths and empty bins

table_13["is_empty"] = False

for contig_i in table_13.contig.unique():
    print(contig_i)
    chr_str = f"chr{int(contig_i)}"

    if contig_i == 23:
        chr_str = "chrX"
    elif contig_i == 24:
        chr_str = "chrY"

    end_count = int(float(end_list_dict[chr_str]))
    print(end_count)

    for bin_num in range(0, end_count + 1):

        bin_title = f"{chr_str}:{float(bin_num)}"

        if bin_title in list(table_13.index):
            print("present!")
        else:
            print(f"Check {bin_num} == {end_count}")
            is_end = bin_num == end_count
            if is_end:
                custom_length = ends_idx.loc[bin_title, "length"]
            else:
                custom_length = 10_000_000

            custom_row = [
                contig_i,
                bin_num,
                is_end,
                custom_length,
                0,
                0,
                0,
                None,
                None,
                True,
            ]

            table_13.loc[bin_title] = custom_row

            print("missing filled in!")


# Median does not need to be explicitly calculated
table_13["total_genes"] = (
    table_13.dominant + table_13.recessive + table_13["non-disease"]
)
table_13["Disease Proportion"] = (table_13.dominant + table_13.recessive) / (
    table_13.dominant + table_13.recessive + table_13["non-disease"]
)
table_13["Dominant Proportion"] = (table_13.dominant) / (
    table_13.dominant + table_13.recessive + table_13["non-disease"]
)
table_13["Recessive Proportion"] = (table_13.recessive) / (
    table_13.dominant + table_13.recessive + table_13["non-disease"]
)
table_13["Non-Disease Proportion"] = (table_13["non-disease"]) / (
    table_13.dominant + table_13.recessive + table_13["non-disease"]
)

# Version of table 13, but only with necessary columns
table_13_export = table_13.sort_values(by=["contig", "bin_num"])[
    [
        "dominant",
        "recessive",
        "non-disease",
        "total_genes",
        "RD Scaled Difference",
        "Non-Disease Scaled Difference",
        "Disease Proportion",
        "Dominant Proportion",
        "Recessive Proportion",
        "Non-Disease Proportion",
        "length",
        "is_end",
        "is_empty",
    ]
]


# Calculate counts per Contig

chr_sum_dict = chr_sum_nongerm_counts.to_dict()

per_chr_counts = {}

for chr_i in chr_sum_nongerm_counts.reset_index(drop=False)["chr"].unique():

    per_chr_counts[chr_i] = {}

    per_chr_counts[chr_i]["contig"] = chr_sum_dict["contig"][chr_i, "Non-Disease"]

    # Dominant
    try:
        per_chr_counts[chr_i]["dominant"] = chr_sum_dict["count"][chr_i, "Dominant"]
    except:
        print("No dominant in: ", chr_i)
        per_chr_counts[chr_i]["dominant"] = 0

    # Recessive
    try:
        per_chr_counts[chr_i]["recessive"] = chr_sum_dict["count"][chr_i, "Recessive"]
    except:
        per_chr_counts[chr_i]["recessive"] = 0

    # Non-Disease
    try:
        per_chr_counts[chr_i]["non-disease"] = chr_sum_dict["count"][
            chr_i, "Non-Disease"
        ]
    except:
        #         print('No dominant in: ',chr_i)
        per_chr_counts[chr_i]["non-disease"] = 0

    try:
        per_chr_counts[chr_i]["recessive"] = chr_sum_dict["count"][chr_i, "Recessive"]
    except:
        per_chr_counts[chr_i]["recessive"] = 0

    r_i = per_chr_counts[chr_i]["recessive"]
    d_i = per_chr_counts[chr_i]["dominant"]
    total_disease = r_i + d_i

    nd_i = per_chr_counts[chr_i]["non-disease"]

    try:
        per_chr_counts[chr_i]["RD Scaled Difference"] = (r_i - d_i) / (total_disease)
    except:
        per_chr_counts[chr_i]["RD Scaled Difference"] = None

    try:
        per_chr_counts[chr_i]["Non-Disease Scaled Difference"] = (
            nd_i - total_disease
        ) / (nd_i + total_disease)
    except:
        per_chr_counts[chr_i]["Non-Disease Scaled Difference"] = None

    # print('Not applicable for: ',chr_i)


table_13_chr = pd.DataFrame.from_dict(per_chr_counts, orient="index").sort_values(
    by=["contig"]
)
table_13_chr.to_csv(f"Table_RS_Fig13_chromosome_{DATE_NAME}.tsv", sep="\t")


# Total Chromosomal Density

# Disease Counts and Makeup per Contig
chr_sum_nongerm = df6_hist[["chr"]].value_counts().sort_index().to_dict()

chr_sum_nongerm_counts = df6_hist[["chr"]].value_counts().sort_index().to_frame()

chr_sum_nongerm_counts["chr_proportion"] = None

for xi, yi in chr_sum_nongerm_counts.iterrows():
    chr_sum = chr_sum_nongerm[(xi[0],)]

    chr_sum_nongerm_counts.loc[xi, "per_1MB"] = yi["count"] / (
        chromosome_lengths[xi[0].split(":")[0].replace("chr", "")] / 1_000_000
    )
    chr_sum_nongerm_counts.loc[xi, "contig"] = int(
        xi[0].split(":")[0].replace("chr", "").replace("X", "23").replace("Y", "24")
    )

with pd.option_context(
    "display.max_rows",
    None,
):
    display(chr_sum_nongerm_counts.sort_values(by="contig"))

chr_sum_nongerm_counts.sort_values(by="contig").to_csv(
    f"Table_RS_total_chromosomal_{DATE_NAME}.tsv", sep="\t"
)


# additional table, by sex-linked status as well
# Summarize results for 3 Eras with Disease Information (without sex linked details)

df5_nonexcl_germy = df5_nonexcl[
    [
        "disease_gene_inheritance_publication",
        "BRAIN",
        "ECTO",
        "MESO",
        "ENDO",
        "TESTIS",
        "OVARY",
    ]
]

ret_autosomal_resolution = return_stats(
    df5_nonexcl_germy, input_ids=["disease_gene_inheritance_publication"]
)
with pd.option_context(
    "display.max_rows",
    None,
):
    display(ret_autosomal_resolution)

ret_autosomal_resolution.to_csv(
    f"Table_RS_FigS_mc_autosomal_dis_counts_{DATE_NAME}.tsv", sep="\t"
)


# Preparing to separately write 10MB bin data out to an Excel table with separate sheets,
# one with empty & end bins, and one without

table_13_export_plotting = table_13_export[
    (table_13_export.is_empty == False) & (table_13_export.is_end == False)
]


with pd.ExcelWriter("marten_10MB_stats_20250905.xlsx") as writer:
    table_13_export_plotting.to_excel(writer, sheet_name="S5E_plotting_10MBstats")
    table_13_export.to_excel(writer, sheet_name="S5E_empty_and_end_10MBstats")
