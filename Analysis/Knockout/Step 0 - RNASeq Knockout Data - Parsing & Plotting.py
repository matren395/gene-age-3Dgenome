#!/usr/bin/env python
# coding: utf-8


# Daniel Marten 2025-05-27
# Import & Parse Cell-Level Knockout RNASeq Data

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from IPython.display import display  # type: ignore

pd.set_option("display.max_columns", None)


# For saving & writing
date_name = "20250623"
do_plot = True


# List of files
rna_seq_base = "~/Downloads/Hap1_RNASeq/"
file_names = os.popen(f"ls {rna_seq_base}*.txt.gz").read().split("\n")[:-1]


# Read all files and place into a list of DataFrames

df_list = []

for xi in file_names:
    name_splits = xi.split("_")
    df_temp = (
        pd.read_csv(f"{xi}", sep="\t", header=None)
        .rename({0: "ENSG", 1: "Count"}, axis=1)
        .set_index("ENSG")
    )
    df_temp["Replicate"] = name_splits[-2]
    df_temp["Condition"] = "_".join(name_splits[4:6])
    df_temp["Triplicate"] = name_splits[6]

    df_list.append(df_temp)


# Merge into one large DataFrame
df_mega = (
    pd.concat(df_list)
    .sort_values(by=["Replicate", "Condition", "Triplicate"])
    .sort_index()
)

df_mega["Condition_Replicate"] = df_mega.Condition + "_" + df_mega.Replicate


# For each condition, compute statistics across biological replicates
# Ignore the technical lane suffix (e.g., _L003 vs _L006)


ensg_counts = (
    df_mega.reset_index(drop=False)
    .groupby(["ENSG", "Condition"], as_index=False)["Count"]
    .describe()
    .set_index("ENSG")
)  # .to_frame()
ensg_counts["coeff_of_variation"] = ensg_counts["std"] / ensg_counts["mean"]


## Set up Translation Dictionaries

translationDict = pd.read_csv(
    r"/Users/marten/Downloads/biomart_export_grch38_ens93_jul2018_20231121.txt",
    sep="\t",
)  # Using GRCh38 mapping
# ENSGs remain stable between versions (e.g., ENSG000567 does not become ENSG000771).
# What changes is the mapping to transcripts and proteins.


translationDict["Ensembl_gene_length"] = (
    translationDict["Gene end (bp)"] - translationDict["Gene start (bp)"]
)

# Note: consider dropping rows with missing values if needed.
translationDict.rename(
    columns={
        "Gene stable ID": "ENSG",
        "Transcript stable ID": "ENST",
        "Protein stable ID": "ENSP",
    },
    inplace=True,
)
# display(translationDict)

translationDict = translationDict[["ENSG", "ENST", "ENSP", "Ensembl_gene_length"]]

# Make different versions of the dictionary indexed on ENSG, ENST, and ENSP.
transEnsg = translationDict.copy()
transEnsg.set_index("ENSG", inplace=True)

transEnst = translationDict.copy()
transEnst.set_index("ENST", inplace=True)

transEnsp = translationDict.copy()
transEnsp.set_index("ENSP", inplace=True)

# 120,755 rows; many are alternative transcripts of the same gene/protein.


# Import mean read count data used throughout the paper, containing OMIM disease status
s3a = pd.read_csv(
    "~/Serious_Work/Batch_Scripts/unannotated-genes-with-all/Marten OMIM Disease Gene - 202503134/Tables/marten_s3a_omimdisease_20250314.tsv",
    sep="\t",
)


# Filter S3A to entries with known ENSGs
s3a_ensg = s3a[~s3a.ENSG.isna()]  # 19,283 OK


unique_conditions = list(df_mega.Condition.unique())


# Filter table with OMIM gene-disease associations and expression data to those with recorded RNA-Seq counts from knockout experiments
# Done by matching gene identifiers from knockout data to those in our original Table S1A
s3a_ensg_choice = s3a_ensg[s3a_ensg.ENSG.isin(ensg_counts.index)]


# Index for future joining and filtering
ensg_multiindex = ensg_counts.reset_index().set_index(["ENSG", "Condition"])


# Pull mean and standard deviation from prior data
for annotation in unique_conditions:
    for measurement in ["mean", "std"]:
        col_name = f"{annotation}_{measurement}"
        print(col_name)
        s3a_ensg_choice[col_name] = [
            ensg_multiindex.loc[(xi, annotation), measurement]
            for xi in s3a_ensg_choice.ENSG
        ]


# Handling for multiple proteins per ENSG
s3a_ensg_choice_expression = s3a_ensg_choice[s3a_ensg_choice.GRCh38_Expression_Plotting]
s3a_ensg_choice_sort = s3a_ensg_choice_expression.sort_values(
    ["ps", "plength"], ascending=[True, False]
)
s3a_ensg_choice_sort_ensg = s3a_ensg_choice_sort[
    ~s3a_ensg_choice_sort.ENSG.duplicated(keep="first")
]


s3a_ensg_choice_sort_ensg["evolutionary_category_3era"] = [
    xi.replace("5-Primate", "3-Chordate").replace("4-Mammal", "3-Chordate")
    for xi in s3a_ensg_choice_sort_ensg["evolutionary_category"]
]


constants = ["ENSG", "evolutionary_category_3era"]
condition_means = [f"{condition_i}_mean" for condition_i in unique_conditions]
condition_stdevs = [f"{condition_i}_std" for condition_i in unique_conditions]
columns_list = constants + condition_means + condition_stdevs


# "Melt" knockout data to fit prior formatting

s3a_ensg_choice_mini = s3a_ensg_choice_sort_ensg[columns_list]
s3a_ensg_choice_minimelt = s3a_ensg_choice_mini.melt(
    id_vars=constants, value_vars=condition_means, var_name="Count"
)

s3a_ensg_choice_minimelt["Count"] = [
    xi.replace("_mean", "") for xi in s3a_ensg_choice_minimelt["Count"]
]


# Palette for plotting
knockout_palette = {
    "DKO_KO": "orange",
    "Wapl_KO": "cyan",
    "Hap1_WT": "red",
    "SCC4_KO": "lime",
}


hue_order = ["Hap1_WT", "SCC4_KO", "Wapl_KO", "DKO_KO"]

# Logic for log10 plotting with zero values (replace zeros with the smallest non-zero value)
smallest_nonzero = s3a_ensg_choice_minimelt[
    s3a_ensg_choice_minimelt.value > 0
].value.min()
print("Smallest non-zero count as: ", smallest_nonzero)
s3a_ensg_choice_minimelt["value_nonzero"] = s3a_ensg_choice_minimelt.value.replace(
    0.0, smallest_nonzero
)
s3a_ensg_choice_minimelt["value_log_nonzero"] = np.log10(
    s3a_ensg_choice_minimelt["value_nonzero"]
)


# Generate plots for:
# Evolutionary category by 3 era
# # - With and without log10 scale

# Plots for relevant setups

for evo_config in ["evolutionary_category_3era"]:
    for log_status in ["value", "value_log_nonzero"]:

        plt.figure(figsize=(12, 8), dpi=80)
        sns.boxplot(
            data=s3a_ensg_choice_minimelt.sort_values(by=evo_config),
            x=evo_config,
            y=log_status,
            hue="Count",
            showfliers=False,
            linewidth=3,
            palette=knockout_palette,
            hue_order=hue_order,
        )
        plt.xlabel("Evolutionary Era", size=32)
        plt.ylabel(
            f"RNASeq Counts{' Log10' if log_status=='value_log_nonzero' else ''}",
            size=32,
        )
        xOnes, xTwos = plt.xticks()

        xfont = 26

        xTwos = ["Ancient", "Metazoan", "Chordate"]

        era_text = "3era"
        era_num = "3"

        # print(plt.ylim())

        plt.ylim((-133.80833333333334, 2809.975))
        if log_status == "value_log_nonzero":
            plt.ylim((-1.0614816416482264, 5.171786966172596))

        plt.xticks(xOnes, xTwos, size=xfont)
        plt.yticks(size=27)
        plt.legend(
            title="Cell Line Status",
            title_fontsize=22,
            fontsize=16,
            bbox_to_anchor=(1.00, 1),
        )
        plt.title(
            f"Knockout Treatment Expression{' Log10' if log_status=='value_log_nonzero' else ''} {era_num} Era",
            size=32,
            pad=15,
        )

        image_name = f"knockout-expression-{'Log10' if log_status=='value_log_nonzero' else ''}{era_num}eras-{date_name}"
        if do_plot:
            plt.savefig(f"{image_name}.png", bbox_inches="tight")
            plt.savefig(f"{image_name}.pdf", bbox_inches="tight")

        plt.show()


for setup in [
    ["evolutionary_category_3era"],
]:
    temp_table = s3a_ensg_choice_sort_ensg.groupby(setup)[condition_means].describe()
    save_str = "_".join(setup)
    print("to: ", save_str)
    display(temp_table)
    temp_table.to_csv(f"marten_summary_counts_{save_str}_{date_name}.tsv", sep="\t")


# Desired: marten_summary_counts_evolutionary_category_3era
