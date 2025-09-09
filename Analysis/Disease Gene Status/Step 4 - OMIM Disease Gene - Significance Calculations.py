#!/usr/bin/env python
# coding: utf-8

"""
Disease Gene Status Significance Analysis Script.

This script performs comprehensive statistical significance testing for disease gene
expression patterns across different evolutionary eras, tissue categories, and inheritance
patterns. It generates significance tables comparing gene expression between different
disease status groups and evolutionary categories using Mann-Whitney U tests with
multiple hypothesis correction.

The script processes OMIM disease gene data and performs pairwise comparisons across:
- Evolutionary eras (3-era classification)
- Tissue categories (Brain, Ectoderm, Mesoderm, Endoderm, Ovary, Testis)
- Disease inheritance patterns (Dominant, Recessive, Non-Disease, etc.)
- Sex-linked vs autosomal gene status

Outputs:
    - `marten_significance_joint_autosomal_all_sep25.tsv`: Complete significance results
    - `Table_Significance_ErrorReWrite.xlsx`: Excel file with organized results by category
"""

# Daniel Marten
# Significance Stats for Disease Genes

import pandas as pd
import scipy
import statsmodels

# Set visualizations
pd.set_option("display.max_columns", None)


def sig_table(
    eras_input,
    tissue_input,
    incoming_df,
    era_col,
    tissue_col,
    verbose=False,
    do_mean=False,
    novel_limitation=True,
):
    """
    Create table of significance calculations using Mann-Whitney U tests.

    Performs pairwise comparisons between different era and tissue combinations
    using Mann-Whitney U tests. Handles data grouping, statistical testing,
    and result compilation into a comprehensive DataFrame.

    Args:
        eras_input (list): List of era categories to compare
        tissue_input (list): List of tissue categories to compare
        incoming_df (pd.DataFrame): Input DataFrame with expression data
        era_col (str): Column name for era categories
        tissue_col (str): Column name for tissue categories
        verbose (bool): Whether to print progress information
        do_mean (bool): Whether to use mean values for grouping
        novel_limitation (bool): Whether to apply minimum count limitations

    Returns:
        pd.DataFrame: Results table with statistical comparisons including
                     p-values, means, medians, standard deviations, and counts
    """

    # Create empty lists that will become the dataframe
    records = 0

    left = []
    right = []
    floatvals = []
    pvals = []

    incoming_df = (
        incoming_df.reset_index(drop=False)
        .set_index([era_col, tissue_col])
        .sort_index()
    )

    for x1_idx in range(len(eras_input)):
        for x2 in eras_input[x1_idx:]:
            x1 = eras_input[x1_idx]
            if verbose:
                print(f"Comparing from: {x1} & {x2}")

            for y1_idx in range(len(tissue_input)):
                for y2 in tissue_input[y1_idx:]:

                    records += 1

                    y1 = tissue_input[y1_idx]
                    if verbose:
                        print(f"comparing ({x1},{y1}) to ({x2},{y2})")
                        print("Count at: ", records)

                    if (x1 == x2) and (y1 == y2):
                        if verbose:
                            print("Duplicate - pass")

                    else:

                        xinput = incoming_df.loc[(x1, y1)]

                        yinput = incoming_df.loc[(x2, y2)]

                        if do_mean == True:
                            if "Name" in xinput.columns:
                                xinput = xinput.groupby("Name")["Mean(Count)"].mean()
                                yinput = yinput.groupby("Name")["Mean(Count)"].mean()
                            else:
                                xinput = xinput.groupby(xinput.index)[
                                    "Mean(Count)"
                                ].mean()
                                yinput = yinput.groupby(yinput.index)[
                                    "Mean(Count)"
                                ].mean()
                        else:
                            xinput = xinput["Mean(Count)"]
                            yinput = yinput["Mean(Count)"]

                        x_mean, x_median, x_std, x_stderr, x_count = (
                            xinput.mean(),
                            xinput.median(),
                            xinput.std(),
                            xinput.sem(),
                            len(xinput),
                        )

                        y_mean, y_median, y_std, y_stderr, y_count = (
                            yinput.mean(),
                            yinput.median(),
                            yinput.std(),
                            yinput.sem(),
                            len(yinput),
                        )

                        left.append(
                            (x1, y1, x_mean, x_median, x_std, x_stderr, x_count)
                        )
                        right.append(
                            (x2, y2, y_mean, y_median, y_std, y_stderr, y_count)
                        )

                        if novel_limitation and any(
                            [xi < 20 for xi in [x_count, y_count]]
                        ):
                            # Skip statistical testing for low count groups
                            statistics_obj_list = [None] * 2

                        else:

                            statistics_obj_list = scipy.stats.mannwhitneyu(
                                x=xinput, y=yinput, nan_policy="omit"
                            )

                        floatvals.append(statistics_obj_list[0])
                        pvals.append(statistics_obj_list[1])

    dfResults = pd.DataFrame()

    dfResults["Left Status"] = [tup[0] for tup in left]
    dfResults["Left Tissue_Group"] = [tup[1] for tup in left]
    dfResults["Left Mean"] = [tup[2] for tup in left]
    dfResults["Left Median"] = [tup[3] for tup in left]
    dfResults["Left STDev"] = [tup[4] for tup in left]
    dfResults["Left STDerr"] = [tup[5] for tup in left]
    dfResults["Left Count"] = [tup[-1] for tup in left]

    dfResults["Right Status"] = [tup[0] for tup in right]
    dfResults["Right Tissue_Group"] = [tup[1] for tup in right]
    dfResults["Right Mean"] = [tup[2] for tup in right]
    dfResults["Right Median"] = [tup[3] for tup in right]
    dfResults["Right STDev"] = [tup[4] for tup in right]
    dfResults["Right STDerr"] = [tup[5] for tup in right]
    dfResults["Right Count"] = [tup[-1] for tup in right]

    dfResults["Floats"] = floatvals
    dfResults["pvals"] = pvals
    dfResults.sort_values(by="pvals", inplace=True)
    dfResults.reset_index(drop=True, inplace=True)

    return dfResults


# Load table containing mean read count data for 6 tissue categories, with disease gene status
df5 = pd.read_csv(
    "Marten OMIM Disease Gene - 202503134/Tables/marten_s3a_omimdisease_20250314.tsv",
    sep="\t",
    index_col="Identifier",
)

# Map values for different groupings

# Remove Sex-Linked Resolution
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

# Keep Sex-Linked split out
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

# Actually perform the mappings
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


# Filter to only include relevant genes for analysis
df5_nonexcl = df5[df5.GRCh38_Expression_Plotting == True]
df5_nonexcl["Autosomal_Status"] = [
    "Sex-Linked" if xi in ["chrX", "chrY"] else "Autosomal" for xi in df5_nonexcl["chr"]
]


# Select relevant columns for analysis
df_meancount = df5_nonexcl[
    [
        "BRAIN",
        "ECTO",
        "MESO",
        "ENDO",
        "OVARY",
        "TESTIS",
        "evolutionary_category_3era",
        "disease_gene_inheritance_merged",
    ]
]


# Select relevant columns for analysis
df_meancount_sexsplit = df5_nonexcl[
    [
        "BRAIN",
        "ECTO",
        "MESO",
        "ENDO",
        "OVARY",
        "TESTIS",
        "evolutionary_category_3era",
        "disease_gene_inheritance_publication",
    ]
]


# Create melted dataframe for analysis
melt_df = pd.melt(
    frame=df_meancount,
    id_vars=["evolutionary_category_3era", "disease_gene_inheritance_merged"],
    var_name="Tissue",
    value_name="Mean(Count)",
    ignore_index=False,
)
melt_df["joint_disease_era"] = [
    f"{xi[0]}_{xi[1]}"
    for xi in zip(
        melt_df.evolutionary_category_3era, melt_df.disease_gene_inheritance_merged
    )
]


# Create melted dataframe for analysis
melt_df_sexsplit = pd.melt(
    frame=df_meancount_sexsplit,
    id_vars=["evolutionary_category_3era", "disease_gene_inheritance_publication"],
    var_name="Tissue",
    value_name="Mean(Count)",
    ignore_index=False,
)
melt_df_sexsplit["joint_disease_era"] = [
    f"{xi[0]}_{xi[1]}"
    for xi in zip(
        melt_df_sexsplit.evolutionary_category_3era,
        melt_df_sexsplit.disease_gene_inheritance_publication,
    )
]


# Define groups for statistical comparisons
# Generate joint groups to enable cross-category comparisons (e.g., Ancient Dominant vs Chordate Non-Disease)

eragroups_1 = sorted(melt_df.evolutionary_category_3era.unique())
tissueCatKeys = sorted(melt_df.Tissue.unique())
joint_groups = sorted(melt_df.joint_disease_era.unique())
disease_groups = sorted(melt_df.disease_gene_inheritance_merged.unique())


disease_groups_sexsplit = sorted(
    melt_df_sexsplit.disease_gene_inheritance_publication.unique()
)
joint_groups_sexsplit = sorted(melt_df_sexsplit.joint_disease_era.unique())


# Create tissue category mapping dictionary
tissue_dict = {}

for xi in tissueCatKeys:
    tissue_dict[xi] = [xi]


# Generate significance table comparing disease status and evolutionary era
df_disease_era = sig_table(
    eras_input=joint_groups,
    tissue_input=tissueCatKeys,
    incoming_df=melt_df,
    era_col="joint_disease_era",
    tissue_col="Tissue",
    do_mean=False,
    verbose=True,
)


df_disease_era_sexsplit = sig_table(
    eras_input=joint_groups_sexsplit,
    tissue_input=tissueCatKeys,
    incoming_df=melt_df_sexsplit,
    era_col="joint_disease_era",
    tissue_col="Tissue",
    do_mean=False,
    verbose=True,
)


# Generate significance table comparing only disease status (without era)
df_disease = sig_table(
    eras_input=disease_groups,
    tissue_input=tissueCatKeys,
    incoming_df=melt_df,
    era_col="disease_gene_inheritance_merged",
    tissue_col="Tissue",
    do_mean=False,
    verbose=True,
)


df_disease_sexsplit = sig_table(
    eras_input=disease_groups_sexsplit,
    tissue_input=tissueCatKeys,
    incoming_df=melt_df_sexsplit,
    era_col="disease_gene_inheritance_publication",
    tissue_col="Tissue",
    do_mean=False,
    verbose=True,
)


# Establish metadata for joining different significance tables
df_disease["includes_era"] = False
df_disease_era["includes_era"] = True

df_disease["autosomal_split"] = False
df_disease_era["autosomal_split"] = False


# _sexsplit
df_disease_sexsplit["includes_era"] = False
df_disease_era_sexsplit["includes_era"] = True

df_disease_sexsplit["autosomal_split"] = True
df_disease_era_sexsplit["autosomal_split"] = True


# Join tables before any statistical corrections, so that we have the proper number of observations
df_disease_joint = pd.concat(
    [df_disease, df_disease_era, df_disease_sexsplit, df_disease_era_sexsplit]
)


# Filter out failed statistical tests (for comparisons, not summary stats)
df_disease_joint_failures = df_disease_joint[df_disease_joint.pvals.isna()]

df_disease_joint = df_disease_joint[~df_disease_joint.pvals.isna()]


fixedpvals_disease_joint = statsmodels.stats.multitest.multipletests(
    df_disease_joint["pvals"], method="fdr_bh"
)[1]

df_disease_joint["Adjusted_Pvals_BH"] = fixedpvals_disease_joint
df_disease_joint["Adjusted_Pvals < 0.05"] = [
    fp < 0.05 for fp in fixedpvals_disease_joint
]

df_disease_joint_failures["Adjusted_Pvals_BH"] = "Not Applicable"
df_disease_joint_failures["Adjusted_Pvals < 0.05"] = "Not Applicable"
df_disease_joint_failures["Floats"] = "Not Applicable"
df_disease_joint_failures["pvals"] = "Not Applicable"
df_disease_joint_failures["pvals"] = "Not Applicable"


df_disease_joint = pd.concat([df_disease_joint, df_disease_joint_failures])


def comparison_key(df):
    """
    Create a unique comparison key for each statistical comparison.

    Generates a string identifier that uniquely identifies each pairwise
    comparison by combining the left and right status and tissue group
    information.

    Args:
        df (pd.DataFrame): DataFrame containing comparison results

    Returns:
        str: Unique comparison identifier in format:
             "LeftStatus_LeftTissue_RightStatus_RightTissue"
    """
    return (
        df["Left Status"]
        + "_"
        + df["Left Tissue_Group"]
        + "_"
        + df["Right Status"]
        + "_"
        + df["Right Tissue_Group"]
    )


df_disease_joint["comparison_key"] = comparison_key(df_disease_joint)
df_disease_joint.set_index("comparison_key", inplace=True)


df_disease_joint.to_csv("marten_significance_joint_autosomal_all_sep25.tsv", sep="\t")


df_disease_joint_era = df_disease_joint[
    (df_disease_joint.includes_era == True)
    & (df_disease_joint.autosomal_split == False)
]

df_disease_joint_simple = df_disease_joint[
    (df_disease_joint.includes_era == False)
    & (df_disease_joint.autosomal_split == False)
]

df_disease_joint_era_sexsplit = df_disease_joint[
    (df_disease_joint.includes_era == True) & (df_disease_joint.autosomal_split == True)
]

df_disease_joint_simple_sexsplit = df_disease_joint[
    (df_disease_joint.includes_era == False)
    & (df_disease_joint.autosomal_split == True)
]


# Read knockout significance data and export combined results

knockout = pd.read_csv("~/ug-gc/marten_knockout_era_significance.tsv", sep="\t")

with pd.ExcelWriter("Table_Significance_ErrorReWrite.xlsx") as writer:
    df_disease_joint_simple.to_excel(writer, sheet_name="Disease")
    df_disease_joint_era.to_excel(writer, sheet_name="Era_Disease")
    df_disease_joint_simple_sexsplit.to_excel(writer, sheet_name="SXLDisease")
    df_disease_joint_era_sexsplit.to_excel(writer, sheet_name="Era_SXLDisease")
    knockout.to_excel(writer, sheet_name="Knockout")
