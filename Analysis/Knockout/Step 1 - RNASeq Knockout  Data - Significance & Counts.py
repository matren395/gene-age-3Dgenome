#!/usr/bin/env python
# coding: utf-8


# Significance Stats for Knockout Analysis

import pandas as pd
import scipy
import statsmodels
from IPython.display import display  # type: ignore


# Display options
pd.set_option("display.max_columns", None)


# Simplified: optionally compute group means.


def sig_table(
    eras_input,
    tissue_input,
    incoming_df,
    era_col,
    verbose=False,
    do_mean=False,
):

    # Create containers that will become the DataFrame.

    dfResults = pd.DataFrame()

    left = []
    right = []
    floatvals = []
    pvals = []
    # removed unused accumulators

    for x1_idx in range(len(eras_input)):
        for x2 in eras_input[x1_idx:]:
            x1 = eras_input[x1_idx]
            # x1 and x2 iterate over all combinations of era labels.
            if verbose:
                print(f"Comparing from: {x1} & {x2}")

            for y1_idx in range(len(tissue_input)):
                for y2 in tissue_input[y1_idx:]:
                    y1 = tissue_input[y1_idx]
                    # print(x1,x2)
                    # print(y1,y2)
                    # y1 and y2 iterate over all combinations of knockout-group labels.
                    if verbose:
                        print(f"comparing ({x1},{y1}) to ({x2},{y2})")
                    if (x1 == x2) and (y1 == y2):
                        if verbose:
                            print("Duplicate - pass")
                    else:

                        xinput = incoming_df[
                            (incoming_df[era_col] == x1) & (incoming_df.Knockout == y1)
                        ]
                        yinput = incoming_df[
                            (incoming_df[era_col] == x2) & (incoming_df.Knockout == y2)
                        ]

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

                        statistics_object_list = scipy.stats.mannwhitneyu(
                            x=xinput, y=yinput, nan_policy="omit"
                        )

                        left.append(
                            (x1, y1, x_mean, x_median, x_std, x_stderr, x_count)
                        )
                        right.append(
                            (x2, y2, y_mean, y_median, y_std, y_stderr, y_count)
                        )
                        floatvals.append(statistics_object_list[0])
                        pvals.append(statistics_object_list[1])

    dfResults["Left Status"] = [tup[0] for tup in left]
    dfResults["Left Knockout_Group"] = [tup[1] for tup in left]
    dfResults["Left Mean"] = [tup[2] for tup in left]
    dfResults["Left Median"] = [tup[3] for tup in left]
    dfResults["Left STDev"] = [tup[4] for tup in left]
    dfResults["Left STDerr"] = [tup[5] for tup in left]
    dfResults["Left Count"] = [tup[-1] for tup in left]

    dfResults["Right Status"] = [tup[0] for tup in right]
    dfResults["Right Knockout_Group"] = [tup[1] for tup in right]
    dfResults["Right Mean"] = [tup[2] for tup in right]
    dfResults["Right Median"] = [tup[3] for tup in right]
    dfResults["Right STDev"] = [tup[4] for tup in right]
    dfResults["Right STDerr"] = [tup[5] for tup in right]
    dfResults["Right Count"] = [tup[-1] for tup in right]

    dfResults["Floats"] = floatvals
    dfResults["pvals"] = pvals
    dfResults.sort_values(by="pvals", inplace=True)
    dfResults.reset_index(drop=True, inplace=True)

    display(dfResults)

    return dfResults


# Table containing mean read count data (with RNASeq Knockout annotation)
df5 = pd.read_csv(
    "Marten Knockout RNASeq 20250529/tables/full_source_genes_with_KOexpression.tsv",
    sep="\t",
    index_col="Identifier",
)


# Filter to genes flagged for GRCh38 expression plotting (relevant subset).
df5_nonexcl = df5[df5.GRCh38_Expression_Plotting == True]


# These serve as the tissue categories (functionally equivalent).
# Still iterate over disease status and evolutionary category.
mean_stats = ["DKO_KO_mean", "Hap1_WT_mean", "Wapl_KO_mean", "SCC4_KO_mean"]


df_meancount = df5_nonexcl[
    [
        "evolutionary_category",
        "evolutionary_category_3era",
        "disease_gene_inheritance_merged",
    ]
    + mean_stats
]


unique_eras = {
    "1-Ancient": "1-Ancient",
    "2-Metazoan": "2-Metazoan",
    "3-Chordate": "3-Chordate_5era",
    "4-Mammal": "4-Mammal_5era",
    "5-Primate": "5-Primate_5era",
}


df_meancount = df_meancount.rename(
    {"evolutionary_category": "evolutionary_category_5era"}, axis=1
)
df_meancount["evolutionary_category_5era"] = [
    unique_eras[xi] for xi in df_meancount.evolutionary_category_5era
]


melt_df = pd.melt(
    frame=df_meancount,
    id_vars=[
        "evolutionary_category_3era",
        "evolutionary_category_5era",
        "disease_gene_inheritance_merged",
    ],
    var_name="Knockout",
    value_name="Mean(Count)",
    ignore_index=False,
)
melt_df["joint_disease_era_3era"] = [
    f"{xi[0]}_{xi[1]}"
    for xi in zip(
        melt_df.evolutionary_category_3era, melt_df.disease_gene_inheritance_merged
    )
]
melt_df["joint_disease_era_5era"] = [
    f"{xi[0]}_{xi[1]}"
    for xi in zip(
        melt_df.evolutionary_category_5era, melt_df.disease_gene_inheritance_merged
    )
]


# Groups and arguments for the statistics method.
eragroups_3era = sorted(melt_df.evolutionary_category_3era.unique())
eragroups_5era = sorted(melt_df.evolutionary_category_5era.unique())

tissueCatKeys = sorted(melt_df.Knockout.unique())

joint_groups_3era = sorted(melt_df.joint_disease_era_3era.unique())
joint_groups_5era = sorted(melt_df.joint_disease_era_5era.unique())

disease_groups = sorted(melt_df.disease_gene_inheritance_merged.unique())


tissue_dict = {}

for xi in tissueCatKeys:
    tissue_dict[xi] = [xi]


# Calculations by knockout group and Evolutionary Category

df_era_3era = sig_table(
    eras_input=eragroups_3era,
    tissue_input=tissueCatKeys,
    incoming_df=melt_df,
    era_col="evolutionary_category_3era",
    do_mean=False,
    verbose=False,
)

# Already mean-aggregated; no need to recompute.
df_era_5era = sig_table(
    eras_input=eragroups_5era,
    tissue_input=tissueCatKeys,
    incoming_df=melt_df,
    era_col="evolutionary_category_5era",
    do_mean=False,
    verbose=False,
)


# Calculations by knockout group and disease status

df_disease = sig_table(
    eras_input=disease_groups,
    tissue_input=tissueCatKeys,
    incoming_df=melt_df,
    era_col="disease_gene_inheritance_merged",
    do_mean=False,
    verbose=True,
)


# Ultra-simplified significance stats — only between treatment groups.
# Fewer iterations: compare only within a single column (treatment groups).


def sig_table_simplified(tissue_input, incoming_df, verbose=False, do_mean=False):

    # Create empty lists that will become the dataframe

    dfResults = pd.DataFrame()

    left = []
    right = []
    floatvals = []
    pvals = []

    for y1_idx in range(len(tissue_input)):
        for y2 in tissue_input[y1_idx:]:
            y1 = tissue_input[y1_idx]
            if verbose:
                print(f"comparing ({y1}) to ({y2})")
            if y1 == y2:
                if verbose:
                    print("Duplicate - pass")
            else:

                xinput = incoming_df[(incoming_df.Knockout == y1)]
                yinput = incoming_df[(incoming_df.Knockout == y2)]

                if do_mean == True:
                    if "Name" in xinput.columns:
                        xinput = xinput.groupby("Name")["Mean(Count)"].mean()
                        yinput = yinput.groupby("Name")["Mean(Count)"].mean()
                    else:
                        xinput = xinput.groupby(xinput.index)["Mean(Count)"].mean()
                        yinput = yinput.groupby(yinput.index)["Mean(Count)"].mean()
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

                statistics_object_list = scipy.stats.mannwhitneyu(
                    x=xinput, y=yinput, nan_policy="omit"
                )

                left.append(
                    ("All_Genes", y1, x_mean, x_median, x_std, x_stderr, x_count)
                )
                right.append(
                    ("All_Genes", y2, y_mean, y_median, y_std, y_stderr, y_count)
                )
                floatvals.append(statistics_object_list[0])
                pvals.append(statistics_object_list[1])

    dfResults["Left Status"] = [tup[0] for tup in left]
    dfResults["Left Knockout_Group"] = [tup[1] for tup in left]
    dfResults["Left Mean"] = [tup[2] for tup in left]
    dfResults["Left Median"] = [tup[3] for tup in left]
    dfResults["Left STDev"] = [tup[4] for tup in left]
    dfResults["Left STDerr"] = [tup[5] for tup in left]
    dfResults["Left Count"] = [tup[-1] for tup in left]

    dfResults["Right Status"] = [tup[0] for tup in right]
    dfResults["Right Knockout_Group"] = [tup[1] for tup in right]
    dfResults["Right Mean"] = [tup[2] for tup in right]
    dfResults["Right Median"] = [tup[3] for tup in right]
    dfResults["Right STDev"] = [tup[4] for tup in right]
    dfResults["Right STDerr"] = [tup[5] for tup in right]
    dfResults["Right Count"] = [tup[-1] for tup in right]

    dfResults["Floats"] = floatvals
    dfResults["pvals"] = pvals
    dfResults.sort_values(by="pvals", inplace=True)
    dfResults.reset_index(drop=True, inplace=True)

    display(dfResults)

    return dfResults


# Evaluate, only comparing among knockout groups
df_knockout = sig_table_simplified(
    tissue_input=tissueCatKeys,
    incoming_df=melt_df,
    do_mean=False,
    verbose=True,
)


with pd.option_context(
    "display.max_rows",
    None,
):
    display(df_disease)


df_disease.value_counts(["Left Status", "Left Count"])


df_disease.value_counts(["Right Status", "Right Count"])


df_disease.value_counts(["Left Status", "Left Mean"])


df_era_3era.value_counts(["Left Status", "Left Count"])


# Let's make one merged, de-duplicated `df_era`.
# From `df_era_5era`, we can exclude those between 1-Ancient and 2-Metazoan.
# Original count is 130 comparisons for 5-era + knockout status.
df_era_5era_slim = df_era_5era[
    (df_era_5era["Left Status"].str.contains("_5era"))
    | (df_era_5era["Right Status"].str.contains("_5era"))
]
display(df_era_5era_slim)  # 108 unique analyses among 5era sets


df_era_3era["Analysis_Group"] = "3_Era"
df_era_5era["Analysis_Group"] = "5_Era_Unique"

df_era = pd.concat([df_era_3era, df_era_5era_slim])

# Combine all results to apply multiple-testing corrections together.

df_disease["Analysis_Group"] = "Disease"
# df_era['Analysis_Group'] = 'Era' # Already added
df_knockout["Analysis_Group"] = "Knockout_Only"

df_disease_joint = pd.concat([df_disease, df_era, df_knockout])


# Corrected p-values from statsmodels.
# fixedpvals_disease = statsmodels.stats.multitest.multipletests(df_disease['pvals'],method='fdr_bh')[1]
# fixedpvals_disease_era = statsmodels.stats.multitest.multipletests(df_disease_era['pvals'],method='fdr_bh')[1]
fixedpvals_disease_joint = statsmodels.stats.multitest.multipletests(
    df_disease_joint["pvals"], method="fdr_bh"
)[1]

df_disease_joint["Adjusted_Pvals_BH"] = fixedpvals_disease_joint
df_disease_joint["Adjusted_Pvals < 0.05"] = [
    fp < 0.05 for fp in fixedpvals_disease_joint
]


def comparison_key(df):
    return (
        df["Left Status"]
        + "_"
        + df["Left Knockout_Group"]
        + "_"
        + df["Right Status"]
        + "_"
        + df["Right Knockout_Group"]
    )


df_disease_joint["comparison_key"] = comparison_key(df_disease_joint)
df_disease_joint.set_index("comparison_key", inplace=True)


df_disease_joint.to_csv("marten_joint_knockout_significance_calculations.tsv", sep="\t")


era_terms = ["1-Ancient", "2-Metazoan", "3-Chordate"]
excluded_terms = [
    "Non-Disease",
    "Dominant",
    "Recessive",
    "4-Mammal_5era",
    "5-Primate_5era",
    "All_Genes",
]

df_disease_paper = df_disease_joint[df_disease_joint.Analysis_Group == "3_Era"]

df_disease_paper.to_csv("marten_knockout_era_significance.tsv", sep="\t")
