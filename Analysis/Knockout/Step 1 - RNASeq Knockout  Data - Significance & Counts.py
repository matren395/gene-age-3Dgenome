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
    knockout_input,
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

            for y1_idx in range(len(knockout_input)):
                for y2 in knockout_input[y1_idx:]:
                    y1 = knockout_input[y1_idx]
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
# Still iterate over statuses.
mean_stats = ["DKO_KO_mean", "Hap1_WT_mean", "Wapl_KO_mean", "SCC4_KO_mean"]


df_meancount = df5_nonexcl[
    [
        "evolutionary_category_3era",
    ]
    + mean_stats
]


melt_df = pd.melt(
    frame=df_meancount,
    id_vars=["evolutionary_category_3era"],
    var_name="Knockout",
    value_name="Mean(Count)",
    ignore_index=False,
)


# Groups and arguments for the statistics method.
eragroups_3era = sorted(melt_df.evolutionary_category_3era.unique())

tissueCatKeys = sorted(melt_df.Knockout.unique())


tissue_dict = {}

for xi in tissueCatKeys:
    tissue_dict[xi] = [xi]


# Calculations by knockout group and Evolutionary Category

df_era_3era = sig_table(
    eras_input=eragroups_3era,
    knockout_input=tissueCatKeys,
    incoming_df=melt_df,
    era_col="evolutionary_category_3era",
    do_mean=False,
    verbose=False,
)

# Corrected p-values from statsmodels.

fixedpvals_joint = statsmodels.stats.multitest.multipletests(
    df_era_3era["pvals"], method="fdr_bh"
)[1]

df_era_3era["Adjusted_Pvals_BH"] = fixedpvals_joint
df_era_3era["Adjusted_Pvals < 0.05"] = [fp < 0.05 for fp in fixedpvals_joint]


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


df_era_3era["comparison_key"] = comparison_key(df_era_3era)
df_era_3era.set_index("comparison_key", inplace=True)


df_era_3era.to_csv("marten_knockout_era_significance.tsv", sep="\t")
