#!/usr/bin/env python
# coding: utf-8

## Daniel Marten
## Distributions & Histograms for Disease Genes

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from cmapPy.pandasGEXpress.parse_gct import parse as tpm_parser

# Display options
pd.set_option("display.max_columns", None)


# Output from Step #0 Notebook
df5 = pd.read_csv(
    "Marten OMIM Disease Gene - 202503134/Tables/marten_s3a_omimdisease_20250314.tsv",
    sep="\t",
    index_col="Identifier",
)


date_name = "_centromere_3DG_20250902"  # variable for versioned plotting


# Enable or disable for saving plots
DO_PLOT = False


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


# Order for plotting (not alphabetical)
priority_map = {
    "Non-Disease": 0,
    "Non-Disease Autosomal": 0,
    "Non-Disease Sex-Linked": 4.5,
    "AD": 3,
    "XLD": 6,
    "AR": 2,
    "XLR": 5,
    "Mixed Autosomal": 7,
    "Mixed Autosomal Digenic": 9,
    "Mixed Sex-Linked": 8,
    "AD/AR MOI": 4,
    "Digenic": 10,
    "EXCLUDED-Unannotated": "EXCLUDED-Unannotated",
    "EXCLUDED-Intergenic Control": "EXCLUDED-Intergenic Control",
}


# Output first 20 colors from viridis palette for use in later palettes
viridis_xx = sns.color_palette("viridis", n_colors=20)


# Build universal palette for all disease work
disease_palette = {
    "Dominant": viridis_xx[-12],
    "Recessive": "pink",
    "AD": viridis_xx[-12],
    "AR": "pink",
    "XLR": viridis_xx[4],
    "XLD": viridis_xx[-4],
    "Non-Disease": "yellow",
    "Non-Disease Autosomal": "yellow",
    "Non-Disease Sex-Linked": "cyan",
}

# Colorblind colors
colorblind_11 = sns.color_palette("colorblind", n_colors=11)

disease_palette["XLR"] = [
    xi / 256 for xi in [80, 70, 164]
]  # Custom modification of XLR for custom color


# Work for plotting mixed statuses
disease_palette["AD/AR MOI"] = colorblind_11[4]

disease_palette["Mixed Autosomal"] = colorblind_11[-4]
disease_palette["Mixed Sex-Linked"] = colorblind_11[-3]
disease_palette["Mixed Autosomal Digenic"] = colorblind_11[-2]

disease_palette["Digenic"] = colorblind_11[1]

disease_palette["Mixed Dominant-Recessive"] = disease_palette["Mixed Autosomal"]
disease_palette["Mixed Dominant-Recessive Digenic"] = disease_palette[
    "Mixed Autosomal Digenic"
]


# Begin custom work for Histograms & Distribution Plots
df6_hist = df5[
    (~df5.annotation.str.contains("Unannotated"))
    & (~df5.annotation.str.contains("orf"))
]

# Create option to graph in only 3 Eras, since 4-Mammal and 5-Primate can be sparse

df6_hist["evolutionary_category_3era"] = [
    xi.replace("5-Primate", "3-Chordate").replace("4-Mammal", "3-Chordate")
    for xi in df6_hist["evolutionary_category"]
]

df6_hist["disease_gene_inheritance"] = [
    "Digenic" if xi == "digenic" else xi for xi in df6_hist["disease_gene_inheritance"]
]


# Establish priority for orders when plotting
df6_hist["priority"] = [priority_map[xi] for xi in df6_hist.disease_gene_inheritance]


# Histogram for disease counts by 3 and 5 eras
# for all disease statuses, only recessive & dominant, & R+D split by sex chr

for hist_config in ["evolutionary_category", "evolutionary_category_3era"]:

    era_str = "5_era"
    if "3era" in hist_config:
        era_str = "3_era"

    for inheritance in [
        "disease_gene_inheritance",
        "disease_gene_inheritance_merged",
        "disease_gene_inheritance_publication",
    ]:

        df6_hist_temp = df6_hist.copy()

        if inheritance != "disease_gene_inheritance":
            df6_hist_temp = df6_hist_temp[
                ~df6_hist_temp["disease_gene_inheritance_publication"].str.contains(
                    "EXCL"
                )
            ]

        # Overall by Disease Gene Status
        plt.figure(figsize=(10, 8), dpi=80)
        sns.histplot(
            data=df6_hist_temp[df6_hist_temp.disease_gene_binary == True].sort_values(
                by=["priority", hist_config, inheritance],
                ascending=[False, True, False],
            ),
            x=hist_config,
            stat="count",
            common_norm=False,
            linewidth=1,
            hue=inheritance,
            palette=disease_palette,
            multiple="stack",
        )

        labels = ["Recessive", "Dominant"]

        if inheritance == "disease_gene_inheritance_publication":
            labels = [
                "Autosomal Recessive",
                "Autosomal Dominant",
                "Sex-Linked Recessive",
                "Sex-Linked Dominant",
            ]
        elif inheritance == "disease_gene_inheritance":
            labels = [
                "Autosomal Recessive",
                "Autosomal Dominant",
                "AD/AR MOI",
                "Sex-Linked Recessive",
                "Sex-Linked Dominant",
                "Mixed Autosomal",
                "Mixed Sex-Linked",
                "Mixed Autosomal Digenic",
                "Digenic",
            ]

        xOnes, xTwos = plt.xticks()

        xTwos = ["Ancient", "Metazoan", "Chordate", "Mammal", "Primate"]

        if era_str == "3_era":
            xTwos = xTwos[:3]

        plt.xticks(xOnes, xTwos, size=20)
        plt.yticks(size=20)
        plt.xlabel("Evolutionary Category", size=32)
        plt.ylabel("Count", size=32)
        plt.legend(
            bbox_to_anchor=(1.01, 1),
            title="Disease Gene Counts",
            labels=labels,
            title_fontsize=20,
            fontsize=20,
        )

        image_name = f"fig4a-diseasegene-hist-{inheritance}-{era_str}-count-{date_name}"
        if DO_PLOT:
            plt.savefig(f"{image_name}.png", bbox_inches="tight")
            plt.savefig(f"{image_name}.pdf", bbox_inches="tight")

        plt.show()


# Histogram for Proportional Makeup by 3&5 Eras
# - for all disease statuses, only recessive & dominant, & R+D split by sex chr

for hist_config in ["evolutionary_category", "evolutionary_category_3era"]:

    era_str = "5_era"
    if "3era" in hist_config:
        era_str = "3_era"

    for inheritance in [
        "disease_gene_inheritance",
        "disease_gene_inheritance_merged",
        "disease_gene_inheritance_publication",
    ]:

        df6_hist_temp = df6_hist.copy()

        if inheritance != "disease_gene_inheritance":
            df6_hist_temp = df6_hist_temp[
                ~df6_hist_temp["disease_gene_inheritance_publication"].str.contains(
                    "EXCL"
                )
            ]

        # Overall by Disease Gene Status
        plt.figure(figsize=(10, 8), dpi=80)
        ax = sns.histplot(
            data=df6_hist_temp[df6_hist_temp.disease_gene_binary == True].sort_values(
                by=["priority", hist_config, inheritance],
                ascending=[False, True, False],
            ),
            x=hist_config,
            stat="count",
            common_norm=True,
            shrink=0.5,
            linewidth=1,
            hue=inheritance,
            palette=disease_palette,
            multiple="fill",
        )

        labels = ["Recessive", "Dominant"]

        if inheritance == "disease_gene_inheritance_publication":
            labels = [
                "Autosomal Recessive",
                "Autosomal Dominant",
                "Sex-Linked Recessive",
                "Sex-Linked Dominant",
            ]
        elif inheritance == "disease_gene_inheritance":
            labels = [
                "Autosomal Recessive",
                "Autosomal Dominant",
                "AD/AR MOI",
                "Sex-Linked Recessive",
                "Sex-Linked Dominant",
                "Mixed Autosomal",
                "Mixed Sex-Linked",
                "Mixed Autosomal Digenic",
                "Digenic",
            ]

        xOnes, xTwos = plt.xticks()

        xTwos = ["Ancient", "Metazoan", "Chordate", "Mammal", "Primate"]

        if era_str == "3_era":
            xTwos = xTwos[:3]

        plt.xticks(xOnes, xTwos, size=20)

        plt.yticks(size=20)
        plt.xlabel("Evolutionary Category", size=36)
        plt.ylabel("Proportion", size=32)
        plt.title(
            "Evolutionary Category Makeup\nby Disease Gene Status", size=32, pad=15
        )

        plt.legend(
            bbox_to_anchor=(1.01, 1),
            title="Disease Gene Status",
            labels=labels,
            title_fontsize=20,
            fontsize=20,
        )

        ax.spines["top"].set_visible(False)

        image_name = (
            f"fig4b-diseasegene-hist-{inheritance}-{era_str}-proportion-{date_name}"
        )
        if DO_PLOT:
            plt.savefig(f"{image_name}.png", bbox_inches="tight")
            plt.savefig(f"{image_name}.pdf", bbox_inches="tight")

        plt.show()


# Moving forward with plotting only by Binary Status, if yes-no
df6_hist = df6_hist[df6_hist.GRCh38_Expression_Plotting]

# Set this as a String instead of Boolean, due to some plotting particularities
df6_hist["disease_gene_binary_str"] = [
    "OMIM-Disease" if xi else "Non-Disease" for xi in df6_hist.disease_gene_binary
]


# Extremely basic palette for only plotting disease or non-disease
simple_palette = {"OMIM-Disease": disease_palette["Dominant"], "Non-Disease": "yellow"}


# Histogram for Counts by Eras (3 & 5) by OMIM Disease & Non-Disease Status

for hist_config in ["evolutionary_category", "evolutionary_category_3era"]:

    era_str = "5_era"
    if "3era" in hist_config:
        era_str = "3_era"

    for inheritance in ["disease_gene_binary_str"]:

        # Overall by Disease Gene Status
        plt.figure(figsize=(10, 8), dpi=80)
        sns.histplot(
            data=df6_hist.sort_values(
                by=[hist_config, inheritance], ascending=[True, False]
            ),
            x=hist_config,
            stat="count",
            common_norm=False,
            linewidth=1,
            hue=inheritance,
            multiple="stack",
            palette=simple_palette,
        )

        labels = ["Non-Disease", "OMIM-Disease"]

        xOnes, xTwos = plt.xticks()

        xTwos = ["Ancient", "Metazoan", "Chordate", "Mammal", "Primate"]

        if era_str == "3_era":
            xTwos = xTwos[:3]

        plt.xticks(xOnes, xTwos, size=20)

        plt.yticks(size=20)
        plt.xlabel("Evolutionary Category", size=36)
        plt.ylabel("Count", size=32)
        plt.legend(
            bbox_to_anchor=(1.01, 1),
            title="Disease Gene Status",
            labels=labels,
            title_fontsize=20,
            fontsize=20,
        )

        image_name = f"fig5a-diseasegene-hist-binary-count-{date_name}"
        if DO_PLOT:
            plt.savefig(f"{image_name}.png", bbox_inches="tight")
            plt.savefig(f"{image_name}.pdf", bbox_inches="tight")

            plt.show()


# Fig5 Hist for b:Proportion - Eras (3 & 5) by OMIM Disease & Non-Disease Status


for hist_config in ["evolutionary_category", "evolutionary_category_3era"]:

    era_str = "5_era"
    if "3era" in hist_config:
        era_str = "3_era"

    for inheritance in ["disease_gene_binary_str"]:

        # Overall by Disease Gene Status
        plt.figure(figsize=(10, 8), dpi=80)
        ax = sns.histplot(
            data=df6_hist.sort_values(
                by=[hist_config, inheritance], ascending=[True, False]
            ),
            x=hist_config,
            stat="percent",
            common_norm=True,
            linewidth=1,
            hue=inheritance,
            multiple="fill",
            shrink=0.5,
            palette=simple_palette,
        )

        ax.spines["top"].set_visible(False)

        labels = ["Non-Disease", "OMIM-Disease"]

        xOnes, xTwos = plt.xticks()

        xTwos = ["Ancient", "Metazoan", "Chordate", "Mammal", "Primate"]

        if era_str == "3_era":
            xTwos = xTwos[:3]

        plt.xticks(xOnes, xTwos, size=22)
        plt.yticks(size=20)
        plt.xlabel("Evolutionary Category", size=36)
        plt.ylabel("Proportion", size=32)
        plt.legend(
            bbox_to_anchor=(1.01, 1),
            title="Disease Gene Status",
            labels=labels,
            title_fontsize=20,
            fontsize=20,
        )

        image_name = f"fig5b-diseasegene-hist-binary-proportion-{date_name}"
        if DO_PLOT:
            plt.savefig(f"{image_name}.png", bbox_inches="tight")
            plt.savefig(f"{image_name}.pdf", bbox_inches="tight")

        plt.show()


# Fig5 Hist for c:Count - Eras (3 & 5) by OMIM Disease (split by R or D MOI) & Non-Disease Status

for hist_config in ["evolutionary_category", "evolutionary_category_3era"]:

    era_str = "5_era"
    if "3era" in hist_config:
        era_str = "3_era"

    for inheritance in ["disease_gene_inheritance_merged"]:

        # Overall by Disease Gene Status
        plt.figure(figsize=(10, 8), dpi=80)
        ax = sns.histplot(
            data=df6_hist.sort_values(by=[hist_config], ascending=[True]),
            x=hist_config,
            stat="count",
            common_norm=False,
            linewidth=1,
            hue=inheritance,
            multiple="stack",
            palette=disease_palette,
            hue_order=["Dominant", "Recessive", "Non-Disease"],
        )

        ax.spines["top"].set_visible(False)

        labels = ["Non-Disease", "Recessive", "Dominant"]

        xOnes, xTwos = plt.xticks()

        xTwos = ["Ancient", "Metazoan", "Chordate", "Mammal", "Primate"]

        if era_str == "3_era":
            xTwos = xTwos[:3]

        plt.xticks(xOnes, xTwos, size=20)
        plt.yticks(size=20)
        plt.xlabel("Evolutionary Category", size=36)
        plt.ylabel("Count", size=32)
        plt.legend(
            bbox_to_anchor=(1.01, 1),
            title="Disease Gene Status",
            labels=labels,
            title_fontsize=20,
            fontsize=20,
        )

        image_name = f"fig5c-diseasegene-hist-diseasestatus-count-revision-{date_name}"
        if True:
            plt.savefig(f"{image_name}.png", bbox_inches="tight")
            plt.savefig(f"{image_name}.pdf", bbox_inches="tight")

        plt.show()


# Histogram for d:Proportion - Eras (3 & 5) by OMIM Disease (split by R or D MOI) & Non-Disease Status

for hist_config in ["evolutionary_category", "evolutionary_category_3era"]:

    era_str = "5_era"
    if "3era" in hist_config:
        era_str = "3_era"

    for inheritance in ["disease_gene_inheritance_merged"]:

        # Overall by Disease Gene Status
        plt.figure(figsize=(10, 8), dpi=80)
        ax = sns.histplot(
            data=df6_hist.sort_values(by=[hist_config], ascending=[True]),
            x=hist_config,
            stat="percent",
            common_norm=True,
            linewidth=1,
            hue=inheritance,
            multiple="fill",
            shrink=0.5,
            palette=disease_palette,
            hue_order=["Dominant", "Recessive", "Non-Disease"],
        )

        ax.spines["top"].set_visible(False)

        labels = ["Non-Disease", "Recessive", "Dominant"]

        xOnes, xTwos = plt.xticks()

        xTwos = ["Ancient", "Metazoan", "Chordate", "Mammal", "Primate"]

        if era_str == "3_era":
            xTwos = xTwos[:3]

        plt.xticks(xOnes, xTwos, size=22)
        plt.yticks(size=20)
        plt.xlabel("Evolutionary Category", size=36)
        plt.ylabel("Proportion", size=32)
        plt.legend(
            bbox_to_anchor=(1.01, 1),
            title="Disease Gene Status",
            labels=labels,
            title_fontsize=20,
            fontsize=20,
        )

        image_name = (
            f"fig5d-diseasegene-hist-diseasestatus-proportion-revision-{date_name}"
        )
        if True:
            plt.savefig(f"{image_name}.png", bbox_inches="tight")
            plt.savefig(f"{image_name}.pdf", bbox_inches="tight")

        plt.show()


# Histogram for e:Count -
# Eras (3 & 5) by OMIM Disease (split by R/D and Sex-Linked Status) & Non-Disease Status

for hist_config in ["evolutionary_category", "evolutionary_category_3era"]:

    era_str = "5_era"
    if "3era" in hist_config:
        era_str = "3_era"

    for inheritance in ["disease_gene_inheritance_publication"]:

        # Overall by Disease Gene Status
        plt.figure(figsize=(10, 8), dpi=80)
        ax = sns.histplot(
            data=df6_hist.sort_values(by=[hist_config], ascending=[True]),
            x=hist_config,
            stat="count",
            common_norm=False,
            linewidth=1,
            hue=inheritance,
            multiple="stack",
            palette=disease_palette,
            hue_order=[
                "XLD",
                "XLR",
                "Non-Disease Sex-Linked",
                "AD",
                "AR",
                "Non-Disease Autosomal",
            ],
        )

        ax.spines["top"].set_visible(False)

        labels = [
            "Non-Disease Autosomal",
            "AR",
            "AD",
            "Non-Disease Sex-Linked",
            "XLR",
            "XLD",
        ]

        xOnes, xTwos = plt.xticks()

        xTwos = ["Ancient", "Metazoan", "Chordate", "Mammal", "Primate"]

        if era_str == "3_era":
            xTwos = xTwos[:3]

        plt.xticks(xOnes, xTwos, size=20)
        plt.yticks(size=20)
        plt.xlabel("Evolutionary Category", size=36)
        plt.ylabel("Count", size=32)
        plt.legend(
            bbox_to_anchor=(1.01, 1),
            title="Disease Gene Status",
            labels=labels,
            title_fontsize=20,
            fontsize=20,
        )

        image_name = f"fig5e-diseasegene-hist-sexlinked-count-revision-{date_name}"
        if DO_PLOT:
            plt.savefig(f"{image_name}.png", bbox_inches="tight")
            plt.savefig(f"{image_name}.pdf", bbox_inches="tight")

        plt.show()


# Fig5 Hist for d:Proportion -
# Eras (3 & 5) by OMIM Disease (split by R/D and Sex-Linked Status) & Non-Disease Status

for hist_config in ["evolutionary_category", "evolutionary_category_3era"]:

    era_str = "5_era"
    if "3era" in hist_config:
        era_str = "3_era"

    for inheritance in ["disease_gene_inheritance_publication"]:

        # Overall by Disease Gene Status
        plt.figure(figsize=(10, 8), dpi=80)
        ax = sns.histplot(
            data=df6_hist.sort_values(by=[hist_config], ascending=[True]),
            x=hist_config,
            stat="proportion",
            common_norm=True,
            linewidth=1,
            hue=inheritance,
            multiple="fill",
            palette=disease_palette,
            shrink=0.8,
            hue_order=[
                "XLD",
                "XLR",
                "Non-Disease Sex-Linked",
                "AD",
                "AR",
                "Non-Disease Autosomal",
            ],
        )

        ax.spines["top"].set_visible(False)

        labels = [
            "Non-Disease Autosomal",
            "AR",
            "AD",
            "Non-Disease Sex-Linked",
            "XLR",
            "XLD",
        ]

        xOnes, xTwos = plt.xticks()

        xTwos = ["Ancient", "Metazoan", "Chordate", "Mammal", "Primate"]

        if era_str == "3_era":
            xTwos = xTwos[:3]

        plt.xticks(xOnes, xTwos, size=20)
        plt.yticks(size=20)
        plt.xlabel("Evolutionary Category", size=36)
        plt.ylabel("Proportion", size=32)
        plt.legend(
            bbox_to_anchor=(1.01, 1),
            title="Disease Gene Status",
            labels=labels,
            title_fontsize=20,
            fontsize=20,
        )

        image_name = f"fig5f-diseasegene-hist-sexlinked-proportion-revision-{date_name}"
        if DO_PLOT:
            plt.savefig(f"{image_name}.png", bbox_inches="tight")
            plt.savefig(f"{image_name}.pdf", bbox_inches="tight")

        plt.show()


# Status for if a gene is on an autosomal chromosome or a sex-linked chromosome
df6_hist["Autosomal_Status"] = [
    "Sex-Linked" if xi in ["chrX", "chrY"] else "Autosomal" for xi in df6_hist["chr"]
]

df6_hist_autosomal = df6_hist[
    ~df6_hist["disease_gene_inheritance_publication"].str.contains("EXCL")
]

# Establish autosomal status as a custom priority we can sort by
df6_hist_autosomal["autosomal_priority"] = [
    0 if xi == "Non-Disease" else 1
    for xi in df6_hist_autosomal.disease_gene_inheritance_merged
]


# Histogram for a:Proportion of Disease Status (including non-disease) by Autosomal vs Sex-Linked

inheritance = "disease_gene_inheritance_merged"
hist_config = "Autosomal_Status"

# Overall proportion by Disease Gene Status
plt.figure(figsize=(10, 8), dpi=80)
ax = sns.histplot(
    data=df6_hist_autosomal.sort_values(
        by=["autosomal_priority", hist_config, inheritance],
        ascending=[False, True, True],
    ),
    x=hist_config,
    stat="percent",
    common_norm=True,
    linewidth=1,
    hue=inheritance,
    multiple="fill",
    palette=disease_palette,
    shrink=0.75,
)

ax.spines["top"].set_visible(False)

labels = ["Non-Disease", "Recessive", "Dominant"]

plt.xticks(size=22)
plt.yticks(size=20)
plt.xlabel("Chromosomal Status", size=36)
plt.ylabel("Proportion", size=32)
plt.legend(
    bbox_to_anchor=(1.01, 1),
    title="Disease Gene Status",
    labels=labels,
    title_fontsize=20,
    fontsize=20,
)

image_name = f"fig6a-diseasegene-hist-autosomal_v_sexlinked-proportion-{date_name}"
if DO_PLOT:
    plt.savefig(f"{image_name}.png", bbox_inches="tight")
    plt.savefig(f"{image_name}.pdf", bbox_inches="tight")

plt.show()


# Histogram for b:Count of Disease Status (including non-disease) by Autosomal vs Sex-Linked

inheritance = "disease_gene_inheritance_merged"
hist_config = "Autosomal_Status"

# Overall by Disease Gene Status
plt.figure(figsize=(10, 8), dpi=80)
ax = sns.histplot(
    data=df6_hist_autosomal.sort_values(
        by=["autosomal_priority", hist_config, inheritance],
        ascending=[False, True, True],
    ),
    x=hist_config,
    stat="count",
    linewidth=1,
    hue=inheritance,
    multiple="stack",
    palette=disease_palette,
    shrink=0.75,
)

ax.spines["top"].set_visible(False)


labels = ["Non-Disease", "Recessive", "Dominant"]

plt.xticks(size=22)
plt.yticks(size=20)
plt.xlabel("Chromosomal Status", size=36)
plt.ylabel("Count", size=32)
plt.legend(
    bbox_to_anchor=(1.01, 1),
    title="Disease Gene Status",
    labels=labels,
    title_fontsize=20,
    fontsize=20,
)
#
image_name = f"fig6b-diseasegene-hist-autosomal_v_sexlinked-count-{date_name}"
if DO_PLOT:
    plt.savefig(f"{image_name}.png", bbox_inches="tight")
    plt.savefig(f"{image_name}.pdf", bbox_inches="tight")


plt.show()


## Histogram Work

# Setting up some count for Chromosomes, since Integer Sorting =/= String Sorting
df6_hist_autosomal["smart_chr"] = [
    int(xi.replace("chr", "").replace("X", "23").replace("Y", "24"))
    for xi in df6_hist_autosomal["chr"]
]


# Plot Disease Gene Status Counts by Chromosome - Count

plt.figure(figsize=(10, 8), dpi=80)
ax = sns.histplot(
    data=df6_hist_autosomal.sort_values(
        by=["smart_chr", "autosomal_priority", hist_config, inheritance],
        ascending=[True, False, True, True],
    ),
    x="chr",
    stat="count",
    linewidth=1,
    hue=inheritance,
    multiple="stack",
    palette=disease_palette,
    shrink=0.75,
)

ax.spines["top"].set_visible(False)


labels = ["Non-Disease", "Recessive", "Dominant"]

plt.xticks(size=16, rotation=-90)
plt.yticks(size=20)
plt.xlabel("Chromosome", size=32)
plt.ylabel("Count", size=32)
plt.legend(
    bbox_to_anchor=(1.01, 1),
    title="Disease Gene Status",
    labels=labels,
    title_fontsize=20,
    fontsize=20,
)

image_name = f"fig7a-diseasegene-hist-chromosome-count-{date_name}"
if DO_PLOT:
    plt.savefig(f"{image_name}.png", bbox_inches="tight")
    plt.savefig(f"{image_name}.pdf", bbox_inches="tight")


plt.show()


# Plot Disease Gene Status by Chromosomeal Makeup - Proportion

plt.figure(figsize=(10, 8), dpi=80)
ax = sns.histplot(
    data=df6_hist_autosomal.sort_values(
        by=["smart_chr", "autosomal_priority", hist_config, inheritance],
        ascending=[True, False, True, True],
    ),
    x="chr",
    stat="percent",
    linewidth=1,
    hue=inheritance,
    multiple="fill",
    common_norm=True,
    palette=disease_palette,
    shrink=0.75,
)

ax.spines["top"].set_visible(False)


labels = ["Non-Disease", "Recessive", "Dominant"]

plt.xticks(size=16, rotation=-90)
plt.yticks(size=20)
plt.xlabel("Chromosome", size=32)
plt.ylabel("Proportion", size=32)
plt.title("Chromosome Make-Up by\nDisease Gene Status", size=32, pad=20)
plt.legend(
    bbox_to_anchor=(1.01, 1),
    title="Disease Gene Status",
    labels=labels,
    title_fontsize=20,
    fontsize=20,
)

image_name = f"fig7b-diseasegene-hist-chromosome-proportion-{date_name}"
if DO_PLOT:
    plt.savefig(f"{image_name}.png", bbox_inches="tight")
    plt.savefig(f"{image_name}.pdf", bbox_inches="tight")

plt.show()


# Plot Disease Gene Status by Chromosomeal Makeup - Proportion - WITH total trend line

plt.figure(figsize=(10, 8), dpi=80)
ax = sns.histplot(
    data=df6_hist_autosomal.sort_values(
        by=["smart_chr", "autosomal_priority", hist_config, inheritance],
        ascending=[True, False, True, True],
    ),
    x="chr",
    stat="percent",
    linewidth=1,
    hue=inheritance,
    multiple="fill",
    common_norm=True,
    palette=disease_palette,
    shrink=0.75,
)

ax.spines["top"].set_visible(False)


labels = ["Non-Disease", "Recessive", "Dominant"]

plt.xticks(size=16, rotation=-90)
plt.yticks(size=20)
plt.xlabel("Chromosome", size=32)
plt.ylabel("Proportion", size=32)
plt.title("Chromosome Make-Up by\nDisease Gene Status", size=32, pad=20)
plt.legend(
    bbox_to_anchor=(1.2, 1),
    title="Disease Gene Status",
    labels=labels,
    title_fontsize=20,
    fontsize=20,
)

chrom_counts = (
    df6_hist_autosomal.value_counts(["chr", "smart_chr"])
    .to_frame()
    .sort_values(by="smart_chr")
    .reset_index()
)

ax2 = plt.twinx()
sns.lineplot(data=chrom_counts, ax=ax2, x="chr", y="count", color="black", linewidth=2)

plt.ylabel("Total Gene Count", size=32)
plt.yticks(size=20)
plt.ylim([0, 2500])

image_name = (
    f"fig7c-diseasegene-hist-chromosome-proportion-trendline-revision-{date_name}"
)
if True:
    plt.savefig(f"{image_name}.png", bbox_inches="tight")
    plt.savefig(f"{image_name}.pdf", bbox_inches="tight")


plt.show()


# BEGIN WORK ON CHROMOSOMAL LOCALIZATION PLOTS


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
    "X": 156040895,
    "Y": 57227415,
}


# De-Index and Bin by 10MB Bins

df6_hist_deindex = df6_hist_autosomal.reset_index(drop=False).sort_values(
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


# Re-Establish Chromosomal Priority
df6_hist_deindex["autosomal_priority"] = [
    0 if xi == "Non-Disease" else 1
    for xi in df6_hist_deindex.disease_gene_inheritance_merged
]


# We do not include End Bins (with length less than 10MB) in our 10MB plots for consistency
# Using a collaborator's binning results as a starting point to establish this

# This contains all of the 10MB bins and coordinates we will be working with, generated from the same data
ten_mb_density = pd.read_excel(
    "/Users/marten/Downloads/for Dan/10MB_bins_start_20250506.xlsx"
)
ten_mb_density["chr"] = None
ten_mb_density["end"] = False

chr_keys = list(chromosome_lengths.keys())

# For each of our bins (10MB across the whole chromosome)
for xi, yi in ten_mb_density.iterrows():

    # If this is the first bin of the chromosome)
    if yi["start"] == 1:
        # Pop the chromosome name from the NCBI list , taking from the front of the list
        new_chr = chr_keys.pop(0)
        new_chr_str = f"chr{new_chr}"

    # Append this with NCBI-proofed keys
    ten_mb_density.loc[xi, "chr"] = new_chr_str

    # if the end of the bin as after the end of the chromosome
    if (yi["stop"] - chromosome_lengths[new_chr]) > 1:
        # Then this is a last bin! and mark it as such
        print("endpoint at: ", xi)
        ten_mb_density.loc[xi, "end"] = True

ends = ten_mb_density[ten_mb_density.end == True]

ends["bin_name"] = [
    f"{xi[0]}:{float(xi[1]//10_000_000)}" for xi in zip(ends.chr, ends.start)
]


# Mark if each gene is in an end bin
df6_hist_deindex["is_end"] = [
    xi in list(ends.bin_name) for xi in df6_hist_deindex["10MB_bin"]
]
df6_hist_deindex.is_end.value_counts()


# Dictionary of each bin and the count of genes which start within it
bin_count_dict = df6_hist_deindex["10MB_bin"].value_counts().to_dict()


# DataFrame of each bin and the count of genes which start within it
bin_count_frame = df6_hist_deindex["10MB_bin"].value_counts().to_frame()


# DataFrame of Age Breakdown within each bin
age_data = (
    df6_hist_deindex[df6_hist_deindex.is_end == False]
    .value_counts(["10MB_bin", "evolutionary_category", "evolutionary_category_3era"])
    .to_frame()
    .reset_index()
)


# Pivot table constructed from prior table, containing the counts for each era within each bin
pivot_bin_age = pd.pivot(
    age_data.reset_index(drop=True),
    index="10MB_bin",
    columns=["evolutionary_category"],
    values="count",
).fillna(0.0)
pivot_bin_age["sum"] = pivot_bin_age.sum(axis=1)

# Counts for anything either Chordate, Mammal, or Primate
pivot_bin_age["3-Chordate-Or-Younger"] = (
    pivot_bin_age["3-Chordate"] + pivot_bin_age["4-Mammal"] + pivot_bin_age["5-Primate"]
)


# Also calculate the proportion of each evolutionary category within each bin
for xi in pivot_bin_age.columns:
    description = f"{xi}_proportion"
    pivot_bin_age[description] = pivot_bin_age[xi] / pivot_bin_age["sum"]


# Similarly, construct and format bin data by Dominant, Recessive, and Non-Disease Status

bin_data = (
    df6_hist_deindex[df6_hist_deindex.is_end == False]
    .value_counts(["10MB_bin", "disease_gene_inheritance_merged"])
    .to_frame()
    .reset_index()
)

pivot_bin = pd.pivot(
    bin_data.reset_index(drop=True),
    index="10MB_bin",
    columns="disease_gene_inheritance_merged",
    values="count",
)

pivot_bin["sum"] = pivot_bin.sum(axis=1)
pivot_bin["dominant_norm"] = pivot_bin["Dominant"] / pivot_bin["sum"]
pivot_bin["recessive_norm"] = pivot_bin["Recessive"] / pivot_bin["sum"]
pivot_bin["nondisease_norm"] = pivot_bin["Non-Disease"] / pivot_bin["sum"]


# Asign simple integer counts for Chromosome and Start position, to be easily sorted.
# as a string, 10,000,000 would be sorted before 9

pivot_bin["smart_chr"] = [
    int(xi.split(":")[0].replace("chr", "").replace("X", "23").replace("Y", "24"))
    for xi in pivot_bin.index
]
pivot_bin["smart_start"] = [float(xi.split(":")[1]) for xi in pivot_bin.index]

pivot_bin = (
    pivot_bin.sort_values(by=["smart_chr", "smart_start"])
    .reset_index(drop=False)
    .set_index("10MB_bin")
)


# Append pivot bin table (with disease information) with table with genetic age information, and drop redundant columns
pivot_bin = pivot_bin.join(pivot_bin_age.drop(["sum", "sum_proportion"], axis=1))


# Similarly, assign a numerical value for the start of each bin within each entry on the table
df6_hist_deindex["smart_start"] = [
    float(xi.split(":")[1]) for xi in df6_hist_deindex["10MB_bin"]
]


# Logic for annotating importing and annotating centromere regions onto 10MB localization graphs

centromere_df = pd.read_csv(
    "~/Downloads/Modeled_regions_for_GRCh38.tsv",
    sep="\t",
    header=0,
    names=["region_name", "chr", "start", "stop", "length"],
)
centromere_df["chr"] = "chr" + centromere_df["chr"]
centromere_df = centromere_df[centromere_df.region_name.str.contains("CEN")]
centromere_df["start_bin_name"] = [
    c_i[0] + ":" + str(float(c_i[1] // 10_000_000))
    for c_i in zip(centromere_df["chr"], centromere_df["start"])
]
centromere_df["end_bin_name"] = [
    c_i[0] + ":" + str(float(c_i[1] // 10_000_000))
    for c_i in zip(centromere_df["chr"], centromere_df["stop"])
]

centromere_list = list(
    set(list(centromere_df.start_bin_name) + list(centromere_df.end_bin_name))
)


# Disease Gene Distribution (Non-Disease, Recessive, Dominant) in 10MB Bins (Excluding Ends)
# WITH additional trendline to show total gene count throughout all bin

plt.figure(figsize=(50, 8), dpi=80)
ax = sns.histplot(
    data=df6_hist_deindex[df6_hist_deindex.is_end == False].sort_values(
        by=[
            "smart_chr",
            "smart_start",
            "autosomal_priority",
            hist_config,
            "disease_gene_inheritance_merged",
        ],
        ascending=[True, True, False, True, True],
    ),
    x="10MB_bin",
    stat="percent",
    linewidth=1,
    hue="disease_gene_inheritance_merged",
    multiple="fill",
    common_norm=True,
    palette=disease_palette,
)

ax.spines["top"].set_visible(False)


labels = ["Non-Disease", "Recessive", "Dominant"]

xOnes, xTwos = plt.xticks()


# Written so as to also be able to include centromere annotations
xTwo_alt = [xi.get_text() for xi in xTwos]
xTwo_contig_unique = []

centromeric_indices = []

for x2_i in range(len(xTwo_alt)):

    x2_i_val = xTwo_alt[x2_i]
    x2_i_contig = x2_i_val.split(":")[0]

    tick_temp = None

    if x2_i_contig in xTwo_contig_unique:
        tick_temp = ""
        # print('Nothing cool...')
    else:
        tick_temp = x2_i_contig
        xTwo_contig_unique.append(x2_i_contig)

    # separate check
    if x2_i_val in centromere_list:
        tick_temp = "-" + tick_temp
        centromeric_indices.append(x2_i)

    xTwo_alt[x2_i] = tick_temp

plt.xticks(xOnes, xTwo_alt, size=26, rotation=-90)
plt.yticks(size=20)
plt.xlim(-1, len(xOnes))

plt.xticks(xOnes, xTwo_alt, size=26, rotation=-90)
plt.yticks(size=20)
plt.xlim(-1, len(xOnes))

plt.xlabel("Chromosome", size=32)
plt.ylabel("Proportion of\nDisease Gene Status", size=32)
plt.title("Disease Gene Distribution\nBy 10MB Bin\nExcluding End Bins", size=32, pad=20)
plt.legend(
    bbox_to_anchor=(1.15, 1.15),
    title="Disease Gene Status",
    labels=labels,
    title_fontsize=20,
    fontsize=20,
)

# Code to plot total counts throughout bins , to give context to bin proportions
ax2 = plt.twinx()
sns.lineplot(data=pivot_bin, color="black", ax=ax2, x="10MB_bin", y="sum", linewidth=2)

plt.ylabel("Total Count", size=32)
plt.yticks(size=20)
plt.ylim([0, 400])

image_name = f"fig8b-hist-trendline-wide-distribution-diseasegenedistribution-10MB-noendbins-{date_name}"
if True:
    plt.savefig(f"{image_name}.png", bbox_inches="tight")
    plt.savefig(f"{image_name}.pdf", bbox_inches="tight")

plt.show()


# Transformations and adjustments before plotting by Evolutioinary Category
era_five_colours["2-Metazoan"] = era_five_colours["2-Metazoa"]

df6_hist_deindex = df6_hist_deindex.sort_values(
    by="evolutionary_category", ascending=True
)


# Some additional and alternative color work
viridis_palette = sns.color_palette("viridis", 10)
muted_palette = sns.color_palette("muted", 10)
colorblind_palette = sns.color_palette("colorblind", 10)

alternative_palette = {
    "1-Ancient": muted_palette[2],
    "2-Metazoan": muted_palette[1],
    "3-Chordate": muted_palette[0],
    "4-Mammal": muted_palette[6],
    "5-Primate": colorblind_palette[-2],
}


# Chromosome length reference

chr_df = pd.DataFrame.from_dict(data=chromosome_lengths, orient="index").reset_index(
    drop=False
)
chr_df = chr_df.rename({0: "length", "index": "Chromosome"}, axis=1)
chr_df["length_mb"] = [xi / 1_000_000 for xi in chr_df["length"]]


# Chromosomal Length
# Derived from Total Length for Human Genome Assembly GRCh38
# https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh38

plt.figure(figsize=(16, 8), dpi=80)

ax = sns.barplot(
    data=chr_df,
    x="Chromosome",
    y="length_mb",
    edgecolor="black",
    width=0.8,
    color="C1",
    linewidth=1.25,
)
plt.xlabel("Chromosome", size=22)
plt.ylabel("Length (Mb)", size=22)

plt.xticks(size=20)
plt.yticks(size=20)

image_name = f"fig14-chromsomal-length-{date_name}"
if DO_PLOT:
    plt.savefig(f"{image_name}.png", bbox_inches="tight")
    plt.savefig(f"{image_name}.pdf", bbox_inches="tight")

plt.show()


# Plot Disease Gene Status by Chromosomeal Makeup - WITH total trend line

plt.figure(figsize=(10, 8), dpi=80)
ax = sns.histplot(
    data=df6_hist_autosomal.sort_values(
        by=["smart_chr", "autosomal_priority", hist_config, inheritance],
        ascending=[True, False, True, True],
    ),
    x="chr",
    stat="percent",
    linewidth=1,
    hue=inheritance,
    multiple="fill",
    common_norm=True,
    palette=disease_palette,
    shrink=0.75,
)

ax.spines["top"].set_visible(False)


labels = ["Non-Disease", "Recessive", "Dominant"]

#         if 'merged' not in inheritance:
#             labels = ['Autosomal Dominant','Autosomal Recessive','X-Linked Dominant', 'X-Linked Recessive']

plt.xticks(size=16, rotation=-90)
plt.yticks(size=20)
plt.xlabel("Chromosome", size=32)
plt.ylabel("Proportion", size=32)
plt.title("Chromosome Make-Up by\nDisease Gene Status", size=32, pad=20)
plt.legend(
    bbox_to_anchor=(1.2, 1),
    title="Disease Gene Status",
    labels=labels,
    title_fontsize=20,
    fontsize=20,
)

chrom_counts = (
    df6_hist_autosomal.value_counts(["chr", "smart_chr"])
    .to_frame()
    .sort_values(by="smart_chr")
    .reset_index()
)

ax2 = plt.twinx()
sns.lineplot(data=chrom_counts, ax=ax2, x="chr", y="count", color="black", linewidth=2)

plt.ylabel("Total Gene Count", size=32)
plt.yticks(size=20)
plt.ylim([0, 2500])

image_name = (
    f"fig7c-diseasegene-hist-chromosome-proportion-trendline-revision-{date_name}"
)
if True:
    plt.savefig(f"{image_name}.png", bbox_inches="tight")
    plt.savefig(f"{image_name}.pdf", bbox_inches="tight")


plt.show()


# Write out annotated version of our data used for above count and histogram plotting
df6_hist_autosomal.to_csv("histogram_working_table.tsv", sep="\t")
