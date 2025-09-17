#!/usr/bin/env python
# coding: utf-8

## Daniel Marten
## Plot Expression by OMIM Disease Gene status

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

date_name = "_20250502"  # variable for versioned plotting


# Enable or disable for saving plots
DO_PLOT = True


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


# Only select eras and six major tissue categories; no metadata included
df6_plotting = df5[~df5.annotation.str.contains("Unannotated")]
df6_plotting = df6_plotting[
    [
        "BRAIN",
        "ECTO",
        "MESO",
        "ENDO",
        "TESTIS",
        "OVARY",
        "evolutionary_category",
        "disease_gene_binary",
        "disease_gene_inheritance",
        "disease_gene_inheritance_merged",
        "disease_gene_inheritance_publication",
    ]
]
df6_plotting = df6_plotting[
    ~df6_plotting.evolutionary_category.str.contains("Intergenic")
]

# Graphing in only three eras
df6_plotting["evolutionary_category"] = [
    xi.replace("5-Primate", "3-Chordate").replace("4-Mammal", "3-Chordate")
    for xi in df6_plotting["evolutionary_category"]
]

# Priority for plotting order
df6_plotting["priority"] = [
    priority_map[xi] for xi in df6_plotting.disease_gene_inheritance
]

# Rename for simpler plotting order - alphanumerical order!
df6 = df6_plotting.rename(
    columns={
        "BRAIN": "1-Brain",
        "ECTO": "2-Ecto",
        "MESO": "3-Meso",
        "ENDO": "4-Endo",
        "TESTIS": "6-Testis",
        "OVARY": "5-Ovary",
    }
)

# Melt into desired DataFrame
df6 = pd.melt(
    frame=df6,
    id_vars=[
        "evolutionary_category",
        "disease_gene_binary",
        "disease_gene_inheritance",
        "disease_gene_inheritance_merged",
        "priority",
        "disease_gene_inheritance_publication",
    ],
    var_name="Germ",
    value_name="Mean Count",
    ignore_index=False,
)


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
# relevant for our own work and checking
disease_palette["AD/AR MOI"] = colorblind_11[4]

disease_palette["Mixed Autosomal"] = colorblind_11[-4]
disease_palette["Mixed Sex-Linked"] = colorblind_11[-3]
disease_palette["Mixed Autosomal Digenic"] = colorblind_11[-2]

disease_palette["Digenic"] = colorblind_11[1]

disease_palette["Mixed Dominant-Recessive"] = disease_palette["Mixed Autosomal"]
disease_palette["Mixed Dominant-Recessive Digenic"] = disease_palette[
    "Mixed Autosomal Digenic"
]


# PLOT SET 1 - For all tissue categories, Mean Count by Era, split by Status
# a - for all Status options (mixed, AD/DR)

for germ_i in list(df6.Germ.unique()):

    germ_text = germ_i.split("-")[-1]

    config = "evolutionary_category"

    df6_temp = df6[df6.Germ == germ_i]
    df6_temp = df6_temp.sort_values(by=["priority", config, "Germ"])

    plt.figure(figsize=(12, 8), dpi=80)
    sns.boxplot(
        data=df6_temp,
        x=config,
        y="Mean Count",
        hue="disease_gene_inheritance",
        showfliers=False,
        linewidth=3,
        palette=disease_palette,
    )
    plt.xlabel("Evolutionary Era", size=32)
    plt.ylabel("Mean Count", size=32)
    xOnes, xTwos = plt.xticks()

    xfont = 26
    if config == "evolutionary_category":
        xTwos = ["Ancient", "Metazoan", "Chordate"]
    #     lablos = ['Brain','Ecto','Meso','Endo','Ovary','Testis']
    plt.xticks(xOnes, xTwos, size=xfont)
    plt.yticks(size=27)
    plt.legend(
        title="Disease Gene", title_fontsize=22, fontsize=22, bbox_to_anchor=(1.00, 1)
    )
    plt.title(
        f"{germ_text} Tissue Mean Count per Disease Gene Status\nby Evolutionary Era",
        size=32,
        pad=15,
    )

    image_name = f"fig1a-germstratified-{germ_text}-diseasegene-meancount-era-allstratifications-{date_name}"
    if DO_PLOT:
        plt.savefig(f"{image_name}.png", bbox_inches="tight")
        plt.savefig(f"{image_name}.pdf", bbox_inches="tight")

    plt.show()

    # PLOT SET 1 - Tissue Category Mean Count by Era, split by Status
    # ^^ b - for only Non-Disease, Recessive, Dominant statuses

    df6_temp = df6[df6.Germ == germ_i]
    df6_temp = df6_temp.sort_values(by=["priority", config, "Germ"])

    df6_temp = df6_temp[
        ~df6_temp.disease_gene_inheritance_publication.str.contains("EXCLUDED")
    ]

    plt.figure(figsize=(12, 8), dpi=80)
    sns.boxplot(
        data=df6_temp,
        x=config,
        y="Mean Count",
        hue="disease_gene_inheritance_merged",
        showfliers=False,
        linewidth=3,
        palette=disease_palette,
    )
    plt.xlabel("Evolutionary Era", size=32)
    plt.ylabel("Mean Count", size=32)
    xOnes, xTwos = plt.xticks()

    xfont = 26
    if config == "evolutionary_category":
        xTwos = ["Ancient", "Metazoan", "Chordate"]
    #     lablos = ['Brain','Ecto','Meso','Endo','Ovary','Testis']
    plt.xticks(xOnes, xTwos, size=xfont)
    plt.yticks(size=27)
    plt.legend(
        title="Disease Gene", title_fontsize=22, fontsize=22, bbox_to_anchor=(1.00, 1)
    )
    plt.title(
        f"{germ_text} Tissue Mean Count per Disease Gene Status\nby Evolutionary Era",
        size=32,
        pad=15,
    )

    image_name = (
        f"fig1b-germstratified-{germ_text}-diseasegene-meancount-era-simple-{date_name}"
    )
    if DO_PLOT:
        plt.savefig(f"{image_name}.png", bbox_inches="tight")
        plt.savefig(f"{image_name}.pdf", bbox_inches="tight")

    plt.show()
    # fig notes:

    # PLOT SET 1 - Tissue Category Mean Count by Era, split by Status
    # ^^ c - with Sex-Linked Statuses split out

    df6_temp = df6[df6.Germ == germ_i]
    df6_temp = df6_temp.sort_values(by=["priority", config, "Germ"])

    df6_temp = df6_temp[
        ~df6_temp.disease_gene_inheritance_publication.str.contains("EXCLUDED")
    ]

    plt.figure(figsize=(12, 8), dpi=80)
    sns.boxplot(
        data=df6_temp,
        x=config,
        y="Mean Count",
        hue="disease_gene_inheritance_publication",
        showfliers=False,
        linewidth=3,
        palette=disease_palette,
    )
    plt.xlabel("Evolutionary Era", size=32)
    plt.ylabel("Mean Count", size=32)
    xOnes, xTwos = plt.xticks()

    xfont = 26
    if config == "evolutionary_category":
        xTwos = ["Ancient", "Metazoan", "Chordate"]
    plt.xticks(xOnes, xTwos, size=xfont)
    plt.yticks(size=27)
    plt.legend(
        title="Disease Gene", title_fontsize=22, fontsize=22, bbox_to_anchor=(1.00, 1)
    )
    plt.title(
        f"{germ_text} Tissue Mean Count per Disease Gene Status\nby Evolutionary Era",
        size=32,
        pad=15,
    )

    image_name = f"fig1c-germstratified-{germ_text}-diseasegene-meancount-era-sexsplit-{date_name}"
    if DO_PLOT:
        plt.savefig(f"{image_name}.png", bbox_inches="tight")
        plt.savefig(f"{image_name}.pdf", bbox_inches="tight")

    plt.show()


# PLOT SET 2 - Mean Counts by Germ Layer, Split by Disease Gene Status
# (flipped version of set #1)
# ^^ a - with all statuses

for config in ["evolutionary_category"]:

    df6_temp = df6.sort_values(by=["priority", config, "Germ"]).copy()
    df6_temp["Germ"] = [germitext.split("-")[-1] for germitext in df6_temp["Germ"]]
    df6_temp[config] = [configii.split("-")[-1] for configii in df6_temp[config]]

    if config == "evolutionary_category":
        palette_i = [yi for xi, yi in era_five_colours.items()]

    plt.figure(figsize=(26, 16), dpi=80)
    sns.boxplot(
        data=df6_temp,
        x="Germ",
        y="Mean Count",
        hue="disease_gene_inheritance",
        showfliers=False,
        linewidth=3,
        palette=palette_i,
    )
    plt.xlabel("Germ Layer", size=46)
    plt.ylabel("Mean Counts", size=46)

    xOnes, xTwos = plt.xticks()

    plt.xticks(size=46)
    plt.yticks(size=46)
    plt.legend(
        bbox_to_anchor=(1.01, 1), title="Disease Gene", title_fontsize=36, fontsize=32
    )
    plt.title(
        "Tissue Mean Counts per Disease Gene Status\nby Germ Layer", size=46, pad=15
    )
    plt.subplots_adjust(bottom=0.2, left=0.1)

    image_name = f"fig2a-diseasegene-germlayer-allstratifications-{date_name}"
    if DO_PLOT:
        plt.savefig(f"{image_name}.png", bbox_inches="tight")
        plt.savefig(f"{image_name}.pdf", bbox_inches="tight")

    plt.show()


# PLOT SET 2 - Mean Counts by Germ Layer, Split by Disease Gene Status
# (flipped version of set #1)
# ^^ b - with only 3 statuses: Non-Disease, Recessive, and Dominant

for config in ["evolutionary_category"]:

    df6_temp = df6.sort_values(by=["priority", config, "Germ"]).copy()
    df6_temp["Germ"] = [germitext.split("-")[-1] for germitext in df6_temp["Germ"]]
    df6_temp[config] = [configii.split("-")[-1] for configii in df6_temp[config]]

    df6_temp = df6_temp[
        ~df6_temp.disease_gene_inheritance_publication.str.contains("EXCL")
    ]

    plt.figure(figsize=(16, 12), dpi=80)
    sns.boxplot(
        data=df6_temp,
        x="Germ",
        y="Mean Count",
        hue="disease_gene_inheritance_merged",
        showfliers=False,
        linewidth=3,
        palette=disease_palette,
    )
    plt.xlabel("Germ Layer", size=46)
    plt.ylabel("Mean Counts", size=46)

    xOnes, xTwos = plt.xticks()

    plt.xticks(size=46)
    plt.yticks(size=46)
    plt.legend(
        bbox_to_anchor=(1.01, 1), title="Disease Gene", title_fontsize=36, fontsize=32
    )
    plt.title(
        "Tissue Mean Counts per Disease Gene Status\nby Germ Layer", size=46, pad=15
    )
    plt.subplots_adjust(bottom=0.2, left=0.1)

    image_name = f"fig2b-diseasegene-germlayer-simple-{date_name}"
    if DO_PLOT:
        plt.savefig(f"{image_name}.png", bbox_inches="tight")
        plt.savefig(f"{image_name}.pdf", bbox_inches="tight")

    plt.show()


# PLOT SET 2 - Mean Counts by Germ Layer, Split by Disease Gene Status
# (flipped version of set #1)
# ^^ c - with split out Sex Linked statuses

for config in ["evolutionary_category"]:

    df6_temp = df6.sort_values(by=["priority", config, "Germ"]).copy()
    df6_temp["Germ"] = [germitext.split("-")[-1] for germitext in df6_temp["Germ"]]
    df6_temp[config] = [configii.split("-")[-1] for configii in df6_temp[config]]

    df6_temp = df6_temp[
        ~df6_temp.disease_gene_inheritance_publication.str.contains("EXCL")
    ]

    plt.figure(figsize=(16, 12), dpi=80)
    sns.boxplot(
        data=df6_temp,
        x="Germ",
        y="Mean Count",
        hue="disease_gene_inheritance_publication",
        showfliers=False,
        linewidth=3,
        palette=disease_palette,
    )
    plt.xlabel("Germ Layer", size=46)
    plt.ylabel("Mean Counts", size=46)

    xOnes, xTwos = plt.xticks()

    plt.xticks(size=46)
    plt.yticks(size=46)
    plt.legend(
        bbox_to_anchor=(1.01, 1), title="Disease Gene", title_fontsize=36, fontsize=32
    )
    plt.title(
        "Tissue Mean Counts per Disease Gene Status\nby Germ Layer", size=46, pad=15
    )
    plt.subplots_adjust(bottom=0.2, left=0.1)

    if config == "evoera5":
        config = "five_eras"

    image_name = f"fig2c-diseasegene-germlayer-sexsplit-{date_name}"
    if DO_PLOT:
        plt.savefig(f"{image_name}.png", bbox_inches="tight")
        plt.savefig(f"{image_name}.pdf", bbox_inches="tight")

    plt.show()


# PLOT SET 3 - Total Brain - a:for all statuses

df6_temp = df6[df6.Germ == "1-Brain"]
df6_temp = df6_temp.sort_values(by=["priority", config, "Germ"])
df6_temp["Germ"] = [germitext.split("-")[-1] for germitext in df6_temp["Germ"]]

plt.figure(figsize=(10, 8), dpi=80)
sns.boxplot(
    data=df6_temp,
    x="disease_gene_inheritance",
    y="Mean Count",
    showfliers=False,
    linewidth=3,
    palette=disease_palette,
)
plt.xlabel("Disease Gene Status", size=32)
plt.ylabel("Mean Count", size=32)
xOnes, xTwos = plt.xticks()

xfont = 20
plt.xticks(xOnes, xTwos, size=xfont, rotation=-90)
plt.yticks(size=27)
plt.title("Brain Mean Count per Disease Gene Status", size=32, pad=15)

image_name = f"fig3a-diseasegene-brainmeancount-noera-allstratifications-{date_name}"
if DO_PLOT:
    plt.savefig(f"{image_name}.png", bbox_inches="tight")
    plt.savefig(f"{image_name}.pdf", bbox_inches="tight")

plt.show()


# PLOT SET 3 - Total Brain - b:for only 3 statuses (Non-Disease, Recessive, and Dominant)

df6_temp = df6[df6.Germ == "1-Brain"]
df6_temp = df6_temp.sort_values(by=["priority", config, "Germ"])
df6_temp["Germ"] = [germitext.split("-")[-1] for germitext in df6_temp["Germ"]]

df6_temp = df6_temp[~df6_temp.disease_gene_inheritance_publication.str.contains("EXCL")]

plt.figure(figsize=(10, 8), dpi=80)
sns.boxplot(
    data=df6_temp,
    x="disease_gene_inheritance_merged",
    y="Mean Count",
    showfliers=False,
    linewidth=3,
    palette=disease_palette,
)
plt.xlabel("Disease Gene Status", size=32)
plt.ylabel("Mean Count", size=32)
xOnes, xTwos = plt.xticks()

xfont = 22
plt.xticks(xOnes, xTwos, size=xfont, rotation=-90)
plt.yticks(size=27)
plt.title("Brain Mean Count per Disease Gene Status", size=32, pad=15)

image_name = f"fig3b-diseasegene-brainmeancount-noera-simple-{date_name}"
if DO_PLOT:
    plt.savefig(f"{image_name}.png", bbox_inches="tight")
    plt.savefig(f"{image_name}.pdf", bbox_inches="tight")

plt.show()


# PLOT SET 3 - Total Brain - c:with sex-linked statuses split out

df6_temp = df6[df6.Germ == "1-Brain"]
df6_temp = df6_temp.sort_values(by=["priority", config, "Germ"])
df6_temp["Germ"] = [germitext.split("-")[-1] for germitext in df6_temp["Germ"]]

df6_temp = df6_temp[~df6_temp.disease_gene_inheritance_publication.str.contains("EXCL")]

plt.figure(figsize=(10, 8), dpi=80)
sns.boxplot(
    data=df6_temp,
    x="disease_gene_inheritance_publication",
    y="Mean Count",
    showfliers=False,
    linewidth=3,
    palette=disease_palette,
)
plt.xlabel("Disease Gene Status", size=32)
plt.ylabel("Mean Count", size=32)
xOnes, xTwos = plt.xticks()

xfont = 22
# if config=='evolutionary_category':
#     xTwos = ['Ancient','Metazoan','Chordate']
# #     lablos = ['Brain','Ecto','Meso','Endo','Ovary','Testis']
for xi in range(len(xTwos)):
    texto = xTwos[xi].get_text()
    xTwos[xi] = texto.replace("Disease", "Disease\n")

plt.xticks(xOnes, xTwos, size=xfont, rotation=-90)
plt.yticks(size=27)
# plt.legend(title='Neural',title_fontsize=22,fontsize=22,bbox_to_anchor=(1.00, 1))
plt.title("Brain Mean Count per Disease Gene Status", size=32, pad=15)

image_name = f"fig3c-diseasegene-brainmeancount-noera-sexsplit-{date_name}"
if DO_PLOT:
    plt.savefig(f"{image_name}.png", bbox_inches="tight")
    plt.savefig(f"{image_name}.pdf", bbox_inches="tight")

plt.show()
