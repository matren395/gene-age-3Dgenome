#!/usr/bin/env python
# coding: utf-8

## Daniel Marten
## Read in OMIM Disease Data, then annotate to existing tables

import math
import os
import re
import statistics
import sys
import warnings
from collections import Counter

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import qtl.io as io
import qtl.norm as norm
import scipy
import scipy.stats as scistats
import seaborn as sns
import sklearn
import statsmodels.api as sm
import statsmodels.stats.weightstats as sm_stats
from cmapPy.pandasGEXpress.parse_gct import parse as tpm_parser
from matplotlib.pyplot import figure
from statsmodels.formula.api import ols as formula_OLS

pd.set_option("display.max_columns", None)


# Import list of OMIM Disease Genes, as provided by Sarah Stenton
disease_list = pd.read_csv(
    "/Users/marten/Downloads/omim_disease_genes.txt",
    sep="\t",
)


# Analysis does not include Phenotypes denoted with ?, [, or {
# Analysis does not include digenic modes of inheritance
# And classifies X-Linked as a type of X-Linked Dominant (for comparisons to Autosomal)
disease_list["Phenotype_Start"] = [xi[0] for xi in disease_list.Phenotypes]
disease_list = disease_list[~disease_list["Phenotype_Start"].isin(["?", "[", "{"])]

disease_list = disease_list[disease_list.moi != "digenic"]

disease_list["moi"] = ["XLD" if xi == "XL" else xi for xi in disease_list["moi"]]


# Group genes by all Modes Of Inheritance, confirm that some genes have multiple
grouped_list = disease_list.groupby("Ensembl Gene ID")["moi"].unique().to_frame()


# Counts for length of unique MOIs
grouped_list["moi_counts"] = [len(xi) for xi in grouped_list["moi"]]
# Filter to only genes with plural for later filtering
grouped_list = grouped_list[grouped_list.moi_counts > 1]


# Sorted set, for consistency and to see what counts we have
grouped_list["moi_sorted_set"] = [list(set(sorted(xi))) for xi in grouped_list.moi]


# Defined Mixed MOI for later classifications as "mixed sex-linked" or so on
grouped_list["mixed_moi"] = None

for xi, yi in grouped_list.iterrows():
    if "Digenic" in yi.moi_sorted_set:
        raise Exception("Previously filtered to no digenic...")
    elif "XLR" in yi.moi_sorted_set:
        grouped_list.loc[xi, "mixed_moi"] = "Mixed Sex-Linked"
    else:
        grouped_list.loc[xi, "mixed_moi"] = "Mixed Autosomal"


# Annotated this Mixed Phenotypes/MOIs onto the Disease Table
disease_list["mixed_phenotypes"] = [
    xi in list(grouped_list.index) for xi in disease_list["Ensembl Gene ID"]
]

disease_list["moi"] = disease_list["moi"].replace("AD/AR", "AD/AR MOI")


# Annotate properly mixed MOIs
for xi, yi in disease_list.iterrows():
    if yi["Ensembl Gene ID"] in list(grouped_list.index):
        disease_list.loc[xi, "moi"] = grouped_list.loc[
            yi["Ensembl Gene ID"], "mixed_moi"
        ]
    else:
        pass


# Only interested in those with Ensembl Gene ID
disease_list = disease_list[~disease_list["Ensembl Gene ID"].isna()]
disease_list = disease_list.set_index("Ensembl Gene ID")
# Binary to be used later - for if a gene is associated with any disease, regardless of moi
disease_list["disease_gene_binary"] = True


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

translationDict.rename(
    columns={
        "Gene stable ID": "ENSG",
        "Transcript stable ID": "ENST",
        "Protein stable ID": "ENSP",
    },
    inplace=True,
)

translationDict = translationDict[["ENSG", "ENST", "ENSP", "Ensembl_gene_length"]]

# Make different dictionary versions indexed on ENSG, ENST, and ENSP.
transEnsg = translationDict.copy()
transEnsg.set_index("ENSG", inplace=True)

transEnst = translationDict.copy()
transEnst.set_index("ENST", inplace=True)

transEnsp = translationDict.copy()
transEnsp.set_index("ENSP", inplace=True)


# Table S1A from the gene-age-3Dgenome publication, as a source to annotate data onto
# Containing information on all genes and their ages, with coordinates in hg19
s1a = pd.read_excel("/Users/marten/Downloads/Table_S1.xlsm", skiprows=4)
s1a_origin = s1a.copy()
s1a = s1a.set_index("Protein_ID").join(transEnsp)


# Filter to only relevant columns for later joining
omim_join = disease_list[
    [
        "Gene/Locus And Other Related Symbols",
        "Gene Name",
        "Approved Gene Symbol",
        "Phenotypes",
        "moi",
    ]
]

# Analysis is not phenotype-specific, only MOI
omim_join = (
    omim_join.groupby(omim_join.index)[
        [
            "Gene/Locus And Other Related Symbols",
            "Gene Name",
            "Approved Gene Symbol",
            "Phenotypes",
            "moi",
        ]
    ]
    .agg(["unique"])
    .droplevel(1, axis=1)
)

omim_join["disease_gene_inheritance"] = [str(xi[0]) for xi in omim_join["moi"]]
omim_join["disease_gene_binary"] = True


# Only join by ENSG
s1a_ENSG = s1a.reset_index(drop=False).set_index("ENSG").join(omim_join)


# Binary is True if a gene is any sort of disease gene
s1a_ENSG["disease_gene_binary"] = [xi == True for xi in s1a_ENSG["disease_gene_binary"]]


# Further annotations - descriptive notes for Non-Disease Genes

s1a_ENSG = s1a_ENSG.reset_index(drop=False)

for xi, yi in s1a_ENSG.iterrows():

    # Set default to Non-Disease, for those not joined to
    if str(yi.moi) in ["None", "NaN", "nan"]:
        s1a_ENSG.loc[xi, "disease_gene_inheritance"] = "Non-Disease"

    # If default and on Sex Chromosome, note as such!
    if (str(yi["Chr"]) in ["chrX", "chrY"]) & (
        s1a_ENSG.loc[xi, "disease_gene_inheritance"] == "Non-Disease"
    ):
        s1a_ENSG.loc[xi, "disease_gene_inheritance"] = "Non-Disease Sex-Linked"

    # For ALL Unannotated Genes, note as Excluded
    if str(yi["Protein_ID"])[:4] != "ENSP":
        s1a_ENSG.loc[xi, "disease_gene_inheritance"] = "EXCLUDED-Unannotated"

s1a_ENSG = s1a_ENSG.set_index("Protein_ID")


s1a_ENSG = s1a_ENSG.drop("moi", axis=1)


# Table from the gene-age-3Dgenome publication, as a source to annotate data onto
# Containing expression information in 6 tissue categories, SeqSim age, and Synteny age
# With coordinates in GRCh38
s3a = pd.read_csv(
    "~/Serious_Work/Batch_Scripts/unannotated-genes-with-all/Tables S3/tables_3a_marten_synteny_20241007.tsv",
    sep="\t",
)
s3a = s3a.set_index("Identifier").join(transEnsp)


# Index by ENSG, join to OMIM table
s3a_ENSG = (
    s3a.reset_index(drop=False)
    .set_index("ENSG")
    .join(omim_join)
    .reset_index(drop=False)
    .set_index("Identifier")
)

s3a_ENSG["disease_gene_binary"] = [xi == True for xi in s3a_ENSG["disease_gene_binary"]]

# Add details about non-disease genes as well

s3a_ENSG_q = s3a_ENSG.reset_index(drop=False)

for xi, yi in s3a_ENSG_q.iterrows():

    #     s3a_ENSG_q[aguwrin]

    # Set default to Non-Disease
    if str(yi.moi) in ["None", "NaN", "nan"]:
        s3a_ENSG_q.loc[xi, "disease_gene_inheritance"] = "Non-Disease"

    # If default and on Sex Chromosome, note as such!
    if (str(yi["chr"]) in ["chrX", "chrY"]) & (
        s3a_ENSG_q.loc[xi, "disease_gene_inheritance"] == "Non-Disease"
    ):
        s3a_ENSG_q.loc[xi, "disease_gene_inheritance"] = "Non-Disease Sex-Linked"

    # For ALL Unannotated Genes, note as Excluded
    if "Unannotated" in str(yi["annotation"]):
        s3a_ENSG_q.loc[xi, "disease_gene_inheritance"] = "EXCLUDED-Unannotated"

    # For ALL Intergenic Controls, note as Excluded
    if str(yi["annotation"]) in ["orf", "norf"]:
        s3a_ENSG_q.loc[xi, "disease_gene_inheritance"] = "EXCLUDED-Intergenic Control"


s3a_ENSG_q = s3a_ENSG_q.drop("moi", axis=1).set_index("Identifier")


# Excluded & Plotted Logic
# For future decisions to not plot Mixed MOIs


s3a_ENSG_q["Excluded"] = [
    "EXCLUDED" in xi for xi in s3a_ENSG_q.disease_gene_inheritance
]

s3a_ENSG_q["GRCh38_Expression_Plotting"] = [
    not any(["AD/AR" in xi, "EXCLUDED" in xi, "Mixed" in xi])
    for xi in s3a_ENSG_q.disease_gene_inheritance
]


s1a_ENSG["Excluded"] = ["EXCLUDED" in xi for xi in s1a_ENSG.disease_gene_inheritance]

s1a_ENSG["GRCh38_Expression_Plotting"] = [
    not any(["AD/AR" in xi, "EXCLUDED" in xi, "Mixed" in xi])
    for xi in s1a_ENSG.disease_gene_inheritance
]


s3a_ENSG_q.to_csv("marten_s3a_omimdisease_20250314.tsv", sep="\t")
s1a_ENSG.to_csv("marten_s1a_omimdisease_20250314.tsv", sep="\t")
