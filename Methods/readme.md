**Code here was used to generate results for later analysis and to prepare files for processing.**

**GRCh38**

Step 1 - Intergenic control BAM files were delivered to us by collaborator Kee-Myoung Nam and downloaded for GRCh38 and split into 5 control groups. The file 'GRCh38 Step 1 - Deduplication & Construction.ipynb' shows the steps to de-duplicate these controls and check for overlap against annotated genes, unannotated genes, and other control sequences. The results were written out.

Step 2 - Controls were removed which did not map over cleanly to GRCh37 chromosomes, as per 'GRCh38 Step 2 - BED Construction with 37_Crossover Failure Removal Check.ipynb'. Additional sequence data is annotated to the tables here as well. 

Step 2+1 - Sufficiently filtered BED files were run through 'bed2gtf.py' and then 'clean_gtf.py' , and then passed into a Terra workspace, where our GTEx RNA-Seq Pipeline was run to calculate the mean RNA expression of each gene or control.

Step 3 - The outputs of the GTEx RNA-Seq pipeline, containing the mean read counts per tissue sample, were first copied out of the Terra workspace and into a GCP bucket for the project, then read into a notebook, normalized via DESeq2’s median of ratios to correct for RNA composition and sequencing depth, and then mean scores calculated per tissue type. Though RNA-Seq calling was done for all sets of intergenic controls, only one set of controls was included in the analysis, and therefore only one was included in normalization. Done in 'GRCh38 Step 3 - Merge Samples and Normalization and Mean Count Outputs.ipynb'. The outputs of this were sent on to /Analysis.

NOTE 1 - Additional unannotated genes were removed after GTEx RNA-Seq calling but before any meaning or normalization steps. These 47 inconsistent unannotated genes are removed in the prior 'Step 3' notebook.

NOTE 2 - Internally and briefly, the five intergenic control sets (split by Igen ORF and Igen Non-ORF) were compared against each other to ensure that expression wasn't discordant (by Igen ORF or Igen Non-ORF type) between control sets.

**GRCh37**

NOTE - To ensure that GRCh37 controls were not discordant with GRCh38 controls, they were also run through the GTEx RNA-Seq Pipeline. The pipeline and workspace is for GTEx v8 in GRCh38, so for these checks the GRCh37 controls were lifted over to GRCh38.

Step 1 - In 'GRCh37 Step 1 - Crossover_38 Deduplication.ipynb', controls were lifted over and similarly to GRCh38 controls checked for overlap/duplication internally and against the exact same set of GRCh38 annotated and unannotated genes. This process also removed any controls that did not cleanly map over to canonical GRCh38 chromosomes. This file produced the GRCh37->38 controls to be passed into the GTEx RNA-Seq pipeline to compare against GRCh38.

Step 2 - In 'GRCh37 Step 2 - Output Construction.ipynb', tables of the GRCh37 Controls (in GRCh37 coordinates) that passed the filtering are produced. Additional sequence data is annotated here as well. 

2+1 - These controls were run through the GTEx RNA-Seq pipeline, normalized in the same way as above, and found to not be discordant with GRCh38 controls.

**Gene Ontology**

3D_Genome_Overlap.py used to generate input lists for gene ontology analysis.
