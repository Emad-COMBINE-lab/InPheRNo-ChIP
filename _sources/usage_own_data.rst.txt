Using InPheRNo-ChIP with Your Own Data
======================================

InPheRNo-ChIP is a computational framework designed to reconstruct phenotype-relevant gene regulatory networks (GRNs) by integrating RNA-seq data, TF-specific ChIP-seq data, and phenotypic labels. 
This guide provides detailed instructions on how to prepare your own data to use with InPheRNo-ChIP, outlining necessary preprocessing steps, tools used, and considerations for optional steps.
We will cover:

- **Necessary data inputs** and their required formats.
- **Data preprocessing steps** for each data type.
- **Optional steps** that may enhance your analysis.
- **Tools and methods** used in our analysis, which you can adopt or modify based on your needs.

.. contents::
   :local:
   :depth: 2


Overview of Required Data In This Guide
---------------------------------------

You will need the following data:

1. **RNA-seq Data**: A gene expression count matrix (genes × samples) containing expression values for both target genes and transcription factors (TFs) across samples with known phenotypic labels.

2. **TF-specific ChIP-seq Data**: ChIP-seq peak files in BED format for the TFs of interest in your study.

3. **Phenotypic Labels**: Labels indicating the phenotype of each sample (e.g., disease status, cell type).

Preparing RNA-seq Data
----------------------
In this section, we will explain how to prepare the RNA-seq data required for InPheRNo-ChIP, specifically focusing on generating the three essential inputs for InPheC_Step1.py, as shown in our GitHub repo's `README.md <https://github.com/Emad-COMBINE-lab/InPheRNo-ChIP/tree/main#readme>`_:

- **Input 1.1**: Normalized counts (e.g., ``RNAseq_expr_voom_normalized.csv``)
- **Input 1.2**: P-values of gene-phenotype associations for significant genes (e.g., ``RNAseq_pval_gene-phenotype_significant.csv``)
- **Input 1.3**: P-values of gene-phenotype associations for all genes (e.g., ``RNAseq_pval_gene-phenotype_allGenes.csv``)

We will walk through the steps needed to prepare these inputs from your RNA-seq and phenotype data.

**Approach**:

1. **Data Acquisition**:

   - If you have raw RNA-seq reads, process them to generate a gene expression count matrix (genes × samples) using standard RNA-seq pipelines.
     
      - Alignment: Use aligners like **STAR** or **HISAT2** to map reads to the reference genome.
      - Quantification: Use tools like **featureCounts** or **HTSeq-count** to obtain raw counts.
   
   - Alternatively, if you already have a gene expression count matrix, you can proceed directly.

2. **Sample Selection and Quality Control**:

   - Convert all gene identifiers to a consistent format (e.g., HGNC symbols) to ensure compatibility across datasets.
   - Assess the quality of your gene expression data.
   - Filter out low-quality samples or genes with low expression levels as appropriate for your analysis.
   - **Optional**: Remove genes located on sex chromosomes if sex-specific effects are not relevant to your study.

3. **Normalization and Batch Effect Removal**:

   - **Normalization**:

     - Normalize the count matrix to account for library size differences and other technical biases.
     - In our analysis, we used the **TMM** (Trimmed Mean of M-values) normalization method implemented in **edgeR**'s `calcNormFactors` function.
     - Alternative normalization methods include **DESeq2**'s median-of-ratios method or **limma**'s `voom` normalization.

   - **Batch Effect Removal**:

     - If your data comes from multiple sources or experiments, consider removing batch effects to minimize unwanted variability.
     - In our analysis, we used **ComBat-seq** from the **sva** package to adjust for batch effects.
     - Alternative tools include **limma**'s `removeBatchEffect` function or **RUVSeq**.

4. **Transformation for Linear Modeling**:

   - Transform the normalized count data to be suitable for linear modeling.
   - In our analysis, we applied the **voom** transformation using **limma**'s `voom` function with quantile normalization. 
   - We then transformed each gene's expression values across samples to follow a standard normal distribution (mean zero, variance one), ensuring comparability across genes.

5. **Differential Expression Analysis**:

   - Perform differential expression (DE) analysis to identify (non-TF) genes associated with your phenotype.
   - In our analysis, we used **edgeR** to perform DE analysis, incorporating factors like batch and lineage in the design matrix.
   - Alternative tools include **DESeq2**.

6. **Prepare Input Files**:

   - Save the voom-transformed and inverse-quantile normalized expression data from step 4 here to a CSV file (e.g., input 1.1 in the GitHub repo's `README.md <https://github.com/Emad-COMBINE-lab/InPheRNo-ChIP/tree/main#readme>`_): ``RNAseq_expr_voom_normalized.csv``).
   - Save the DE results for significant genes from step 5 here in a CSV file (e.g., input 1.2: ``RNAseq_pval_gene-phenotype_significant.csv``).
   - Save the DE results for all genes (not just significant ones) from step 5 here in a CSV file  (e.g., input 1.3: ``RNAseq_pval_gene-phenotype_allGenes.csv``)
   - All three files should be placed in ``./Data/RNA_seq/``

**Notes**:

- **Confounding Variables**: Include factors such as batch, gender, or other experimental conditions in your design matrix to control for their effects during DE analysis.
- **TFs in Expression Data**: Ensure you have expression data for the TFs you plan to include in your GRN.

Preparing ChIP-seq Data
-----------------------

In this section, we will explain how to prepare ChIP-seq data.

**Approach**:

1. **Data Acquisition**:

   - Obtain ChIP-seq peak files (BED format) for your TFs of interest.
   - If you have raw ChIP-seq reads, process them to generate peak files using standard ChIP-seq analysis pipelines. 
   - Tools for processing raw reads include aligners like **Bowtie2** or **BWA** and peak-calling software such as **MACS2**.
   - Alternatively, use publicly available ChIP-seq peak datasets relevant to your study.

2. **Sample Filtering**:

   - Convert all gene identifiers to a consistent format (e.g., HGNC symbols) to ensure compatibility across datasets.
   - Remove low-quality samples and those not matching your criteria.
   - In our analysis, we filtered out samples with low quality or incorrect annotations based on metadata review.

3. **Combining Replicates**:

   - If you have biological replicates, combine peaks from replicates to identify reproducible binding sites.
   - In our analysis, we used the **IDR (Irreproducible Discovery Rate)** framework from `ENCODE <https://www.encodeproject.org/software/idr/>`_ to combine peaks from replicates and select reproducible binding sites.

4. **Filtering Peaks**:

   - Remove peaks in blacklisted regions to avoid artifacts and regions prone to mapping errors.
   - In our analysis, we used **BEDTools** to exclude blacklisted regions, such as the ENCODE blacklist regions.

5. **Assigning Peaks to Genes**:

   - Assign peaks to genes based on proximity to transcription start sites (TSS) or other genomic features.
   - In our analysis, we used **T-Gene** from MEME suite (`link <https://meme-suite.org/meme/doc/tgene.html?man_type=web>`_) to associate peaks with potential target genes.

6. **Prepare Input Files**:

   - For each TF, extract p-values or scores for each TF-gene pair from the peak assignment results and output to a .BED file.
   - Place processed BED files in the ``./Data/ChIP_GSExxx`` directory (replace GSExxx with the appropriate identifier).


