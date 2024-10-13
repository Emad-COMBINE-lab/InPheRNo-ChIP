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


Overview of Required Data
-------------------------

To run InPheRNo-ChIP, you will need the following:

1. **RNA-seq Data**: A gene expression count matrix (genes × samples) containing normalized expression values for both target genes and transcription factors (TFs) across samples with known phenotypic labels.

2. **TF-specific ChIP-seq Data**: ChIP-seq peak files in BED format for the TFs of interest in your study.

3. **Phenotypic Labels**: Labels indicating the phenotype of each sample (e.g., disease status, cell type).

Preparing RNA-seq Data
----------------------


**Requirement**: A gene expression count matrix (genes × samples) containing normalized expression values.

**Approach**:

1. **Data Acquisition**:

   - Obtain a gene expression count matrix for your samples. This can be from your own experiments or publicly available datasets.
   - If you have raw RNA-seq reads, process them to generate the count matrix using standard RNA-seq pipelines. This typically involves read alignment and quantification, but in this guide we expect the input to be a count matrix.
   - Tools for processing raw reads include aligners like **STAR** or **HISAT2** and quantification tools like **featureCounts** or **HTSeq-count**.

2. **Quality Control**:

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

4. **Differential Expression Analysis**:

   - Perform differential expression (DE) analysis to identify genes associated with your phenotype.
   - In our analysis, we used **edgeR** to perform DE analysis, incorporating factors like batch and lineage in the design matrix.
   - Alternative tools include **DESeq2** and **limma-voom**.
   - **Optional**: Exclude transcription factor genes from the list of differentially expressed genes if they will serve as regulators in your GRN.

5. **Transformation for Linear Modeling**:

   - Transform the normalized count data to be suitable for linear modeling.
   - We applied the **voom** transformation using **limma**'s `voom` function with quantile normalization.
   - This transformation estimates the mean-variance relationship and provides weights for linear modeling.

6. **Prepare Input Files**:

   - Save the list of differentially expressed genes and their associated p-values.
   - Ensure that gene identifiers are consistent across all datasets (e.g., using HGNC symbols).

**Notes**:

- **Confounding Variables**: Include factors such as batch, gender, or other experimental conditions in your design matrix to control for their effects during DE analysis.
- **TFs in Expression Data**: Ensure you have expression data for the TFs you plan to include in your GRN.

Preparing ChIP-seq Data
-----------------------

**Requirement**: ChIP-seq peak files in BED format for your TFs of interest.

**Approach**:

1. **Data Acquisition**:

   - Obtain ChIP-seq peak files (BED format) for your TFs of interest.
   - If you have raw ChIP-seq reads, process them to generate peak files using standard ChIP-seq analysis pipelines. 
   - Tools for processing raw reads include aligners like **Bowtie2** or **BWA** and peak-calling software such as **MACS2**.
   - Alternatively, use publicly available ChIP-seq peak datasets relevant to your study.

2. **Sample Filtering**:

   - Map TFs in your sample metadata to universal gene names to ensure consistency across datasets.
   - Remove low-quality samples and those not matching your criteria.
   - In our analysis, we filtered out samples with low quality or incorrect annotations based on metadata review.

3. **Combining Replicates**:

   - If you have biological replicates, combine peaks from replicates to identify reproducible binding sites.
   - In our analysis, we used the **IDR (Irreproducible Discovery Rate)** framework to combine peaks from replicates and select reproducible binding sites.
   - Alternative tools include **PePr**, **SICER**, or the **ENCODE Consistency Model**.

4. **Filtering Peaks**:

   - Remove peaks in blacklisted regions to avoid artifacts and regions prone to mapping errors.
   - In our analysis, we used **BEDTools** to exclude blacklisted regions, such as the ENCODE blacklist regions.

5. **Assigning Peaks to Genes**:

   - Assign peaks to genes based on proximity to transcription start sites (TSS) or other genomic features.
   - In our analysis, we used **T-Gene** to associate peaks with potential target genes.

6. **Prepare Input Files**:

   - Extract p-values or scores for each TF-gene pair from the peak assignment results.
   - Ensure that gene identifiers are consistent with those used in other data types.


