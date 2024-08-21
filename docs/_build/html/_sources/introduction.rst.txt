Introduction
=============================

InPheRNo-ChIP: Inference of Phenotype-relevant Regulatory Networks with ChIP-seq Integration
--------------------------------------------------------------------------------------------

InPheRNo-ChIP is a computational framework designed to reconstruct phenotype-relevant gene regulatory networks (GRNs) by integrating RNA-seq, TF-specific ChIP-seq, and phenotypic labels. This tool employs a probabilistic graphical model to model the simultaneous effects of transcription factors on target genes, facilitating the identification of regulatory interactions that are crucial for specific phenotypic outcomes. By leveraging multi-modal data, InPheRNo-ChIP enhances the accuracy of GRN predictions, making it particularly useful in developmental biology and disease modeling.

Framework Overview
------------------

**InPheRNo-ChIP** main steps, details can be found in the accompanying research paper:

1. **RNA-seq Data Processing**: Estimating p-values for TF-gene interactions and gene-phenotype associations using three RNA-seq datasets. 
2. **ChIP-seq Data Processing**: Deriving p-values for TF-gene interactions from ChIP-seq data.
3. **PGM Integration**: Merging calculated p-values within a PGM, introducing binary variables to represent TF-gene pairs.
4. **MCMC Sampling**: Estimating posterior probabilities for these variables to form an initial regulatory graph.
5. **Normalization and Filtering**: Refining the initial graph to produce a precise and phenotype-relevant GRN.

Key Features
------------

- **Comprehensive Data Integration**: **InPheRNo-ChIP** extends the foundational InPheRNo algorithm (available on `GitHub <https://github.com/KnowEnG/InPheRNo>`_) by integrating various datasets, including RNA-seq and ChIP-seq data, from human embryonic and endoderm cell lines.
- **Probabilistic Graphical Model (PGM)**: The framework employs a PGM to systematically model the influence of transcription factors on target genes.
- **Advanced Gene-Filtering and Normalization**: Incorporating advanced techniques, **InPheRNo-ChIP** refines gene regulatory network (GRN) inference, enhancing the accuracy and relevance of its findings.
- **Phenotype-Relevant GRN Inference**: The tool is adept at constructing GRNs that are directly relevant to phenotypic variations, particularly focusing on key endoderm markers.

Intended Audience
-----------------

**InPheRNo-ChIP** is tailored for:

- *Developmental Biologists*: Researchers studying embryogenesis, cellular differentiation, and developmental processes will find InPheRNo-ChIP particularly useful for understanding gene regulation during these critical stages.
- *Bioinformaticians and Computational Biologists*: Professionals who specialize in analyzing biological data, especially those with a focus on high-throughput sequencing data, will utilize InPheRNo-ChIP for in-depth analysis and interpretation.
- and much more!
