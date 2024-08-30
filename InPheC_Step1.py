"""
This script is part of the InPheRNo-ChIP project and is responsible for process and analyze RNA-seq and ChIP-seq data. 
It includes functionality for:

    - Parsing command-line arguments to configure paths and options for data processing.
    - Loading ChIP-seq and RNA-seq data from specified directories: ./Data/
    - Processing and filtering data to ensure consistency among gene and TF sets.
    - Performing Elastic Net regression on gene and TF expression data.

Usage:

    - Run the script from the command line "python InPheC_Step1.py" directly.
    - For detailed argument options, use the help option: `python InPheC_Step1.py --help`

Note: Ensure that all dependencies are installed and the Python environment is correctly set up for running this script.
"""

import argparse
import glob
import logging
import os
import sys
import warnings

import numpy as np
import pandas as pd
import scipy.stats as ss
import statsmodels.api as sm
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import ElasticNetCV
from tqdm import tqdm

from utils_logo import print_logo


class Step1ArgumentParser:
    """
    Handles command-line arguments for the InPheRNo-ChIP project, setting up and parsing them.
    
    :param argparse.ArgumentParser parser: Configures the command-line arguments.
    :param argparse.Namespace args: Stores the values of parsed arguments.
    """

    def __init__(self):
        """Initializes the ArgumentParser with a description of the script's purpose and setups the arguments."""
        self.parser = argparse.ArgumentParser(
            description="Configure paths and options for RNA-seq and ChIP-seq data processing."
        )
        self._add_arguments()

    def _add_arguments(self):
        """
        Adds command-line options for the script, specifying directories, files, and other parameters
        required for processing RNA-seq and ChIP-seq data.
        """
        # Directories
        self.parser.add_argument(
            "-id_r",
            "--in_dir_rna",
            default="./Data/RNA_seq",
            help="Path to GEx data directory.",
        )
        self.parser.add_argument(
            "-id_c",
            "--in_dir_chip",
            default="./Data/ChIP_GSE61475_G",
            help="Path to ChIP-seq data directory.",
        )
        self.parser.add_argument(
            "-od",
            "--out_dir",
            default="./step1/",
            help="Path to output directory for step 1.",
        )

        # Input Files
        self.parser.add_argument(
            "-it",
            "--input_tf_addr",
            default="./Data/mapped_hmTFs_legitTFs-HumanTFDB.csv",
            help="CSV file containing list of tfs.",
        )
        self.parser.add_argument(
            "-ie",
            "--input_expression",
            default="RNAseq_expr_sample_gse981-gse361-gse371_voom-zScore.csv",
            help="CSV file with gene and tf expression data.",
        )
        self.parser.add_argument(
            "-igp",
            "--input_gene_phenotype_interest",
            default="RNAseq_pval_gene-phenotype_gse361-gse981-gse371_logFC2_FDR0.01.csv",
            help="CSV file with p-values of gene-phenotype.",
        )

        # Options and Parameters
        self.parser.add_argument(
            "-loi",
            "--lineage_of_interest",
            nargs="+",
            default=["h64", "dEN"],
            help="Prefix for lineages of interest in ChIP-seq filenames.",
        )
        self.parser.add_argument(
            "-cpoi",
            "--chip_pvalue_of_interest",
            default="Distance_P_Value",
            help="P-value of interest, this corresponds to the t-Gene output column names.",
        )

        self.parser.add_argument(
            "-num_coefs",
            "--max_num_coefs",
            default=22,
            type=int,
            help="maximum number of TFs to be selected by Elastic Net, \
                we set to 22 bc all 22 TFs are relevant to the differential process. \
                This number allows the model to consider all TFs while still guarding against overfitting, \
                given the sample size of 28",
        )

        # Output Files
        self.parser.add_argument(
            "-ogp",
            "--output_gene_phenotype",
            default="RNAseq_pval_gene-pheno_intersect_chip_",
            help="Output file with gene-phenotype p-values.",
        )
        self.parser.add_argument(
            "-tgt",
            "--output_gene_tf",
            default="RNAseq_pval_gene-tf_intersect_chip_",
            help="Output file with gene-tf p-values.",
        )
        self.parser.add_argument(
            "-tgpk",
            "--output_chip_tf_gene_peaks",
            default="ChIPseq_pval_gene-tf_minQ_",
            help="Output file prefix with ChIP-seq cleaned data, which contains the minimum p-value across peaks for each lineage",
        )

    def parse_args(self):
        """Parses the command-line arguments."""
        return self.parser.parse_args()


class DataLoader:
    """
    Loads and prepares RNA-seq and ChIP-seq data from specified directories.

    :param chip_seq_dir: Path to the ChIP-seq data directory.
    :type chip_seq_dir: str
    :param rna_seq_dir: Path to the RNA-seq data directory.
    :type rna_seq_dir: str
    """

    def __init__(self, chip_seq_dir, rna_seq_dir):
        """
        Initializes the DataLoader with specified directories for ChIP-seq and RNA-seq data.
        """
        self.chip_seq_dir = chip_seq_dir
        self.rna_seq_dir = rna_seq_dir

    def load_chip_data(self, lineage_of_interest):
        """
        Loads and merges ChIP-seq data files for specified lineages.

        :param lineage_of_interest: Lineages to include in the analysis.
        :type lineage_of_interest: list
        :return: Combined DataFrame of ChIP-seq data across specified lineages.
        :rtype: pandas.DataFrame
        """
        combined_df = pd.DataFrame()
        for lineage in lineage_of_interest:
            file_nms = sorted(glob.glob(os.path.join(self.chip_seq_dir, f"{lineage}*.bed.tsv")))
            for file_nm in file_nms:
                df = pd.read_csv(file_nm, sep="\t")
                df["lineage"] = lineage
                tf_name = os.path.basename(file_nm).split("_")[1]  # Extracting TF name (second element after split)
                df["TF"] = tf_name
                combined_df = pd.concat([combined_df, df], ignore_index=True)
        return combined_df

    def load_rna_data(self, f_exp, f_gp):
        """
        Loads RNA-seq expression data and associated gene-phenotype p-values from specified files.

        :param f_exp: Filename of the RNA-seq expression data file.
        :type f_exp: str
        :param f_gp: Filename of the gene-phenotype p-values data file.
        :type f_gp: str
        :return: Tuple containing expression data and gene-phenotype p-values.
        :rtype: (pandas.DataFrame, pandas.DataFrame)
        """
        expr = pd.read_csv(
            os.path.join(self.rna_seq_dir, f_exp), index_col=0, delimiter="," if f_exp.endswith(".csv") else "\t"
        )
        gene_pheno_pval = pd.read_csv(
            os.path.join(self.rna_seq_dir, f_gp), index_col=0, delimiter="," if f_gp.endswith(".csv") else "\t"
        )
        return expr, gene_pheno_pval

    def load_tf_list(self, f_tf):
        """
        Loads a list of transcription factors from a specified CSV file.

        :param f_tf: Path to the CSV file containing the list of transcription factors.
        :type f_tf: str
        :return: DataFrame containing the list of transcription factors.
        :rtype: pandas.DataFrame
        """
        tf_list = pd.read_csv(f_tf)
        return tf_list


class DataProcessor:
    """
    Processes and filters RNA-seq and ChIP-seq data to ensure consistency among gene and TF sets.
    This class is crucial for preparing the data for subsequent analysis, including filtering based on shared genes and TFs.
    """

    def __init__(self, chip_data, GEx_gene_tf, gene_phenotype, hm_tf_df, chip_pvalue_of_interest):
        """
        Initializes the DataProcessor with various data sources for processing.

        :param chip_data: ChIP-seq data as a pandas DataFrame.
        :type chip_data: pandas.DataFrame
        :param GEx_gene_tf: Gene and TF expression data as a pandas DataFrame.
        :type GEx_gene_tf: pandas.DataFrame
        :param gene_phenotype: Gene-phenotype data as a pandas DataFrame.
        :type gene_phenotype: pandas.DataFrame
        :param hm_tf_df: Human transcription factors data as a pandas DataFrame.
        :type hm_tf_df: pandas.DataFrame
        :param chip_pval_flag: Column name in ChIP-seq data for the p-value of interest.
        :type chip_pval_flag: str
        """
        self.chip_data = chip_data
        self.GEx_gene_tf = GEx_gene_tf
        self.gene_phenotype = gene_phenotype
        self.hm_tf_list = set(hm_tf_df.Symbol.tolist())
        self.gex_genes_tfs = set(list(self.GEx_gene_tf.index.values))
        self.chip_pval_flag = chip_pvalue_of_interest

    def _find_common_genes(self):  # for DE genes, ChIP genes, and GEx genes
        """
        Identifies common genes across RNA-seq and ChIP-seq datasets.

        :return: A set of gene IDs common across all datasets.
        :rtype: set
        """
        chip_genes = set(self.chip_data["Gene_ID"])  # 14982 genes
        de_genes = set(self.gene_phenotype.index.values)  # 1745 DE genes
        intersect = chip_genes & de_genes & self.gex_genes_tfs  # 1624 DE-ChIP-voom genes
        print(
            f" >> Originally, # DE genes: {len(de_genes)}, # ChIP genes: {len(chip_genes)}, # GEx genes: {len(self.gex_genes_tfs)}"
        )
        print(f" >> after intersection, # DE-ChIP genes: {len(intersect)}")
        return intersect

    def _find_common_tfs(self):
        """
        Identifies common transcription factors across RNA-seq and ChIP-seq datasets.

        :return: A set of TF identifiers common to all datasets.
        :rtype: set
        """
        chip_tfs = set(self.chip_data["TF"])
        return chip_tfs & self.hm_tf_list & self.gex_genes_tfs

    def filter_data(self):
        """
        Filters RNA-seq and ChIP-seq data to ensure consistency among gene and transcription factor (TF) sets,
        and prepares the data for subsequent analysis steps. This method processes the data by:
        - Identifying common genes and TFs across different datasets.
        - Filtering ChIP-seq data to keep rows that have genes and TFs present in both lineages.
        - Filtering to keep the minimal p-value across peaks for each (TF, gene, lineage) combination.
        - Preparing expression data matrices for genes and TFs for Elastic Net regression.

        :return: A tuple containing:
            - pandas.DataFrame: Filtered ChIP-seq data with minimal p-values across peaks, reset index.
            - pandas.DataFrame: Gene expression data for shared genes.
            - pandas.DataFrame: TF expression data for shared TFs.
            - pandas.DataFrame: Gene-phenotype data for shared genes.
        :rtype: tuple
        """
        shared_genes = self._find_common_genes()  # 1624 DE-ChIP-Pheno genes
        shared_tfs = self._find_common_tfs()

        print(f" Now filtering ChIP-seq data..")
        # filter 1: filter chip to keep rows that have gene on both lineages
        tmp = self.chip_data[self.chip_data["Gene_ID"].isin(shared_genes) & self.chip_data["TF"].isin(shared_tfs)]
        # filter 2: filter chip to keep rows that have TF on both lineages
        filtered_chip = tmp.groupby(["Gene_ID", "TF"]).filter(
            lambda x: x["lineage"].nunique() > 1
        )  # 1405 DE-ChIP-voom genes

        # filter 3: filter chip to keep (TF, gene, lineage) with lowest p-value across PEAKS
        filtered_chip_lowest_pval_across_peaks = (
            filtered_chip.sort_values(self.chip_pval_flag).groupby(["TF", "Gene_ID", "lineage"]).head(1)
        )

        chip_genes = set(filtered_chip_lowest_pval_across_peaks["Gene_ID"].tolist())
        chip_tfs = set(filtered_chip_lowest_pval_across_peaks["TF"].tolist())

        # [GEX data] split to gene and tf expression, prepare for Elastic Net
        expr_gene = self.GEx_gene_tf.loc[list(chip_genes)]  # (1405, 28)
        expr_tf = self.GEx_gene_tf.loc[list(chip_tfs)]  # (22, 28)
        filtered_gp = self.gene_phenotype.loc[list(chip_genes)]  # (1405, 5)

        # find shared genes and tfs after filtering
        refined_shared_genes = set(expr_gene.index.tolist()) & set(chip_genes)
        refined_shared_tfs = set(expr_tf.index.tolist()) & set(chip_tfs)

        print(f" >> (refined) # shared genes in ChIP-seq, DE, and GEx: {len(refined_shared_genes)}")
        print(f" >> (refined) # shared TFs in ChIP-seq, DE, and GEx: {len(refined_shared_tfs)}")

        # assert
        assert expr_gene.shape[0] == len(chip_genes), "Number of genes in expr_gene does not match chip_genes"
        assert expr_tf.shape[0] == len(chip_tfs), "Number of genes in expr_tf does not match chip_tfs"

        return filtered_chip_lowest_pval_across_peaks.reset_index(drop=True), expr_gene, expr_tf, filtered_gp


class ElasticNetProcessor:
    """
    Conducts Elastic Net regression to identify relationships between gene expressions and transcription factors
    within the InPheRNo-ChIP project. It dynamically adjusts model parameters to optimize fit and manage convergence.

    :param expr_gene: DataFrame containing gene expression data.
    :type expr_gene: pandas.DataFrame
    :param expr_tf: DataFrame containing TF expression data.
    :type expr_tf: pandas.DataFrame
    :param max_num_coefs: Maximum allowed non-zero coefficients in the model, provides a control over model complexity.
    :type max_num_coefs: int, optional
    """

    def __init__(self, expr_gene, expr_tf, max_num_coefs=None):
        """
        Initializes the processor with gene and TF expression data sets, and sets up the Elastic Net model parameters.

        :param expr_gene: Gene expression data as a pandas DataFrame.
        :param expr_tf: TF expression data as a pandas DataFrame.
        :param max_num_coefs: Optional; sets a cap on the number of non-zero coefficients to prevent overfitting.
        """
        self.expr_gene = expr_gene
        self.expr_tf = expr_tf
        self.l1_rat = 0.5
        self.eps = sys.float_info.min  # 3e-308
        self.max_num_coefs = max_num_coefs
        self.ls_genes = list(self.expr_gene.index.values)
        self.ls_tfs = list(self.expr_tf.index.values)

        assert self.max_num_coefs is not None, "max_num_coefs must be specified"

    def perform_elastic_net(self):
        """ 
        Performs Elastic Net regression across all genes, using transcription factor expressions to predict gene activity.

        :return: DataFrame with p-values indicating the significance of associations between genes and TFs.
        :rtype: pandas.DataFrame

        :note: This function is adapted from the InPheRNo paper with updated eps and n_alphas values to improve model performance and stability.
        """
        print(f"Performing Elastic Net regression on gene and TF expression data..")
        X_features = self.expr_tf.values.T
        pvalue_gt_array = (-1) * np.ones((len(self.ls_genes), len(self.ls_tfs)))
        convergence_issues = {}

        for i, gene in enumerate(self.ls_genes):
            logging.info(f"Pvalue_gene_tf: {i}")
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always", ConvergenceWarning)

                y = self.expr_gene.iloc[i].values
                EN_model = ElasticNetCV(l1_ratio=self.l1_rat, max_iter=100000)

                #### make sure that number of nonzero coefs do not exceed max_num_coefs
                alphas1, coefs1, _ = EN_model.path(
                    X_features, y, eps=0.001, n_alphas=10000
                )  # this setting gives # non-zero coefs closed to 22.
                num_coefs = np.sum(coefs1 != 0, axis=0)
                # print(num_coefs)
                # print(num_coefs[num_coefs <= max_num_coefs][-1])
                rep_EN = 0
                if num_coefs[-1] < self.max_num_coefs:
                    EN_coef = coefs1[:, -1]
                    selected_ind = np.array(range(len(self.ls_tfs)))[EN_coef != 0]
                else:
                    while (
                        (num_coefs[0] != num_coefs[-1])
                        and (max(num_coefs[num_coefs <= self.max_num_coefs]) != self.max_num_coefs)
                        and (rep_EN < 10)
                    ):
                        rep_EN += 1
                        alpha_min = alphas1[(num_coefs <= self.max_num_coefs)][-1]
                        alpha_max = alphas1[(num_coefs > self.max_num_coefs)][0]
                        alphas3 = np.linspace(alpha_min, alpha_max, 10)
                        alphas1, coefs1, _ = EN_model.path(X_features, y, alphas=alphas3)
                        num_coefs = np.sum(coefs1 != 0, axis=0)
                        # print(num_coefs)
                        # print(num_coefs[num_coefs <= max_num_coefs][-1])
                    print("repeat", rep_EN)
                    if num_coefs[0] == num_coefs[-1]:
                        EN_coef = coefs1[:, 0]
                        selected_ind = np.array(range(len(self.ls_tfs)))[EN_coef != 0]
                    else:
                        EN_coef = coefs1[:, len(num_coefs[num_coefs <= self.max_num_coefs]) - 1]
                        selected_ind = np.array(range(len(self.ls_tfs)))[EN_coef != 0]

                print(len(selected_ind))
                # n_selected_tf = len(selected_ind)
                X_features_new = X_features[:, selected_ind]

                model = sm.OLS(y, X_features_new)
                results = model.fit()

                ts_b = results.tvalues
                if any(item.category == ConvergenceWarning for item in w):
                    convergence_issues[gene] = 1  # Set to 1 to indicate a warning for this gene

            #####prcise estimation of p-vals
            N = np.shape(X_features)[0] - 1  # degrees of freedom in ttest for regression
            pvals_precise = []
            x1_bar = 0  # mean of the first distribution
            n1 = (N + 2) // 2  # num obs first distribution
            n2 = N + 2 - n1  # num obs second distribution
            x2_std = 0
            t_precise = []
            for j in range(len(ts_b)):
                T = ts_b[j]  # T statistic
                x2_bar = -np.sign(T)  # mean of the second distribution
                x1_std = 1 / T * np.sqrt((n1 + n2 - 2) / (n1 - 1) / (1 / n1 + 1 / n2))  # std of first distribution
                (t, p) = ss.ttest_ind_from_stats(x1_bar, x1_std, n1, x2_bar, x2_std, n2)
                if p < self.eps:
                    p = self.eps
                pvals_precise.append(p)
                t_precise.append(t)

            pvalue_gt_array[i, selected_ind] = pvals_precise

        pvalue_gt_df = pd.DataFrame(
            pvalue_gt_array,
            index=self.ls_genes,
            columns=self.ls_tfs,
        )
        return pvalue_gt_df


def save_output_files(args, chip_gt, rna_gt, rna_gp):
    """
    Saves processed data to specified output files as part of the InPheRNo-ChIP project. This function handles the
    creation of output directories and manages the naming conventions for output files based on the provided command-line arguments.

    :param args: Parsed command-line arguments containing output file paths and names.
    :type args: argparse.Namespace
    :param chip_gt: Processed ChIP-seq data to be saved.
    :type chip_gt: pandas.DataFrame
    :param rna_gt: Processed RNA-seq gene-TF data to be saved.
    :type rna_gt: pandas.DataFrame
    :param rna_gp: Processed RNA-seq gene-phenotype data to be saved.
    :type rna_gp: pandas.DataFrame
    """
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)
    addr_out_pvalue_gp = os.path.join(
        args.out_dir,
        f"{args.output_gene_phenotype}{'_'.join(args.lineage_of_interest)}.csv",
    )
    rna_gp.to_csv(addr_out_pvalue_gp, index=True)
    addr_out_pvalue_gt = os.path.join(
        args.out_dir,
        f"{args.output_gene_tf}{'_'.join(args.lineage_of_interest)}.csv",
    )
    rna_gt.to_csv(addr_out_pvalue_gt, index=True)
    addr_out_chip_gt = os.path.join(
        args.out_dir,
        f"{args.output_chip_tf_gene_peaks}{'_'.join(args.lineage_of_interest)}.csv",
    )
    chip_gt.to_csv(addr_out_chip_gt, index=False)
    print(f"Saved 3 output files to {args.out_dir}.")


def step1_main():
    """
    Main function for the first step in the InPheRNo-ChIP project pipeline. Orchestrates the execution of data loading,
    processing, Elastic Net regression, and saving the output files. This function acts as the entry point when running the script.
    
    This function follows these steps:
    - 1. Parses command-line arguments.
    - 2. Loads data from specified directories.
    - 3. Processes data to align gene and TF datasets.
    - 4. Performs Elastic Net regression to analyze gene-TF relationships.
    - 5. Saves the processed data to output files.
    
    :raises SystemExit: If the dimensions of gene-phenotype data and results from Elastic Net regression do not match, indicating a potential issue in data alignment or processing.
    """
    parser = Step1ArgumentParser()
    args = parser.parse_args()

    # Load data
    data_loader = DataLoader(args.in_dir_chip, args.in_dir_rna)
    chip_both_lineage = data_loader.load_chip_data(args.lineage_of_interest)
    expr_all, gene_phenotype_interest = data_loader.load_rna_data(
        args.input_expression, args.input_gene_phenotype_interest
    )
    tf_list = data_loader.load_tf_list(args.input_tf_addr)

    # Filter RNA-seq and ChIP-seq data based on common genes and TFs
    processor = DataProcessor(
        chip_both_lineage, expr_all, gene_phenotype_interest, tf_list, args.chip_pvalue_of_interest
    )
    chip_gt, tmp_expr_gene, tmp_expr_tf, rna_gp = processor.filter_data()  # 1405 DE-ChIP-Voom Genes

    # Perform Elastic Net on GENE and TF expression
    en_processor = ElasticNetProcessor(tmp_expr_gene, tmp_expr_tf, args.max_num_coefs)
    rna_gt = en_processor.perform_elastic_net()

    # Save results
    if (rna_gp.index != rna_gt.index).sum() > 0:
        sys.exit("[GEx] Dimensions do not match!")
    save_output_files(args, chip_gt, rna_gt, rna_gp)


if __name__ == "__main__":
    print_logo("Running step 1")
    step1_main()
