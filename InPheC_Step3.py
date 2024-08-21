"""
InPheC_Step3.py

This script is an integral part of the InPheRNo-ChIP project's step 3 process. 

Key Operations:

    - Combines multiple batches
    - Averages posterior probabilities across different chains from step 2's PGM outputs.
    - Applies gene-wise filters to the averaged data.
    - Performs min-max normalization and form the final phenotype-relevant GRN.
    
Usage:

    - Run the script from the command line "python InPheC_Step3.py" directly.
    - For detailed argument options, use the help option: `python InPheC_Step3.py --help`

Note:
Ensure the Python environment is properly set up and all dependencies are installed before running this script.
"""
import argparse
import logging
import os
import platform
import random
import re
import sys
import time
import warnings
from multiprocessing import Pool
from typing import Any, Dict, Tuple, Union

import cloudpickle
import numpy as np
import pandas as pd
from tqdm import tqdm

from utils_logo import print_logo


class Step3ArgumentParser:
    """
    Parses command-line arguments for the script.

    This class configures and retrieves command-line options, returning the parsed arguments as an argparse.Namespace object. 
    It sets default values for input and output directories, the number of batches, and the number of burn-in iterations for posterior processing.

    :return: An argparse.Namespace object containing all the command-line arguments.
    :rtype: argparse.Namespace

    :Example:

    args = ArgumentParser.parse()
    print(args.in_dir)  # Prints the input directory specified in the command-line
    """
    @staticmethod
    def parse() -> argparse.Namespace:
        """
        Sets up the command-line arguments and parses them.

        :return: Parsed arguments from the command line.
        :rtype: argparse.Namespace
        """
        parser = argparse.ArgumentParser()
        parser.add_argument("-id", "--in_dir", default="./step2")
        parser.add_argument("-od", "--out_dir", default="./step3")
        parser.add_argument(
            "-nba", "--n_batches", type=int, default=13, help="number of the batches in the input folder."
        )
        parser.add_argument(
            "-nb",
            "--n_burn",
            type=int,
            default=400,
            help="number of iterations to burn, \
            it's not uncommon to discard the first 10 percent of the post-tuning draws as burn-in",
        )
        return parser.parse_args()


class FileNotFound(Exception):
    """
    Custom exception for handling file-not-found errors.

    This exception is raised when a required file is not found in the specified directory, 
    extending the base Exception class.

    :param str message: The error message to display.
    """
    pass


class ParameterExtractor:
    """
    Extracts parameters from file names in the specified directory.

    Designed to read file names in a given directory and extract key parameters such as T_prior, n_c, n_i, and n_t, which are essential for processing steps in the pipeline.

    :param in_dir: Input directory where files are located.
    :type in_dir: str

    :return: A dictionary containing extracted parameters with keys like 'sample_fn', 'T_prior', 'n_c', 'n_i', and 'n_t'.
    :rtype: Dict[str, Union[str, float]]

    :raises FileNotFound: If no relevant pickle file is found in the directory.
    """
    def __init__(self, in_dir: str):
        """
        Initializes the ParameterExtractor with the specified input directory.

        :param in_dir: Directory containing the files to extract parameters from.
        :type in_dir: str
        """
        self.in_dir = in_dir

    def extract(self) -> Dict[str, Union[str, float]]:
        """
        Extracts parameters from the filenames in the input directory.

        :return: A dictionary containing the extracted parameters.
        :rtype: Dict[str, Union[str, float]]
        :raises FileNotFound: If no relevant pickle file is found in the directory.
        """
        sample_file = next((f for f in os.listdir(self.in_dir) if f.startswith("RnCnP_") and f.endswith(".pkl")), None)
        if sample_file is None:
            raise FileNotFound("No pickle file found in directory.")

        return {
            "sample_fn": sample_file,
            **dict(zip(["T_prior", "n_c", "n_i", "n_t"], map(float, re.findall(r"[-+]?\d*\.\d+|\d+", sample_file)))),
        }


class PosteriorCombiner:
    """
    Combines posterior distributions from multiple batches and chains.

    Processes each batch, extracting and combining posterior distributions into a single DataFrame from pickle files.

    :param in_dir: Input directory containing the pickle files.
    :type in_dir: str
    :param n_batches: Number of batches to process.
    :type n_batches: int
    :param n_burn: Number of burn-in iterations.
    :type n_burn: int
    :param common_params: Common parameters used across the batches.
    :type common_params: dict

    :return: Tuple containing paths to the combined posterior file and the in2pgm data.
    :rtype: tuple
    """

    def __init__(self, in_dir: str, n_batches: int, n_burn: int, common_params: dict):
        """
        Initializes the PosteriorCombiner with the necessary parameters.

        :param in_dir: Directory containing the input pickle files.
        :param n_batches: Number of batches to process.
        :param n_burn: Number of burn-in iterations to discard.
        :param common_params: Common parameters used across the batches.
        """
        self.in_dir = in_dir
        self.n_batches = n_batches
        self.n_burn = n_burn
        self.common_params = common_params

    def _ordinal(self, n: int) -> str:
        """
        Converts an integer to its ordinal representation (e.g., 1st, 2nd, 3rd).

        :param n: Number to convert.
        :type n: int

        :return: Ordinal string of the number.
        :rtype: str
        """
        suffix = ["th", "st", "nd", "rd", "th"][min(n % 10, 4)]
        if 11 <= (n % 100) <= 13:
            suffix = "th"
        return str(n) + suffix

    def _process_chain(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Processes each batch within the chain, extracting and combining posterior distributions.

        :return: Tuple of DataFrames containing the recovered posteriors and the genes in the PGM.
        :rtype: Tuple[pd.DataFrame, pd.DataFrame]
        :raises FileNotFound: If a pickle file is not found in the directory for a specific batch.
        """
        logging.info("Processing all chains together")

        T_recovered_allGenes_ls = []
        genes_in2pgm = []

        # Inner progress bar for tracking progress within the chain
        for i_b in tqdm(range(self.n_batches), desc="Processing Batches", leave=False):
            logging.debug(f"\t batch id = {i_b}")  # Debug message

            # Locate the pickle
            file_suffix = f"batch{i_b}.pkl"
            f = next(
                (
                    filename
                    for filename in os.listdir(self.in_dir)
                    if os.path.isfile(os.path.join(self.in_dir, filename)) and filename.endswith(file_suffix)
                ),
                None,
            )
            logging.info(f)
            if f is None:
                raise FileNotFound(f" Oops! No pickle file found in directory for batch {i_b}, please check.")
            else:
                with open(os.path.join(self.in_dir, f), "rb") as file:
                    genes_onebatch = cloudpickle.load(file)

                for gene, values in genes_onebatch["trace"].items():
                    T_samples = values.posterior[
                        "T"
                    ]  # (4, 1000, 1, 23, 1): (n_chian, n_samples, n_genes, n_tfs, n_lineage)
                    T_burn = T_samples.sel(draw=slice(self.n_burn, None))  # (4, 900, 1, 23, 1)

                    # Pooling across all chains
                    T_pooled = np.mean(
                        T_burn, axis=(0, 1)
                    )  # Averaging over chains and draws (not the same as concatenation of chains!)

                    jj_tfs = genes_onebatch["info_per_batch"][gene]["TF"]
                    T_recovered_jj = pd.DataFrame(T_pooled.values.reshape(1, -1), index=[gene], columns=list(jj_tfs))
                    T_recovered_allGenes_ls.append(T_recovered_jj)

            # load genes_onebatch into genes_in2pgm
            genes_in2pgm.append(pd.concat(genes_onebatch["info_per_batch"].values()))

        return pd.concat(T_recovered_allGenes_ls, ignore_index=False), pd.concat(genes_in2pgm)

    def combine(self, out_dir: str):
        """
        Combines the posteriors and saves the results to specified output files.

        :param out_dir: Directory where the output files will be saved.
        :type out_dir: str

        :return: Tuple of paths to the combined in2pgm and posterior files.
        :rtype: Tuple[str, str]
        """
        # Part 1: combine posteriors
        result, inputs_to_pgm = self._process_chain()

        logging.info(f"Completed combining chains: Saving to ./step3/step2-in2pgm.pkl")

        # output: save the genes_in2pgm data outside the loop, since the in2pgm is the same for all processes/batches.
        fout_in2pgm = os.path.join(out_dir, f"step2-in2pgm.pkl")
        with open(fout_in2pgm, "wb") as handle:
            cloudpickle.dump(inputs_to_pgm, handle)

        # output: pickle the combined posteriors
        fout_posterior = os.path.join(
            out_dir, re.sub(r"(_noburn).*", f"_nb{self.n_burn}.pkl", self.common_params["sample_fn"])
        )
        with open(fout_posterior, "wb") as handle:
            cloudpickle.dump(result, handle)

        return fout_posterior


class DataProcessor:
    """
    Processes and analyzes the combined posterior data from probabilistic graphical models.

    This class loads and processes the in2pgm and posterior data, combines them, and applies filtering and normalization 
    to produce a final gene regulatory network (GRN).

    :param fn1: Path to the first input file containing in2pgm data.
    :type fn1: str
    :param fn2: Path to the second input file containing posterior data.
    :type fn2: str
    """

    def __init__(self, fn: str):
        self.fn = fn
        self.pgmout = self._load_data()

    def _load_data(self):
        """
        Loads data from the specified files.

        :return: Tuple of DataFrames containing the in2pgm data and posterior data.
        :rtype: Tuple[pd.DataFrame, pd.DataFrame]
        :raises ValueError: If either of the files is empty.
        """
        # Check if files are empty
        if os.path.getsize(self.fn) == 0:
            raise ValueError(f"Empty file detected. Please delete files in ./step3/tmp and re-run step3.")
        # load
        with open(self.fn, "rb") as file:
            data = cloudpickle.load(file)
        return data

    def combine_in2pgm_outposterior(self) -> pd.DataFrame:
        """
        Combines the in2pgm and posterior data into a single DataFrame.

        This function merges the two datasets on the TF and Gene columns, effectively creating an edge list for the GRN.

        :return: Combined DataFrame with the in2pgm and posterior data.
        :rtype: pd.DataFrame
        """
        # convert the pivot table to an edge list
        melted_df = self.pgmout.reset_index().melt(
            id_vars="index", value_vars=self.pgmout.columns, var_name="TF", value_name="raw_posterior_mean"
        )
        melted_df = melted_df.rename(columns={"index": "Gene"})

        # combine data1 and data2
        combined_df = pd.merge(self.in2pgm, melted_df, on=["TF", "Gene"], how="inner")
        return combined_df

    def _calculate_max_differential_metrics(self, case: pd.Series, a_col: str) -> dict:
        """
        Calculates differential metrics to identify significant thresholds for filtering.

        This function computes differences and ratios between consecutive elements in a sorted series 
        to determine the most significant threshold.

        :param case: Series of data points to analyze.
        :type case: pd.Series
        :param a_col: Column identifier in the DataFrame to calculate metrics.
        :type a_col: str

        :return: Dictionary of calculated metrics.
        :rtype: Dict[str, float]
        :raises ValueError: If the computed metrics do not agree on a single cutoff point.
        """
        y = case[a_col].values
        info = {
            "Diff: yi-yj": [i - j for i, j in zip(y[:-1], y[1:])],
            "ratio_1: (yi-yj) / y_max": [(i - j) / max(y) for i, j in zip(y[:-1], y[1:])],
            "ratio_4: (yi-yj) / (y_max-y_ref)": [(i - j) / (max(y) - min(y)) for i, j in zip(y[:-1], y[1:])],
        }
        maxx_dict = pd.DataFrame.from_dict(info).idxmax().to_dict()
        maxLoc = set(maxx_dict.values())
        if len(maxLoc) > 1:
            raise ValueError("3 metric do not agree with each other")
        return {"metric-Diff/R1/R4": y[maxx_dict["ratio_4: (yi-yj) / (y_max-y_ref)"] + 1]}

    def _apply_normalization_to_each_gene(self, df, coi) -> pd.DataFrame:
        """
        Applies min-max normalization to each gene in the DataFrame.

        :param df: DataFrame containing raw scores to normalize.
        :type df: pd.DataFrame
        :param coi: Column of interest in the DataFrame to normalize.
        :type coi: str

        :return: DataFrame with normalization applied.
        :rtype: pd.DataFrame
        """
        df["min_max_normalized"] = df.groupby("Gene")[coi].transform(lambda x: (x - x.min()) / (x.max() - x.min()))
        return df

    def _apply_filters_to_each_gene(self, df, coi, Tp_val) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Applies filtering based on differential metrics and threshold values to each gene.

        This function determines filtering criteria and then applies these filters to identify significant gene-TF interactions.

        :param df: DataFrame to apply filters.
        :type df: pd.DataFrame
        :param coi: Column of interest where filtering is applied.
        :type coi: str
        :param Tp_val: Threshold prior value used in filtering.
        :type Tp_val: float

        :return: Tuple containing the filtered DataFrame and a DataFrame with filter criteria.
        :rtype: Tuple[pd.DataFrame, pd.DataFrame]
        """
        # ---- find filtering criterion ----
        # filter 1: calculates differences and ratios between consecutive elements in a sorted series
        filter_criterion = self._calculate_max_differential_metrics(df, coi)

        # filter 2: get value of T_prior
        filter_criterion["T_prior"] = Tp_val

        # filter 3: get Tij values of 3 controls in each gene
        ref_tfs = df[df["TF"].str.startswith("ref_")]
        ref_tf_values = dict(zip(ref_tfs["TF"], ref_tfs[coi]))
        filter_criterion.update(ref_tf_values)

        # find the maximum value of the filter criterion
        max_value = max(filter_criterion.values())

        # ---- apply filters ----
        mask_df = df.copy()
        mask_df["final_importance_score"] = mask_df[mask_df[coi] > max_value][coi]
        # mask the value < max to 0
        # mask_df["final_importance_score"] = mask_df[coi].where(df[coi] > max_value, 0)

        filter_criterion["max_val"] = max_value

        return mask_df, pd.DataFrame([filter_criterion])  # List wrapping ensures a single-row DataFrame

    def form_final_grn(self, col_raw_score, Tp_val) -> pd.DataFrame:
        """
        Forms the final gene regulatory network (GRN) based on normalized and filtered data.

        :param col_raw_score: Column name of the raw scores used for final processing.
        :type col_raw_score: str
        :param Tp_val: Threshold prior value used in forming the GRN.
        :type Tp_val: float

        :return: Tuple containing the final GRN DataFrame and gene-specific thresholds.
        :rtype: pd.DataFrame
        """
        processed_dfs = []
        
        # Group by 'Gene' and apply normalization and filtering to each group
        long_format_df = self.pgmout.reset_index().melt(id_vars='index', var_name='TF', value_name=col_raw_score).rename(columns={'index': 'Gene'})
        for gene, group_df in tqdm(long_format_df.groupby("Gene")):
            # sort the group by 'raw_posterior_mean' in descending order
            group_df = group_df.sort_values(by=col_raw_score, ascending=False)

            norm_df = self._apply_normalization_to_each_gene(group_df, col_raw_score)
            flt_norm_df, _ = self._apply_filters_to_each_gene(norm_df, "min_max_normalized", Tp_val)

            # Append the processed group to the list
            processed_dfs.append(flt_norm_df)

        # Combine all processed groups into a single df
        tmp_df = pd.concat(processed_dfs, ignore_index=True)
        
        # Sparse network
        formatted_df = tmp_df[~tmp_df['TF'].str.startswith('ref_')].fillna('').pivot(index='Gene', columns='TF', values='final_importance_score').apply(pd.to_numeric, errors='coerce')
        return formatted_df


class Utils:
    @staticmethod
    def configure_logging(log_dir: str) -> None:
        """
        Configures the logging system for the script, directing log output to a file in a specified directory.

        :param log_dir: Directory where the log file will be stored.
        :type log_dir: str
        """
        os.makedirs(log_dir, exist_ok=True)
        logging.basicConfig(
            level=logging.INFO,  # Set the log level
            format="%(asctime)s - %(levelname)s - %(message)s",  # Set the log format
            filename=os.path.join(log_dir, "step3.log"),  # Optional: Specify a file to write logs to
            filemode="a",  # 'w' for overwrite, 'a' for append
        )

    @staticmethod
    def check_pkl_file_count(in_dir: str, expected_count: int):
        """
        Checks if the number of .pkl files in the input directory matches the expected number.

        :param in_dir: Directory to check for .pkl files.
        :type in_dir: str
        :param expected_count: The number of .pkl files expected in the directory.
        :type expected_count: int
        :raises ValueError: If the number of .pkl files does not match the expected count.
        """
        pkl_files = [f for f in os.listdir(in_dir) if f.endswith(".pkl")]
        if len(pkl_files) != expected_count:
            raise ValueError(
                f"Number of .pkl files found ({len(pkl_files)}) does not match the expected number ({expected_count}). "
                "Please ensure that step 2 was completed successfully, all batches were processed, "
                "and the `-nba` argument is correctly set according to the batch size and the total number of genes analyzed in step 2."
            )

    @staticmethod
    def print_python_details():
        """
        Prints detailed information about the Python environment running the script.
        Outputs various Python-related system information, including version, compiler, and other system details.
        """
        # Print the rest of the details in original language
        print("\nThank you for using InPheRNo-ChIP! Python Environment Details:")
        print("----------------------------")
        print("Python version:", platform.python_version())
        print("Compiler:", platform.python_compiler())
        print("Build:", platform.python_build())
        print("System:", platform.system())
        print("Release:", platform.release())
        print("Machine:", platform.machine())
        print("Processor:", platform.processor())
        print("sys.version:", sys.version.replace("\n", " "))
        print("-------------Bye-------------")


def step3_main() -> None:
    """
    Executes the main functionality of the Step 3 script in the InPheRNo-ChIP project. This function manages the workflow
    that combines and processes posterior distributions from probabilistic graphical models (PGMs) generated in Step 2.

    The main steps include:
    - Configuring logging for debugging and information tracking.
    - Parsing command line arguments to get directories and processing parameters.
    - Checking the count of .pkl files to ensure all expected files are present.
    - Creating necessary directories for outputs.
    - Extracting and logging parameters from file names.
    - Checking for and combining existing processed data or processing new data.
    - Merging, filtering, and normalizing data to form the final Gene Regulatory Network (GRN).
    - Saving the final outputs in an Excel file with multiple sheets for easy review.

    :raises FileNotFoundError: If required files are not found in the specified directories.
    :raises ValueError: If processing parameters are incorrect or if files are empty.
    """
    print_logo("Running step3")
    
    # Use the Utils class for configuring logging
    Utils.configure_logging("./logs")

    args = Step3ArgumentParser.parse()
    # Check if the number of .pkl files matches the number of batches
    Utils.check_pkl_file_count(args.in_dir, args.n_batches)
    
    # output directory
    os.makedirs(args.out_dir, exist_ok=True)
    intermediate_folder = os.path.join(args.out_dir, "tmp")
    os.makedirs(intermediate_folder, exist_ok=True)
    
    # Initialize the process timer
    start_time = time.time()

    common_params = ParameterExtractor(in_dir=args.in_dir).extract()
    logging.info(
        f"Running [step3A.py], input folder: {args.in_dir}; \nn_chain = {common_params['n_c']}, n_batch = {args.n_batches}"
    )

    # ----- part a. combine posteriors -----
    combiner = PosteriorCombiner(
        in_dir=args.in_dir, n_batches=args.n_batches, n_burn=args.n_burn, common_params=common_params
    )
    faddr_posterior = combiner.combine(out_dir=intermediate_folder)
    logging.info(
        f"Step 3 part a - posteriors combination is completed, TIME:{time.time() - start_time}. \n Now applying filtering and normalization..."
    )

    # ----- part b. filter and normalize -----
    postprocess = DataProcessor(faddr_posterior)
    sort_by = "raw_posterior_mean"
    out_grn = postprocess.form_final_grn(sort_by, common_params["T_prior"])

    # save to excel with 2 tabs
    base_name = os.path.basename(faddr_posterior)  # Extracts 'posterior_file.pkl'
    fout = os.path.join(args.out_dir, f"{base_name.replace('.pkl', '_minMax.csv')}")
    out_grn.to_csv(fout)

    # Use the Utils class to print Python environment details at the end
    Utils.print_python_details()


if __name__ == "__main__":
    step3_main()
