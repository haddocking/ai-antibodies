"""Main CLI."""
import argparse
import logging
import sys
from pathlib import Path
import os
import time
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
from aiabs.modules.analysis import (
    read_caprievals,
    generate_tables,
#    create_comparison_graphs
)

logging.basicConfig(filename="aiabsanalysis")
log = logging.getLogger("aiabsanalysis")

ch = logging.StreamHandler()
formatter = logging.Formatter(
    " [%(asctime)s %(module)s:L%(lineno)d %(levelname)s] %(message)s"
)
ch.setFormatter(formatter)
log.addHandler(ch)


argument_parser = argparse.ArgumentParser()

argument_parser.add_argument(
    "--ref_folder",
    help=""
    )

argument_parser.add_argument(
    "--capri_string",
    help="",
    default="caprieval"
    )

argument_parser.add_argument(
    "--exclude_pdbs",
    help="list of comma-separated pdb to exclude",
    default="7kpj,7kn4"
    )

def load_args(arguments):
    """
    Load argument parser.

    Parameters
    ----------
    arguments : argparse.ArgumentParser
        Argument parser.

    Returns
    -------
    cmd : argparse.Namespace
        Parsed command-line arguments.

    """
    return arguments.parse_args()


def cli(arguments, main_func):
    """
    Command-line interface entry point.

    Parameters
    ----------
    arguments : argparse.ArgumentParser
        Argument parser.
    main_func : function
        Main function.

    """
    cmd = load_args(arguments)
    main_func(**vars(cmd))


def maincli():
    """Execute main client."""
    cli(argument_parser, main)


def main(
    ref_folder,
    capri_string,
    exclude_pdbs,
):
    """Main function."""
    log.setLevel("INFO")
    start_time = time.time()
    log.info(f"CWD {os.getcwd()}")
   
    log.info(f"ref_folder {ref_folder}")

    # creating analysis folders
    apath = Path(ref_folder, "analysis")
    if not apath.exists():
        os.mkdir(apath)

    to_exclude = []
    if exclude_pdbs:
        to_exclude = exclude_pdbs.split(",")
    log.info(f"pdbs to exclude from the analysis {to_exclude}")

    # 1 initialising antibodies dictionary
    #antibodies = [el for el in os.listdir(ref_folder) if len(el) == 10 and el not in to_exclude]
    antibodies = [el for el in os.listdir(ref_folder) if (len(el) == 10 or len(el) == 4) and el not in to_exclude]
    log.info(f"{len(antibodies)} antibodies in folder {ref_folder}: {antibodies}")
    
    capri_ss = {}
    capri_clt = {}
    for pdb in antibodies:
        log.info(f"processing pdb {pdb}")
        run_dirs = [fl for fl in os.listdir(Path(ref_folder, pdb)) if fl.startswith("run")]
        #log.info(f"run_dirs are {run_dirs}")
        for run in run_dirs:
            if capri_string == "flexref_analysis":
                rel_path = Path(ref_folder, pdb, run, "analysis")
            else:
                rel_path = Path(ref_folder, pdb, run)
            capri_folders = [fold.split("_")[0] for fold in os.listdir(rel_path) if fold.endswith(capri_string)]
            #log.info(f"capri_folders are {capri_folders}")
            for capri in capri_folders:
                if capri not in capri_ss.keys():
                    capri_ss[capri] = {}
                df, df_clt = read_caprievals(capri, rel_path, capri_string)
                if df is None:
                    log.warning(f"None single structure for {pdb} {capri}")
                else:
                    if pdb not in capri_ss[capri]:
                        capri_ss[capri][pdb] = {}
                    capri_ss[capri][pdb][run] = df
                if df_clt is None or np.unique(df_clt["cluster_rank"]) == ["-"]:
                    log.debug(f"None clustering for {pdb} {capri}")
                else:
                    if capri not in capri_clt.keys():
                        capri_clt[capri] = {}
                    if pdb not in capri_clt[capri]:
                        capri_clt[capri][pdb] = {}
                    capri_clt[capri][pdb][run] = df_clt
    log.info(f"capri reading took {round((time.time() - start_time), 3)} seconds")
    # aggregated tables
    generate_tables(capri_ss, capri_clt, ref_folder)

    log.info(f"aiabs-analysis succesful for input folder {ref_folder}")
    elap_time = round((time.time() - start_time), 3)
    log.info(f"aiabs-analysis run took {elap_time} seconds")

if __name__ == "__main__":
    sys.exit(maincli())
