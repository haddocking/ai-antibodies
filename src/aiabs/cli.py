"""Main CLI."""
import argparse
import logging
import sys
from pathlib import Path
import os
import time
import pandas as pd

from aiabs.modules.pdb import generate_fasta, generate_target, preprocess_ab_file, preprocess_ant_file, get_resids, preprocess_ensemble_file
from aiabs.modules.file import setup_directory, prepare_jobs
from aiabs.modules.restraints import (
    check_align,
    get_restrs_from_pairs,
    write_restraints,
    write_zdock_restraints,
    get_restrs_from_list,
    write_split_tbl,
    write_unambig_tbl,
    join_tbls,
    standardize_unambig_tbls,
    get_solvent_exposed_residues,
    )
from aiabs.modules.cluster import cluster_antibodies

logging.basicConfig(filename="aiabslog")
log = logging.getLogger("aiabslog")

ch = logging.StreamHandler()
formatter = logging.Formatter(
    " [%(asctime)s %(module)s:L%(lineno)d %(levelname)s] %(message)s"
)
ch.setFormatter(formatter)
log.addHandler(ch)


argument_parser = argparse.ArgumentParser()

argument_parser.add_argument(
    "pdb",
    help=""
    )

argument_parser.add_argument(
    "--input_dir",
    help=""
)

argument_parser.add_argument(
    "--output_dir",
    help="absolute path to the output directory"
    )

argument_parser.add_argument(
    "--flag_str",
    help="flag to specify if one non-standard structure must be specified in the ensemble. (example: ig:mdref,ab:omm",
    default=None
    )

argument_parser.add_argument(
    "--ex_files_path",
    help="if specified, the path to the example_files folder",
    default="example_files"
    )

argument_parser.add_argument(
    "--act_act_path",
    help="if specified, the path to the generate_act_act.sh file",
    default="generate_act_act.sh"
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
    pdb,
    input_dir,
    output_dir,
    flag_str,
    ex_files_path,
    act_act_path
):
    """Main function."""
    log.setLevel("DEBUG")
    start_time = time.time()
    input_files = [Path(input_dir, fname) for fname in os.listdir(input_dir)]

    # checking restraint sh file
    act_act_path = Path(act_act_path)
    if not act_act_path.exists():
        raise Exception(f"act_act file {act_act_path} does not exist")
    
    if len(input_files) == 0:
        raise Exception("no input files available")
    pdb_files = [fpath for fpath in input_files if fpath.name.endswith("pdb")]

    flag_dict = None
    if flag_str:
        flag_dict = {}
        for el in flag_str.split(","):
            splt = el.split(":")
            keyword = splt[0]
            refinement_method = splt[1]
            flag_dict[keyword] = refinement_method
    
    # 1 setting up directory
    pdb_dir, pdb_datadir, pdb_alndir = setup_directory(output_dir, pdb)
    
    # 2 preprocess antibodies
    ab_keywords = ["ABlooper", "ABodyBuilder2", "AF2", "IgFold"]
    inp_ab_files, out_ab_files, ab_fastas = {}, {}, {}

    for keyword in ab_keywords:
        log.info(f"preprocessing {keyword} file")
        
        inp_ab_files[keyword] = [el for el in pdb_files if el.name.startswith(keyword)][0]
        out_ab_files[keyword] = preprocess_ab_file(inp_ab_files[keyword],
                                                pdb_datadir,
                                                f"{pdb}_{keyword}_haddock-ready.pdb")
        ab_fastas[keyword] = generate_fasta(inp_ab_files[keyword], pdb_datadir)
    
    # ABodyBuilder2 input : create also an ensemble
    abb2_files = os.listdir(Path(input_dir, f"ABB2_ensemble_models_{pdb}"))
    inp_abb2_files = [Path(input_dir, f"ABB2_ensemble_models_{pdb}", el) for el in abb2_files]
    log.info(f"abb2_files available {inp_abb2_files}")
    abb2_fastas, proc_abb2_files = {}, {}
    for abb2_fl in sorted(inp_abb2_files):
        rank = abb2_fl.stem.split("_")[0]
        proc_abb2_files[rank] = preprocess_ab_file(abb2_fl,
                                                pdb_datadir,
                                                f"{pdb}_abb2{rank}_haddock-ready.pdb")
        abb2_fastas[rank] = generate_fasta(abb2_fl, pdb_datadir)
    log.info(f"proc_abb2_files {proc_abb2_files}")
    # abb2 ensemble creation
    ensemble_fname = preprocess_ensemble_file(proc_abb2_files,
                                              pdb_datadir,
                                              f"{pdb}_ABodyBuilder2ens_haddock-ready.pdb",
                                              flag_dict)
    
    # Alphafold input: create also an ensemble
    af_files = os.listdir(Path(input_dir, f"AF2_ensemble_models_{pdb}"))
    inp_af_files = [Path(input_dir, f"AF2_ensemble_models_{pdb}", el) for el in af_files]
    log.info(f"af_files available {inp_af_files}")
    af_fastas, proc_af_files = {}, {}
    for af_fl in sorted(inp_af_files):
        rank = f'afrank{af_fl.stem.split("rank_")[-1]}'
        proc_af_files[rank] = preprocess_ab_file(af_fl,
                                                pdb_datadir,
                                                f"{pdb}_{rank}_haddock-ready.pdb")
        af_fastas[rank] = generate_fasta(af_fl, pdb_datadir)
    log.info(f"proc_af_files {proc_af_files}")
    # af ensemble creation
    ensemble_fname = preprocess_ensemble_file(proc_af_files,
                                                pdb_datadir,
                                                f"{pdb}_AF2ens_haddock-ready.pdb",
                                                flag_dict)
    
    # IgFold input: create also an ensemble
    ig_files = os.listdir(Path(input_dir, f"IgFold_ensemble_models_{pdb}"))
    inp_ig_files = [Path(input_dir, f"IgFold_ensemble_models_{pdb}", el) for el in ig_files]
    log.info(f"ig_files available {inp_ig_files}")
    ig_fastas, proc_ig_files = {}, {}
    for ig_fl in sorted(inp_ig_files):
        rank = f'igrank{ig_fl.stem.split("rank_")[-1]}'
        proc_ig_files[rank] = preprocess_ab_file(ig_fl,
                                                pdb_datadir,
                                                f"{pdb}_{rank}_haddock-ready.pdb")
        ig_fastas[rank] = generate_fasta(ig_fl, pdb_datadir)
    log.info(f"proc_ig_files {proc_ig_files}")
    # ig ensemble creation
    ensemble_fname = preprocess_ensemble_file(proc_ig_files,
                                              pdb_datadir,
                                              f"{pdb}_IgFoldens_haddock-ready.pdb",
                                              flag_dict)

    # 3 preprocess antigen
    inp_ant_file = [el for el in pdb_files if el.name.endswith("antigen.pdb")][0]
    log.info(f"preprocessing antigen file {inp_ant_file}")
    out_ant_file, antigen_residues = preprocess_ant_file(inp_ant_file, pdb_datadir, f"{pdb}_antigen_haddock-ready.pdb")
    log.info(f"antigen file is {out_ant_file}")
    # af2 antigen
    inp_ant_af2_file = [el for el in pdb_files if el.name.endswith("_antigen_model.pdb")][0]
    out_ant_af2_file, antigen_af2_residues = preprocess_ant_file(inp_ant_af2_file, pdb_datadir, f"{pdb}_af2_antigen_haddock-ready.pdb")
    log.info(f"af2 antigen file is {out_ant_af2_file}")
    # fastas for antigen
    ant_fastas = {}
    ant_fastas["antigen"] = generate_fasta(inp_ant_file, pdb_datadir)
    ant_fastas["af2_antigen"] = generate_fasta(inp_ant_af2_file, pdb_datadir)


    # 4 preprocess target
    log.info(f"expected target file name = {pdb}_true_complex.pdb")
    log.info(f"looking for it in pdb_files {[el.name for el in pdb_files]}")
    target_file = [el for el in pdb_files if el.name == f"{pdb}_true_complex.pdb"]
    if len(target_file) == 0:
        target_file = [el for el in pdb_files if el.name.endswith(f"true_complex.pdb")]
    log.info(f"target_file {target_file}")
    # extracting chains
    c_list_ab = Path(input_dir, f"{pdb}_residue_constraints_antibody.csv")
    c_list_df = pd.read_csv(c_list_ab)
    
    chains = c_list_df["chain"].unique()
    chains = [chains[0], chains[1]] # heavy chain should be first I imagine
    log.info(f"antibody heavy-light chain {chains}")
    target_ab_fname, target_ant_fname, target_fname, antibody_bound_resids = generate_target(target_file[0], pdb_datadir, pdb, chains)
    print(f"antibody_bound_resids {antibody_bound_resids}")
    log.info(f"target_ab_fname {target_ab_fname}")
    ab_fastas["target"] = generate_fasta(target_ab_fname, pdb_datadir)
    # 5 check alignments and residue ids
    dict.update(ab_fastas, abb2_fastas)
    # add af and igfold sequences
    dict.update(ab_fastas, af_fastas)
    dict.update(ab_fastas, ig_fastas)

    check_aligned = check_align(ab_fastas, pdb_alndir)

    # 5.b check antigen alignments
    check_aligned_antigen = check_align(ant_fastas, pdb_alndir, pdb)

    antibodies_resids = {}
    for inp_ab in inp_ab_files.keys():
        antibodies_resids[inp_ab] = get_resids(inp_ab_files[inp_ab])
    # getting resids from abb2 structures
    for n in range(len(inp_abb2_files)):
        key = f"abb2_{n}"
        antibodies_resids[key] = get_resids(inp_abb2_files[n])
    log.info(f"antibodies_resids {antibodies_resids.keys()}")
    for ab1 in antibodies_resids.keys():
        for ab2 in antibodies_resids.keys():
            if (ab1 != ab2):
                residues_ids_i = [el[1:] for el in antibodies_resids[ab1]]
                residues_ids_j = [el[1:] for el in antibodies_resids[ab2]]
                if antibodies_resids[ab1] != antibodies_resids[ab2]:
                    if residues_ids_i == residues_ids_j:
                        log.warning(f"resids mismatch between {ab1} and {ab2}: residue ids are the same, but the chains are different")
                    else:
                        raise Exception(f"resids mismatch between {ab1} and {ab2}: residue ids are different")
    # antigen residues
    ant_filename = [el for el in pdb_files if el.name.endswith("antigen.pdb")][0]
    #antigen_residues = get_resids(ant_filename)
    print(f"antigen_residues {antigen_residues}")
    # AF2 antigen residues
    ant_af2_filename = [el for el in pdb_files if el.name.endswith("AF2_antigen_model.pdb")][0]
    #antigen_af2_residues = get_resids(ant_af2_filename)

    # 6 coupled restraints
    c_pairs_file = Path(input_dir, f"{pdb}_constraint_pairs.txt")
    constr_ab, constr_ant = get_restrs_from_pairs(c_pairs_file,
                                                  antibodies_resids["IgFold"],
                                                  antigen_residues,
                                                  chain_mapping={chains[0]: "H", chains[1]: "L"})
    restr_line_filename = Path(pdb_datadir, f"{pdb}_restraint_line.txt")
    tbl_filename = Path(pdb_datadir, f"{pdb}_ambig_Para_Epi.tbl")

    write_restraints({"A": constr_ab, "B": constr_ant}, restr_line_filename, tbl_filename, act_act_path)

    # bound antibody restraints
    constr_ab_bound, constr_ant_bound = get_restrs_from_pairs(c_pairs_file,
                                                              antibody_bound_resids,
                                                              antigen_residues,
                                                              chain_mapping={chains[0]: chains[0], chains[1]: chains[1]})
    print(f"constr_ant_bound {constr_ant_bound}")
    restr_line_bound_filename = Path(pdb_datadir, f"{pdb}_bound_restraint_line.txt")
    tbl_filename_bound = Path(pdb_datadir, f"{pdb}_ambig_bound_Para_Epi.tbl")
    write_restraints({"A": constr_ab_bound, "B": constr_ant_bound}, restr_line_bound_filename, tbl_filename_bound, act_act_path)
                                                              
    # 6.b ZDOCK coupled restraints
    zdock_ab_restr_filename = Path(pdb_datadir, f"{pdb}_zdock_ab_Para_Epi.csv")
    zdock_ant_restr_filename = Path(pdb_datadir, f"{pdb}_zdock_ant_Para_Epi.csv")
    write_zdock_restraints({"A": constr_ab}, zdock_ab_restr_filename)
    write_zdock_restraints({"B": constr_ant}, zdock_ant_restr_filename)

    # 6.c af2 coupled restraints
    af2_c_pairs_file = Path(input_dir, f"{pdb}_af2_constraint_pairs.txt")
    constr_ab_af2, constr_ant_af2 = get_restrs_from_pairs(af2_c_pairs_file,
                                                  antibodies_resids["IgFold"],
                                                  antigen_af2_residues,
                                                  chain_mapping={chains[0]: "H", chains[1]: "L"})
    # 22/5/2024: unfortunately this is not true anymore
    #assert constr_ab == constr_ab_af2
    restr_line_af2_filename = Path(pdb_datadir, f"{pdb}_af2_restraint_line.txt")
    tbl_filename_af2 = Path(pdb_datadir, f"{pdb}_ambig_af2_Para_Epi.tbl")
    write_restraints({"A": constr_ab_af2, "B": constr_ant_af2}, restr_line_af2_filename, tbl_filename_af2, act_act_path)

    # 6.d ZDOCK af2 coupled restraints (no need to regenerate those of the antibody)
    zdock_af2ant_restr_filename = Path(pdb_datadir, f"{pdb}_zdock_af2ant_Para_Epi.csv")
    write_zdock_restraints({"B": constr_ant_af2}, zdock_af2ant_restr_filename)

    # 7 decoupled restraints
    c_list_ab = Path(input_dir, f"{pdb}_residue_constraints_antibody.csv")
    c_list_ant = Path(input_dir, f"{pdb}_residue_constraints_antigen.csv")
    c_list_ant_af2 = Path(input_dir, f"{pdb}_af2_residue_constraints_antigen.csv")
    constr_ab_dec = get_restrs_from_list(c_list_ab, antibodies_resids["IgFold"], chain_mapping={chains[0]: "H", chains[1]: "L"})
    constr_ant_dec = get_restrs_from_list(c_list_ant, antigen_residues)
    constr_ant_af2_dec = get_restrs_from_list(c_list_ant_af2, antigen_af2_residues)
    
    restr_line_dec_filename = Path(pdb_datadir, f"{pdb}_decoupled_restraint_line.txt")
    tbl_filename_dec = Path(pdb_datadir, f"{pdb}_ambig_CDR_EpiVag_act-act.tbl")
    write_restraints({"A": constr_ab_dec, "B": constr_ant_dec}, restr_line_dec_filename, tbl_filename_dec, act_act_path)

    new_tbl_filename = Path(pdb_datadir, f"{pdb}_ambig_CDR_EpiVag.tbl")
    splt_string = f"!{os.linesep}! HADDOCK AIR restraints for 2nd partner"
    write_split_tbl(tbl_filename_dec, splt_string, new_tbl_filename)

    # bound antibody restraints
    constr_ab_bound_dec = get_restrs_from_list(c_list_ab, antibody_bound_resids, chain_mapping={chains[0]: chains[0], chains[1]: chains[1]})

    restr_line_bound_dec_filename = Path(pdb_datadir, f"{pdb}_bound_decoupled_restraint_line.txt")
    tbl_filename_bound_dec = Path(pdb_datadir, f"{pdb}_ambig_bound_CDR_EpiVag_act-act.tbl")
    write_restraints({"A": constr_ab_bound_dec, "B": constr_ant_dec}, restr_line_bound_dec_filename, tbl_filename_bound_dec, act_act_path)

    new_tbl_filename_bound = Path(pdb_datadir, f"{pdb}_ambig_bound_CDR_EpiVag.tbl")
    splt_string = f"!{os.linesep}! HADDOCK AIR restraints for 2nd partner"
    write_split_tbl(tbl_filename_bound_dec, splt_string, new_tbl_filename_bound)

    # 7.b ZDOCK decoupled restraints
    zdock_ab_restr_filename_dec = Path(pdb_datadir, f"{pdb}_zdock_ab_CDR_EpiVag.csv")
    zdock_ant_restr_filename_dec = Path(pdb_datadir, f"{pdb}_zdock_ant_CDR_EpiVag.csv")
    zdock_af2ant_restr_filename_dec = Path(pdb_datadir, f"{pdb}_zdock_af2ant_CDR_EpiVag.csv")
    write_zdock_restraints({"A": constr_ab_dec}, zdock_ab_restr_filename_dec)
    write_zdock_restraints({"B": constr_ant_dec}, zdock_ant_restr_filename_dec)
    write_zdock_restraints({"B": constr_ant_af2_dec}, zdock_af2ant_restr_filename_dec)

    # bonus: CDR (active)- surface residues (passive) restraints
    #af2_exposed_residues = get_solvent_exposed_residues(out_ant_af2_file, 0.4)
    #log.info(f"exposed residues on af2 structure are {af2_exposed_residues}")
    #
    #uniq_epi_af2_resids = list(set(constr_ant_af2))
    #overlap = [el for el in uniq_epi_af2_resids if el in af2_exposed_residues]
    #log.info(f"overlap residues with Para-Epi {overlap}, a fraction of {len(overlap)/len(uniq_epi_af2_resids)}")
    #
    #tbl_filename_cdr_asa = Path(pdb_datadir, f"{pdb}_ambig_CDR_ASA_act-pas.tbl")
    #restr_line_asa_filename = Path(pdb_datadir, f"{pdb}_cdr-asa_restraint_line.txt")
    #write_restraints({"A": constr_ab_dec, "B": af2_exposed_residues}, restr_line_asa_filename, tbl_filename_cdr_asa, act_act_path)
    #write_split_tbl(tbl_filename_cdr_asa, splt_string, tbl_filename_cdr_asa)

    # af2 decoupled restraints
    restr_line_af2_dec_filename = Path(pdb_datadir, f"{pdb}_af2_decoupled_restraint_line.txt")
    tbl_af2_filename_dec = Path(pdb_datadir, f"{pdb}_ambig_af2_CDR_EpiVag_act-act.tbl")
    write_restraints({"A": constr_ab_dec, "B": constr_ant_af2_dec}, restr_line_af2_dec_filename, tbl_af2_filename_dec, act_act_path)

    new_tbl_af2_filename = Path(pdb_datadir, f"{pdb}_ambig_af2_CDR_EpiVag.tbl")
    splt_string = f"!{os.linesep}! HADDOCK AIR restraints for 2nd partner"
    write_split_tbl(tbl_af2_filename_dec, splt_string, new_tbl_af2_filename)

    # 8 unambig restraints
    unambig_ant_fname = Path(pdb_datadir, f"{pdb}_unambig_ant.tbl")
    write_unambig_tbl(out_ant_file, unambig_ant_fname)
    unambig_ant_af2_fname = Path(pdb_datadir, f"{pdb}_unambig_af2_ant.tbl") # the two tbls can be different
    write_unambig_tbl(out_ant_af2_file, unambig_ant_af2_fname)

    unambig_files_std, unambig_files_af2 = {}, {}
    for key in out_ab_files.keys():
        unambig_ab_fname =  Path(pdb_datadir, f"{pdb}_unambig_{key}_abody.tbl")
        unambig_fname = Path(pdb_datadir, f"{pdb}_unambig_{key}.tbl")
        unambig_af2_fname = Path(pdb_datadir, f"{pdb}_unambig_af2_{key}.tbl")
        write_unambig_tbl(out_ab_files[key], unambig_ab_fname)
        # join antibody and antigen tbls
        join_tbls([unambig_ab_fname, unambig_ant_fname], unambig_fname)
        unambig_files_std[key] = unambig_fname
        # join antibody and af2-antigen tbls
        join_tbls([unambig_ab_fname, unambig_ant_af2_fname], unambig_af2_fname)
        unambig_files_af2[f"{key}"] = unambig_af2_fname
    
    # 8.b unambig restraints for the bound antibody
    unambig_ab_bound_fname =  Path(pdb_datadir, f"{pdb}_unambig_bound_abody.tbl")
    unambig_bound_fname = Path(pdb_datadir, f"{pdb}_unambig_bound.tbl")
    unambig_bound_af2_fname = Path(pdb_datadir, f"{pdb}_unambig_af2_bound.tbl")
    write_unambig_tbl(target_ab_fname, unambig_ab_bound_fname)
    # join antibody and antigen tbls
    join_tbls([unambig_ab_bound_fname, unambig_ant_fname], unambig_bound_fname)
    #unambig_files_std[key] = unambig_bound_fname
    # join antibody and af2-antigen tbls
    join_tbls([unambig_ab_bound_fname, unambig_ant_af2_fname], unambig_bound_af2_fname)
    #unambig_files_af2[f"{key}"] = unambig_bound_af2_fname
        
    # 9 ensemble creation
    ensemble_fname = preprocess_ensemble_file(out_ab_files,
                                              pdb_datadir,
                                              f"{pdb}_ensemble.pdb",
                                              flag_dict)
    # 9.b no AF2 ensemble creation
    out_ab_files.pop("AF2")
    ensemble_fname = preprocess_ensemble_file(out_ab_files,
                                              pdb_datadir,
                                              f"{pdb}_noAF2_ensemble.pdb",
                                              flag_dict=None) # the antibodies have already been substituted at this point

    # 9.c clustered ensemble creation
    # merge af2 and igfold ensembles
    overall_ensemble_files = {}
    # update dictionary
    dict.update(overall_ensemble_files, proc_ig_files)
    dict.update(overall_ensemble_files, proc_af_files)
    dict.update(overall_ensemble_files, proc_abb2_files)
    # with ablooper
    overall_ensemble_files["ABlooper"] = out_ab_files["ABlooper"]

    cl_ensemble_files = cluster_antibodies(overall_ensemble_files, constr_ab_dec)
    
    ensemble_fname = preprocess_ensemble_file(cl_ensemble_files,
                                                pdb_datadir,
                                                f"{pdb}_clust_ensemble.pdb",
                                                flag_dict)

    # 10 copying and sedding cfg and job files.
    ex_files_path = Path(ex_files_path)
    prepare_jobs(pdb, output_dir = pdb_dir, flag_dict = flag_dict, ex_files_path=ex_files_path)

    # 11 join different unambig tbl files (take the mean distance)
    standardize_unambig_tbls(unambig_files_std.values(), Path(pdb_datadir, f"{pdb}_unambig_ens.tbl"), "A")
    standardize_unambig_tbls(unambig_files_af2.values(), Path(pdb_datadir, f"{pdb}_unambig_af2_ens.tbl"), "A")

    log.info(f"aiabs succesful for pdb {pdb}")
    elap_time = round((time.time() - start_time), 3)
    log.info(f"aiabs run took {elap_time} seconds")

if __name__ == "__main__":
    sys.exit(maincli())
