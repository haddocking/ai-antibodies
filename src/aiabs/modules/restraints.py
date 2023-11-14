import os
import sys
from Bio import Align
from pathlib import Path
from Bio.Align import substitution_matrices
import numpy as np
import subprocess
import shlex
import logging
import numpy as np
import MDAnalysis as mda
import freesasa as fs


log = logging.getLogger("aiabslog")

def get_solvent_exposed_residues(pdb_struct, threshold):
    """
    Gets the list residues with relative freesasa higher than threshold.
    """
    mdu = mda.Universe(pdb_struct)
    resids = list(mdu.residues.resids)
    struct = fs.Structure(str(pdb_struct))
    result = fs.calc(struct)
    s = result.residueAreas()

    # do calculations
    exp_residues = []
    for n in resids:
        rel_sasa = s["B"][str(n)].relativeTotal
        if rel_sasa > threshold:
            exp_residues.append(n)
    return exp_residues

def check_fastas(seq_1, seq_2, aln_fname):
    """align sequences and gives identity in output."""
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    
    alns = aligner.align(seq_1, seq_2)
    top_aln = alns[0]

    with open(aln_fname, "w") as fh:
        fh.write(str(top_aln))

    identity = (
                str(top_aln).count("|") / float(min(len(seq_1), len(seq_2)))
                ) * 100

    return identity

def run_subproc(cmd):
    """run subprocess."""
    p = subprocess.run(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = p.stdout.decode("utf-8").split(os.linesep)
    return out

def write_tbl(content, tbl_filename):
    """write content to tbl."""
    with open(tbl_filename, "w") as wfile:
        for ln in content:
            wfile.write(f"{ln}{os.linesep}")

def check_align(dict_of_fastas, output_dir, pdb=None):
    """check that all the aligments share 100% identity."""
    sequences_dict = {}
    for keyword in dict_of_fastas.keys():
        fasta_string = open(dict_of_fastas[keyword]).read().split("\n")[1:]
        sequences_dict[keyword] = "".join(fasta_string)

    alignments_list = []
    for ab1 in sequences_dict.keys():
        for ab2 in sequences_dict.keys():
            if (ab1 != ab2):
                key = f"{ab1}-{ab2}"
                if key not in alignments_list:
                    alignments_list.append(key)
                    aln_fname = Path(output_dir, f"{key}_blosum62.aln")
                    if len(sequences_dict[ab1]) != len(sequences_dict[ab2]):
                        if pdb is None:
                            raise Exception(f"sequences for key {key} have different length")
                        elif pdb in ["7k7h"]:
                            log.warning(f"sequences for key {key} have different length")
                    identity = check_fastas(sequences_dict[ab1], sequences_dict[ab2], aln_fname)
                    if identity != 100.0:
                        raise Exception(f"alignment failed for key {key}")
    return True

def get_restrs_from_pairs(cpairs_file, ab_res, ant_res):
    """Get restraints from residue-residue pairs."""
    constr_ab, constr_ant = [], []
    with open(cpairs_file, "r") as pairs_fl:
        for line in pairs_fl:
            if line.startswith("#") == False:
                contacts = line.split(":")
                # contacts[0]
                chain_ab = contacts[0].split(",")[0]
                if chain_ab not in ["H","L"]:
                    raise Exception(f" chain {chain_ab} not compatible with antibody")
                resnum = contacts[0].split(",")[1]
                res_id = chain_ab + "-" + resnum
                constr_ab.append(ab_res.index(res_id) + 1)
                # contacts[1]
                chain_antigen = contacts[1].split(",")[0]
                resnum = contacts[1].split(",")[1]
                res_id = chain_antigen + "-" + resnum
                constr_ant.append(ant_res.index(res_id) + 1)
    return constr_ab, constr_ant

def get_restrs_from_list(c_list, residues):
    """Get restraints from list."""
    constr = []
    with open(c_list, "r") as ab_fl:
        for line in ab_fl:
            if line.startswith("chain") == False: # ignore header
                chain = line.split(",")[0]
                resnum = line.split(",")[1]
                res_id = chain + "-" + resnum
                if res_id in constr:
                    raise Exception(f"potential duplicate {res_id} in constr_ab_list")
                else:
                    constr.append(residues.index(res_id) + 1)
    return constr

def write_restraints(dict_of_resids, restr_line_filename, tbl_filename, act_act_path):
    """generate restraints."""
    list_restrs = []
    for key in dict_of_resids.keys():
        list_restrs.append(f"{key} " + "+".join([str(el) for el in np.unique(dict_of_resids[key])]))
    restr_line = " ".join(list_restrs)
    
    log.info(f"overall restraint_line\n{restr_line}")

    with open(restr_line_filename, "w") as wfile:
        wfile.write(restr_line)

    cmd = f"{act_act_path} {restr_line}"
    out = run_subproc(cmd)
    write_tbl(out, tbl_filename)

def write_split_tbl(tbl_filename, splt_string, new_tbl_filename):
    first_partner_string = open(tbl_filename, "r").read()
    content = first_partner_string.split(splt_string)[0].split(os.linesep)
    write_tbl(content, new_tbl_filename)

def write_unambig_tbl(inp_pdb_f, tbl_filename):
    """write unambiguous restraints."""
    cmd = f"python2 /trinity/login/abonvin/haddock_git/haddock-tools/restrain_bodies.py {inp_pdb_f}"
    unambig_restr = run_subproc(cmd)
    write_tbl(unambig_restr, tbl_filename)

def join_tbls(list_of_tbls, tbl_filename):
    full_string = ""
    for tbl in list_of_tbls:
        full_string += open(tbl, "r").read()
    with open(tbl_filename, "w") as wfile:
        wfile.write(full_string)

def standardize_unambig_tbls(list_of_tbls, tbl_filename, chain):
    """
    when you have multiple unambig.tbl and you want to standardize the distances
    
    WARNING : if files in the ensemble are considerably different this might lead to errors
    """
    unambig_restrs = {} # dictionary of unambig restraints. assumption: every unambig file must contain the same restrained atoms
    antigen_restrs = [] # list of antigen restraints. they shoold be equivalent in all unambig files
    unambig_lines = {}
    full_string = ""
    for fl in list_of_tbls:
        lines = open(fl).read().strip().split("\n")
        for ln in lines:
            if ln != "":
                splt_line = ln.split()
                if splt_line[2] == chain:
                    restr_key = f"{splt_line[5]}{splt_line[2]}-{splt_line[13]}{splt_line[10]}"
                    if restr_key in unambig_restrs.keys():
                        unambig_restrs[restr_key].append(float(splt_line[17]))
                    else:
                        unambig_restrs[restr_key] = [float(splt_line[17])]
                        unambig_lines[restr_key] = ln
                else: # this collects antigen restraints
                    restr_key = f"{splt_line[5]}{splt_line[2]}-{splt_line[13]}{splt_line[10]}"
                    if restr_key not in antigen_restrs:
                        full_string += f"{ln}{os.linesep}"
                        antigen_restrs.append(restr_key)
    #
    for restr_key in unambig_restrs.keys():
        mean_value = np.mean(unambig_restrs[restr_key])
        deviation = [val - mean_value for val in unambig_restrs[restr_key]]
        if any(np.array(deviation) > 2.0):
            log.info(f"\nwarning : unambigous key {restr_key} has deviations > 2.0 wrt mean value")
            log.info(f"deviations from mean value are {deviation}")
        orig_line_splt = unambig_lines[restr_key].split()
        orig_line_splt[17] = f"{mean_value:.3f}"
        output_line = f"{' '.join(orig_line_splt)}{os.linesep}"
        full_string += output_line
    
    # writing
    with open(tbl_filename, "w") as wfile:
        wfile.write(full_string)

def write_zdock_restraints(dict_of_resids, zdock_restr_filename):
    """generate ZDOCK restraints."""
    with open(zdock_restr_filename, "w") as wfile:
        wfile.write(f"chain,number,residue{os.linesep}")
        for key in dict_of_resids.keys():
            for resid in np.unique(dict_of_resids[key]):
                wfile.write(f"{key},{resid},dummy{os.linesep}")