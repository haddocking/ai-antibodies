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
import re
import random
import itertools


log = logging.getLogger("aiabslog")

def calc_euclidean(i, j):
    return ((j[0]-i[0])**2 + (j[1]-i[1])**2 + (j[2]-i[2])**2)**0.5

def read_structure(pdbf, exclude=None):
    """
    Reads a PDB file and returns a list of parsed atoms
    """
    _atoms = {'CA', 'P'}  # Alpha-Carbon (Prot), Backbone Phosphorous (DNA)
    _altloc = {' ', 'A'}

    if not exclude:
        exclude = set()
    else:
        exclude = set(exclude)

    res_list = []
    with open(pdbf, 'r') as pdb_handle:
        for line in pdb_handle:
            field = line[0:4]
            if field != 'ATOM':
                continue

            aname = line[12:16].strip()
            chain = line[21] if line[21].strip() else line[72:76].strip()  # chain ID or segID
            if chain not in exclude and aname in _atoms and line[16] in _altloc:
                resi = int(line[22:26])
                coords = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                res_list.append((chain, resi, aname, coords))

    if not res_list:
        logging.critical('[!] PDB File seems empty or no CA/P atoms found: {0}'.format(pdbf))
        sys.exit(1)

    return res_list


def get_bodies(atom_lst, prot_threshold=4.0, dna_threshold=7.5):
    """
    Determines gaps in an atom list following simple distance based criteria.
    Returns continuous fragments.
    """

    bodies = []

    threshold = {'CA': prot_threshold, 'P': dna_threshold}

    body_start = 0
    i = None
    for i, atom in enumerate(atom_lst[1:], start=1):
        p_atom = atom_lst[i-1]

        chain, resi, aname, xyz = atom
        p_chain, p_resi, p_aname, p_xyz = p_atom

        if (chain == p_chain) and (aname == p_aname):  # Internal Gap
            d_xyz = calc_euclidean(xyz, p_xyz)
            if d_xyz >= threshold[aname]:
                logging.debug('[+++] (Internal) Body: {0}:{1}'.format(body_start, i-1))
                bodies.append((body_start, i-1))
                body_start = i  # Set new beginning

        elif (chain != p_chain) or (aname != p_aname):  # Different molecules/types
            logging.debug('[+++] Body: {0}:{1}'.format(body_start, i-1))
            bodies.append((body_start, i-1))
            body_start = i  # Set new beginning

    if not bodies:  # Single continuous molecule
        bodies.append((0, len(atom_lst)))
    else:
        logging.debug('[+++] Body: {0}:{1}'.format(body_start, i))
        bodies.append((body_start, i))  # Last body

    logging.info('[++] Found {0} bodies'.format(len(bodies)))

    return bodies


def build_restraints(bodies):
    """
    Generates distance restraints to maintain the relative
    orientation of the different bodies during the simulations.

    Generates two unique restraints per pair of bodies.

    Each restraint is created using two random atoms on each body
    and using their exact euclidean distance as target distance.
    """

    def pick_residues(body, max_trials=10):
        # Pick two random residues in each body
        # Make sure they are far apart from each other
        n_trials = 0
        while 1:
            try:
                res_i, res_ii = random.sample(body, 2)
            except ValueError:
                # Likely, sample size is 1
                logging.warning('[!] One-sized body found. This may lead to problems..')
                return body[0], body[0]

            logging.debug('[+++] Trial {0}: {1} & {2}'.format(n_trials, res_i, res_ii))
            if abs(res_i - res_ii) > 3:
                logging.info('[++] Picked residues {0} & {1}'.format(res_i, res_ii))
                return res_i, res_ii
            n_trials += 1
            if n_trials == max_trials:
                msg = '[!] Could not pick two unique distant residues in body after {0} tries'
                logging.info(msg.format(max_trials))
                return res_i, res_ii

    restraints = []

    n_bodies = range(len(bodies))
    combinations = itertools.combinations(n_bodies, 2)

    for pair_bodies in combinations:
        body_i, body_j = pair_bodies
        logging.debug('[+++] Restraining body {0} to body {1}'.format(body_i, body_j))

        st_body_i, en_body_i = bodies[body_i]
        st_body_j, en_body_j = bodies[body_j]
        res_i, res_ii = pick_residues(range(st_body_i, en_body_i+1))
        res_j, res_jj = pick_residues(range(st_body_j, en_body_j+1))

        logging.info('[++] Created restraint: {0}:{1} <--> {2}:{3}'.format(body_i, res_i, body_j, res_j))
        restraints.append((res_i, res_j))
        logging.info('[++] Created restraint: {0}:{1} <--> {2}:{3}'.format(body_i, res_ii, body_j, res_jj))
        restraints.append((res_ii, res_jj))

    return restraints


def generate_tbl(atom_lst, restraints):
    """
    Makes a list of TBL-formatted restraints.

    Parameters
    ----------
    atom_lst : list
        List of atoms in the form (chain, resi, aname, coords)
    
    restraints : list
        List of restraints in the form (res_i, res_j)
    """
    tbl_content = []
    for r in restraints:
        i, j = r
        atom_i, atom_j = atom_lst[i], atom_lst[j]
        dist_ij = calc_euclidean(atom_i[3], atom_j[3])

        tbl = "assign (segid {0[0]} and resi {0[1]} and name {0[2]}) ".format(atom_i)
        tbl += "(segid {0[0]} and resi {0[1]} and name {0[2]}) ".format(atom_j)
        tbl += "{0:3.3f} 0.0 0.0".format(dist_ij)
        tbl_content.append(tbl)
    return tbl_content


def check_parenthesis(file):
    open_parenthesis = re.compile('[(]')
    close_parenthesis = re.compile('[)]')
    quotation_marks = re.compile('[\"]')
    opened = 0
    closed = 0
    quote = 0
    for _match in open_parenthesis.finditer(file):
        opened += 1
    for _match in close_parenthesis.finditer(file):
        closed += 1
    for _match in quotation_marks.finditer(file):
        quote += 1
    if opened != closed:
        raise Exception("Problem with TBL file parentheses ({:d} opening for {:d} "
                        "closing parentheses)".format(opened, closed))
    if quote % 2 != 0:
        raise Exception("Problem with TBL file, odd number of quotation marks "
                        "({:d} quotation marks)".format(quote))

def restrain_bodies(structure, exclude=None, verbose=0):  # noqa: E501
    """Create distance restraints to lock several chains together.
    
    Parameters
    ----------
    structure : str
        The PDB structure to be restrained.

    exclude : str
        Chains to exclude from the calculation.
    
    verbose : int
        Tune verbosity of the output.
    """
    if verbose== 1:
        logging.basicConfig(level=logging.INFO)
    elif verbose > 1:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.WARNING)

    # Main logic
    atom_lst = read_structure(structure, exclude=exclude)
    bodies = get_bodies(atom_lst)
    restraints = build_restraints(bodies)
    tbl_content = generate_tbl(atom_lst, restraints)
    return tbl_content

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
                    len_ab1 = len(sequences_dict[ab1])
                    len_ab2 = len(sequences_dict[ab2])
                    if len_ab1 != len_ab2:
                        if pdb is None:
                            exc_message = f"sequences for key {key} have different lengths ({len_ab1} vs {len_ab2})"
                            if "target" in key:
                                log.warning(exc_message)
                            else:
                                raise Exception(exc_message)
                        elif pdb in ["7k7h"]:
                            log.warning(f"sequences for key {key} have different length")
                    identity = check_fastas(sequences_dict[ab1], sequences_dict[ab2], aln_fname)
                    if identity != 100.0:
                        raise Exception(f"alignment failed for key {key}, identity is {identity}")
    return True

def get_restrs_from_pairs(cpairs_file, ab_res, ant_res, chain_mapping):
    """Get restraints from residue-residue pairs."""
    constr_ab, constr_ant = [], []
    nwarns = 0
    print(f"ab_res {ab_res}")
    with open(cpairs_file, "r") as pairs_fl:
        for line in pairs_fl:
            if line.startswith("#") == False:
                contacts = line.split(":")
                # contacts[0]
                chain_ab = contacts[0].split(",")[0]
                # 22/5/2024: this is not true anymore
                #if chain_ab not in ["H","L"]:
                #    raise Exception(f" chain {chain_ab} not compatible with antibody")
                resnum = contacts[0].split(",")[1]
                res_id = chain_mapping[chain_ab] + "-" + resnum
                # 26/5/2024: need to check the last character of the res_id
                # if it is a string, then we should not count it
                res_id_splt = res_id.split("-")[1]
                last_character = res_id_splt[-1]
                if last_character.isalpha():
                    resid = int(res_id_splt[:-1])
                    #continue
                else:
                    resid = int(res_id_splt)
                # now the check
                if res_id not in ab_res and resid > 130:
                    log.warning(f"residue {res_id} not in ab_res")
                    nwarns += 1
                    continue
                constr_ab.append(ab_res.index(res_id) + 1)
                # contacts[1]
                chain_antigen = contacts[1].split(",")[0]
                resnum = contacts[1].split(",")[1]
                res_id = chain_antigen + "-" + resnum
                constr_ant.append(ant_res.index(res_id) + 1)
                print(f"line {line} -> {constr_ab[-1]} {constr_ant[-1]}")
    # if everything went smoothly the number of warnings should be very small
    log.info(f"number of warnings in get_restrs_from_pairs {nwarns}")
    if nwarns/len(constr_ab) > 0.1:
        raise Exception(f"too many warnings {nwarns} in get_restrs_from_pairs")
    return constr_ab, constr_ant

def get_restrs_from_list(c_list, residues, chain_mapping=None):
    """Get restraints from list."""
    constr = []
    with open(c_list, "r") as ab_fl:
        for line in ab_fl:
            if line.startswith("chain") == False: # ignore header
                chain = line.split(",")[0]
                resnum = line.split(",")[1]
                if chain_mapping:
                    res_id = chain_mapping[chain] + "-" + resnum
                else:
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
    #cmd = f"python2 /trinity/login/abonvin/haddock_git/haddock-tools/restrain_bodies.py {inp_pdb_f}"
    #unambig_restr = run_subproc(cmd)
    unambig_restr = restrain_bodies(inp_pdb_f)
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