from pdbtools.pdb_tofasta import pdb_to_fasta
from pathlib import Path
import shutil
from pdbtools.pdb_splitchain import split_chain
from pdbtools.pdb_merge import concatenate_files
from pdbtools.pdb_tidy import tidy_pdbfile
from pdbtools.pdb_chain import alter_chain
from pdbtools.pdb_reres import renumber_residues
from pdbtools.pdb_reatom import renumber_atom_serials
from pdbtools.pdb_delhetatm import remove_hetatm
from pdbtools.pdb_mkensemble import make_ensemble
import logging

log = logging.getLogger("aiabslog")

def generate_fasta(inp_pdb_f, output_dir):
    
    log.info(f"file {inp_pdb_f}")
    fasta_fname = Path(output_dir, f"{inp_pdb_f.stem}.fasta")
    with open(inp_pdb_f, "r") as pdb_fh:
        with open(fasta_fname, "w") as wfile:
            fasta = pdb_to_fasta(pdb_fh, multi=None)
            wfile.write(''.join(fasta))

    return fasta_fname

def remove_occ_pdb(inp_filename, out_filename):
    # removing occupancies
    with open(out_filename, "w") as wfile:
        with open(inp_filename, "r") as rfile:
            for line in rfile:
                wfile.write(f"{line[:26]} {line[27:]}")
    return

def tidy_pdb(inp_pdb_f):
    """
    Tidy PDB file.

    Parameters
    ----------
    inp_pdb_f : Path
        Path to PDB file.

    Returns
    -------
    Path
        Path to PDB file.
    """
    # log.debug("Tidying PDB file")
    out_pdb_fname = Path(f"{inp_pdb_f.stem}-tidy.pdb")
    with open(inp_pdb_f, "r") as pdb_fh:
        with open(out_pdb_fname, "w") as f:
            for line in tidy_pdbfile(pdb_fh):
                f.write(line)
    return out_pdb_fname


def reatom_pdb(inp_pdb_f, atomid=1):
    out_pdb_fname = Path(f"{inp_pdb_f.stem}-reres{atomid}.pdb")
    with open(inp_pdb_f, "r") as pdb_fh:
        with open(out_pdb_fname, "w") as f:
            for line in renumber_atom_serials(pdb_fh, atomid):
                f.write(line)
    return out_pdb_fname

def reres_pdb(inp_pdb_f, resid=1):

    out_pdb_fname = Path(f"{inp_pdb_f.stem}-reres{resid}.pdb")
    with open(inp_pdb_f, "r") as pdb_fh:
        with open(out_pdb_fname, "w") as f:
            for line in renumber_residues(pdb_fh, resid):
                f.write(line)
    return out_pdb_fname

def chain_pdb(inp_pdb_f, chain_id):
    out_pdb_fname = Path(f"{inp_pdb_f.stem}-chain{chain_id}.pdb")
    with open(inp_pdb_f, "r") as pdb_fh:
        with open(out_pdb_fname, "w") as f:
            for line in alter_chain(pdb_fh, chain_id):
                f.write(line)
    return out_pdb_fname

def merge_pdb(inp_pdb_fs):
    out_pdb_fname = Path(f"{inp_pdb_fs[0].stem}-merge.pdb")
    pdb_fhs = [open(pdb, "r") for pdb in inp_pdb_fs]
    merged_pdb = concatenate_files(pdb_fhs)
    with open(out_pdb_fname, "w") as f:
            for line in merged_pdb:
                f.write(line)
    return out_pdb_fname

def mkensemble_pdb(inp_pdb_fs):
    out_pdb_fname = Path(f"ensemble.pdb")
    ensemble_pdb = make_ensemble(inp_pdb_fs)
    with open(out_pdb_fname, "w") as f:
            for line in ensemble_pdb:
                f.write(line)
    return out_pdb_fname

def split_pdb(inp_pdb_f):
    with open(inp_pdb_f, "r") as pdb_fh:
        split_chain(pdb_fh)


def delhetatom_pdb(inp_pdb_f):
    out_pdb_fname = Path(f"{inp_pdb_f.stem}-delhetatom.pdb")
    with open(inp_pdb_f, "r") as pdb_fh:
        with open(out_pdb_fname, "w") as f:
            for line in remove_hetatm(pdb_fh):
                f.write(line)
    return out_pdb_fname
    
def generate_target(target_file, output_dir, pdb):
    split_pdb(target_file)
    # moving away antibody files
    ab_files = list(Path(".").glob("*true_complex_[H-L]*"))
    ab_files = sorted(ab_files)
    log.info(f"ab_files {ab_files}")

    # merging antibody chains
    merged_pdb = merge_pdb(ab_files)
    chained_pdb = chain_pdb(merged_pdb, "A")
    reresnum_pdb = reres_pdb(chained_pdb, 1)
    tidyed_pdb = tidy_pdb(reresnum_pdb)
    target_ab_fname = Path(output_dir, f"{pdb}_targetab.pdb")
    remove_occ_pdb(tidyed_pdb, target_ab_fname)
    #shutil.move(tidyed_pdb, target_ab_fname)

    for fl in ab_files:
        fl.unlink()
    merged_pdb.unlink()
    chained_pdb.unlink()
    reresnum_pdb.unlink()
    tidyed_pdb.unlink()
    
    # now the antigen!
    ant_files = list(Path(".").glob("*true_complex*"))
    ant_files = sorted(ant_files)
    log.info(f"antigen files {ant_files}")
    merged_pdb = merge_pdb(ant_files)
    chained_pdb = chain_pdb(merged_pdb, "B")
    reresnum_pdb = reres_pdb(chained_pdb, 1)
    tidyed_pdb = tidy_pdb(reresnum_pdb)
    target_ant_fname = Path(output_dir, f"{pdb}_targetant.pdb")
    remove_occ_pdb(tidyed_pdb, target_ant_fname)
    #shutil.move(tidyed_pdb, target_ant_fname)
    for fl in ant_files:
        fl.unlink()
    merged_pdb.unlink()
    chained_pdb.unlink()
    reresnum_pdb.unlink()
    tidyed_pdb.unlink()
    # merging them together
    merged_pdb = merge_pdb([target_ab_fname, target_ant_fname])
    reatomnum_pdb = reatom_pdb(merged_pdb, 1)
    tidyed_pdb = tidy_pdb(reatomnum_pdb)
    target_fname = Path(output_dir, f"{pdb}_target.pdb")
    #shutil.move(tidyed_pdb, target_fname)
    remove_occ_pdb(tidyed_pdb, target_fname)
    merged_pdb.unlink()
    reatomnum_pdb.unlink()
    tidyed_pdb.unlink()

    return target_ab_fname, target_ant_fname, target_fname


def preprocess_ab_file(inp_pdb_f, output_dir, out_fname):
    delheted_pdb = delhetatom_pdb(inp_pdb_f)
    chained_pdb = chain_pdb(delheted_pdb, "A")
    reresnum_pdb = reres_pdb(chained_pdb, 1)
    reatomnum_pdb = reatom_pdb(reresnum_pdb, 1)
    tidyed_pdb = tidy_pdb(reatomnum_pdb)
    out_fname = Path(output_dir, out_fname)
    # removing occupancies
    #with open(out_fname, "w") as wfile:
    #    with open(tidyed_pdb, "r") as rfile:
    #        for line in rfile:
    #            wfile.write(f"{line[:26]} {line[27:]}")
    remove_occ_pdb(tidyed_pdb, out_fname)
    
    # unlinking
    delheted_pdb.unlink()
    chained_pdb.unlink()
    reresnum_pdb.unlink()
    reatomnum_pdb.unlink()
    tidyed_pdb.unlink()
    return out_fname

def preprocess_ant_file(inp_pdb_f, output_dir, out_fname):
    delheted_pdb = delhetatom_pdb(inp_pdb_f)
    chained_pdb = chain_pdb(delheted_pdb, "B")
    reresnum_pdb = reres_pdb(chained_pdb, 1)
    reatomnum_pdb = reatom_pdb(reresnum_pdb, 1)
    tidyed_pdb = tidy_pdb(reatomnum_pdb)
    out_fname = Path(output_dir, out_fname)
    shutil.move(tidyed_pdb, out_fname)
    
    # unlinking
    delheted_pdb.unlink()
    chained_pdb.unlink()
    reresnum_pdb.unlink()
    reatomnum_pdb.unlink()

    return out_fname

def preprocess_ensemble_file(inp_pdb_dict, output_dir, out_fname, flag_dict=None):
    """
    Parameters 

    inp_pdb_dict : dict
        dictionary of input pdbs
    output_dir : Path or str
        output string
    out_fname : str
        ensemble filename
    flag_dict : dict or None
        
    Returns
    -------
    out_fname : Path or str
        ensemble output filename
    """
    if flag_dict:
        for keyword in flag_dict:
            log.info(f"using {flag_dict[keyword]} for {keyword} pdb")
            if keyword not in inp_pdb_dict.keys():
                log.warning(f"keyword {keyword} not in {inp_pdb_dict.keys()}")
                continue
            # remove key from input pdbs
            ori_name = inp_pdb_dict[keyword]
            new_name = Path(ori_name.parent, f"{ori_name.stem}_{flag_dict[keyword]}{ori_name.suffix}")
            log.info(f"new {keyword} filename {new_name}")
            inp_pdb_dict[keyword] = new_name

    inp_pdb_fs = inp_pdb_dict.values()

    ensemble_pdb = mkensemble_pdb(inp_pdb_fs)
    tidyed_pdb = tidy_pdb(ensemble_pdb)

    out_fname = Path(output_dir, out_fname)
    shutil.move(tidyed_pdb, out_fname)

    # unlinking
    ensemble_pdb.unlink()
    
    return out_fname

def get_resids(inp_pdb_f):
    pdb_resids = []
    with open(inp_pdb_f, "r") as ab_file:
        for ln in ab_file:
            if ln[13:15] == "CA":
                res_id = ln[21] + "-" + ln[22:27].strip()
                pdb_resids.append(res_id)
    return pdb_resids