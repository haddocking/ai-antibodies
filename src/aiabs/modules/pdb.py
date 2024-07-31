from pdbtools.pdb_tofasta import pdb_to_fasta
from pathlib import Path
import shutil
from pdbtools.pdb_splitchain import split_chain
from pdbtools.pdb_merge import concatenate_files
from pdbtools.pdb_tidy import tidy_pdbfile
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

def _remove_altloc(lines):
    # the altloc ID is removed in processed altloc lines
    for line_num, line in lines:
        yield (line_num, line[:16] + " " + line[17:])


def _flush(register, option, others):
    """
    Processes the collected atoms according to the selaltloc option.
    """
    lines_to_yield = []
    select_by_occupancy = option is None

    atom_lines = ("ATOM", "HETATM")

    # anisou lines are treated specially
    anisou_lines = ("ANISOU",)

    for resnum, atomnames in register.items():
        for atomname, altlocs in atomnames.items():
            if select_by_occupancy:
                # gathers all alternative locations for the atom
                all_lines = []
                for altloc, lines in altlocs.items():
                    all_lines.extend(lines)

                # identifies the highest occupancy combining dictionary
                # and sorting
                new = {}
                for line_number, line in all_lines:
                    if line.startswith(atom_lines):
                        occupancy_number = line[54:60]
                        list_ = new.setdefault(occupancy_number, [])
                        list_.append((line_number, line))

                    # assumes ANISOU succeed the respective ATOM line
                    elif line.startswith(anisou_lines):
                        list_.append((line_number, line))

                # sort keys by occupancy
                keys_ = sorted(
                    new.keys(), key=lambda x: float(x.strip()), reverse=True
                )

                these_atom_lines = new[keys_[0]]

                # always yield the first line
                lines_to_yield.extend(_remove_altloc(these_atom_lines[0:1]))

                del all_lines, new

            # selected by option:
            else:
                if option in altlocs:
                    # selects the option, that's it
                    lines_to_yield.extend(_remove_altloc(altlocs[option]))

                else:
                    # if the option does not exist, add all altlocs
                    for altloc, lines in altlocs.items():
                        lines_to_yield.extend(lines)

    # add comments
    lines_to_yield.extend(others)

    # lines are sorted to the line number so that the output is sorted
    # the same way as in the input PDB
    lines_to_yield.sort(key=lambda x: x[0])

    # the line number is ignored, only the line is yield
    for line_number, line in lines_to_yield:
        yield line


def select_by_occupancy(fhandle, option=None):
    """
    Selects altloc labels for the entire PDB file.

    Parameters
    ----------
    fhandle : an iterable giving the PDB file line-by-line.

    option : str or `None`.
        The alternative location identifier to select. By default
        (`None`) selects the alternative location with highest
        occupancy. In this case, if the different alternative locations
        have the same occupancy, selects the one that comes first.
        Selecting by highest occupancy removes all altloc labels for all
        atoms. Provide an option (e.g. 'A') to select only atoms with
        altloc label `A`. If you select `A` and an atom has conformers
        with altlocs `B` and `C`, both B and C will be kept in the
        output. Despite not an official format, many times alternative
        locations are identified by a blank character ' ' (space), and a
        [A-Z] character.  In these cases, to select the alternative
        location identified by a blank character give `option=' '`.

    Returns
    -------
    generator
        A generator object. To exhaust the generator, that is, to
        process the PDB file (or PDB lines), convert it to a list.

        >>> from pdbtools.pdb_selaltloc import run
        >>> with('input.pdb', 'r') as fin:
        >>>     processed_lines = list(run(fin))

        For more example see:

        >>> import pdbtools
        >>> help(pdbtools)
    """
    records = ("ATOM", "HETATM", "ANISOU")
    terminators = ("TER", "END", "CONECT", "END", "ENDMDL", "MODEL")

    # register atom information
    register = dict()

    # register comment lines
    others = []

    # register current chain
    chain = None
    prev_chain = None

    # keep record of the line number. This will be used to sort lines
    # after selecting the desired alternative location
    nline = 0

    # the loop will collect information on the different atoms
    # throughout the PDB file until a new chain or any terminal line is
    # found. At that point, the collected information is flushed because
    # all altlocs for that block have been defined.
    for line in fhandle:
        nline += 1

        if line.startswith(records):
            # here resnum + insertion code are taken to identify
            # different residues
            resnum = line[22:27]
            atomname = line[12:16]
            altloc = line[16]
            chain = line[21:22]

            # flush lines because we enter a new chain
            if chain != prev_chain:
                # the "yield from" statement is avoided to keep
                # compatibility with Python 2.7
                for _line in _flush(register, option, others):
                    yield _line

                # Python 2.7 compatibility. Do not use .clear() method
                # restart help variables
                del register, others
                register, others = dict(), []

            # organizes information hierarchically
            resnum_d = register.setdefault(resnum, {})
            atomname_d = resnum_d.setdefault(atomname, {})
            altloc_d = atomname_d.setdefault(altloc, [])

            # adds info to dictionary
            altloc_d.append((nline, line))

        # flush information because we reached the end of a block
        elif line.startswith(terminators):
            for _line in _flush(register, option, others):
                yield _line

            del register, others
            register, others = dict(), []

            yield line  # yield the current line after flush

        else:
            # append comments to flush list
            # The reason to add comments to a list instead of yielding
            # them directly is to cover the possibility of having
            # comments in the middle of the PDB file. Obviously is this
            # extremely unlikely. But just in case...
            others.append((nline, line))

        prev_chain = chain

    # at the end of the PDB, flush the remaining lines
    for _line in _flush(register, option, others):
        yield _line

def occ_pdb(inp_pdb_f):
    """
    Select residues with highest occupancy.

    Parameters
    ----------
    inp_pdb_f : Path
        Path to PDB file.

    Returns
    -------
    out_pdb_fname : Path
        Path to PDB file.
    """
    # log.debug("Selecting residues with highest occupancy")
    out_pdb_fname = Path(f"{inp_pdb_f.stem}-occ.pdb")
    with open(inp_pdb_f, "r") as pdb_fh:
        with open(out_pdb_fname, "w") as f:
            for line in select_by_occupancy(pdb_fh):
                f.write(line)
    return out_pdb_fname


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
            for line in pdb_fh:
                f.write(f"{line[:21]}{chain_id}{line[22:]}")
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
    
def generate_target(target_file, output_dir, pdb, chains):
    split_pdb(target_file)
    # moving away antibody files
    print(f"chains {chains}")
    ab_files = [
        Path(f"{pdb}_true_complex_{chains[0]}.pdb"),
        Path(f"{pdb}_true_complex_{chains[1]}.pdb"),
    ]
    # 16 may 2024
    log.info(f"ab_files {ab_files}")
    #antibody_bound_resids = get_resids(ab_files[0])#.extend(get_resids(ab_files[1]))
    
    # add to this object the ab_files[1] residues
    #antibody_bound_resids.extend(get_resids(ab_files[1]))

    # merging antibody chains
    occ_pdb_first = occ_pdb(ab_files[0])
    occ_pdb_second = occ_pdb(ab_files[1])
    merged_pdb = merge_pdb([occ_pdb_first, occ_pdb_second])
    antibody_bound_resids = get_resids(merged_pdb)
    print(f"antibody_bound_resids {antibody_bound_resids}")
    #import os
    #print(f"merged_pdb {merged_pdb} ls {os.listdir()}")
    chained_pdb = chain_pdb(merged_pdb, "A")
    tidyed_pdbch = tidy_pdb(chained_pdb)
    reresnum_pdb = reres_pdb(tidyed_pdbch, 1)
    tidyed_pdb = tidy_pdb(reresnum_pdb)
    target_ab_fname = Path(output_dir, f"{pdb}_targetab.pdb")
    remove_occ_pdb(tidyed_pdb, target_ab_fname)
    # shutil.move(chained_pdb, "chain.pdb")
    # shutil.move(reresnum_pdb, "reres.pdb")
    # shutil.move(tidyed_pdb, "tidy.pdb")
    # shutil.move(tidyed_pdbch, "tidych.pdb")
    # shutil.move(occ_pdb_first, "occ1.pdb")
    # shutil.move(occ_pdb_second, "occ2.pdb")
    # shutil.move(merged_pdb, "merged.pdb")

    for fl in ab_files:
        fl.unlink()
    occ_pdb_first.unlink()
    occ_pdb_second.unlink()
    merged_pdb.unlink()
    chained_pdb.unlink()
    reresnum_pdb.unlink()
    tidyed_pdb.unlink()
    tidyed_pdbch.unlink()
    
    # now the antigen!
    ant_files = list(Path(".").glob("*true_complex*"))
    print(f"ant_files {ant_files}")
    ant_files = sorted(ant_files)
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

    return target_ab_fname, target_ant_fname, target_fname, antibody_bound_resids


def preprocess_ab_file(inp_pdb_f, output_dir, out_fname):
    delheted_pdb = delhetatom_pdb(inp_pdb_f)
    chained_pdb = chain_pdb(delheted_pdb, "A")
    reresnum_pdb = reres_pdb(chained_pdb, 1)
    reatomnum_pdb = reatom_pdb(reresnum_pdb, 1)
    tidyed_pdb = tidy_pdb(reatomnum_pdb)
    out_fname = Path(output_dir, out_fname)
    # removing occupancies
    remove_occ_pdb(tidyed_pdb, out_fname)
    
    # unlinking
    delheted_pdb.unlink()
    chained_pdb.unlink()
    reresnum_pdb.unlink()
    reatomnum_pdb.unlink()
    tidyed_pdb.unlink()
    return out_fname

def preprocess_ant_file(inp_pdb_f, output_dir, out_fname):
    #occ_pdb = Path(f"{inp_pdb_f.name}-noocc.pdb")
    noocc_pdb = occ_pdb(inp_pdb_f)
    delheted_pdb = delhetatom_pdb(noocc_pdb)
    # 28/5/2024: we need to get the residues that are in the antigen
    ant_residues = get_resids(delheted_pdb)
    chained_pdb = chain_pdb(delheted_pdb, "B")
    reresnum_pdb = reres_pdb(chained_pdb, 1)
    reatomnum_pdb = reatom_pdb(reresnum_pdb, 1)
    tidyed_pdb = tidy_pdb(reatomnum_pdb)
    out_fname = Path(output_dir, out_fname)
    shutil.move(tidyed_pdb, out_fname)
    
    # unlinking
    noocc_pdb.unlink()
    delheted_pdb.unlink()
    chained_pdb.unlink()
    reresnum_pdb.unlink()
    reatomnum_pdb.unlink()

    return out_fname, ant_residues

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
