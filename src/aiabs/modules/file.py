from pathlib import Path
import os
import shutil
import logging
log = logging.getLogger("aiabslog")

def setup_directory(output_dir, pdb):
    """sets up the output dir"""
    tpath = Path(output_dir, pdb)
    if not tpath.exists():
        os.mkdir(tpath)
    data_path = Path(tpath, "data")
    if not data_path.exists():
        os.mkdir(data_path)
    aln_path = Path(tpath, "alignments")
    if not aln_path.exists():
        os.mkdir(aln_path)
    return tpath, data_path, aln_path

def copy_files(files_list, output_dir):
    copied_files = []
    for fl in files_list:
        target_path = Path(output_dir, fl.name)
        shutil.copy(fl, target_path)
        copied_files.append(target_path)
    
    return copied_files

def sed_files(files_list, to_sed, to_replace):
    for fl in files_list:
        fl_string = open(fl, "r").read()
        new_string = fl_string.replace(to_sed, to_replace)
        with open(fl, "w") as wfile:
            wfile.write(new_string)

def prepare_jobs(pdb, output_dir, flag_dict = None, ex_files_path = Path(".")):
    log.info(f"preparing jobs for {pdb}")
    log.info(f"retrieving example files from {ex_files_path}")
    cfg_files = list(ex_files_path.glob("*cfg"))
    job_files = list(ex_files_path.glob("*job"))
    if len(cfg_files) == 0 or len(job_files) == 0:
        raise Exception("no cfg or job files, something is wrong")
    target_cfg_files = copy_files(cfg_files, output_dir = output_dir)
    target_job_files = copy_files(job_files, output_dir = output_dir)

    sed_files(target_cfg_files, to_sed="PDB", to_replace=pdb)
    sed_files(target_job_files, to_sed="PDB", to_replace=pdb)
    sed_files(target_job_files, to_sed="DATADIR", to_replace=str(output_dir))

    # correcting cfg files
    flag_map = {
        "ABodyBuilder": "abodybuilder",
        "IgFold": "igfold"
    }
    if flag_dict:
        for keyword in flag_dict.keys():
            subset_cfg_files = [fl for fl in target_cfg_files if flag_map[keyword] in fl.name]
            sed_files(subset_cfg_files,
                      to_sed=f"{pdb}_{keyword}_haddock-ready.pdb",
                      to_replace=f"{pdb}_{keyword}_haddock-ready_{flag_dict[keyword]}.pdb")
