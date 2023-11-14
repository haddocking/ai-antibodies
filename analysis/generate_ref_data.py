import pandas as pd
from pathlib import Path
import numpy as np
import MDAnalysis as mda
from functions import get_epitope_and_paratope

#comparison to model accuracy
ref_data = {}
data = []
with open(Path("data", "cdr_rmsd_all_models_updated.csv")) as rfile:
    for ln in rfile:
        graffe = ln.split("{")[1:]
        for ab in graffe:
            splt_comma = ab.split(",")[:-1]
            cur_dict = {}
            for el in splt_comma:
                keyval = el.split(":")
                key = keyval[0].replace("'", "").replace(" ", "").replace("\\", "").replace("}", "").replace('"', "")
                value = keyval[1].replace("'", "").replace(" ", "").replace("\\", "").replace("}", "").replace('"', "")
                cur_dict[key] = value
            pdb = cur_dict["pdb"]
            print(pdb, cur_dict)
            if pdb not in ["7kpj", "7kn4"]:
                if cur_dict["method"] != "ABB":
                    data_list = [cur_dict["pdb"], cur_dict["method"], cur_dict["cdrh1"], cur_dict["cdrh2"], cur_dict["cdrh3"], cur_dict["cdrl1"], cur_dict["cdrl2"], cur_dict["cdrl3"]]
                    data.append(data_list)
print("saving reference data")
ref_data_df = pd.DataFrame(data)
ref_data_df.columns = ["pdb", "method", "cdrh1", "cdrh2", "cdrh3", "cdrl1", "cdrl2", "cdrl3"]
ref_data_df.to_csv(Path("data", "ref_data.csv"),float_format='%.3f',sep=" ")

# plddt categories
af2_values = pd.read_csv(Path("data", "AF2_antigen_rmsd.csv"), sep=",")
assert af2_values.shape[0] == 79
print(f"mean full antigen RMSD {np.mean(af2_values['rmsd'])} median {np.median(af2_values['rmsd'])}")
print(f"mean epitope RMSD {np.mean(af2_values['rmsd_epitope'])} median {np.median(af2_values['rmsd_epitope'])}")
print(f"number of structures with epi_rmsd < 2.0 {np.sum(af2_values['rmsd_epitope'] < 2.0)}")
print(f"number of structures with overall RMSD < 2.0 {np.sum(af2_values['rmsd'] < 2.0)}")


bench_folder = Path("..", "benchmark_haddock_10_Jan_2023")
ave_plddt_list = []
ave_plddt_histo_epi = {"bad": [], "good" : [], "high": []}
ave_plddt_histo_full = {"bad": [], "good" : [], "high": []}
for pdb in np.unique(af2_values["pdb"]):
    full_rmsd = af2_values[af2_values["pdb"] == pdb]["rmsd"].values[0]
    epi_rmsd = af2_values[af2_values["pdb"] == pdb]["rmsd_epitope"].values[0]
    if full_rmsd > 5.0:
        cat_full = "bad"
    elif full_rmsd < 2.0:
        cat_full = "high"
    else:
        cat_full = "good"
    if epi_rmsd > 5.0:
        cat_epi = "bad"
    elif epi_rmsd < 2.0:
        cat_epi = "high"
    else:
        cat_epi = "good"
    # now plddt
    constr_file = Path(bench_folder, pdb, f"{pdb}_af2_constraint_pairs.txt")
    af2_file = Path(bench_folder, pdb, f"{pdb}_AF2_antigen_model.pdb")
    # getting epitope and paratope
    epitope, paratope = get_epitope_and_paratope(constr_file)
    # processing epitope
    chain = list(epitope.keys())[0]
    mdu = mda.Universe(af2_file)
    selection_string = f"chainID {chain} and name CA and resid {' '.join([str(el) for el in epitope[chain]])}"
    full_selection = f"chainID {chain} and name CA"
    epitope_atoms = mdu.select_atoms(selection_string)
    full_atoms = mdu.select_atoms(full_selection)
    # getting average plddt
    ave_plddt_epi = np.mean(epitope_atoms.tempfactors)
    ave_plddt_full = np.mean(full_atoms.tempfactors)
    if ave_plddt_epi < 70.0:
        cat_epi_plddt = "bad"
    elif ave_plddt_epi > 90.0:
        cat_epi_plddt = "high"
    else:
        cat_epi_plddt = "good"
    
    if ave_plddt_full < 70.0:
        cat_full_plddt = "bad"
    elif ave_plddt_full > 90.0:
        cat_full_plddt = "high"
    else:
        cat_full_plddt = "good"

    ave_plddt_list.append([pdb, ave_plddt_epi, ave_plddt_full, cat_epi, cat_full, cat_epi_plddt, cat_full_plddt])
    
    print(f"pdb {pdb} rmsd_epi {epi_rmsd} cat_epi {cat_epi} ave_plddt_epi {ave_plddt_epi}, cat_epi_plddt {cat_epi_plddt}")

ave_plddt_df = pd.DataFrame(ave_plddt_list)
print(ave_plddt_df.shape)

ave_plddt_df.columns = ["pdb", "epi-ave-plddt", "full-ave-plddt", "epi-cat", "full-cat", "epi-cat-plddt", "full-cat-plddt"]
print(print(f"good epitope structures {ave_plddt_df.where(ave_plddt_df['epi-cat'] == 'good').dropna()}"))

af2_values_full = pd.merge(af2_values, ave_plddt_df, on='pdb')
af2_values_full.to_csv(Path("data", "AF2_antigen_rmsd_plddt.csv"), sep="\t")

# plddt categories paratope
ave_plddt_paratope_list = []
ave_plddt_paratope_histo = {"bad": [], "good" : [], "high": []}

ave_abb2conf_paratope_list = []
ave_abb2conf_paratope_histo = {"bad": [], "good" : [], "high": []}

for pdb in np.unique(af2_values["pdb"]):
    if pdb not in ["7kpj", "7kn4"]:
        constr_file = Path(bench_folder, pdb, f"{pdb}_af2_constraint_pairs.txt")
        # loading files
        af2_file = Path(bench_folder, pdb, f"AF2_{pdb}_antibody_model_imgt.pdb")
        abb2_file = Path(bench_folder, pdb, f"ABodyBuilder2_{pdb}_antibody_model.pdb")
        mdu = mda.Universe(af2_file)
        mdu_abb2 = mda.Universe(abb2_file)
        # get epitope and paratope
        epitope, paratope = get_epitope_and_paratope(constr_file)
        print(f"paratope {paratope}")
        # define selection
        sel_strings = []
        chains = list(paratope.keys())
        for chain in chains:
            selection_string = f"(chainID {chain} and name CA and resid {' '.join([str(el) for el in paratope[chain]])})"
            sel_strings.append(selection_string)
        selection_string = " or ".join(sel_strings)
        print("full sel_string", selection_string)
        # select atoms
        paratope_atoms = mdu.select_atoms(selection_string)
        paratope_atoms_abb2 = mdu_abb2.select_atoms(selection_string)
        # assert same number of atoms
        len_paratope = sum([len(chain) for chain in paratope.values()])
        assert len_paratope == len(paratope_atoms)
        assert len_paratope == len(paratope_atoms_abb2)

        ave_plddt = np.mean(paratope_atoms.tempfactors)
        ave_abb2conf = np.mean(paratope_atoms_abb2.tempfactors)
        
        ave_plddt_paratope_list.append([pdb, ave_plddt])
        ave_abb2conf_paratope_list.append([pdb, ave_abb2conf])
        if ave_plddt < 70.0:
            cat = "bad"
        elif ave_plddt > 90.0:
            cat = "high"
        else:
            cat = "good"
        ave_plddt_paratope_histo[cat].append(pdb)
        
        if ave_abb2conf > 1.0:
            cat_abb2 = "bad"
        elif ave_abb2conf < 0.5:
            cat_abb2 = "high"
        else:
            cat_abb2 = "good"
        ave_abb2conf_paratope_histo[cat_abb2].append(pdb)
        
        print(f"pdb {pdb} ave_plddt {ave_plddt}, cat {cat}, ave_abb2conf {ave_abb2conf}, cat_abb2 {cat_abb2}")
    
# save data
ave_plddt_paratope_df = pd.DataFrame(ave_plddt_paratope_list)
ave_plddt_paratope_df.columns = ["pdb", "para-ave-plddt"]
# subset of ref_data: the af2 subset
ref_data_af2 = ref_data_df.loc[ref_data_df["method"] == "AF2"]
# merging
af2_values_paratope_full = pd.merge(ref_data_af2, ave_plddt_paratope_df, on='pdb')
af2_values_paratope_full.to_csv(Path("data", "AF2_antibody_rmsd_plddt.csv"), sep="\t")

# abb2
ave_abb2conf_paratope_df = pd.DataFrame(ave_abb2conf_paratope_list)
ave_abb2conf_paratope_df.columns = ["pdb", "para-ave-conf"]
# subset of ref_data: the af2 subset
ref_data_abb2 = ref_data_df.loc[ref_data_df["method"] == "ABB2"]
# merging
abb2_values_paratope_full = pd.merge(ref_data_abb2, ave_abb2conf_paratope_df, on='pdb')
abb2_values_paratope_full.to_csv(Path("data", "ABB2_antibody_rmsd_conf.csv"), sep="\t")