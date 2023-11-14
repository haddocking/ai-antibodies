import numpy as np
import pandas as pd
from pathlib import Path

pdbs_to_exclude = ["6xsw", "7cj2", "7e9b", "7m3n", "7n3c", "7n3d", "7o52", "7or9"]

def load_data():
    """
    Load all the necessary data for the analysis.
    """
    df_ss = pd.read_csv(Path("data", "df_ss.tsv"), sep=" ")
    df_ss_flexref = pd.read_csv(Path("data", "df_ss_flexref.tsv"), sep=" ")
    # bound data
    df_ss_bound = pd.read_csv(Path("data", "df_ss_bound.tsv"), sep=" ")
    df_ss_bound_flexref = pd.read_csv(Path("data", "df_ss_bound_flexref.tsv"), sep=" ")
    # rigid emref data
    df_ss_rigid_emref = pd.read_csv(Path("data", "df_ss_rigid-emref.tsv"), sep=" ")

    # ref data
    ref_data = pd.read_csv(Path("data", "ref_data.csv"), sep=" ")
    # zdock data
    zdock_ss = pd.read_csv(Path("data", "df_ss_zdock.tsv"), sep=" ")

    # remove pdbs
    df_ss = df_ss.loc[~df_ss["pdb"].isin(pdbs_to_exclude)]
    df_ss_flexref = df_ss_flexref.loc[~df_ss_flexref["pdb"].isin(pdbs_to_exclude)]
    df_ss_bound = df_ss_bound.loc[~df_ss_bound["pdb"].isin(pdbs_to_exclude)]
    df_ss_bound_flexref = df_ss_bound_flexref.loc[~df_ss_bound_flexref["pdb"].isin(pdbs_to_exclude)]
    df_ss_rigid_emref = df_ss_rigid_emref.loc[~df_ss_rigid_emref["pdb"].isin(pdbs_to_exclude)]
    ref_data = ref_data.loc[~ref_data["pdb"].isin(pdbs_to_exclude)]
    zdock_ss = zdock_ss.loc[~zdock_ss["pdb"].isin(pdbs_to_exclude)]

    # check that the number of runs is correct
    tot_runs = np.unique(df_ss["pdb"]).shape[0] # should be 79
    print(f"total number of runs {tot_runs}")
    assert tot_runs == 71
    zdock_tot_runs = np.unique(zdock_ss["pdb"]).shape[0] # should be 79
    print(f"total number of ZDOCK runs {zdock_tot_runs}")
    assert zdock_tot_runs == 71

    # rigidbody dfs
    rigidbody_capri = df_ss.loc[df_ss["capri_step"].isin([2,3])]
    rigidbody_capri_bound = df_ss_bound.loc[df_ss_bound["capri_step"].isin([2,3])]

    # emref dfs
    emref_capri7_list = ["run-ens-Para-Epi-mpi-196-clt", "run-ens-CDR-EpiVag-mpi-196-clt", 
                         "run-ens-CDR-EpiVag-AA-mpi-196-clt", "run-af2ens-Para-Epi-mpi-196-clt",
                         "run-af2ens-CDR-EpiVag-mpi-196-clt", "run-af2ens-CDR-EpiVag-AA-mpi-196-clt"]
    emref_capri6_list = ["run-ens-Para-Epi-mpi-196-48", "run-ens-CDR-EpiVag-mpi-196-48",
                         "run-ens-CDR-EpiVag-AA-mpi-196-48", "run-af2ens-Para-Epi-mpi-196-48",
                         "run-af2ens-CDR-EpiVag-mpi-196-48", "run-af2ens-CDR-EpiVag-AA-mpi-196-48"]

    emref_capri5_list = [el for el in np.unique(df_ss["run"]) if el not in 
                         emref_capri7_list and el not in emref_capri6_list]

    emref_capri_5 = df_ss.loc[df_ss["run"].isin(emref_capri5_list)]
    emref_capri_5 = emref_capri_5.loc[emref_capri_5["capri_step"].isin([5])]

    emref_capri_7 = df_ss.loc[df_ss["run"].isin(emref_capri7_list)]
    emref_capri_7 = emref_capri_7.loc[emref_capri_7["capri_step"].isin([7])]

    emref_capri_6 = df_ss.loc[df_ss["run"].isin(emref_capri6_list)]
    emref_capri_6 = emref_capri_6.loc[emref_capri_6["capri_step"].isin([6])]

    emref_capri = pd.concat([emref_capri_5, emref_capri_6, emref_capri_7])

    emref_capri_bound = df_ss_bound.loc[df_ss_bound["capri_step"].isin([5])]
    # emref after rigid
    emref_capri_rigid = df_ss_rigid_emref.loc[df_ss_rigid_emref["capri_step"].isin([4])]

    return rigidbody_capri, rigidbody_capri_bound, emref_capri, emref_capri_bound, df_ss_flexref, df_ss_bound_flexref, zdock_ss, emref_capri_rigid

def load_clt_data():
    df_clt = pd.read_csv(Path("data", "df_clt.tsv"), sep=" ")
    df_clt_bound = pd.read_csv(Path("data", "df_clt_bound.tsv"), sep=" ")
    df_clt_loose = pd.read_csv(Path("data", "df_clt_loose.tsv"), sep=" ")
    df_clt_bound_loose = pd.read_csv(Path("data", "df_clt_bound_loose.tsv"), sep=" ")

    # remove pdbs to exclude
    df_clt = df_clt.loc[~df_clt["pdb"].isin(pdbs_to_exclude)]
    df_clt_bound = df_clt_bound.loc[~df_clt_bound["pdb"].isin(pdbs_to_exclude)]
    df_clt_loose = df_clt_loose.loc[~df_clt_loose["pdb"].isin(pdbs_to_exclude)]
    df_clt_bound_loose = df_clt_bound_loose.loc[~df_clt_bound_loose["pdb"].isin(pdbs_to_exclude)]
    return df_clt, df_clt_bound, df_clt_loose, df_clt_bound_loose

def create_bound_gap_dictionary(capri_df, acc_key="acc"):
    """
    Creates a dictionary with the number of models in the top N for each run.

    Parameters
    ----------
    capri_df : pandas.DataFrame
        DataFrame with the capri values of a specific stage.
    
    Returns
    -------
    bound_gap : dict
        Dictionary with the number of models in the top N for each run.
    """
    print(f"acc key is {acc_key}")
    af2_runs = [el for el in np.unique(capri_df["run"]) if el.startswith("run-af2") and not el.startswith("run-af2-")]
    tops = [1, 5, 10, 20, 48]
    tops_dict = {el : 0 for el in tops}
    af2_bound_gap = {"Para-Epi" : {}, "CDR-EpiVag" : {}, "CDR-EpiVag-AA": {}}
    # loop over af2_runs
    for run in af2_runs:
        if "Para-Epi" in run:
            af2_bound_gap["Para-Epi"][run] = tops_dict.copy()
        elif "CDR-EpiVag-AA" in run:
            af2_bound_gap["CDR-EpiVag-AA"][run] = tops_dict.copy()
        else:
            af2_bound_gap["CDR-EpiVag"][run] = tops_dict.copy()
    # fill the dictionary
    for top in tops:
        for run in af2_runs:
            topN = capri_df.loc[capri_df["run"] == run][["pdb", f"{acc_key}_top{top}"]]
            #
            splt_run = run.split("-")
            splt_run[1] = splt_run[1][3:]
            bound_run = "-".join(splt_run)
            #print(f"bound_run for run {run} is {bound_run}")
            bound_topN = capri_df.loc[capri_df["run"] == bound_run][["pdb", f"{acc_key}_top{top}"]]
            bound_mod = bound_topN.loc[bound_topN[f"{acc_key}_top{top}"] > 0].shape[0]
            af2_mod = topN.loc[topN[f"{acc_key}_top{top}"] > 0].shape[0]

            if "Para-Epi" in run:
                af2_bound_gap["Para-Epi"][run][top] = [bound_mod, af2_mod]
            elif "CDR-EpiVag-AA" in run:
                af2_bound_gap["CDR-EpiVag-AA"][run][top] = [bound_mod, af2_mod]
            else:
                af2_bound_gap["CDR-EpiVag"][run][top] = [bound_mod, af2_mod]
    # return the filled dictionary
    return af2_bound_gap

def create_bound_bound_dictionary(capri_bound, acc_key="acc"):
    """
    
    """
    print(f"acc key is {acc_key}")
    tops = [1, 5, 10, 20, 48]
    tops_dict = {el : 0 for el in tops}
    bound_bound = {"Para-Epi" : tops_dict.copy(), "CDR-EpiVag" : tops_dict.copy(), "CDR-EpiVag-AA": tops_dict.copy()}
    bound_bound_runs = np.unique(capri_bound["run"])
    for top in tops:
        for run in bound_bound_runs:
            topN = capri_bound.loc[capri_bound["run"] == run][["pdb", f"{acc_key}_top{top}"]]
            bound_bound_mod = topN.loc[topN[f"{acc_key}_top{top}"] > 0].shape[0]
            if "Para-Epi" in run:
                bound_bound["Para-Epi"][top] = bound_bound_mod
            elif "CDR-EpiVag-AA" in run:
                bound_bound["CDR-EpiVag-AA"][top] = bound_bound_mod
            else:
                bound_bound["CDR-EpiVag"][top] = bound_bound_mod
    return bound_bound

def create_bound_bound_clt_dictionary(capri_clt_bound):
    """
    
    """
    tops = [1, 2, 3, 4, 5]
    tops_dict = {el : 0 for el in tops}
    bound_bound_clt = {"Para-Epi" : tops_dict.copy(), "CDR-EpiVag" : tops_dict.copy(), "CDR-EpiVag-AA": tops_dict.copy()}
    bound_bound_runs = np.unique(capri_clt_bound["run"])
    for top in tops:
        for run in bound_bound_runs:
            topN = capri_clt_bound.loc[capri_clt_bound["run"] == run][["pdb", f"acc_top{top}"]]
            bound_bound_mod = topN.loc[topN[f"acc_top{top}"] > 0].shape[0]
            if "Para-Epi" in run:
                bound_bound_clt["Para-Epi"][top] = bound_bound_mod
            elif "CDR-EpiVag-AA" in run:
                bound_bound_clt["CDR-EpiVag-AA"][top] = bound_bound_mod
            else:
                bound_bound_clt["CDR-EpiVag"][top] = bound_bound_mod
    return bound_bound_clt

def create_zdock_dictionary(zdock_ss):
    """
    Creates a ZDOCK dictionary with the number of models in the top N for each run.
    """
    zdock_bound_runs = ["run_capri_zdock_ab-CDR-EpiVag", "run_capri_zdock_abl-CDR-EpiVag",
                    "run_capri_zdock_af2-CDR-EpiVag", "run_capri_zdock_ig-CDR-EpiVag",
                    "run_capri_zdock_ab-Para-Epi", "run_capri_zdock_abl-Para-Epi",
                    "run_capri_zdock_af2-Para-Epi", "run_capri_zdock_ig-Para-Epi"
                   ]
    zdock_af2_runs = ["run_capri_zdock_af2ab-CDR-EpiVag", "run_capri_zdock_af2abl-CDR-EpiVag",
                  "run_capri_zdock_af2af2-CDR-EpiVag", "run_capri_zdock_af2ig-CDR-EpiVag",
                  "run_capri_zdock_af2ab-Para-Epi", "run_capri_zdock_af2abl-Para-Epi",
                  "run_capri_zdock_af2af2-Para-Epi", "run_capri_zdock_af2ig-Para-Epi"]

    tops = [1, 5, 10, 20, 48]
    tops_dict = {el : 0 for el in tops}
    zdock_data = {"Para-Epi" : {}, "CDR-EpiVag" : {}}
    for el in ["Para-Epi" , "CDR-EpiVag"]:
        for run in zdock_bound_runs:
            run_name = run.split("_")[-1]
            if "Para-Epi" in run:
                zdock_data["Para-Epi"][run_name] = tops_dict.copy()
            elif "CDR-EpiVag" in run:
                zdock_data["CDR-EpiVag"][run_name] = tops_dict.copy()

    for top in tops:
        for n in range(len(zdock_bound_runs)):
            zdock_bound_run = zdock_bound_runs[n]
            zdock_af2_run = zdock_af2_runs[n]
            topN = zdock_ss.loc[zdock_ss["run"] == zdock_bound_run][["pdb", f"acc_top{top}"]]
            topN_af2 = zdock_ss.loc[zdock_ss["run"] == zdock_af2_run][["pdb", f"acc_top{top}"]]
            # filling the dict
            bound_mod = topN.loc[topN[f"acc_top{top}"] > 0].shape[0]
            af2_mod = topN_af2.loc[topN_af2[f"acc_top{top}"] > 0].shape[0]
            run_name = zdock_bound_run.split("_")[-1]
            if "Para-Epi" in run_name:
                zdock_data["Para-Epi"][run_name][top] = [bound_mod, af2_mod]
            else:
                zdock_data["CDR-EpiVag"][run_name][top] = [bound_mod, af2_mod]
    return zdock_data

def create_bound_gap_clt_dictionary(capri_clt):
    """
    """
    af2_runs = [el for el in np.unique(capri_clt["run"]) if el.startswith("run-af2") and not el.startswith("run-af2-")]
    tops = [1, 2, 3]
    tops_dict = {el : 0 for el in tops}
    af2_bound_gap_clt = {"Para-Epi" : {}, "CDR-EpiVag" : {}, "CDR-EpiVag-AA": {}}
    for run in af2_runs:
        if "Para-Epi" in run:
            af2_bound_gap_clt["Para-Epi"][run] = tops_dict.copy()
        elif "CDR-EpiVag-AA" in run:
            af2_bound_gap_clt["CDR-EpiVag-AA"][run] = tops_dict.copy()
        else:
            af2_bound_gap_clt["CDR-EpiVag"][run] = tops_dict.copy()

    for top in tops:
        for run in af2_runs:
            topN = capri_clt.loc[capri_clt["run"] == run][["pdb", f"acc_top{top}"]]
            #
            splt_run = run.split("-")
            splt_run[1] = splt_run[1][3:]
            bound_run = "-".join(splt_run)
            bound_topN = capri_clt.loc[capri_clt["run"] == bound_run][["pdb", f"acc_top{top}"]]
            bound_mod = bound_topN.loc[bound_topN[f"acc_top{top}"] > 0].shape[0]
            af2_mod = topN.loc[topN[f"acc_top{top}"] > 0].shape[0]

            if "Para-Epi" in run:
                af2_bound_gap_clt["Para-Epi"][run][top] = [bound_mod, af2_mod]
            elif "CDR-EpiVag-AA" in run:
                af2_bound_gap_clt["CDR-EpiVag-AA"][run][top] = [bound_mod, af2_mod]
            else:
                af2_bound_gap_clt["CDR-EpiVag"][run][top] = [bound_mod, af2_mod]
    return af2_bound_gap_clt

def get_epitope_and_paratope(constr_file):
    epitope_residues = {}
    paratope_residues = {}
    with open(constr_file, "r") as rfile:
        for ln in rfile:
            if not ln.startswith("#"):
                splt_ln_ab = ln.split(":")[0].split(",")
                splt_ln_ant = ln.split(":")[1].split(",")
                chain_ab = splt_ln_ab[0]
                resid_ab = splt_ln_ab[1]
                chain_ant = splt_ln_ant[0]
                resid_ant = int(splt_ln_ant[1])
                if chain_ant not in epitope_residues:
                    epitope_residues[chain_ant] = []
                if chain_ab not in paratope_residues:
                    paratope_residues[chain_ab] = []
                epitope_residues[chain_ant].append(resid_ant)
                paratope_residues[chain_ab].append(resid_ab)
    for ch in epitope_residues.keys():
        epitope_residues[ch] = list(np.unique(epitope_residues[ch]))
    for ch in paratope_residues.keys():
        paratope_residues[ch] = list(np.unique(paratope_residues[ch]))
    return epitope_residues, paratope_residues