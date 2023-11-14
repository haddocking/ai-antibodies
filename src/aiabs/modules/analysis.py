import os
import pandas as pd
import numpy as np
import logging
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from pathlib import Path

log = logging.getLogger("aiabsanalysis")

def read_capri_table(capri_filename, comment="#"):
    """
    Read capri table with pandas.
    
    Parameters
    ----------
    capri_filename : str or Path
        capri single structure filename
    comment : str
        the string used to denote a commented line in capri tables
    Returns
    -------
    capri_df : pandas DataFrame
        dataframe of capri values
    """
    capri_df = pd.read_csv(
        capri_filename,
        sep="\t",
        comment=comment)
    return capri_df

def read_caprievals(module_id, rel_path, capri_string="caprieval"):
    capri_dir = str(module_id) + f"_{capri_string}"
    df = None
    df_clt = None
    if capri_dir not in os.listdir(rel_path):
        log.warning("directory not found")
    else:
        capri_path = os.path.join(rel_path, capri_dir)
        if "capri_ss.tsv" in os.listdir(capri_path):
            df = read_capri_table(os.path.join(rel_path, capri_dir) + "/capri_ss.tsv")
        if "capri_clt.tsv" in os.listdir(capri_path):
            df_clt = read_capri_table(os.path.join(rel_path, capri_dir) + "/capri_clt.tsv")
    return df, df_clt

def acceptable_clusters_loose(caprieval_ss):
    """
    if one of the first four models of a cluster is correct, then the cluster is correct
    """
    if "model-cluster-ranking" in caprieval_ss.columns:
        keyword = "model-cluster-ranking"
    else:
        keyword = "self.model-cluster-ranking"
    npu = np.unique(caprieval_ss[keyword])
    if npu == "-":
        #print("not enough models in caprieval_ss")
        return None, None, None
    capri_df_clust = caprieval_ss.where((caprieval_ss[keyword] <= 3)).dropna()
    #print(f"capri_df_clust {capri_df_clust}")
    #cl_ids = np.unique(capri_df_clust["cluster-ranking"])
    #(f"cluster ids: {cl_ids}")
    capri_acc = []
    capri_medium = []
    capri_high = []
    for n in range(5): # top 5 clusters
        #capri_df_clust = capri_df.loc[capri_df["cluster-ranking"] <= n+1 & capri_df["model-cluster-ranking"] <= 4]
        capri_df_clust_subset = capri_df_clust.where((capri_df_clust["cluster-ranking"] <= n+1)).dropna()
        uni_cl = np.unique(capri_df_clust_subset["cluster-ranking"])
        acc_cl = 0
        med_cl = 0
        high_cl = 0
        for cl in uni_cl:
            cl_df = capri_df_clust_subset.where(capri_df_clust_subset["cluster-ranking"] == cl).dropna()
            cl_acc_len = len(acceptable_models_capri(cl_df))
            cl_med_len = len(medium_models_capri(cl_df))
            cl_high_len = len(high_models_capri(cl_df))
            if cl_acc_len > 0:
                acc_cl += 1
            if cl_med_len > 0:
                med_cl += 1
            if cl_high_len > 0:
                high_cl += 1
        capri_acc.append(acc_cl)
        capri_medium.append(med_cl)
        capri_high.append(high_cl)
    return capri_acc, capri_medium, capri_high

def acceptable_models_capri(caprieval_ss, threshold=100):
    red_df = caprieval_ss.iloc[:threshold,:]
    npw = np.where(red_df["fnat"] > 0.1)
    npw_lrmsd = np.where(red_df["lrmsd"].iloc[npw[0]] < 10.0)
    npw_irmsd = np.where(red_df["irmsd"].iloc[npw[0]] < 4.0)
    np_conc = np.concatenate((npw[0][npw_lrmsd[0]], npw[0][npw_irmsd[0]]))
    return np.unique(np_conc)+1
    

def medium_models_capri(caprieval_ss,threshold=100):
    red_df = caprieval_ss.iloc[:threshold,:]
    npw = np.where(red_df["fnat"] > 0.3)
    npw_lrmsd = np.where(red_df["lrmsd"].iloc[npw[0]] < 5.0)
    npw_irmsd = np.where(red_df["irmsd"].iloc[npw[0]] < 2.0)
    np_conc = np.concatenate((npw[0][npw_lrmsd[0]], npw[0][npw_irmsd[0]]))
    return np.unique(np_conc)+1


def high_models_capri(caprieval_ss,threshold=100):
    red_df = caprieval_ss.iloc[:threshold,:]
    npw = np.where(red_df["fnat"] > 0.5)
    npw_lrmsd = np.where(red_df["lrmsd"].iloc[npw[0]] < 1.0)
    npw_irmsd = np.where(red_df["irmsd"].iloc[npw[0]] < 1.0)
    np_conc = np.concatenate((npw[0][npw_lrmsd[0]], npw[0][npw_irmsd[0]]))
    return np.unique(np_conc)+1


def get_colors_top(list_of_values,acc,mediums,highs):
    colors = ["w","w","w","w","w","w"]
    for i in range(len(list_of_values)):
        values = range(list_of_values[i]+1)
        if any(val in acc for val in values):
            colors[i] = "aquamarine"
    # mediums
    for i in range(len(list_of_values)):
        values = range(list_of_values[i]+1)
        if any(val in mediums for val in values):
            colors[i] = "lightgreen"
    # highs
    for i in range(len(list_of_values)):
        values = range(list_of_values[i]+1)
        if any(val in highs for val in values):
            colors[i] = "darkgreen"
    return colors

def get_acc_models(capri_ss, clt_loose=False):
    """
    get single structure acceptable models
    """
    a_models, m_models, h_models = {}, {}, {}
    a_clusters_loose, m_clusters_loose, h_clusters_loose = {}, {}, {}
    
    for capri_step in capri_ss.keys():
        log.info(f"capri_step {capri_step}")
        if capri_step not in a_models.keys():
            a_models[capri_step], m_models[capri_step], h_models[capri_step] = {}, {}, {}
            if clt_loose:
                a_clusters_loose[capri_step], m_clusters_loose[capri_step], h_clusters_loose[capri_step] = {}, {}, {}
        for pdb in capri_ss[capri_step].keys():
            if pdb not in a_models[capri_step].keys():
                a_models[capri_step][pdb], m_models[capri_step][pdb], h_models[capri_step][pdb] = {}, {}, {}
                if clt_loose:
                    a_clusters_loose[capri_step][pdb], m_clusters_loose[capri_step][pdb], h_clusters_loose[capri_step][pdb] = {}, {}, {}
            for run in capri_ss[capri_step][pdb]:
                a_models[capri_step][pdb][run] = acceptable_models_capri(capri_ss[capri_step][pdb][run])
                m_models[capri_step][pdb][run] = medium_models_capri(capri_ss[capri_step][pdb][run])
                h_models[capri_step][pdb][run] = high_models_capri(capri_ss[capri_step][pdb][run])
                if clt_loose:
                    a_clusters_loose[capri_step][pdb][run], m_clusters_loose[capri_step][pdb][run], h_clusters_loose[capri_step][pdb][run] = acceptable_clusters_loose(capri_ss[capri_step][pdb][run])
                    #print(f"run {run} {a_clusters_loose[capri_step][pdb][run]} acceptable, {m_clusters_loose[capri_step][pdb][run]} medium, {h_clusters_loose[capri_step][pdb][run]} high")
    return a_models, m_models, h_models, a_clusters_loose, m_clusters_loose, h_clusters_loose

def retrieve_cell_text(capri_ss, key, capri_step, columns, run_dict, a_models, m_models, h_models):
    cell_text = {"acc": [], "med" : [], "high" : []}
    top_rank = [int(el.lstrip("top")) for el in columns[:-1]]
    for run in sorted(list(capri_ss[key][capri_step]["7bbj"].keys())):
        log.info(f"capri_step {capri_step} run {run}{os.linesep}")
        colors, no_acc = [], []
        for pdb in capri_ss[key][capri_step].keys():
            try:
                new_colors = get_colors_top(top_rank,
                                            a_models[key][capri_step][pdb][run],
                                            m_models[key][capri_step][pdb][run],
                                            h_models[key][capri_step][pdb][run])
                colors.append(new_colors)
                if new_colors[-2] == "w":
                    no_acc.append(pdb)
            except Exception as e:
                raise Exception(f"exception {e} occurred for {key} {capri_step} {pdb} {run}")
        if no_acc != []:
            if "Para-Epi" in run:
                log.warning(f"some pdbs have no acceptable models: {no_acc}")
        # colors
        df_colors = pd.DataFrame(colors)
        df_colors.columns = columns
        acc_row = []
        med_row = []
        hig_row = []
        for col in df_colors.columns[:-1]:
            base_dict = {"aquamarine" : 0, "darkgreen" : 0, "w" : 0, "lightgreen" : 0}
            npu = np.unique(df_colors[col], return_counts = True)
            for l in range(len(npu[0])):
                base_dict[npu[0][l]] = npu[1][l]
            nacc = base_dict["aquamarine"] + base_dict["darkgreen"] + base_dict["lightgreen"]
            nmed = base_dict["darkgreen"] + base_dict["lightgreen"]
            nhig = base_dict["darkgreen"]
            ntot = base_dict["aquamarine"] + base_dict["darkgreen"] + base_dict["lightgreen"] + base_dict["w"]
            acc_row.append(round(nacc/ntot, 2))
            med_row.append(round(nmed/ntot, 2))
            hig_row.append(round(nhig/ntot, 2))
        run_splt = run.split("-")
        run_name = f"{run_dict[run_splt[1]]}-{run_splt[2]}-{run_splt[3]}"
        if run_splt[4] == "AA":
            run_name += "-AA"
        if run_splt[1] == "ens":
            run_name += f"-{run_splt[-1]}"
        acc_row.append(run_name)
        med_row.append(run_name)
        hig_row.append(run_name)
        cell_text["acc"].append(acc_row)
        cell_text["med"].append(med_row)
        cell_text["high"].append(hig_row)
    return cell_text

def create_table(cell_text, columns, table_name, table_title):
    n_cols = len(columns)
    n_rows = len(cell_text)
    #log.info(f"n_cols {n_cols} n_rows {n_rows}")
    fig, ax = plt.subplots(figsize=(12,12))
    ax.axis('tight')
    ax.axis('off')
    colWidths = [0.1,0.1,0.1,0.1,0.1,0.3]
    the_table = ax.table(cellText=cell_text, colLabels=columns, loc='center', colWidths=colWidths)
    the_table.set_fontsize(40)
    # best rendering happens when scale is 1,4 for 6 columns and 9 rows
    scaling_factor = 0.9*n_rows/n_cols
    the_table.scale(1,scaling_factor)
    plt.title(table_title, fontsize=24)
    plt.savefig(table_name, dpi=200)
    plt.close()

def create_comparison_graphs(capri_dict, ref_data, output_file_basename):
    """
    compare each quantity to the reference data
    """
    data = {}
    df_dict = {
        "ab" : "ABB2",
        "abl" : "ABlooper",
        "ig" : "IgFold",
        "af2" : "AF2"
    }
    color_dict = {
        "ab" : plt.cm.tab10.colors[0],
        "abl" : plt.cm.tab10.colors[1],
        "ig" : plt.cm.tab10.colors[2],
        "af2" : plt.cm.tab10.colors[3],
        "abens" : plt.cm.tab10.colors[4],
        "ens" : "black",
        "ensnoaf2" : "gray",
        "af2ab" : plt.cm.tab10.colors[5],
        "af2abl" : plt.cm.tab10.colors[6],
        "af2ig" : plt.cm.tab10.colors[7],
        "af2af2" : plt.cm.tab10.colors[8],
        "af2abens" : plt.cm.tab10.colors[9],
        "af2ens" : "cyan",
        "af2ensnoaf2" : "brown",
        "bound" : "darkblue"
    }
    for pdb in capri_dict.keys():
        #pdb_line = np.where(ref_data["pdb"] == pdb)[0][0]
        #pdb_row = ref_df.iloc[pdb_line, :]
        for run in capri_dict[pdb].keys():
            run_family = "-".join(run.split("-")[2:5])
            if run_family not in data.keys():
                data[run_family] = []
            
            run_ab = run.split("-")[1]
            #print(f"run_ab {run_ab}")
            if run_ab in df_dict.keys():
                cdrh1 = ref_data[pdb][df_dict[run_ab]]["cdrh1"]
                cdrh2 = ref_data[pdb][df_dict[run_ab]]["cdrh2"]
                cdrh3 = ref_data[pdb][df_dict[run_ab]]["cdrh3"]
                cdrl1 = ref_data[pdb][df_dict[run_ab]]["cdrl1"]
                cdrl2 = ref_data[pdb][df_dict[run_ab]]["cdrl2"]
                cdrl3 = ref_data[pdb][df_dict[run_ab]]["cdrl3"]
            elif run_ab in ["ens", "af2ens"]:
                cdrh1 = np.mean([ref_data[pdb]["AF2"]["cdrh1"], ref_data[pdb]["ABB2"]["cdrh1"], ref_data[pdb]["IgFold"]["cdrh1"], ref_data[pdb]["ABlooper"]["cdrh1"]])
                cdrh2 = np.mean([ref_data[pdb]["AF2"]["cdrh2"], ref_data[pdb]["ABB2"]["cdrh2"], ref_data[pdb]["IgFold"]["cdrh2"], ref_data[pdb]["ABlooper"]["cdrh2"]])
                cdrh3 = np.mean([ref_data[pdb]["AF2"]["cdrh3"], ref_data[pdb]["ABB2"]["cdrh3"], ref_data[pdb]["IgFold"]["cdrh3"], ref_data[pdb]["ABlooper"]["cdrh3"]])
                cdrl1 = np.mean([ref_data[pdb]["AF2"]["cdrl1"], ref_data[pdb]["ABB2"]["cdrl1"], ref_data[pdb]["IgFold"]["cdrl1"], ref_data[pdb]["ABlooper"]["cdrl1"]])
                cdrl2 = np.mean([ref_data[pdb]["AF2"]["cdrl2"], ref_data[pdb]["ABB2"]["cdrl2"], ref_data[pdb]["IgFold"]["cdrl2"], ref_data[pdb]["ABlooper"]["cdrl2"]])
                cdrl3 = np.mean([ref_data[pdb]["AF2"]["cdrl3"], ref_data[pdb]["ABB2"]["cdrl3"], ref_data[pdb]["IgFold"]["cdrl3"], ref_data[pdb]["ABlooper"]["cdrl3"]])
            elif run_ab in ["ensnoaf2", "af2ensnoaf2"]:
                cdrh1 = np.mean([ref_data[pdb]["ABB2"]["cdrh1"], ref_data[pdb]["IgFold"]["cdrh1"], ref_data[pdb]["ABlooper"]["cdrh1"]])
                cdrh2 = np.mean([ref_data[pdb]["ABB2"]["cdrh2"], ref_data[pdb]["IgFold"]["cdrh2"], ref_data[pdb]["ABlooper"]["cdrh2"]])
                cdrh3 = np.mean([ref_data[pdb]["ABB2"]["cdrh3"], ref_data[pdb]["IgFold"]["cdrh3"], ref_data[pdb]["ABlooper"]["cdrh3"]])
                cdrl1 = np.mean([ref_data[pdb]["ABB2"]["cdrl1"], ref_data[pdb]["IgFold"]["cdrl1"], ref_data[pdb]["ABlooper"]["cdrl1"]])
                cdrl2 = np.mean([ref_data[pdb]["ABB2"]["cdrl2"], ref_data[pdb]["IgFold"]["cdrl2"], ref_data[pdb]["ABlooper"]["cdrl2"]])
                cdrl3 = np.mean([ref_data[pdb]["ABB2"]["cdrl3"], ref_data[pdb]["IgFold"]["cdrl3"], ref_data[pdb]["ABlooper"]["cdrl3"]])
            elif run_ab in ["abens", "af2abens"]:
                cdrh1 = np.mean([ref_data[pdb]["ABB2_ensemble_0"]["cdrh1"], ref_data[pdb]["ABB2_ensemble_1"]["cdrh1"], ref_data[pdb]["ABB2_ensemble_2"]["cdrh1"], ref_data[pdb]["ABB2_ensemble_3"]["cdrh1"]])
                cdrh2 = np.mean([ref_data[pdb]["ABB2_ensemble_0"]["cdrh2"], ref_data[pdb]["ABB2_ensemble_1"]["cdrh2"], ref_data[pdb]["ABB2_ensemble_2"]["cdrh2"], ref_data[pdb]["ABB2_ensemble_3"]["cdrh2"]])
                cdrh3 = np.mean([ref_data[pdb]["ABB2_ensemble_0"]["cdrh3"], ref_data[pdb]["ABB2_ensemble_1"]["cdrh3"], ref_data[pdb]["ABB2_ensemble_2"]["cdrh3"], ref_data[pdb]["ABB2_ensemble_3"]["cdrh3"]])
                cdrl1 = np.mean([ref_data[pdb]["ABB2_ensemble_0"]["cdrl1"], ref_data[pdb]["ABB2_ensemble_1"]["cdrl1"], ref_data[pdb]["ABB2_ensemble_2"]["cdrl1"], ref_data[pdb]["ABB2_ensemble_3"]["cdrl1"]])
                cdrl2 = np.mean([ref_data[pdb]["ABB2_ensemble_0"]["cdrl2"], ref_data[pdb]["ABB2_ensemble_1"]["cdrl2"], ref_data[pdb]["ABB2_ensemble_2"]["cdrl2"], ref_data[pdb]["ABB2_ensemble_3"]["cdrl2"]])
                cdrl3 = np.mean([ref_data[pdb]["ABB2_ensemble_0"]["cdrl3"], ref_data[pdb]["ABB2_ensemble_1"]["cdrl3"], ref_data[pdb]["ABB2_ensemble_2"]["cdrl3"], ref_data[pdb]["ABB2_ensemble_3"]["cdrl3"]])
            # computing ave quantities
            ave_dockq = np.mean(capri_dict[pdb][run]["dockq"])
            ave_fnat = np.mean(capri_dict[pdb][run]["fnat"])
            ave_irmsd = np.mean(capri_dict[pdb][run]["irmsd"])
            ave_lrmsd = np.mean(capri_dict[pdb][run]["lrmsd"])

            run_ab_color = color_dict[run_ab]
            data[run_family].append([pdb, run_ab, run_ab_color,
                                     cdrh1, cdrh2, cdrh3, cdrl1, cdrl2, cdrl3, 
                                     ave_dockq, ave_fnat, ave_irmsd, ave_lrmsd])
    for run_family in data.keys():
        if data[run_family] != []: #remove this check
            data[run_family] = pd.DataFrame(data[run_family])
            data[run_family].columns = ["pdb", "run", "run_ab_color",
                                    "cdrh1", "cdrh2", "cdrh3", "cdrl1", "cdrl2", "cdrl3",
                                    "ave_dockq", "ave_fnat", "ave_irmsd", "ave_lrmsd"]
    
    # plotting observables
    for obs in ["ave_dockq", "ave_fnat", "ave_irmsd", "ave_lrmsd"]:
        for comp_obs in ["cdrh1", "cdrh2", "cdrh3", "cdrl1", "cdrl2", "cdrl3"]:
            for run_family in data.keys():
                output_fname = f"{output_file_basename}-{run_family}-{obs}-{comp_obs}.png"
                fig, ax = plt.subplots()
                for i, dff in data[run_family].groupby("run_ab_color"):
                    try:
                        corr_coeff = pearsonr(dff[comp_obs], dff[obs])
                    except Exception as e:
                        log.info(f"corr_coeff {run_family}-{obs}-{comp_obs} failed: {e}")
                    #log.info(f"corr_coeff {run_family}-{obs}-{comp_obs} = {corr_coeff}")
                    scatter = ax.scatter(dff[comp_obs], dff[obs],
                                         color=dff["run_ab_color"].iloc[0], label = f"{dff['run'].iloc[0]} ({round(corr_coeff[0], 2)})")
                plt.title(output_fname.split("/")[-1])
                plt.ylabel(obs)
                plt.xlabel(comp_obs)
                plt.legend()
                if obs in ["ave_dockq", "ave_fnat"]:
                    plt.ylim((0.0,1.0))
                plt.savefig(output_fname)
                plt.close()

    return data

def generate_tables(capri_ss, capri_clt, ref_folder):
    """
    generate top models table.

    Parameters
    ----------
    capri_ss : dict
        dictionary with the capri ss data
    capri_clt : dict
        dictionary with the capri clt data
    ref_folder : str
        path to the reference folder
    """
    # get the acceptable, medium, and high-quality models.
    # these are dictionaries that have the same structure as the capri_ss and capri_clt ones
    a_models_capri, m_models_capri, h_models_capri, a_models_clt_loose, m_models_clt_loose, h_models_clt_loose = get_acc_models(capri_ss, clt_loose=True)
    a_clusters_capri, m_clusters_capri, h_clusters_capri, a, b, c = get_acc_models(capri_clt)
    #print(f"a_models_clt_loose: {a_models_clt_loose}")
    #print(f"a_models_capri: {a_models_capri}")
    # baseline caprieval tables
    run_dict = {
        "ig" : "IgFold",
        "ab" : "ABBuilder",
        "abl" : "ABLooper",
        "af2" : "AF2",
        "ens" : "ens",
        "ensnoaf2" : "ensnoaf2"
    }
    acc_nums = [1,5,10,20,48]
    acc_column_names = [f"acc_top{n}" for n in acc_nums]
    med_column_names = [f"med_top{n}" for n in acc_nums]
    high_column_names = [f"high_top{n}" for n in acc_nums]
    column_names = ["capri_step", "pdb", "run"]
    column_names.extend(acc_column_names)
    column_names.extend(med_column_names)
    column_names.extend(high_column_names)

    data = []
    for capri_step in list(capri_ss.keys()):
        for pdb in capri_ss[capri_step].keys():
            for run in capri_ss[capri_step][pdb].keys():
                #print(f"step {capri_step} pdb {pdb} run {run} a_models_capri for {pdb} {a_models_capri[capri_step][pdb][run]}")
                acc_top = {}
                med_top = {}
                high_top = {}
                for num in acc_nums:
                    top_acc = len([el for el in a_models_capri[capri_step][pdb][run] if el <= num])
                    top_med = len([el for el in m_models_capri[capri_step][pdb][run] if el <= num])
                    top_high = len([el for el in h_models_capri[capri_step][pdb][run] if el <= num])
                    acc_top[num] = top_acc
                    med_top[num] = top_med
                    high_top[num] = top_high
                data_list = [capri_step, pdb, run]
                data_list.extend(list(acc_top.values()))
                data_list.extend(list(med_top.values()))
                data_list.extend(list(high_top.values()))
                data.append(data_list)

    df_overall = pd.DataFrame(data)
    #print(f"df_overall: {df_overall}")
    df_overall.columns = column_names
    df_overall.to_csv(Path(ref_folder, "df_ss.tsv"), sep=" ", index=None)

    acc_nums = [1,2,3,4,5]
    acc_column_names = [f"acc_top{n}" for n in acc_nums]
    med_column_names = [f"med_top{n}" for n in acc_nums]
    high_column_names = [f"high_top{n}" for n in acc_nums]
    column_names = ["capri_step", "pdb", "run"]
    column_names.extend(acc_column_names)
    column_names.extend(med_column_names)
    column_names.extend(high_column_names)

    cl_data = []
    for capri_step in list(capri_clt.keys()):
        for pdb in capri_clt[capri_step].keys():
            for run in capri_clt[capri_step][pdb].keys():
                #print(f"step {capri_step} pdb {pdb} run {run} a_clusters_capri for {pdb} {a_clusters_capri[capri_step][pdb][run]}")
                acc_top = {}
                med_top = {}
                high_top = {}
                for num in acc_nums:
                    top_acc = len([el for el in a_clusters_capri[capri_step][pdb][run] if el <= num])
                    top_med = len([el for el in m_clusters_capri[capri_step][pdb][run] if el <= num])
                    top_high = len([el for el in h_clusters_capri[capri_step][pdb][run] if el <= num])
                    acc_top[num] = top_acc
                    med_top[num] = top_med
                    high_top[num] = top_high
                cl_data_list = [capri_step, pdb, run]
                cl_data_list.extend(list(acc_top.values()))
                cl_data_list.extend(list(med_top.values()))
                cl_data_list.extend(list(high_top.values()))
                cl_data.append(cl_data_list)

    df_clusters_overall = pd.DataFrame(cl_data)
    df_clusters_overall.columns = column_names
    df_clusters_overall.to_csv(Path(ref_folder, "df_clt.tsv"), sep=" ", index=None)
    
    # cluster loose definition
    acc_nums = [1,2,3,4,5]
    # column names are the same, no need to change them
    cl_loose_data = []
    for capri_step in list(capri_clt.keys()):
        for pdb in capri_clt[capri_step].keys():
            for run in capri_clt[capri_step][pdb].keys():
                #print(f"step {capri_step} pdb {pdb} run {run} a_clusters_capri for {pdb} {a_clusters_capri[capri_step][pdb][run]}")
                acc_top = {}
                med_top = {}
                high_top = {}
                for num in acc_nums:
                    top_acc = a_models_clt_loose[capri_step][pdb][run][num-1]
                    top_med = m_models_clt_loose[capri_step][pdb][run][num-1]
                    top_high = h_models_clt_loose[capri_step][pdb][run][num-1]
                    acc_top[num] = top_acc
                    med_top[num] = top_med
                    high_top[num] = top_high
                cl_loose_data_list = [capri_step, pdb, run]
                cl_loose_data_list.extend(list(acc_top.values()))
                cl_loose_data_list.extend(list(med_top.values()))
                cl_loose_data_list.extend(list(high_top.values()))
                cl_loose_data.append(cl_loose_data_list)

    df_clusters_loose_overall = pd.DataFrame(cl_loose_data)
    df_clusters_loose_overall.columns = column_names
    df_clusters_loose_overall.to_csv(Path(ref_folder, "df_clt_loose.tsv"), sep=" ", index=None)
