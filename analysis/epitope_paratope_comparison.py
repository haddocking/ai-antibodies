import pandas as pd
from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np
from functions import load_data
plt.rcParams["font.family"] = "Helvetica"

cat_dict = {"Para-Epi": "Para-Epi",
            "CDR-EpiVag": "CDR-VagueEpi" ,
            "CDR-EpiVag-AA" : "CDR-VagueEpi-AA"}

# LOAD DATA
af2_values_full = pd.read_csv(Path("data", "AF2_antigen_rmsd_plddt_multi_regions.csv"))
af2_values = pd.read_csv(Path("data", "AF2_antigen_rmsd_plddt.csv"), sep="\t")
print(f"af2_values {af2_values}")
filter_pdbs = ['7kpj','7kn4', "6xsw", "7cj2", "7e9b", "7m3n", "7n3c", "7n3d", "7o52", "7or9"]

af2_values = af2_values.query('pdb not in @filter_pdbs').reset_index()
af2_values_full = af2_values_full.query('pdb not in @filter_pdbs').reset_index()
rigidbody_capri, rigidbody_capri_bound, emref_capri, emref_capri_bound, df_ss_flexref, df_ss_bound_flexref, zdock_ss, emref_rigid_capri = load_data()


# PLOT CORRELATION ON TWO PLOT
fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
plt.subplots_adjust(wspace=0.2, hspace=0.1)
plt.suptitle("pLDDT-RMSD correlation", fontsize=24, verticalalignment="center")

axs[0].scatter(af2_values_full["rmsd_para_epi"], af2_values_full["plddt_ave_para_epi"])
axs[0].set_xlabel("epitope RMSD ($\\AA$)", fontsize=24)
axs[0].set_ylabel("pLDDT", fontsize=24)
axs[0].set_ylim((50,101))
axs[0].set_xlim((-0.5,20.0))
axs[1].scatter(af2_values_full["rmsd_full"], af2_values_full["plddt_ave_full"])
axs[1].set_xlabel("full antigen RMSD ($\\AA$)", fontsize=24)
axs[1].set_ylim((50,101))
axs[1].set_xlim((-0.5,20.0))
# ticks
axs[1].set_yticks([])
axs[0].tick_params(axis='both', which='major', labelsize=16)
axs[1].tick_params(axis='both', which='major', labelsize=16)

for pdb in list(af2_values_full["pdb"]):
    df_subset = af2_values_full.loc[af2_values_full["pdb"] == pdb]
    epi_rmsd = df_subset["rmsd_para_epi"].iloc[0]
    full_rmsd = df_subset["rmsd_full"].iloc[0]
    epi_plddt = df_subset["plddt_ave_para_epi"].iloc[0]
    full_plddt = df_subset["plddt_ave_full"].iloc[0]
    if epi_rmsd > 5.0 or epi_plddt < 70.0:
        axs[0].text(x=epi_rmsd, y=epi_plddt+0.3, s=pdb.upper(), fontsize=14)
    if full_rmsd > 5.0 or full_plddt < 75.0:
        if pdb == "7mrz":
            xloc = full_rmsd - 2.0
            yloc = full_plddt + 1.0
        else:
            xloc = full_rmsd
            yloc = full_plddt + 0.3
        axs[1].text(x=xloc, y=yloc, s=pdb.upper(), fontsize=14)
        
        
plt.tight_layout()
plt.savefig(Path("figures", "epitope_paratope_comparison", f"plddt-rmsd.png"), dpi=300, bbox_inches="tight")

# PLOT CORRELATION ON ONE PLOT
fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(6, 6))
plt.subplots_adjust(wspace=0.2, hspace=0.1)
plt.suptitle("pLDDT-RMSD correlation (AF2 antigen)", fontsize=24, verticalalignment="center")

axs.scatter(af2_values_full["rmsd_para_epi"], af2_values_full["plddt_ave_para_epi"])
axs.set_xlabel("epitope RMSD ($\\AA$)", fontsize=24)
axs.set_ylabel("epitope pLDDT", fontsize=24)
axs.set_ylim((50,101))
axs.set_xlim((-0.2,6.5))

# ticks
axs.tick_params(axis='both', which='major', labelsize=16)

for pdb in list(af2_values_full["pdb"]):
    df_subset = af2_values_full.loc[af2_values_full["pdb"] == pdb]
    epi_rmsd = df_subset["rmsd_para_epi"].iloc[0]
    full_rmsd = df_subset["rmsd_full"].iloc[0]
    epi_plddt = df_subset["plddt_ave_para_epi"].iloc[0]
    full_plddt = df_subset["plddt_ave_full"].iloc[0]
    if epi_rmsd > 5.0 or epi_plddt < 70.0:
        axs.text(x=epi_rmsd-0.2, y=epi_plddt+0.75, s=pdb.upper(), fontsize=14)    

print(f"correlation coefficient for af2 epi-ave-plddt = {np.corrcoef(af2_values_full['rmsd_para_epi'], af2_values_full['plddt_ave_para_epi'])[0,1]}")
plt.tight_layout()
plt.savefig(Path("figures", "epitope_paratope_comparison", f"plddt-rmsd-epionly.png"), dpi=300, bbox_inches="tight")


# EPITOPE RMSD PLOT
# generating ave_af2rmsdepi_histo dictionary
ave_af2rmsdepi_histo = {}
ave_af2rmsdepi_histo["high"] = af2_values.where(af2_values["epi-cat"] == "high").dropna()["pdb"]
ave_af2rmsdepi_histo["good"] = af2_values.where(af2_values["epi-cat"] == "good").dropna()["pdb"]
ave_af2rmsdepi_histo["bad"] = af2_values.where(af2_values["epi-cat"] == "bad").dropna()["pdb"]
print(f"ave_af2rmsdepi_histo = {ave_af2rmsdepi_histo}")
tops = [1, 5, 10, 20, 48]
tops_dict = {el : 0 for el in tops}

af2_rmsdepi = {"Para-Epi" : {}, "CDR-EpiVag" : {}, "CDR-EpiVag-AA": {}}

cat_runs = {"Para-Epi" : {}, "CDR-EpiVag" : {}, "CDR-EpiVag-AA": {}}
for cat in cat_runs.keys():
    cat_runs[cat] = [f"run-af2ab-{cat}-mpi-50-50", f"run-af2abens-{cat}-mpi-50-50", f"run-af2abl-{cat}-mpi-50-50",
            f"run-af2af2-{cat}-mpi-50-50", f"run-af2ens-{cat}-mpi-50-50", f"run-af2ens-{cat}-mpi-196-48",
            f"run-af2ens-{cat}-mpi-196-clt", f"run-af2ensnoaf2-{cat}-mpi-50-50", f"run-af2ig-{cat}-mpi-50-50"]

for cat in ["Para-Epi" , "CDR-EpiVag", "CDR-EpiVag-AA"]:
    af2_rmsdepi[cat] = {}
    for stage in ["rigidbody", "emref", "flexref"]:
        af2_rmsdepi[cat][stage] = {}
        for quality in ["bad", "good", "high"]:
            af2_rmsdepi[cat][stage][quality] = {}
            for run in cat_runs[cat]:
                af2_rmsdepi[cat][stage][quality][run] = tops_dict.copy()

for cat in cat_runs.keys():
    for top in tops:
        for run in cat_runs[cat]:
            # rigidbody
            run_df_rigidbody = rigidbody_capri.loc[rigidbody_capri["run"] == run]
            
            bad_pdbs_epirmsd = run_df_rigidbody.loc[run_df_rigidbody["pdb"].isin(ave_af2rmsdepi_histo["bad"])]
            good_pdbs_epirmsd = run_df_rigidbody.loc[run_df_rigidbody["pdb"].isin(ave_af2rmsdepi_histo["good"])]
            high_pdbs_epirmsd = run_df_rigidbody.loc[run_df_rigidbody["pdb"].isin(ave_af2rmsdepi_histo["high"])]
            
            af2_mod_bad_rigid_epirmsd = bad_pdbs_epirmsd.loc[bad_pdbs_epirmsd[f"acc_top{top}"] > 0].shape[0]/bad_pdbs_epirmsd.shape[0]
            af2_mod_good_rigid_epirmsd = good_pdbs_epirmsd.loc[good_pdbs_epirmsd[f"acc_top{top}"] > 0].shape[0]/good_pdbs_epirmsd.shape[0]
            af2_mod_high_rigid_epirmsd = high_pdbs_epirmsd.loc[high_pdbs_epirmsd[f"acc_top{top}"] > 0].shape[0]/high_pdbs_epirmsd.shape[0]
            
            # emref
            run_df_emref = emref_capri.loc[emref_capri["run"] == run]
            
            bad_pdbs_epirmsd = run_df_emref.loc[run_df_emref["pdb"].isin(ave_af2rmsdepi_histo["bad"])]
            good_pdbs_epirmsd = run_df_emref.loc[run_df_emref["pdb"].isin(ave_af2rmsdepi_histo["good"])]
            high_pdbs_epirmsd = run_df_emref.loc[run_df_emref["pdb"].isin(ave_af2rmsdepi_histo["high"])]
            
            af2_mod_bad_emref_epirmsd = bad_pdbs_epirmsd.loc[bad_pdbs_epirmsd[f"acc_top{top}"] > 0].shape[0]/bad_pdbs_epirmsd.shape[0]
            af2_mod_good_emref_epirmsd = good_pdbs_epirmsd.loc[good_pdbs_epirmsd[f"acc_top{top}"] > 0].shape[0]/good_pdbs_epirmsd.shape[0]
            af2_mod_high_emref_epirmsd = high_pdbs_epirmsd.loc[high_pdbs_epirmsd[f"acc_top{top}"] > 0].shape[0]/high_pdbs_epirmsd.shape[0]
            # flexref
            run_df_flexref = df_ss_flexref.loc[df_ss_flexref["run"] == run]
            
            bad_pdbs_epirmsd = run_df_flexref.loc[run_df_flexref["pdb"].isin(ave_af2rmsdepi_histo["bad"])]
            good_pdbs_epirmsd = run_df_flexref.loc[run_df_flexref["pdb"].isin(ave_af2rmsdepi_histo["good"])]
            high_pdbs_epirmsd = run_df_flexref.loc[run_df_flexref["pdb"].isin(ave_af2rmsdepi_histo["high"])]
            
            af2_mod_bad_flexref_epirmsd = bad_pdbs_epirmsd.loc[bad_pdbs_epirmsd[f"acc_top{top}"] > 0].shape[0]/bad_pdbs_epirmsd.shape[0]
            af2_mod_good_flexref_epirmsd = good_pdbs_epirmsd.loc[good_pdbs_epirmsd[f"acc_top{top}"] > 0].shape[0]/good_pdbs_epirmsd.shape[0]
            af2_mod_high_flexref_epirmsd = high_pdbs_epirmsd.loc[high_pdbs_epirmsd[f"acc_top{top}"] > 0].shape[0]/high_pdbs_epirmsd.shape[0]
            
            af2_rmsdepi[cat]["rigidbody"]["bad"][run][top] = af2_mod_bad_rigid_epirmsd
            af2_rmsdepi[cat]["rigidbody"]["good"][run][top] = af2_mod_good_rigid_epirmsd
            af2_rmsdepi[cat]["rigidbody"]["high"][run][top] = af2_mod_high_rigid_epirmsd
            af2_rmsdepi[cat]["emref"]["bad"][run][top] = af2_mod_bad_emref_epirmsd
            af2_rmsdepi[cat]["emref"]["good"][run][top] = af2_mod_good_emref_epirmsd
            af2_rmsdepi[cat]["emref"]["high"][run][top] = af2_mod_high_emref_epirmsd
            af2_rmsdepi[cat]["flexref"]["bad"][run][top] = af2_mod_bad_flexref_epirmsd
            af2_rmsdepi[cat]["flexref"]["good"][run][top] = af2_mod_good_flexref_epirmsd
            af2_rmsdepi[cat]["flexref"]["high"][run][top] = af2_mod_high_flexref_epirmsd

# plot
plt.rcParams["font.family"] = "Helvetica"
fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(15, 12))
plt.subplots_adjust(wspace=0.1, hspace=0.1, bottom=0.25, top=0.925)
plt.suptitle("Epitope RMSD and SR", fontsize=24, verticalalignment="center")
col = 0
colors = plt.cm.Set2.colors

# rigidbody
for cat in ["Para-Epi", "CDR-EpiVag-AA", "CDR-EpiVag"]:
    protocols = ["ab" ,"abens", "abl", "af2", "ens196_48", "ens196_clt", "ens48_48", "ensnoaf", "ig"]
    runs = [f"run-af2ab-{cat}-mpi-50-50", f"run-af2abens-{cat}-mpi-50-50", f"run-af2abl-{cat}-mpi-50-50",
            f"run-af2af2-{cat}-mpi-50-50", f"run-af2ig-{cat}-mpi-50-50", f"run-af2ens-{cat}-mpi-196-48", 
            f"run-af2ens-{cat}-mpi-196-clt", f"run-af2ens-{cat}-mpi-50-50", f"run-af2ensnoaf2-{cat}-mpi-50-50"]
    mod_dict_1 = {
        "bad" : [af2_rmsdepi[cat]["rigidbody"]["bad"][run][1] for run in runs],
        "good": [af2_rmsdepi[cat]["rigidbody"]["good"][run][1] for run in runs],
        "high": [af2_rmsdepi[cat]["rigidbody"]["high"][run][1] for run in runs],
    }
    
    mod_dict_10 = {
        "bad" : [af2_rmsdepi[cat]["rigidbody"]["bad"][run][10] for run in runs],
        "good": [af2_rmsdepi[cat]["rigidbody"]["good"][run][10] for run in runs],
        "high": [af2_rmsdepi[cat]["rigidbody"]["high"][run][10] for run in runs],
    }
        
    n=len(mod_dict_1["bad"])
    r = np.arange(n)
    width = 0.3
    
    axs[0][col].bar(r, mod_dict_1["good"], color = colors[1],
                width = width, edgecolor = 'black',
                label='GOOD')
    axs[0][col].bar(r + width, mod_dict_1["high"], color = colors[2],
                width = width, edgecolor = 'black',
                label='HIGH')
    
    # TOP 10 blurred
    axs[0][col].bar(r, mod_dict_10["good"], color = colors[1], alpha=0.4,
                width = width, #edgecolor = 'black',
                label='GOOD')
    axs[0][col].bar(r + width, mod_dict_10["high"], color = colors[2], alpha=0.4,
                width = width, #edgecolor = 'black',
                label='HIGH')
    axs[0][col].set_ylim((0,1.01))
    if col > 0:
        axs[0][col].set_yticks([])
    else:
        axs[0][col].set_ylabel("Rigid-body Success Rate", fontsize=16)  
    axs[0][col].set_xticks([])
    col += 1

col = 0
for cat in ["Para-Epi", "CDR-EpiVag-AA", "CDR-EpiVag"]:
    runs = [f"run-af2ab-{cat}-mpi-50-50", f"run-af2abens-{cat}-mpi-50-50", f"run-af2abl-{cat}-mpi-50-50",
            f"run-af2af2-{cat}-mpi-50-50", f"run-af2ig-{cat}-mpi-50-50", f"run-af2ens-{cat}-mpi-196-48", 
            f"run-af2ens-{cat}-mpi-196-clt", f"run-af2ens-{cat}-mpi-50-50", f"run-af2ensnoaf2-{cat}-mpi-50-50", ]
    xticks_lower = []
    for run in runs:
        run_name = run.split('-')[1][3:]
        if run_name == "ab":
            run_name = "abb"
        if run_name == "abens":
            run_name = "abbe"
        if run_name == "ens":
            run_name += f"{run.split('-')[-2]}-{run.split('-')[-1]}"
        if run_name == "ens50-50":
            run_name = "ens"
        xticks_lower.append(run_name)
    xticks = [xt.upper() for xt in xticks_lower]
    #print(f"xticks {xticks}")
    mod_dict_1_emref = {
        "bad" : [af2_rmsdepi[cat]["emref"]["bad"][run][1] for run in runs],
        "good": [af2_rmsdepi[cat]["emref"]["good"][run][1] for run in runs],
        "high": [af2_rmsdepi[cat]["emref"]["high"][run][1] for run in runs],
    }
    #
    
    mod_dict_10_emref = {
        "bad" : [af2_rmsdepi[cat]["emref"]["bad"][run][10] for run in runs],
        "good": [af2_rmsdepi[cat]["emref"]["good"][run][10] for run in runs],
        "high": [af2_rmsdepi[cat]["emref"]["high"][run][10] for run in runs],
    }
    print(mod_dict_10_emref["bad"])
    axs[1][col].bar(r, mod_dict_1_emref["good"], color = colors[1],
                width = width, edgecolor = 'black',
                label='GOOD T1')
    axs[1][col].bar(r + width, mod_dict_1_emref["high"], color = colors[2],
                width = width, edgecolor = 'black',
                label='HIGH T1')
    
    # TOP 10 blurred
    axs[1][col].bar(r, mod_dict_10_emref["good"], color = colors[1], alpha=0.4,
                width = width, #edgecolor = 'black',
                label='GOOD T10')
    axs[1][col].bar(r + width, mod_dict_10_emref["high"], color = colors[2], alpha=0.4,
                width = width, #edgecolor = 'black',
                label='HIGH T10')
    axs[1][col].set_ylim((0,1.01))
    if col > 0:
        axs[1][col].set_yticks([])
    else:
        axs[1][col].set_ylabel("Emref Success Rate", fontsize=16)
    axs[1][col].set_xticks(r + width/2, xticks, fontsize=16, rotation = 90)
    col+=1

axs[0][0].set_title(cat_dict["Para-Epi"], fontsize=16)
axs[0][1].set_title(cat_dict["CDR-EpiVag-AA"], fontsize=16)
axs[0][2].set_title(cat_dict["CDR-EpiVag"], fontsize=16)

axs[1][1].set_xlabel("Docking protocol", fontsize=24, labelpad=20)

y_values = [0.2, 0.4, 0.6, 0.8]
for y in y_values:
    axs[0][0].axhline(y=y, color='black', linestyle='-', alpha=0.1, zorder=0)
    axs[0][1].axhline(y=y, color='black', linestyle='-', alpha=0.1, zorder=0)
    axs[0][2].axhline(y=y, color='black', linestyle='-', alpha=0.1, zorder=0)
    axs[1][0].axhline(y=y, color='black', linestyle='-', alpha=0.1, zorder=0)
    axs[1][1].axhline(y=y, color='black', linestyle='-', alpha=0.1, zorder=0)
    axs[1][2].axhline(y=y, color='black', linestyle='-', alpha=0.1, zorder=0)

lines_labels = [axs[1][1].get_legend_handles_labels()]
lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]

fig.legend(lines, labels, ncol=6, fontsize=16, loc="lower center", bbox_to_anchor = (0, 0.02, 1, 1))

plt.savefig(Path("figures", "epitope_paratope_comparison", "rigid_emref_1-10_rmsdepi.png"), dpi=300, bbox_inches="tight")

# same plot but only two columns
fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(10, 12))
plt.subplots_adjust(wspace=0.1, hspace=0.1, bottom=0.25, top=0.925)
plt.suptitle("Epitope RMSD and SR", fontsize=24, verticalalignment="center")
col = 0
colors = plt.cm.Set2.colors

# rigidbody
for cat in ["Para-Epi", "CDR-EpiVag-AA"]:
    print(f"\n{cat}\n")
    protocols = ["ab" ,"abens", "abl", "af2", "ens196_48", "ens196_clt", "ens48_48", "ensnoaf", "ig"]
    runs = [f"run-af2ab-{cat}-mpi-50-50", f"run-af2abens-{cat}-mpi-50-50", f"run-af2abl-{cat}-mpi-50-50",
            f"run-af2af2-{cat}-mpi-50-50", f"run-af2ig-{cat}-mpi-50-50", f"run-af2ens-{cat}-mpi-196-48",
            f"run-af2ens-{cat}-mpi-196-clt", f"run-af2ens-{cat}-mpi-50-50", f"run-af2ensnoaf2-{cat}-mpi-50-50"]
    mod_dict_1 = {
        "bad" : [af2_rmsdepi[cat]["rigidbody"]["bad"][run][1] for run in runs],
        "good": [af2_rmsdepi[cat]["rigidbody"]["good"][run][1] for run in runs],
        "high": [af2_rmsdepi[cat]["rigidbody"]["high"][run][1] for run in runs],
    }
    
    mod_dict_10 = {
        "bad" : [af2_rmsdepi[cat]["rigidbody"]["bad"][run][10] for run in runs],
        "good": [af2_rmsdepi[cat]["rigidbody"]["good"][run][10] for run in runs],
        "high": [af2_rmsdepi[cat]["rigidbody"]["high"][run][10] for run in runs],
    }
        
    n=len(mod_dict_1["bad"])
    r = np.arange(n)
    width = 0.3
    
    print(f"rigidbody mod_dict_1 {mod_dict_1}")

    axs[0][col].bar(r, mod_dict_1["good"], color = colors[1],
                width = width, edgecolor = 'black',
                label='GOOD')
    axs[0][col].bar(r + width, mod_dict_1["high"], color = colors[2],
                width = width, edgecolor = 'black',
                label='HIGH')
    
    # TOP 10 blurred
    axs[0][col].bar(r, mod_dict_10["good"], color = colors[1], alpha=0.4,
                width = width, #edgecolor = 'black',
                label='GOOD')
    axs[0][col].bar(r + width, mod_dict_10["high"], color = colors[2], alpha=0.4,
                width = width, #edgecolor = 'black',
                label='HIGH')
    axs[0][col].set_ylim((0,1.01))
    if col > 0:
        axs[0][col].set_yticks([])
    else:
        axs[0][col].set_ylabel("Rigid-body Success Rate", fontsize=16)  
    axs[0][col].set_xticks([])
    col += 1

col = 0
for cat in ["Para-Epi", "CDR-EpiVag-AA"]:
    print(f"\n{cat}\n")
    runs = [f"run-af2ab-{cat}-mpi-50-50", f"run-af2abens-{cat}-mpi-50-50", f"run-af2abl-{cat}-mpi-50-50",
            f"run-af2af2-{cat}-mpi-50-50", f"run-af2ig-{cat}-mpi-50-50", f"run-af2ens-{cat}-mpi-196-48", 
            f"run-af2ens-{cat}-mpi-196-clt", f"run-af2ens-{cat}-mpi-50-50", f"run-af2ensnoaf2-{cat}-mpi-50-50", ]
    xticks_lower = []
    for run in runs:
        run_name = run.split('-')[1][3:]
        if run_name == "ab":
            run_name = "abb"
        if run_name == "abens":
            run_name = "abbe"
        if run_name == "ens":
            run_name += f"{run.split('-')[-2]}-{run.split('-')[-1]}"
        if run_name == "ens50-50":
            run_name = "ens"
        xticks_lower.append(run_name)
    xticks = [xt.upper() for xt in xticks_lower]
    #print(f"xticks {xticks}")
    mod_dict_1_emref = {
        "bad" : [af2_rmsdepi[cat]["emref"]["bad"][run][1] for run in runs],
        "good": [af2_rmsdepi[cat]["emref"]["good"][run][1] for run in runs],
        "high": [af2_rmsdepi[cat]["emref"]["high"][run][1] for run in runs],
    }
    #
    
    mod_dict_10_emref = {
        "bad" : [af2_rmsdepi[cat]["emref"]["bad"][run][10] for run in runs],
        "good": [af2_rmsdepi[cat]["emref"]["good"][run][10] for run in runs],
        "high": [af2_rmsdepi[cat]["emref"]["high"][run][10] for run in runs],
    }
    print(f"mod_dict_1_emref[high] {mod_dict_1_emref['high']}")
    axs[1][col].bar(r, mod_dict_1_emref["good"], color = colors[1],
                width = width, edgecolor = 'black',
                label='GOOD T1')
    axs[1][col].bar(r + width, mod_dict_1_emref["high"], color = colors[2],
                width = width, edgecolor = 'black',
                label='HIGH T1')
    
    # TOP 10 blurred
    print(f"plotting {mod_dict_10_emref['bad']} for cat {cat}")
    axs[1][col].bar(r, mod_dict_10_emref["good"], color = colors[1], alpha=0.4,
                width = width, #edgecolor = 'black',
                label='GOOD T10')
    axs[1][col].bar(r + width, mod_dict_10_emref["high"], color = colors[2], alpha=0.4,
                width = width, #edgecolor = 'black',
                label='HIGH T10')
    axs[1][col].set_ylim((0,1.01))
    if col > 0:
        axs[1][col].set_yticks([])
    else:
        axs[1][col].set_ylabel("Emref Success Rate", fontsize=16)
    axs[1][col].set_xticks(r + width/2, xticks, fontsize=16, rotation = 90)
    col+=1

axs[0][0].set_title(cat_dict["Para-Epi"], fontsize=16)
axs[0][1].set_title(cat_dict["CDR-EpiVag-AA"], fontsize=16)
#axs[0][2].set_title(cat_dict["CDR-EpiVag"], fontsize=16)

# adding letters
axs[0][0].text(0.02, 0.93, "a)", transform=axs[0][0].transAxes, size=15, weight='bold')
axs[0][1].text(0.02, 0.93, "b)", transform=axs[0][1].transAxes, size=15, weight='bold')
axs[1][0].text(0.02, 0.93, "c)", transform=axs[1][0].transAxes, size=15, weight='bold')
axs[1][1].text(0.02, 0.93, "d)", transform=axs[1][1].transAxes, size=15, weight='bold')

for y in y_values:
    axs[0][0].axhline(y=y, color='black', linestyle='-', alpha=0.1, zorder=0)
    axs[0][1].axhline(y=y, color='black', linestyle='-', alpha=0.1, zorder=0)
    axs[1][0].axhline(y=y, color='black', linestyle='-', alpha=0.1, zorder=0)
    axs[1][1].axhline(y=y, color='black', linestyle='-', alpha=0.1, zorder=0)

fig.text(0.5, 0.1, "Docking protocol", fontsize=24, ha='center')
lines_labels = [axs[1][1].get_legend_handles_labels()]
lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]

fig.legend(lines, labels, ncol=6, fontsize=16, loc="lower center", bbox_to_anchor = (0, 0.02, 1, 1))

plt.savefig(Path("figures", "epitope_paratope_comparison", "rigid_emref_1-10_rmsdepi_2cols.png"), dpi=300, bbox_inches="tight")


#<<<<<<< HEAD
## PARATOPE RMSD-pLDDT PLOT
#af2_para_plddt = pd.read_csv(Path("data", "AF2_antibody_rmsd_plddt.csv"), sep="\t")
#print(f"af2 paratope pLDDT shape {af2_para_plddt.shape}")
#values_para_rmsd = pd.read_csv(Path("data", "paratope_rmsds.csv"), sep=",")
#af2_para_rmsd = values_para_rmsd.loc[values_para_rmsd["model"] == "AF2"]
#print(f"af2 paratope RMSD shape {af2_para_rmsd.shape}")
##print(af2_para_rmsd["pdb"].values)
##print(af2_para_plddt["pdb"].values)
#assert np.array_equal(af2_para_rmsd["pdb"].values, af2_para_plddt["pdb"].values)
#=======
# PARATOPE RMSD-PLDDT PLOT
af2_para_plddt = pd.read_csv(Path("data", "AF2_antibody_rmsd_plddt_multi_regions.csv"))
af2_para_plddt = af2_para_plddt.query('pdb not in @filter_pdbs').reset_index()
print(af2_para_plddt.shape)

#>>>>>>> 00357db2d0784ad01207e40cb1f2716d18e44946

fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(6, 6))
plt.subplots_adjust(wspace=0.2, hspace=0.1)
plt.suptitle("pLDDT-RMSD correlation", fontsize=24, verticalalignment="center")

print("correlation coefficient for af2 plddt_ave_para_epi = ", np.corrcoef(af2_para_plddt["rmsd_para_epi"], af2_para_plddt["plddt_ave_para_epi"])[0,1])

axs.scatter(af2_para_plddt["rmsd_para_epi"], af2_para_plddt["plddt_ave_para_epi"],color="lightgreen")
axs.set_xlabel("paratope RMSD ($\\AA$)", fontsize=24)
axs.set_ylabel("paratope pLDDT", fontsize=24)
axs.set_ylim((70,101))
axs.set_xlim((-0.1,6.0))

# ticks
axs.tick_params(axis='both', which='major', labelsize=16)

#<<<<<<< HEAD
#for pdb in list(af2_para_rmsd["pdb"]):
#    df_subset = af2_para_rmsd.loc[af2_para_rmsd["pdb"] == pdb]
#    para_rmsd = df_subset["rmsd_paratope"].iloc[0]
#    df_plddt = af2_para_plddt.loc[af2_para_plddt["pdb"] == pdb]
#    para_plddt = df_plddt["para-ave-plddt"].iloc[0]
#    if para_rmsd > 4.0 or para_plddt < 80.0:
#        if pdb in ["7n3c", "7kf0", "7qnw", "7msq"]:
#            axs.text(x=para_rmsd-0.45, y=para_plddt+0.5, s=pdb.upper(), fontsize=14)
#=======
for pdb in list(af2_para_plddt["pdb"]):
    df_subset = af2_para_plddt.loc[af2_para_plddt["pdb"] == pdb]
    para_rmsd = df_subset["rmsd_para_epi"].iloc[0]
    para_plddt = df_subset["plddt_ave_para_epi"].iloc[0]
    if para_rmsd > 4.0 or para_plddt < 75.0:
        if pdb in ["7f7e", "7rfb", "7msq"]:
            axs.text(x=para_rmsd-0.4, y=para_plddt+0.25, s=pdb.upper(), fontsize=12)
#>>>>>>> 00357db2d0784ad01207e40cb1f2716d18e44946
        else:
            axs.text(x=para_rmsd+0.05, y=para_plddt+0.25, s=pdb.upper(), fontsize=12)

plt.tight_layout()
plt.savefig(Path("figures", "epitope_paratope_comparison", f"plddt-rmsd-paraonly.png"), dpi=300, bbox_inches="tight")

# PARATOPE RMSD-CONF SCORE PLOT
abb2_para_conf = pd.read_csv(Path("data", "ABB2_antibody_rmsd_conf.csv"))
abb2_para_conf = abb2_para_conf.query('pdb not in @filter_pdbs').reset_index()
print(f"abb2_para_conf {abb2_para_conf}")
#<<<<<<< HEAD
#assert list(abb2_para_rmsd["pdb"]) == list(abb2_para_conf["pdb"]) # pdbs must be ordered
#assert abb2_para_rmsd.shape[0] == 79
#assert abb2_para_conf.shape[0] == 79

# plot
fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(6, 6))
plt.subplots_adjust(wspace=0.2, hspace=0.1)
plt.suptitle("RMSPE-RMSD correlation", fontsize=24, verticalalignment="center")
print("correlation coefficient for abb2 para_ave_conf = ", np.corrcoef(abb2_para_conf["rmsd_paratope"], abb2_para_conf["para_ave_conf"])[0,1])
axs.scatter(abb2_para_conf["rmsd_paratope"], abb2_para_conf["para_ave_conf"],color="darkgreen")
axs.set_xlabel("paratope RMSD ($\\AA$)", fontsize=24)
axs.set_ylabel("paratope RMSPE", fontsize=24)
axs.set_ylim((-0.1,2.5))
axs.set_xlim((-0.1,7.0))

for pdb in list(abb2_para_conf["pdb"]):
    df_subset = abb2_para_conf.loc[abb2_para_conf["pdb"] == pdb]
    para_rmsd = df_subset["rmsd_paratope"].iloc[0]
#<<<<<<< HEAD
#    df_conf = abb2_para_conf.loc[abb2_para_conf["pdb"] == pdb]
#    para_conf = df_conf["para-ave-conf"].iloc[0]
#    if para_rmsd > 4.0 or para_conf > 1.0:
#        if pdb in ["7msq", "7bnv"]:
#            axs.text(x=para_rmsd-0.45, y=para_conf+0.035, s=pdb.upper(), fontsize=14)
#        elif pdb == "7n4j":
#            axs.text(x=para_rmsd-0.2, y=para_conf+0.035, s=pdb.upper(), fontsize=14)
#=======
    para_conf = df_subset["para_ave_conf"].iloc[0]
    if para_rmsd > 4.0 or para_conf > 1.5:
        if pdb == "7msq":
            axs.text(x=para_rmsd-0.4, y=para_conf+0.07, s=pdb.upper(), fontsize=14)
        elif pdb == "7pr0":
            axs.text(x=para_rmsd-0.1, y=para_conf+0.025, s=pdb.upper(), fontsize=14)
#>>>>>>> 00357db2d0784ad01207e40cb1f2716d18e44946
        else:
            axs.text(x=para_rmsd, y=para_conf+0.025, s=pdb.upper(), fontsize=14)

# ticks
axs.tick_params(axis='both', which='major', labelsize=16)

plt.tight_layout()
plt.savefig(Path("figures", "epitope_paratope_comparison", f"confscore-rmsd-paraonly.png"), dpi=300, bbox_inches="tight")
