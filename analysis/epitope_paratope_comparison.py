import pandas as pd
from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np
from functions import get_sorted_runs, load_data

plt.rcParams["font.family"] = "Helvetica"

cat_dict = {"Para-Epi": "Para-Epi",
            "CDR-EpiVag": "CDR-VagueEpi" ,
            "CDR-EpiVag-AA" : "CDR-VagueEpi-AA"}
XTICKS = ["ABB", "ABL", "AF2", "IG", "ABBE", "AF2E", "IGE", "ENS", "ENSNOAF2", "CLE", "ENS196-48", "ENS196-CLT"]

# LOAD DATA
af2_values = pd.read_csv(Path("data", "AF2_antigen_rmsd_plddt.csv"), sep=",")
print(f"af2_values {af2_values}")

rigidbody_capri, rigidbody_capri_bound, emref_capri, emref_capri_bound, df_ss_flexref, df_ss_bound_flexref, zdock_ss, emref_rigid_capri, af2multimer_ss = load_data()


# PLOT CORRELATION ON TWO PLOT
fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
plt.subplots_adjust(wspace=0.2, hspace=0.1)
plt.suptitle("pLDDT-RMSD correlation", fontsize=24, verticalalignment="center")

axs[0].scatter(af2_values["epi-rmsd"], af2_values["epi-plddt"])
axs[0].set_xlabel("epitope RMSD ($\\AA$)", fontsize=24)
axs[0].set_ylabel("pLDDT", fontsize=24)
axs[0].set_ylim((50,101))
axs[0].set_xlim((-0.5,20.0))
axs[1].scatter(af2_values["antigen-rmsd"], af2_values["antigen-plddt"])
axs[1].set_xlabel("full antigen RMSD ($\\AA$)", fontsize=24)
axs[1].set_ylim((50,101))
axs[1].set_xlim((-0.5,20.0))
# ticks
axs[1].set_yticks([])
axs[0].tick_params(axis='both', which='major', labelsize=16)
axs[1].tick_params(axis='both', which='major', labelsize=16)

for pdb in list(af2_values["pdb"]):
    df_subset = af2_values.loc[af2_values["pdb"] == pdb]
    epi_rmsd = df_subset["epi-rmsd"].iloc[0]
    full_rmsd = df_subset["antigen-rmsd"].iloc[0]
    epi_plddt = df_subset["epi-plddt"].iloc[0]
    full_plddt = df_subset["antigen-plddt"].iloc[0]
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

axs.scatter(af2_values["epi-rmsd"], af2_values["epi-plddt"])
axs.set_xlabel("epitope RMSD ($\\AA$)", fontsize=24)
axs.set_ylabel("epitope pLDDT", fontsize=24)
axs.set_ylim((50,101))
axs.set_xlim((-0.2,6.5))

# ticks
axs.tick_params(axis='both', which='major', labelsize=16)

for pdb in list(af2_values["pdb"]):
    df_subset = af2_values.loc[af2_values["pdb"] == pdb]
    epi_rmsd = df_subset["epi-rmsd"].iloc[0]
    full_rmsd = df_subset["antigen-rmsd"].iloc[0]
    epi_plddt = df_subset["epi-plddt"].iloc[0]
    full_plddt = df_subset["antigen-plddt"].iloc[0]
    if epi_rmsd > 5.0 or epi_plddt < 70.0:
        axs.text(x=epi_rmsd-0.2, y=epi_plddt+0.75, s=pdb.upper(), fontsize=14)    

print(f"correlation coefficient for af2 epi-ave-plddt = {np.corrcoef(af2_values['epi-rmsd'], af2_values['epi-plddt'])[0,1]}")
plt.tight_layout()
plt.savefig(Path("figures", "epitope_paratope_comparison", f"plddt-rmsd-epionly.png"), dpi=300, bbox_inches="tight")


# EPITOPE RMSD PLOT
# generating ave_af2rmsdepi_histo dictionary
ave_af2rmsdepi_histo = {}
ave_af2rmsdepi_histo["high"] = af2_values.where(af2_values["cat_epi_plddt"] == "high").dropna()["pdb"]
ave_af2rmsdepi_histo["good"] = af2_values.where(af2_values["cat_epi_plddt"] == "good").dropna()["pdb"]
ave_af2rmsdepi_histo["bad"] = af2_values.where(af2_values["cat_epi_plddt"] == "bad").dropna()["pdb"]
print(f"ave_af2rmsdepi_histo = {ave_af2rmsdepi_histo}")
tops = [1, 5, 10, 20, 48]
tops_dict = {el : 0 for el in tops}

af2_rmsdepi = {"Para-Epi" : {}, "CDR-EpiVag" : {}, "CDR-EpiVag-AA": {}}

cat_runs = {"Para-Epi" : {}, "CDR-EpiVag" : {}, "CDR-EpiVag-AA": {}}
for cat in cat_runs.keys():
    cat_runs[cat] = get_sorted_runs(cat)

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
    runs = get_sorted_runs(cat)
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
    runs = get_sorted_runs(cat)
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
    axs[1][col].set_xticks(r + width/2, XTICKS, fontsize=16, rotation = 90)
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
    runs = get_sorted_runs(cat)
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
    runs = get_sorted_runs(cat)
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
    print(f"mod_dict_1_emref[good] {mod_dict_1_emref['good']}")
    print(f"mod_dict_1_emref[high] {mod_dict_1_emref['high']}")
    print(f"mod_dict_10_emref[good] {mod_dict_10_emref['good']}")
    print(f"mod_dict_10_emref[high] {mod_dict_10_emref['high']}")
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
    axs[1][col].set_xticks(r + width/2, XTICKS, fontsize=16, rotation = 90)
    col+=1

axs[0][0].set_title(cat_dict["Para-Epi"], fontsize=16)
axs[0][1].set_title(cat_dict["CDR-EpiVag-AA"], fontsize=16)

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

#=======
# PARATOPE RMSD-PLDDT PLOT
paratope_data = pd.read_csv(Path("data", "paratope_h3_rmsds.csv"))
af2_para_plddt = paratope_data.loc[paratope_data["model"] == "AF2"]
assert af2_para_plddt.shape[0] == 82
print(f"af2_para_plddt {af2_para_plddt}")

fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(6, 6))
plt.subplots_adjust(wspace=0.2, hspace=0.1)
plt.suptitle("pLDDT-RMSD correlation", fontsize=24, verticalalignment="center")

print("correlation coefficient for af2 para-plddt = ", np.corrcoef(af2_para_plddt["paratope_rmsd"], af2_para_plddt["paratope_conf"])[0,1])

axs.scatter(af2_para_plddt["paratope_rmsd"], af2_para_plddt["paratope_conf"],color="lightgreen")
axs.set_xlabel("paratope RMSD ($\\AA$)", fontsize=24)
axs.set_ylabel("paratope pLDDT", fontsize=24)
axs.set_ylim((70,101))
axs.set_xlim((-0.1,6.0))

# ticks
axs.tick_params(axis='both', which='major', labelsize=16)

for pdb in list(af2_para_plddt["pdb"]):
    df_subset = af2_para_plddt.loc[af2_para_plddt["pdb"] == pdb]
    para_rmsd = df_subset["paratope_rmsd"].iloc[0]
    para_plddt = df_subset["paratope_conf"].iloc[0]
    if para_rmsd > 4.0 or para_plddt < 75.0:
        axs.text(x=para_rmsd-0.4, y=para_plddt+0.25, s=pdb.upper(), fontsize=12)

plt.tight_layout()
plt.savefig(Path("figures", "epitope_paratope_comparison", f"plddt-rmsd-paraonly.png"), dpi=300, bbox_inches="tight")

# PARATOPE RMSD-CONF SCORE PLOT
abb2_para_conf = paratope_data.loc[paratope_data["model"] == "ABodyBuilder2"]

# plot
fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(6, 6))
plt.subplots_adjust(wspace=0.2, hspace=0.1)
plt.suptitle("RMSPE-RMSD correlation", fontsize=24, verticalalignment="center")
print("correlation coefficient for abb2 para_ave_conf = ", np.corrcoef(abb2_para_conf["paratope_rmsd"], abb2_para_conf["paratope_conf"])[0,1])
axs.scatter(abb2_para_conf["paratope_rmsd"], abb2_para_conf["paratope_conf"],color="darkgreen")
axs.set_xlabel("paratope RMSD ($\\AA$)", fontsize=24)
axs.set_ylabel("paratope RMSPE", fontsize=24)
axs.set_ylim((-0.1,2.5))
axs.set_xlim((-0.1,8.0))

for pdb in list(abb2_para_conf["pdb"]):
    df_subset = abb2_para_conf.loc[abb2_para_conf["pdb"] == pdb]
    para_rmsd = df_subset["paratope_rmsd"].iloc[0]
    para_conf = df_subset["paratope_conf"].iloc[0]
    if para_rmsd > 4.0 or para_conf > 1.5:
        if pdb == "7l7e_C0-D0":
            axs.text(x=para_rmsd-1.7, y=para_conf+0.025, s=pdb.upper(), fontsize=14)
        elif pdb == "7bnv_H0-L0":
            axs.text(x=para_rmsd-1.5, y=para_conf+0.025, s=pdb.upper(), fontsize=14)
        elif pdb == "7mzm_H0-L0":
            axs.text(x=para_rmsd+0.025, y=para_conf-0.05, s=pdb.upper(), fontsize=14)
        else:
            axs.text(x=para_rmsd, y=para_conf+0.025, s=pdb.upper(), fontsize=14)

# ticks
axs.tick_params(axis='both', which='major', labelsize=16)

plt.tight_layout()
plt.savefig(Path("figures", "epitope_paratope_comparison", f"confscore-rmsd-paraonly.png"), dpi=300, bbox_inches="tight")
