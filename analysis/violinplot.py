import pandas as pd
from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np
pdbs_to_exclude = ["6xsw", "7cj2", "7e9b", "7m3n", "7n3c", "7n3d", "7o52", "7or9"]
# CDR H3
unb_data = pd.read_csv(Path("data", "unbound_rmsd.csv"), sep="\t")
ref_data = pd.read_csv(Path("data", "ref_data.csv"), sep=" ")
ref_data = ref_data.where(~ref_data["pdb"].isin(pdbs_to_exclude)).dropna()
assert np.unique(ref_data["pdb"]).shape[0] == 71

cdr_h3_abb = ref_data.where(ref_data["method"] == "ABB2").dropna()["cdrh3"]
cdr_h3_abl = ref_data.where(ref_data["method"] == "ABlooper").dropna()["cdrh3"]
cdr_h3_af2 = ref_data.where(ref_data["method"] == "AF2").dropna()["cdrh3"]
cdr_h3_ig = ref_data.where(ref_data["method"] == "IgFold").dropna()["cdrh3"]
# stats
median_abb = np.median(cdr_h3_abb)
median_abl = np.median(cdr_h3_abl)
median_af2 = np.median(cdr_h3_af2)
median_ig = np.median(cdr_h3_ig)
mean_abb = np.mean(cdr_h3_abb)
mean_abl = np.mean(cdr_h3_abl)
mean_af2 = np.mean(cdr_h3_af2)
mean_ig = np.mean(cdr_h3_ig)
std_abb = np.std(cdr_h3_abb)
std_abl = np.std(cdr_h3_abl)
std_af2 = np.std(cdr_h3_af2)
std_ig = np.std(cdr_h3_ig)
print(f"cdrh3 median abb {median_abb} abl {median_abl} af2 {median_af2} ig {median_ig}")
print(f"cdrh3 mean abb {mean_abb} +- {std_abb} abl {mean_abl} +- {std_abl} af2 {mean_af2} +- {std_af2} ig {mean_ig} +- {std_ig}")
# extracting abens best model
abens_list = ["ABB2_ensemble_0", "ABB2_ensemble_1", "ABB2_ensemble_2", "ABB2_ensemble_3"]
cdr_h3_abens = []
for pdb in np.unique(ref_data["pdb"]):
    abens_df = ref_data.where(ref_data["pdb"] == pdb).dropna().where(ref_data["method"].isin(abens_list)).dropna()
    cdr_h3_abens.append(min(abens_df["cdrh3"]))
# extracting ens best model, namely best model in the 'ABB', 'ABL', 'AF2', "IG"
cdr_h3_ens = []
for pdb in np.unique(ref_data["pdb"]):
    ens_df = ref_data.where(ref_data["pdb"] == pdb).dropna().where(ref_data["method"].isin(["ABB2", "ABlooper", "AF2", "IgFold"])).dropna()
    cdr_h3_ens.append(min(ens_df["cdrh3"]))

fig, ax = plt.subplots(1, 1, figsize=(9, 9))
# Create a list of the data to be plotted
data = [cdr_h3_abb, cdr_h3_abl, cdr_h3_af2, cdr_h3_ig, cdr_h3_abens, cdr_h3_ens]
# Set the colors for the violins based on the category
colors = ['Blue', 'Cyan', 'Purple', "Orange", "Red", "Darkgreen"]
# Create the violin plot
plots = ax.violinplot(data, showmedians=True, showextrema=False, widths=1)

# Set the color of the violin patches
for pc, color in zip(plots['bodies'], colors):
    pc.set_facecolor(color)

# Set the color of the median lines
plots['cmedians'].set_colors("black")
# write medians
medians = [np.median(d) for d in data]
for i, median in enumerate(medians):
    ax.text(i + 1, median+0.1, f'${median:.2f}\ \AA$', ha='center', va='bottom', color='black', fontsize=14)

# Set the labels
ax.set_xticks([1, 2, 3, 4, 5, 6], labels=['ABB', 'ABL', 'AF2', "IG", "ABBE", "ENS"], fontsize=16)
ax.tick_params(axis='y', labelsize=16)

ax.set_ylabel("CDRH3 RMSD ($\AA$)",fontsize=16)
#ax.set_xlabel("Method",fontsize=16)
plt.savefig(Path("figures", "violinplots", "cdrh3_violinplot.png"), dpi=400, bbox_inches="tight")

# unbound H3
fig, ax = plt.subplots(1, 1, figsize=(9, 9))
# Create a list of the data to be plotted
data_h3_unb = [cdr_h3_abb, cdr_h3_abl, cdr_h3_af2, cdr_h3_ig, cdr_h3_abens, cdr_h3_ens, unb_data["h3_rmsd"]]
# Set the colors for the violins based on the category
colors = ['Blue', 'Cyan', 'Purple', "Orange", "Red", "Darkgreen", "Gray"]
# Create the violin plot
plots = ax.violinplot(data_h3_unb, showmedians=True, showextrema=False, widths=1)

# Set the color of the violin patches
for pc, color in zip(plots['bodies'], colors):
    pc.set_facecolor(color)

# Set the color of the median lines
plots['cmedians'].set_colors("black")
# write medians
medians = [np.median(d) for d in data_h3_unb]
for i, median in enumerate(medians):
    ax.text(i + 1, median+0.1, f'${median:.2f}\ \AA$', ha='center', va='bottom', color='black', fontsize=14)

# Set the labels
ax.set_xticks([1, 2, 3, 4, 5, 6, 7], labels=['ABB', 'ABL', 'AF2', "IG", "ABBE", "ENS", "UNB"], fontsize=16)
ax.tick_params(axis='y', labelsize=16)

ax.set_ylabel("CDRH3 RMSD ($\AA$)",fontsize=16)
#ax.set_xlabel("Method",fontsize=16)
plt.savefig(Path("figures", "violinplots", "cdrh3_violinplot_unb.png"), dpi=400, bbox_inches="tight")

# PARATOPE
abens_list_para = ["ABB2_ensemble_rank0", "ABB2_ensemble_rank1", "ABB2_ensemble_rank2", "ABB2_ensemble_rank3"]
paratope_rmsd = pd.read_csv(Path("data", "paratope_rmsds.csv"), sep=",")
paratope_rmsd = paratope_rmsd.where(~paratope_rmsd["pdb"].isin(pdbs_to_exclude)).dropna()
assert np.unique(paratope_rmsd["pdb"]).shape[0] == 71
para_rmsd_abb = paratope_rmsd.where(paratope_rmsd["model"] == "ABodyBuilder2").dropna()["rmsd_paratope"]
para_rmsd_abl = paratope_rmsd.where(paratope_rmsd["model"] == "ABlooper").dropna()["rmsd_paratope"]
para_rmsd_af2 = paratope_rmsd.where(paratope_rmsd["model"] == "AF2").dropna()["rmsd_paratope"]
para_rmsd_ig = paratope_rmsd.where(paratope_rmsd["model"] == "IgFold").dropna()["rmsd_paratope"]
# extracting stats
para_median_abb = np.median(para_rmsd_abb)
para_median_abl = np.median(para_rmsd_abl)
para_median_af2 = np.median(para_rmsd_af2)
para_median_ig = np.median(para_rmsd_ig)
para_mean_abb = np.mean(para_rmsd_abb)
para_mean_abl = np.mean(para_rmsd_abl)
para_mean_af2 = np.mean(para_rmsd_af2)
para_mean_ig = np.mean(para_rmsd_ig)
para_std_abb = np.std(para_rmsd_abb)
para_std_abl = np.std(para_rmsd_abl)
para_std_af2 = np.std(para_rmsd_af2)
para_std_ig = np.std(para_rmsd_ig)

print(f"ABodyBuilder2 paratope: median {para_median_abb:.2f} mean {para_mean_abb:.2f} std {para_std_abb:.2f}")
print(f"ABlooper paratope: median {para_median_abl:.2f} mean {para_mean_abl:.2f} std {para_std_abl:.2f}")
print(f"AF2 paratope: median {para_median_af2:.2f} mean {para_mean_af2:.2f} std {para_std_af2:.2f}")
print(f"IgFold paratope: median {para_median_ig:.2f} mean {para_mean_ig:.2f} std {para_std_ig:.2f}")

# extracting ens best model, namely best model in the 'ABB', 'ABL', 'AF2', "IG"
para_rmsd_ens = []
for pdb in np.unique(paratope_rmsd["pdb"]):
    ens_df = paratope_rmsd.where(paratope_rmsd["pdb"] == pdb).dropna().where(paratope_rmsd["model"].isin(["ABodyBuilder2", "ABlooper", "AF2", "IgFold"])).dropna()
    para_rmsd_ens.append(min(ens_df["rmsd_paratope"]))
para_rmsd_abens = []
for pdb in np.unique(ref_data["pdb"]):
    abens_df = paratope_rmsd.where(paratope_rmsd["pdb"] == pdb).dropna().where(paratope_rmsd["model"].isin(abens_list_para)).dropna()
    para_rmsd_abens.append(min(abens_df["rmsd_paratope"]))
#print(para_rmsd_ens)

fig, ax = plt.subplots(1, 1, figsize=(9, 9))
# Create a list of the data to be plotted
data = [para_rmsd_abb, para_rmsd_abl, para_rmsd_af2, para_rmsd_ig, para_rmsd_abens, para_rmsd_ens]
# Set the colors for the violins based on the category
colors = ['Blue', 'Cyan', 'Purple', "Orange", "Red", "Darkgreen"]
# Create the violin plot
plots = ax.violinplot(data, showmedians=True, showextrema=False, widths=1)

# Set the color of the violin patches
for pc, color in zip(plots['bodies'], colors):
    pc.set_facecolor(color)

# Set the color of the median lines
plots['cmedians'].set_colors("black")
# write medians
medians = [np.median(d) for d in data]
for i, median in enumerate(medians):
    ax.text(i + 1, median+0.1, f'${median:.2f}\ \AA$', ha='center', va='bottom', color='black', fontsize=14)

# Set the labels
ax.set_xticks([1, 2, 3, 4, 5, 6], labels=['ABB', 'ABL', 'AF2', "IG", "ABBE", "ENS"], fontsize=16)
ax.tick_params(axis='y', labelsize=16)

ax.set_ylabel("Paratope RMSD ($\AA$)",fontsize=16)
#ax.set_xlabel("Method",fontsize=16)
plt.savefig(Path("figures", "violinplots", "paratope_violinplot.png"), dpi=400, bbox_inches="tight")

# adding unbound data
#para_rmsd_unb = [0.4509972259629457,0.4390521344154917,0.7788158223885456 ,0.46914223178365955,
#1.5749394164436636,1.8495151759984525,0.4081729543513871,0.9956042358221723,
#0.5241969637819592,0.7288005318925351,1.7681836908830881,0.39291736179706477,
#0.5958420766833828,0.5037569828402556,0.47233065317342027,0.5391822131740583]

fig, ax = plt.subplots(1, 1, figsize=(9, 9))
# Create a list of the data to be plotted
data_para_unb = [para_rmsd_abb, para_rmsd_abl, para_rmsd_af2, para_rmsd_ig, para_rmsd_abens, para_rmsd_ens, unb_data["paratope_rmsd"]]
# Set the colors for the violins based on the category
colors = ['Blue', 'Cyan', 'Purple', "Orange", "Red", "Darkgreen", "Gray"]
# Create the violin plot
plots = ax.violinplot(data_para_unb, showmedians=True, showextrema=False, widths=1)

# Set the color of the violin patches
for pc, color in zip(plots['bodies'], colors):
    pc.set_facecolor(color)

# Set the color of the median lines
plots['cmedians'].set_colors("black")
# write medians
medians = [np.median(d) for d in data_para_unb]
for i, median in enumerate(medians):
    ax.text(i + 1, median+0.1, f'${median:.2f}\ \AA$', ha='center', va='bottom', color='black', fontsize=14)

# Set the labels
ax.set_xticks([1, 2, 3, 4, 5, 6, 7], labels=['ABB', 'ABL', 'AF2', "IG", "ABBE", "ENS", "UNB"], fontsize=16)
ax.tick_params(axis='y', labelsize=16)

ax.set_ylabel("Paratope RMSD ($\AA$)",fontsize=16)
#ax.set_xlabel("Method",fontsize=16)
plt.savefig(Path("figures", "violinplots", "paratope_violinplot_unb.png"), dpi=400, bbox_inches="tight")

# joined unbound plot
fig, ax = plt.subplots(1, 3, figsize=(18, 9), width_ratios=[7, 7, 4])
plt.subplots_adjust(wspace=0.2, hspace=0.1)
plots = ax[0].violinplot(data_h3_unb, showmedians=True, showextrema=False, widths=1)

# Set the color of the violin patches
for pc, color in zip(plots['bodies'], colors):
    pc.set_facecolor(color)
plots['cmedians'].set_colors("black")

plots = ax[1].violinplot(data_para_unb, showmedians=True, showextrema=False, widths=1)
plots['cmedians'].set_colors("black")
# Set the color of the violin patches
for pc, color in zip(plots['bodies'], colors):
    pc.set_facecolor(color)

# antigen violinplot
antigen_rmsd = pd.read_csv(Path("data", "AF2_antigen_rmsd.csv"), sep=",")
# exclude the pdbs
antigen_rmsd = antigen_rmsd.where(~antigen_rmsd["pdb"].isin(pdbs_to_exclude)).dropna()
assert antigen_rmsd.shape[0] == 71
# print ratios for the antigen
ratio_20_epitope = np.where(antigen_rmsd["rmsd_epitope"]<2.0)[0].shape[0]/antigen_rmsd.shape[0]
ratios_20_antrmsd = np.where(antigen_rmsd["rmsd"]<2.0)[0].shape[0]/antigen_rmsd.shape[0]
print(f"ratios_20_epitope {ratio_20_epitope}")
print(f"ratios_20_antrmsd {ratios_20_antrmsd}")

data_antigen = [antigen_rmsd["rmsd"], antigen_rmsd["rmsd_epitope"]]
plots = ax[2].violinplot(data_antigen, showmedians=True, showextrema=False, widths=1)
print(f"antigen mean {np.mean(antigen_rmsd['rmsd'])} +- {np.std(antigen_rmsd['rmsd'])}")
print(f"antigen mean epitope {np.mean(antigen_rmsd['rmsd_epitope'])} +- {np.std(antigen_rmsd['rmsd_epitope'])}")
# Set the color of the median lines
plots['cmedians'].set_colors("black")
for pc, color in zip(plots['bodies'], colors):
    pc.set_facecolor("Purple")
# write medians
medians_h3 = [np.median(d) for d in data_h3_unb]
medians_para = [np.median(d) for d in data_para_unb]
medians_antigen = [np.median(d) for d in data_antigen]

ratios_30 = [np.where(np.array(d)<3.0)[0].shape[0]/len(d) for d in data_para_unb]
ratios_25 = [np.where(np.array(d)<2.5)[0].shape[0]/len(d) for d in data_para_unb]
print(f"ratios {ratios_30}")
print(f"ratios {ratios_25}")

ax[1].axhline(3.0, color="black", linestyle="--")
ax[1].axhline(2.5, color="black", linestyle="--")

for i, median in enumerate(medians_h3):
    ax[0].text(i + 1, median+0.1, f'${median:.2f}\ \AA$', ha='center', va='bottom', color='black', fontsize=12)
for i, median in enumerate(medians_para):
    ax[1].text(i + 1, median+0.05, f'${median:.2f}\ \AA$', ha='center', va='bottom', color='black', fontsize=12)
for i, median in enumerate(medians_antigen):
    ax[2].text(i + 1, median+0.1, f'${median:.2f}\ \AA$', ha='center', va='bottom', color='black', fontsize=12)

# percentage of structures with RMSD < 3.0
for i, ratio in enumerate(ratios_30):
    ax[1].text(i+1, 3.05, f'{ratios_30[i]*100:.0f}%', ha='center', va='bottom', color='black', fontsize=12)
for i, ratio in enumerate(ratios_25):
    ax[1].text(i+1, 2.55, f'{ratios_25[i]*100:.0f}%', ha='center', va='bottom', color='black', fontsize=12)

# Set the labels
ax[0].set_xticks([1, 2, 3, 4, 5, 6, 7], labels=['ABB', 'ABL', 'AF2', "IG", "ABBE", "ENS", "UNB"], fontsize=16)
ax[0].tick_params(axis='y', labelsize=16)
ax[1].set_xticks([1, 2, 3, 4, 5, 6, 7], labels=['ABB', 'ABL', 'AF2', "IG", "ABBE", "ENS", "UNB"], fontsize=16)
ax[1].tick_params(axis='y', labelsize=16)
ax[2].set_xticks([1, 2], labels=['Full', 'Epitope'], fontsize=16)
ax[2].tick_params(axis='y', labelsize=16)

ax[0].set_ylabel("CDRH3 RMSD ($\AA$)",fontsize=16)
ax[1].set_ylabel("Paratope RMSD ($\AA$)",fontsize=16)
ax[2].set_ylabel("Antigen RMSD ($\AA$)",fontsize=16)
ax[0].set_ylim((0.0,10.5))
ax[1].set_ylim((0.0,7.0))
ax[2].set_ylim((0.0,18.0))

plt.savefig(Path("figures", "violinplots", "violinplot_full.png"), dpi=400, bbox_inches="tight")

# separate antigen violinplot
fig, ax = plt.subplots(1, 1, figsize=(4, 9))
plots = ax.violinplot(data_antigen, showmedians=True, showextrema=False, widths=1)
# Set the color of the median lines
plots['cmedians'].set_colors("black")
for pc, color in zip(plots['bodies'], colors):
    pc.set_facecolor("Purple")
for i, median in enumerate(medians_antigen):
    ax.text(i + 1, median+0.1, f'${median:.2f}\ \AA$', ha='center', va='bottom', color='black', fontsize=12)

ax.set_xticks([1, 2], labels=['Full', 'Epitope'], fontsize=16)
ax.tick_params(axis='y', labelsize=16)
ax.set_ylabel("Antigen RMSD ($\AA$)",fontsize=16)
ax.set_ylim((0.0,10.5))
ax.set_ylim((0.0,7.0))
ax.set_ylim((0.0,18.0))
plt.savefig(Path("figures", "violinplots", "violinplot_antigen.png"), dpi=400, bbox_inches="tight")

