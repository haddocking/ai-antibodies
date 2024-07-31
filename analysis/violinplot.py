import pandas as pd
from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np
NPDBS = 82

XTICK_ANTIBODY = ['ABB', 'ABL', 'AF2', "IG", "ABBE", "AFE", "IGE", "ENS", "CLE"]
# CDR H3
unb_data = pd.read_csv(Path("data", "unbound_rmsd.csv"), sep="\t")
ref_data = pd.read_csv(Path("data", "paratope_h3_rmsds.csv"), sep=",")
assert np.unique(ref_data["pdb"]).shape[0] == NPDBS

cdr_h3_abb = ref_data.where(ref_data["model"] == "ABodyBuilder2").dropna()["h3_rmsd"]
cdr_h3_abl = ref_data.where(ref_data["model"] == "ABlooper").dropna()["h3_rmsd"]
cdr_h3_af2 = ref_data.where(ref_data["model"] == "AF2").dropna()["h3_rmsd"]
cdr_h3_ig = ref_data.where(ref_data["model"] == "IgFold").dropna()["h3_rmsd"]
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
print(f"cdrh3 mean abb {mean_abb:.3f} +- {std_abb:.3f} abl {mean_abl:.3f} +- {std_abl:.3f} af2 {mean_af2:.3f} +- {std_af2:.3f} ig {mean_ig:.3f} +- {std_ig:.3f}")
# extracting abens best model
abens_list = ["ABB2_ensemble_rank0", "ABB2_ensemble_rank1", "ABB2_ensemble_rank2", "ABB2_ensemble_rank3"]
cdr_h3_abens = []
for pdb in np.unique(ref_data["pdb"]):
    abens_df = ref_data.where(ref_data["pdb"] == pdb).dropna().where(ref_data["model"].isin(abens_list)).dropna()
    cdr_h3_abens.append(min(abens_df["h3_rmsd"]))
print(f"ABBE cdrh3 {np.median(cdr_h3_abens):.2f} +- {np.std(cdr_h3_abens):.2f}")
afens_list = ["AF2_ensemble_rank0", "AF2_ensemble_rank1", "AF2_ensemble_rank2", "AF2_ensemble_rank3"]
cdr_h3_afens = []
for pdb in np.unique(ref_data["pdb"]):
    afens_df = ref_data.where(ref_data["pdb"] == pdb).dropna().where(ref_data["model"].isin(afens_list)).dropna()
    cdr_h3_afens.append(min(afens_df["h3_rmsd"]))
print(f"AFENS cdrh3 {np.median(cdr_h3_afens):.2f} +- {np.std(cdr_h3_afens):.2f}")
# igens
igens_list = ["IG_ensemble_rank0", "IG_ensemble_rank1", "IG_ensemble_rank2", "IG_ensemble_rank3"]
cdr_h3_igens = []
for pdb in np.unique(ref_data["pdb"]):
    igens_df = ref_data.where(ref_data["pdb"] == pdb).dropna().where(ref_data["model"].isin(igens_list)).dropna()
    cdr_h3_igens.append(min(igens_df["h3_rmsd"]))
print(f"IGENS cdrh3 {np.median(cdr_h3_igens):.2f} +- {np.std(cdr_h3_igens):.2f}")
# extracting ens best model, namely best model in the 'ABB', 'ABL', 'AF2', "IG"
cdr_h3_ens = []
for pdb in np.unique(ref_data["pdb"]):
    ens_df = ref_data.where(ref_data["pdb"] == pdb).dropna().where(ref_data["model"].isin(["ABodyBuilder2", "ABlooper", "AF2", "IgFold"])).dropna()
    cdr_h3_ens.append(min(ens_df["h3_rmsd"]))
# clens
cdr_h3_clens = ref_data.where(ref_data["model"] == "clens").dropna()["h3_rmsd"]
print(f"clens cdrh3 {np.median(cdr_h3_clens):.2f} +- {np.std(cdr_h3_clens):.2f}")

fig, ax = plt.subplots(1, 1, figsize=(9, 9))
# Create a list of the data to be plotted
data = [cdr_h3_abb, cdr_h3_abl, cdr_h3_af2, cdr_h3_ig, cdr_h3_abens, cdr_h3_afens, cdr_h3_igens, cdr_h3_ens]
# Set the colors for the violins based on the category
colors = ['Blue', 'Cyan', 'Purple', "Orange", "Red", "Darkgreen", "Gray", "coral", "Lightgreen"]
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
ax.set_xticks([1, 2, 3, 4, 5, 6, 7, 8, 9], labels=XTICK_ANTIBODY, fontsize=16)
ax.tick_params(axis='y', labelsize=16)

ax.set_ylabel("CDRH3 RMSD ($\AA$)",fontsize=16)
#ax.set_xlabel("Method",fontsize=16)
plt.savefig(Path("figures", "violinplots", "cdrh3_violinplot.png"), dpi=400, bbox_inches="tight")

# unbound H3
fig, ax = plt.subplots(1, 1, figsize=(9, 9))
# Create a list of the data to be plotted
data_h3_unb = [cdr_h3_abb, cdr_h3_abl, cdr_h3_af2, cdr_h3_ig, cdr_h3_abens, cdr_h3_afens, cdr_h3_igens, cdr_h3_ens, cdr_h3_clens, unb_data["h3_rmsd"]]
# Set the colors for the violins based on the category
colors.extend(["Gray"])
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
ax.set_xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], labels=XTICK_ANTIBODY + ["UNB"], fontsize=16)
ax.tick_params(axis='y', labelsize=16)

ax.set_ylabel("CDRH3 RMSD ($\AA$)",fontsize=16)
#ax.set_xlabel("Method",fontsize=16)
plt.savefig(Path("figures", "violinplots", "cdrh3_violinplot_unb.png"), dpi=400, bbox_inches="tight")

# PARATOPE
para_rmsd_abb = ref_data.where(ref_data["model"] == "ABodyBuilder2").dropna()["paratope_rmsd"]
para_rmsd_abl = ref_data.where(ref_data["model"] == "ABlooper").dropna()["paratope_rmsd"]
para_rmsd_af2 = ref_data.where(ref_data["model"] == "AF2").dropna()["paratope_rmsd"]
para_rmsd_ig = ref_data.where(ref_data["model"] == "IgFold").dropna()["paratope_rmsd"]
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
for pdb in np.unique(ref_data["pdb"]):
    ens_df = ref_data.where(ref_data["pdb"] == pdb).dropna().where(ref_data["model"].isin(["ABodyBuilder2", "ABlooper", "AF2", "IgFold"])).dropna()
    para_rmsd_ens.append(min(ens_df["paratope_rmsd"]))
para_rmsd_abens = []
for pdb in np.unique(ref_data["pdb"]):
    abens_df = ref_data.where(ref_data["pdb"] == pdb).dropna().where(ref_data["model"].isin(abens_list)).dropna()
    para_rmsd_abens.append(min(abens_df["paratope_rmsd"]))
para_rmsd_afens = []
for pdb in np.unique(ref_data["pdb"]):
    afens_df = ref_data.where(ref_data["pdb"] == pdb).dropna().where(ref_data["model"].isin(afens_list)).dropna()
    para_rmsd_afens.append(min(afens_df["paratope_rmsd"]))
para_rmsd_igens = []
for pdb in np.unique(ref_data["pdb"]):
    igens_df = ref_data.where(ref_data["pdb"] == pdb).dropna().where(ref_data["model"].isin(igens_list)).dropna()
    para_rmsd_igens.append(min(igens_df["paratope_rmsd"]))

print(f"ABBE paratope {np.median(para_rmsd_abens):.2f} +- {np.std(para_rmsd_abens):.2f}")
print(f"AFENS paratope {np.median(para_rmsd_afens):.2f} +- {np.std(para_rmsd_afens):.2f}")
print(f"IGENS paratope {np.median(para_rmsd_igens):.2f} +- {np.std(para_rmsd_igens):.2f}")
print(f"ENS paratope {np.median(para_rmsd_ens):.2f} +- {np.std(para_rmsd_ens):.2f}")
para_rmsd_clens = ref_data.where(ref_data["model"] == "clens").dropna()["paratope_rmsd"]
para_mean_clens = np.mean(para_rmsd_clens)
para_median_clens = np.median(para_rmsd_clens)
para_std_clens = np.std(para_rmsd_clens)
print(f"clens paratope: median {para_median_clens:.2f} mean {para_mean_clens:.2f} std {para_std_clens:.2f}")

fig, ax = plt.subplots(1, 1, figsize=(9, 9))
# Create a list of the data to be plotted
data = [para_rmsd_abb, para_rmsd_abl, para_rmsd_af2, para_rmsd_ig, para_rmsd_abens, para_rmsd_afens, para_rmsd_igens, para_rmsd_ens, para_rmsd_clens]
# Set the colors for the violins based on the category
colors = ['Blue', 'Cyan', 'Purple', "Orange", "Red", "Darkgreen", "Gray", "coral", "Lightgreen"]
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
ax.set_xticks([1, 2, 3, 4, 5, 6, 7, 8, 9], labels=XTICK_ANTIBODY, fontsize=16)
ax.tick_params(axis='y', labelsize=16)

ax.set_ylabel("Paratope RMSD ($\AA$)",fontsize=16)
#ax.set_xlabel("Method",fontsize=16)
plt.savefig(Path("figures", "violinplots", "paratope_violinplot.png"), dpi=400, bbox_inches="tight")

fig, ax = plt.subplots(1, 1, figsize=(9, 9))
# Create a list of the data to be plotted
data_para_unb = [para_rmsd_abb, para_rmsd_abl, para_rmsd_af2, para_rmsd_ig, para_rmsd_abens, para_rmsd_afens, para_rmsd_igens, para_rmsd_ens, para_rmsd_clens, unb_data["paratope_rmsd"]]
# Set the colors for the violins based on the category
colors.extend(["Gray"])
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
ax.set_xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], labels=XTICK_ANTIBODY + ["UNB"], fontsize=16)
ax.tick_params(axis='y', labelsize=16)

ax.set_ylabel("Paratope RMSD ($\AA$)",fontsize=16)
#ax.set_xlabel("Method",fontsize=16)
plt.savefig(Path("figures", "violinplots", "paratope_violinplot_unb.png"), dpi=400, bbox_inches="tight")

# joined unbound plot
fig, ax = plt.subplots(1, 3, figsize=(21, 9), width_ratios=[8.5, 8.5, 4])
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
antigen_rmsd = pd.read_csv(Path("data", "AF2_antigen_rmsd_plddt.csv"), sep=",")
# exclude the pdbs
assert antigen_rmsd.shape[0] == NPDBS
# print ratios for the antigen
lower_20_epitope = np.where(antigen_rmsd["epi-rmsd"]<2.0)[0].shape[0]
#ratios_20_antrmsd = np.where(antigen_rmsd["rmsd"]<2.0)[0].shape[0]/antigen_rmsd.shape[0]
print(f"lower_20_epitope {lower_20_epitope} ratio {lower_20_epitope/antigen_rmsd.shape[0]:.3f}")
lower_50_epitope = np.where(antigen_rmsd["epi-rmsd"]<5.0)[0].shape[0]
print(f"lower_50_epitope {lower_50_epitope} ratio {lower_50_epitope/antigen_rmsd.shape[0]:.3f}")
print(f"lower_50_epitope pdbs {np.where(antigen_rmsd['epi-rmsd']<5.0)}")
lower_20_full = np.where(antigen_rmsd["antigen-rmsd"]<2.0)[0].shape[0]
print(f"lower_20_full {lower_20_full} ratio {lower_20_full/antigen_rmsd.shape[0]:.3f}")

#data_antigen = [antigen_rmsd["rmsd"], antigen_rmsd["rmsd_epitope"]]
data_antigen = [antigen_rmsd["antigen-rmsd"], antigen_rmsd["epi-rmsd"]]
plots = ax[2].violinplot(data_antigen, showmedians=True, showextrema=False, widths=1)
print(f"antigen mean full {np.mean(antigen_rmsd['antigen-rmsd'])} +- {np.std(antigen_rmsd['antigen-rmsd'])}")
print(f"antigen mean epitope {np.mean(antigen_rmsd['epi-rmsd'])} +- {np.std(antigen_rmsd['epi-rmsd'])}")
# medians
print(f"antigen median full {np.median(antigen_rmsd['antigen-rmsd'])}")
print(f"antigen median epitope {np.median(antigen_rmsd['epi-rmsd'])}")
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
ax[0].set_xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], labels=XTICK_ANTIBODY + ["UNB"], fontsize=16)
ax[0].tick_params(axis='y', labelsize=16)
ax[1].set_xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], labels=XTICK_ANTIBODY+ ["UNB"], fontsize=16)
ax[1].tick_params(axis='y', labelsize=16)
ax[2].set_xticks([1, 2], labels=['Full', 'Epitope'], fontsize=16)
#ax[2].set_xticks([1], labels=['Epitope'], fontsize=16)
ax[2].tick_params(axis='y', labelsize=16)

ax[0].set_ylabel("CDRH3 RMSD ($\AA$)",fontsize=16)
ax[1].set_ylabel("Paratope RMSD ($\AA$)",fontsize=16)
ax[2].set_ylabel("Antigen RMSD ($\AA$)",fontsize=16)
ax[0].set_ylim((0.0,10.5))
ax[1].set_ylim((0.0,7.0))
ax[2].set_ylim((0.0,18.0))

plt.savefig(Path("figures", "violinplots", "violinplot_full_unb.png"), dpi=400, bbox_inches="tight")

#
# SAME BUT WITHOUT UNB DATA
#
fig, ax = plt.subplots(1, 3, figsize=(21, 9), width_ratios=[8.5, 8.5, 4])
plt.subplots_adjust(wspace=0.2, hspace=0.1)
plots = ax[0].violinplot(data_h3_unb[:-1], showmedians=True, showextrema=False, widths=1)

# Set the color of the violin patches
for pc, color in zip(plots['bodies'], colors):
    pc.set_facecolor(color)
plots['cmedians'].set_colors("black")

plots = ax[1].violinplot(data_para_unb[:-1], showmedians=True, showextrema=False, widths=1)
plots['cmedians'].set_colors("black")
# Set the color of the violin patches
for pc, color in zip(plots['bodies'], colors):
    pc.set_facecolor(color)

data_antigen = [antigen_rmsd["antigen-rmsd"], antigen_rmsd["epi-rmsd"]]
plots = ax[2].violinplot(data_antigen, showmedians=True, showextrema=False, widths=1)
print(f"antigen full mean {np.mean(antigen_rmsd['antigen-rmsd'])} +- {np.std(antigen_rmsd['antigen-rmsd'])}")
print(f"antigen mean epitope {np.mean(antigen_rmsd['epi-rmsd'])} +- {np.std(antigen_rmsd['epi-rmsd'])}")
# Set the color of the median lines
plots['cmedians'].set_colors("black")
for pc, color in zip(plots['bodies'], colors):
    pc.set_facecolor("Purple")
# write medians
medians_h3 = [np.median(d) for d in data_h3_unb[:-1]]
medians_para = [np.median(d) for d in data_para_unb[:-1]]
medians_antigen = [np.median(d) for d in data_antigen]

ratios_30 = [np.where(np.array(d)<3.0)[0].shape[0]/len(d) for d in data_para_unb[:-1]]
ratios_25 = [np.where(np.array(d)<2.5)[0].shape[0]/len(d) for d in data_para_unb[:-1]]
print(f"ratios {ratios_30}")
print(f"ratios {ratios_25}")

ax[1].axhline(3.0, color="black", linestyle="--")
ax[1].axhline(2.5, color="black", linestyle="--")

for i, median in enumerate(medians_h3):
    ax[0].text(i + 1, median+0.1, f'${median:.2f}\ \AA$', ha='center', va='bottom', color='black', fontsize=12)
for i, median in enumerate(medians_para):
    if i == 1:
        ax[1].text(i + 1, median-0.25, f'${median:.2f}\ \AA$', ha='center', va='bottom', color='black', fontsize=12)
    else:
        ax[1].text(i + 1, median+0.05, f'${median:.2f}\ \AA$', ha='center', va='bottom', color='black', fontsize=12)
for i, median in enumerate(medians_antigen):
    ax[2].text(i + 1, median+0.1, f'${median:.2f}\ \AA$', ha='center', va='bottom', color='black', fontsize=12)

# percentage of structures with RMSD < 3.0
for i, ratio in enumerate(ratios_30):
    ax[1].text(i+1, 3.05, f'{ratios_30[i]*100:.0f}%', ha='center', va='bottom', color='black', fontsize=12)
for i, ratio in enumerate(ratios_25):
    ax[1].text(i+1, 2.55, f'{ratios_25[i]*100:.0f}%', ha='center', va='bottom', color='black', fontsize=12)

# Set the labels
ax[0].set_xticks([1, 2, 3, 4, 5, 6, 7, 8, 9], labels=XTICK_ANTIBODY, fontsize=16)
ax[0].tick_params(axis='y', labelsize=16)
ax[1].set_xticks([1, 2, 3, 4, 5, 6, 7, 8, 9], labels=XTICK_ANTIBODY, fontsize=16)
ax[1].tick_params(axis='y', labelsize=16)
ax[2].set_xticks([1, 2], labels=['Full', 'Epitope'], fontsize=16)
#ax[2].set_xticks([1], labels=['Epitope'], fontsize=16)
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

