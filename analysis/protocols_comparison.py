import pandas as pd
from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np
from functions import load_data, create_bound_gap_dictionary, create_bound_bound_dictionary, get_sorted_runs
NPDBS = 82
plt.rcParams["font.family"] = "Helvetica"

cat_dict = {"Para-Epi": "Para-Epi",
            "CDR-EpiVag": "CDR-VagueEpi" ,
            "CDR-EpiVag-AA" : "CDR-VagueEpi-AA"}

# LOAD DATA
rigidbody_capri, rigidbody_capri_bound, emref_capri, emref_capri_bound, zdock_ss, emref_rigid_capri, af2multimer_ss = load_data()
#rigidbody_capri, rigidbody_capri_bound, emref_capri, emref_capri_bound, df_ss_flexref, df_ss_bound_flexref = load_data()
tot_runs = np.unique(rigidbody_capri["pdb"]).shape[0] # should be 79
print(f"total number of runs {tot_runs}")
assert tot_runs == NPDBS
# extract rigidbody, flexref and emref data
acc_key="acc"
af2_bound_gap_rigid = create_bound_gap_dictionary(rigidbody_capri, acc_key=acc_key)
af2_bound_gap_emref = create_bound_gap_dictionary(emref_capri, acc_key=acc_key)
bound_bound_rigid = create_bound_bound_dictionary(rigidbody_capri_bound, acc_key=acc_key)
bound_bound_emref = create_bound_bound_dictionary(emref_capri_bound, acc_key=acc_key)
# comparison with af2 multimer
af2_multimer_dict = {}
for a_key in ["acc", "med", "high"]:
    af2_multimer_dict[a_key] = {}
    for top in [1,10]:
        column = f"{a_key}_top{top}"
        # where this column is not 0
        succ_runs = af2multimer_ss.loc[af2multimer_ss[column] != 0]
        af2_multimer_dict[a_key][top] = succ_runs.shape[0]/tot_runs

print(f"af2multimer ACC {af2_multimer_dict}")

# PLOT
print("Creating stacked barplot for paper")
fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(15, 9))
plt.subplots_adjust(wspace=0.1, hspace=0.1)
plt.suptitle("HADDOCK3 performances", fontsize=24, verticalalignment="center")
col = 0

# rigidbody
for cat in ["Para-Epi", "CDR-EpiVag-AA", "CDR-EpiVag"]:
    print(f"\ncat {cat}\n")
    sorted_runs = get_sorted_runs(cat)
    # top 1
    mod_b_1 = [af2_bound_gap_rigid[cat][run][1][0] for run in sorted_runs]
    mod_af2_1 = [af2_bound_gap_rigid[cat][run][1][1] for run in sorted_runs]
    mod_b_frac_1 = [el/tot_runs for el in mod_b_1]
    mod_af2_frac_1 = [el/tot_runs for el in mod_af2_1]
    # top 10
    mod_b_10 = [af2_bound_gap_rigid[cat][run][10][0] for run in sorted_runs]
    mod_af2_10 = [af2_bound_gap_rigid[cat][run][10][1] for run in sorted_runs]
    mod_b_frac_10 = [el/tot_runs for el in mod_b_10]
    mod_af2_frac_10 = [el/tot_runs for el in mod_af2_10]
    print(f"mod_b_frac_1 {[round(el, 3) for el in mod_b_frac_1]}\n mod_af2_frac_1 {[round(el, 3) for el in mod_af2_frac_1]}")
    print(f"mod_b_frac_10 {[round(el, 3) for el in mod_b_frac_10]}\n mod_af2_frac_10 {[round(el, 3) for el in mod_af2_frac_10]}")
    # parameters
    n=len(mod_b_frac_1)
    r = np.arange(n)
    width = 0.25
    
    #bars
    axs[0][col].bar(r, mod_b_frac_1, color = 'b',
            width = width, edgecolor = 'black',
            label='Bound T1')
    axs[0][col].bar(r + width, mod_af2_frac_1, color = 'orange',
            width = width, edgecolor = 'black',
            label='AF2 T1')
    # TOP 10 blurred
    axs[0][col].bar(r, mod_b_frac_10, color = 'b', alpha=0.4,
            width = width, edgecolor = 'black',
            label='Bound T2')
    axs[0][col].bar(r + width, mod_af2_frac_10, color = 'orange', alpha=0.4,
            width = width, edgecolor = 'black',
            label='AF2 T10')
    
    # adding bound bound results
    axs[0][col].bar(r[-1]+1,
                    bound_bound_rigid[cat][1]/tot_runs,
                    color="r",
                    width = width,
                    label="T1 BB")
    axs[0][col].bar(r[-1]+1,
                    bound_bound_rigid[cat][10]/tot_runs,
                    color="r",
                    alpha=0.4,
                    width = width,
                    label="T10 BB")
    
    # adding Alphafold2 multimer results
    axs[0][col].bar(r[-1]+2,
                        af2_multimer_dict["acc"][1],
                        color="g",
                        width = width,
                        label="AF2 Multimer T1")
    axs[0][col].bar(r[-1]+2,
                        af2_multimer_dict["acc"][10],
                        color="g",
                        alpha=0.4,
                        width = width,
                        label="AF2 Multimer T10")
                        
    axs[0][col].set_ylim((0,1.01))
    if col > 0:
        axs[0][col].set_yticks([])
    else:
        axs[0][col].set_ylabel("Rigid-body Success Rate", fontsize=16)
        
    axs[0][col].set_xticks([])
    col += 1
    print(f"Adding bound_bound_rigid results for {cat} (acc_key = {acc_key}): t1 {bound_bound_rigid[cat][1]/tot_runs:.3f} t10 {bound_bound_rigid[cat][10]/tot_runs:.3f}")


print(f"\n\nnow emref data\n\n")

col = 0
for cat in ["Para-Epi", "CDR-EpiVag-AA", "CDR-EpiVag"]:
    print(f"\ncat {cat}\n")
    xticks = ["ABB", "ABL", "AF2", "IG", "ABBE", "AF2E", "IGE", "ENS", "ENSNOAF2", "CLE", "ENS196-48", "ENS196-CLT"]
    xticks.append("BOUND-BOUND")
    xticks.append("AF2-MULTIMER")
    print(f"xticks {xticks}")
    
    sorted_runs = get_sorted_runs(cat)

    mod_b_1 = [af2_bound_gap_emref[cat][run][1][0] for run in sorted_runs]
    mod_af2_1 = [af2_bound_gap_emref[cat][run][1][1] for run in sorted_runs]
    mod_b_frac_1 = [el/tot_runs for el in mod_b_1]
    mod_af2_frac_1 = [el/tot_runs for el in mod_af2_1]
    # top 10
    mod_b_10 = [af2_bound_gap_emref[cat][run][10][0] for run in sorted_runs]
    mod_af2_10 = [af2_bound_gap_emref[cat][run][10][1] for run in sorted_runs]
    mod_b_frac_10 = [el/tot_runs for el in mod_b_10]
    mod_af2_frac_10 = [el/tot_runs for el in mod_af2_10]

    print(f"mod_b_frac_1 {[round(el, 3) for el in mod_b_frac_1]} mod_af2_frac_1 {[round(el, 3) for el in mod_af2_frac_1]}")
    print(f"mod_b_frac_10 {[round(el, 3) for el in mod_b_frac_10]} mod_af2_frac_10 {[round(el, 3) for el in mod_af2_frac_10]}")

    axs[1][col].bar(r, mod_b_frac_1, color = 'b',
            width = width, edgecolor = 'black',
            label='Bound T1')
    axs[1][col].bar(r + width, mod_af2_frac_1, color = 'orange',
            width = width, edgecolor = 'black',
            label='AF2 T1')
    # TOP 10 blurred
    axs[1][col].bar(r, mod_b_frac_10, color = 'b', alpha=0.4,
            width = width, edgecolor = 'black',
            label='Bound T10')
    axs[1][col].bar(r + width, mod_af2_frac_10, color = 'orange', alpha=0.4,
            width = width, edgecolor = 'black',
            label='AF2 T10')
    if col > 0:
        axs[1][col].set_yticks([])
    else:
        axs[1][col].set_ylabel("Emref Success Rate", fontsize=16)
        
    # adding bound bound results
    axs[1][col].bar(r[-1]+1,
                    bound_bound_emref[cat][1]/tot_runs,
                    color="r",
                    width = width,
                    label="T1 BB")
    axs[1][col].bar(r[-1]+1,
                    bound_bound_emref[cat][10]/tot_runs,
                    color="r",
                    alpha=0.4,
                    width = width,
                    label="T10 BB")
    
    # adding Alphafold2 multimer results
    axs[1][col].bar(r[-1]+2,
                        af2_multimer_dict["acc"][1],
                        color="g",
                        width = width,
                        label="AF2 Multimer T1")
    axs[1][col].bar(r[-1]+2,
                        af2_multimer_dict["acc"][10],
                        color="g",
                        alpha=0.4,
                        width = width,
                        label="AF2 Multimer T10")
    ##
    axs[1][col].set_ylim((0,1.01))
    new_r = np.arange(n)
    new_r = np.append(new_r, n - width/2)
    new_r = np.append(new_r, n+1 - width/2)
    axs[1][col].set_xticks(new_r + width/2, xticks, fontsize=16, rotation = 90)
    print(f"Adding bound_bound_emref results for {cat} (acc_key = {acc_key}): t1 {bound_bound_emref[cat][1]/tot_runs:.3f} t10 {bound_bound_emref[cat][10]/tot_runs:.3f}")

    
    col += 1
axs[0][0].set_title(cat_dict["Para-Epi"], fontsize=16)
axs[0][2].set_title(cat_dict["CDR-EpiVag"], fontsize=16)
axs[0][1].set_title(cat_dict["CDR-EpiVag-AA"], fontsize=16)

y_values = [0.2, 0.4, 0.6, 0.8]
for y in y_values:
    axs[0][0].axhline(y=y, color='black', linestyle='-', alpha=0.1, zorder=0)
    axs[1][0].axhline(y=y, color='black', linestyle='-', alpha=0.1, zorder=0)
    axs[0][1].axhline(y=y, color='black', linestyle='-', alpha=0.1, zorder=0)
    axs[1][1].axhline(y=y, color='black', linestyle='-', alpha=0.1, zorder=0)
    axs[0][2].axhline(y=y, color='black', linestyle='-', alpha=0.1, zorder=0)
    axs[1][2].axhline(y=y, color='black', linestyle='-', alpha=0.1, zorder=0)

# inserting letters
axs[0][0].text(0.005, 0.93, "a)", transform=axs[0][0].transAxes, size=14, weight='bold')
axs[0][1].text(0.005, 0.93, "b)", transform=axs[0][1].transAxes, size=14, weight='bold')
axs[0][2].text(0.005, 0.93, "c)", transform=axs[0][2].transAxes, size=14, weight='bold')
axs[1][0].text(0.005, 0.93, "d)", transform=axs[1][0].transAxes, size=14, weight='bold')
axs[1][1].text(0.005, 0.93, "e)", transform=axs[1][1].transAxes, size=14, weight='bold')
axs[1][2].text(0.005, 0.93, "f)", transform=axs[1][2].transAxes, size=14, weight='bold')

axs[1][1].set_xlabel("Docking protocol", fontsize=24, labelpad=20)
lines_labels = [axs[1][1].get_legend_handles_labels()]
lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]

fig.legend(lines, labels, ncol=4, fontsize=16, loc="lower center", bbox_to_anchor = (0,-0.1, 1, 1))
plt.tight_layout()
plt.savefig(Path("figures", "af2_bound_comparison", f"rigid_emref_1-10_fracbound.png"), dpi=300, bbox_inches="tight")
plt.close()

# 2x2 plot for presentation
print("Creating stacked barplot for presentation")
fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(10, 9))
plt.subplots_adjust(wspace=0.1, hspace=0.1)
    
#plt.suptitle("HADDOCK3 performances", fontsize=24, verticalalignment="center")
col = 0
# rigidbody
for cat in ["Para-Epi", "CDR-EpiVag-AA"]:
    #sorted_runs = [f'run-af2ab-{cat}-mpi-50-50', f'run-af2abens-{cat}-mpi-50-50', f'run-af2abl-{cat}-mpi-50-50',
    #               f'run-af2af2-{cat}-mpi-50-50', f'run-af2ig-{cat}-mpi-50-50', f'run-af2ens-{cat}-mpi-50-50', 
    #               f'run-af2ensnoaf2-{cat}-mpi-50-50', f'run-af2ens-{cat}-mpi-196-48', f'run-af2ens-{cat}-mpi-196-clt']
    sorted_runs = get_sorted_runs(cat)

    # top 1
    mod_b_1 = [af2_bound_gap_rigid[cat][run][1][0] for run in sorted_runs]
    mod_af2_1 = [af2_bound_gap_rigid[cat][run][1][1] for run in sorted_runs]
    mod_b_frac_1 = [el/tot_runs for el in mod_b_1]
    mod_af2_frac_1 = [el/tot_runs for el in mod_af2_1]
    # top 10
    mod_b_10 = [af2_bound_gap_rigid[cat][run][10][0] for run in sorted_runs]
    mod_af2_10 = [af2_bound_gap_rigid[cat][run][10][1] for run in sorted_runs]
    mod_b_frac_10 = [el/tot_runs for el in mod_b_10]
    mod_af2_frac_10 = [el/tot_runs for el in mod_af2_10]
    
    # parameters
    n=len(mod_b_frac_1)
    r = np.arange(n)
    width = 0.25
    
    #bars
    axs[0][col].bar(r, mod_b_frac_1, color = 'b',
            width = width, edgecolor = 'black',
            label='Bound T1')
    axs[0][col].bar(r + width, mod_af2_frac_1, color = 'orange',
            width = width, edgecolor = 'black',
            label='AF2 T1')
    # TOP 10 blurred
    axs[0][col].bar(r, mod_b_frac_10, color = 'b', alpha=0.4,
            width = width, edgecolor = 'black',
            label='Bound T2')
    axs[0][col].bar(r + width, mod_af2_frac_10, color = 'orange', alpha=0.4,
            width = width, edgecolor = 'black',
            label='AF2 T10')
    
    # adding bound bound results
    print(f"Adding bound_bound_rigid results for {cat} (acc_key = {acc_key}): t1 {bound_bound_rigid[cat][1]/tot_runs} t10 {bound_bound_rigid[cat][10]/tot_runs}")
    axs[0][col].bar(r[-1]+1,
                    bound_bound_rigid[cat][1]/tot_runs,
                    color="r",
                    width = width,
                    label="T1 BB")
    axs[0][col].bar(r[-1]+1,
                    bound_bound_rigid[cat][10]/tot_runs,
                    color="r",
                    alpha=0.4,
                    width = width,
                    label="T10 BB")
    
    # adding Alphafold2 multimer results
    axs[0][col].bar(r[-1]+2,
                        af2_multimer_dict["acc"][1],
                        color="g",
                        width = width,
                        label="AF2 Multimer T1")
    axs[0][col].bar(r[-1]+2,
                        af2_multimer_dict["acc"][10],
                        color="g",
                        alpha=0.4,
                        width = width,
                        label="AF2 Multimer T10")
    
    axs[0][col].set_ylim((0,1.01))
    if col > 0:
        axs[0][col].set_yticks([])
    else:
        axs[0][col].set_ylabel("Rigid-body Success Rate", fontsize=16)
        
    axs[0][col].set_xticks([])
    col += 1

col = 0
for cat in ["Para-Epi", "CDR-EpiVag-AA"]:
    #xticks.append("BOUND-BOUND")
    print(f"xticks {xticks}")
    sorted_runs = get_sorted_runs(cat)
    # top 1    
    mod_b_1 = [af2_bound_gap_emref[cat][run][1][0] for run in sorted_runs]
    mod_af2_1 = [af2_bound_gap_emref[cat][run][1][1] for run in sorted_runs]
    mod_b_frac_1 = [el/tot_runs for el in mod_b_1]
    mod_af2_frac_1 = [el/tot_runs for el in mod_af2_1]
    # top 10
    mod_b_10 = [af2_bound_gap_emref[cat][run][10][0] for run in sorted_runs]
    mod_af2_10 = [af2_bound_gap_emref[cat][run][10][1] for run in sorted_runs]
    mod_b_frac_10 = [el/tot_runs for el in mod_b_10]
    mod_af2_frac_10 = [el/tot_runs for el in mod_af2_10]

    axs[1][col].bar(r, mod_b_frac_1, color = 'b',
            width = width, edgecolor = 'black',
            label='Bound T1')
    axs[1][col].bar(r + width, mod_af2_frac_1, color = 'orange',
            width = width, edgecolor = 'black',
            label='AF2 T1')
    # TOP 10 blurred
    axs[1][col].bar(r, mod_b_frac_10, color = 'b', alpha=0.4,
            width = width, edgecolor = 'black',
            label='Bound T10')
    axs[1][col].bar(r + width, mod_af2_frac_10, color = 'orange', alpha=0.4,
            width = width, edgecolor = 'black',
            label='AF2 T10')
    if col > 0:
        axs[1][col].set_yticks([])
    else:
        axs[1][col].set_ylabel("Emref Success Rate", fontsize=16)
    
    print(f"Adding bound_bound_emref results for {cat} (acc_key = {acc_key}): t1 {bound_bound_emref[cat][1]/tot_runs:.3f} t10 {bound_bound_emref[cat][10]/tot_runs:.3f}")
    # adding bound bound results
    axs[1][col].bar(r[-1]+1,
                    bound_bound_emref[cat][1]/tot_runs,
                    color="r",
                    width = width,
                    label="T1 BB")
    axs[1][col].bar(r[-1]+1,
                    bound_bound_emref[cat][10]/tot_runs,
                    color="r",
                    alpha=0.4,
                    width = width,
                    label="T10 BB")
    
    axs[1][col].bar(r[-1]+2,
                        af2_multimer_dict["acc"][1],
                        color="g",
                        width = width,
                        label="AF2 Multimer T1")
    axs[1][col].bar(r[-1]+2,
                        af2_multimer_dict["acc"][10],
                        color="g",
                        alpha=0.4,
                        width = width,
                        label="AF2 Multimer T10")
    ##
    axs[1][col].set_ylim((0,1.01))
    new_r = np.arange(n)
    new_r = np.append(new_r, n - width/2)
    new_r = np.append(new_r, n+1 - width/2)
    axs[1][col].set_xticks(new_r + width/2, xticks, fontsize=16, rotation = 90)
    
    col += 1
axs[0][0].set_title(cat_dict["Para-Epi"], fontsize=16)
axs[0][1].set_title(cat_dict["CDR-EpiVag"], fontsize=16) # people don't want to see the AA

fig.text(0.5, -0.02, "Docking protocol", fontsize=24, ha='center')
lines_labels = [axs[1][1].get_legend_handles_labels()]
lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]

y_values = [0.2, 0.4, 0.6, 0.8]
for y in y_values:
    axs[0][0].axhline(y=y, color='black', linestyle='-', alpha=0.1, zorder=0)
    axs[1][0].axhline(y=y, color='black', linestyle='-', alpha=0.1, zorder=0)
    axs[0][1].axhline(y=y, color='black', linestyle='-', alpha=0.1, zorder=0)
    axs[1][1].axhline(y=y, color='black', linestyle='-', alpha=0.1, zorder=0)
#plt.grid(color='black', linestyle='-', linewidth=20, axis="y", which="both")

fig.legend(lines, labels, ncol=4, fontsize=16, loc="lower center", bbox_to_anchor = (0,-0.125, 1, 1))
plt.tight_layout()
plt.savefig(Path("figures", "af2_bound_comparison", f"rigid_emref_1-10_fracbound_2cols.png"), dpi=300, bbox_inches="tight")
plt.close()

# medium quality poses
print("Creating stacked barplots for medium and high quality poses")
acc_key="med"
af2_bound_gap_rigid_med = create_bound_gap_dictionary(rigidbody_capri, acc_key=acc_key)
af2_bound_gap_emref_med = create_bound_gap_dictionary(emref_capri, acc_key=acc_key)
bound_bound_rigid_med = create_bound_bound_dictionary(rigidbody_capri_bound, acc_key=acc_key)
bound_bound_emref_med = create_bound_bound_dictionary(emref_capri_bound, acc_key=acc_key)

fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(10, 9))
plt.subplots_adjust(wspace=0.1, hspace=0.1)

col = 0
# rigidbody
for cat in ["Para-Epi", "CDR-EpiVag-AA"]:
    print(f"\ncat is {cat}\n")
    sorted_runs = get_sorted_runs(cat)

    # top 1
    mod_b_1 = [af2_bound_gap_rigid_med[cat][run][1][0] for run in sorted_runs]
    mod_af2_1 = [af2_bound_gap_rigid_med[cat][run][1][1] for run in sorted_runs]
    mod_b_frac_1 = [el/tot_runs for el in mod_b_1]
    mod_af2_frac_1 = [el/tot_runs for el in mod_af2_1]
    # top 10
    mod_b_10 = [af2_bound_gap_rigid_med[cat][run][10][0] for run in sorted_runs]
    mod_af2_10 = [af2_bound_gap_rigid_med[cat][run][10][1] for run in sorted_runs]
    mod_b_frac_10 = [el/tot_runs for el in mod_b_10]
    mod_af2_frac_10 = [el/tot_runs for el in mod_af2_10]

    print(f"Rigid_med_bound_top1 {[round(el, 3) for el in mod_b_frac_1]}")
    print(f"Rigid_med_af2_top1 {[round(el, 3) for el in mod_af2_frac_1]}")
    print(f"Rigid_med_bound_top10 {[round(el, 3) for el in mod_b_frac_10]}")
    print(f"Rigid_med_af2_top10 {[round(el, 3) for el in mod_af2_frac_10]}")
    
    # parameters
    n=len(mod_b_frac_1)
    r = np.arange(n)
    width = 0.25
    
    #bars
    axs[0][col].bar(r, mod_b_frac_1, color = 'b',
            width = width, edgecolor = 'black',
            label='Bound T1')
    axs[0][col].bar(r + width, mod_af2_frac_1, color = 'orange',
            width = width, edgecolor = 'black',
            label='AF2 T1')
    # TOP 10 blurred
    axs[0][col].bar(r, mod_b_frac_10, color = 'b', alpha=0.4,
            width = width, edgecolor = 'black',
            label='Bound T2')
    axs[0][col].bar(r + width, mod_af2_frac_10, color = 'orange', alpha=0.4,
            width = width, edgecolor = 'black',
            label='AF2 T10')
    
    # adding bound bound results
    print(f"Adding bound_bound_rigid results for {cat} (acc_key = {acc_key}):: t1 {bound_bound_rigid_med[cat][1]/tot_runs:.2f} t10 {bound_bound_rigid_med[cat][10]/tot_runs:.2f}")
    axs[0][col].bar(r[-1]+1,
                    bound_bound_rigid_med[cat][1]/tot_runs,
                    color="r",
                    width = width,
                    label="T1 BB")
    axs[0][col].bar(r[-1]+1,
                    bound_bound_rigid_med[cat][10]/tot_runs,
                    color="r",
                    alpha=0.4,
                    width = width,
                    label="T10 BB")
    
    # adding Alphafold2 multimer results
    axs[0][col].bar(r[-1]+2,
                        af2_multimer_dict["med"][1],
                        color="g",
                        width = width,
                        label="AF2 Multimer T1")
    axs[0][col].bar(r[-1]+2,
                        af2_multimer_dict["med"][10],
                        color="g",
                        alpha=0.4,
                        width = width,
                        label="AF2 Multimer T10")
    
    axs[0][col].set_ylim((0,1.01))
    if col > 0:
        axs[0][col].set_yticks([])
    else:
        axs[0][col].set_ylabel("Rigid-body Success Rate", fontsize=16)
        
    axs[0][col].set_xticks([])
    col += 1

col = 0
for cat in ["Para-Epi", "CDR-EpiVag-AA"]:
    print(f"\ncat is {cat}\n")
    xticks = ["ABB", "ABL", "AF2", "IG", "ABBE", "AF2E", "IGE", "ENS", "ENSNOAF2", "CLE", "ENS196-48", "ENS196-CLT"]
    xticks.append("BOUND-BOUND")
    xticks.append("AF2-MULTIMER")
    print(f"xticks {xticks}")
    sorted_runs = get_sorted_runs(cat)
    # top 1    
    mod_b_1 = [af2_bound_gap_emref_med[cat][run][1][0] for run in sorted_runs]
    mod_af2_1 = [af2_bound_gap_emref_med[cat][run][1][1] for run in sorted_runs]
    mod_b_frac_1 = [el/tot_runs for el in mod_b_1]
    mod_af2_frac_1 = [el/tot_runs for el in mod_af2_1]
    # top 10
    mod_b_10 = [af2_bound_gap_emref_med[cat][run][10][0] for run in sorted_runs]
    mod_af2_10 = [af2_bound_gap_emref_med[cat][run][10][1] for run in sorted_runs]
    mod_b_frac_10 = [el/tot_runs for el in mod_b_10]
    mod_af2_frac_10 = [el/tot_runs for el in mod_af2_10]

    print(f"Emref_med_bound_top1 {[round(el, 3) for el in mod_b_frac_1]}")
    print(f"Emref_med_af2_top1 {[round(el, 3) for el in mod_af2_frac_1]}")
    print(f"Emref_med_bound_top10 {[round(el, 3) for el in mod_b_frac_10]}")
    print(f"Emref_med_af2_top10 {[round(el, 3) for el in mod_af2_frac_10]}")
    

    axs[1][col].bar(r, mod_b_frac_1, color = 'b',
            width = width, edgecolor = 'black',
            label='Bound T1')
    axs[1][col].bar(r + width, mod_af2_frac_1, color = 'orange',
            width = width, edgecolor = 'black',
            label='AF2 T1')
    # TOP 10 blurred
    axs[1][col].bar(r, mod_b_frac_10, color = 'b', alpha=0.4,
            width = width, edgecolor = 'black',
            label='Bound T10')
    axs[1][col].bar(r + width, mod_af2_frac_10, color = 'orange', alpha=0.4,
            width = width, edgecolor = 'black',
            label='AF2 T10')
    if col > 0:
        axs[1][col].set_yticks([])
    else:
        axs[1][col].set_ylabel("Emref Success Rate", fontsize=16)
    
    print(f"Adding bound_bound_emref results for {cat} (acc_key = {acc_key}):: t1 {bound_bound_emref_med[cat][1]/tot_runs:.2f} t10 {bound_bound_emref_med[cat][10]/tot_runs:.2f}")
    # adding bound bound results
    axs[1][col].bar(r[-1]+1,
                    bound_bound_emref_med[cat][1]/tot_runs,
                    color="r",
                    width = width,
                    label="T1 BB")
    axs[1][col].bar(r[-1]+1,
                    bound_bound_emref_med[cat][10]/tot_runs,
                    color="r",
                    alpha=0.4,
                    width = width,
                    label="T10 BB")
    
    # adding Alphafold2 multimer results
    axs[1][col].bar(r[-1]+2,
                        af2_multimer_dict["med"][1],
                        color="g",
                        width = width,
                        label="AF2 Multimer T1")
    axs[1][col].bar(r[-1]+2,
                        af2_multimer_dict["med"][10],
                        color="g",
                        alpha=0.4,
                        width = width,
                        label="AF2 Multimer T10")
    ##
    axs[1][col].set_ylim((0,1.01))
    new_r = np.arange(n)
    new_r = np.append(new_r, n - width/2)
    new_r = np.append(new_r, n+1 - width/2)
    axs[1][col].set_xticks(new_r + width/2, xticks, fontsize=16, rotation = 90)
    
    col += 1
axs[0][0].set_title(cat_dict["Para-Epi"], fontsize=16)
axs[0][1].set_title(cat_dict["CDR-EpiVag-AA"], fontsize=16) # people don't want to see the AA

# inserting letters
axs[0][0].text(0.01, 0.93, "a)", transform=axs[0][0].transAxes, size=15, weight='bold')
axs[0][1].text(0.01, 0.93, "b)", transform=axs[0][1].transAxes, size=15, weight='bold')
axs[1][0].text(0.01, 0.93, "c)", transform=axs[1][0].transAxes, size=15, weight='bold')
axs[1][1].text(0.01, 0.93, "d)", transform=axs[1][1].transAxes, size=15, weight='bold')

fig.text(0.5, -0.02, "Docking protocol", fontsize=24, ha='center')
lines_labels = [axs[1][1].get_legend_handles_labels()]
lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]

y_values = [0.2, 0.4, 0.6, 0.8]
for y in y_values:
    axs[0][0].axhline(y=y, color='black', linestyle='-', alpha=0.1, zorder=0)
    axs[1][0].axhline(y=y, color='black', linestyle='-', alpha=0.1, zorder=0)
    axs[0][1].axhline(y=y, color='black', linestyle='-', alpha=0.1, zorder=0)
    axs[1][1].axhline(y=y, color='black', linestyle='-', alpha=0.1, zorder=0)
#plt.grid(color='black', linestyle='-', linewidth=20, axis="y", which="both")

fig.legend(lines, labels, ncol=4, fontsize=16, loc="lower center", bbox_to_anchor = (0,-0.125, 1, 1))
plt.tight_layout()
plt.savefig(Path("figures", "af2_bound_comparison", f"rigid_emref_1-10_fracbound_2cols_med.png"), dpi=300, bbox_inches="tight")
plt.close()

# high quality poses
acc_key = "high"
af2_bound_gap_rigid_high = create_bound_gap_dictionary(rigidbody_capri, acc_key=acc_key)
af2_bound_gap_emref_high = create_bound_gap_dictionary(emref_capri, acc_key=acc_key)
bound_bound_rigid_high = create_bound_bound_dictionary(rigidbody_capri_bound, acc_key=acc_key)
bound_bound_emref_high = create_bound_bound_dictionary(emref_capri_bound, acc_key=acc_key)
fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(10, 9))
plt.subplots_adjust(wspace=0.1, hspace=0.1)

col = 0
# rigidbody
for cat in ["Para-Epi", "CDR-EpiVag-AA"]:
    print(f"\ncat is {cat}\n")
    sorted_runs = get_sorted_runs(cat)

    # top 1
    mod_b_1 = [af2_bound_gap_rigid_high[cat][run][1][0] for run in sorted_runs]
    mod_af2_1 = [af2_bound_gap_rigid_high[cat][run][1][1] for run in sorted_runs]
    mod_b_frac_1 = [el/tot_runs for el in mod_b_1]
    mod_af2_frac_1 = [el/tot_runs for el in mod_af2_1]
    # top 10
    mod_b_10 = [af2_bound_gap_rigid_high[cat][run][10][0] for run in sorted_runs]
    mod_af2_10 = [af2_bound_gap_rigid_high[cat][run][10][1] for run in sorted_runs]
    mod_b_frac_10 = [el/tot_runs for el in mod_b_10]
    mod_af2_frac_10 = [el/tot_runs for el in mod_af2_10]
    # top 48
    mod_b_48 = [af2_bound_gap_rigid_high[cat][run][48][0] for run in sorted_runs]
    mod_af2_48 = [af2_bound_gap_rigid_high[cat][run][48][1] for run in sorted_runs]
    mod_b_frac_48 = [el/tot_runs for el in mod_b_48]
    mod_af2_frac_48 = [el/tot_runs for el in mod_af2_48]
    
    print(f"mod_b_frac_10: {[round(el, 3) for el in mod_b_frac_10]}")
    # parameters
    n=len(mod_b_frac_1)
    r = np.arange(n)
    width = 0.25
    
    #bars
    axs[0][col].bar(r, mod_b_frac_1, color = 'b',
            width = width, edgecolor = 'black',
            label='Bound T1')
    axs[0][col].bar(r + width, mod_af2_frac_1, color = 'orange',
            width = width, edgecolor = 'black',
            label='AF2 T1')
    # TOP 10 blurred
    axs[0][col].bar(r, mod_b_frac_10, color = 'b', alpha=0.4,
            width = width, edgecolor = 'black',
            label='Bound T2')
    axs[0][col].bar(r + width, mod_af2_frac_10, color = 'orange', alpha=0.4,
            width = width, edgecolor = 'black',
            label='AF2 T10')
        # TOP 48 blurred
    axs[0][col].bar(r, mod_b_frac_48, color = 'b', alpha=0.2,
            width = width, edgecolor = 'black',
            label='Bound T48')
    axs[0][col].bar(r + width, mod_af2_frac_48, color = 'orange', alpha=0.2,
            width = width, edgecolor = 'black',
            label='AF2 T48')
    
    # adding bound bound results
    print(f"Adding bound_bound_rigid results for {cat} (acc_key = {acc_key}): t1 {bound_bound_rigid_high[cat][1]/tot_runs} t10 {bound_bound_rigid_high[cat][10]/tot_runs}")
    axs[0][col].bar(r[-1]+1,
                    bound_bound_rigid_high[cat][1]/tot_runs,
                    color="r",
                    width = width,
                    label="T1 BB")
    axs[0][col].bar(r[-1]+1,
                    bound_bound_rigid_high[cat][10]/tot_runs,
                    color="r",
                    alpha=0.4,
                    width = width,
                    label="T10 BB")
    axs[0][col].bar(r[-1]+1,
                    bound_bound_rigid_high[cat][48]/tot_runs,
                    color="r",
                    alpha=0.2,
                    width = width,
                    label="T48 BB")
    # adding Alphafold2 multimer results
    axs[0][col].bar(r[-1]+2,
                af2_multimer_dict["high"][1],
                color="g",
                width = width,
                label="AF2 Multimer T1")
    axs[0][col].bar(r[-1]+2,
                af2_multimer_dict["high"][10],
                color="g",
                alpha=0.4,
                width = width,
                label="AF2 Multimer T10")
    
    axs[0][col].set_ylim((0,1.01))
    if col > 0:
        axs[0][col].set_yticks([])
    else:
        axs[0][col].set_ylabel("Rigid-body Success Rate", fontsize=16)
        
    axs[0][col].set_xticks([])
    col += 1

col = 0
for cat in ["Para-Epi", "CDR-EpiVag-AA"]:
    print(f"\ncat is {cat}\n")
    xticks = ["ABB", "ABL", "AF2", "IG", "ABBE", "AF2E", "IGE", "ENS", "ENSNOAF2", "CLE", "ENS196-48", "ENS196-CLT"]
    xticks.append("BOUND-BOUND")
    xticks.append("AF2-MULTIMER")
    print(f"xticks {xticks}")
    sorted_runs = get_sorted_runs(cat)
    # top 1    
    mod_b_1 = [af2_bound_gap_emref_high[cat][run][1][0] for run in sorted_runs]
    mod_af2_1 = [af2_bound_gap_emref_high[cat][run][1][1] for run in sorted_runs]
    mod_b_frac_1 = [el/tot_runs for el in mod_b_1]
    mod_af2_frac_1 = [el/tot_runs for el in mod_af2_1]
    # top 10
    mod_b_10 = [af2_bound_gap_emref_high[cat][run][10][0] for run in sorted_runs]
    mod_af2_10 = [af2_bound_gap_emref_high[cat][run][10][1] for run in sorted_runs]
    mod_b_frac_10 = [el/tot_runs for el in mod_b_10]
    mod_af2_frac_10 = [el/tot_runs for el in mod_af2_10]
    # top 48
    mod_b_48 = [af2_bound_gap_emref_high[cat][run][48][0] for run in sorted_runs]
    mod_af2_48 = [af2_bound_gap_emref_high[cat][run][48][1] for run in sorted_runs]
    mod_b_frac_48 = [el/tot_runs for el in mod_b_48]
    mod_af2_frac_48 = [el/tot_runs for el in mod_af2_48]

    print(f"Emref_high_bound_top1 {[round(el, 3) for el in mod_b_frac_1]}")
    print(f"Emref_high_af2_top1 {[round(el, 3) for el in mod_af2_frac_1]}")
    print(f"Emref_high_bound_top10 {[round(el, 3) for el in mod_b_frac_10]}")
    print(f"Emref_high_af2_top10 {[round(el, 3) for el in mod_af2_frac_10]}")
          
    axs[1][col].bar(r, mod_b_frac_1, color = 'b',
            width = width, edgecolor = 'black',
            label='Bound T1')
    axs[1][col].bar(r + width, mod_af2_frac_1, color = 'orange',
            width = width, edgecolor = 'black',
            label='AF2 T1')
    # TOP 10 blurred
    axs[1][col].bar(r, mod_b_frac_10, color = 'b', alpha=0.4,
            width = width, edgecolor = 'black',
            label='Bound T10')
    axs[1][col].bar(r + width, mod_af2_frac_10, color = 'orange', alpha=0.4,
            width = width, edgecolor = 'black',
            label='AF2 T10')
    # TOP 48 more blurred
    axs[1][col].bar(r, mod_b_frac_48, color = 'b', alpha=0.2,
        width = width, edgecolor = 'black',
        label='Bound T48')
    axs[1][col].bar(r + width, mod_af2_frac_48, color = 'orange', alpha=0.2,
        width = width, edgecolor = 'black',
        label='AF2 T48')

    if col > 0:
        axs[1][col].set_yticks([])
    else:
        axs[1][col].set_ylabel("Emref Success Rate", fontsize=16)
    
    print(f"Adding bound_bound_emref results for {cat} (acc_key = {acc_key}): t1 {bound_bound_emref_high[cat][1]/tot_runs:.2f} t10 {bound_bound_emref_high[cat][10]/tot_runs:.2f}")
    # adding bound bound results
    axs[1][col].bar(r[-1]+1,
                    bound_bound_emref_high[cat][1]/tot_runs,
                    color="r",
                    width = width,
                    label="T1 BB")
    axs[1][col].bar(r[-1]+1,
                    bound_bound_emref_high[cat][10]/tot_runs,
                    color="r",
                    alpha=0.4,
                    width = width,
                    label="T10 BB")
    axs[1][col].bar(r[-1]+1,
                    bound_bound_emref_high[cat][48]/tot_runs,
                    color="r",
                    alpha=0.2,
                    width = width,
                    label="T48 BB")
    
    # adding Alphafold2 multimer results
    axs[1][col].bar(r[-1]+2,
                af2_multimer_dict["high"][1],
                color="g",
                width = width,
                label="AF2 Multimer T1")
    axs[1][col].bar(r[-1]+2,
                af2_multimer_dict["high"][10],
                color="g",
                alpha=0.4,
                width = width,
                label="AF2 Multimer T10")
    ##
    axs[1][col].set_ylim((0,1.01))
    new_r = np.arange(n)
    new_r = np.append(new_r, n - width/2)
    new_r = np.append(new_r, n+1 - width/2)
    axs[1][col].set_xticks(new_r + width/2, xticks, fontsize=16, rotation = 90)
    
    col += 1
axs[0][0].set_title(cat_dict["Para-Epi"], fontsize=16)
axs[0][1].set_title(cat_dict["CDR-EpiVag-AA"], fontsize=16) # people don't want to see the AA

# inserting letters
axs[0][0].text(0.01, 0.93, "a)", transform=axs[0][0].transAxes, size=15, weight='bold')
axs[0][1].text(0.01, 0.93, "b)", transform=axs[0][1].transAxes, size=15, weight='bold')
axs[1][0].text(0.01, 0.93, "c)", transform=axs[1][0].transAxes, size=15, weight='bold')
axs[1][1].text(0.01, 0.93, "d)", transform=axs[1][1].transAxes, size=15, weight='bold')

fig.text(0.5, -0.02, "Docking protocol", fontsize=24, ha='center')
lines_labels = [axs[1][1].get_legend_handles_labels()]
lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]

y_values = [0.2, 0.4, 0.6, 0.8]
for y in y_values:
    axs[0][0].axhline(y=y, color='black', linestyle='-', alpha=0.1, zorder=0)
    axs[1][0].axhline(y=y, color='black', linestyle='-', alpha=0.1, zorder=0)
    axs[0][1].axhline(y=y, color='black', linestyle='-', alpha=0.1, zorder=0)
    axs[1][1].axhline(y=y, color='black', linestyle='-', alpha=0.1, zorder=0)
#plt.grid(color='black', linestyle='-', linewidth=20, axis="y", which="both")

fig.legend(lines, labels, ncol=4, fontsize=16, loc="lower center", bbox_to_anchor = (0,-0.165, 1, 1))
plt.tight_layout()
plt.savefig(Path("figures", "af2_bound_comparison", f"rigid_emref_1-10_fracbound_2cols_high.png"), dpi=300, bbox_inches="tight")
plt.close()
