from functions import get_sorted_runs, load_data, create_bound_gap_dictionary
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
NPDBS = 82
cat_dict = {"Para-Epi": "Para-Epi",
            "CDR-EpiVag": "CDR-VagueEpi" ,
            "CDR-EpiVag-AA" : "CDR-VagueEpi-AA"}

rigidbody_capri, rigidbody_capri_bound, emref_capri, emref_capri_bound, df_ss_flexref, df_ss_bound_flexref, zdock_ss, emref_rigid_capri, af2multimer_ss = load_data()

af2_bound_gap_rigid = create_bound_gap_dictionary(rigidbody_capri, acc_key="acc")
af2_bound_gap_rigid_emref = create_bound_gap_dictionary(emref_rigid_capri, acc_key="acc")

print(f"af2_bound_gap_rigid_emref {af2_bound_gap_rigid_emref}")

tot_runs = np.unique(emref_rigid_capri["pdb"]).shape[0] # should be 79
assert tot_runs == NPDBS

print("CREATING STACKED BARPLOT FOR PAPER")
plt.rcParams["font.family"] = "Helvetica"

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(15, 8))
plt.subplots_adjust(wspace=0.1, hspace=0.1, bottom=0.25, top=0.925)
plt.suptitle("HADDOCK3-rigid vs emref", fontsize=24, verticalalignment="center")
col = 0
colors = plt.cm.tab20b.colors

width = 0.2
# rigidbody
for cat in ["Para-Epi", "CDR-EpiVag-AA"]:
    #xticks = ["ABB", "ABBE", "ABL", "AF2", "IG", "ENS", "ENSNOAF2", "ENS196-48", "ENS196-CLT"]
    xticks = ["ABB", "ABL", "AF2", "IG", "ABBE", "AF2E", "IGE", "ENS", "ENSNOAF2", "CLE", "ENS196-48", "ENS196-CLT"]
    sorted_runs = get_sorted_runs(cat)
    
    # haddock data
    mod_b_haddock_1 = [af2_bound_gap_rigid[cat][run][1][0] for run in sorted_runs]
    mod_af2_haddock_1 = [af2_bound_gap_rigid[cat][run][1][1] for run in sorted_runs]
    mod_b_haddock_10 = [af2_bound_gap_rigid[cat][run][10][0] for run in sorted_runs]
    mod_af2_haddock_10 = [af2_bound_gap_rigid[cat][run][10][1] for run in sorted_runs]
    
    # fractions
    mod_b_haddock_frac_1 = [el/tot_runs for el in mod_b_haddock_1]
    mod_af2_haddock_frac_1 = [el/tot_runs for el in mod_af2_haddock_1]
    mod_b_haddock_frac_10 = [el/tot_runs for el in mod_b_haddock_10]
    mod_af2_haddock_frac_10 = [el/tot_runs for el in mod_af2_haddock_10]
    #print(mod_b_haddock_frac_1)
    
    mod_b_rigid_emref_1 = [af2_bound_gap_rigid_emref[cat][run][1][0] for run in sorted_runs]
    mod_af2_rigid_emref_1 = [af2_bound_gap_rigid_emref[cat][run][1][1] for run in sorted_runs]
    mod_b_rigid_emref_10 = [af2_bound_gap_rigid_emref[cat][run][10][0] for run in sorted_runs]
    mod_af2_rigid_emref_10 = [af2_bound_gap_rigid_emref[cat][run][10][1] for run in sorted_runs]
    
    #rigid_emref fractions
    mod_b_rigid_emref_frac_1 = [el/tot_runs for el in mod_b_rigid_emref_1]
    mod_af2_rigid_emref_frac_1 = [el/tot_runs for el in mod_af2_rigid_emref_1]
    mod_b_rigid_emref_frac_10 = [el/tot_runs for el in mod_b_rigid_emref_10]
    mod_af2_rigid_emref_frac_10 = [el/tot_runs for el in mod_af2_rigid_emref_10]

    n=len(mod_b_haddock_frac_1)
    r = np.arange(n)

    axs[col].bar(r, mod_b_haddock_frac_1, color = colors[2],
                width = width, edgecolor = 'black',
                label='Rigid Bound')
    axs[col].bar(r + width, mod_b_rigid_emref_frac_1, color = colors[0],
                width = width, edgecolor = 'black',
                label='Emref Bound')
    axs[col].bar(r + width*2, mod_af2_haddock_frac_1, color = colors[6],
                width = width, edgecolor = 'black',
                label='Rigid AF2')
    axs[col].bar(r + width*3, mod_af2_rigid_emref_frac_1, color = colors[4],
                width = width, edgecolor = 'black',
                label='Emref AF2')
    # TOP 10 blurred
    axs[col].bar(r, mod_b_haddock_frac_10, color = colors[2], alpha=0.4,
                width = width, edgecolor = 'black',
                label='Rigid Bound T10')
    axs[col].bar(r + width, mod_b_rigid_emref_frac_10, color = colors[0], alpha=0.4,
                width = width, edgecolor = 'black',
                label='Emref BoundT10')
    axs[col].bar(r + width*2, mod_af2_haddock_frac_10, color = colors[6], alpha=0.4,
                width = width, edgecolor = 'black',
                label='Rigid AF2 T10')
    axs[col].bar(r + width*3, mod_af2_rigid_emref_frac_10, color = colors[4], alpha=0.4,
                width = width, edgecolor = 'black',
                label='Emref AF2 T10')

    # print comparison of performance
    print(f"Rigid Bound {cat} 1: {np.mean(mod_b_haddock_frac_1)}")
    print(f"Rigid+Emref Bound {cat} 1: {np.mean(mod_b_rigid_emref_frac_1)}")
    print(f"Rigid AF2 {cat} 1: {np.mean(mod_af2_haddock_frac_1)}")
    print(f"Rigid+Emref AF2 {cat} 1: {np.mean(mod_af2_rigid_emref_frac_1)}")
    print(f"Rigid Bound {cat} 10: {np.mean(mod_b_haddock_frac_10)}")
    print(f"Rigid+Emref Bound {cat} 10: {np.mean(mod_b_rigid_emref_frac_10)}")
    print(f"Rigid AF2 {cat} 10: {np.mean(mod_af2_haddock_frac_10)}")
    print(f"Rigid+Emref AF2 {cat} 10: {np.mean(mod_af2_rigid_emref_frac_10)}")
    # deltas between emref and rigid
    print(f"Delta Bound {cat} 1: {np.mean(mod_b_rigid_emref_frac_1) - np.mean(mod_b_haddock_frac_1)}")
    print(f"Delta AF2 {cat} 1: {np.mean(mod_af2_rigid_emref_frac_1) - np.mean(mod_af2_haddock_frac_1)}")
    print(f"Delta Bound {cat} 10: {np.mean(mod_b_rigid_emref_frac_10) - np.mean(mod_b_haddock_frac_10)}")
    print(f"Delta AF2 {cat} 10: {np.mean(mod_af2_rigid_emref_frac_10) - np.mean(mod_af2_haddock_frac_10)}")

    # print overall deltas
    print(f"Delta Bound {cat} 1: {np.mean(mod_b_rigid_emref_frac_1):.2f} - {np.mean(mod_b_haddock_frac_1):.2f} = {(np.mean(mod_b_rigid_emref_frac_1) - np.mean(mod_b_haddock_frac_1)):.2f}")
    print(f"Delta Bound {cat} 10: {np.mean(mod_b_rigid_emref_frac_10):.2f} - {np.mean(mod_b_haddock_frac_10):.2f} = {(np.mean(mod_b_rigid_emref_frac_10) - np.mean(mod_b_haddock_frac_10)):.2f}")

    # now loop over deltas to extract improvements
    for k in range(len(xticks)):
        delta_bound_1 = mod_b_rigid_emref_frac_1[k] - mod_b_haddock_frac_1[k]
        delta_af2_1 = mod_af2_rigid_emref_frac_1[k] - mod_af2_haddock_frac_1[k]
        delta_bound_10 = mod_b_rigid_emref_frac_10[k] - mod_b_haddock_frac_10[k]
        delta_af2_10 = mod_af2_rigid_emref_frac_10[k] - mod_af2_haddock_frac_10[k]
        print(f"{xticks[k]} Delta Bound {cat} 1 {xticks[k]}: {delta_bound_1:.2f}")
        print(f"{xticks[k]} Delta AF2 {cat} 1 {xticks[k]}: {delta_af2_1:.2f}")
        print(f"{xticks[k]} Delta Bound {cat} 10 {xticks[k]}: {delta_bound_10:.2f}")
        print(f"{xticks[k]} Delta AF2 {cat} 10 {xticks[k]}: {delta_af2_10:.2f}")
    axs[col].set_ylim((0,1.01))
    if col > 0:
        axs[col].set_yticks([])
    else:
        axs[col].set_ylabel("Docking Success Rate", fontsize=16)
    #
    #axs[col].set_xticks([])
    axs[col].set_ylim((0,1.01))
    new_r = np.arange(n)
    #new_r = np.append(new_r, n - 3/2*width)
    axs[col].set_xticks(new_r + 3/2*width, xticks, fontsize=16, rotation = 90)
    col += 1

axs[0].set_title(cat_dict["Para-Epi"], fontsize=16)
#axs[2].set_title(cat_dict["CDR-EpiVag"], fontsize=16)
axs[1].set_title(cat_dict["CDR-EpiVag-AA"], fontsize=16)

y_values = [0.2, 0.4, 0.6, 0.8]
for y in y_values:
    axs[0].axhline(y=y, color='black', linestyle='-', alpha=0.1, zorder=0)
    axs[1].axhline(y=y, color='black', linestyle='-', alpha=0.1, zorder=0)

fig.text(0.5, -0.02, "Docking protocol", fontsize=24, ha='center')
#axs[1][1].set_xlabel("Docking protocol", fontsize=24, labelpad=20)
lines_labels = [axs[1].get_legend_handles_labels()]
lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]

fig.legend(lines, labels, ncol=4, fontsize=16, loc="lower center", bbox_to_anchor = (0,-0.15, 1, 1))
plt.tight_layout()
plt.savefig(Path("figures", "af2_bound_comparison", "rigid_vs_rigid-emref_1-10_fracbound.png"), dpi=300, bbox_inches="tight")