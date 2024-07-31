import pandas as pd
from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np
from functions import get_sorted_runs, load_clt_data, create_bound_gap_clt_dictionary, create_bound_bound_clt_dictionary
NPDBS = 82
plt.rcParams["font.family"] = "Helvetica"

cat_dict = {"Para-Epi": "Para-Epi",
            "CDR-EpiVag": "CDR-VagueEpi" ,
            "CDR-EpiVag-AA" : "CDR-VagueEpi-AA"}

# LOAD DATA
df_clt, df_clt_bound, df_clt_loose, df_clt_bound_loose = load_clt_data()
tot_runs = np.unique(df_clt_loose["pdb"]).shape[0] # should be 83
print(f"total number of runs {tot_runs}")
assert tot_runs == NPDBS

# CREATE DICTIONARIES
af2_bound_gap_clt = create_bound_gap_clt_dictionary(df_clt)
af2_bound_gap_clt_loose = create_bound_gap_clt_dictionary(df_clt_loose)
bound_bound_clt = create_bound_bound_clt_dictionary(df_clt_bound)
bound_bound_clt_loose = create_bound_bound_clt_dictionary(df_clt_bound_loose)
#print(af2_bound_gap_clt)
#print(af2_bound_gap_clt_loose)
print(F"plotting clustering data")
# PLOT
# standard df_clt based plot

# loose clustering SR definition
def stacked_clustering_barplot(bound_gap_clt, bound_bound_clt, fname):
    """
    Create a stacked barplot of the clustering results.

    Parameters
    ----------
    bound_gap_clt : dict
        dictionary of bound-gap clustering results
    fname : str
        filename of the plot
    
    Returns
    -------
    None
    """
    print(f"creating plot {fname}")
    fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(15, 6))
    plt.subplots_adjust(wspace=0.1, hspace=0.1)
    plt.suptitle("HADDOCK3 clustering performances", fontsize=24, verticalalignment="center")
    col = 0
    for cat in ["Para-Epi", "CDR-EpiVag-AA", "CDR-EpiVag"]:
        print("category", cat)
        #xticks = ["ABB", "ABBE", "ABL", "AF2", "IG", "ENS", "ENSNOAF2", "ENS196-48", "ENS196-CLT"]
        xticks = ["ABB", "ABL", "AF2", "IG", "ABBE", "AF2E", "IGE", "ENS", "ENSNOAF2", "CLE", "ENS196-48", "ENS196-CLT"]
        xticks.append("BOUND-BOUND")
        print(f"xticks {xticks}")
        sorted_runs = get_sorted_runs(cat)
        print(f"sorted_runs {sorted_runs}")

        mod_b_1 = [bound_gap_clt[cat][run][1][0] for run in sorted_runs]
        mod_af2_1 = [bound_gap_clt[cat][run][1][1] for run in sorted_runs]
        mod_b_frac_1 = [el/tot_runs for el in mod_b_1]
        mod_af2_frac_1 = [el/tot_runs for el in mod_af2_1]
        # top 2
        mod_b_2 = [bound_gap_clt[cat][run][2][0] for run in sorted_runs]
        mod_af2_2 = [bound_gap_clt[cat][run][2][1] for run in sorted_runs]
        mod_b_frac_2 = [el/tot_runs for el in mod_b_2]
        mod_af2_frac_2 = [el/tot_runs for el in mod_af2_2]
        # top 3
        mod_b_3 = [bound_gap_clt[cat][run][3][0] for run in sorted_runs]
        mod_af2_3 = [bound_gap_clt[cat][run][3][1] for run in sorted_runs]
        mod_b_frac_3 = [el/tot_runs for el in mod_b_3]
        mod_af2_frac_3 = [el/tot_runs for el in mod_af2_3]
        
        print(f"mod_b_frac_1 {[round(el, 3) for el in mod_b_frac_1]}")
        print(f"mod_af2_frac_1 {[round(el, 3) for el in mod_af2_frac_1]}")
        print(f"mod_b_frac_2 {[round(el, 3) for el in mod_b_frac_2]}")
        print(f"mod_af2_frac_2 {[round(el, 3) for el in mod_af2_frac_2]}")
        print(f"mod_b_frac_3 {[round(el, 3) for el in mod_b_frac_3]}")
        print(f"mod_af2_frac_3 {[round(el, 3) for el in mod_af2_frac_3]}")
        
        # parameters
        n=len(mod_b_frac_1)
        r = np.arange(n)
        width = 0.25

        #bars
        axs[col].bar(r, mod_b_frac_1, color = 'b',
                width = width, edgecolor = 'black',
                label='Bound T1')
        axs[col].bar(r + width, mod_af2_frac_1, color = 'orange',
                width = width, edgecolor = 'black',
                label='AF2 T1')
        # TOP 2 blurred
        axs[col].bar(r, mod_b_frac_2, color = 'b', alpha=0.6,
                width = width, edgecolor = 'black',
                label='Bound T2')
        axs[col].bar(r + width, mod_af2_frac_2, color = 'orange', alpha=0.6,
                width = width, edgecolor = 'black',
                label='AF2 T2')
        # TOP 3 more blurred
        axs[col].bar(r, mod_b_frac_3, color = 'b', alpha=0.3,
                width = width, edgecolor = 'black',
                label='Bound T3')
        axs[col].bar(r + width, mod_af2_frac_3, color = 'orange', alpha=0.3,
                width = width, edgecolor = 'black',
                label='AF2 T3')

        # adding bound bound results
        axs[col].bar(r[-1]+1,
                       bound_bound_clt[cat][1]/tot_runs,
                       color="r",
                       width = width,
                       label="T1 BB")
        axs[col].bar(r[-1]+1,
                       bound_bound_clt[cat][2]/tot_runs,
                       color="r",
                       alpha=0.6,
                       width = width,
                       label="T2 BB")
        axs[col].bar(r[-1]+1,
                       bound_bound_clt[cat][3]/tot_runs,
                       color="r",
                       alpha=0.2,
                       width = width,
                       label="T3 BB")

        axs[col].set_ylim((0,1.01))
        if col > 0:
            axs[col].set_yticks([])
        else:
            axs[col].set_ylabel("Clustering Success Rate", fontsize=16)

        new_r = np.arange(n)
        new_r = np.append(new_r, n - width/2)
        axs[col].set_xticks(new_r + width/2, xticks, fontsize=16, rotation = 90)
        col += 1

    axs[0].set_title(cat_dict["Para-Epi"], fontsize=16)
    axs[2].set_title(cat_dict["CDR-EpiVag"], fontsize=16)
    axs[1].set_title(cat_dict["CDR-EpiVag-AA"], fontsize=16)

    y_values = [0.2, 0.4, 0.6, 0.8]
    for y in y_values:
        axs[0].axhline(y=y, color='black', linestyle='-', alpha=0.1, zorder=0)
        axs[1].axhline(y=y, color='black', linestyle='-', alpha=0.1, zorder=0)
        axs[2].axhline(y=y, color='black', linestyle='-', alpha=0.1, zorder=0)

    axs[1].set_xlabel("Docking protocol", fontsize=24, labelpad=20)
    lines_labels = [axs[1].get_legend_handles_labels()]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    fig.legend(lines, labels, ncol=9, fontsize=15, loc="lower center", bbox_to_anchor = (0,-0.1, 1, 1))
    plt.tight_layout()

    plt.savefig(Path("figures", "clustering", fname), dpi=300, bbox_inches="tight")

stacked_clustering_barplot(af2_bound_gap_clt, bound_bound_clt, fname="clustering_1-3_fracbound.png")
print(f"\n\nLOOSE DEFINITION\n\n")
stacked_clustering_barplot(af2_bound_gap_clt_loose, bound_bound_clt_loose, fname="clustering_loose_1-3_fracbound.png")
