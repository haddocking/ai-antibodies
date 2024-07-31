import pandas as pd
from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np
from functions import load_data, create_bound_gap_dictionary, create_bound_bound_dictionary, create_zdock_dictionary
plt.rcParams["font.family"] = "Helvetica"
NPDBS = 82
cat_dict = {"Para-Epi": "Para-Epi",
            "CDR-EpiVag": "CDR-VagueEpi" ,
            "CDR-EpiVag-AA" : "CDR-VagueEpi-AA"}

# LOAD DATA
rigidbody_capri, rigidbody_capri_bound, emref_capri, emref_capri_bound, df_ss_flexref, df_ss_bound_flexref, zdock_ss, emref_rigid_capri, af2multimer_ss = load_data()
tot_runs = np.unique(zdock_ss["pdb"]).shape[0] # should be 83
print(f"total number of runs {tot_runs}")
assert tot_runs == NPDBS
# extract rigidbody, flexref and emref data
af2_bound_gap_rigid = create_bound_gap_dictionary(rigidbody_capri)
af2_bound_gap_emref = create_bound_gap_dictionary(emref_capri)
bound_bound_rigid = create_bound_bound_dictionary(rigidbody_capri_bound)
bound_bound_emref = create_bound_bound_dictionary(emref_capri_bound)
zdock_data = create_zdock_dictionary(zdock_ss)

# PLOT
print("CREATING RIGIDBODY PLOTS")
tops = [1, 5, 10, 20, 48]

for cat in ["Para-Epi", "CDR-EpiVag"]:
    haddock_rb_af2_runs = [f"run-af2ab-{cat}-mpi-50-50", f"run-af2abl-{cat}-mpi-50-50",
                             f"run-af2af2-{cat}-mpi-50-50", f"run-af2ig-{cat}-mpi-50-50"]
    haddock_rb_af2_aa_runs = [f"run-af2ab-{f'{cat}-AA'}-mpi-50-50", f"run-af2abl-{f'{cat}-AA'}-mpi-50-50",
                             f"run-af2af2-{f'{cat}-AA'}-mpi-50-50", f"run-af2ig-{f'{cat}-AA'}-mpi-50-50"]
    for top in tops:
        plt.figure(figsize=(15, 9), dpi=250)
        plt.ylim((0,81))
        if cat == "CDR-EpiVag":
            mod_b_haddock_aa = [af2_bound_gap_rigid[f"{cat}-AA"][run][top][0] for run in haddock_rb_af2_aa_runs]
            mod_af2_haddock_aa = [af2_bound_gap_rigid[f"{cat}-AA"][run][top][1] for run in haddock_rb_af2_aa_runs]
            mod_b_haddock_aa_frac = [el/tot_runs for el in mod_b_haddock_aa]
            mod_af2_haddock_aa_frac = [el/tot_runs for el in mod_af2_haddock_aa]
            
        mod_b_haddock = [af2_bound_gap_rigid[cat][run][top][0] for run in haddock_rb_af2_runs]
        mod_af2_haddock = [af2_bound_gap_rigid[cat][run][top][1] for run in haddock_rb_af2_runs]
        mod_b_zdock = [zdock_data[cat][run][top][0] for run in zdock_data[cat]]
        mod_af2_zdock = [zdock_data[cat][run][top][1] for run in zdock_data[cat]]
        # fractions
        mod_b_haddock_frac = [el/tot_runs for el in mod_b_haddock]
        mod_af2_haddock_frac = [el/tot_runs for el in mod_af2_haddock]
        mod_b_zdock_frac = [el/tot_runs for el in mod_b_zdock]
        mod_af2_zdock_frac = [el/tot_runs for el in mod_af2_zdock]
        
        xticks = []
        for run in zdock_data[cat]:
            run_name = run.split('-')[0]
            xticks.append(run_name)
        
        n=len(mod_b_haddock)
        r = np.arange(n)
        width = 0.2
          
        plt.bar(r, mod_b_haddock, color = 'b',
                width = width, edgecolor = 'black',
                label='HADDOCK Bound')
        
        plt.bar(r + width, mod_af2_haddock, color = 'g',
                width = width, edgecolor = 'black',
                label='HADDOCK AF2')
        
        plt.bar(r + width*2, mod_b_zdock, color = 'lightblue',
                width = width, edgecolor = 'black',
                label='ZDOCK Bound')
        
        plt.bar(r + width*3, mod_af2_zdock, color = 'lightgreen',
                width = width, edgecolor = 'black',
                label='ZDOCK AF2')
        
        plt.axhline(y = bound_bound_rigid[cat][top], color = 'r', linestyle = '-', label="HADDOCK bound-bound")
        
        plt.xlabel("Run", fontsize = 20)
        plt.ylabel(f"Acceptable pdbs (top{top})", fontsize = 20)
        
        plt.title(f"ZDOCK vs HADDOCK (rigidbody) {cat_dict[cat]} scenario)",fontsize=24)
        plt.xticks(r + 1.5*width, xticks, fontsize=16)
        plt.yticks(fontsize=16)
        plt.legend(loc=(1.01, 0.75), ncol=1, borderaxespad=0., fontsize=20)
        plt.tight_layout()
        
        plt.savefig(Path("figures", "zdock_comparison", f"zdock_rigidbody_{cat}_top{top}.png"), dpi=300)
        
        plt.close()
        #plt.show()
        
        # fraction plot
        plt.figure(figsize=(15, 9))
        plt.ylim((0,1.01))
        plt.bar(r, mod_b_haddock_frac, color = 'b',
                width = width, edgecolor = 'black',
                label='HADDOCK Bound')
        plt.bar(r + width, mod_af2_haddock_frac, color = 'g',
                width = width, edgecolor = 'black',
                label='HADDOCK AF2')
        plt.bar(r + width*2, mod_b_zdock_frac, color = 'lightblue',
                width = width, edgecolor = 'black',
                label='ZDOCK Bound')
        plt.bar(r + width*3, mod_af2_zdock_frac, color = 'lightgreen',
                width = width, edgecolor = 'black',
                label='ZDOCK AF2')
        
        plt.axhline(y = bound_bound_rigid[cat][top]/tot_runs, color = 'r', linestyle = '-', label="HADDOCK bound-bound")
        plt.xlabel("Docking protocol", fontsize = 20)
        plt.ylabel(f"Fraction of acceptable pdbs (top{top})", fontsize = 20)
        plt.title(f"ZDOCK vs HADDOCK (rigidbody) {cat_dict[cat]} scenario)",fontsize=24)
        #
        plt.xticks(r + 1.5*width, xticks, fontsize=16)
        plt.yticks(fontsize=16)
        plt.legend(loc=(1.01, 0.75), ncol=1, borderaxespad=0., fontsize= 20)
        plt.tight_layout()
        plt.savefig(Path("figures", "zdock_comparison", f"zdock_rigidbody_{cat}_top{top}_frac.png"), dpi=300)
        plt.close()
        
        # CDR-EpiVag-AA scenario
        if cat == "CDR-EpiVag":
            plt.figure(figsize=(15, 9), dpi=250)
            plt.ylim((0,81))
            plt.bar(r, mod_b_haddock_aa, color = 'b',
                width = width, edgecolor = 'black',
                label='HADDOCK Bound')
        
            plt.bar(r + width, mod_af2_haddock_aa, color = 'g',
                width = width, edgecolor = 'black',
                label='HADDOCK AF2')
        
            plt.bar(r + width*2, mod_b_zdock, color = 'lightblue',
                width = width, edgecolor = 'black',
                label='ZDOCK Bound')
        
            plt.bar(r + width*3, mod_af2_zdock, color = 'lightgreen',
                width = width, edgecolor = 'black',
                label='ZDOCK AF2')
            plt.axhline(y = bound_bound_rigid[f"{cat}-AA"][top], color = 'r', linestyle = '-', label="HADDOCK bound-bound")
            
            plt.xlabel("Run", fontsize = 20)
            plt.ylabel(f"Acceptable pdbs (top{top})", fontsize = 20)
            plt.title(f"ZDOCK vs HADDOCK (rigidbody) {cat_dict[f'{cat}-AA']} scenario)",fontsize=24)
            #
            plt.xticks(r + 1.5*width, xticks, fontsize=16)
            plt.yticks(fontsize=16)
            plt.legend(loc=(1.01, 0.85), ncol=1, borderaxespad=0., fontsize=20)
            plt.tight_layout()
            plt.savefig(Path("figures", "zdock_comparison", f"zdock_rigidbody_{f'{cat}-AA'}_top{top}.png"), dpi=300)
            plt.close()
            #plt.show()
            # fraction plot
            plt.figure(figsize=(15, 9))
            plt.ylim((0,1.01))
            plt.bar(r, mod_b_haddock_aa_frac, color = 'b',
                    width = width, edgecolor = 'black',
                    label='HADDOCK Bound')
            plt.bar(r + width, mod_af2_haddock_aa_frac, color = 'g',
                    width = width, edgecolor = 'black',
                    label='HADDOCK AF2')
            plt.bar(r + width*2, mod_b_zdock_frac, color = 'lightblue',
                    width = width, edgecolor = 'black',
                    label='ZDOCK Bound')
            plt.bar(r + width*3, mod_af2_zdock_frac, color = 'lightgreen',
                    width = width, edgecolor = 'black',
                    label='ZDOCK AF2')
            
            plt.axhline(y = bound_bound_rigid[f'{cat}-AA'][top]/tot_runs, color = 'r', linestyle = '-', label="HADDOCK bound-bound")
            plt.xlabel("Docking protocol", fontsize = 20)
            plt.ylabel(f"Fraction of acceptable pdbs (top{top})", fontsize = 20)
            plt.title(f"ZDOCK vs HADDOCK (rigidbody) {cat_dict[f'{cat}-AA']} scenario)",fontsize=24)
            #
            plt.xticks(r + 1.5*width, xticks, fontsize=16)
            plt.yticks(fontsize=16)
            plt.legend(loc=(1.01, 0.75), ncol=1, borderaxespad=0., fontsize= 20)
            plt.tight_layout()
            plt.savefig(Path("figures", "zdock_comparison", f"zdock_rigidbody_{f'{cat}-AA'}_top{top}_frac.png"), dpi=300)
            
            #plt.show()
            plt.close()

# now emref data
print("CREATING EMREF PLOTS")
for cat in ["Para-Epi", "CDR-EpiVag"]:
    haddock_rb_af2_runs = [f"run-af2ab-{cat}-mpi-50-50", f"run-af2abl-{cat}-mpi-50-50",
                             f"run-af2af2-{cat}-mpi-50-50", f"run-af2ig-{cat}-mpi-50-50"]
    haddock_rb_af2_aa_runs = [f"run-af2ab-{f'{cat}-AA'}-mpi-50-50", f"run-af2abl-{f'{cat}-AA'}-mpi-50-50",
                             f"run-af2af2-{f'{cat}-AA'}-mpi-50-50", f"run-af2ig-{f'{cat}-AA'}-mpi-50-50"]
    for top in tops:
        print(top, cat)
        plt.figure(figsize=(15, 9), dpi=250)
        plt.ylim((0,81))
        if cat == "CDR-EpiVag":
            mod_b_haddock_aa = [af2_bound_gap_emref[f"{cat}-AA"][run][top][0] for run in haddock_rb_af2_aa_runs]
            mod_af2_haddock_aa = [af2_bound_gap_emref[f"{cat}-AA"][run][top][1] for run in haddock_rb_af2_aa_runs]
            mod_b_haddock_aa_frac = [el/tot_runs for el in mod_b_haddock_aa]
            mod_af2_haddock_aa_frac = [el/tot_runs for el in mod_af2_haddock_aa]
            
        mod_b_haddock = [af2_bound_gap_emref[cat][run][top][0] for run in haddock_rb_af2_runs]
        mod_af2_haddock = [af2_bound_gap_emref[cat][run][top][1] for run in haddock_rb_af2_runs]
        mod_b_zdock = [zdock_data[cat][run][top][0] for run in zdock_data[cat]]
        mod_af2_zdock = [zdock_data[cat][run][top][1] for run in zdock_data[cat]]
        # fractions
        mod_b_haddock_frac = [el/tot_runs for el in mod_b_haddock]
        mod_af2_haddock_frac = [el/tot_runs for el in mod_af2_haddock]
        mod_b_zdock_frac = [el/tot_runs for el in mod_b_zdock]
        mod_af2_zdock_frac = [el/tot_runs for el in mod_af2_zdock]
        
        xticks = []
        for run in zdock_data[cat]:
            run_name = run.split('-')[0]
            xticks.append(run_name)
        
        n=len(mod_b_haddock)
        r = np.arange(n)
        width = 0.2
          
        plt.bar(r, mod_b_haddock, color = 'b',
                width = width, edgecolor = 'black',
                label='HADDOCK Bound')
        
        plt.bar(r + width, mod_af2_haddock, color = 'g',
                width = width, edgecolor = 'black',
                label='HADDOCK AF2')
        
        plt.bar(r + width*2, mod_b_zdock, color = 'lightblue',
                width = width, edgecolor = 'black',
                label='ZDOCK Bound')
        
        plt.bar(r + width*3, mod_af2_zdock, color = 'lightgreen',
                width = width, edgecolor = 'black',
                label='ZDOCK AF2')
        
        plt.axhline(y = bound_bound_emref[cat][top], color = 'r', linestyle = '-', label="HADDOCK bound-bound")
        
        plt.xlabel("Run", fontsize = 20)
        plt.ylabel(f"Acceptable pdbs (top{top})", fontsize = 20)
        
        plt.title(f"ZDOCK vs HADDOCK (emref) {cat_dict[cat]} scenario)",fontsize=24)
        plt.xticks(r + 1.5*width, xticks, fontsize=16)
        plt.yticks(fontsize=16)
        plt.legend(loc=(1.01, 0.75), ncol=1, borderaxespad=0., fontsize=20)
        plt.tight_layout()
        
        plt.savefig(Path("figures", "zdock_comparison", f"zdock_emref_{cat}_top{top}.png"), dpi=300)
        plt.close()
        #plt.show()
        # fraction plot
        plt.figure(figsize=(15, 9))
        plt.ylim((0,1.01))
        plt.bar(r, mod_b_haddock_frac, color = 'b',
                width = width, edgecolor = 'black',
                label='HADDOCK Bound')
        plt.bar(r + width, mod_af2_haddock_frac, color = 'g',
                width = width, edgecolor = 'black',
                label='HADDOCK AF2')
        plt.bar(r + width*2, mod_b_zdock_frac, color = 'lightblue',
                width = width, edgecolor = 'black',
                label='ZDOCK Bound')
        plt.bar(r + width*3, mod_af2_zdock_frac, color = 'lightgreen',
                width = width, edgecolor = 'black',
                label='ZDOCK AF2')
        
        plt.axhline(y = bound_bound_emref[cat][top]/tot_runs, color = 'r', linestyle = '-', label="HADDOCK bound-bound")
        plt.xlabel("Docking protocol", fontsize = 20)
        plt.ylabel(f"Fraction of acceptable pdbs (top{top})", fontsize = 20)
        plt.title(f"ZDOCK vs HADDOCK (emref) {cat_dict[cat]} scenario)",fontsize=24)
        #
        plt.xticks(r + 1.5*width, xticks, fontsize=16)
        plt.yticks(fontsize=16)
        plt.legend(loc=(1.01, 0.75), ncol=1, borderaxespad=0., fontsize= 20)
        plt.tight_layout()
        plt.savefig(Path("figures", "zdock_comparison", f"zdock_emref_{cat}_top{top}_frac.png"), dpi=300)
        
        #plt.show()
        plt.close()
        # CDR-EpiVag-AA scenario
        if cat == "CDR-EpiVag":
            plt.figure(figsize=(15, 9), dpi=250)
            plt.ylim((0,81))
            plt.bar(r, mod_b_haddock_aa, color = 'b',
                width = width, edgecolor = 'black',
                label='HADDOCK Bound')
        
            plt.bar(r + width, mod_af2_haddock_aa, color = 'g',
                width = width, edgecolor = 'black',
                label='HADDOCK AF2')
        
            plt.bar(r + width*2, mod_b_zdock, color = 'lightblue',
                width = width, edgecolor = 'black',
                label='ZDOCK Bound')
        
            plt.bar(r + width*3, mod_af2_zdock, color = 'lightgreen',
                width = width, edgecolor = 'black',
                label='ZDOCK AF2')
            plt.axhline(y = bound_bound_emref[f"{cat}-AA"][top], color = 'r', linestyle = '-', label="HADDOCK bound-bound")
            
            plt.xlabel("Run", fontsize = 20)
            plt.ylabel(f"Acceptable pdbs (top{top})", fontsize = 20)
            plt.title(f"ZDOCK vs HADDOCK (emref) {cat_dict[f'{cat}-AA']} scenario)",fontsize=24)
            #
            plt.xticks(r + 1.5*width, xticks, fontsize=16)
            plt.yticks(fontsize=16)
            plt.legend(loc=(1.01, 0.85), ncol=1, borderaxespad=0., fontsize=20)
            plt.tight_layout()
            plt.savefig(Path("figures", "zdock_comparison", f"zdock_emref_{f'{cat}-AA'}_top{top}.png"), dpi=300)
            plt.close()
            #plt.show()
            # fraction plot
            plt.figure(figsize=(15, 9))
            plt.ylim((0,1.01))
            plt.bar(r, mod_b_haddock_aa_frac, color = 'b',
                    width = width, edgecolor = 'black',
                    label='HADDOCK Bound')
            plt.bar(r + width, mod_af2_haddock_aa_frac, color = 'g',
                    width = width, edgecolor = 'black',
                    label='HADDOCK AF2')
            plt.bar(r + width*2, mod_b_zdock_frac, color = 'lightblue',
                    width = width, edgecolor = 'black',
                    label='ZDOCK Bound')
            plt.bar(r + width*3, mod_af2_zdock_frac, color = 'lightgreen',
                    width = width, edgecolor = 'black',
                    label='ZDOCK AF2')
            
            plt.axhline(y = bound_bound_emref[f'{cat}-AA'][top]/tot_runs, color = 'r', linestyle = '-', label="HADDOCK bound-bound")
            plt.xlabel("Docking protocol", fontsize = 20)
            plt.ylabel(f"Fraction of acceptable pdbs (top{top})", fontsize = 20)
            plt.title(f"ZDOCK vs HADDOCK (emref) {cat_dict[f'{cat}-AA']} scenario)",fontsize=24)
            #
            plt.xticks(r + 1.5*width, xticks, fontsize=16)
            plt.yticks(fontsize=16)
            plt.legend(loc=(1.01, 0.75), ncol=1, borderaxespad=0., fontsize= 20)
            plt.tight_layout()
            plt.savefig(Path("figures", "zdock_comparison", f"zdock_emref_{f'{cat}-AA'}_top{top}_frac.png"), dpi=300)
            
            #plt.show()
            plt.close()

print("CREATING STACKED BARPLOT FOR PAPER")
plt.rcParams["font.family"] = "Helvetica"

fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(10, 12))
plt.subplots_adjust(wspace=0.1, hspace=0.1, bottom=0.25, top=0.925)
plt.suptitle("HADDOCK3-ZDOCK comparison", fontsize=24, verticalalignment="center")
col = 0
colors = plt.cm.Paired.colors
print(colors)

# rigidbody
for cat in ["Para-Epi", "CDR-EpiVag-AA"]:
    # haddock runs
    haddock_rb_af2_runs = [f"run-af2ab-{cat}-mpi-50-50", f"run-af2abl-{cat}-mpi-50-50",
                             f"run-af2af2-{cat}-mpi-50-50", f"run-af2ig-{cat}-mpi-50-50"]
    haddock_rb_af2_aa_runs = [f"run-af2ab-{f'{cat}-AA'}-mpi-50-50", f"run-af2abl-{f'{cat}-AA'}-mpi-50-50",
                             f"run-af2af2-{f'{cat}-AA'}-mpi-50-50", f"run-af2ig-{f'{cat}-AA'}-mpi-50-50"]
    
    # haddock data
    mod_b_haddock_1 = [af2_bound_gap_rigid[cat][run][1][0] for run in haddock_rb_af2_runs]
    mod_af2_haddock_1 = [af2_bound_gap_rigid[cat][run][1][1] for run in haddock_rb_af2_runs]
    mod_b_haddock_10 = [af2_bound_gap_rigid[cat][run][10][0] for run in haddock_rb_af2_runs]
    mod_af2_haddock_10 = [af2_bound_gap_rigid[cat][run][10][1] for run in haddock_rb_af2_runs]
    
    # fractions
    mod_b_haddock_frac_1 = [el/tot_runs for el in mod_b_haddock_1]
    mod_af2_haddock_frac_1 = [el/tot_runs for el in mod_af2_haddock_1]
    mod_b_haddock_frac_10 = [el/tot_runs for el in mod_b_haddock_10]
    mod_af2_haddock_frac_10 = [el/tot_runs for el in mod_af2_haddock_10]
    #print(mod_b_haddock_frac_1)
    # zdock data
    if cat == "CDR-EpiVag-AA":
        zdock_cat = "CDR-EpiVag"
    else:
        zdock_cat = cat
    mod_b_zdock_1 = [zdock_data[zdock_cat][run][1][0] for run in zdock_data[zdock_cat]]
    mod_af2_zdock_1 = [zdock_data[zdock_cat][run][1][1] for run in zdock_data[zdock_cat]]
    mod_b_zdock_10 = [zdock_data[zdock_cat][run][10][0] for run in zdock_data[zdock_cat]]
    mod_af2_zdock_10 = [zdock_data[zdock_cat][run][10][1] for run in zdock_data[zdock_cat]]
    
    #zdock fractions
    mod_b_zdock_frac_1 = [el/tot_runs for el in mod_b_zdock_1]
    mod_af2_zdock_frac_1 = [el/tot_runs for el in mod_af2_zdock_1]
    mod_b_zdock_frac_10 = [el/tot_runs for el in mod_b_zdock_10]
    mod_af2_zdock_frac_10 = [el/tot_runs for el in mod_af2_zdock_10]

    axs[0][col].bar(r, mod_b_haddock_frac_1, color = colors[1],
                width = width, edgecolor = 'black',
                label='HADDOCK Bound')
    axs[0][col].bar(r + width, mod_af2_haddock_frac_1, color = colors[3],
                width = width, edgecolor = 'black',
                label='HADDOCK AF2')
    axs[0][col].bar(r + width*2, mod_b_zdock_frac_1, color = colors[7],
                width = width, edgecolor = 'black',
                label='ZDOCK Bound')
    axs[0][col].bar(r + width*3, mod_af2_zdock_frac_1, color = colors[9],
                width = width, edgecolor = 'black',
                label='ZDOCK AF2')
    # TOP 10 blurred
    axs[0][col].bar(r, mod_b_haddock_frac_10, color = colors[1], alpha=0.4,
                width = width, edgecolor = 'black',
                label='HADDOCK Bound T10')
    axs[0][col].bar(r + width, mod_af2_haddock_frac_10, color = colors[3], alpha=0.4,
                width = width, edgecolor = 'black',
                label='HADDOCK AF2 T10')
    axs[0][col].bar(r + width*2, mod_b_zdock_frac_10, color = colors[7], alpha=0.4,
                width = width, edgecolor = 'black',
                label='ZDOCK BoundT10')
    axs[0][col].bar(r + width*3, mod_af2_zdock_frac_10, color = colors[9], alpha=0.4,
                width = width, edgecolor = 'black',
                label='ZDOCK AF2 T10')
    axs[0][col].set_ylim((0,1.01))
    if col > 0:
        axs[0][col].set_yticks([])
    else:
        axs[0][col].set_ylabel("Rigid-body Success Rate", fontsize=16)
        
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
    
    axs[0][col].set_xticks([])
    col += 1

col = 0
for cat in ["Para-Epi", "CDR-EpiVag-AA"]:
    # haddock runs
    haddock_rb_af2_runs = [f"run-af2ab-{cat}-mpi-50-50", f"run-af2abl-{cat}-mpi-50-50",
                             f"run-af2af2-{cat}-mpi-50-50", f"run-af2ig-{cat}-mpi-50-50"]
    haddock_rb_af2_aa_runs = [f"run-af2ab-{f'{cat}-AA'}-mpi-50-50", f"run-af2abl-{f'{cat}-AA'}-mpi-50-50",
                             f"run-af2af2-{f'{cat}-AA'}-mpi-50-50", f"run-af2ig-{f'{cat}-AA'}-mpi-50-50"]
    # xticks
    xticks_lower = []
    for run in haddock_rb_af2_runs:        
        run_name = run.split('-')[1][3:]
        if run_name == "ab":
            run_name = "abb"
        xticks_lower.append(run_name)
    xticks = [xt.upper() for xt in xticks_lower]
    xticks.append("BB")
    #print(f"xticks {xticks}")
    
    # haddock data
    mod_b_haddock_1_emref = [af2_bound_gap_emref[cat][run][1][0] for run in haddock_rb_af2_runs]
    mod_af2_haddock_1_emref = [af2_bound_gap_emref[cat][run][1][1] for run in haddock_rb_af2_runs]
    mod_b_haddock_10_emref = [af2_bound_gap_emref[cat][run][10][0] for run in haddock_rb_af2_runs]
    mod_af2_haddock_10_emref = [af2_bound_gap_emref[cat][run][10][1] for run in haddock_rb_af2_runs]
    
    # fractions
    mod_b_haddock_frac_1_emref = [el/tot_runs for el in mod_b_haddock_1_emref]
    mod_af2_haddock_frac_1_emref = [el/tot_runs for el in mod_af2_haddock_1_emref]
    mod_b_haddock_frac_10_emref = [el/tot_runs for el in mod_b_haddock_10_emref]
    mod_af2_haddock_frac_10_emref = [el/tot_runs for el in mod_af2_haddock_10_emref]
    
    # zdock data
    if cat == "CDR-EpiVag-AA":
        zdock_cat = "CDR-EpiVag"
    else:
        zdock_cat = cat
    mod_b_zdock_1 = [zdock_data[zdock_cat][run][1][0] for run in zdock_data[zdock_cat]]
    mod_af2_zdock_1 = [zdock_data[zdock_cat][run][1][1] for run in zdock_data[zdock_cat]]
    mod_b_zdock_10 = [zdock_data[zdock_cat][run][10][0] for run in zdock_data[zdock_cat]]
    mod_af2_zdock_10 = [zdock_data[zdock_cat][run][10][1] for run in zdock_data[zdock_cat]]
    
    #zdock fractions
    mod_b_zdock_frac_1 = [el/tot_runs for el in mod_b_zdock_1]
    mod_af2_zdock_frac_1 = [el/tot_runs for el in mod_af2_zdock_1]
    mod_b_zdock_frac_10 = [el/tot_runs for el in mod_b_zdock_10]
    mod_af2_zdock_frac_10 = [el/tot_runs for el in mod_af2_zdock_10]
    
    axs[1][col].bar(r, mod_b_haddock_frac_1_emref, color = colors[1],
                width = width, edgecolor = 'black',
                label='HADDOCK Bound T1')
    axs[1][col].bar(r + width, mod_af2_haddock_frac_1_emref, color = colors[3],
                width = width, edgecolor = 'black',
                label='HADDOCK AF2 T1')
    axs[1][col].bar(r + width*2, mod_b_zdock_frac_1, color = colors[7],
                width = width, edgecolor = 'black',
                label='ZDOCK Bound T1')
    axs[1][col].bar(r + width*3, mod_af2_zdock_frac_1, color = colors[9],
                width = width, edgecolor = 'black',
                label='ZDOCK AF2 T1')
    # TOP 10 blurred
    axs[1][col].bar(r, mod_b_haddock_frac_10_emref, color = colors[1], alpha=0.4,
                width = width, edgecolor = 'black',
                label='HADDOCK Bound T10')
    axs[1][col].bar(r + width, mod_af2_haddock_frac_10_emref, color = colors[3], alpha=0.4,
                width = width, edgecolor = 'black',
                label='HADDOCK AF2 T10')
    axs[1][col].bar(r + width*2, mod_b_zdock_frac_10, color = colors[7], alpha=0.4,
                width = width, edgecolor = 'black',
                label='ZDOCK Bound T10')
    axs[1][col].bar(r + width*3, mod_af2_zdock_frac_10, color = colors[9], alpha=0.4,
                width = width, edgecolor = 'black',
                label='ZDOCK AF2 T10')

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
    
    axs[1][col].set_ylim((0,1.01))
    new_r = np.arange(n)
    new_r = np.append(new_r, n - 3/2*width)
    axs[1][col].set_xticks(new_r + 3/2*width, xticks, fontsize=16, rotation = 90)
    col+=1

axs[0][0].set_title(cat_dict["Para-Epi"], fontsize=16)
#axs[0][2].set_title(cat_dict["CDR-EpiVag"], fontsize=16)
axs[0][1].set_title(cat_dict["CDR-EpiVag-AA"], fontsize=16)
fig.text(0.5, 0.17, "Docking protocol", fontsize=24, ha='center')
#axs[1][1].set_xlabel("Docking protocol", fontsize=24, labelpad=20)
lines_labels = [axs[1][1].get_legend_handles_labels()]
lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]

fig.legend(lines, labels, ncol=2, fontsize=16, loc="lower center")

plt.savefig(Path("figures", "zdock_comparison", f"zdock_rigid_emref_1-10_fracbound.png"), dpi=300)

# here's the 3-figure plot
# fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(15, 12))
# plt.subplots_adjust(wspace=0.1, hspace=0.1, bottom=0.25, top=0.925)
# plt.suptitle("HADDOCK3-ZDOCK comparison", fontsize=24, verticalalignment="center")
# col = 0
# colors = plt.cm.Paired.colors
# print(colors)

# # rigidbody
# for cat in ["Para-Epi", "CDR-EpiVag-AA", "CDR-EpiVag"]:
#     # haddock runs
#     haddock_rb_af2_runs = [f"run-af2ab-{cat}-mpi-50-50", f"run-af2abl-{cat}-mpi-50-50",
#                              f"run-af2af2-{cat}-mpi-50-50", f"run-af2ig-{cat}-mpi-50-50"]
#     haddock_rb_af2_aa_runs = [f"run-af2ab-{f'{cat}-AA'}-mpi-50-50", f"run-af2abl-{f'{cat}-AA'}-mpi-50-50",
#                              f"run-af2af2-{f'{cat}-AA'}-mpi-50-50", f"run-af2ig-{f'{cat}-AA'}-mpi-50-50"]
    
#     # haddock data
#     mod_b_haddock_1 = [af2_bound_gap_rigid[cat][run][1][0] for run in haddock_rb_af2_runs]
#     mod_af2_haddock_1 = [af2_bound_gap_rigid[cat][run][1][1] for run in haddock_rb_af2_runs]
#     mod_b_haddock_10 = [af2_bound_gap_rigid[cat][run][10][0] for run in haddock_rb_af2_runs]
#     mod_af2_haddock_10 = [af2_bound_gap_rigid[cat][run][10][1] for run in haddock_rb_af2_runs]
    
#     # fractions
#     mod_b_haddock_frac_1 = [el/tot_runs for el in mod_b_haddock_1]
#     mod_af2_haddock_frac_1 = [el/tot_runs for el in mod_af2_haddock_1]
#     mod_b_haddock_frac_10 = [el/tot_runs for el in mod_b_haddock_10]
#     mod_af2_haddock_frac_10 = [el/tot_runs for el in mod_af2_haddock_10]
#     #print(mod_b_haddock_frac_1)
#     # zdock data
#     if cat == "CDR-EpiVag-AA":
#         zdock_cat = "CDR-EpiVag"
#     else:
#         zdock_cat = cat
#     mod_b_zdock_1 = [zdock_data[zdock_cat][run][1][0] for run in zdock_data[zdock_cat]]
#     mod_af2_zdock_1 = [zdock_data[zdock_cat][run][1][1] for run in zdock_data[zdock_cat]]
#     mod_b_zdock_10 = [zdock_data[zdock_cat][run][10][0] for run in zdock_data[zdock_cat]]
#     mod_af2_zdock_10 = [zdock_data[zdock_cat][run][10][1] for run in zdock_data[zdock_cat]]
    
#     #zdock fractions
#     mod_b_zdock_frac_1 = [el/tot_runs for el in mod_b_zdock_1]
#     mod_af2_zdock_frac_1 = [el/tot_runs for el in mod_af2_zdock_1]
#     mod_b_zdock_frac_10 = [el/tot_runs for el in mod_b_zdock_10]
#     mod_af2_zdock_frac_10 = [el/tot_runs for el in mod_af2_zdock_10]

#     axs[0][col].bar(r, mod_b_haddock_frac_1, color = colors[1],
#                 width = width, edgecolor = 'black',
#                 label='HADDOCK Bound')
#     axs[0][col].bar(r + width, mod_af2_haddock_frac_1, color = colors[3],
#                 width = width, edgecolor = 'black',
#                 label='HADDOCK AF2')
#     axs[0][col].bar(r + width*2, mod_b_zdock_frac_1, color = colors[7],
#                 width = width, edgecolor = 'black',
#                 label='ZDOCK Bound')
#     axs[0][col].bar(r + width*3, mod_af2_zdock_frac_1, color = colors[9],
#                 width = width, edgecolor = 'black',
#                 label='ZDOCK AF2')
#     # TOP 10 blurred
#     axs[0][col].bar(r, mod_b_haddock_frac_10, color = colors[1], alpha=0.4,
#                 width = width, edgecolor = 'black',
#                 label='HADDOCK Bound T10')
#     axs[0][col].bar(r + width, mod_af2_haddock_frac_10, color = colors[3], alpha=0.4,
#                 width = width, edgecolor = 'black',
#                 label='HADDOCK AF2 T10')
#     axs[0][col].bar(r + width*2, mod_b_zdock_frac_10, color = colors[7], alpha=0.4,
#                 width = width, edgecolor = 'black',
#                 label='ZDOCK BoundT10')
#     axs[0][col].bar(r + width*3, mod_af2_zdock_frac_10, color = colors[9], alpha=0.4,
#                 width = width, edgecolor = 'black',
#                 label='ZDOCK AF2 T10')
#     axs[0][col].set_ylim((0,1.01))
#     if col > 0:
#         axs[0][col].set_yticks([])
#     else:
#         axs[0][col].set_ylabel("Rigid-body Success Rate", fontsize=16)
        
#     # adding bound bound results
#     axs[0][col].bar(r[-1]+1,
#                     bound_bound_rigid[cat][1]/tot_runs,
#                     color="r",
#                     width = width,
#                     label="T1 BB")
#     axs[0][col].bar(r[-1]+1,
#                     bound_bound_rigid[cat][10]/tot_runs,
#                     color="r",
#                     alpha=0.4,
#                     width = width,
#                     label="T10 BB")
    
#     axs[0][col].set_xticks([])
#     col += 1

# col = 0
# for cat in ["Para-Epi", "CDR-EpiVag-AA", "CDR-EpiVag"]:
#     # haddock runs
#     haddock_rb_af2_runs = [f"run-af2ab-{cat}-mpi-50-50", f"run-af2abl-{cat}-mpi-50-50",
#                              f"run-af2af2-{cat}-mpi-50-50", f"run-af2ig-{cat}-mpi-50-50"]
#     haddock_rb_af2_aa_runs = [f"run-af2ab-{f'{cat}-AA'}-mpi-50-50", f"run-af2abl-{f'{cat}-AA'}-mpi-50-50",
#                              f"run-af2af2-{f'{cat}-AA'}-mpi-50-50", f"run-af2ig-{f'{cat}-AA'}-mpi-50-50"]
#     # xticks
#     xticks_lower = []
#     for run in haddock_rb_af2_runs:        
#         run_name = run.split('-')[1][3:]
#         if run_name == "ab":
#             run_name = "abb"
#         xticks_lower.append(run_name)
#     xticks = [xt.upper() for xt in xticks_lower]
#     xticks.append("BB")
#     #print(f"xticks {xticks}")
    
#     # haddock data
#     mod_b_haddock_1_emref = [af2_bound_gap_emref[cat][run][1][0] for run in haddock_rb_af2_runs]
#     mod_af2_haddock_1_emref = [af2_bound_gap_emref[cat][run][1][1] for run in haddock_rb_af2_runs]
#     mod_b_haddock_10_emref = [af2_bound_gap_emref[cat][run][10][0] for run in haddock_rb_af2_runs]
#     mod_af2_haddock_10_emref = [af2_bound_gap_emref[cat][run][10][1] for run in haddock_rb_af2_runs]
    
#     # fractions
#     mod_b_haddock_frac_1_emref = [el/tot_runs for el in mod_b_haddock_1_emref]
#     mod_af2_haddock_frac_1_emref = [el/tot_runs for el in mod_af2_haddock_1_emref]
#     mod_b_haddock_frac_10_emref = [el/tot_runs for el in mod_b_haddock_10_emref]
#     mod_af2_haddock_frac_10_emref = [el/tot_runs for el in mod_af2_haddock_10_emref]
    
#     # zdock data
#     if cat == "CDR-EpiVag-AA":
#         zdock_cat = "CDR-EpiVag"
#     else:
#         zdock_cat = cat
#     mod_b_zdock_1 = [zdock_data[zdock_cat][run][1][0] for run in zdock_data[zdock_cat]]
#     mod_af2_zdock_1 = [zdock_data[zdock_cat][run][1][1] for run in zdock_data[zdock_cat]]
#     mod_b_zdock_10 = [zdock_data[zdock_cat][run][10][0] for run in zdock_data[zdock_cat]]
#     mod_af2_zdock_10 = [zdock_data[zdock_cat][run][10][1] for run in zdock_data[zdock_cat]]
    
#     #zdock fractions
#     mod_b_zdock_frac_1 = [el/tot_runs for el in mod_b_zdock_1]
#     mod_af2_zdock_frac_1 = [el/tot_runs for el in mod_af2_zdock_1]
#     mod_b_zdock_frac_10 = [el/tot_runs for el in mod_b_zdock_10]
#     mod_af2_zdock_frac_10 = [el/tot_runs for el in mod_af2_zdock_10]
    
#     axs[1][col].bar(r, mod_b_haddock_frac_1_emref, color = colors[1],
#                 width = width, edgecolor = 'black',
#                 label='HADDOCK Bound T1')
#     axs[1][col].bar(r + width, mod_af2_haddock_frac_1_emref, color = colors[3],
#                 width = width, edgecolor = 'black',
#                 label='HADDOCK AF2 T1')
#     axs[1][col].bar(r + width*2, mod_b_zdock_frac_1, color = colors[7],
#                 width = width, edgecolor = 'black',
#                 label='ZDOCK Bound T1')
#     axs[1][col].bar(r + width*3, mod_af2_zdock_frac_1, color = colors[9],
#                 width = width, edgecolor = 'black',
#                 label='ZDOCK AF2 T1')
#     # TOP 10 blurred
#     axs[1][col].bar(r, mod_b_haddock_frac_10_emref, color = colors[1], alpha=0.4,
#                 width = width, edgecolor = 'black',
#                 label='HADDOCK Bound T10')
#     axs[1][col].bar(r + width, mod_af2_haddock_frac_10_emref, color = colors[3], alpha=0.4,
#                 width = width, edgecolor = 'black',
#                 label='HADDOCK AF2 T10')
#     axs[1][col].bar(r + width*2, mod_b_zdock_frac_10, color = colors[7], alpha=0.4,
#                 width = width, edgecolor = 'black',
#                 label='ZDOCK Bound T10')
#     axs[1][col].bar(r + width*3, mod_af2_zdock_frac_10, color = colors[9], alpha=0.4,
#                 width = width, edgecolor = 'black',
#                 label='ZDOCK AF2 T10')

#     if col > 0:
#         axs[1][col].set_yticks([])
#     else:
#         axs[1][col].set_ylabel("Emref Success Rate", fontsize=16)
        
#     # adding bound bound results
#     axs[1][col].bar(r[-1]+1,
#                     bound_bound_emref[cat][1]/tot_runs,
#                     color="r",
#                     width = width,
#                     label="T1 BB")
#     axs[1][col].bar(r[-1]+1,
#                     bound_bound_emref[cat][10]/tot_runs,
#                     color="r",
#                     alpha=0.4,
#                     width = width,
#                     label="T10 BB")
    
#     new_r = np.arange(n)
#     new_r = np.append(new_r, n - 3/2*width)
#     axs[1][col].set_xticks(new_r + 3/2*width, xticks, fontsize=16, rotation = 90)
#     col+=1

# axs[0][0].set_title(cat_dict["Para-Epi"], fontsize=16)
# axs[0][2].set_title(cat_dict["CDR-EpiVag"], fontsize=16)
# axs[0][1].set_title(cat_dict["CDR-EpiVag-AA"], fontsize=16)

# axs[1][1].set_xlabel("Docking protocol", fontsize=24, labelpad=20)
# lines_labels = [axs[1][1].get_legend_handles_labels()]
# lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
# #fig.set_tight_layout(True)

# fig.legend(lines, labels, ncol=5, fontsize=16, loc="lower center", bbox_to_anchor = (0, 0.075, 1, 1))

# plt.savefig(Path("figures", "zdock_comparison", f"zdock_rigid_emref_1-10_fracbound.png"), dpi=300)