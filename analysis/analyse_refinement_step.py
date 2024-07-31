import numpy as np
from functions import load_data, create_bound_gap_dictionary, create_bound_bound_dictionary
NPDBS = 82
print("Para-Epi scenario")
xticks = ["ABB", "ABBE", "ABL", "AF2", "IG", "ENS", "ENSNOAF2", "ENS196-48", "ENS196-CLT"]

# LOAD DATA
rigidbody_capri, rigidbody_capri_bound, emref_capri, emref_capri_bound, zdock_ss, emref_rigid_capri, af2multimer_ss = load_data()
tot_runs = np.unique(rigidbody_capri["pdb"]).shape[0] # should be 79
print(f"total number of runs {tot_runs}")

assert tot_runs == NPDBS
# extract rigidbody, flexref and emref data
acc_key="acc"
cat = "Para-Epi"
sorted_runs = [f'run-af2ab-{cat}-mpi-50-50', f'run-af2abens-{cat}-mpi-50-50', f'run-af2abl-{cat}-mpi-50-50',
                f'run-af2af2-{cat}-mpi-50-50', f'run-af2ig-{cat}-mpi-50-50', f'run-af2ens-{cat}-mpi-50-50', 
                f'run-af2ensnoaf2-{cat}-mpi-50-50', f'run-af2ens-{cat}-mpi-196-48', f'run-af2ens-{cat}-mpi-196-clt']

af2_bound_gap_rigid = create_bound_gap_dictionary(rigidbody_capri, acc_key=acc_key)
af2_bound_gap_emref = create_bound_gap_dictionary(emref_capri, acc_key=acc_key)
bound_bound_rigid = create_bound_bound_dictionary(rigidbody_capri_bound, acc_key=acc_key)
bound_bound_emref = create_bound_bound_dictionary(emref_capri_bound, acc_key=acc_key)

delta_bound_top1, delta_bound_top10 = [], []
delta_af2_top1, delta_af2_top10 = [], []
for n in range(len(xticks)):
    delta_bound_top1.append(af2_bound_gap_emref[cat][sorted_runs[n]][1][0]/tot_runs - af2_bound_gap_rigid[cat][sorted_runs[n]][1][0]/tot_runs)
    delta_bound_top10.append(af2_bound_gap_emref[cat][sorted_runs[n]][10][0]/tot_runs - af2_bound_gap_rigid[cat][sorted_runs[n]][10][0]/tot_runs)
    delta_af2_top1.append(af2_bound_gap_emref[cat][sorted_runs[n]][1][1]/tot_runs - af2_bound_gap_rigid[cat][sorted_runs[n]][1][1]/tot_runs)
    delta_af2_top10.append(af2_bound_gap_emref[cat][sorted_runs[n]][10][1]/tot_runs - af2_bound_gap_rigid[cat][sorted_runs[n]][10][1]/tot_runs)

print(f"BOUND: avg Delta_top1 {np.mean(delta_bound_top1):.2f} & avg Delta_top10 {np.mean(delta_bound_top10):.2f} \\\\")
print(f"AF2: avg Delta_top1 {np.mean(delta_af2_top1):.2f} & avg Delta_top10 {np.mean(delta_af2_top10):.2f} \\\\")

# CDR-EpiVag-AA
print("\nCDR-EpiVag-AA scenario\n")
cat = "CDR-EpiVag-AA"
sorted_runs = [f'run-af2ab-{cat}-mpi-50-50', f'run-af2abens-{cat}-mpi-50-50', f'run-af2abl-{cat}-mpi-50-50',
                f'run-af2af2-{cat}-mpi-50-50', f'run-af2ig-{cat}-mpi-50-50', f'run-af2ens-{cat}-mpi-50-50', 
                f'run-af2ensnoaf2-{cat}-mpi-50-50', f'run-af2ens-{cat}-mpi-196-48', f'run-af2ens-{cat}-mpi-196-clt']

delta_cdr_bound_top1, delta_cdr_bound_top10 = [], []
delta_cdr_af2_top1, delta_cdr_af2_top10 = [], []
for n in range(len(xticks)):
    delta_cdr_bound_top1.append(af2_bound_gap_emref[cat][sorted_runs[n]][1][0]/tot_runs - af2_bound_gap_rigid[cat][sorted_runs[n]][1][0]/tot_runs)
    delta_cdr_bound_top10.append(af2_bound_gap_emref[cat][sorted_runs[n]][10][0]/tot_runs - af2_bound_gap_rigid[cat][sorted_runs[n]][10][0]/tot_runs)
    delta_cdr_af2_top1.append(af2_bound_gap_emref[cat][sorted_runs[n]][1][1]/tot_runs - af2_bound_gap_rigid[cat][sorted_runs[n]][1][1]/tot_runs)
    delta_cdr_af2_top10.append(af2_bound_gap_emref[cat][sorted_runs[n]][10][1]/tot_runs - af2_bound_gap_rigid[cat][sorted_runs[n]][10][1]/tot_runs)
    
print(f"BOUND: avg Delta_top1 {np.mean(delta_cdr_bound_top1):.2f} & avg Delta_top10 {np.mean(delta_cdr_bound_top10):.2f} \\\\")
print(f"AF2: avg Delta_top1 {np.mean(delta_cdr_af2_top1):.2f} & avg Delta_top10 {np.mean(delta_cdr_af2_top10):.2f} \\\\")

# now with the medium-quality poses
# Para-Epi
print("\n\nmedium-quality poses")
print("\nPara-Epi scenario\n")
acc_key="med"
cat = "Para-Epi"
sorted_runs = [f'run-af2ab-{cat}-mpi-50-50', f'run-af2abens-{cat}-mpi-50-50', f'run-af2abl-{cat}-mpi-50-50',
                f'run-af2af2-{cat}-mpi-50-50', f'run-af2ig-{cat}-mpi-50-50', f'run-af2ens-{cat}-mpi-50-50', 
                f'run-af2ensnoaf2-{cat}-mpi-50-50', f'run-af2ens-{cat}-mpi-196-48', f'run-af2ens-{cat}-mpi-196-clt']

af2_bound_gap_rigid = create_bound_gap_dictionary(rigidbody_capri, acc_key=acc_key)
af2_bound_gap_emref = create_bound_gap_dictionary(emref_capri, acc_key=acc_key)
bound_bound_rigid = create_bound_bound_dictionary(rigidbody_capri_bound, acc_key=acc_key)
bound_bound_emref = create_bound_bound_dictionary(emref_capri_bound, acc_key=acc_key)

delta_bound_top1, delta_bound_top10 = [], []
delta_af2_top1, delta_af2_top10 = [], []
for n in range(len(xticks)):
    delta_bound_top1.append(af2_bound_gap_emref[cat][sorted_runs[n]][1][0]/tot_runs - af2_bound_gap_rigid[cat][sorted_runs[n]][1][0]/tot_runs)
    delta_bound_top10.append(af2_bound_gap_emref[cat][sorted_runs[n]][10][0]/tot_runs - af2_bound_gap_rigid[cat][sorted_runs[n]][10][0]/tot_runs)
    delta_af2_top1.append(af2_bound_gap_emref[cat][sorted_runs[n]][1][1]/tot_runs - af2_bound_gap_rigid[cat][sorted_runs[n]][1][1]/tot_runs)
    delta_af2_top10.append(af2_bound_gap_emref[cat][sorted_runs[n]][10][1]/tot_runs - af2_bound_gap_rigid[cat][sorted_runs[n]][10][1]/tot_runs)

print(f"BOUND: avg Delta_top1 {np.mean(delta_bound_top1):.2f} & avg Delta_top10 {np.mean(delta_bound_top10):.2f} \\\\")
print(f"AF2: avg Delta_top1 {np.mean(delta_af2_top1):.2f} & avg Delta_top10 {np.mean(delta_af2_top10):.2f} \\\\")

#CDR-EpiVag-AA
print("\nCDR-EpiVag-AA scenario\n")
cat = "CDR-EpiVag-AA"
sorted_runs = [f'run-af2ab-{cat}-mpi-50-50', f'run-af2abens-{cat}-mpi-50-50', f'run-af2abl-{cat}-mpi-50-50',
                f'run-af2af2-{cat}-mpi-50-50', f'run-af2ig-{cat}-mpi-50-50', f'run-af2ens-{cat}-mpi-50-50', 
                f'run-af2ensnoaf2-{cat}-mpi-50-50', f'run-af2ens-{cat}-mpi-196-48', f'run-af2ens-{cat}-mpi-196-clt']


delta_cdr_bound_top1, delta_cdr_bound_top10 = [], []
delta_cdr_af2_top1, delta_cdr_af2_top10 = [], []
for n in range(len(xticks)):
    delta_cdr_bound_top1.append(af2_bound_gap_emref[cat][sorted_runs[n]][1][0]/tot_runs - af2_bound_gap_rigid[cat][sorted_runs[n]][1][0]/tot_runs)
    delta_cdr_bound_top10.append(af2_bound_gap_emref[cat][sorted_runs[n]][10][0]/tot_runs - af2_bound_gap_rigid[cat][sorted_runs[n]][10][0]/tot_runs)
    delta_cdr_af2_top1.append(af2_bound_gap_emref[cat][sorted_runs[n]][1][1]/tot_runs - af2_bound_gap_rigid[cat][sorted_runs[n]][1][1]/tot_runs)
    delta_cdr_af2_top10.append(af2_bound_gap_emref[cat][sorted_runs[n]][10][1]/tot_runs - af2_bound_gap_rigid[cat][sorted_runs[n]][10][1]/tot_runs)

print(f"BOUND: avg Delta_top1 {np.mean(delta_cdr_bound_top1):.2f} & avg Delta_top10 {np.mean(delta_cdr_bound_top10):.2f} \\\\")
print(f"AF2: avg Delta_top1 {np.mean(delta_cdr_af2_top1):.2f} & avg Delta_top10 {np.mean(delta_cdr_af2_top10):.2f} \\\\")