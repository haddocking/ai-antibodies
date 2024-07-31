import MDAnalysis as mda
from MDAnalysis.analysis import align
from pathlib import Path

import MDAnalysis.analysis.rms
import numpy as np
from sklearn.cluster import AgglomerativeClustering
from scipy.spatial.distance import squareform
import os
import logging

log = logging.getLogger("aiabslog")

def get_clustering_dict(clusters, ligands):
    """
    Gets dictionary of clusters.

    Parameters
    ----------
    clusters : list
        list of cluster IDs
    ligands : list
        names of the ligands

    Returns
    -------
    cl_dict : dict
        dictionary of clustered interfaces
        *example* {
            1 : [
            'interface_1', 'interface_3'
                ] ,
            2 : [
            'interface_2'
            ],
            ...
            }
    """
    cl_dict = {}
    # loop over clusters
    for cl in range(len(clusters)):
        if clusters[cl] not in cl_dict.keys():
            cl_dict[clusters[cl]] = [ligands[cl]]
        else:
            cl_dict[clusters[cl]].append(ligands[cl])
    #log.info(f"Cluster dictionary {cl_dict}")
    return cl_dict

def cond_index(i: int, j: int, n: int) -> float:
    """
    Get the condensed index from two matrix indexes.

    Parameters
    ----------
    i : int
        Index of the first element.
    j : int
        Index of the second element.
    n : int
        Number of observations.
    """
    return n * (n - 1) / 2 - (n - i) * (n - i - 1) / 2 + j - i - 1



def get_cluster_center(npw, rmsd_matrix, n_obs):
    """
    Get the cluster centers.

    Parameters
    ----------
    npw: np.ndarray
        Indexes of the cluster over cluster_list array
    n_obs : int
        Number of overall observations (models).
    rmsd_matrix : np.ndarray
        RMSD matrix.

    Returns
    -------
    cluster_center : int
        Index of cluster center
    """
    intra_cl_distances = {el: 0.0 for el in npw}

    # iterating over the elements of the cluster
    for m_idx in range(len(npw)):
        npws = npw[m_idx + 1:]
        pairs = [int(cond_index(npw[m_idx], npw_el, n_obs)) for npw_el in npws]
        for pair_idx in range(len(pairs)):
            intra_cl_distances[npw[m_idx]] += rmsd_matrix[pairs[pair_idx]]
            intra_cl_distances[npws[pair_idx]] += rmsd_matrix[pairs[pair_idx]]
    cluster_center = min(intra_cl_distances, key=intra_cl_distances.get)  # type: ignore  # noqa : E501
    return cluster_center

def write_dict(input_dict, out_filename, keyword, sep=" "):
    """
    Writes dictionary to file.

    Parameters
    ----------
    input_dict : dict
        input dictionary
    out_filename : str or Path
        name of the output filename
    keyword : str
        keyword to be used before each entry
    """
    #log.info(f"Writing {keyword} information to file {out_filename}")
    with open(out_filename, "w") as wfile:
        for key in input_dict.keys():
            cl_string = sep.join([str(el) for el in input_dict[key]])
            wfile.write(f"{keyword} {key} -> " + cl_string + os.linesep)


def calculate_quantities(model_uni, ref_uni, loop_residues, heavy_end=118):

    # overall rmsd
    sel = f"name CA and protein and not resid {heavy_end}"
    model_bkb = model_uni.select_atoms(sel)
    ref_bkb = ref_uni.select_atoms(sel)
    #if model_bkb.n_atoms != ref_bkb.n_atoms:
    #    for n in range(model_bkb.n_atoms):
    #        print(f"{model_bkb[n].resid}-{model_bkb[n].name}, {ref_bkb[n].resid}-{ref_bkb[n].name}")
    rmsd_full = MDAnalysis.analysis.rms.rmsd(model_uni.select_atoms(sel).positions, ref_uni.select_atoms(sel).positions, superposition=True)
    
    model_residues = model_uni.select_atoms(f"protein and name CA and chainID A and not resid {heavy_end}").resids
    ref_residues = ref_uni.select_atoms(f"protein and name CA and chainID A and not resid {heavy_end}").resids
    
    mod_framework_hresids = [model_residues[n] for n in range(len(model_residues)) if n not in loop_residues]
    
    ref_framework_hresids = [ref_residues[n] for n in range(len(ref_residues)) if n not in loop_residues]
    
    mod_framework_hatoms = model_uni.select_atoms(f"protein and name CA and chainID A and resid {' '.join([str(el) for el in mod_framework_hresids])}")
    ref_framework_hatoms = ref_uni.select_atoms(f"protein and name CA and chainID A and resid {' '.join([str(el) for el in ref_framework_hresids])}")
    
    mod_com = model_uni.atoms.center_of_mass()
    ref_com = ref_uni.atoms.center_of_mass()
    
    R, rmsd_framework = align.rotation_matrix(mod_framework_hatoms.positions - mod_framework_hatoms.center_of_mass(), ref_framework_hatoms.positions - ref_framework_hatoms.center_of_mass())
    
    #loop
    mod_loop_sel = f"not name H* and backbone and chainID A and ("
    loop_resids = [f" resid {loop_residues[n]} " for n in range(len(loop_residues))]
    mod_loop_sel += " or ".join(loop_resids)
    mod_loop_sel += " )"
    mod_h1_atoms = model_uni.select_atoms(mod_loop_sel)
    ref_h1_atoms = ref_uni.select_atoms(mod_loop_sel)

    # loop rmsd
    mod_loop_com_coords = mod_h1_atoms.select_atoms('not name H* and backbone').center_of_mass()
    ref_loop_com_coords = ref_h1_atoms.select_atoms('not name H* and backbone').center_of_mass()
    mod_h1_atoms.atoms.translate(-mod_loop_com_coords)
    mod_h1_atoms.atoms.rotate(R)
    mod_h1_atoms.atoms.translate(ref_loop_com_coords)
    rmsd_loop = mda.analysis.rms.rmsd(mod_h1_atoms.positions, ref_h1_atoms.positions, superposition=False)
    
    return rmsd_full, rmsd_framework, rmsd_loop


def cluster_antibodies(pdb_dict, loop_resids):
    
    pdb_keys = list(pdb_dict.keys())
    npdbs = len(pdb_keys)
    log.info(f"{npdbs} pdbs to cluster")
    npairs = npdbs * (npdbs - 1) // 2
    loop_rmsd_matrix = np.zeros(npairs)
    knt = 0
    for i in range(len(pdb_keys)):
        pdb_file = pdb_dict[pdb_keys[i]]
        for j in range(i+1, npdbs):
            pdb_file2 = pdb_dict[pdb_keys[j]]
            model_uni = mda.Universe(pdb_file)
            ref_uni = mda.Universe(pdb_file2)
            rmsd_full, rmsd_framework, rmsd_loop = calculate_quantities(model_uni, ref_uni, loop_resids)
            loop_rmsd_matrix[knt] = rmsd_loop
            knt += 1
    log.info(f"average loop rmsd {np.mean(loop_rmsd_matrix):.2f} std {np.std(loop_rmsd_matrix):.2f}")
    clustering = AgglomerativeClustering(n_clusters=4,
                                     metric="precomputed",
                                     linkage="complete").fit(squareform(loop_rmsd_matrix))
    # plot dendrogram
    cluster_centers_data = []

    clusters = clustering.labels_
    cl_dict = get_clustering_dict(clusters, pdb_keys)
    
    assert len(np.unique(clusters)) == 4
    # now for each cluster we want to calculate the dispersion of the h3 rmsd vs the reference
    for cl in cl_dict:
        log.info(f"cluster {cl} has {len(cl_dict[cl])} elements")
        npw = np.where(clusters == cl)[0]
        cluster_center_idx = get_cluster_center(npw, loop_rmsd_matrix, npdbs)
        cluster_center = pdb_keys[cluster_center_idx]

        cluster_centers_data.append(cluster_center)
    ensemble_dict = {}
    for el in cluster_centers_data:
        ensemble_dict[el] = pdb_dict[el]
    write_dict(cl_dict, Path(pdb_dict[el].parent, f"clustering_complete_{npdbs}_4.txt"), "h3_cluster")
    return ensemble_dict
            