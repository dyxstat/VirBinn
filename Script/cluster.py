from sklearn.metrics import silhouette_score
import numpy as np
import pandas as pd
import scipy.sparse as scisp
import igraph as ig
import leidenalg
import os
import logging
from Script.utils import gen_bins

logger = logging.getLogger(__name__)

# build graph and cluster by silhouette-driven Leiden
def bin(fasta, path, precompute = False, output_prefix = 'virbinn'):
    combined = scisp.load_npz(os.path.join(path, "tmp", 'combined.npz'))
    viral_ci = pd.read_csv(os.path.join(path, 'viral_contig_info.csv'), header=0).values
    
    vcount = combined.shape[0]
    src, tgt, wei = combined.row, combined.col, combined.data
    mask = src > tgt
    
    edges = list(zip(src[mask], tgt[mask]))
    wei = wei[mask]
    g = ig.Graph(vcount, edges)
    
    combined = combined.toarray().astype(float)
    if precompute:
        logger.info("Use precomputed matrix")
        combined /= combined.max()
        combined = 1-combined

    sil_scores = []
    resolutions = list(range(2, 50))
    for r in resolutions:
        part = leidenalg.find_partition(
            g,
            leidenalg.RBConfigurationVertexPartition,
            weights=wei,
            resolution_parameter=r,
            n_iterations=-1
        )
        labels = np.zeros(vcount, dtype=int)
        
        for cid, cluster in enumerate(part):
            for contig in cluster:
                labels[contig] = cid
            
        if len(np.unique(labels)) == vcount:
            break
        
        if precompute:          
            sil_scores.append(silhouette_score(combined, labels), metric="precomputed")
        else:
            sil_scores.append(silhouette_score(combined, labels))
        
    best_r = resolutions[int(np.argmax(sil_scores))]
    final_part = leidenalg.find_partition(
        g,
        leidenalg.RBConfigurationVertexPartition,
        weights=wei,
        resolution_parameter=best_r,
        n_iterations=-1
    )
    logger.info("Chosen resolution %d => %d viral bins",
                best_r, len(final_part))

    # write cluster file
    assigned = []
    cluster_path = os.path.join(path, f"{output_prefix}_clusters.txt")
    with open(cluster_path, "w") as out:
        for cid, cluster in enumerate(final_part):
            for idx in cluster:
                contig_name = viral_ci[idx, 0]
                assigned.append(contig_name)
                out.write(f"{contig_name}\tgroup{cid}\n")

        # any unassigned
        un_idx = len(final_part)
        for name in viral_ci[: , 0]:
            if name not in assigned:
                out.write(f"{name}\tgroup{un_idx}\n")
                un_idx += 1

    logger.info("Cluster assignments written to %s", cluster_path)
    
    logger.info("Generating bins from %s", cluster_path)
    
    bin_dir = os.path.join(path, f"{output_prefix}_VIRAL_BIN")
    
    gen_bins(fasta, cluster_path, bin_dir)
    