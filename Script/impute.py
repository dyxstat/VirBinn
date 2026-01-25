# Script/impute_module.py

import os
import pandas as pd
import numpy as np
import scipy.sparse as scisp
import logging

from Script.utils import (
    imputation,
    top_percent_edges
)

logger = logging.getLogger(__name__)


def impute(viral_file, path, discard_viral=50, discard_host=80):
    """
    Perform the imputation and generate a cluster assignment file.
    """
    
    # load
    viral_names = pd.read_csv(viral_file, sep='\t', header=None)[0].tolist()
    ci = pd.read_csv(os.path.join(path, 'contig_info.csv'), header=0).values
    mat = scisp.load_npz(os.path.join(path, 'Raw_contact_matrix.npz')).tocoo()

    # normalize by sqrt(diag_i * diag_j)
    diag = mat.diagonal() + 1
    sqrt_diag = np.sqrt(diag)
    diag_max = diag.max()
    rows, cols, data = mat.row, mat.col, mat.data
    new_data = np.empty_like(data, dtype=float)
    for k in range(len(data)):
        denom = (sqrt_diag[rows[k]] * sqrt_diag[cols[k]] / diag_max)
        new_data[k] = data[k] / denom if denom > 0 else 0.0
    mat = scisp.coo_matrix((new_data, (rows, cols)), shape=mat.shape)
    mat = mat.tolil()
    mat.setdiag(0)

    # split indices
    viral_idx, host_idx = [], []
    for i, row in enumerate(ci):
        if row[0] in viral_names:
            viral_idx.append(i)
        else:
            host_idx.append(i)
    order = viral_idx + host_idx
    
    
    mat = mat.tocsr()
    viral_mat = mat[viral_idx, :][:, viral_idx]
    viral_ci = ci[viral_idx , :]
    
    mat = mat[order, :][:, order]
    ci = ci[order, :]
    
    n_virus = len(viral_idx)

    # impute
    im_v = imputation(viral_mat, n_virus, 0.5, discard_viral, masked=False)
    im_h = imputation(mat, n_virus, 0.5, discard_host, masked=True)
    for M in (im_v, im_h):
        M.setdiag(0)


    raw_count = viral_mat.tocoo().data.size // 2

    fv = top_percent_edges(im_v, raw_count)
    fh = top_percent_edges(im_h, raw_count)

    # combine and symmetrize
    vb = scisp.coo_matrix(
        (np.ones_like(viral_mat.tocoo().data),
         (viral_mat.tocoo().row, viral_mat.tocoo().col)),
        shape=viral_mat.shape
    )
    tmp = vb.maximum(fv)
    combined = tmp.toarray() + fh.toarray()
    combined = scisp.coo_matrix(combined)
    logger.info("Number of edges in (viral, tmp, host, combined): %s, %s, %s, %s",
                fv.tocoo().data.size // 2,
                tmp.tocoo().data.size // 2,
                fh.tocoo().data.size // 2,
                combined.tocoo().data.size // 2)
    
    # Save matrices as sparse .npz files
    scisp.save_npz(os.path.join(path, 'tmp', 'fv.npz'), fv)
    scisp.save_npz(os.path.join(path, 'tmp', 'fh.npz'), fh)
    scisp.save_npz(os.path.join(path, 'tmp', 'tmp.npz'), tmp)
    scisp.save_npz(os.path.join(path, 'tmp', 'viral_mat.npz'), viral_mat)
    scisp.save_npz(os.path.join(path, 'tmp', 'combined.npz'), combined)

    with open(os.path.join(path, 'viral_contig_info.csv'),'w') as out:
        out.write(str('Viral contig name')+ ','  + str('Number of restriction sites')+ ',' +str('Contig length') + '\n')
        for row in viral_ci:
            out.write(str(row[0])+ ',' + str(row[1])+ ',' + str(row[2]))
            out.write('\n')
    