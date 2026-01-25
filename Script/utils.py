import time
import math
import os
import io
import gzip
import subprocess
import pickle
import bz2
import sys
import numpy as np
import scipy.sparse as scisp
from scipy.sparse.linalg import norm


def save_object(file_name, obj):
    """
    Serialize an object to a file with gzip compression. .gz will automatically be
    added if missing.

    :param file_name: output file name
    :param obj: object to serialize
    """
    with open_output(file_name, compress='gzip') as out_h:
        pickle.dump(obj, out_h)


def load_object(file_name):
    """
    Deserialize an object from a file with automatic support for compression.

    :param file_name: input file name
    :return: deserialzied object
    """
    with open_input(file_name) as in_h:
        return pickle.load(in_h)


def open_input(file_name):
    """
    Open a text file for input. The filename is used to indicate if it has been
    compressed. Recognising gzip and bz2.

    :param file_name: the name of the input file
    :return: open file handle, possibly wrapped in a decompressor
    """
    suffix = file_name.split('.')[-1].lower()
    if suffix == 'bz2':
        return bz2.BZ2File(file_name, 'r')
    elif suffix == 'gz':
        return gzip.GzipFile(file_name, 'r')
    else:
        return open(file_name, 'r')


def open_output(file_name, append=False, compress=None, gzlevel=6):
    """
    Open a text stream for reading or writing. Compression can be enabled
    with either 'bzip2' or 'gzip'. Additional option for gzip compression
    level. Compressed filenames are only appended with suffix if not included.

    :param file_name: file name of output
    :param append: append to any existing file
    :param compress: gzip, bzip2
    :param gzlevel: gzip level (default 6)
    :return:
    """

    mode = 'w' if not append else 'w+'

    if compress == 'bzip2':
        if not file_name.endswith('.bz2'):
            file_name += '.bz2'
        # bz2 missing method to be wrapped by BufferedWriter. Just directly
        # supply a buffer size
        return bz2.BZ2File(file_name, mode, buffering=65536)
    elif compress == 'gzip':
        if not file_name.endswith('.gz'):
            file_name += '.gz'
        return io.BufferedWriter(gzip.GzipFile(file_name, mode, compresslevel=gzlevel))
    else:
        return io.BufferedWriter(io.FileIO(file_name, mode))


def make_dir(path, exist_ok=False):
    """
    Convenience method for making directories with a standard logic.
    An exception is raised when the specified path exists and is not a directory.
    :param path: target path to create
    :param exist_ok: if true, an existing directory is ok. Existing files will still cause an exception
    """
    if not os.path.exists(path):
        os.mkdir(path)
    elif not exist_ok:
        raise IOError('output directory already exists!')
    elif os.path.isfile(path):
        raise IOError('output path already exists and is a file!')


def app_path(subdir, filename):
    """
    Return path to named executable in a subdirectory of the running application

    :param subdir: subdirectory of application path
    :param filename: name of file
    :return: absolute path
    """
    return os.path.join(sys.path[0], subdir, filename)



def count_fasta_sequences(file_name):
    """
    Estimate the number of fasta sequences in a file by counting headers. Decompression is automatically attempted
    for files ending in .gz. Counting and decompression is by why of subprocess calls to grep and gzip. Uncompressed
    files are also handled. This is about 8 times faster than parsing a file with BioPython and 6 times faster
    than reading all lines in Python.

    :param file_name: the fasta file to inspect
    :return: the estimated number of records
    """
    if file_name.endswith('.gz'):
        proc_uncomp = subprocess.Popen(['gzip', '-cd', file_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        proc_read = subprocess.Popen(['grep', r'^>'], stdin=proc_uncomp.stdout, stdout=subprocess.PIPE)
    else:
        proc_read = subprocess.Popen(['grep', r'^>', file_name], stdout=subprocess.PIPE)

    n = 0
    for _ in proc_read.stdout:
        n += 1
    return n

def floor_to_2dec(x):
    return math.floor(x * 100) / 100

def calc_sparsity(matrix):
    row, col = matrix.shape
    sparsity = matrix.nnz / row / col
    return sparsity


def subtract_and_zero_diag(matrix):
    new_matrix = 1 - matrix.copy()
    np.fill_diagonal(new_matrix, 0)
    return new_matrix    
 

def mask_marker_entries_to_zero(adj, num_marker_contig):
    """
    Set entries in a sparse matrix to zero where both row and column indices
    are less than num_marker_contig, but keep the original matrix shape.

    Parameters:
        adj (scipy.sparse.spmatrix): Input sparse matrix.
        num_marker_contig (int): Threshold for masking.

    Returns:
        scipy.sparse.coo_matrix: Masked sparse matrix in COO format.
    """
    # Convert to COO format
    adj = adj.tocoo()

    # Mask: entries where both row and col < num_marker_contig
    mask = ~((adj.row <= num_marker_contig) & (adj.col <= num_marker_contig))

    # Keep only unmasked entries
    new_row = adj.row[mask]
    new_col = adj.col[mask]
    new_data = adj.data[mask]

    # Return matrix with original shape
    return scisp.coo_matrix((new_data, (new_row, new_col)), shape=adj.shape)


def random_walk(P, rp, perc, tol, num_marker_contig, masked=False):
    #masked represents whether remove the first num_marker_contig entries
    _record = 0
    _start_time = time.time()
    
    if masked == True:
        P = mask_marker_entries_to_zero(P, num_marker_contig)
    
    row = np.arange(num_marker_contig)
    col = np.arange(num_marker_contig)
    data = np.ones(num_marker_contig)
    
    I = scisp.coo_matrix((data, (row, col)), shape=(num_marker_contig, P.shape[0]), dtype=np.float32)
    Q = I.copy()
    delta_previous = 0
    
    for i in range(500):
        Q_new = (1 - rp) * Q.dot(P) + rp * I
        delta = norm(Q - Q_new)
        Q = Q_new.copy()
        
        if i >= 1:
            min_cutoff = np.percentile(Q.tocoo().data, perc)
            if min_cutoff > 0:
                s_before = calc_sparsity(Q)
                Q = Q.multiply(Q > min_cutoff)
                s_after = calc_sparsity(Q)
            print('Iteration{}: sparsity before is {}; sparsity after is {}; delta is {} \n'.format(i , s_before, s_after, delta))
            
        if delta < tol:
            _end_time = time.time()
            print('Random walk ends because of reaching the tolerance.')
            break

        if np.abs(delta - delta_previous) < 0.001:
            _record += 1
            
        delta_previous = delta
        if _record == 5:
            print('Random walk ends because of early stop.')
            break
    _end_time = time.time()
    Q = Q.tocsr()[np.arange(num_marker_contig) , :]
    Q = Q.tocsc()[: , np.arange(num_marker_contig)]
    print('Inputaion takes {} time and {} steps; loss is {}; sparsity of imputed matrix is {}'.format(_end_time-_start_time, i+1, delta, calc_sparsity(Q)))
    return Q


def imputation(hic_matrix, num_viral_contig, rwr_rp, rwr_perc, masked=False):
    '''
    Impute normcc-normalized Hi-C contact matrix using randow walk with restart
    '''
    ###Only impute part of matrix instead of the whole matrix

    A = hic_matrix.copy()
    A = A - scisp.diags(A.diagonal())
    B = A + scisp.diags((A.sum(axis=0).A.ravel() == 0).astype(int))
    d = scisp.diags(1 / B.sum(axis=0).A.ravel())
    P = d.dot(B).astype(np.float32)
    Q = random_walk(P , rwr_rp , rwr_perc , 0.01 , num_viral_contig, masked=masked)

    E = Q.copy()
    E += E.T
    d = E.sum(axis=0).A.ravel()
    d[d == 0] = 1
    b = scisp.diags(1 / np.sqrt(d))
    E = b.dot(E).dot(b)
    
    return E.tolil()


def combine_top_edges_with_prior(adj, X, keep_ratio=0.3):
    # 确保 COO 格式
    adj = adj.tocoo()
    X = X.tocoo()

    shape = adj.shape

    # Step 1: 获取 X 中所有边的 index（只看上三角避免重复）
    mask_X = X.row < X.col
    X_pairs = set(zip(X.row[mask_X], X.col[mask_X]))

    # Step 2: 获取 adj 中上三角非 X 的边
    mask_adj = (adj.row < adj.col)
    row = adj.row[mask_adj]
    col = adj.col[mask_adj]
    data = adj.data[mask_adj]

    # 过滤掉 X 中已包含的边
    edge_candidates = []
    for i, (r, c, w) in enumerate(zip(row, col, data)):
        if (r, c) not in X_pairs:
            edge_candidates.append((r, c, w))

    # Step 3: 选出 top 30% 边
    edge_candidates.sort(key=lambda x: x[2], reverse=True)
    k = int(len(edge_candidates) * keep_ratio)
    selected_edges = edge_candidates[:k]

    # Step 4: 合并 X 中边 和 选出的边（确保对称性）
    combined_rows = []
    combined_cols = []

    # 加入 X 中边
    for r, c in X_pairs:
        if r!= c:
            combined_rows += [r, c]
            combined_cols += [c, r]

    # 加入 top edges
    for r, c, _ in selected_edges:
        combined_rows += [r, c]
        combined_cols += [c, r]

    # 所有值设为 1
    data_ones = np.ones(len(combined_rows), dtype=int)

    combined = scisp.coo_matrix((data_ones, (combined_rows, combined_cols)), shape=shape)

    return combined



def top_percent_edges(adj, num_retained_edges=3000):
    # Ensure COO format
    adj = adj.tocoo()
    shape = adj.shape

    # Only take upper triangle to avoid duplicates
    mask = adj.row < adj.col
    row = adj.row[mask]
    col = adj.col[mask]
    data = adj.data[mask]

    # Create list of (i, j, weight)
    edges = list(zip(row, col, data))

    # Sort by weight descending and keep top 10%
    edges.sort(key=lambda x: x[2], reverse=True)
    k = int(min(num_retained_edges, len(edges)))
    selected = edges[:k]

    # Symmetrize and set all values to 1
    final_rows = []
    final_cols = []

    for r, c, _ in selected:
        final_rows += [r, c]
        final_cols += [c, r]

    data_ones = np.ones(len(final_rows), dtype=int)
    combined = scisp.coo_matrix((data_ones, (final_rows, final_cols)), shape=shape)

    return combined


def gen_bins(fastafile,resultfile,outputdir):
    # read fasta file
    sequences={}
    if fastafile.endswith("gz"):
        with gzip.open(fastafile,'r') as f:
            for line in f:
                line=str(line,encoding="utf-8")
                if line.startswith(">"):
                    if " " in line:
                        seq,others=line.split(' ', 1)
                        sequences[seq] = ""
                    else :
                        seq=line.rstrip("\n")
                        sequences[seq] = ""
                else:
                    sequences[seq] += line.rstrip("\n")
    else:
        with open(fastafile,'r') as f:
            for line in f:
                if line.startswith(">"):
                    if " " in line:
                        seq,others=line.split(' ', 1)
                        sequences[seq] = ""
                    else :
                        seq=line.rstrip("\n")
                        sequences[seq] = ""
                else:
                    sequences[seq] += line.rstrip("\n")
    dic={}
    with open(resultfile,"r") as f:
        for line in f:
            contig_name,cluster_name=line.strip().split('\t')
            try:
                dic[cluster_name].append(contig_name)
            except:
                dic[cluster_name]=[]
                dic[cluster_name].append(contig_name)
    print("Writing bins in \t{}".format(outputdir))
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    
    bin_name=0
    for _,cluster in dic.items():
        if bin_name < 10:
            bin = 'VIRAL_BIN'+ '000' + str(bin_name) + '.fa'
        elif bin_name >= 10 and bin_name < 100:
            bin = 'VIRAL_BIN'+ '00' + str(bin_name) + '.fa'
        elif bin_name >= 100 and bin_name < 1000:
            bin = 'VIRAL_BIN'+ '0' + str(bin_name) + '.fa'
        else:
            bin = 'VIRAL_BIN'+str(bin_name) + '.fa'
        binfile=os.path.join(outputdir,"{}".format(bin))
        with open(binfile,"w") as f:
            for contig_name in cluster:
                contig_name=">"+contig_name
                try:
                    sequence=sequences[contig_name]
                except:
                    continue
                f.write(contig_name+"\n")
                f.write(sequence+"\n")
        bin_name+=1