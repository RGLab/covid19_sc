import scipy.sparse as sparse
import h5py as h5py
from scipy.sparse import vstack
import numpy as np

merge_path = '/fh/fast/gottardo_r/ytian_working/covid19_datasets/h5seurat_reprocessed/covid19_datasets.h5seurat'
big_path = '/fh/fast/gottardo_r/ytian_working/covid19_datasets/h5seurat_reprocessed/bigSeu_tomerge.h5seurat'
stephenson_path = '/fh/fast/gottardo_r/ytian_working/covid19_datasets/h5seurat_reprocessed/stephenson_2021_tomerge.h5seurat'

# Pull out sparse matrix and convert to 64bit int
def get_matrix(h5seurat_path, matrix): 
    f = h5py.File(h5seurat_path, 'r')
    matrix = f['assays/SCT/' + matrix]
    data = matrix['data']
    indices = matrix['indices']
    indptr = matrix['indptr']
    mx = sparse.csr_matrix((data, indices, indptr))
    mx.indices = mx.indices.astype(np.dtype('int64'))
    mx.indptr = mx.indptr.astype(np.dtype('int64'))
    return mx

# concat two csr matrices
def concatenate_csr_matrices_by_row(matrix1, matrix2):
    new_data = np.concatenate((matrix1.data, matrix2.data))
    new_indices = np.concatenate((matrix1.indices, matrix2.indices))
    new_ind_ptr = matrix2.indptr + len(matrix1.data)
    new_ind_ptr = new_ind_ptr[1:]
    new_ind_ptr = np.concatenate((matrix1.indptr, new_ind_ptr))
    return sparse.csr_matrix((new_data, new_indices, new_ind_ptr))

# Features are the same between the two, so grab them from either. 
features = h5py.File(big_path, 'r')['assays/SCT/features'][...]

# write result to new assay in merged file
hfile = h5py.File(merge_path)
sct = hfile.create_group("assays/SCT")
sct.create_dataset("features", data = features)
sct.attrs.create("key", ['sct_'])
sct.attrs.create("s4class", ['Seurat:SCTAssay'])

# Data
data_big = get_matrix(big_path, "data")
data_ste = get_matrix(stephenson_path, "data")
data_merged = concatenate_csr_matrices_by_row(data_big, data_ste)

sct_data = sct.create_group("data")
sct_data.create_dataset("data", data = data_merged.data)
sct_data.create_dataset("indices", data = data_merged.indices)
sct_data.create_dataset("indptr", data = data_merged.indptr)
sct_data.attrs.create("dims", data=[11205, data_merged.shape[0]])

# Counts
counts_big = get_matrix(big_path, "counts")
counts_ste = get_matrix(stephenson_path, "counts")
counts_merged = concatenate_csr_matrices_by_row(counts_big, counts_ste)

sct_counts = sct.create_group("counts")
sct_counts.create_dataset("data", data = counts_merged.data)
sct_counts.create_dataset("indices", data = counts_merged.indices)
sct_counts.create_dataset("indptr", data = counts_merged.indptr)
sct_counts.attrs.create("dims", data=[11205, counts_merged.shape[0]])

# now, we can set "SCT" as the active assay
hfile.attrs.create('active.assay', ['SCT'])

hfile.close()

