import ctypes
import torch
from scipy.sparse import coo_matrix
import os
import numpy as np

fusedmm_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'libmyfusedmm_shared.so')
# print(fusedmm_path)
# lib = ctypes.cdll.LoadLibrary('./libmyfusedmm_shared.so')
lib = ctypes.cdll.LoadLibrary(fusedmm_path)

lib.SpMM.argtypes = [ctypes.c_longlong, ctypes.c_longlong, ctypes.c_longlong, ctypes.c_float, ctypes.c_longlong, ctypes.c_longlong, ctypes.c_longlong, ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_longlong), ctypes.POINTER(ctypes.c_longlong), ctypes.POINTER(ctypes.c_longlong), ctypes.POINTER(ctypes.c_float), ctypes.c_longlong, ctypes.POINTER(ctypes.c_float), ctypes.c_longlong, ctypes.c_float, ctypes.POINTER(ctypes.c_float), ctypes.c_longlong] 
lib.SpMM.restype = ctypes.c_void_p
lib.mytest_csr.argtypes = [ctypes.c_char, ctypes.c_longlong, ctypes.c_longlong, ctypes.c_longlong, ctypes.c_float, ctypes.c_longlong, ctypes.c_longlong, ctypes.c_longlong, ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_longlong), ctypes.POINTER(ctypes.c_longlong), ctypes.POINTER(ctypes.c_longlong), ctypes.POINTER(ctypes.c_float), ctypes.c_longlong, ctypes.POINTER(ctypes.c_float), ctypes.c_longlong, ctypes.c_float, ctypes.POINTER(ctypes.c_float), ctypes.c_longlong] 
lib.mytest_csr.restype = ctypes.c_void_p
kernel_type = {
    't-dist': ord('t'),
    'fr': ord('f'),
    'sigmoid': ord('s'),
    'spmm': ord('m'),
    'gcn': ord('g')
}

def spmm_fusedmm_v1(A_rowptr, A_colid, A_vals, B_dense):
    alpha = 1; beta = 0; tkern = 'm';  # Don't change

    # Matrix Dimensions: 
    M = A_rowptr.size(0) -1 
    N = B_dense.size(0)
    K = B_dense.size(1)

    # Sparse Matrix properties
    S_rows = M
    S_cols = N
    S_nnz = A_vals.size(0)

    
    S_rowptr = ctypes.cast(A_rowptr.data_ptr(), ctypes.POINTER(ctypes.c_longlong))
    S_rowptr_plus_1 = ctypes.cast(A_rowptr[1:].data_ptr(), ctypes.POINTER(ctypes.c_longlong))
    S_colids = ctypes.cast(A_colid.data_ptr(), ctypes.POINTER(ctypes.c_longlong))
    S_values = ctypes.cast(A_vals.data_ptr(), ctypes.POINTER(ctypes.c_float))

    lda = ldb = ldc = K         # Note: ***ATL_Cachelen is not integrated yet***

    #Sp = M x N
    szA = 1     # Unused for SpMM
    szB = N * K # Elements of B dense matrix
    szC = M * K # Elements of resultant C dense matrix

    a_tensor = torch.zeros(1)
    a = ctypes.cast(a_tensor.data_ptr(), ctypes.POINTER(ctypes.c_float))
    
    b_flat = torch.flatten(B_dense)
    b = ctypes.cast(b_flat.data_ptr(), ctypes.POINTER(ctypes.c_float))
    
    c_tensor = torch.zeros(szC)
    c = ctypes.cast(c_tensor.data_ptr(), ctypes.POINTER(ctypes.c_float))    # C Matrix, initally empty, to store result matrix

    # lib.SpMM(M, N, K, alpha, S_nnz, S_rows, S_cols, S_values, S_colids, S_rowptr, S_rowptr_plus_1, a, lda, b, ldb, beta, c, ldc)
    lib.mytest_csr(kernel_type['spmm'], M, N, K, alpha, S_nnz, S_rows, S_cols, S_values, S_colids, S_rowptr, S_rowptr_plus_1, a, lda, b, ldb, beta, c, ldc)
    
    return c_tensor.reshape(M, K)

def spmm_fusedmm_v2(idx, val, m, n, matrix):
    csr = coo_matrix((val, idx), shape=(m,n), copy=False).tocsr(False)
    return spmm_fusedmm_v1(torch.as_tensor(csr.indptr, dtype=torch.int64), torch.as_tensor(csr.indices, dtype=torch.int64), val, matrix)


# # v-1 test: (CogDL Interface)
##-----------------------------------
# a_rowptr = torch.tensor([0, 3, 6], dtype=torch.int64)
# a_colid = torch.tensor([0, 1, 2, 0, 1, 2], dtype=torch.int64)
# a_vals = torch.tensor([1, 2, 3, 4, 5, 6], dtype=torch.float32)
# b_dense = torch.tensor([[10, 11], [20, 21], [30, 31]], dtype=torch.float32)
# c = spmm_fusedmm_v1(a_rowptr, a_colid, a_vals, b_dense)
# print(c)

## v-2 test: (Torch Sparse Interface)
##-----------------------------------
# index = torch.tensor([[0, 0, 1, 2, 2],
#                       [0, 2, 1, 0, 1]])
# value = torch.Tensor([1, 2, 4, 1, 3])
# matrix = torch.Tensor([[1, 4], [2, 5], [3, 6]])
# c = spmm_fusedmm_v2(index, value, 3, 3, matrix)
# print(c)

