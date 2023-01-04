import ctypes
import torch

lib = ctypes.cdll.LoadLibrary('./libmyfusedmm_shared.so')
lib.SpMM.argtypes = [ctypes.c_longlong, ctypes.c_longlong, ctypes.c_longlong, ctypes.c_float, ctypes.c_longlong, ctypes.c_longlong, ctypes.c_longlong, ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_longlong), ctypes.POINTER(ctypes.c_longlong), ctypes.POINTER(ctypes.c_longlong), ctypes.POINTER(ctypes.c_float), ctypes.c_longlong, ctypes.POINTER(ctypes.c_float), ctypes.c_longlong, ctypes.c_float, ctypes.POINTER(ctypes.c_float), ctypes.c_longlong] 
lib.SpMM.restype = ctypes.c_void_p

def performDummySpMM():
    # lib.performDummySpMM()

    # Peforms this mm: https://miro.medium.com/max/1400/1*YGcMQSr0ge_DGn96WnEkZw.png

    # SpMM Function=============
    # Performs MM: S X B = C
    # S Dim = M x N, B Dim = N X K, C Dim = M x K

    lib.SpMM.argtypes = [ctypes.c_longlong, ctypes.c_longlong, ctypes.c_longlong, ctypes.c_float, ctypes.c_longlong, ctypes.c_longlong, ctypes.c_longlong, ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_longlong), ctypes.POINTER(ctypes.c_longlong), ctypes.POINTER(ctypes.c_longlong), ctypes.POINTER(ctypes.c_float), ctypes.c_longlong, ctypes.POINTER(ctypes.c_float), ctypes.c_longlong, ctypes.c_float, ctypes.POINTER(ctypes.c_float), ctypes.c_longlong] 
    lib.SpMM.restype = ctypes.c_void_p


    alpha = 1; beta = 0; tkern = 'm';  # Don't change

    # Matrix Dimensions: 
    M = 2; N = 3; K = 2

    # Sparse Matrix properties
    S_rows = 2 ; S_cols = 2; S_nnz = 6

    S_rowidx = [0, 3, 6]
    S_rowptr = (ctypes.c_longlong * 3)(*S_rowidx)
    S_rowptr_plus_1 = (ctypes.c_longlong * 2)(*S_rowidx[1:])
    S_colids = (ctypes.c_longlong * 6)(*[0, 1, 2, 0, 1, 2])
    S_values = (ctypes.c_float * 6)(*[1, 2, 3, 4, 5, 6])


    lda = ldb = ldc = K         # Note: ***ATL_Cachelen is not integrated yet***

    szA = 1     # Unused for SpMM
    szB = N * K # Elements of B dense matrix
    szC = M * K # Elements of resultant C dense matrix



    a = (ctypes.c_float * szA)(*[0 for i in range(szA)])    # Unused
    b_flat = torch.flatten(b_dense)
    b = ctypes.cast(b_flat.data_ptr(), ctypes.POINTER(ctypes.c_float))
    # print(b[0])
    # b = (ctypes.c_float * szB)(*[10, 11, 20, 21, 30, 31])   # B Matrix, Dense
    c = (ctypes.c_float * szC)(*[0 for i in range(szC)])    # C Matrix, initally empty, to store result matrix

    lib.SpMM(M, N, K, alpha, S_nnz, S_rows, S_cols, S_values, S_colids, S_rowptr, S_rowptr_plus_1, a, lda, b, ldb, beta, c, ldc)

    print("Result Matrix C: \n")
    for i in range(szC):
        print(c[i], end='\t')
        if ((i + 1) % M == 0):
            print("")

def csr_spmm_cpu(A_rowptr, A_colid, A_vals, B_dense):
    alpha = 1; beta = 0; tkern = 'm';  # Don't change

    # Matrix Dimensions: 
    M = A_rowptr.size(0) -1 
    N = B_dense.size(0)
    K = B_dense.size(1)

    # Sparse Matrix properties
    S_rows = M
    S_cols = N
    S_nnz = A_vals.size(0)

    # S_rowidx = [0, 3, 6]
    S_rowptr = ctypes.cast(A_rowptr.data_ptr(), ctypes.POINTER(ctypes.c_longlong))
    S_rowptr_plus_1 = ctypes.cast(A_rowptr[1:].data_ptr(), ctypes.POINTER(ctypes.c_longlong))
    S_colids = ctypes.cast(A_colid.data_ptr(), ctypes.POINTER(ctypes.c_longlong))
    S_values = ctypes.cast(A_vals.data_ptr(), ctypes.POINTER(ctypes.c_float))

    lda = ldb = ldc = K         # Note: ***ATL_Cachelen is not integrated yet***

    szA = 1     # Unused for SpMM
    szB = N * K # Elements of B dense matrix
    szC = M * K # Elements of resultant C dense matrix

    a_tensor = torch.zeros(1)
    a = ctypes.cast(a_tensor.data_ptr(), ctypes.POINTER(ctypes.c_float))
    
    # a = (ctypes.c_float * szA)(*[0 for i in range(szA)])    # Unused
    # b_tensor = torch.tensor([10, 20, 20, 21, 30, 31], dtype=torch.float32)
    b_flat = torch.flatten(B_dense)
    b = ctypes.cast(b_flat.data_ptr(), ctypes.POINTER(ctypes.c_float))
    # b = ctypes.cast(B_dense.data_ptr(), ctypes.POINTER(ctypes.c_float))
    # print(b_flat)
    # b = (ctypes.c_float * szB)(*[10, 11, 20, 21, 30, 31])   # B Matrix, Dense
    
    c_tensor = torch.zeros(szC)
    c = ctypes.cast(c_tensor.data_ptr(), ctypes.POINTER(ctypes.c_float))    # C Matrix, initally empty, to store result matrix

    lib.SpMM(M, N, K, alpha, S_nnz, S_rows, S_cols, S_values, S_colids, S_rowptr, S_rowptr_plus_1, a, lda, b, ldb, beta, c, ldc)

    return c_tensor.reshape(M, K)

a_rowptr = torch.tensor([0, 3, 6], dtype=torch.int64)
a_colid = torch.tensor([0, 1, 2, 0, 1, 2], dtype=torch.int64)
a_vals = torch.tensor([1, 2, 3, 4, 5, 6], dtype=torch.float32)
b_dense = torch.tensor([[10, 11], [20, 21], [30, 31]], dtype=torch.float32)

c = csr_spmm_cpu(a_rowptr, a_colid, a_vals, b_dense)

print(c)

# performDummySpMM()

# a_rowptr = torch.tensor([])
# a_col = torch.tensor([])
# a_vals = torch.tensor([])
# b_dense = torch.tensor([10, 11, 20, 21, 30, 31], dtype=torch.float32)

# print(b_dense)

# performDummySpMM()
# csr_spmm_cpu(A_rowptr : torch.tensor, A_colind : torch.tensor, A_csrVal : torch.tensor, B : torch.tensor):
#     pass


# ptr = torch.flatten(b_dense).data_ptr()
# print(ptr[0])
# import numpy as np
# a = np.array([10, 20, 20, 21, 30, 31],dtype=np.float32)
# b = torch.tensor(a , dtype=torch.float32)
# p = ctypes.cast(b.data_ptr() , ctypes.POINTER(ctypes.c_float))
# print(p+1)

# # np.set_printoptions(precision=10)
# a = np.array([10, 20, 30],dtype=np.float32)

# t = torch.tensor(a , dtype=torch.float32)
# ptr = ctypes.cast(t.data_ptr(), ctypes.POINTER(ctypes.c_float))
# print(f'{a}\n{t}\n{ptr[0]}')
# # print(np.ndarray(1, dtype=np.uint32,buffer=a ))

