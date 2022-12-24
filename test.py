import ctypes
import numpy

lib = ctypes.cdll.LoadLibrary('./shared_object/C/libmyfusedmm_shared.so')


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
b = (ctypes.c_float * szB)(*[10, 11, 20, 21, 30, 31])   # B Matrix, Dense
c = (ctypes.c_float * szC)(*[0 for i in range(szC)])    # C Matrix, initally empty, to store result matrix

lib.SpMM(M, N, K, alpha, S_nnz, S_rows, S_cols, S_values, S_colids, S_rowptr, S_rowptr_plus_1, a, lda, b, ldb, beta, c, ldc)

print("Result Matrix C: \n")
for i in range(szC):
    print(c[i], end='\t')
    if ((i + 1) % M == 0):
        print("")

# S_rowptr = numpy.array([0, 3, 6], dtype=numpy.int64).ctypes.data_as(ctypes.POINTER(ctypes.c_longlong))
# S_rowptr1 = numpy.array([3, 6], dtype=numpy.int64).ctypes.data_as(ctypes.POINTER(ctypes.c_longlong))
# S_colids = numpy.array([0, 1, 2, 0, 1, 2], dtype=numpy.int64).ctypes.data_as(ctypes.POINTER(ctypes.c_longlong))

#xx
# S_values = numpy.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], dtype=numpy.float16).ctypes.data_as(ctypes.POINTER(ctypes.c_float))
# S_values = numpy.array([2 for _ in range(6)], dtype=numpy.float16).ctypes.data_as(ctypes.POINTER(ctypes.c_float))
# xx

# S_values = numpy.array([0 for i in range(6)], dtype=numpy.float16).ctypes.data_as(ctypes.POINTER(ctypes.c_float))
# for i in range(6):
#     S_values[i]  = i + 1


# a = numpy.array([0 for i in range(szA)], dtype=numpy.float16).ctypes.data_as(ctypes.POINTER(ctypes.c_float))
# b = numpy.array([0 for i in range(szB)], dtype=numpy.float16).ctypes.data_as(ctypes.POINTER(ctypes.c_float))
# c = numpy.array([0 for i in range(szC)], dtype=numpy.float16).ctypes.data_as(ctypes.POINTER(ctypes.c_float))

# b[0] = 10.0
# b[1] = 11.0
# b[2] = 20.0
# b[3] = 21.0
# b[4] = 30.0
# b[5] = 31.0


# for i in range(6):
#     print(S_values[i])

