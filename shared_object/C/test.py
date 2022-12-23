import ctypes
import numpy

lib = ctypes.cdll.LoadLibrary('./libmyfusedmm_shared.so')

lib.performDummySpMM()


# M = 2
# N = 3
# K = 2

# alpha = 1
# beta = 0
# tkern = 'm'

# S_rows = 2
# S_cols = 2
# S_nnz = 6

# S_rowptr = numpy.array([0, 3, 6], dtype=numpy.int64)
# S_colids = numpy.array([0, 1, 2, 0, 1, 2], dtype=numpy.int64)
# S_values = numpy.array([1, 2, 3, 4, 5, 6], dtype=numpy.float16)

# lda = ldb = ldc = K

# szA = szB = szC = 16

# a = numpy.array([0 for i in range(szA)], dtype=numpy.float16)
# b = numpy.array([0 for i in range(szB)], dtype=numpy.float16)
# c = numpy.array([0 for i in range(szC)], dtype=numpy.float16)

# b[0] = 10
# b[1] = 11
# b[2] = 20
# b[3] = 21
# b[4] = 30
# b[5] = 31


# lib.SpMM(M, N, K, alpha, S_nnz, S_rows, S_cols, S_values.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), S_colids.ctypes.data_as(ctypes.POINTER(ctypes.c_long)), S_rowptr.ctypes.data_as(ctypes.POINTER(ctypes.c_long)), S_rowptr.ctypes.data_as(ctypes.POINTER(ctypes.c_long)), a.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), lda, b.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), ldb, beta, c.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), ldc)

# SpMM(M, N, K, alpha, S_nnz, S_rows, S_cols, S_values, 
# S_colids, S_rowptr, S_rowptr+1, a, lda, b, ldb, 
# beta, c, ldc);   

