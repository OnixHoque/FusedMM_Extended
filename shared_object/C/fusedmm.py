import ctypes
import os

dll = ctypes.CDLL('./libmyfusedmm_shared_kernels.so', mode=os.RTLD_LAZY)
# lib = ctypes.cdll.LoadLibrary('./libmyfusedmm_shared.so', mode=1)
