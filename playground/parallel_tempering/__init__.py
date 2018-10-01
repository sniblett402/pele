# sn402: If you don't have access to an MPI library, these two lines can be removed.
# MPI is NOT required to use bhpt, and is not required to reproduce the results in this dataset.
from _base_mpi_ptmc import _MPI_Parallel_Tempering
from mpi_ptmc import *
