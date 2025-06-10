from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

my_values = np.zeros(3)

for i in range (3):
    if rank == 0:
        my_values[i] = i +1 

print('Prima: ', my_values[0], ' ', my_values[1], ' ', my_values[2], 'per il processo: ', rank, '\n')

comm.Bcast(my_values, root=0)

print('Dopo: ', my_values[0], ' ', my_values[1], ' ', my_values[2], 'per il processo: ', rank, '\n')
