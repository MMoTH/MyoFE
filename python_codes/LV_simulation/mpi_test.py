from mpi4py import MPI

def square(n):
    return n*n
comm = MPI.COMM_WORLD

rank_id = comm.Get_rank()

size = comm.Get_size()
print('Rank %0.0f is called'%rank_id)

if rank_id == 0:
    a=2
    for i in range(1,size):
        a_2 = square(a)
        comm.send(a_2,dest=i)

        print('Rank 0 has sent a_2 to rank %0.0f'%i)
else:
    a_2 = comm.recv(source = 0)
    print('Rank %0.0f has recived a_2 from rank 0'%rank_id)


def square(n):
    return n*n