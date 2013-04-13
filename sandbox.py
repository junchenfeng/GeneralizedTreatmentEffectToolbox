
#import  cProfile
#import sampleSimulation as ss
#import localLinearRegression as LLR

#Y,D,X,Z = ss.underiv_data(100)


#cProfile.run('LLR.optimalBandwidthSelection(Y,X[:,0:1])')

#test1 = LLR.LOORiskEstFast(Y, X[:, 0:1], 1)
#test2 = LLR.LOORiskEstSlow(Y, X[:, 0:1], 1)

#print test1, test2

#from mpi4py import MPI
#import numpy

#comm = MPI.COMM_WORLD
#rank = comm.Get_rank()

#if rank == 0:
#    data = numpy.arange(1000, dtype = 'i')
#    comm.Send([data, MPI.INT], dest = 1, tag = 77)
#elif rank == 1:
#    data = numpy.empty(1000, dtype = 'i')
#    comm.Recv([data, MPI.INT], source = 0, tag = 77)


import sampleSimulation as ss


import localLinearRegression as LLR
import numpy
Y,D,X,Z = ss.underiv_data(500)
h = numpy.arange(0.05,2,0.05)

LLR.optimalBandwidthSelection(Y, X[:,0:1], h)


# Retrieves the result calculated by job1
# The value of job1() is the same as sum_primes(100)
# If the job has not been finished yet, execution will wait here until result is available
