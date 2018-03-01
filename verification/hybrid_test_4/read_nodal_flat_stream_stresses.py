from numpy import *
import os

#
set_printoptions(threshold='nan')
#
stream_fname = "./wns0000030_stream"
#
fileobj = open(stream_fname, mode='rb')
num_cols = 26
num_nodes = 1568
num_vals_to_read = num_nodes * num_cols
data = fromfile(file=fileobj,count=num_vals_to_read,
       dtype=float64).reshape(num_nodes, num_cols)
x0 = data[:,0]
x1 = data[:,1]
x2 = data[:,2]
x3 = data[:,3]
x4 = data[:,4]
x5 = data[:,5]
x6 = data[:,6]
x21 = data[:,21]
print  " %15.7e  %15.7e  %15.7e"%(
     x0[num_nodes-1] ,x1[num_nodes-1], x2[num_nodes-1])

exit(0)

