from numpy import *
from scipy.integrate import trapz, cumtrapz
import matplotlib
matplotlib.use("PDF")  # non-interactive plot making
import matplotlib.pyplot as plt
import os

# ----------------------------------------------------------------------------
#
#    begin main
#
# ----------------------------------------------------------------------------
#
set_printoptions(threshold='nan')
#
directory = './'

#
#        
#
eta = 0.5
NA = (eta/2.0)*(1.0+eta)
NB = 1.0-eta*eta
NC = -eta*0.5*(1.0-eta)


print ".. NA: ", NA
print ".. NB: ", NB
print ".. NC: ", NC
print "... sum: ", NA+NB+NC
#
s= 0.0
t= 0.0
N1= 0.25*(1.0-s)*(1.0-t)*(-s-t-1.0)
N2= 0.25*(1.0+s)*(1.0-t)*(s-t-1.0)
N3= 0.25*(1.0+s)*(1.0+t)*(s+t-1.0)
N4= 0.25*(1.0-s)*(1.0+t)*(-s+t-1.0)

N5= 0.5*(1.0-t)*(1.0+s)*(1.0-s)
N6= 0.5*(1.0-t)*(1.0+t)*(1.0+s)
N7= 0.5*(1.0+t)*(1.0+s)*(1.0-s)
N8= 0.5*(1.0-s)*(1.0+t)*(1.0-t)

print ".. N1: ", N1
print ".. N2: ", N2
print ".. N3: ", N3
print ".. N4: ", N4
print ".. N5: ", N5
print ".. N6: ", N6
print ".. N7: ", N7
print ".. N8: ", N8

print "... sum: ", N1+N2+N3+N4+N5+N6+N7+N8


exit(0)
