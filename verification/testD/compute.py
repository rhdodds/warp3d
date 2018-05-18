import math
import numpy as np
from scipy.integrate import trapz, cumtrapz
import matplotlib
matplotlib.use("PDF")  # non-interactive plot making
import matplotlib.pyplot as plt
import os

ym=30000.
nu=0.3
alpha=10.0e-06
theta_edge= 0.6*166.667

sig_theta_e = ym * alpha * theta_edge / (1.0-nu)

print "sig_theta_e: ", sig_theta_e


J=0.2587E-01

print "Norm J: ", J/(sig_theta_e*sig_theta_e*0.6/ym)

KI=(ym*J/.91)**0.5

print "Norm KI: ", KI/sig_theta_e/(3.14159*0.6)**0.5
