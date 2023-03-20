
March 20, 2023

R. Dodds

Python3 program to show impact of using the PCHIP method 
(Piecewise Cubic Hermite Interpolation Polynomials) to convert
segmental stress-strain curves into curves with continuous
plastic tangent modulus (H') at each point on the curve.

The method was derived originally for general 1-D data sets where
conventional spline fits lead to large over/under interpolations
in data sets.

Wikipedia article on Monotone Cubic Interpolation included in
directory as wikipedia_article.pdf

Original paper in SIAM J. Sci. Stat. Computing included as pchip_paper.pdf.
That algorithm is used in this Python code and in WARP3D. 
Fritsch-Butland with Brodlie modification for non-uniformly
spaced strain values. 

The derivative for point 1 is computed from the one-sided, weighted 
harmonic mean value suggested in C. Moler’s 2004 book 
Numerical Computing with MATLAB and also used in SciPy.

This could also be used for the last point but not yet included in this
code.

The Python code here and SciPy now match.


This Code:
==========

This is a Python3 program. Imports standard packages: numpy, scipy
It also uses the supporting Python3 file Dodds_plot_support.py included
in the directory.

Output:
=======

Plot: stress_strain_curve.pdf shows
stress vs. plastic strain for (1) data file points, 
(2) the WARP3D implementation of PCHIP, (3) the SciPY 
implementation of PCHIP.

Plot: hprime.pdf shows
the slope of stress vs. plastic strain curve generally referred
to as H' (hprime). Derived from WARP3D pchip curve and the
SciPy pchip curve. These 2 curves should nearly always 
be very much identical.

Usage:
======

Create one or more (text) input files with a segmental stress-strain
curve. These are total stress vs. total strain values. The file format is:

<Young's modulus>
<number of data points>
0.0 0.0
<strain value> <stress value>
...
...

Example stress-strain curves are included in this directory.
The code is agnostic about units. Just use consistent units for the
modulus and stress values.

The file "wiggly_curve" illustrates how the method handles curves with
dips.

Execution: (assume $ is command line prompt in terminal)
=========

The code prompts for input.

$python3 driver.py    ==>  assumes python3 has NumPy and SciPy pkgs. 

> file w/ curve points (use quit to stop, CR for same file): n10_EPRI_expanded_curve
> min,max stress; min,max plastic strain for plot (4 values, no commas): 0 100 0 1
> min,max range for hprime in hprime vs. plastic strain plot (2 values): 0 5000
.... last curve values from file:  4.0042769 128.3082
.... Building pchip slopes to use at pts on user curve
.... Building WARP3D pchip interpolated points to plot
.... Building SciPy pchip interpolated points to plot
.... Starting plot stress-strain
.... PDF for plot: stress_strain_curve.pdf
.... Starting plot hprime
.... PDF for plot: hprime.pdf

The input loop continues to simplify exploration of different X-Y axis
plot limits. 


The WARP3D version and the SciPy curves in both plots should have
v. small differences.

** Advice:
==========

The PCHIP algorithm of the input segmental curve may be deceivingly
smooth. The PCHIP curve is always C0 and C1 continuous.

Check the hprime.pdf plot. This is slope of the stress vs. plastic
strain curve fit. It is C0 continuous but not C1 continuous. Quite small
perturbations of the segmental points may trigger small humps increasing
H' when a continuously decreasing H' is often desired. Points on the segmental 
curve may require tweaking/nudging to achieve no "surprises" in the H' values.





