# warp3d
Open Source Static and Dynamic Nonlinear Analysis of Solids:  Windows, Linux and Mac OS X Computers

WARP3D is an open-source, production-quality code under continuing development 
to meet the challenges of large-scale, 3-D solid simulations for focused 
investigations on fatigue and fracture mechanisms and behavior in metallic 
components/structures. WARP3D provides an alternative computational resource 
for this narrower class of simulations in comparison with the expansive, 
general-purpose commercial codes (e.g. Abaqus, Ansys) and evolving families 
of the non-public, multi-physics driven codes at the national labs (e.g. 
Sierra, Grizzly, Moose, Diablo, …). The WARP3D source code, ready-to-run 
executables, extensive technical and user-documentation, verification and 
example problem suites, post-processors and documented workflows are 
provided under the University of Illinois/NCSA, Open Source License
(allows free unrestricted use, modification, redistribution, commercialization).

Several unique features of the code that support challenging simulations to 
understand fatigue and fracture processes include:

•	Comprehensive library of solid and interface-cohesive elements all supporting 
    finite-deformation behavior 
    
•	An array of constitutive models for metals including temperature, strain-rate, 
    creep and finite-deformation effects; extensive models for crystal plasticity simulations of 
    fcc, bcc, bcc48, hcp6, hcp18, and a single-slip systems with multiple options for 
    slip-rate relationships; 3D nonlinear cohesive constitutive relationships with 
    mode I-II-III interactions, boundary cavitation/slip models, nonlocal effects 
    on cohesive behavior and finite interface separations-rotations. All implemented 
    in a robust finite deformation formulation based on decompositions of the 
    deformation gradients F.
    
•	Computation of the 3D J-integral including the combined effects of residual 
    strains/stresses, crack face load-ing, thermal loading, inertia, functionally 
    graded and anisotropic materials with arbitrary orientations. The newest capabilities 
    (fall 2018) support thermo-mechanical process simulations with extensive 
    plastic deformations (e.g. bead-by-bead weld simulation) prior to insertion 
    of one or more cracks – this unique and complete J-integral formulation 
    in the code provides path (domain) independent values under these severe 
    conditions. Interaction integrals to compute KI, KII , KIII  and T-stress
    components for linear elastic solutions.
•	Options to grow cracks during a 3D simulation include: (1) node release 
    via manual and/or automatic track-ing of CTOD/CTOA along a front, 
    (2) manual or automatic tracking of element damage/degradation and 
    extinction (Gurson-Tvergaard, SMCS, …), (3) 3D triangular/quadrilateral 
    interface elements with degrading, mixed-mode local/nonlocal cohesive 
    behavior. Crack growth drivers adaptively govern the global solution 
    processes to avoid truncation of highly local damage evolution along crack fronts. 
•	General 3D mesh-tying, rigid body contact, extensive library of model loading capabilities

•	An exceptionally robust, globally implicit-iterative solver with multi-level 
    adaptive sub stepping combined with line-search and extrapolation. 
    Employs the industry leading, high-performance Pardiso equation solver 
    (free in Intel’s MKL) allowing models with millions of nodes/elements 
    on single workstations and clusters, with LLNLs hypre solver also used on clusters.
    
•	Plug-compatible interfaces with Abaqus UMAT and  UEXTERNALDB to incorporate 
    existing, highly-specialized constitutive models and external data access.
    
•	Multiple output schemes including support for ParaView.

•	Web sites: warp3d.net (overview, all documentation, ready-to-run 
    executables, complete system in .zip package); github.com/rhdodds/warp3d (developmental
    repository with all latest updates)
    
•	Supported platforms: Windows 10/8/7 (parallel via OpenMP), Linux: Redhat, 
    Ubuntu, … (parallel via OpenMP threads and MPI+OpenMP for clusters), 
    Mac OSX (parallel via OpenMP). Linux and Mac OSX versions may be re-built 
    to include user modifications with gfortran or Intel Fortran.
    Rebuilds on Windows require Intel Fortran. No compilers are 
    required to download and use the ready-to-run executables on Windows, Linux and Mac OSX. 

WARP3D is supported on two web sites

1. See http://www.warp3d.net for general and more detailed information, 
full system downloads, ready-to-run executables, documentation and 
current notes/issues.

2. This web site on GitHub is used by system developers to update 
and maintain the software and documentation including all source 
code, development tools, verification suites.

Most WARP3D users will not need to use this web site (GitHub).
If you want to have the most up-to-date source code, documentation, 
etc. you may clone this repository. This will not include 
ready-to-run executables.

GitHub will download a zip file for the entire system or you 
can clone this repository.

For cloning, you will need to have the (free) git software for 
distributed source control management tool. git runs on all platforms.
To create a clone of the current development (but working) WARP3D 
system version, use the following command (in a bash shell on Linux and 
OSX). Choose any name you want for the directory to hold the clone (below uses warp3d)

% git clone https://github.com/rhdodds/warp3d  \<name for your copy of repository\>


Follow the step-by-step instructions to rebuild WARP3D 
executable(s) on a platform provided in the files 
README_for_xxxx (where xxxx is the platform).
