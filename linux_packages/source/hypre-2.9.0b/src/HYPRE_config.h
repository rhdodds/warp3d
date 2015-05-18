/* HYPRE_config.h.  Generated from HYPRE_config.h.in by configure.  */
/*BHEADER**********************************************************************
 * Copyright (c) 2008,  Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * This file is part of HYPRE.  See file COPYRIGHT for details.
 *
 * HYPRE is free software; you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * $Revision: 1.21 $
 ***********************************************************************EHEADER*/


/* config/HYPRE_config.h.in.  Generated from configure.in by autoheader.  */

/* Release name */
#define HYPRE_RELEASE_NAME "hypre"

/* Version number */
#define HYPRE_RELEASE_VERSION "2.9.0b"

/* Date of release */
#define HYPRE_RELEASE_DATE "2012/10/30"

/* Time of release */
#define HYPRE_RELEASE_TIME "00:00:00"

/* Bug reports */
#define HYPRE_RELEASE_BUGS "hypre-support@llnl.gov"

/* Define to 1 for Solaris. */
/* #undef HYPRE_SOLARIS */

/* Define to 1 for Linux on platforms running any version of CHAOS */
/* #undef HYPRE_LINUX_CHAOS */

/* Define to 1 for Linux platforms */
#define HYPRE_LINUX 1

/* Define to 1 for Alpha platforms */
/* #undef HYPRE_ALPHA */

/* Define to 1 for RS6000 platforms */
/* #undef HYPRE_RS6000 */

/* Define to 1 for IRIX64 platforms */
/* #undef HYPRE_IRIX64 */

/* Define to 1 if using long long int for HYPRE_Int */
/* #undef HYPRE_BIGINT */

/* Define to 1 if an MPI library is found */
#define HYPRE_HAVE_MPI 1

/* Define to 1 if the routine MPI_Comm_f2c is found */
/* #undef HYPRE_HAVE_MPI_COMM_F2C */

/* Disable MPI, enable serial codes */
/* #undef HYPRE_SEQUENTIAL */

/* Using HYPRE timing routines */
/* #undef HYPRE_TIMING */

/* Using dxml for BLAS */
/* #undef HYPRE_USING_DXML */

/* Using essl for BLAS */
/* #undef HYPRE_USING_ESSL */

/* Using internal Hypre routines */
/* #undef HYPRE_USING_HYPRE_BLAS */

/* Using internal Hypre routines */
/* #undef HYPRE_USING_HYPRE_LAPACK */

/* No global partitioning being used */
/* #undef HYPRE_NO_GLOBAL_PARTITION */

/* Print HYPRE errors */
/* #undef HYPRE_PRINT_ERRORS */

/* Enable OpenMP support */
#define HYPRE_USING_OPENMP 1

/* Define as follows to set the Fortran name mangling scheme:
 * 0 = unspecified
 * 1 = no underscores
 * 2 = one underscore
 * 3 = two underscores
 * 4 = caps, no underscores */
#define HYPRE_FMANGLE 0

/* Define as in HYPRE_FMANGLE to set the BLAS name mangling scheme */
#define HYPRE_FMANGLE_BLAS 0

/* Define as in HYPRE_FMANGLE to set the LAPACK name mangling scheme */
#define HYPRE_FMANGLE_LAPACK 0

/* Define to a macro mangling the given C identifier (in lower and upper
 * case), which must not contain underscores, for linking with Fortran. */
#define HYPRE_F77_FUNC(name,NAME) name ## _

/* As HYPRE_F77_FUNC, but for C identifiers containing underscores. */
#define HYPRE_F77_FUNC_(name,NAME) name ## _
