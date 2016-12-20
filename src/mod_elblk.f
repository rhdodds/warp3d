c     ****************************************************************
c     *                                                              *
c     *                     f-90 data modules:                       *
c     *                        elblk_data                            *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 03/10/04                   *
c     *                                                              *
c     *     define the data structures for a block of elements       *
c     *     used in various places. this replaces common.elblk       *
c     *                                                              *
c     ****************************************************************
c
c
c
      module elblk_data
c     -----------------
c
      include 'param_def'
c
c                stresses, strains, histories, etc. for a block
c
      double precision
     & ddtse(mxvl,nstr,mxgp),
     & urcs_blk_n(mxvl,nstrs,mxgp),
     & urcs_blk_n1(mxvl,nstrs,mxgp), rot_blk_n1(mxvl,9,mxgp)
      save rtse, ddtse, urcs_blk_n, urcs_blk_n1, rot_blk_n1
c
      double precision,
     &       dimension(:,:,:), save, allocatable ::
     &        elem_hist, elem_hist1
c
c                internal forces, transformations for block of elements
c
      double precision
     & elestr(mxvl,mxoupr,mxoupt)
      save elestr
c
      integer, save :: blk_size_hist, blk_size_gp
c
      intrinsic allocated
c
      end module
