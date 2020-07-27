c
      module patwarp_data
c
      parameter  ( maxprp = 2000 )
      parameter  ( maxmtl = 2000 )
      parameter  ( maxset = 200 )
      parameter  ( maxcon = 100 )
      parameter  ( maxtrn = 2000 )
      parameter  ( maxetp = 5000 )
      parameter  ( mxnmgp = 30000 )
      parameter  ( mxnmbl = 2*(8000000/128) )
      parameter  ( mxblsz = 512 )
      parameter  ( maxprocs = 10000 )
      parameter  ( maxaccgrps = 10000 )
c
      integer  numnod, numele, nummtl, numprp, intact, maxinc,
     &         maxadv, maxnls, maxndt, maxnvl, maxeld, maxdld,
     &         termin, termot, nfile, ofile, cd_file, inbl_file,
     &         ct_file, ofile_store, comdummy, blksz, nmpc, maxmpc,
     &         cdf_len, inbl_len, ctf_len, nxtels, numgrp, nelblk,
     &         two16, nxtnls, nxtndt, msgnod, msgele, nxtnvl,
     &         eincpt, eadvpt, nxteld, nxtdld, block_method

      integer  ncset(6), ynanswers(9), nlodtb(maxset,5),
     &         econfg(maxetp), ettypt(maxetp), etty_tail_pt(maxetp),
     &         elblks(0:mxblsz,mxnmbl), blkproc(mxnmbl),
     &         grpnum(mxnmgp), grphead(mxnmgp), grptail(mxnmgp)
c
      integer, allocatable, dimension (:) :: nodcon, nodndf, nodcfg,
     &                                       nodcid, nodpcn,
     &                                       nod_trn_list,
     &                                       elepid, elecid, eletyp,
     &                                       elennd, elecfg, eleipt,
     &                                       elenad, eleadp, eletran,
     &                                       eletrpoin, elreno, eleord,
     &                                       elegnl, elem_access,
     &                                       elreon, eleproc, grplst,
     &                                       eleinc, eletrinc,
     &                                       ndcndt, nodval, nodptr,
     &                                       mpctrm
c
      integer, allocatable, dimension (:,:) :: cset, etypls, ndcnls,
     &                                         eloads, deload,
     &                                         mpcnod, mpcdof
c
      real, allocatable, dimension (:) :: eleadv, rndcndt, rnodval
c
      real, allocatable, dimension (:,:) :: eleang
c
      double precision, allocatable, dimension (:,:) :: coord, mpcmlt
c
      character*8   time
      character*12  date, versn
      character*80  neunam, outnam, ordnam, coord_file, const_file,
     &              incid_block_file, ustitl, elem_print_file,
     &              nptfle, resfile, prefix_name
c
      character*16  etypes(maxetp)
c
      character*1, allocatable, dimension (:) :: nodtyp
c
      logical  debug, prnode, princi, prcons, prnlod, prntmp, prnmpc,
     &         prcset, scalar, prdistl, parallel, sep_file_control,
     &         made_transition
c
c
      end module patwarp_data
