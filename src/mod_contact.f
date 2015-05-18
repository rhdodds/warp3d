c     ****************************************************************
c     *                                                              *
c     *                    f-90 module contact                       *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 5/10/04                    *
c     *                                                              *
c     *     define the variables and data structures to support      *
c     *     rigid body contact                                       *
c     *                                                              *
c     ****************************************************************
c
c
c
      module contact
c
      parameter ( maxcontact = 20, maxcntprm = 10 )
c
      integer, save :: num_contact
      integer, save, allocatable :: contact_cause(:,:)
      integer, save :: contact_shape(maxcontact)
c
      logical, save :: use_contact
c     
#dbl      double precision, save ::
#sgl      real, save ::
     &	   cplane_vec (3,2,maxcontact),
     &     cshape_norm (3,maxcontact),
     &     cshape_pnt (3,maxcontact),
     &     cshape_rate (3,maxcontact),
     &     cshape_param (maxcntprm,maxcontact),
     &     contact_depth (maxcontact),
     &     contact_stiff (maxcontact), 
     &     contact_fric (maxcontact)
c
#dbl      double precision, allocatable, save ::
#sgl      real, allocatable, save ::
     &     contact_force(:)
c
      end module


