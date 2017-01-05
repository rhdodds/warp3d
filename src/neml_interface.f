      interface
            function create_nemlmodel(fname, mname, ier) bind(C)
                  use iso_c_binding
                  implicit none
                  type(c_ptr) :: create_nemlmodel
                  character(kind=c_char) :: fname(*)
                  character(kind=c_char) :: mname(*)
                  integer :: ier
            end function

            subroutine destroy_nemlmodel(model, ier) bind(C)
                  use iso_c_binding
                  implicit none
                  type(c_ptr), value :: model
                  integer :: ier
            end subroutine

            function nstore_nemlmodel(model) bind(C)
                  use iso_c_binding
                  implicit none
                  integer :: nstore_nemlmodel
                  type(c_ptr), value :: model
            end function

            subroutine init_store_nemlmodel(model, store, ier) bind(C)
                  use iso_c_binding
                  implicit none
                  type(c_ptr), value :: model
                  double precision, intent(out), dimension(*) :: store
                  integer, intent(out) :: ier
            end subroutine

            subroutine update_sd_nemlmodel(model, e_np1, e_n, Temp_np1,
     &                  Temp_n, time_np1, time_n, s_np1, s_n,
     &                  h_np1, h_n,
     &                  A_np1, ier) bind(C)
                  use iso_c_binding
                  implicit none
                  type(c_ptr), value :: model
                  
                  double precision, intent(in), dimension(6) ::
     &                  e_np1, e_n, s_n
                  double precision, intent(out), dimension(6) ::
     &                  s_np1
                 double precision, intent(out), dimension(6,6) ::
     &                  A_np1
                  double precision, intent(in), dimension(*) ::
     &                  h_n
                  double precision, intent(out), dimension(*) ::
     &                  h_np1
                  double precision, intent(in), value ::
     &                  Temp_np1, Temp_n, time_np1, time_n
                  integer, intent(out) :: ier

            end subroutine


      end interface
