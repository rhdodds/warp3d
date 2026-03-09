      module mm10_work_def
c      
      use iso_Fortran_env
      use mm10_defs, only : crystal_props, crystal_state
      implicit none
      include 'param_def'
c
c              declaration of local arrays for mm10_solver
c
      type :: mm10_working_data
        type(crystal_props) :: props
        type(crystal_state) :: np1, np0
        double precision, dimension(max_uhard) :: vec1,vec2
        double precision, dimension(max_uhard,max_uhard) :: arr1,arr2
        complex(kind=real64), dimension(max_uhard) :: ivec1, ivec2
        integer gaspt, solvfnc
        double precision, dimension(max_uhard) :: x2
        logical :: locdebug
      end type
c
      end module mm10_work_def
