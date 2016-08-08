!<@brief The system state vector containing all carbon and
!!hydrological, vegetation and soil pools.
!>
!>@details Defines the system state vector type and creates an instance
!>an array of this type called ssv which is then used throughout the
!>code. Also created is an instance of the type called ssv_temp.
!>
!>@author Mark Lomas
!>@date July 2016
!>@param sldkfj sldkfj
!>@param sldkfj sldkfj
!>@todo
module compartments
  use real_precision, only: dp
  use dims,           only: max_cohorts, max_lai_comps, max_root_comps, &
 max_stem_comps, max_suma_comps

  type ValAge
    real(dp) :: val
    integer :: age
  end type ValAge

  type StemComps
    integer :: no
    type(ValAge), dimension(max_stem_comps) :: c
    real(dp) :: tot
!    contains
!      function total()
!
!      end function total
  end type

  type RootComps
    integer :: no
    type(ValAge), dimension(max_root_comps) :: c
    real(dp) :: tot
!    contains
!      function total()
!
!      end function total
  end type

  type LaiComps
    integer :: no
    type(ValAge), dimension(max_lai_comps) :: c
    real(dp) :: tot
!    contains
!      function total()
!
!      end function total
  end type

  type SumaComps
    integer :: no
    type(ValAge), dimension(max_suma_comps) :: c
    real(dp) :: tot
!    contains
!      function total()
!
!      end function total
  end type

end module compartments



module system_state
use real_precision, only: dp
use dims,           only: max_cohorts, max_lai_comps, max_root_comps, &
 max_stem_comps, max_suma_comps
use compartments

private :: dp, max_cohorts, max_lai_comps, max_root_comps, &
 max_stem_comps, max_suma_comps

public :: ssv, ssv_temp

type SystemState
!----------------------------------------------------------------------*
  real(dp), dimension(12,2) :: assj
  real(dp), dimension(12)   :: assv
!----------------------------------------------------------------------*
  real(dp)                                :: age
  real(dp), dimension(2)               :: bio
  real(dp)                                :: cov
  real(dp)                                :: ppm
  real(dp)                                :: hgt
!----------------------------------------------------------------------*
  real(dp), dimension(4)               :: soil_h2o
  real(dp)                                :: snow
  real(dp)                                :: l_snow
!----------------------------------------------------------------------*
  real(dp), dimension(8)               :: c
  real(dp), dimension(8)               :: n
  real(dp), dimension(3)               :: minn
  real(dp), dimension(3)               :: nppstore
  real(dp)                                :: slc
  real(dp)                                :: rlc
  real(dp)                                :: sln
  real(dp)                                :: rln
!----------------------------------------------------------------------*
  real(dp), dimension(30)              :: sm_trig
  integer                             :: bb
  integer                             :: ss
  integer                             :: bbgs
  integer                             :: dsbb
  integer                             :: chill
  integer                             :: dschill
!----------------------------------------------------------------------*
  type(LaiComps)  :: lai
  type(StemComps) :: stem
  type(RootComps) :: root
  type(SumaComps) :: suma
!----------------------------------------------------------------------*
  real(dp)                                :: stemfr
!----------------------------------------------------------------------*
  real(dp)                                :: npp
  real(dp)                                :: nps
  real(dp)                                :: npl
  real(dp)                                :: npr
  real(dp)                                :: evp
!----------------------------------------------------------------------*
  real(dp)                                :: dslc
  real(dp)                                :: drlc
  real(dp)                                :: dsln
  real(dp)                                :: drln
!----------------------------------------------------------------------*
end type

type(SystemState) :: ssv(max_cohorts), ssv_temp



end module system_state






