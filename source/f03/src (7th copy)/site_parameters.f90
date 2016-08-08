module site_parameters

use real_precision, only: dp
use dims,           only: max_cohorts, max_pftps

private :: dp, max_cohorts, max_pftps

public :: ssp

type SiteParameters
  integer                  :: day
  integer                  :: mnth
  integer                  :: year
  integer                  :: iyear
  real(dp)                 :: lat
  real(dp)                 :: lon
  real(dp)                 :: latres
  real(dp)                 :: lonres
  integer                  :: nft
  integer                  :: lai
  integer                  :: cohort
  integer                  :: jday
  real(dp)                 :: soil_depth
  real(dp)                 :: sand
  real(dp)                 :: silt
  real(dp)                 :: clay
  real(dp)                 :: bulk
  real(dp)                 :: orgc
  real(dp)                 :: wilt
  real(dp)                 :: field
  real(dp)                 :: sat
  real(dp), dimension(200) :: tmem
  real(dp), dimension(4)   :: new_soil_h2o
  real(dp)                 :: new_snow
  real(dp)                 :: new_l_snow
  real(dp), dimension(8)   :: new_c
  real(dp), dimension(8)   :: new_n
  real(dp), dimension(3)   :: new_minn
  real(dp)                 :: new_slc
  real(dp)                 :: new_rlc
  real(dp)                 :: new_sln
  real(dp)                 :: new_rln
  real(dp)                 :: new_cov
  integer                  :: cohorts
  integer, dimension(max_cohorts,max_cohorts) :: co2ftmap
  real(dp), dimension(max_pftps)              :: ftcov
end type

type (SiteParameters) ssp

end module site_parameters
