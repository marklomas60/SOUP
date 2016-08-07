module misc_values

use real_precision

type MiscValues
  real(dp) :: mv_leafmol
  real(dp) :: mv_soilw
  real(dp) :: mv_soil2g
  real(dp) :: mv_respref
end type

type(MiscValues) msv

end module misc_values

