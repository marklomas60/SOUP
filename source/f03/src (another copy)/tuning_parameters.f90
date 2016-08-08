module tuning_parameters

use real_precision, only: dp

private :: dp

public :: tgp

type TuningParameters
  real(dp) :: p_kscale
  real(dp) :: p_stmin
  real(dp) :: p_et
  real(dp) :: p_bs
  real(dp) :: p_roff
  real(dp) :: p_pet
  real(dp) :: p_fprob
  real(dp) :: p_rootfr
  real(dp) :: p_stemfr
  real(dp) :: p_resp
  real(dp) :: p_smbbc
  real(dp) :: p_sla
  real(dp) :: p_nu1
  real(dp) :: p_nu2
  real(dp) :: p_nu3
  real(dp) :: p_nu4
  real(dp) :: p_nleaf
  real(dp) :: p_dresp
  real(dp) :: p_vm
  real(dp) :: p_kgw
  real(dp) :: p_v1
  real(dp) :: p_v2
  real(dp) :: p_v3
  real(dp) :: p_v4
  real(dp) :: p_j1
  real(dp) :: p_j2
  real(dp) :: p_j3
  real(dp) :: p_j4
  real(dp) :: p_roff2
  real(dp) :: p_city_dep
  real(dp) :: p_opt
  real(dp) :: p_laimem
end type

type (TuningParameters) tgp

end module tuning_parameters
