!> Precision 'dp' given by 'selected_real_kind(P,R).
!! where P is the precision ie. number of significant digits
!! and R is the range +- 10**R
!! standard single precision (32 bit), P=7, R = 38
!! double precision (64 bit),          P=15, R=308

module real_precision

  integer, parameter :: dp = selected_real_kind(15, 308)
!  integer, parameter :: dp = selected_real_kind(7, 38)

end module real_precision
