module hydrology_methods

use real_precision
use system_state
use site_parameters
use tuning_parameters
use pft_parameters

implicit none

contains

!**********************************************************************!
!                                                                      !
!                       hydrology :: hydrology_methods                 !
!                       ------------------------------                 !
!                                                                      !
! subroutine hydrology(adp,sfc,sw,sswc,awl,kd,kx,eemm,etmm,pet2,prc,   !
! s1in,t,rlai,evap,tran,roff,interc,evbs,f2,f3,ft)                     !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Inclusion of century water model (which is a 'bucket model)
!! using doly evaporation and interception calculations.
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine hydrology(adp,sfc,sw,sswc,awl,kd,kx,eemm,etmm,pet2,prc, &
 s1in,t,rlai,evap,tran,roff,interc,evbs,f2,f3,ft)
!**********************************************************************!
real(dp) :: adp(4),sfc(4),sw(4),sswc(4),awl(4),kd,kx,tf(13)
real(dp) :: ladp(4),lsfc(4),lsw(4),lsswc(4),s1in
real(dp) :: eemm,etmm,pet2,dp2,t,rlai,evap,tran,roff,f2,f3
real(dp) :: bst,kf,interc,ms,fs,rwc(4),w(4),rem,f1,ds,st,ws,sl,sf,sd
real(dp) :: fav,ev,ss,sums,evbs,pet3,ans1,tfscale,prc
integer :: lai,iter,ft,i
!----------------------------------------------------------------------!

iter = 30

dp2 = prc

!----------------------------------------------------------------------!
! SET parameter for THROUGHFALL.                                       !
!----------------------------------------------------------------------!
tf(1) = 1.0
tf(2) = 0.95
tf(3) = 0.935
tf(4) = 0.92
tf(5) = 0.905
tf(6) = 0.89
tf(7) = 0.875
tf(8) = 0.86
tf(9) = 0.845
tf(10)= 0.83
tf(11)= 0.815
tf(12)= 0.80
tf(13)= 0.785

tfscale = 2.75
do i=1,13
  tf(i) = 1.0 - tfscale*(1.0 - tf(i))
enddo

ds = 0.0
st = 0.0

kf = 1.0 - kd
rem = rlai - int(rlai)
lai = int(rlai) + 1

pet3 = pet2
pet2 = 3.0*pet2

!----------------------------------------------------------------------!
! Adjustment for the 'CITY' functional type.                           !
!----------------------------------------------------------------------!
if (pft(ft)%itag==2) then
  do i=1,4
    ladp(i) = adp(i)*tgp%p_city_dep/ssp%soil_depth/10.0
    lsfc(i) = sfc(i)*tgp%p_city_dep/ssp%soil_depth/10.0
    lsw(i) = sw(i)*tgp%p_city_dep/ssp%soil_depth/10.0
    lsswc(i) = sswc(i)*tgp%p_city_dep/ssp%soil_depth/10.0
  enddo
else
  do i=1,4
    ladp(i) = adp(i)
    lsfc(i) = sfc(i)
    lsw(i) = sw(i)
    lsswc(i) = sswc(i)
  enddo
endif

if (t<0.0) then
  s1in = 0.0
else
  s1in = dp2
endif

!----------------------------------------------------------------------!
if (t<0.0) then
!----------------------------------------------------------------------!
! Snow calculation (mm day-1).                                         !
!----------------------------------------------------------------------!
  ssv(ft)%snow = ssv(ft)%snow + dp2
  evap = 0.0
  interc = 0.0
else
!----------------------------------------------------------------------!
! Interception water loss (evap mm day-1).                             !
!----------------------------------------------------------------------!
  if (rlai>0) then
    interc = dp2*(1.0 - (tf(lai) + rem*(tf(lai+1) - tf(lai))))
    interc = interc*min(1.0,(0.5 + t/32.0))
    evap = eemm
    if (evap>pet2) evap = pet2
    if (evap>interc)  evap = interc
    pet2 = pet2 - evap
    dp2 = dp2 - evap
    interc = evap
  else
    evap = 0.0
    interc = 0.0
  endif

  if (ssv(ft)%snow>0.0) then
    ssv(ft)%l_snow = ssv(ft)%l_snow + dp2
  else
!----------------------------------------------------------------------!
! Soil water input 0-150mm (s1).                                       !
!----------------------------------------------------------------------!
    ssv(ft)%soil_h2o(1) = ssv(ft)%soil_h2o(1) + dp2
  endif
endif

if (ssv(ft)%snow>0.0) then
!----------------------------------------------------------------------!
! Snow sublimination.                                                  !
!----------------------------------------------------------------------!
  sl = 0.85*pet2/50.0
  if (sl>ssv(ft)%snow) then
    sl = ssv(ft)%snow
    pet2 = pet2 - sl
    ssv(ft)%snow = 0.0
  else
    ssv(ft)%snow = ssv(ft)%snow - sl
    pet2 = pet2 - sl
!----------------------------------------------------------------------!
! Liquid snow input (lsn - mm day-1).                                  !
!----------------------------------------------------------------------!
    ms = 0.0
    if (t>0.0) then
      ms = t*40.0/real(iter)
      if (ms>ssv(ft)%snow) then
        ms = ssv(ft)%snow
      endif
      ssv(ft)%snow = ssv(ft)%snow - ms
      ssv(ft)%l_snow = ssv(ft)%l_snow + ms
    endif
  endif
else
  sl = 0.0
endif

!----------------------------------------------------------------------!
! Soil water input for first layer.                                    !
!----------------------------------------------------------------------!
fs = 0.0
if ((ssv(ft)%l_snow>0.05*(ssv(ft)%l_snow + ssv(ft)%snow)).and.(t>0.0)) &
 then
  fs = ssv(ft)%l_snow - 0.05*(ssv(ft)%l_snow + ssv(ft)%snow)
  ssv(ft)%l_snow = 0.05*(ssv(ft)%l_snow + ssv(ft)%snow)
  ssv(ft)%soil_h2o(1) = ssv(ft)%soil_h2o(1) + fs
  s1in = s1in + fs
endif

!----------------------------------------------------------------------!
! Soil water 150-300mm (ssv(ft)%soil_h2o(2)).                                           !
!----------------------------------------------------------------------!
f1 = 0.0
if (ssv(ft)%soil_h2o(1)>lsfc(1)) then
  f1 = ssv(ft)%soil_h2o(1) - lsfc(1)
  ssv(ft)%soil_h2o(1) = lsfc(1)
endif
ssv(ft)%soil_h2o(2) = ssv(ft)%soil_h2o(2) + f1

!----------------------------------------------------------------------!
! Soil water 150-300mm (ssv(ft)%soil_h2o(2)).                                           !
!----------------------------------------------------------------------!
f2 = 0.0
if (ssv(ft)%soil_h2o(2)>lsfc(2)) then
  f2 = ssv(ft)%soil_h2o(2) - lsfc(2)
  ssv(ft)%soil_h2o(2) = lsfc(2)
endif
ssv(ft)%soil_h2o(3) = ssv(ft)%soil_h2o(3) + f2

!----------------------------------------------------------------------!
! Soil water 300-450mm (ssv(ft)%soil_h2o(3)).                                           !
!----------------------------------------------------------------------!
f3 = 0.0
if (ssv(ft)%soil_h2o(3)>lsfc(3)) then
  f3 = ssv(ft)%soil_h2o(3) - lsfc(3)
  ssv(ft)%soil_h2o(3) = lsfc(3)
endif
ssv(ft)%soil_h2o(4) = ssv(ft)%soil_h2o(4) + f3

!----------------------------------------------------------------------!
! Fill up to saturated water content from bottom up.                   !
!----------------------------------------------------------------------!
roff = 0.0
ans1 = lsfc(4) + tgp%p_roff2*(lsswc(4) - lsfc(4))
if (ssv(ft)%soil_h2o(4)>ans1) then
  ssv(ft)%soil_h2o(3) = ssv(ft)%soil_h2o(3) + ssv(ft)%soil_h2o(4) - ans1
  ssv(ft)%soil_h2o(4) = ans1
  ans1 = lsfc(3) + tgp%p_roff2*(lsswc(3) - lsfc(3))
  if (ssv(ft)%soil_h2o(3)>ans1) then
    ssv(ft)%soil_h2o(2) = ssv(ft)%soil_h2o(2) + ssv(ft)%soil_h2o(3) - &
 ans1
    ssv(ft)%soil_h2o(3) = ans1
    ans1 = lsfc(2) + tgp%p_roff2*(lsswc(2) - lsfc(2))
    if (ssv(ft)%soil_h2o(2)>ans1) then
      ssv(ft)%soil_h2o(1) = ssv(ft)%soil_h2o(1) + &
 ssv(ft)%soil_h2o(2) - ans1
      ssv(ft)%soil_h2o(2) = ans1
      ans1 = lsfc(1) + tgp%p_roff2*(lsswc(1) - lsfc(1))
      if (ssv(ft)%soil_h2o(1)>ans1) then
        roff = ssv(ft)%soil_h2o(1) - ans1
        ssv(ft)%soil_h2o(1) = ans1
      endif
    endif
  endif
endif

!----------------------------------------------------------------------!
! Deep storage (ds - mm day-1) and stream (st - mm day-1).             !
!----------------------------------------------------------------------!
sf = 0.0
sd = 0.0
if (ssv(ft)%soil_h2o(4)>lsfc(4)) then
  sf = kf*(ssv(ft)%soil_h2o(4) - lsfc(4))*tgp%p_roff
  sd = kd*(ssv(ft)%soil_h2o(4) - lsfc(4))*tgp%p_roff
  roff = roff + sf + sd
!      ELSE
!        roff = 0.0
endif
ss = ds*kx
ds = ds + sd - ss
st = st + sf + ss
ssv(ft)%soil_h2o(4) = ssv(ft)%soil_h2o(4) - sf - sd

!----------------------------------------------------------------------!
! Calculation of transpiration (tran - mm day-1).                      !
!----------------------------------------------------------------------!
fav = 0.0

if (ssv(ft)%soil_h2o(1)>lsw(1)) then
  fav = fav + ssv(ft)%soil_h2o(1) - lsw(1)
  rwc(1) = (ssv(ft)%soil_h2o(1) - lsw(1))/(lsfc(1)-lsw(1))
  rwc(1) = 0.0
else
  rwc(1) = 0.0
endif

if (ssv(ft)%soil_h2o(2)>lsw(2)) then
  fav = fav + ssv(ft)%soil_h2o(2) - lsw(2)
  rwc(2) = (ssv(ft)%soil_h2o(2) - lsw(2))/(lsfc(2)-lsw(2))
else
  rwc(2) = 0.0
endif

if (ssv(ft)%soil_h2o(3)>lsw(3)) then
  fav = fav + ssv(ft)%soil_h2o(3) - lsw(3)
  rwc(3) = (ssv(ft)%soil_h2o(3) - lsw(3))/(lsfc(3)-lsw(3))
else
  rwc(3) = 0.0
endif

if (ssv(ft)%soil_h2o(4)>lsw(4)) then
  fav = fav + ssv(ft)%soil_h2o(4) - lsw(4)
  rwc(4) = (ssv(ft)%soil_h2o(4) - lsw(4))/(lsfc(4) - lsw(4))
else
  rwc(4) = 0.0
endif

w(1) = rwc(1)*awl(1)*ladp(1)
w(2) = rwc(2)*awl(2)*ladp(2)
w(3) = rwc(3)*awl(3)*ladp(3)
w(4) = rwc(4)*awl(4)*ladp(4)
ws = w(1) + w(2) + w(3) + w(4)

tran = etmm
if (tran>pet2)  tran = pet2
pet2 = pet2 - tran

if (ws>1e-6) then
  ssv(ft)%soil_h2o(1) = ssv(ft)%soil_h2o(1) - tran*w(1)/ws
  ssv(ft)%soil_h2o(2) = ssv(ft)%soil_h2o(2) - tran*w(2)/ws
  ssv(ft)%soil_h2o(3) = ssv(ft)%soil_h2o(3) - tran*w(3)/ws
  ssv(ft)%soil_h2o(4) = ssv(ft)%soil_h2o(4) - tran*w(4)/ws
else
  sums = ssv(ft)%soil_h2o(2) + ssv(ft)%soil_h2o(3) + ssv(ft)%soil_h2o(4)
  if (sums>1e-6) then
     ssv(ft)%soil_h2o(1) = ssv(ft)%soil_h2o(1) - &
 ssv(ft)%soil_h2o(1)*tran/sums
     ssv(ft)%soil_h2o(2) = ssv(ft)%soil_h2o(2) - &
 ssv(ft)%soil_h2o(2)*tran/sums
     ssv(ft)%soil_h2o(3) = ssv(ft)%soil_h2o(3) - &
 ssv(ft)%soil_h2o(3)*tran/sums
     ssv(ft)%soil_h2o(4) = ssv(ft)%soil_h2o(4) - &
 ssv(ft)%soil_h2o(4)*tran/sums
! changed by Ghislain 6/10/03
  else
     ssv(ft)%soil_h2o(2) = 0.0
     ssv(ft)%soil_h2o(3) = 0.0
     ssv(ft)%soil_h2o(4) = 0.0
     tran = 0.0
  endif
endif

if (ssv(ft)%soil_h2o(1)<0.0) then
  ssv(ft)%soil_h2o(2) = ssv(ft)%soil_h2o(2) + ssv(ft)%soil_h2o(1)
  ssv(ft)%soil_h2o(1) = 0.0
endif
if (ssv(ft)%soil_h2o(2)<0.0) then
  ssv(ft)%soil_h2o(3) = ssv(ft)%soil_h2o(3) + ssv(ft)%soil_h2o(2)
  ssv(ft)%soil_h2o(2) = 0.0
endif
if (ssv(ft)%soil_h2o(3)<0.0) then
  ssv(ft)%soil_h2o(4) = ssv(ft)%soil_h2o(4) + ssv(ft)%soil_h2o(3)
  ssv(ft)%soil_h2o(3) = 0.0
endif

if (ssv(ft)%soil_h2o(4)<0.0) then
  ssv(ft)%soil_h2o(3) = ssv(ft)%soil_h2o(3) + ssv(ft)%soil_h2o(4)
  ssv(ft)%soil_h2o(4) = 0.0
  if (ssv(ft)%soil_h2o(3)<0.0) then
    ssv(ft)%soil_h2o(2) = ssv(ft)%soil_h2o(2) + ssv(ft)%soil_h2o(3)
    ssv(ft)%soil_h2o(3) = 0.0
    if (ssv(ft)%soil_h2o(2)<0.0) then
      ssv(ft)%soil_h2o(1) = ssv(ft)%soil_h2o(1) + ssv(ft)%soil_h2o(2)
      ssv(ft)%soil_h2o(2) = 0.0
      if (ssv(ft)%soil_h2o(1)<0.0) then
        tran = tran + ssv(ft)%soil_h2o(1)
        if (tran<0.0)  tran = 0.0
        ssv(ft)%soil_h2o(1) = 0.0
        write(*,*) 'No water at all.'
      endif
    endif
  endif        
endif

!----------------------------------------------------------------------!
! Calculation of bare soil evaporation 'evbs' (mm day-1).              !
!----------------------------------------------------------------------!
!        ev = (rwc(1) - 0.25)/0.75
!        ev = (rwc(1) - 0.01)/0.75
ev = (ssv(ft)%soil_h2o(2) - lsw(2))/(lsfc(2) - lsw(2))
if (ev<0.0)  ev = 0.0
bst = 0.0
if (t>0) bst = (t/16.0)
evbs = ev*tgp%p_bs*0.33*pet3*1.3*bst

if (evbs>pet2)  evbs = pet2
pet2 = pet2 - evbs
!      IF (evbs.GT.ssv(ft)%soil_h2o(1)) THEN
!        IF (evbs.GT.(ssv(ft)%soil_h2o(1)+ssv(ft)%soil_h2o(2))) &
! evbs = ssv(ft)%soil_h2o(1) + ssv(ft)%soil_h2o(2)
!        ssv(ft)%soil_h2o(2) = ssv(ft)%soil_h2o(2) - &
! evbs + ssv(ft)%soil_h2o(1)
!        ssv(ft)%soil_h2o(1) = 0.0
!      ELSE
!        ssv(ft)%soil_h2o(1) = ssv(ft)%soil_h2o(1) - evbs
!      ENDIF
if (evbs>ssv(ft)%soil_h2o(1)) then
  evbs = ssv(ft)%soil_h2o(1)
  ssv(ft)%soil_h2o(1) = 0.0
else
  ssv(ft)%soil_h2o(1) = ssv(ft)%soil_h2o(1) - evbs
endif

evap = evap + sl + evbs

do i=1,29
  ssv(ssp%cohort)%sm_trig(31-i) = ssv(ssp%cohort)%sm_trig(30-i)
enddo
ssv(ssp%cohort)%sm_trig(1) = (s1in - 0.0*evbs)

!----------------------------------------------------------------------!
end subroutine hydrology
!----------------------------------------------------------------------!

!**********************************************************************!
end module hydrology_methods
!**********************************************************************!
