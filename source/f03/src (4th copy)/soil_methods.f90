module soil_methods

use real_precision
use state_methods
use misc_values
use tuning_parameters

implicit none

contains

!**********************************************************************!
!                                                                      !
!                       soil_dynamics2 :: soil_methods                 !
!                       ------------------------------                 !
!                                                                      !
! SUBROUTINE soil_dynamics2(pet,prc,tmp,flow1,flow2,nfix,              !
! nci,sresp,lch,ca,site,year,yr0,yrf,speedc,ft,check_closure)          !
!                                                                      !
!**********************************************************************!
subroutine soil_dynamics2(pet,prc,tmp,flow1,flow2,nfix, &
 nci,sresp,lch,ca,site,year,yr0,yrf,speedc,ft,check_closure)
!**********************************************************************!
real(dp) :: ts,tc,pet,prc,tmp,h2o,flow1,flow2,slnr,rlnr,scl,rcl,ca, &
 nfix,sresp,lch,nci(4),ar,fms,fmr,tm,cal,cap,cas,csp,csa,kr(8),ml,al, &
 c(8),c0(8),n0(8),n(8),slc,rlc,sln,scn,t0,t1,fl(14),minn(3),xft, &
 total_carbon,old_total_carbon,sum
real(dp), dimension(max_cohorts) :: mpet,mprc,mtmp,mh2o,mflow1,mflow2
integer :: ij,yr0,yrf,year,site,ft,i,sddel(max_cohorts),co
logical :: speedc,check_closure
save mpet,mprc,mtmp,mh2o,mflow1,mflow2,sddel
!----------------------------------------------------------------------!

co = ssp%cohort

if ((ssp%iyear == 1).and.(ssp%mnth == 1).and.(ssp%day == 1)) then
  mpet(co) = 0.0
  mprc(co) = 0.0
  mtmp(co) = 0.0
  mh2o(co) = 0.0
  mflow1(co) = 0.0
  mflow2(co) = 0.0
  sddel(co) = 0
  sresp = 0.0
  lch = 0.0
endif

mpet(co) = mpet(co) + pet
mprc(co) = mprc(co) + prc
mtmp(co) = mtmp(co) + tmp
mh2o(co) = mh2o(co) + msv%mv_soilw
mflow1(co) = mflow1(co) + flow1
mflow2(co) = mflow2(co) + flow2
sddel(co) = sddel(co) + 1

if (ssp%day == 30) then
  mpet(co) = mpet(co)/real(sddel(co))*30.0
  mprc(co) = mprc(co)/real(sddel(co))*30.0
  mtmp(co) = mtmp(co)/real(sddel(co))
  mh2o(co) = mh2o(co)/real(sddel(co))
  mflow1(co) = mflow1(co)/real(sddel(co))*30.0
  mflow2(co) = mflow2(co)/real(sddel(co))*30.0

  call SOIL_DYNAMICS(mpet(co),mprc(co),mtmp(co),mh2o(co),mflow1(co),&
 mflow2(co),nfix,nci,sresp,lch,ca,site,year,yr0,yrf,speedc,ft,&
 check_closure)
 sresp = sresp/30.0
 lch = lch/30.0

  mpet(co) = 0.0
  mprc(co) = 0.0
  mtmp(co) = 0.0
  mh2o(co) = 0.0
  mflow1(co) = 0.0
  mflow2(co) = 0.0
  sddel(co) = 0
endif

!**********************************************************************!
end subroutine soil_dynamics2
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                       soil_dynamics :: soil_methods                  !
!                       ------------------------------                 !
!                                                                      !
! SUBROUTINE soil_dynamics(avmnpet,avmnppt,avmnt,h2o,flow1,flow2,nfix, !
! nci,srespm,lchm,ca,site,year,yr0,yrf,speedc,ft,check_closure)        !
!                                                                      !
!**********************************************************************!
subroutine soil_dynamics(avmnpet,avmnppt,avmnt,h2o,flow1,flow2,nfix, &
 nci,srespm,lchm,ca,site,year,yr0,yrf,speedc,ft,check_closure)
!**********************************************************************!
real(dp) :: ts,tc,avmnpet,avmnppt,avmnt,h2o,flow1,flow2,slnr,rlnr,scl,&
 rcl,ca,nfix,srespm,lchm,nci(4),ar,fms,fmr,tm,cal,cap,cas,csp,csa,&
 kr(8),ml,al,c(8),c0(8),n0(8),n(8),slc,rlc,sln,scn,t0,t1,fl(14),&
 minn(3),xft,total_carbon,old_total_carbon,sum
integer :: ij,yr0,yrf,year,site,ft,i
logical :: speedc,check_closure
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Check carbon closure.
!----------------------------------------------------------------------!
  if (check_closure) call SUM_CARBON(old_total_carbon,.false.)
!----------------------------------------------------------------------!

ts = ssp%sand
tc = ssp%clay

slc = ssv(ft)%slc*12.0/ssv(ft)%cov
rlc = ssv(ft)%rlc*12.0/ssv(ft)%cov
sln = ssv(ft)%sln*12.0/ssv(ft)%cov

ssv(ft)%slc = 0.0
ssv(ft)%rlc = 0.0
ssv(ft)%sln = 0.0

sum = 0.0
do i=1,8
  sum = sum + ssv(ft)%c(i)
  c0(i) = ssv(ft)%c(i)
  n0(i) = ssv(ft)%n(i)
enddo
sum = sum + slc/12.0 + rlc/12.0

do i=1,3
  minn(i) = ssv(ft)%minn(i)
enddo

!----------------------------------------------------------------------!
! Calculate lignin to nitrogen ratios based on the C02 concentration.  !
!----------------------------------------------------------------------!
slnr = 0.18*ca + 2.7
rlnr = 0.18*ca + 2.7

scl = 0.12   ! Fractions of surface carbon that are lignin.
rcl = 0.35   ! Fractions of root carbon that are lignin.

if (slc*sln>1e-6) then
  scn = slc/sln ! Ratio of surface carbon to nitrogen.
else
  scn = 10.0
endif

!----------------------------------------------------------------------!
! Compute monthly accumulation of values required by century.          !
!----------------------------------------------------------------------!
call SETCENPAR(ts,tc,avmnpet,avmnppt,avmnt,h2o,flow2,slnr, &
 rlnr,ar,fms,fmr,scl,rcl,tm,xft,cal,cap,cas,csp,csa,kr)

!----------------------------------------------------------------------!
! Using century model compute the new carbon and nitrogen pools.       !
!----------------------------------------------------------------------!

t0 = 0.0
t1 = 1.0/12.0
srespm = 0.0
lchm = 0.0

if (c0(1)+c0(2)+c0(3)+c0(4)+c0(5)+c0(6)+c0(7)+c0(8)>1000.0) then
  call CDYN(c0,c,t0,t1,slc,rlc,scl,rcl,kr,csa,cas,csp,cap, &
 cal,fms,fmr,xft,srespm,lchm,year,yr0,yrf,speedc)

  call FLOWS(c0,c,t0,t1,fl,slc,rlc,scl,rcl,fms,fmr,csa,csp,cas,cap)

  call NDYN(n0,n,minn,t0,t1,fl,c,flow1,flow2,ts,ml,al,scn,nfix,cal,nci,site)

  if (srespm<0.0) then
    write(*,*) 'Error: soil respiration negative'
    write(*,'(''c : '',8f8.0)') c0
    write(*,'(''slc rlc scl rcl : '',4f8.3)') slc,rlc,scl,rcl
    stop
  endif

else

  do ij=1,8
    c(ij) = c0(ij)
    n(ij) = n0(ij)
  enddo

  c(1) = c(1) + (1.0 - fms)*slc/12.0
  c(2) = c(2) + (1.0 - fmr)*rlc/12.0
  c(5) = c(5) + fms*slc/12.0
  c(6) = c(6) + fmr*rlc/12.0

  srespm = 0.0
  lchm = 0.0

endif

!CALL SOILCLOSS(c0,soilcl,scl,kr,rcl,xft,cal)

do i=1,8
  sum = sum - c(i)
  ssv(ft)%c(i) = c(i)
  ssv(ft)%n(i) = n(i)
enddo

do i=1,3
  ssv(ft)%minn(i) = minn(i)
enddo

!----------------------------------------------------------------------!
! Check carbon closure.
!----------------------------------------------------------------------!
if (check_closure) then
  call SUM_CARBON(total_carbon,.false.)
  if (abs(total_carbon-old_total_carbon+(srespm+lchm)*ssv(ft)%cov) > &
 1.0e-1) then
    write(*,*) 'Breach of carbon closure in DOLYMONTH:',&
 total_carbon+(srespm+lchm)*ssv(ft)%cov-old_total_carbon,' g/m^2.'
    write(*,*) &
 total_carbon,srespm*ssv(ft)%cov,lchm*ssv(ft)%cov,old_total_carbon
    stop
  endif
endif

!**********************************************************************!
end subroutine soil_dynamics
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                       fcap :: soil_methods                           !
!                       --------------------                           !
!                                                                      !
! FCAP computes the field capacities and the wilting points of the     !
! soil layers, from the soil sand and clay contents.                   !
!                                                                      !
! subroutine fcap(ts,tc,adp,sfc,sw)                                    !
!                                                                      !
!**********************************************************************!
subroutine fcap(ts,tc,adp,sfc,sw)
real(dp) :: ts,tc,adp(4),sfc(4),sw(4),t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,&
 t11,t12,t13,t14,t15,x,y,class(12)
!**********************************************************************!
integer sumx
!----------------------------------------------------------------------!

x=ts*100.00
y=tc*100.00

t1 = x - 44.32
t2 = y - 39.33
t3 = y + x - 59.55
if ((t3>=0.0).and.(t2>=0.0).and.(t1<=0.0)) then
  sumx = 1
else
  t4 = y - 34.83
  if ((t1>=0.0).and.(t4>=0.0)) then
  sumx = 5
  else
    if (t2>=0.0) then
      sumx = 2
    else
      t10 = y - 19.78
      t15 = y + x - 71.78
      if ((t15>=0.0).and.(t10>=0.0).and.(t1>=0.0)) then
        sumx = 6
      else
        t5 = y - 26.97
        t6 = x - 21.59
        if ((t5>=0.0).and.(t6>=0.0)) then
          sumx = 4
        else
          if (t5>=0.0) then
            sumx = 3
          else
            t7 = y + x - 50.00
            t11 = x - 52.27
            t12 = y - 7.87
            if ((t7>=0.0).and.(t12>=0.0).and.(t11<=0.0)) then
              sumx = 7
            else
              t13 = y - 0.98*x + 69.04
              if ((t7>=0.0).and.(t13>=0.0)) then
                sumx = 12
              else
                t14 = y - 2.31*x + 197.92
                if ((t13<=0.0).and.(t14>=0.0)) then
                  sumx = 10
                else
                  if (t14<=0.0) then
                    sumx = 11
                  else
                    t8 = y + x - 20.00
                    t9 = y - 10.00
                    if ((t9>=0.0).or.(t8>=0.0)) then
                      sumx = 8
                    else
                      sumx = 9
                    endif
                  endif
                endif
              endif
            endif
          endif
        endif
      endif
    endif
  endif
endif

class(1) = 0.468      ! Clay
class(2) = 0.468      ! Silty Clay
class(3) = 0.464      ! Silty Clay Loam
class(4) = 0.465      ! Clay Loam
class(5) = 0.406      ! Sandy Clay
class(6) = 0.404      ! Sandy Clay Loam
class(7) = 0.439      ! Loam
class(8) = 0.476      ! Silty loam
class(9) = 0.476      ! Silt
class(10) = 0.421     ! Loamy Sand
class(11) = 0.339     ! Sand
class(12) = 0.434     ! Sandy Loam

sfc(1) = class(sumx)*adp(1)*1.0
sfc(2) = class(sumx)*(adp(2) - adp(1))*1.0
sfc(3) = class(sumx)*(adp(3) - adp(2))*1.0
sfc(4) = class(sumx)*(adp(4) - adp(3))*1.0
sw(1) = 10.0
sw(2) = 10.0
sw(3) = 10.0
sw(4) = 10.0

!**********************************************************************!
end subroutine fcap
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                        setcenpar :: soil_methods                     !
!                        -------------------------                     !
!                                                                      !
! SETCENPAR sets the century parameters.                               !
!                                                                      !
! subroutine setcenpar(ts,tc,pet,ppt,tmp,h2o,flow2,slnr,rlnr,ar,fms,   !
! fmr,scl,rcl,tm,ft,cal,cap,cas,csp,csa,kr)                            !
!                                                                      !
!                                                                      !
!**********************************************************************!
subroutine setcenpar(ts,tc,pet,ppt,tmp,h2o,flow2,slnr,rlnr,ar,fms, &
 fmr,scl,rcl,tm,ft,cal,cap,cas,csp,csa,kr)
!**********************************************************************!
real(dp) :: ts,tc,pet,ppt,tmp,h2o,flow2,slnr,rlnr,ar,fms,fmr,scl,rcl,&
 tm,ft,cal,cap,cas,csp,csa,kr(8),smois,lcs,lcr,t
!----------------------------------------------------------------------!
! ppt and pet are monthly totals

if (pet>0.0) then
  smois = (h2o + ppt*1.0)/pet*0.2
else
  smois = 1000.0
endif
ar =-eow(smois)*eot(tmp)
if (tmp<0.0) ar = 0.0

fms = 0.85 - 0.018*slnr
fmr = 0.85 - 0.018*rlnr

if (fms<0.001)  fms = 0.001
if (fmr<0.001)  fmr = 0.001

lcs = exp(-3.0*scl)
lcr = exp(-3.0*rcl)

t = 1.0 - ts/100.0
tm = (1.0 - 0.75*t)

ft = 0.85 - 0.68*t
cal = (flow2/18.0)*(0.01 + 0.04*ts/100.0)
cap = 0.003 + 0.032*tc/100.0
cas = (1.0 - cap - cal - ft)

csp = 0.003 - 0.009*tc/100.0
csa = 1.0 - csp - 0.55

kr(1) = 3.9*ar*lcs*tgp%p_kscale
kr(2) = 4.8*ar*lcr*tgp%p_kscale
kr(3) = 7.3*ar*tm*tgp%p_kscale
kr(4) = 6.0*ar*tgp%p_kscale
kr(5) = 14.8*ar*tgp%p_kscale
kr(6) = 18.5*ar*tgp%p_kscale
kr(7) = 0.2*ar*tgp%p_kscale
kr(8) = 0.0045*ar*tgp%p_kscale

!**********************************************************************!
end subroutine setcenpar
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                       soilcloss :: soil_methods                      !
!                       -------------------------                      !
!                                                                      !
! SOILCLOSS calculates carbon losses in the soil due to resperation    !
! and leaching, given the soil carbon pool values.                     !
!                                                                      !
! subroutine soilcloss(c,tot,scl,kr,rcl,ft,cal)                        !
!                                                                      !
!**********************************************************************!
subroutine soilcloss(c,tot,scl,kr,rcl,ft,cal)
!**********************************************************************!
implicit none
real(dp) :: c(8),tot,scl,kr(8),rcl,ft,cal
!----------------------------------------------------------------------!

tot = -0.3*scl*kr(1)*c(1) - 0.6*(1.0 - scl)*kr(1)*c(1) - &
 0.55*(1 - rcl)*kr(2)*c(2) - 0.3*rcl*kr(2)*c(2) - &
 0.6*kr(5)*c(5) - 0.55*kr(6)*c(6) - 0.6*kr(4)*c(4) - &
 0.55*kr(7)*c(7) - 0.55*kr(8)*c(8) - ft*kr(3)*c(3) - cal*kr(3)*c(3)

!**********************************************************************!
end subroutine soilcloss
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                        cdyn :: soil_methods                          !
!                        --------------------                          !
!                                                                      !
! Uses the NAG routine 2EBF to solve the differential equaitons        !
! given by the carbon dynamics.                                        !
!                                                                      !
! subroutine cdyn(c0,c1,t0,t1,slc,rlc,scl,rcl,kr,csa,cas,csp,cap,      !
! cal,fms,fmr,ft,sresp,lch,year,yr0,yrf,speedc)                        !
!                                                                      !
!**********************************************************************!
subroutine cdyn(c0,c1,t0,t1,slc,rlc,scl,rcl,kr,csa,cas,csp,cap, &
 cal,fms,fmr,ft,sresp,lch,year,yr0,yrf,speedc)
!**********************************************************************!
real(dp) :: c0(8),c1(8),t0,t1,slc,rlc,scl,rcl,kr(8),csa,cas,csp,cap, &
 cal,fms,fmr,sresp,lch,t01,t11,cp(10),cjac(10,10),ft,sc(8)
real(dp) :: cslc,crlc,cscl,crcl,ckr(8),ccsa,ccas,ccsp,ccap,ccal,cfms,&
 cfmr,cft
integer :: i,year,yrf,yr0,j
logical :: speedc
common /cder/cslc,crlc,cscl,crcl,ckr,ccsa,ccas,ccsp,ccap,ccal,cfms,&
 cfmr,cft
!----------------------------------------------------------------------!

cslc = slc
crlc = rlc
cscl = scl
crcl = rcl
ccsa = csa
ccas = cas
ccsp = csp
ccap = cap
ccal = cal
cfms = fms
cfmr = fmr
cft = ft
do i=1,8
  c1(i)=c0(i)
  ckr(i) = kr(i)
enddo

t01 = t0
t11 = t1

call DER(c1,cp)
call PDER(cjac)

if (speedc) then
  sc(7) = 10.0*1.0
  sc(8) = 500.0*1.0
  cp(7) = cp(7)*sc(7)
  cp(8) = cp(8)*sc(8)
  do i=1,10
    do j=7,8
      cjac(j,i) = cjac(j,i)*(sc(j) - (sc(j) - 1.0)/ &
 real(yrf - yr0 + 1)*real(year - yr0 + 1))
    enddo
  enddo
endif

call diffeq(t01,t11,c1,cp,cjac,sresp,lch)

!**********************************************************************!
end subroutine cdyn
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                        diffeq :: soil_methods                        !
!                        ----------------------                        !
!                                                                      !
! Gives the Jacobean corresponding to the differential carbon          !
! equations.                                                           !
!                                                                      !
! subroutine diffeq(t01,t11,c1,cp,cjac,sresp,lch)                      !
!                                                                      !
!**********************************************************************!
subroutine diffeq(t01,t11,c1,cp,cjac,sresp,lch)
!**********************************************************************!
real(dp) :: t01,t11,c1(8),cp(10),cjac(10,10),secder(10),delt,sresp,lch
integer  :: i,j
!----------------------------------------------------------------------!
 
delt = t11 - t01
 
do i=1,8
  secder(i) = 0.0
  do j=1,8
    secder(i) = secder(i) + cjac(i,j)*cp(j)
  enddo
  c1(i) = c1(i) + cp(i)*delt + secder(i)*delt**2.0/2.0
enddo
 
secder(9) = 0.0
secder(10) = 0.0
do j=1,10
  secder(9) = secder(9) + cjac(9,j)*cp(j)
  secder(10) = secder(10) + cjac(10,j)*cp(j)
enddo
sresp = sresp + cp(9)*delt + secder(9)*delt**2.0/2.0
lch = lch + cp(10)*delt + secder(10)*delt**2.0/2.0

!**********************************************************************!
end subroutine diffeq
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                         der :: soil_methods                          !
!                         -------------------                          !
!                                                                      !
! Returns the carbon derivative calculated from the carbon state.      !
!                                                                      !
! subroutine der(c,f)                                                  !
!                                                                      !
!**********************************************************************!
subroutine der(c,f)
!**********************************************************************!
implicit none
real(dp) :: c(8),f(10)
real(dp) :: cslc,crlc,cscl,crcl,ckr(8),ccsa,ccas,ccsp,ccap,ccal,cfms,&
 cfmr,cft
common /cder/cslc,crlc,cscl,crcl,ckr,ccsa,ccas,ccsp,ccap,ccal,cfms,&
 cfmr,cft
!----------------------------------------------------------------------!

f(1) = cslc*(1.0 - cfms) + ckr(1)*c(1)
f(2) = crlc*(1.0 - cfmr) + ckr(2)*c(2)
f(3) =-0.45*(1.0 - crcl)*ckr(2)*c(2) - 0.45*ckr(6)*c(6) - &
 0.45*ckr(8)*c(8) - ckr(7)*c(7)*ccsa + ckr(3)*c(3)
f(4) =-0.4*(1.0 - cscl)*ckr(1)*c(1) - 0.4*ckr(5)*c(5) + &
 ckr(4)*c(4)
f(5) = cslc*cfms + ckr(5)*c(5)
f(6) = crlc*cfmr + ckr(6)*c(6)
f(7) =-0.7*cscl*ckr(1)*c(1) - 0.4*ckr(4)*c(4) - &
 0.7*crcl*ckr(2)*c(2) - ckr(3)*c(3)*ccas + ckr(7)*c(7)
f(8) =-ckr(3)*c(3)*ccap - ckr(7)*c(7)*ccsp + ckr(8)*c(8)
f(9) =-(0.6-0.3*cscl)*ckr(1)*c(1) - (0.55-0.25*crcl)*ckr(2)*c(2) - &
 cft*ckr(3)*c(3) - 0.6*ckr(4)*c(4) - 0.6*ckr(5)*c(5) - &
 0.55*ckr(6)*c(6) - 0.55*ckr(7)*c(7) - 0.55*ckr(8)*c(8) 
f(10) = -ckr(3)*c(3)*ccal

!**********************************************************************!
end subroutine der
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                          pder :: soil_methods                        !
!                          --------------------                        !
!                                                                      !
! PDER gives the Jacobean corresponding to the differential carbon     !
! equations.                                                           !
!                                                                      !
! subroutine pder(jac)                                                 !
!                                                                      !
!**********************************************************************!
subroutine pder(jac)
!**********************************************************************!
real(dp) :: jac(10,10)
integer i,j
real(dp) :: cslc,crlc,cscl,crcl,ckr(8),ccsa,ccas,ccsp,ccap,ccal,cfms,&
 cfmr,cft
common /cder/cslc,crlc,cscl,crcl,ckr,ccsa,ccas,ccsp,ccap,ccal,cfms,&
 cfmr,cft
!----------------------------------------------------------------------!

do i=1,10
  do j=1,10
    jac(i,j) = 0.0
  enddo
enddo
jac(1,1) = ckr(1)
jac(2,2) = ckr(2)
jac(3,2) =-0.45*(1.0 - crcl)*ckr(2)
jac(3,3) = ckr(3)
jac(3,6) =-0.45*ckr(6)
jac(3,7) =-ccsa*ckr(7)
jac(3,8) =-0.45*ckr(8)
jac(4,1) =-0.4*(1.0 - cscl)*ckr(1)
jac(4,4) = ckr(4)
jac(4,5) =-0.4*ckr(5)
jac(5,5) = ckr(5)
jac(6,6) = ckr(6)
jac(7,1) =-0.7*cscl*ckr(1)
jac(7,2) =-0.7*crcl*ckr(2)
jac(7,3) =-ccas*ckr(3)
jac(7,4) =-0.4*ckr(4)
jac(7,7) = ckr(7)
jac(8,3) =-ccap*ckr(3)
jac(8,7) =-ccsp*ckr(7)
jac(8,8) = ckr(8)
jac(9,1) =-(0.6 - 0.3*cscl)*ckr(1)
jac(9,2) =-(0.55 - 0.25*crcl)*ckr(2)
jac(9,3) =-cft*ckr(3) 
jac(9,4) =-0.6*ckr(4) 
jac(9,5) =-0.6*ckr(5) 
jac(9,6) =-0.55*ckr(6) 
jac(9,7) =-0.55*ckr(7) 
jac(9,8) =-0.55*ckr(8) 
jac(9,9) = 0.0 
jac(10,3) = -ckr(3)*ccal

!**********************************************************************!
end subroutine pder
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                        clsasi :: soil_methods                        !
!                        ----------------------                        !
!                                                                      !
! This subroutine returns the persentages of soil clay sand and silt.  !
! These are estimated from the soil classification on the islscp data  !
! file 'soltex'. This results in the folowing classifications:         !
!                                                                      !
!        ISLSCP Index      Clay  Sand  Silt       Comments             !
!              1             7    80    13        loamy sand           !
!              2            12    62    26        sandy loam           !
!              3            18    42    40        loam                 !
!              4            27    63    10        sandy clay loam      !
!              5            30    35    35        clay loam            !
!              6             0     0     0        ice                  !
!              7            18    42    40        loam                 !
!              0             0     0     0        ocean                !
!                                                                      !
! subroutine clsasi(tex,tc,ts,tsi)                                     !
!                                                                      !
!**********************************************************************!
subroutine clsasi(tex,tc,ts,tsi)
!**********************************************************************!
real(dp) :: tc,ts,tsi
integer :: tex
!----------------------------------------------------------------------!

if ((tex==0).or.(tex==6)) then
  tc = 0.0
  ts = 0.0
  tsi = 0.0
else
  if (tex==1) then
    tc = 7.0
    ts = 80.0
  elseif (tex==2) then
    tc = 12.0
    ts = 62.0
  elseif (tex==3) then
    tc = 18.0
    ts = 42.0
  elseif (tex==4) then
    tc = 27.0
    ts = 63.0
  elseif (tex==5) then
    tc = 30.0
    ts = 35.0
  elseif (tex==7) then
    tc = 18.0
    ts = 42.0
  endif
  tsi = 100.0 - tc -ts
endif

!**********************************************************************!
end subroutine clsasi
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                         ndyn :: soil_methods                         !
!                         --------------------                         !
!                                                                      !
! Calculate nitrogen dynamics from the carbon dynamics and the carbon  !
! flows calculated using FLOWS.                                        !
!                                                                      !
! subroutine ndyn(n0,n,minn,t0,t1,fl,c,flow1,flow2,ts,ml,al,scn,nfix,  !
! cal,nci,site)                                                        !
!                                                                      !
!----------------------------------------------------------------------!
subroutine ndyn(n0,n,minn,t0,t1,fl,c,flow1,flow2,ts,ml,al,scn,nfix,cal,&
 nci,site)
!**********************************************************************!
real(dp) :: n0(8),n(8),minn(3),t0,t1,fl(14),c(8),flow1,flow2,ts,ml,al,&
 scn,nfix,cal,nci(4),ncsm,nca,ncs,ncp,soiln,plantn,delt,x,totn, &
 ratminn,f1,f2,nin,az,bz,cz,dz,ez,fz,gz,hz,iz
integer i,site
!----------------------------------------------------------------------!

delt = t1 - t0
soiln = minn(3)
plantn = 100.0/(scn + 1.0)
!...check whether 'scn' is correct, 'scn' = surface carbon/surface
! nitrogen, should this be litter carbon/litter nitrogen.

nca = 1.0/ncrata(soiln)

do i=1,8
  n(i)=n0(i)
enddo

!----------------------------------------------------------------------!
! Find mineral nitrogen leaching due to the water flows between the    !
! soil layers.                                                         !
!----------------------------------------------------------------------!
!      ml = minn(2)*flow2*delt*(0.2 + 0.7*ts/100.0)/18.0
!      minn(2) = minn(2) - ml
!      IF (minn(1).GT.0.0) THEN
!        temp = minn(1)*flow1*delt*(0.2 + 0.7*ts/100.0)/18.0
!        minn(1) = minn(1) - temp
!        minn(2) = minn(2) + temp
!      ENDIF

!      minn(1) = minn(1) + nfix*delt + (1.0 - nup)*(sln + rln)*delt

if (flow1<1.0e-10) flow1 = 0.0
if (flow2<1.0e-10) flow2 = 0.0

f1 = flow1*(0.2 + 0.7*ts/100.0)/18.0
f2 = flow2*(0.2 + 0.7*ts/100.0)/18.0
nin = nfix
if (abs(f1-f2)<1.0e-6) then
  f2 = f2 - 0.001
endif

if (flow1>0.0) then
  bz = nin/f1
  az = minn(1) - bz
  minn(1) = az*exp(-f1*delt) + bz
  if (flow2>0.0) then
    dz = f1*az/(f2 - f1)
    ez = f1*bz/f2
    cz = minn(2) - dz - ez
    minn(2) = cz*exp(-f2*delt) + dz*exp(-f1*delt) + ez
    fz = -cz
    gz = -dz*f2/f1
    hz = f2*ez
    iz = -fz - gz
    ml = fz*exp(-f2*delt) + gz*exp(-f1*delt) + hz*delt + iz
  else
    cz = -az
    dz = f1*bz
    ez = minn(2) - cz
    minn(2) = cz*exp(-f1*delt) + dz*delt + ez
    fz = -cz*f2/f1
    gz = f2*dz/2.0
    hz = f2*ez
    iz =-fz
    ml = fz*exp(-f1*delt) + gz*delt**2.0 + hz*delt + iz
  endif
else
  if (flow2>0.0) then
    minn(1) = nfix*delt + minn(1)
    cz = minn(2)
    minn(2) = cz*exp(-f2*delt)
    ml = -cz*exp(-f2*delt) + cz
  else
    minn(1) = nfix*delt + minn(1)
    minn(2) = minn(2)
  endif
endif

if (.not.(minn(2)<10000.0))  &
 write(11,'(''minn(2) > 10000, '',i5,7f12.4)') &
 site,minn(2),cz,f2,dz,f1,ez,delt

if (minn(3)<2.0) then
  al = fl(12)*cal/(15.0 - 6.0*minn(3))
else
  al = fl(12)*cal/3.0
endif

!----------------------------------------------------------------------!
! Calculate total nitrogen and then find how much of that is mineral.  !
!----------------------------------------------------------------------!
totn = minn(1)+minn(2)+n(1)+n(2)+n(3)+n(4)+n(5)+n(6)+n(7)+n(8)-al

if (.not.(minn(1)<10000.0))  &
 write(11,'(''minn(1) > 10000, '',i5,4f12.4)') &
 site,flow1,flow2,minn(1),minn(2)

call NBAL(c,totn,x,scn,nci)
if (.not.(minn(1)<10000.0)) &
 write(11,'(''minn(1) > 10000, '',i5,4f12.4)') &
 site,flow1,flow2,minn(1),minn(2)

!----------------------------------------------------------------------!
! Adjust all nitrogen pools accordingly.                               !
!----------------------------------------------------------------------!
if (minn(2)>0.0) then
  ratminn = minn(1)/minn(2)
  minn(1) = x/(1.0 + 1.0/ratminn)
  minn(2) = x/(1.0 + ratminn)
  minn(3) = x
else
  minn(1) = x
  minn(2) = 0.0
  minn(3) = x
endif

if (.not.(minn(1)<10000.0)) &
 write(11,'(''minn(1) > 10000, '',i5,6f12.4)') &
 site,flow1,flow2,minn(1),minn(2),x,ratminn

ncsm = 1.0/ncratsm(plantn)
nca = 1.0/ncrata(x)
ncs = 1.0/ncrats(x)
ncp = 1.0/ncratp(x)

n(1) = c(1)*nci(1)
n(5) = c(5)*nci(2)
n(2) = c(2)*nci(3)
n(6) = c(6)*nci(4)
n(4) = ncsm*c(4)
n(8) = ncp*c(8)
n(7) = ncs*c(7)
n(3) = nca*c(3)

!**********************************************************************!
end subroutine ndyn
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                        nbal :: soil_methods                          !
!                        --------------------                          !
!                                                                      !
! Computes the mineral nitrogen pool, given the value of the           !
! total nitrogen content of the soil. To do this a quartic equation    !
! needs to be solved, this is done using the Newton Raphson method.    !
!                                                                      !
! subroutine nbal(c,tn,x,scn,nci)                                      !
!                                                                      !
!**********************************************************************!
subroutine nbal(c,tn,x,scn,nci)
!**********************************************************************!
real(dp) :: c(8),tn,x,scn,nci(4),const,ncsm,plantn,fx0,fx05,fx1,fx15,&
 fx2,b,at
integer :: i
!----------------------------------------------------------------------!

plantn = 100.0/(scn + 1.0)
ncsm = 1.0/ncratsm(plantn)

const = tn - c(1)*nci(1) - c(5)*nci(2) - c(2)*nci(3) - c(6)*nci(4) - &
 c(4)*ncsm

!      cx1 =-7.2
!      cx2 = 102.0 + 7.2*const
!      cx3 = 4.8*c(8) + 1.8*c(7) + 1.2*c(3) - 450.0 - 102.0*const
!      cx4 =-36.0*c(8) - 16.5*c(7) - 14.0*c(3) + 600.0 + 450.0*const
!      cx5 = 60.0*c(8) + 30.0*c(7)  + 40.0*c(3) - 600.0*const

x=0.0
fx0 = 60.0*c(8) - 36.0*c(8)*x + 4.8*c(8)*x**2 + 30.0*c(7) - &
 16.5*c(7)*x +1.8*c(7)*x**2 + 40.0*c(3) - 14.0*c(3)*x + &
 1.2*c(3)*x**2 + 600.0*x - 450.0*x**2 + 102.0*x**3 - 7.2*x**4 - &
 600.0*const + 450.0*const*x - 102.0*const*x**2 + 7.2*const*x**3
x=0.5
 fx05 = 60.0*c(8) - 36.0*c(8)*x + 4.8*c(8)*x**2 + 30.0*c(7) - &
 16.5*c(7)*x +1.8*c(7)*x**2 + 40.0*c(3) - 14.0*c(3)*x + &
 1.2*c(3)*x**2 + 600.0*x - 450.0*x**2 + 102.0*x**3 - 7.2*x**4 - &
 600.0*const + 450.0*const*x - 102.0*const*x**2 + 7.2*const*x**3
x=1.0
 fx1 = 60.0*c(8) - 36.0*c(8)*x + 4.8*c(8)*x**2 + 30.0*c(7) - &
 16.5*c(7)*x +1.8*c(7)*x**2 + 40.0*c(3) - 14.0*c(3)*x + &
 1.2*c(3)*x**2 + 600.0*x - 450.0*x**2 + 102.0*x**3 - 7.2*x**4 - &
 600.0*const + 450.0*const*x - 102.0*const*x**2 + 7.2*const*x**3
x=1.5
 fx15 = 60.0*c(8) - 36.0*c(8)*x + 4.8*c(8)*x**2 + 30.0*c(7) - &
 16.5*c(7)*x +1.8*c(7)*x**2 + 40.0*c(3) - 14.0*c(3)*x + &
 1.2*c(3)*x**2 + 600.0*x - 450.0*x**2 + 102.0*x**3 - 7.2*x**4 - &
 600.0*const + 450.0*const*x - 102.0*const*x**2 + 7.2*const*x**3
x=2.0
 fx2 = 60.0*c(8) - 36.0*c(8)*x + 4.8*c(8)*x**2 + 30.0*c(7) - &
 16.5*c(7)*x +1.8*c(7)*x**2 + 40.0*c(3) - 14.0*c(3)*x + &
 1.2*c(3)*x**2 + 600.0*x - 450.0*x**2 + 102.0*x**3 - 7.2*x**4 - &
 600.0*const + 450.0*const*x - 102.0*const*x**2 + 7.2*const*x**3

if ((fx0<0.0).and.(fx2>0.0)) then
  if (fx15<0.0) then
    x = 0.5*fx15/(fx15 - fx2) + 1.5
  elseif (fx1<0.0) then
    x = 0.5*fx1/(fx1 - fx15) + 1.0
  elseif (fx05<0.0) then
    x = 0.5*fx05/(fx05 - fx1) + 0.5
  else
    x = 0.5*fx0/(fx0 - fx05)
  endif
  do i=1,4
    at = 60.0*c(8) - 36.0*c(8)*x + 4.8*c(8)*x**2 + &
 30.0*c(7) - 16.5*c(7)*x +1.8*c(7)*x**2 + 40.0*c(3) - &
 14.0*c(3)*x + 1.2*c(3)*x**2 + 600.0*x - 450.0*x**2 + &
 102.0*x**3 - 7.2*x**4 - 600.0*const + 450.0*const*x - &
 102.0*const*x**2 + 7.2*const*x**3
  b =-36.0*c(8) + 9.6*c(8)*x - 16.5*c(7) + 3.6*c(7)*x - &
 14.0*c(3) + 2.4*c(3)*x + 600.0 - 900.0*x + 306.0*x**2 - &
 28.8*x**3 + 450.0*const - 204.0*const*x + 21.6*const*x**2
  x = x - at/b
  enddo
else if ((fx0>0.0).and.(fx2>0.0)) then
  x = 0.0
else if ((fx0<0.0).and.(fx2<0.0)) then
  x = const - c(3)/3.0 - c(7)/12.0 - c(8)/7.0
endif

!**********************************************************************!
end subroutine nbal
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                        flows :: soil_methods                         !
!                        ---------------------                         !
!                                                                      !
! Computes the fourteen carbon flows associated with the carbon        !
! dynamics.                                                            !
!                                                                      !
! subroutine flows(c0,c,t0,t1,fl,slc,rlc,scl,rcl,fms,fmr,csa,csp,cas,  !
! cap)                                                                 !
!                                                                      !
!**********************************************************************!
subroutine flows(c0,c,t0,t1,fl,slc,rlc,scl,rcl,fms,fmr,csa,csp,cas,cap)
!**********************************************************************!
real(dp) :: c0(8),c(8),t0,t1,fl(14),slc,rlc,scl,rcl,fms,fmr,csa,csp,&
 cas,cap,delt
!----------------------------------------------------------------------!

delt = t1 - t0

fl(1) = slc*(1.0 - fms)*delt
fl(2) = slc*fms*delt
fl(7) = rlc*(1.0 - fmr)*delt
fl(8) = rlc*fmr*delt
fl(3) = (c0(1) - c(1) + fl(1))*scl
fl(4) = (c0(1) - c(1) + fl(1))*(1.0 - scl)
fl(5) = c0(5) - c(5) + fl(2)
fl(9) = (c0(2) - c(2) + fl(7))*rcl
fl(10) = (c0(2) - c(2) + fl(7))*(1.0 - rcl)
fl(11) = c0(6) - c(6) + fl(8)
fl(6) = c0(4) - c(4) + 0.4*fl(4) + 0.4*fl(5)

!----------------------------------------------------------------------!
! Solve three simultaneous equations to find fl(12),fl(13) and fl(14). !
!----------------------------------------------------------------------!
fl(12) = (c0(3) - c(3) + 0.45*fl(10) + 0.45*fl(11) + &
 csa*(c0(7) - c(7) + 0.7*fl(3) + 0.4*fl(6) + 0.7*fl(9)) + &
 0.45*(c0(8) - c(8) + csp*(c0(7) - c(7) + 0.7*fl(3) + &
 0.4*fl(6) + 0.7*fl(9))))/(1.0 - csa*cas - 0.45*csp*cas - 0.45*cap)
fl(13) = c0(7) - c(7) + 0.7*fl(3) + 0.4*fl(6) + 0.7*fl(9) + cas*fl(12)
fl(14) = c0(8) - c(8) + csp*fl(13) + cap*fl(12)

!     WRITE(*,*) 'c1',c(1)-c0(1),fl(1)-fl(3)-fl(4)
!     WRITE(*,*) 'c2',c(2)-c0(2),fl(7)-fl(9)-fl(10)
!     WRITE(*,*) 'c3',c(3)-c0(3),fl(13)*csa+fl(10)*.45+fl(11)*.45+
!   & fl(14)*.45-fl(12)
!     WRITE(*,*) 'c4',c(4)-c0(4),.4*fl(4)+.4*fl(5)-fl(6)
!     WRITE(*,*) 'c5',c(5)-c0(5),fl(2)-fl(5)
!     WRITE(*,*) 'c6',c(6)-c0(6),fl(8)-fl(11)
!     WRITE(*,*) 'c7',c(7)-c0(7),.7*fl(3)+.4*fl(6)+.7*fl(9)+
!   & cas*fl(12)-fl(13)
!     WRITE(*,*) 'c8',c(8)-c0(8),csp*fl(13)+cap*fl(12)-fl(14)

!**********************************************************************!
end subroutine flows
!**********************************************************************!





!**********************************************************************!
!                          FUNCTIONS                                   !
!**********************************************************************!

!**********************************************************************!
!                                                                      !
!                        ncrata :: soil_methods                        !
!                        ----------------------                        !
!                                                                      !
! Calculates N:C for the active organic matter pool.                   !
!                                                                      !
! function ncrata(x)                                                   !
!                                                                      !
!**********************************************************************!
real(dp) function ncrata(x)
!**********************************************************************!
real(dp) :: x
!----------------------------------------------------------------------!

if (x>2.0) then
  ncrata = 3.0
else
  ncrata = 15.0 - 6.0*x
endif
if (x<0.0) then
  ncrata = 15.0
endif

!**********************************************************************!
end function ncrata
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                         ncrats :: soil_methods                       !
!                         ----------------------                       !
!                                                                      !
! Calculates N:C for the slow organic matter pool.            !
!                                                                      !
! function ncrats(x)                                                   !
!                                                                      !
!**********************************************************************!
real(dp) function ncrats(x)
!----------------------------------------------------------------------!
real(dp) :: x
!----------------------------------------------------------------------!

if (x>2.0) then
  ncrats = 12.0
else
  ncrats = 20.0 - 4.0*x
endif
if (x<0.0) then
  ncrats = 20.0
endif

!**********************************************************************!
end function ncrats
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                       ncratp :: soil_methods                         !
!                       ----------------------                         !
!                                                                      !
! Calculates N:C for the passive organic matter pool.                  !
!                                                                      !
! real(dp) function ncratp(x)                                          !
!                                                                      !
!**********************************************************************!
real(dp) function ncratp(x)
!**********************************************************************!
real(dp) :: x
!----------------------------------------------------------------------!

if (x>2.0) then
  ncratp = 7.0
else
  ncratp = 10.0 - 1.5*x
endif
if (x<0.0) then
  ncratp = 10.0
endif

!**********************************************************************!
end function ncratp
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                       ncratsm :: soil_methods                        !
!                       -----------------------                        !
!                                                                      !
! Calculates N:C for the surface microbe pool.                         !
!                                                                      !
! real(dp) function ncratsm(x)                                         !
!                                                                      !
!**********************************************************************!
real(dp) function ncratsm(x)
!**********************************************************************!
real(dp) :: x
!----------------------------------------------------------------------!

if (x>2.0) then
  ncratsm = 10.0
else
  ncratsm = 20.0 - 5.0*x
endif
if (x<0.0) then
  ncratsm = 20.0
endif

!**********************************************************************!
end function ncratsm
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                       eot :: soil_methods                            !
!                       -------------------                            !
!                                                                      !
! Calculates the effect of temperature used to calculate the           !
! decay rate multiplier 'a'.                                           !
!                                                                      !
! real(dp) function eot(x)                                             !
!                                                                      !
!**********************************************************************!
real(dp) function eot(x)
!**********************************************************************!
real(dp) :: x
!----------------------------------------------------------------------!

if (x<0.0) then
  eot = 0.04
elseif (x<10.0) then
  eot = (0.16 - 0.03)*(x - 0.0)/(10.0 - 0.0) + 0.03
elseif (x<20.0) then
  eot = (0.50 - 0.16)*(x - 10.0)/(20.0 - 10.0) + 0.16
elseif (x<30.0) then
  eot = (0.88 - 0.50)*(x - 20.0)/(30.0 - 20.0) + 0.50
elseif (x<40.0) then
  eot = (0.96 - 0.88)*(x - 30.0)/(40.0 - 30.0) + 0.88
elseif (x<47.0) then
  eot = (0.86 - 0.96)*(x - 40.0)/(47.0 - 40.0) + 0.96
else
  eot = 0.0
endif

!**********************************************************************!
end function eot
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                       eot1 :: soil_methods                           !
!                       --------------------                           !
!                                                                      !
! Calculates the effect of temperature used to calculate the           !
! decay rate multiplier 'a'.                                           !
!                                                                      !
! real(dp) function eot1(x)                                            !
!                                                                      !
!**********************************************************************!
real(dp) function eot1(x)
!**********************************************************************!
real(dp) :: x,a,b
!----------------------------------------------------------------------!

a = (45.0 - x)/10.0
b = exp(0.076*(1.0 - a**2.63))
eot1 = (a**0.2)*b

!**********************************************************************!
end function eot1
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                       eow :: soil_methods                            !
!                       -------------------                            !
!                                                                      !
! Calculates the effect of water used to calculate the decay           !
! rate multiplier 'a'.                                                 !
!                                                                      !
! real(dp) function eow(x)                                             !
!                                                                      !
!**********************************************************************!
real(dp) function eow(x)
!**********************************************************************!
real(dp) :: x
!----------------------------------------------------------------------!

if (x<0.2) then
  eow = (0.21 - 0.09)*(x - 0.0)/(0.2 - 0.0) + 0.09
elseif (x<0.4) then
  eow = (0.48 - 0.21)*(x - 0.2)/(0.4 - 0.2) + 0.21
elseif (x<0.6) then
  eow = (0.82 - 0.48)*(x - 0.4)/(0.6 - 0.4) + 0.48
elseif (x<0.8) then
  eow = (0.97 - 0.82)*(x - 0.6)/(0.8 - 0.6) + 0.82
elseif (x<1.0) then
  eow = (0.97 - 0.97)*(x -0.8)/(1.0 - 0.8) + 0.97
elseif (x<2.0) then
  eow = (0.63 - 0.97)*(x - 1.0)/(2.0 - 1.0) + 0.97
else
  eow = (0.29 - 0.63)*(x - 2.0)/(3.0 - 2.0) + 0.63
endif

if ((eow<0.0).and.(x>1.0)) then
  eow = 0.1
endif

!**********************************************************************!
end function eow
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                       eow1 :: soil_methods                           !
!                       --------------------                           !
!                                                                      !
! Calculates the effect of water used to calculate the decay           !
! rate multiplier 'a'.                                                 !
!                                                                      !
! real(dp) function eow1(x)                                            !
!                                                                      !
!**********************************************************************!
real(dp) function eow1(x)
!**********************************************************************!
real(dp) :: x
!----------------------------------------------------------------------!

if (x<1.5) then
  eow1 = 1.0/(1.0 + 30.0*exp(-8.5*x))
else
  eow1 = 1.0 - 0.7*(x - 1.5)/1.5
endif

if (eow1<0.0001)  eow1 = 0.0001

!**********************************************************************!
end function eow1
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                        wsparam :: soil_methods                       !
!                        -----------------------                       !
!                                                                      !
! Sets some system parameters for both the water and soils             !
! models.                                                              !
!                                                                      !
! subroutine wsparam(l_b_and_c,nupc,awl,kd,kx,nci,nfix,adp,topsl,sfc,  !
! sw,sswc)                                                             !
!                                                                      !
!**********************************************************************!
subroutine wsparam(l_b_and_c,nupc,awl,kd,kx,nci,nfix,adp,topsl,sfc,sw,&
 sswc)
!**********************************************************************!
real(dp) :: nupc,awl(4),kd,kx,nci(4),nfix,adp(4),sfc(4),sw(4),sswc(4),&
 topsl,bc_wres,bc_wqs,ans1,bc_a,bc_l
logical :: l_b_and_c
!----------------------------------------------------------------------!

awl(1) = 0.05 ! Relative root density of the ith layer.
awl(2) = 0.05 !                  "
awl(3) = 0.03 !                  "
awl(4) = 0.01 !                  "
kd = 0.5    ! Fraction of excess water flowing to deep storage.
kx = 0.2    !                   "                 stream.
nci(1) = 0.04       ! NC ratio of litter pools.
nci(2) = 0.06666667 !          "
nci(3) = 0.04       !          "
nci(4) = 0.06666667 !          "
nupc = 1.0   ! Amount of nitrogen uptake ratio to litter.
nfix = 0.5   ! Soil deposition plus N fixation (g/m2/y).

if (l_b_and_c) then

!----------------------------------------------------------------------!
!                        Brookes Corey                                 !
!----------------------------------------------------------------------!
bc_wres = 0.0
bc_wqs = 1.0-ssp%bulk/2.65
ans1 = -0.475662 +0.005165*ssp%sand + 0.002066*ssp%silt - &
 0.023327*ssp%clay - 0.040628*ssp%bulk + 0.006824*ssp%orgc - &
 0.000136162*ssp%sand**2 - 0.0000954716*ssp%silt**2 + &
 0.000298891*ssp%clay**2 - 0.0637*ssp%bulk**2 - &
 0.00031679*ssp%orgc**2 + 0.0000010388*ssp%sand**3 + &
 0.0000010299*ssp%silt**3 - 0.0000013251*ssp%clay**3

bc_l = 10**ans1

ans1 = 0.4104 + 0.002684*ssp%sand + 0.006354*ssp%silt + &
 0.17766*ssp%bulk + 0.00013855*ssp%clay**2

bc_a = 1/10**ans1

ssp%wilt = bc_wres + (bc_wqs - bc_wres)*(15300.0*bc_a)**(-bc_l)
ssp%field = bc_wres + (bc_wqs - bc_wres)*(51.0*bc_a)**(-bc_l)
ssp%sat = 1.0-ssp%bulk/2.65

endif

!----------------------------------------------------------------------!
! Convert soil layer height into mm of water.                          !
!----------------------------------------------------------------------!
adp(1) = topsl*10.0
adp(2) = (ssp%soil_depth/3.0 - topsl)*10.0
adp(3) = ssp%soil_depth*10.0/3.0
adp(4) = ssp%soil_depth*10.0/3.0
sfc(1) = ssp%field*adp(1)
!      sfc(1) = adp(1)
sfc(2) = ssp%field*adp(2)
sfc(3) = ssp%field*adp(3)
sfc(4) = ssp%field*adp(4)
sw(1) = ssp%wilt*adp(1)
sw(2) = ssp%wilt*adp(2)
sw(3) = ssp%wilt*adp(3)
sw(4) = ssp%wilt*adp(4)
sswc(1) = ssp%sat*adp(1)
sswc(2) = ssp%sat*adp(2)
sswc(3) = ssp%sat*adp(3)
sswc(4) = ssp%sat*adp(4)


!**********************************************************************!
end subroutine wsparam
!**********************************************************************!

end module soil_methods
