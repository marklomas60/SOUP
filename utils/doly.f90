!> @brief Collection of subroutines to read in data.
!! @details
!! @author Mark Lomas
!! @date July 2016

module doly

use real_precision
use phenology_methods
use productivity_methods
use misc_values
use site_parameters
use tuning_parameters

implicit none

contains

!**********************************************************************!
!                                                                      !
!                          dolyday :: doly                             !
!                          ---------------                             !
!                                                                      !
!           (D)ynamic gl(O)ba(L) phytogeograph(Y) model                !
! SUBROUTINE DOLYDAY(tmp,prc,hum,ca,soilc,soiln,minn,adp,sfc,sw,sswc,  !
! awl,kd,kx,daygpp,resp_l,rlai,evap,tran,roff,interc,evbs,flow1,flow2, !
! year,mnth,day,pet,ht,thty_dys,ft,lmor_sc,nleaf,leaflitter,hrs,q,     !
! qdirect,qdiff,fpr,canga,gsn,rn,check_closure)                        !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Read internal parameters from "param.dat" file, and io
!! parameters from "misc_params.dat".
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine dolyday(tmp,prc,hum,ca,soilc,soiln,minn,adp,sfc,sw,sswc, &
 awl,kd,kx,daygpp,resp_l,rlai,evap,tran,roff,interc,evbs,flow1,flow2, &
 year,mnth,day,pet,ht,thty_dys,ft,lmor_sc,nleaf,leaflitter,hrs,q, &
 qdirect,qdiff,fpr,canga,gsn,rn,check_closure)
!**********************************************************************!
real(dp), parameter :: oi = 21000.0
real(dp) :: bb0,bbmax,bblim,sslim,tmp,prc,hum,ca,maxevap, &
 lairat,interc,evbs,t,suma,sum,laiinc, &
 soilc,soiln,rh,tk,rd,rn,pet,petmm,petv,q,qdiff,sla, &
 qdirect,hrs,vpd,canga,lam,rho,s,gam,amx,gsn,dp2,pet2,wtfc, &
 etwt,etmm,ftagh,soilw,maxc,soil2g,can2a,daygpp,daynpp, &
 gsum,sapresp,nppstorx,leafmol,ht,tf(13), &
 eemm,amax,respref,can2g,cangs,svp,rlai, &
 wtwp,p,nppstor2,canres,et,ee,flow1,flow2,tran,rem, &
 adp(4),evap,f2,sfc(4),sw(4),awl(4),kd, &
 kx,f3,minn,roff,windspeed,leaflit,resp_l, &
 npp_eff,eff_dec,eff_age,eff,nppstore,stemnpp, &
 rootnpp,leafnpp,tfscale,sswc(4),yield,yld,resp,fpr, &
 s1in,stemfr,lmor_sc(3600),nleaf,leaflitter,old_total_carbon, &
 total_carbon
integer leafls,stemls,rootls,bbm,ssm,sss,ftphen,c3c4,thty_dys,ft, &
 mnth,i,iter,ndsum(12),lai,day,year,bb,bbgs,ftdth,ss,dsbb, &
 chill,dschill,co
logical veg,check_closure
!----------------------------------------------------------------------!

co = ssp%cohort

!----------------------------------------------------------------------!
! Check carbon closure.
!----------------------------------------------------------------------!
if ((check_closure).and.(pft(co)%sla > 0.0)) then
  old_total_carbon = ssv(co)%lai%tot*12.0/pft(co)%sla/18.0 + &
 ssv(co)%nppstore(1) + ssv(co)%stem%tot + ssv(co)%root%tot + ssv(co)%bio(1) + &
 ssv(co)%bio(2)
endif
!----------------------------------------------------------------------!

wtwp     = ssp%wilt
wtfc     = ssp%field
sla      = pft(ft)%sla
c3c4     = pft(ft)%c3c4
ftphen   = pft(ft)%phen
ftagh    = pft(ft)%crop
ftdth    = pft(ft)%d2h
leafls   = pft(ft)%lls
stemls   = pft(ft)%sls
rootls   = pft(ft)%rls
bbm      = pft(ft)%bbmem
bb0      = pft(ft)%bb0
bbmax    = pft(ft)%bbmax
bblim    = pft(ft)%bblim
ssm      = pft(ft)%senm
sss      = pft(ft)%sens
sslim    = pft(ft)%senlim
lairat   = pft(ft)%lrat

leaflitter = 0.0

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

!----------------------------------------------------------------------!
! Data input for number of days for each month of the year 'ndays',    !
! the mid julian day of each month 'ndsum'.                            !
!----------------------------------------------------------------------!
data ndsum/16,46,75,106,136,167,197,228,259,289,320,350/

p = 101325.0

rlai = ssv(ssp%cohort)%lai%tot
rem = rlai - int(rlai)
lai = int(rlai) + 1

!----------------------------------------------------------------------!
! INITIALISATIONS.                                                     !
!----------------------------------------------------------------------!
! Night-time respiration ? CHECK
rd = 0.82e-6
! maxc, moisture response parameter
maxc = 690.0/622.6
! Leaf molecular weight.
if (sla>0.0) then
  leafmol = 1.0/(sla*18.0)
else
  leafmol =0.0
endif

! initialisation a faire a chaque fois
if (pft(ft)%phen == 0) then
   veg = .false.
else
   veg = .true.
endif

!----------------------------------------------------------------------!
! Water model variables and rain generater.                            !
!----------------------------------------------------------------------!
iter = 30

! unsused anymore... so not important if set to this value every day. 
! Better: should be remove from nppcalc
amx = -1000.0
amax = -1000.0

!----------------------------------------------------------------------!
soilw = ssv(ft)%soil_h2o(1) + ssv(ft)%soil_h2o(2) + &
 ssv(ft)%soil_h2o(3) + ssv(ft)%soil_h2o(4)
soil2g = soilw/(ssp%soil_depth*10.0)

if (soil2g>1.0) soil2g = 1.0

!----------------------------------------------------------------------!
! Bare ground IF statement                                             !
!----------------------------------------------------------------------!
if (.not.(veg)) rlai = 0.0
!----------------------------------------------------------------------!

t = tmp
tk = 273.0 + t
rh = hum
!      abmin = tmp(mnth,day)*1.29772 - 19.5362
if (rh>95.0)  rh=95.0
if (rh<30.0)  rh=30.0

! Sapwood respiration N.B. 8.3144 J/K/mol= Universal gas constant.
! sapresp = exp(21.6 - 5.367e4/(8.3144*tk))
sapresp = 0.4*exp(0.06*t)*0.1
sapresp = 0.4*exp(0.02*t)*0.35

! Converted to a monthly figure (N.B. 24 hour day).
! sapresp = (sapresp*3600.0*24.0*30.0)/1000000.0
! Use daily value
sapresp = (sapresp*3600.0*24.0)/1000000.0

! Totalled for the year.
if (soil2g>ssp%wilt) then
  respref = sapresp*tgp%p_resp*((soil2g - ssp%wilt)/&
 min(1.0,(wtfc - ssp%wilt)))**tgp%p_kgw
else
  respref = 0.0
endif

if (t<0.0) then
  respref = 0.0
endif

t = tmp
lam = 2500.0 - (2.367*t)
s = 48.7*exp(0.0532*t)
rn = 0.96*(q*1000000.0/4.0 + 208.0 + 6.0*t)
rn = rn*0.52
pet = (1.26*s*rn)/(s + 66.0)*tgp%p_pet
petmm = (pet*3600.0)/(lam*1000.0)*10.0
if (petmm>0.0) then
  petv = petmm/1.0
else
  petv = 0.0
endif
pet = petv

!----------------------------------------------------------------------!
! canga=k^2 u / (log[(z-d)/z0])^2
! k=von Karman constant. k=0.41
! z=reference height
! d=zero plane displacement
! z0=roughness length
!----------------------------------------------------------------------!
windspeed= 5.0 ! in m/s
canga = 0.168*windspeed/log((200.0 - 0.7*ht)/(0.1*ht))**2

!----------------------------------------------------------------------!

npp_eff = 0.0
sum = 0.0
eff_dec = 0.75
eff_age = 90.0
!DO i=1,leafls
!  IF (i.LE.eff_age) THEN
!    eff = 1.0
!  ELSE
!    eff = eff_dec + real(leafls - i)*(1.0 - eff_dec)/(real(leafls) - &
! eff_age + 1.0)
!  ENDIF
!  npp_eff = npp_eff +  leafv(i)*eff
!  sum = sum + leafv(i)
!ENDDO

!IF (abs(sum-rlai).gt.0.001) THEN
!  WRITE(*,*) 'leafv not = rlai ',SUM,RLAI,mnth,day,nppstore
!  STOP
!ENDIF
!IF (sum.GT.0.0) THEN
!  npp_eff = npp_eff/sum
!ELSE
!  npp_eff = 0.0
!ENDIF
npp_eff = 1.0

if (veg) then
  call NPPCALC(npp_eff,c3c4,maxc,soilc,soiln,minn,soil2g,wtwp,wtfc, &
 rd,rlai,t,rh,ca,oi,rn,qdirect,qdiff,can2a,can2g,canres,suma, &
 amx,amax,gsum,hrs,canga/1.3,p,day,nleaf,fpr)

!      write(*,*) 'o',canga,can2a,can2g,canres,suma,gsum

!----------------------------------------------------------------------!
! Calculate the canopy conductance.                                    !
!----------------------------------------------------------------------!
  cangs = 8.3144*tk/p*can2g
  cangs = 1.6*cangs             ! convert CO2 gs into H2O gs
  gsn = cangs
!----------------------------------------------------------------------!
else
  can2a = 0.0
  can2g = 0.0
  canres = 0.
  cangs = 0.0
  gsn = 0.0
  suma = 0.0
endif

call SUMA_ADD(suma)

dp2 = prc
pet2 = petv

pet = eemm
pet2 = pet

!----------------------------------------------------------------------!
! Set switch 'stseas' to find when a new growing season has started.   *
! 'stseas = 1' for the first day of a new growing season.              !
!----------------------------------------------------------------------!

soil2g = soilw/(ssp%soil_depth*10.0)

!----------------------------------------------------------------------!
!                         DAILY npp                                    !
!----------------------------------------------------------------------!
daynpp = can2a*3600.0*hrs/1000000.0 - canres*3600.0*(24.0-hrs)/1000000.0
daynpp = daynpp*12.0
resp_l = canres*3600.0*(24.0-hrs)*12.0/1000000.0
daygpp= can2a*3600.0*hrs*12.0/1000000.0

!----------------------------------------------------------------------!
! End of the DAILY LOOP.                                               !
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Check carbon closure.
!----------------------------------------------------------------------!
if ((check_closure).and.(pft(co)%sla > 0.0)) then
  total_carbon = ssv(co)%lai%tot*12.0/pft(co)%sla/18.0 + &
 ssv(co)%nppstore(1) + ssv(co)%stem%tot + ssv(co)%root%tot + ssv(co)%bio(1) + &
 ssv(co)%bio(2)
  if (abs(total_carbon-old_total_carbon) >  1.0e-3) then
    write(*,*) 'Breach of carbon closure in DOLYDAY:', &
 total_carbon-old_total_carbon,' g/m^2.'
    stop
  endif
endif

end subroutine dolyday





!**********************************************************************!
!                                                                      !
!                     evapotranspiration :: doly                       !
!                     --------------------------                       !
!                                                                      !
! SUBROUTINE evapotranspiration(t,rh,rn,canga,gsn,hrs,eemm,etmm)       !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Read internal parameters from "param.dat" file, and io
!! parameters from "misc_params.dat".
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine evapotranspiration(t,rh,rn,canga,gsn,hrs,eemm,etmm)
!**********************************************************************!
real(dp) :: t,rh,rn,canga,gsn,eemm,etmm,et,svp,vpd,lam,rho,s,gam,hrs, &
 maxevap,etwt,ee

!----------------------------------------------------------------------!
! Penman-Monteith equation for evapotranspiration.                     !
! Units of ET = W/m2  (CHECK) N.B. Conductances in moles.              !
!----------------------------------------------------------------------!
svp = 6.108*exp((17.269*t)/(237.3 + t))*100.0
vpd = (1.0 - rh/100.0)*svp
lam = 2500.0 - 2.367*t
rho = 1288.4 - 4.103*t
s = 48.7*exp(0.0532*t)
gam = 101325.0*1.012/(0.622*lam)

if ((ssv(ssp%cohort)%lai%tot>0.1).and.(msv%mv_soil2g>ssp%wilt)) then
  et = (s*rn + rho*1.012*canga*vpd)/(s + gam*(1.0 + canga/gsn))*tgp%p_et
! watch dog ajoute par ghislain 
  if (et<0.0) then
    et=0.0
  endif

  etwt = (et*3600.0*hrs)/lam
  etmm = etwt/1000.0
  if (etmm>0.0) then
    if (etmm/hrs>maxevap)  maxevap = etmm/hrs
  endif
else
  et = 0.0
  etwt = 0.0
  etmm = 0.0
endif

ee = (s*rn + rho*1.012*canga*vpd)/(s + gam)
eemm = (ee*3600.0*hrs)/(lam*1000.0)

! added by Ghislain 20/10/03
if (ee<0.0) then
! WRITE(*,*) 'ee is negativ ee=',ee,rn
  ee=0.0
  eemm=0.0
endif

end subroutine evapotranspiration



end module doly













