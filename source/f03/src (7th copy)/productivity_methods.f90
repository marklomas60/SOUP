module productivity_methods

use real_precision
use dims
use sunshade
use system_state
use pft_parameters
use site_parameters
use tuning_parameters

implicit none

contains

!**********************************************************************!
!                                                                      !
!                   nppcalc :: productivity_methods                    !
!                   -------------------------------                    !
!                                                                      !
! SUBROUTINE NPPCALC(npp_eff,c3,maxc,soilc,soiln,minn,soil2g,wtwp,     !
! wtfc,rd,rlai,t,rh,ca,oi,rn,qdirect,qdiff,can2a,can2g,canres,suma,amx,!
! amax,gsum,hrs,canga,p,day,nleaf_sum,fpr)                             !
!                                                                      !
!**********************************************************************!
subroutine nppcalc(npp_eff,c3,maxc,soilc,soiln,minn,soil2g,wtwp, &
 wtfc,rd,rlai,t,rh,ca,oi,rn,qdirect,qdiff,can2a,can2g,canres,suma,amx,&
 amax,gsum,hrs,canga,p,day,nleaf_sum,fpr)
!**********************************************************************!
real(dp) :: suma,sumd,rlai,soilc,soil2g,wtfc,nmult,npp_eff,soiln,y1,y0,&
 x1,x0, mmult,minn,nup,t,tk,up,kc,ko,tau,sum,nleaf_sum,rem,can(12),&
 vm(12),jmx(12),oi,canres,dresp(12),drespt(12),upt(12),nupw,vmx(12),rh,&
 c1,c2,canga,p,ga,nleaf,q,qdiff,qdirect,cs,coeff,can2a,can2g,rn,wtwp,&
 kg,maxc,rht,tpav,tppcv,tpgsv,ca,rd,tpaj,tppcj,tpgsj,hrs,a(12),gs(12),&
 ci(12),tpac4,tppcc4,tpgsc4,xvmax,asunlit,ashade,pcsunlit,pcshade, &
 gssunlit,gsshade,jsunlit(12),jshade(12),qsunlit(12),qshade(12), &
 fsunlit(12),fshade(12),ax,amax,amx,gsum,lyr,qt,t154,xdvy,px,oitau1, &
 oitau2,gsmin,x,y,dv,kcoiko,qt1,qt2,vmxp1,jp1,fpr,apar
! light-limited assimilation rate (j) and irradiance (q) for sunlit
! and shade
! fraction of sunlit and shade
integer :: i,lai,k,c3,day
real(dp) :: soilalbedo,leafalbedo,kbeam,kdiff,m,kbeamstar,canopyalbedo,&
 albedodiff,albedobeam
!----------------------------------------------------------------------!

apar=0.0

rem = rlai - int(rlai)
lai = int(rlai) + 1
tk = 273.0 + t

!----------------------------------------------------------------------!
! Set some parameters used in GOUDRIANS LAW routine (quicker to do it
! here).
!----------------------------------------------------------------------!
! soilalbedo: value from Wang
! leafalbedo: leaf transmittance and reflectance are equal. Then albedo
! is the double kbeam: extinction coefficent of the direct beam
! kdiff: extinction coefficient of diffuse radiations
! I found somewhere an approximate relation: kbeam=0.75*kdiff
! ... hum, I don't know if albedodiff is REAL*8ly different or not...
! so equal
!----------------------------------------------------------------------!
soilalbedo = 0.15
leafalbedo = 0.15
kbeam = 0.5
kdiff = 0.66
m = sqrt(1.0-leafalbedo)
kbeamstar = m*kbeam
canopyalbedo=2.0*kbeam/(kbeam+kdiff)*(1.0-m)/(1.0+m)
albedobeam=canopyalbedo+(soilalbedo-canopyalbedo)*exp(-kbeam*rlai)
albedodiff=albedobeam
!----------------------------------------------------------------------!
qt = 3.0
qt1 = qt**(t/10.0)/qt**(2.5)
qt2 = 100.0*qt**(30.0/10.0)/qt**(2.5)

suma = 0.0
sumd = 0.0
nleaf_sum = 0.0

kc = exp(35.8 - 80.5/(0.00831*tk))
ko = exp(9.6 - 14.51/(0.00831*tk))*1000.0
tau = exp(-3.949 + 28.99/(0.00831*tk))

! changed by Ghislain 08/10/03      IF (rlai.GT.0.1) THEN
q = qdiff + qdirect
if ((rlai>0.1).and.(q>0.0)) then

  if (soil2g>wtwp) then
    kg = maxc*((soil2g - wtwp)/(wtfc - wtwp))**tgp%p_kgw
  else
    kg = 0.0
  endif
  if (kg>maxc)  kg = maxc
  vmxp1 = npp_eff*kg**tgp%p_nu3

! Nitrogen uptake CHECK.
! nup = 120.0*(exp(-0.00008*soilc))
  nupw = tgp%p_nu1*(exp(tgp%p_nu2*soilc))*kg**tgp%p_nu3
  if (nupw<0.0) nupw = 0.0

! Nitrogen multiplier.
  nmult = soiln*tgp%p_nu4
  if (nmult>=1.0)  nmult = 1.0

  y1 = 1.3
  y0 = 1.0
  x1 = 50.0
  x0 = 10.0
  mmult = minn*(y1 - y0)/(x1 - x0) + (y0*x1 - x0*y1)/(x1 - x0)
  if (mmult>y1)  mmult = y1
  if (mmult<y0)  mmult = y0

  nup = nupw*nmult*mmult
  nup = nupw*nmult
  up = nup

  if (up<0.0)  up = 0.0

!----------------------------------------------------------------------!
! Total the Beer's Law values for each lai up to the current one.      !
!----------------------------------------------------------------------!
  coeff =-0.5
  sum = 0.0 
  do i=1,lai-1
    k = i - 1
    sum = sum + exp(coeff*real(k))
  enddo
  sum = sum + exp(coeff*real(lai - 1))*rem

!----------------------------------------------------------------------!
  do i=1,lai
    lyr = real(i)

    can(i) = exp(coeff*(lyr - 1.0))/sum

    upt(i) = up*can(i)

    nleaf = (upt(i))*tgp%p_nleaf
    if (i<lai) then
      nleaf_sum = nleaf_sum + nleaf
    else
      nleaf_sum = nleaf_sum + nleaf*rem
    endif

    dresp(i) = nleaf*tgp%p_dresp
    drespt(i) = dresp(i)

    if (kg>1.0e-10) then

!        vm(i) = nleaf*tgp%p_vm
!        jm(i) = (29.1 + 1.64*vm(i))

!        IF (jm(i).lt.0.0) jm(i) = 0.0
!        vm(i) = vm(i)/1000000.0
!        jm(i) = jm(i)/1000000.0

!        xx = (-(2.0*tgp%p_j2-150.0*tgp%p_j3)-(
!     &(2.0*tgp%p_j2-150.0*tgp%p_j3)**2.0-12.0*tgp%p_j3*
!     &(tgp%p_j1-50.0*tgp%p_j2+1875.0*tgp%p_j3))**0.5)/6.0/tgp%p_j1

!        IF (t.GT.xx) THEN
!          jmx(i) = jm(i)*(1.0 + tgp%p_j1*(t-25.0) +
!     &tgp%p_j2*(t - 25.0)**2 + tgp%p_j3*(t - 25.0)**3)
!        ELSE
!          jmx(i) = jm(i)*(1.0 + tgp%p_j1*(xx-25.0) +
!     &tgp%p_j2*(xx - 25.0)**2 + tgp%p_j3*(xx - 25.0)**3)
!        ENDIF

!        jmx(i) = jm(i)*(0.4 + t*0.6/30.0)
!        IF (jmx(i).GT.jm(i))  jmx(i) = jm(i)

!        yy = (-(2.0*tgp%p_v2-150.0*tgp%p_v3)-(
!     &(2.0*tgp%p_v2-150.0*tgp%p_v3)**2.0-12.0*tgp%p_v3*
!     &(tgp%p_v1-50.0*tgp%p_v2+1875.0*tgp%p_v3))**0.5)/6.0/tgp%p_v
!        IF (t.GT.yy) THEN
!          vmx(i) = vm(i)*(1.0 + tgp%p_v1*(t - 25.0) +
!     &tgp%p_v2*(t - 25.0)**2 + tgp%p_v3*(t - 25.0)**3)
!        ELSE
!          vmx(i) = vm(i)*(1.0 + tgp%p_v1*(yy - 25.0) +
!     &tgp%p_v2*(yy - 25.0)**2 + tgp%p_v3*(yy - 25.0)**3)
!        ENDIF

!        vmx(i) = vm(i)*(0.4 + t*0.6/30.0)
!        IF (vmx(i).GT.vm(i))  vmx(i) = vm(i)

!        xx = (0.25 + (1.0-0.25)/25.0*t)**1.5

      vm(i) = nleaf*tgp%p_vm
      vmx(i)=vm(i)*qt1
      if (t>30.0) vmx(i)=qt2
      vmx(i) = vmx(i)*vmxp1

! vmx(i) = vmx(i)*(1.0+0.25*(1.0-ca/35.0))

      jmx(i) = (29.1 + 1.64*vmx(i))
      vmx(i) = vmx(i)*1.0e-6
      jmx(i) = jmx(i)*1.0e-6

      rd = vmx(i)*0.02

! change Ghislain
!        irc = q*exp(-0.5*lyr)

!        CALL BEERSLAW(lyr-0.5,rlai,q,0.0,
!     &       fsunlit(i),qsunlit(i), fshade(i),qshade(i))

      call GOUDRIAANSLAW2(lyr-0.5,rlai,qdirect,qdiff,fsunlit(i),&
 qsunlit(i),fshade(i),qshade(i),soilalbedo,leafalbedo,kbeam,kdiff,m,&
 kbeamstar,canopyalbedo,albedobeam,albedodiff)

!        j(i) = 0.24*irc/(1.0 + (0.24**2)*(irc**2)/(jmx(i)**2))
!     &**0.5

!     calculate separatly the sunlit and shaded light limited 
!     assimilation rate
  jp1 = jmx(i)**2
  jsunlit(i)=0.24*qsunlit(i)/(1.0+0.0576*(qsunlit(i)**2)/jp1)**0.5
  jshade(i)=0.24*qshade(i)/(1.0+0.0576*(qshade(i)**2)/jp1)**0.5

!  jsunlit(i)=0.24*qsunlit(i)/(1.0+(0.24**2)*(qsunlit(i)**2)/&
! (jmx(i)**2))**0.5
!  jshade(i)=0.24*qshade(i)/(1.0+(0.24**2)*(qshade(i)**2)/&
! (jmx(i)**2))**0.5

    else

      call GOUDRIAANSLAW2(lyr-0.5,rlai,qdirect,qdiff,fsunlit(i),&
 qsunlit(i),fshade(i),qshade(i),soilalbedo,leafalbedo,kbeam,kdiff,m,&
 kbeamstar,canopyalbedo,albedobeam,albedodiff)

      vmx(i) = 0.0
      jsunlit(i) = 0.0
      jshade(i) = 0.0

    endif

  enddo

!----------------------------------------------------------------------!
!     fAPAR calculation added here from gh nppcalc.f amf251105
!----------------------------------------------------------------------!
  call GOUDRIAANSLAW2(rlai,rlai,qdirect,qdiff,fsunlit(lai+1),&
 qsunlit(lai+1),fshade(lai+1),qshade(i),soilalbedo,leafalbedo,kbeam,&
 kdiff,m,kbeamstar,canopyalbedo,albedobeam,albedodiff)

  apar = apar + rem*(fsunlit(lai+1)*qsunlit(lai+1)+&
 fshade(lai+1)*qshade(lai+1))
  fpr=apar/(qdiff+qdirect)
      
!----------------------------------------------------------------------!
! removed by Ghislain 24/11/03
! the output xvmx, xjmx and xj are not used anymore
!      CALL XVJMAX(xvmx,xjmx,xj,t,q,sum,nup,oi,xdresp)
!----------------------------------------------------------------------!
  canres = 0.0
  do i=1,lai
    if (i<lai) then
      canres = canres + dresp(i)
    else
      canres = canres + dresp(i)*rem
    endif
  enddo

else
  do i=1,lai
    can(i) = 0.0
    upt(i) = 0.0
    dresp(i) = 0.0
    drespt(i) = dresp(i)
    vm(i) = 0.0
    jmx(i) = 0.0
    vmx(i) = 0.0

    jshade(i) = 0.0
    jsunlit(i) = 0.0
    fsunlit(i) = 0.0
    fshade(i) = 0.0
    qsunlit(i) = 0.0
    qshade(i) = 0.0
  enddo
  canres = 0.0
endif

c1 = 142.4 - 4.8*t
if (c1>80.0)  c1 = 80.0
if (c1<8.0)  c1 = 8.0
c2 = 12.7 - 0.207*t
if (c2>10.0)  c2 = 10.0
if (c2<6.9)  c2 = 6.9

!----------------------------------------------------------------------!
! Assimilation calculations using subroutines ASSVMAX and ASSJ.        !
!----------------------------------------------------------------------!
can2a = 0.0
can2g = 0.0

t154 = 1.54*t
px = 1.0e-6*p
oitau1 = oi/tau
oitau2 = 0.5*oi/tau
gsmin = 0.005
x = 3.0 + 0.2*t
if (x<0.25) x = 0.25
y = 1.5
dv = 0.6108*exp(17.269*t/(237.3 + t))*(1.0 - rh/100.0)
xdvy = x/(1.0 + dv/y)
kcoiko = kc*(1.0 + oi/ko)

do i=1,lai

  ssp%lai = i
  dresp(i) = drespt(i)

  if (rlai>0.1) then
    if ((t>0.0).and.(rn>0.0).and.(soil2g>wtwp).and.&
 (fshade(i)+fsunlit(i)>0.1)) then
      ga = canga/(8.3144*tk/p)/rlai

!----------------------------------------------------------------------!
! Only recompute assimilation every few days to speed things up, but   !
! make sure it is calculated on the first day of any growing season.   !
!----------------------------------------------------------------------!
      rht = rh

      if (c3==1) then
        call ASSVMAX3(tpav,tppcv,cs,tpgsv,oi,ca,vmx(i),rh,kg,ko,kc,&
 tau,rd,ga,t,p,t154,xdvy,px,oitau1,oitau2,gsmin,kcoiko)
!        CALL ASSVMAX(tpav,tppcv,cs,tpgsv,oi,ca,vmx(i),rh,kg,ko,kc,tau,&
! rd,ga,t,p)
!        print*,tpav
!     make the sum of assimilation, internal concentration and stomatal
!     conductance.
        a(i)=0.0
        ci(i)=0.0
        gs(i)=0.0

        if (fshade(i)>0.01) then
          call ASSJ3(tpaj,tppcj,tpgsj,oi,ca,rh,kg,tau,rd,ga,t,p,&
 jshade(i),t154,xdvy,px,oitau1,oitau2,gsmin)
!          CALL ASSJ(tpaj,tppcj,tpgsj,oi,ca,rh,kg,tau,rd,ga,t,p,&
! jshade(i))
!          print*,'fsh ',tpaj
          call LIMITATIONPROCESS(tpav,tppcv,tpgsv,tpaj,tppcj,tpgsj,&
 ashade,pcshade,gsshade)
          a(i)=a(i)+fshade(i)*ashade
          ci(i)=ci(i)+fshade(i)*pcshade
          gs(i)=gs(i)+fshade(i)*gsshade
        endif

        if (fsunlit(i)>0.01) then
          call ASSJ4(tpaj,tppcj,tpgsj,oi,ca,rh,kg,tau,rd,ga,t,p,&
 jsunlit(i),t154,xdvy,px,oitau1,oitau2,gsmin)
!          CALL ASSJ(tpaj,tppcj,tpgsj,oi,ca,rh,kg,tau,rd,ga,t,p,&
! jsunlit(i))
!          print*,'fsun ',tpaj
          call LIMITATIONPROCESS(tpav,tppcv,tpgsv,tpaj,tppcj,tpgsj,&
 asunlit,pcsunlit,gssunlit)
          a(i)=a(i)+fsunlit(i)*asunlit
          ci(i)=ci(i)+fsunlit(i)*pcsunlit
          gs(i)=gs(i)+fsunlit(i)*gssunlit
        endif

      else

        tpac4=0.0
        tppcc4=0.0
        tpgsc4=0.0

        if (fshade(i)>0.01) then
          call ASSC4(ashade,pcshade,gsshade,ca,xvmax,rht,t,qshade(i),&
 p,upt(i))
          tpac4=tpac4+fshade(i)*ashade
          tppcc4=tppcc4+fshade(i)*pcshade
          tpgsc4=tpgsc4+fshade(i)*gsshade
        endif
        if (fsunlit(i)>0.01) then
          call ASSC4(asunlit,pcsunlit,gssunlit,ca,xvmax,rht,t,&
 qsunlit(i),p,upt(i))
          tpac4=tpac4+fsunlit(i)*asunlit
          tppcc4=tppcc4+fsunlit(i)*pcsunlit
          tpgsc4=tpgsc4+fsunlit(i)*gssunlit
        endif

        a(i) = tpac4
        gs(i) = tpgsc4
        ci(i) = tppcc4

        if (gs(i)<0.0)  gs(i) = 0.0
        gs(i) = gs(i)/1000.0

      endif
!----------------------------------------------------------------------!
! End of assimilation calculation 'IF' statement.                      !
!----------------------------------------------------------------------!
    else
      a(i) = 0.0
      gs(i) = 0.005
      ci(i) = ca
    endif
!----------------------------------------------------------------------!
! End of 't & rn > 0' 'IF' statement.                                  !
!----------------------------------------------------------------------!
  else
    a(i) = 0.0
    gs(i) = 0.0
    dresp(i) = 0.0
    ci(i) = ca
  endif
!----------------------------------------------------------------------!
! End of 'leaf' 'IF' statement.                                        !
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! At the start of a season recompute dark resperation.                 !
!----------------------------------------------------------------------!
!        IF (stseas.EQ.1)  dresp(i) = drespt(i)

!----------------------------------------------------------------------!
! Cumulate daily assimilation of the last lai layer, 'suma'.           !
!----------------------------------------------------------------------!
  if (i==lai) then
    if (lai>1) then
      suma = suma + (rem*a(lai) + (1.0 - rem)*a(lai - 1))* &
 3600.0*hrs - (rem*dresp(lai) + (1.0 - rem)*dresp(lai - &
 1))*3600.0*(24.0 - hrs)
      sumd = sumd + (rem*dresp(lai) + (1.0 - rem)*dresp(lai - &
 1))*3600.0*(24.0 - hrs)
    else
      suma = suma + rem*a(lai)*3600.0*hrs - rem*dresp(lai)* &
 3600.0*(24.0 - hrs)
      sumd = sumd + rem*dresp(lai)*3600.0*(24.0 - hrs)
      suma = 0.0
    endif
  endif
!----------------------------------------------------------------------!

  if (i<lai) then
    can2a = can2a + a(i)
    can2g = can2g + gs(i)
  else
    can2a = can2a + a(i)*rem
    can2g = can2g + gs(i)*rem
  endif

enddo
!----------------------------------------------------------------------!
! End of LAI loop.                                                     !
!----------------------------------------------------------------------!

if (lai>1) then
  gsum = gs(lai - 1)*(1.0 - rem) + gs(lai)*rem
else
  gsum = gs(1)
endif

if ((rlai>0.0).and.(t>0.0).and.(rn>0.0)) then
  if (mod(day-1,1)==0) then

!          CALL XAMAX(ax,t,oi,xvmx,rd,xj,xdresp,ci(1))
! Temporary change by Ghislain
     ax=0
! end if change
    if (rlai<-1.0) then
      if (ax*rem>amax)  then
        amax = a(1)*rem
      endif
      if (ax*rem>amx)  amx = ax*rem
    else
      if (ax>amax)  then
        amax = a(1)
      endif
      if (ax>amx)  amx = ax
    endif

  endif
endif

suma = suma/1000000.0

!**********************************************************************!
end subroutine nppcalc
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!              limitationprocess :: productivity_methods               !
!              -----------------------------------------               !
!                                                                      !
!     determin which process light, or Rubisco is the limiting         !
!     factor and return the assimilation, partial presure and          !
!     stomatal conductance corresponding to the limiting process       !
!                                                                      !
! SUBROUTINE limitationprocess(av,pcv,gsv,aj,pcj,gsj,a,pc,gs)          !
!                                                                      !
!**********************************************************************!
subroutine limitationprocess(av,pcv,gsv,aj,pcj,gsj,a,pc,gs)
!**********************************************************************!
real(dp) :: av,pcv,gsv,aj,pcj,gsj,a,pc,gs
!----------------------------------------------------------------------!

if (av<aj) then
   a = av
   gs = gsv
   pc = pcv
else
   a = aj
   gs = gsj
   pc = pcj
endif

!**********************************************************************!
end subroutine limitationprocess
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                   assvmax :: productivity_methods                    !
!                   -------------------------------                    !
! ASSVMAX computes the assimilation rate by forming a cubic equation   !
! from the equations defining assimilation rate, stomatal conductance, !
! CO2 partial pressure (assimilation diffusion equaiton) and the       !
! carboxylation rate. The cubic equation is solved and the largest     !
! real root is assumed to be the assimilation rate. The variables here !
! correspond to the Doly paper, in the main program they are such that !
!        po,pa,vmax,g0,g1,rh,ko,kc,kg,tau,rd                           !
!     = oi,ca,vmx(i),c1,c2,rh,ko,kc,kg,tau,rd.                         !
!                                                                      !
! SUBROUTINE assvmax(a,ci,cs,gs,oi,ca,vmx,rh,kg,ko,kc,tau,rd,ga,t,p)   !
!                                                                      !
!**********************************************************************!
subroutine assvmax(a,ci,cs,gs,oi,ca,vmx,rh,kg,ko,kc,tau,rd,ga,t,p)
!**********************************************************************!
real(dp) :: a,ci,gs,oi,ca,vmx,rh,ko,kc,tau,rd,fa0,faf,gam,step, &
 p,ga,x,y,dv,t,a0,cs,af,fa,kg,ax,bx,cx,gsmin
integer :: i
!----------------------------------------------------------------------!

gsmin = 0.005

x = 3.0 + 0.2*t
if (x<0.25) x = 0.25
y = 1.5

dv = 0.6108*exp(17.269*t/(237.3 + t))*(1.0 - rh/100.0)

if (5>30) then
do i=0,200
  a = real(i)/10.0
  cs = ca - a*1.0e-6*p/ga
  gs = gsmin + (x*a/(1.0 + dv/y)/(cs*10.0 - 1.54*t))*kg
  ci = cs - a*1.0e-6*p/gs
  gam = 1.0 - 0.5*oi/tau/ci
  fa0 = a - (gam*vmx*ci/(ci + kc*(1.0 + oi/ko)) - rd)*1.0e6
  write(*,'(7f10.3)') a,fa0,cs,gs,ci,gam,(x*a/(1.0 + dv/y)/&
 (cs*10.0 - 1.54*t))
enddo
endif

step = 1.0
a = 0.0
cs = ca - a*1.0e-6*p/ga
gs = gsmin + (x*a/(1.0 + dv/y)/(cs*10.0 - 1.54*t))*kg
ci = cs - a*1.0e-6*p/gs
gam = 1.0 - 0.5*oi/tau/ci
fa0 = a - (gam*vmx*ci/(ci + kc*(1.0 + oi/ko)) - rd)*1.0e6

a0 = a

50    continue
  a = a + step
  cs = ca - a*1.0e-6*p/ga
  gs = gsmin + (x*a/(1.0 + dv/y)/(cs*10.0 - 1.54*t))*kg
  ci = cs - a*1.0e-6*p/gs
  gam = 1.0 - 0.5*oi/tau/ci
  faf = a - (gam*vmx*ci/(ci + kc*(1.0 + oi/ko)) - rd)*1.0e6
if (faf<0.0) then
  fa0 = faf
  a0 = a
  goto 50
endif
af = a

if (af>100) then
  a = 0.0
  cs = ca - a*1.0e-6*p/ga
  gs = gsmin + (x*a/(1.0 + dv/y)/(cs*10.0 - 1.54*t))*kg
  ci = cs - a*1.0e-6*p/gs
  gam = 1.0 - 0.5*oi/tau/ci
else
  a = (a0 + af)/2.0
  cs = ca - a*1.0e-6*p/ga
  gs = gsmin + (x*a/(1.0 + dv/y)/(cs*10.0 - 1.54*t))*kg
  ci = cs - a*1.0e-6*p/gs
  gam = 1.0 - 0.5*oi/tau/ci
  fa = a - (gam*vmx*ci/(ci + kc*(1.0 + oi/ko)) - rd)*1.0e6

  bx = ((a+a0)*(faf-fa)/(af-a)-(af+a)*(fa-fa0)/(a-a0))/(a0-af)
  ax = (faf-fa-bx*(af-a))/(af**2-a**2)
  cx = fa-bx*a-ax*a**2

  if (abs(ax)>0.0) then
    if (bx**2-4.0*ax*cx>0.0) then
      a = (-bx+(bx**2-4.0*ax*cx)**0.5)/(2.0*ax)
  if (a>af)  a = (-bx-(bx**2-4.0*ax*cx)**0.5)/(2.0*ax)
    else
      a = 0.0
    endif
  else
    if (abs(bx)>0.0) then
      a =-cx/bx
    else
      a = 0.0
    endif
  endif

endif

  cs = ca - a*1.0e-6*p/ga
  gs = gsmin + (x*a/(1.0 + dv/y)/(cs*10.0 - 1.54*t))*kg
  ci = cs - a*1.0e-6*p/gs
  gam = 1.0 - 0.5*oi/tau/ci

!**********************************************************************!
end subroutine assvmax
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                   assvmax3 :: productivity_methods                   !
!                   --------------------------------                   !
! ! ASSVMAX computes the assimilation rate by forming a cubic equation !
! from the equations defining assimilation rate, stomatal conductance, !
! CO2 partial pressure (assimilation diffusion equaiton) and the       !
! carboxylation rate. The cubic equation is solved and the largest     !
! real root is assumed to be the assimilation rate. The variables here !
! correspond to the Doly paper, in the main program they are such that !
!        po,pa,vmax,g0,g1,rh,ko,kc,kg,tau,rd                           !
!     = oi,ca,vmx(i),c1,c2,rh,ko,kc,kg,tau,rd.                         !
!                                                                      !
! SUBROUTINE assvmax3(a,ci,cs,gs,oi,ca,vmx,rh,kg,ko,kc,tau,rd,ga,t,p,  !
! t154,xdvy,px,oitau1,oitau2,gsmin,kcoiko)                             !
!                                                                      !
!**********************************************************************!
subroutine assvmax3(a,ci,cs,gs,oi,ca,vmx,rh,kg,ko,kc,tau,rd,ga,t,p,&
 t154,xdvy,px,oitau1,oitau2,gsmin,kcoiko)
!**********************************************************************!
real(dp) :: a,ci,gs,oi,ca,vmx,rh,ko,kc,tau,rd,fa0,faf,gam,step, &
 p,ga,x,y,dv,t,a0,cs,af,fa,kg,ax,bx,cx,gsmin,t154,xdvy,px,oitau1,&
 oitau2,kcoiko
!----------------------------------------------------------------------!

if (ssv(ssp%cohort)%assv(ssp%lai) > 0) then
  a = ssv(ssp%cohort)%assv(ssp%lai)
  cs = ca - a*px/ga
  gs = gsmin + (xdvy*a/(cs*10.0 - t154))*kg
  ci = cs - a*px/gs
  gam = 1.0 - oitau2/ci
  a = (gam*vmx*ci/(ci + kcoiko) - rd)*1.0e6
!  print*,'it ',ssv(ssp%cohort)%assv(ssp%lai),a
else
  call ASSVMAX(a,ci,cs,gs,oi,ca,vmx,rh,kg,ko,kc,tau,rd,ga,t,p)
endif
ssv(ssp%cohort)%assv(ssp%lai) = a

!**********************************************************************!
end subroutine assvmax3
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                   assj :: productivity_methods                       !
!                   ----------------------------                       !
!                                                                      !
!! ASSVMAX computes the assimilation rate by forming a cubic equation  !
! from the equations defining assimilation rate, stomatal conductance, !
! CO2 partial pressure (assimilation diffusion equaiton) and the       !
! carboxylation rate. The cubic equation is solved and the largest     !
! real root is assumed to be the assimilation rate. The variables here !
! correspond to the Doly paper, in the main program they are such that !
!        po,pa,vmax,g0,g1,rh,ko,kc,kg,tau,rd                           !
!     = oi,ca,vmx(i),c1,c2,rh,ko,kc,kg,tau,rd.                         !
!                                                                      !
! SUBROUTINE assj(a,ci,gs,oi,ca,rh,kg,tau,rd,ga,t,p,j)                 !
!                                                                      !
!**********************************************************************!
subroutine assj(a,ci,gs,oi,ca,rh,kg,tau,rd,ga,t,p,j)
!**********************************************************************!
real(dp) :: a,ci,gs,oi,ca,rh,tau,rd,j,gam,fa0,faf,a0,af,step,p,ga,x,y,&
 dv,t,cs,fa,kg,ax,bx,cx,gsmin
!----------------------------------------------------------------------!

gsmin = 0.005

x = 3.0 + 0.2*t
if (x<0.25) x = 0.25
y = 1.5

dv = 0.6108*exp(17.269*t/(237.3 + t))*(1.0 - rh/100.0)

step = 1.0
a = 0.0
cs = ca - a*1.0e-6*p/ga
gs = gsmin + (x*a/(1.0 + dv/y)/(cs*10.0 - 1.54*t))*kg
ci = cs - a*1.0e-6*p/gs
gam = 1.0 - 0.5*oi/tau/ci
fa0 = a - (gam*j*ci/4.0/(ci + oi/tau) - rd)*1.0e6

a0 = a

50    continue
  a = a + step
  cs = ca - a*1.0e-6*p/ga
  gs = gsmin + (x*a/(1.0 + dv/y)/(cs*10.0 - 1.54*t))*kg
  ci = cs - a*1.0e-6*p/gs
  gam = 1.0 - 0.5*oi/tau/ci
  faf = a - (gam*j*ci/4.0/(ci + oi/tau) - rd)*1.0e6
if (faf<0.0) then
  fa0 = faf
  a0 = a
  goto 50
endif
af = a

if (af>100) then
  a = 0.0
  cs = ca - a*1.0e-6*p/ga
  gs = gsmin + (x*a/(1.0 + dv/y)/(cs*10.0 - 1.54*t))*kg
  ci = cs - a*1.0e-6*p/gs
  gam = 1.0 - 0.5*oi/tau/ci
else

  a = (a0 + af)/2.0
  cs = ca - a*1.0e-6*p/ga
  gs = gsmin + (x*a/(1.0 + dv/y)/(cs*10.0 - 1.54*t))*kg
  ci = cs - a*1.0e-6*p/gs
  gam = 1.0 - 0.5*oi/tau/ci
  fa = a - (gam*j*ci/4.0/(ci + oi/tau) - rd)*1.0e6

  bx = ((a+a0)*(faf-fa)/(af-a)-(af+a)*(fa-fa0)/(a-a0))/(a0-af)
  ax = (faf-fa-bx*(af-a))/(af**2-a**2)
  cx = fa-bx*a-ax*a**2

  if (abs(ax)>0.0) then
    if (bx**2-4.0*ax*cx>0.0) then
      a = (-bx+(bx**2-4.0*ax*cx)**0.5)/(2.0*ax)
  if (a>af)  a = (-bx-(bx**2-4.0*ax*cx)**0.5)/(2.0*ax)
    else
      a = 0.0
    endif
  else
    if (abs(bx)>0.0) then
      a =-cx/bx
    else
      a = 0.0
    endif
  endif

endif

cs = ca - a*1.0e-6*p/ga
gs = gsmin + (x*a/(1.0 + dv/y)/(cs*10.0 - 1.54*t))*kg
ci = cs - a*1.0e-6*p/gs
gam = 1.0 - 0.5*oi/tau/ci

!**********************************************************************!
end subroutine assj
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                      assj3 :: productivity_methods                   !
!                      -----------------------------                   !
!                                                                      !
! ASSVMAX computes the assimilation rate by forming a cubic equation   !
! from the equations defining assimilation rate, stomatal conductance, !
! CO2 partial pressure (assimilation diffusion equaiton) and the       !
! carboxylation rate. The cubic equation is solved and the largest     !
! real root is assumed to be the assimilation rate. The variables here !
! correspond to the Doly paper, in the main program they are such that !
!        po,pa,vmax,g0,g1,rh,ko,kc,kg,tau,rd                           !
!     = oi,ca,vmx(i),c1,c2,rh,ko,kc,kg,tau,rd.                         !
!                                                                      !
! SUBROUTINE assj3(a,ci,gs,oi,ca,rh,kg,tau,rd,ga,t,p,j,t154,xdvy,px,   !
! oitau1,oitau2,gsmin)                                                                     !
!                                                                      !
!**********************************************************************!
subroutine assj3(a,ci,gs,oi,ca,rh,kg,tau,rd,ga,t,p,j,t154,xdvy,px,&
 oitau1,oitau2,gsmin)
!**********************************************************************!
real(dp) :: a,ci,gs,oi,ca,rh,tau,rd,j,gam,fa0,faf,a0,af,step,p,ga,x,y,&
 dv,t,cs,fa,kg,ax,bx,cx,gsmin,t154,xdvy,px,oitau1,oitau2
!----------------------------------------------------------------------!

if (ssv(ssp%cohort)%assj(ssp%lai,1) > 0) then
  a = ssv(ssp%cohort)%assj(ssp%lai,1)
  cs = ca - a*px/ga
  gs = gsmin + (xdvy*a/(cs*10.0 - t154))*kg
  ci = cs - a*px/gs
  gam = 1.0 - oitau2/ci
  a = (gam*j*ci*0.25/(ci + oitau1) - rd)*1.0e6
!  print*,'it ',ssv(ssp%cohort)%assj(ssp%lai,1),a
else
  call ASSJ(a,ci,gs,oi,ca,rh,kg,tau,rd,ga,t,p,j)
endif
ssv(ssp%cohort)%assj(ssp%lai,1) = a

!**********************************************************************!
end subroutine assj3
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                     assj4 :: productivity_methods                    !
!                     -----------------------------                    !
!                                                                      !
! ASSVMAX computes the assimilation rate by forming a cubic equation   !
! from the equations defining assimilation rate, stomatal conductance, !
! CO2 partial pressure (assimilation diffusion equaiton) and the       !
! carboxylation rate. The cubic equation is solved and the largest     !
! real root is assumed to be the assimilation rate. The variables here !
! correspond to the Doly paper, in the main program they are such that !
!        po,pa,vmax,g0,g1,rh,ko,kc,kg,tau,rd                           !
!     = oi,ca,vmx(i),c1,c2,rh,ko,kc,kg,tau,rd.                         !
!                                                                      !
! SUBROUTINE assj4(a,ci,gs,oi,ca,rh,kg,tau,rd,ga,t,p,j,t154,xdvy,px,&  !
! oitau1,oitau2,gsmin)                                                 !
!                                                                      !
!**********************************************************************!
subroutine assj4(a,ci,gs,oi,ca,rh,kg,tau,rd,ga,t,p,j,t154,xdvy,px,&
 oitau1,oitau2,gsmin)
!**********************************************************************!
real(dp) :: a,ci,gs,oi,ca,rh,tau,rd,j,gam,fa0,faf,a0,af,step,p,ga,x,y,dv,t,cs, &
 fa,kg,ax,bx,cx,gsmin,t154,xdvy,px,oitau1,oitau2
!----------------------------------------------------------------------!

if (ssv(ssp%cohort)%assj(ssp%lai,2) > 0) then
  a = ssv(ssp%cohort)%assj(ssp%lai,2)
  cs = ca - a*px/ga
  gs = gsmin + (xdvy*a/(cs*10.0 - t154))*kg
  ci = cs - a*px/gs
  gam = 1.0 - oitau2/ci
  a = (gam*j*ci*0.25/(ci + oitau1) - rd)*1.0e6
!  print*,'it ',ssv(ssp%cohort)%assj(ssp%lai,2),a
else
  call ASSJ(a,ci,gs,oi,ca,rh,kg,tau,rd,ga,t,p,j)
endif
ssv(ssp%cohort)%assj(ssp%lai,2) = a

!**********************************************************************!
end subroutine assj4
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                  assvmax2 :: productivity_methods                    !
!                  --------------------------------                    !
!                                                                      !
! ASSVMAX computes the assimilation rate by forming a cubic equation   !
! from the equations defining assimilation rate, stomatal conductance, !
! CO2 partial pressure (assimilation diffusion equaiton) and the       !
! carboxylation rate. The cubic equation is solved and the largest     !
! real root is assumed to be the assimilation rate. The variables here !
! correspond to the Doly paper, in the main program they are such that !
!        po,pa,vmax,g0,g1,rh,ko,kc,kg,tau,rd                           !
!     = oi,ca,vmx(i),c1,c2,rh,ko,kc,kg,tau,rd.                         !
!                                                                      !
! SUBROUTINE assvmax2(a,w,pc,gs,po,pa,vmax,g0,g1,rh,ko,kc,kg,tau,&     !
! rd,acheck)                                                           !
!                                                                      !
!**********************************************************************!
subroutine assvmax2(a,w,pc,gs,po,pa,vmax,g0,g1,rh,ko,kc,kg,tau,rd,&
 acheck)
!**********************************************************************!
real(dp) :: a,w,pc,gs,po,pa,vmax,g0,g1,rh,ko,kc,kg,tau,rd,acheck,r2q3, &
 ax,bx,cx,dx,ex,fx,gx,hx,p,q,r,qt,rt,th,a1,a2,a3,pi,at,bt,sc,arg1
!----------------------------------------------------------------------!

pi = atan(1.0)*4.0
sc = 1000000.0
!----------------------------------------------------------------------!
! 'a2,b2,c2,d2' are such that 'w=(a2a + b2)/(c2a + d2)'.               !
!----------------------------------------------------------------------!
ax = ko*pa*vmax*(kg*g1*rh - 160.0)
bx = kg*g0*pa**2.0*vmax*ko
cx =-160.0*pa*ko + kc*kg*ko*g1*rh + pa*kg*ko*g1*rh + kc*kg*po*g1*rh
dx = pa**2.0*kg*ko*g0 + kc*kg*po*g0*pa + kc*kg*ko*g0*pa

!----------------------------------------------------------------------!
! 'ex,fx,gx,hx' are such that '0.5po/(taupc)-1=(exa+fx)/(gxa+hx)'.     !
!----------------------------------------------------------------------!
ex =-tau*pa*kg*g1*rh + 160.0*tau*pa + 0.5*po*kg*g1*rh
fx =-tau*pa**2*kg*g0 + 0.5*po*kg*g0*pa
gx = tau*pa*(kg*g1*rh - 160.0)
hx = tau*pa**2.0*kg*g0

!----------------------------------------------------------------------!
! Cubic coefficients 'p,q,r'.                                          !
!----------------------------------------------------------------------!
p = sc*rd + (sc*ax*ex + cx*hx + dx*gx)/(cx*gx)
q = (sc*rd*cx*hx + sc*rd*dx*gx + sc*ax*fx + sc*bx*ex + dx*hx)/(cx*gx)
r = (sc*rd*dx*hx + sc*bx*fx)/(cx*gx)

!----------------------------------------------------------------------!
! Sove the cubic equaiton to find 'a'.                                 !
!----------------------------------------------------------------------!
Qt = (p**2 - 3.0*q)/9.0
Rt = (2.0*p**3 - 9.0*p*q + 27.0*r)/54.0

r2q3 = Rt**2 - Qt**3

if (r2q3<0.0) then
  arg1 = Rt/Qt**1.5
  if (arg1>1.0)  arg1 = 1.0
  if (arg1<-1.0)  arg1 = -1.0
  th = acos(arg1)
  a1 =-2.0*Qt**0.5*cos(th/3.0) - p/3.0
  a2 =-2.0*Qt**0.5*cos((th + 2.0*pi)/3.0) - p/3.0
  a3 =-2.0*Qt**0.5*cos((th + 4.0*pi)/3.0) - p/3.0
  a = a1
  if (a2>a)  a = a2
  if (a3>a)  a = a3
else
  at =-Rt/abs(Rt)*(abs(Rt) + r2q3**0.5)**(1.0/3.0)
  if (abs(at)>1e-6) then
    bt = Qt/at
  else
    bt = 0.0
  endif
  a1 = at + bt - p/3.0
!        a2 =-0.5*(at + bt) - p/3.0 + i*3.0**0.5*(at - bt)/2.0
!        a3 =-0.5*(at + bt) - p/3.0 - i*3.0**0.5*(at - bt)/2.0
  a = a1
endif

!----------------------------------------------------------------------!
! Compute 'gs,pc,w' corresponding to 'a'.                              !
!----------------------------------------------------------------------!
gs = (g0 + g1*a*rh/pa)*kg
pc = pa - a*160.0/gs
w = vmax*pc/(pc + kc*(1 + po/ko))
acheck = (w*(1.0 - 0.5*po/(tau*pc)) - rd)*sc

!**********************************************************************!
end subroutine assvmax2
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                   assj2 :: productivity_methods                      !
!                   -----------------------------                      !
!                                                                      !
! ASSJ computes the assimilation rate by forming a cubic equation      !
! from the equations defining assimilation rate, stomatal conductance, !
! CO2 partial pressure (assimilation diffusion equaiton) and the       !
! carboxylation rate. The cubic equation is solved and the largest     !
! real root is assumed to be the assimilation rate. The variables here !
! correspond to the Doly paper, in the main program they are such that !
!        po,pa,j,g0,g1,rh,kg,tau,rd                                    !
!     = oi,ca,j(i),c1,c2,rh,kg,tau,rd                                  !
!                                                                      !
! SUBROUTINE assj2(a,w,pc,gs,po,pa,j,g0,g1,rh,kg,tau,rd,acheck)        !
!                                                                      !
!**********************************************************************!
subroutine assj2(a,w,pc,gs,po,pa,j,g0,g1,rh,kg,tau,rd,acheck)
!**********************************************************************!
real(dp) :: a,w,pc,gs,po,pa,j,g0,g1,rh,kg,tau,rd,acheck,r2q3,arg1, &
 ax,bx,cx,dx,ex,fx,gx,hx,p,q,r,qt,rt,th,a1,a2,a3,pi,at,bt,sc
!----------------------------------------------------------------------!

pi = atan(1.0)*4.0
sc = 1000000.0

!----------------------------------------------------------------------!
! 'a2,b2,c2,d2' are such that 'w=(a2a + b2)/(c2a + d2)'.               !
!----------------------------------------------------------------------!
ax = j*tau*pa*(kg*g1*rh - 160.0)
bx = j*tau*pa*kg*g0*pa
cx = 4.0*pa*kg*tau*g1*rh - 640.0*pa*tau + 4.0*po*kg*g1*rh
dx = 4.0*pa**2.0*kg*tau*g0 + 4.0*po*kg*g0*pa

!----------------------------------------------------------------------!
! 'ex,fx,gx,hx' are such that '0.5po/(taupc)-1=(exa+fx)/(gxa+hx)'.     *
!----------------------------------------------------------------------!
ex =-tau*pa*kg*g1*rh + 160.0*tau*pa + 0.5*po*kg*g1*rh
fx =-tau*pa**2.0*kg*g0 + 0.5*po*kg*g0*pa
gx = tau*pa*(kg*g1*rh - 160.0)
hx = tau*pa**2.0*kg*g0

!----------------------------------------------------------------------!
! Cubic coefficients 'p,q,r'.                                          !
!----------------------------------------------------------------------!
p = sc*rd + (sc*ax*ex + cx*hx + dx*gx)/(cx*gx)
q = (sc*rd*cx*hx + sc*rd*dx*gx + sc*ax*fx + sc*bx*ex + dx*hx)/(cx*gx)
r = (sc*rd*dx*hx + sc*bx*fx)/(cx*gx)

!----------------------------------------------------------------------!
! Sove the cubic equaiton to find 'a'.                                 !
!----------------------------------------------------------------------!
Qt = (p**2 - 3.0*q)/9.0
Rt = (2.0*p**3 - 9.0*p*q + 27.0*r)/54.0

r2q3 = Rt**2 - Qt**3

if (r2q3<0.0) then
  arg1 = Rt/Qt**1.5
  if (arg1>1.0)  arg1 = 1.0
  if (arg1<-1.0)  arg1 = -1.0
  th = acos(arg1)
  a1 =-2.0*Qt**0.5*cos(th/3.0) - p/3.0
  a2 =-2.0*Qt**0.5*cos((th + 2.0*pi)/3.0) - p/3.0
  a3 =-2.0*Qt**0.5*cos((th + 4.0*pi)/3.0) - p/3.0
  a = a1
  if (a2>a)  a = a2
  if (a3>a)  a = a3
else
  at =-Rt/abs(Rt)*(abs(Rt) + r2q3**0.5)**(1.0/3.0)
  if (abs(at)>1e-6) then
    bt = Qt/at
  else
    bt = 0.0
  endif
  a1 = at + bt - p/3.0
!        a2 =-0.5*(at + bt) - p/3.0 + i*3.0**0.5*(at - bt)/2.0
!        a3 =-0.5*(at + bt) - p/3.0 - i*3.0**0.5*(at - bt)/2.0
  a = a1
endif

!----------------------------------------------------------------------!
! Compute 'gs,pc,w' corresponding to 'a'.                              !
!----------------------------------------------------------------------!
gs = (g0 + g1*a*rh/pa)*kg
pc = pa - a*160.0/gs
w = j*pc/4.0/(pc + po/tau)
acheck = (w*(1.0 - 0.5*po/(tau*pc)) - rd)*sc

!**********************************************************************!
end subroutine assj2
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                   xamax :: productivity_methods                      !
!                   -----------------------------                      !
!                                                                      !
! This subroutine calculates the value of 'amax' at the maximum        !
! iradience as opposed to the average irradiance in the main program.  !
!                                                                      !
! SUBROUTINE xamax(a,xt,oi,vmx,rd,j,dresp,pc)                          !
!                                                                      !
!**********************************************************************!
subroutine xamax(a,xt,oi,vmx,rd,j,dresp,pc)
!**********************************************************************!
real(dp) :: a,xt,oi,vmx,rd,j,dresp,pc,t,tk,c1,c2,kc,ko,tau,tpav,tpwv,&
 tpaj,tpwj,sc
!       tppcj,tppcv,tpgsv,tpgsj,acheck
!----------------------------------------------------------------------!

t = 2.64 + 1.14*xt
tk = 273.0 + t
sc = 1000000.0

c1 = 142.4 - 4.8*t
if (c1>80.0)  c1 = 80.0
if (c1<8.0)  c1 = 8.0
c2 = 12.7 - 0.207*t
if (c2>10.0)  c2 = 10.0
if (c2<6.9)  c2 = 6.9
kc = exp(35.8 - 80.5/(0.00831*tk))
ko = exp(9.6 - 14.51/(0.00831*tk))*1000.0
tau = exp(-3.949 + 28.99/(0.00831*tk))

!      rht = rh
!      IF (kg*c2*rh.LT.170.0)  rht = 170.0/(kg*c2)
!      CALL ASSVMAX(tpav,tpwv,tppcv,tpgsv,oi,ca,vmx,c1,c2,
!     &rht,ko,kc,kg,tau,rd,acheck)
!      CALL ASSJ(tpaj,tpwj,tppcj,tpgsj,oi,ca,j,c1,c2,rht,
!     &kg,tau,rd,acheck)

tpwv = vmx*pc/(pc + kc*(1.0 + oi/ko))
tpav = (tpwv*(1.0 - 0.5*oi/(tau*pc)) - rd)*sc

tpwj = j*pc/(4.0*(pc + oi/tau))
tpaj = (tpwj*(1.0 - 0.5*oi/(tau*pc)) - rd)*sc

if (tpav<tpaj) then
  A = tpav
else
  A = tpaj
endif
if (A<-dresp)  A =-dresp

!**********************************************************************!
end subroutine xamax
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                   assc4 :: productivity_methods                      !
!                   -----------------------------                      !
! ASS4 computes the assimilation rate by forming a cubic equation      !
! from the equations defining assimilation rate, stomatal conductance, !
! CO2 partial pressure (assimilation diffusion equaiton) and the       !
! carboxylation rate. The cubic equation is solved and the largest     !
! real root is assumed to be the assimilation rate. This is then       !
! compared with assimilation limited to light dependant and vmax       !
! dependant assimilation and the minimum one is taken. The variables   !
! here correspond to the Doly paper, in the main program in the        !
! folowing mannor.                                                     !
!        po,pa,j,g0,g1,rh,kg,tau,rd                                    !
!     = oi,ca,j(i),c1,c2,rh,kg,tau,rd                                  !
!                                                                      !
! SUBROUTINE ASSC4(a,pc,gs,pa,vmax,rh,t,qg,p,up)                       !
!                                                                      !
!**********************************************************************!
subroutine ASSC4(a,pc,gs,pa,vmax,rh,t,qg,p,up)
!**********************************************************************!
real(dp) :: a,pc,gs,pa,vmax,rh,t,qg,up,p,vmq,absx,qrub,f,jqq,dv,amxt,&
 ac4t
!----------------------------------------------------------------------!

absx = 0.8
qrub = 0.11
f = 0.6
vmq = 2.0

dv = 0.6108*exp(17.269*t/(237.3 + t))*(1.0 - rh/100.0)/101.325
! kilo pascals to mol/mol = divide by 101.325

pc = pa*(0.5 - (1.6*dv*p/2500.0/pa)**0.5)

vmax = up*360.0*vmq**((t - 25.0)/10.0)/(360.0 + up)
vmax = vmax/((1 + exp(0.3*(10.0 - t)))*(0.8 + exp(0.14*(t - 36.0))))
vmax = vmax*1.0e-6*1.25

jqq = absx*qrub*f*qg*1.0e+6
vmq = 2.0

amxt = up*190.0/(360.0 + up)
ac4t = amxt*vmq**((t - 25.0)/10.0)/(1.0 + &
 exp(0.3*(10.0 - t)))/(0.8 + exp(0.14*(t - 36.0)))
ac4t = up/2.0*vmq**((t - 25.0)/10.0)/(1.0 +  &
 exp(0.3*(15.0 - t)))/(0.8 + exp(0.14*(t - 33.0)))


if (ac4t<jqq) then
  a = ac4t
else
  a = jqq
endif

gs = a*1.6*p/(pa - pc)*1.0e-3

!**********************************************************************!
end subroutine ASSC4
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                   assc42 :: productivity_methods                     !
!                   ------------------------------                     !
!                                                                      !
! ASS4 computes the assimilation rate by forming a cubic equation      !
! from the equations defining assimilation rate, stomatal conductance, !
! CO2 partial pressure (assimilation diffusion equaiton) and the       !
! carboxylation rate. The cubic equation is solved and the largest     !
! real root is assumed to be the assimilation rate. This is then       !
! compared with assimilation limited to light dependant and vmax       !
! dependant assimilation and the minimum one is taken. The variables   !
! here correspond to the Doly paper, in the main program in the        !
! folowing mannor.                                                     !
!        po,pa,j,g0,g1,rh,kg,tau,rd                                    !
!     = oi,ca,j(i),c1,c2,rh,kg,tau,rd                                  !
!                                                                      !
! SUBROUTINE assc42(a,w,pc,gs,po,pa,vmax,g0,g1,rh,kg,tau,rd,t,qg,      !
! acheck,up)                                                           !
!                                                                      !
!**********************************************************************!
subroutine assc42(a,w,pc,gs,po,pa,vmax,g0,g1,rh,kg,tau,rd,t,qg,acheck,&
 up)
!**********************************************************************!
real(dp) :: a,w,pc,gs,po,pa,vmax,g0,g1,rh,kg,tau,rd,t,qg,acheck,up, &
 p,q,r,qt,rt,th,a1,a2,a3,pi,at,bt,sc,kp,kmq,pg,l,vmq,absx, &
 qrub,f,a0,b0,b1,c1,d1,e1,f1,b2,c2,d2,b3,c3,d3,wq,wv, &
 a4,b4,c4,d4,f2,aq,av,a2t,r2q3,arg1
!----------------------------------------------------------------------!

pi = atan(1.0)*4.0
sc = 1000000.0

kp = 0.7
kmq = 2.0
pg = 101325.0
l = 0.000005
vmq = 2.0
absx = 0.8
qrub = 0.11
f = 0.6

vmax = up*360.0*vmq**((t - 25.0)/10.0)/(360.0 + up)
vmax = vmax/((1 + exp(0.3*(10.0 - t)))*(0.8 + exp(0.14*(t - 36.0))))
vmax = vmax*1.0e-6*1.25

!----------------------------------------------------------------------!
! 'a' and 'b' are such that 'jc=a pc  - b'.                            !
!----------------------------------------------------------------------!
a0 = kp*kmq**((t - 25.0)/10.0)/pg
b0 = l*kmq**((t - 25.0)/10.0)/pg

!----------------------------------------------------------------------!
! Reduce equations to a cubic in 'a'.                                  !
!----------------------------------------------------------------------!
a1 =-b0*sc
b1 = a0*sc
c1 =-0.5*po
d1 = tau
e1 =-tau*rd*sc
f1 = tau

a2 = g0*kg*(pa**2.0)
b2 = g0*kg*pa
c2 =-pa*160.0 + g1*rh*kg*pa
d2 = g1*rh*kg

a3 = a1*b2 + a2*b1
b3 = a1*d2 + c2*b1
c3 = c1*b2 + a2*d1
d3 = c1*d2 + c2*d1

f2 =-f1

a4 = f2*c2*d2
b4 = f2*(a2*d2 + b2*c2) + d3*b3 + e1*c2*d2
c4 = f2*a2*b2 + c3*b3 + d3*a3 + e1*(a2*d2 + c2*b2)
d4 = c3*a3 + e1*a2*b2

!----------------------------------------------------------------------!
! Cubic coefficients 'p,q,r'.                                          !
!----------------------------------------------------------------------!
p = b4/a4
q = c4/a4
r = d4/a4

a2t = a2
!----------------------------------------------------------------------!
! Sove the cubic equaiton to find 'a'.                                 !
!----------------------------------------------------------------------!
Qt = (p**2 - 3.0*q)/9.0
Rt = (2.0*p**3 - 9.0*p*q + 27.0*r)/54.0

r2q3 = Rt**2 - Qt**3

if (r2q3<0.0) then
  arg1 = Rt/Qt**1.5
  if (arg1>1.0)  arg1 = 1.0
  if (arg1<-1.0)  arg1 = -1.0
  th = acos(arg1)
  a1 =-2.0*Qt**0.5*cos(th/3.0) - p/3.0
  a2 =-2.0*Qt**0.5*cos((th + 2.0*pi)/3.0) - p/3.0
  a3 =-2.0*Qt**0.5*cos((th + 4.0*pi)/3.0) - p/3.0
  a = a1
  if (a2>a)  a = a2
  if (a3>a)  a = a3
else
  at =-Rt/abs(Rt)*(abs(Rt) + r2q3**0.5)**(1.0/3.0)
  if (abs(at)>1e-6) then
    bt = Qt/at
  else
    bt = 0.0
  endif
  a1 = at + bt - p/3.0
!        a2 =-0.5*(at + bt) - p/3.0 + i*3.0**0.5*(at - bt)/2.0
!        a3 =-0.5*(at + bt) - p/3.0 - i*3.0**0.5*(at - bt)/2.0
  a = a1
endif

gs = (g0 + g1*a*rh/pa)*kg
pc = pa - a*160.0/gs
w = a0*pc - b0
acheck = (w*(1.0 - 0.5*po/(tau*pc)) - rd)*sc

a2 = a2t

!----------------------------------------------------------------------!
! Compute light dependent assimilation rate.                           !
!----------------------------------------------------------------------!

wq = absx*qrub*f*qg*1.0
p = tau*c2
q = tau*a2 + wq*0.5*po*sc*d2 - tau*c2*(wq*sc - rd*sc)
r = wq*0.5*po*sc*b2 - tau*a2*(wq*sc - rd*sc)
aq = (-q + (abs(q)**2.0 - 4.0*p*r)**0.5)/(2.0*p)

!----------------------------------------------------------------------!
! Compute 'vmax' dependent assimilation rate.                          !
!----------------------------------------------------------------------!

wv = vmax*vmq**((t - 25.0)/10.0)/((1.0 + exp(0.3* &
 (10.0 - t)))*(0.8 + exp(0.14*(t - 36.0))))

q = tau*a2 + wv*0.5*po*sc*d2 - tau*c2*(wv*sc - rd*sc)
r = wv*0.5*po*sc*b2 - tau*a2*(wv*sc - rd*sc)
av = (-q + (abs(q)**2.0 - 4.0*p*r)**0.5)/(2.0*p)

!----------------------------------------------------------------------!
! Find limiting assimilation value.                                    !
!----------------------------------------------------------------------!
if (aq<a) then
  a = aq
  w = wq
  gs = (g0 + g1*a*rh/pa)*kg
  pc = pa - a*160.0/gs
endif

if (av<a) then
  a = av
  w = wv
  gs = (g0 + g1*a*rh/pa)*kg
  pc = pa - a*160.0/gs
endif

!**********************************************************************!
end subroutine assc42
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                   xvjmax :: productivity_methods                     !
!                   ------------------------------                     !
!                                                                      !
! This subroutine calculates the value of 'amax' at the maximum        !
! iradience as opposed to the average irradiance in the main program.  !
!                                                                      !
! SUBROUTINE xvjmax(vmx,jmx,j,xt,xq,sum,nup,oi,dresp)                  !
!                                                                      !
!**********************************************************************!
subroutine xvjmax(vmx,jmx,j,xt,xq,sum,nup,oi,dresp)
!**********************************************************************!
real(dp) :: vmx,jmx,j,xt,xq,sum,nup,oi,tk,shapex,maxx,conv,cc,ea,num,&
 denom,up,kc,ko,tau,can,upt,am,vm,jm,irc,q,t,aa,bb,dresp
!----------------------------------------------------------------------!

t = 2.64 + 1.14*xt
tk = 273.0 + t
q = 2.0*xq

nup = 1.5*nup

shapex = 40.8 + 0.01*t - t**2*0.002
maxx = 0.738 - 0.002*t
if (nup<=0.0) then
  write(*,*) 'Problem in XVJMAX nup=',nup
  stop
endif
conv = 97.4115 - 2.5039*log(nup)

cc = 36.0*exp(-t*0.14) + 20.0
ea = 81000.0*exp(-t*0.14) + 50000.0
num = exp(shapex - conv/(0.00831*tk))
denom = 1.0 + exp((maxx*tk - 205.9)/(0.00831*tk))
up = num/denom

kc = exp(35.8 - 80.5/(0.00831*tk))
ko = exp(9.6 - 14.51/(0.00831*tk))*1000.0
tau = exp(-3.949 + 28.99/(0.00831*tk))
aa = 24.5/(24.5 + kc + kc*oi/ko)
bb = 0.5/(24.5 + kc + kc*oi/ko)


can = 1.0/sum
upt = up*can
am = upt*190.0/(360.0 + upt)

dresp = exp(cc-(ea/(8.3144*tk)))*upt/50.0

vm =-1.0*(am + 0.82)/(-1.0*aa + bb*oi/tau)
jm = 29.1 + 1.64*vm
vm = vm/1000000.0
jm = jm/1000000.0

if (t>6.34365) then
  jmx = jm*(1.0 + 0.0409*(t-25.0)- 1.54e-3*(t - 25.0)**2 - &
 9.42e-5*(t - 25.0)**3)
else
  jmx = jm*0.312633
endif

if (t>9.51718) then
  vmx = vm*(1.0 + 0.0505*(t - 25.0)- 0.248e-3*(t - 25.0)**2 - &
 8.09e-5*(t-25.0)**3)
else
  vmx = vm*0.458928
endif

! Mark we have to discuss about this calculation for irc here...
irc = q*exp(-0.5)
j = 0.24*irc/(1.0 + (0.24**2)*(irc**2)/(jmx**2))**0.5

!**********************************************************************!
end subroutine xvjmax
!**********************************************************************!

end module productivity_methods
