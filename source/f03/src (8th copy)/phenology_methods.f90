!> @brief Methods which control the phenology, carbon allocation, and
!! the aging and updating of the active carbon stem, root, leaf
!! compartments. Along with suma which is also compartmentalised in a similar manner.
!! @details The active carbon pools of stem root and leaf, as well as suma
!! are compartmentalised into bins of an integral number of days long eg.
!! 'lai_comp_length' in the 'dims' module.
!!
!! @todo Basicly the compartments are behaving like a linked list but
!! using more memory. More memoty efficient to use linked lists.
!!
!! @section Compartments
!! Compartments are created as and when required. Each compartment
!! is a structure containing a real variable 'val' and an integer value
!! 'age'. The variable 'age' is set equal to the Julian day that the
!! compartment was created, and remains constant at that value.
!! Once created a compartment will accumulate any carbon for the following
!! number of days given by 'var_comp_length'.
module phenology_methods



use misc_values
use system_state
use pft_parameters
use site_parameters
use dims
use tuning_parameters
use misc_parameters
implicit none

contains

!**********************************************************************!
!                   phenology :: phenological_methods                  !
!                   ---------------------------------                  !
!                                                                      !
! Used for grasses and crops. Lai being controlled by the available    !
! npp storage.                                                         !
!                                                                      !
! SUBROUTINE phenology(yield,laiinc)                                   !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
subroutine phenology(yield,laiinc)
!**********************************************************************!
real(dp) :: yield,laiinc
integer :: co
!----------------------------------------------------------------------!

co = ssp%cohort

if (pft(co)%phen/=0) then
  if (pft(co)%phen==1) then
    call PHENOLOGY1(yield,laiinc)
  elseif (pft(co)%phen==2) then
    call PHENOLOGY2(yield,laiinc)
  else
    write(*,*) 'No phenology defined for ',pft(co)%phen,pft(co)%tag
  stop
  endif
endif

end subroutine phenology





!**********************************************************************!
!                   phenology1 :: phenological_methods                 !
!                   ----------------------------------                 !
!                                                                      !
! Used for trees. Lai being controlled by a weighted memory of npp     !
! storage over the previous 5 years.                                   !
!                                                                      !
! SUBROUTINE phenology1(yield,laiinc)                                    !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
!**********************************************************************!
subroutine phenology1(yield,laiinc)
!**********************************************************************!
real(dp), parameter :: maxlai = 10.9
real(dp) :: rlai,lairat,laiinc,wtwp,wtfc,bb0,bbmax,bblim,sslim,bbsum, &
 yield,smtrig,tsuma,stemfr,maint,ftagh
integer :: lai,leafls,mnth,day,ij,i,bb,gs,bbgs,sssum
integer :: bb2bbmin,bb2bbmax,bbm,ssm,sss,ss,dsbb
integer :: ftdth,harvest,chill,dschill,co
!----------------------------------------------------------------------!

co = ssp%cohort

wtwp = ssp%wilt
day  = ssp%day
mnth = ssp%mnth

wtfc     = ssp%field
ftagh    = pft(co)%crop
ftdth    = pft(co)%d2h
stemfr   = ssv(co)%stemfr
bb       = ssv(co)%bb
ss       = ssv(co)%ss
bbgs     = ssv(co)%bbgs
dsbb     = ssv(co)%dsbb
chill    = ssv(co)%chill
dschill  = ssv(co)%dschill
leafls   = pft(co)%lls
bbm      = pft(co)%bbmem
bb0      = pft(co)%bb0
bbmax    = pft(co)%bbmax
bblim    = pft(co)%bblim
ssm      = pft(co)%senm
sss      = pft(co)%sens
sslim    = pft(co)%senlim
lairat   = pft(co)%lrat

bb2bbmin = 285
bb2bbmax = 390
gs = 60

rlai = ssv(co)%lai%tot

yield = 0.0
harvest = 0

ss = 0

!----------------------------------------------------------------------!
! Check for chilling.                                                  !
!----------------------------------------------------------------------!
if (chill==0) then
  bbsum = 0.0
  do i=1,20
    if (ssp%tmem(i)<-5.0) &
  bbsum = bbsum + max(-10.0,ssp%tmem(i)+5.0)
  enddo
  if (bbsum<-100) then
    chill = 1
    dschill = 1
  endif
endif
if (chill==1) then
  dschill = dschill + 1
endif
if (dschill>260) then
  chill = 0
  dschill = 0
endif

!----------------------------------------------------------------------!
! Bubburst, if no budburst set, and sufficient soil moisture, then     *
! check for  budburst.                                                 !
!----------------------------------------------------------------------!
if (((bb==0).and.(msv%mv_soil2g>wtwp+0.25*(wtfc-wtwp))).or. &
 ((dsbb>bb2bbmax).and.(msv%mv_soil2g>wtwp+0.1*(wtfc-wtwp)))) then

  smtrig = 0.0
  do i=1,30
    smtrig = smtrig + ssv(co)%sm_trig(i)
  enddo

  if ((smtrig>30.0).or.(dsbb>bb2bbmax)) then

!----------------------------------------------------------------------!
! Check for budburst using degree days.                                !
!----------------------------------------------------------------------!
    bbsum = 0.0
    do i=1,bbm
      if (ssp%tmem(i)>bb0)  bbsum = bbsum + min(bbmax,ssp%tmem(i)-bb0)
    enddo

    if ((real(bbsum)>=real(bblim)*exp(-0.01*real(dschill))) &
 .or.(dsbb>bb2bbmax)) then
!----------------------------------------------------------------------!
! Adjust proportion of gpp going into stem production based on suma.   !
! This is essentially the LAI control.                                 !
!----------------------------------------------------------------------!
!      print*,stemfr
      tsuma = ssv(co)%suma%tot
!      print*,'tsuma ',tsuma,stemfr,ssv(co)%nppstore(1)
      maint = max(1.0,(real(leafls)/360.0)*1.0)
      tsuma = tsuma - msv%mv_leafmol*1.25/maint*tgp%p_opt
      if (tsuma>msv%mv_leafmol*1.25/maint*tgp%p_opt) tsuma = msv%mv_leafmol*1.25/maint*tgp%p_opt
      if (tsuma<-msv%mv_leafmol*1.25/maint*tgp%p_opt) tsuma = -msv%mv_leafmol*1.25/maint*tgp%p_opt
      stemfr = stemfr + tsuma*tgp%p_laimem*12.0

      if (stemfr<120.0) stemfr = 120.0

!----------------------------------------------------------------------!
! Budburst occurance.                                                  !
!----------------------------------------------------------------------!
      bb = (mnth-1)*30 + day
      bbgs = 0
      dsbb = 0

      if (stemfr<0.8*ssv(co)%nppstore(1)) then
        ssv(co)%nppstore(3) = ssv(co)%nppstore(1) - stemfr
      else
        if (stemfr<0.75*ssv(co)%nppstore(1)) then
          ssv(co)%nppstore(3) = ssv(co)%nppstore(1) - stemfr
        else
          ssv(co)%nppstore(3) = ssv(co)%nppstore(1)*0.25
        endif
        stemfr = stemfr*0.8
      endif
      laiinc = (ssv(co)%nppstore(1) - 0.0*ssv(co)%nppstore(3))/msv%mv_leafmol/1.25/12.0
      ssv(co)%nppstore(2) = ssv(co)%nppstore(1)
!      print*,'laiinc nppstore ',laiinc,ssv(co)%nppstore(1)
    endif
  endif
endif

!----------------------------------------------------------------------!
! Compute length of growing season, and set to zero when equal to      !
! growing season.                                                      !
!----------------------------------------------------------------------!
if (bb>0)  bbgs = bbgs + 1
if (bbgs-gs>bb2bbmin) then
  bb = 0
  bbgs = 0
endif

if (dsbb < 500) dsbb = dsbb + 1

!----------------------------------------------------------------------!
! Set LAI increase.                                                    !
!----------------------------------------------------------------------!
if ((bb>0).and.(bbgs<gs).and.(ssv(co)%nppstore(1)>0.0)) then
  laiinc = lairat*(ssv(co)%nppstore(1) - 0.0*ssv(co)%nppstore(3))/msv%mv_leafmol/1.0/12.0*2.0
  if (rlai+laiinc>maxlai)  laiinc = maxlai - rlai
  if (rlai+laiinc>11.5)  laiinc = 11.5 - rlai
  if ((rlai>0).and.(ssv(co)%nppstore(1)<0.0)) laiinc = 0.0
else
  laiinc = 0.0
endif

!----------------------------------------------------------------------!
! Senescence, if rlai is greater than zero, compute senescence.        !
!----------------------------------------------------------------------!
if (rlai>1.0e-6) then
  if (abs(ftagh)<1.0e-6) then
    if (msv%mv_soil2g<wtwp*0.0) then
!----------------------------------------------------------------------!
! Senescence event due to soil moisture.                               !
!----------------------------------------------------------------------!
      laiinc = -rlai
      ss = day + (mnth - 1)*30
    elseif (bbgs>100) then
!----------------------------------------------------------------------!
! Check for senescence, senescence occurs there are 'sss' days colder  *
! than 'sslim' out of the last 'ssm' days.                             !
!----------------------------------------------------------------------!
      sssum = 0
      do i=1,ssm
        if (ssp%tmem(i)<sslim)  sssum = sssum + 1
      enddo
      if (sssum>=sss) then
!----------------------------------------------------------------------!
! Senescence event due to temperature.                                 !
!----------------------------------------------------------------------!
        laiinc =-rlai
        ss = day + (mnth - 1)*30
      endif
    endif 
  else
!----------------------------------------------------------------------!
! Crop senescence
!----------------------------------------------------------------------!
    if ((bbgs>leafls).or.(bbgs>ftdth)) then
      harvest = 1
      laiinc =-rlai
      ss = day + (mnth - 1)*30
    endif
  endif
endif

!----------------------------------------------------------------------!
ssv(co)%stemfr      = stemfr
ssv(co)%bb          = bb
ssv(co)%ss          = ss
ssv(co)%bbgs        = bbgs
ssv(co)%dsbb        = dsbb
ssv(co)%chill       = chill
ssv(co)%dschill     = dschill

end subroutine phenology1





!**********************************************************************!
!                   phenology2 :: phenological_methods                 !
!                   ----------------------------------                 !
!                                                                      !
! Used for trees. Lai being controlled by a weighted memory of npp     !
! storage over the previous 5 years.                                   !
!                                                                      !
! SUBROUTINE phenology2(yield,laiinc)                                  !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
!**********************************************************************!
subroutine phenology2(yield,laiinc)
!**********************************************************************!
real(dp), parameter :: maxlai = 10.9
real(dp) :: rlai,lairat,wtwp,laiinc,wtfc,bb0,bbmax,bblim,sslim,bbsum,yield,smtrig, &
 tsuma,stemfr,maint
integer :: lai,leafls,ij,i,bb,gs,bbgs,sssum,bb2bbmin,bb2bbmax,bbm,ssm, &
 sss,ss,dsbb,chill,dschill,co,day,mnth
!----------------------------------------------------------------------!

co = ssp%cohort

wtwp = ssp%wilt
day = ssp%day
mnth = ssp%mnth

wtfc     = ssp%field
stemfr   = ssv(co)%stemfr
bb       = ssv(co)%bb
ss       = ssv(co)%ss
bbgs     = ssv(co)%bbgs
dsbb     = ssv(co)%dsbb
chill    = ssv(co)%chill
dschill  = ssv(co)%dschill
leafls   = pft(co)%lls
bbm      = pft(co)%bbmem
bb0      = pft(co)%bb0
bbmax    = pft(co)%bbmax
bblim    = pft(co)%bblim
ssm      = pft(co)%senm
sss      = pft(co)%sens
sslim    = pft(co)%senlim
lairat   = pft(co)%lrat

bb2bbmin = 315
bb2bbmax = 375
gs = 30

rlai = ssv(co)%lai%tot

yield = 0.0

ss = 0

!----------------------------------------------------------------------!
! Check for chilling.                                                  !
!----------------------------------------------------------------------!
if (chill==0) then
  bbsum = 0.0
  do i=1,20
    if (ssp%tmem(i)<-5.0) &
  bbsum = bbsum + max(-10.0,ssp%tmem(i)+5.0)
  enddo
  if (bbsum<-100) then
    chill = 1
    dschill = 1
  endif
endif
if (chill==1) then
  dschill = dschill + 1
endif
if (dschill>260) then
  chill = 0
  dschill = 0
endif

!if (bb > 0) stop

!----------------------------------------------------------------------!
! Bubburst, if no budburst set, and sufficient soil moisture, then     *
! check for  budburst.                                                 !
!----------------------------------------------------------------------!
if (((bb==0).and.(msv%mv_soil2g>wtwp+0.5*(wtfc-wtwp))).or. &
 ((dsbb>bb2bbmax).and.(msv%mv_soil2g>wtwp+0.1*(wtfc-wtwp)))) then

  smtrig = 0.0
  do i=1,30
    smtrig = smtrig + ssv(co)%sm_trig(i)
  enddo

  if ((smtrig>30.0).or.(dsbb>bb2bbmax)) then

!----------------------------------------------------------------------!
! Check for budburst using degree days.                                !
!----------------------------------------------------------------------!
    bbsum = 0.0
    do i=1,bbm
      if (ssp%tmem(i)>bb0)  bbsum = bbsum + min(bbmax,ssp%tmem(i)-bb0)
    enddo

    if ((real(bbsum)>=real(bblim)*exp(-0.01*real(dschill))) &
 .or.(dsbb>bb2bbmax)) then
!----------------------------------------------------------------------!
! Adjust proportion of gpp going into stem production based on suma.   *
! This is essentially the LAI control.                                 !
!----------------------------------------------------------------------!
      tsuma = ssv(co)%suma%tot
!      print*,'tsuma ',tsuma,stemfr,ssv(co)%nppstore(1)
      maint = max(1.0,(real(leafls)/360.0)*1.0)
      tsuma = tsuma - msv%mv_leafmol*1.25/maint*tgp%p_opt
      if (tsuma>msv%mv_leafmol*1.25/maint*tgp%p_opt) tsuma = msv%mv_leafmol*1.25/maint*tgp%p_opt
      if (tsuma<-msv%mv_leafmol*1.25/maint*tgp%p_opt) tsuma = -msv%mv_leafmol*1.25/maint*tgp%p_opt
      stemfr = stemfr + tsuma*tgp%p_laimem*12.0

      if (stemfr<120.0) stemfr = 120.0

!----------------------------------------------------------------------!
! Budburst occurance.                                                 !
!----------------------------------------------------------------------!
!      print*,mnth,day,' budburst'
      bb = (mnth-1)*30 + day
      bbgs = 0
      dsbb = 0

      if (stemfr<0.75*ssv(co)%nppstore(1)) then
        ssv(co)%nppstore(3) = ssv(co)%nppstore(1) - stemfr
      else
        if (stemfr<0.75*ssv(co)%nppstore(1)) then
          ssv(co)%nppstore(3) = ssv(co)%nppstore(1) - stemfr
        else
          ssv(co)%nppstore(3) = ssv(co)%nppstore(1)*0.25
        endif
        stemfr = stemfr*0.95
      endif
      laiinc = (ssv(co)%nppstore(1) - ssv(co)%nppstore(3))/msv%mv_leafmol/1.25/12.0
      ssv(co)%nppstore(2) = ssv(co)%nppstore(1)
!      print*,'laiinc ',laiinc,tsuma,stemfr
    endif
  endif
endif

!----------------------------------------------------------------------!
! Compute length of growing season, and set to zero when equal to      *
! growing season.                                                      !
!----------------------------------------------------------------------!
if (bb>0)  bbgs = bbgs + 1
if (bbgs-gs>bb2bbmin) then
  bb = 0
  bbgs = 0
endif

if (dsbb < 500) dsbb = dsbb + 1

!----------------------------------------------------------------------!
! Set LAI increase.                                                    !
!----------------------------------------------------------------------!
if ((bb>0).and.(bbgs<gs).and.(ssv(co)%nppstore(1)>1.0)) then
  laiinc = lairat*(ssv(co)%nppstore(2) - ssv(co)%nppstore(3))/msv%mv_leafmol/1.25/12.0
  if (rlai+laiinc>maxlai)  laiinc = maxlai - rlai
  if (rlai+laiinc>11.5)  laiinc = 11.5 - rlai
  if ((rlai>0).and.(ssv(co)%nppstore(1)<60.0)) laiinc = 0.0
else
  laiinc = 0.0
endif

!----------------------------------------------------------------------!
! Senescence, if rlai is greater than zero, compute senescence.        !
!----------------------------------------------------------------------!
if (rlai>1.0e-6) then
  if (msv%mv_soil2g<wtwp*0.0) then
!----------------------------------------------------------------------!
! Senescence event due to soil moisture.                               !
!----------------------------------------------------------------------!
    laiinc = -rlai
    ss = day + (mnth - 1)*30
  elseif (bbgs>100) then
!----------------------------------------------------------------------!
! Check for senescence, senescence occurs there are 'sss' days colder  *
! than 'sslim' out of the last 'ssm' days.                             !
!----------------------------------------------------------------------!
    sssum = 0
    do i=1,ssm
      if (ssp%tmem(i)<sslim)  sssum = sssum + 1
    enddo
    if (sssum>=sss) then
!----------------------------------------------------------------------!
! Senescence event due to temperature.                                 !
!----------------------------------------------------------------------!
      laiinc =-rlai
      ss = day + (mnth - 1)*30
    endif
  endif
endif
!----------------------------------------------------------------------!

ssv(co)%stemfr      = stemfr
ssv(co)%bb          = bb
ssv(co)%ss          = ss
ssv(co)%bbgs        = bbgs
ssv(co)%dsbb        = dsbb
ssv(co)%chill       = chill
ssv(co)%dschill     = dschill

end subroutine phenology2






!**********************************************************************!
!                                                                      !
!                    allocation :: phenological_methods                !
!                    ----------------------------------                !
! SUBROUTINE allocation(laiinc,daygpp,resp_l,lmor_sc,resp,leaflit,&    !
! stemnpp,rootnpp,resp_s,resp_r,resp_m,check_closure)                  !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Allocate GPP between storage, stems, roots and leaves.
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine allocation(laiinc,daygpp,resp_l,lmor_sc,resp,leaflit,&
 stemnpp,rootnpp,resp_s,resp_r,resp_m,check_closure)
!**********************************************************************!
real(dp) :: laiinc,leaflit,resp,ans,yy,stemnpp,total_carbon,old_total_carbon, &
 root_fixed,stem_fixed,lmor_sc(3600),daynpp,rootnpp,resp_s,resp_r,resp_l, &
 daygpp,resp_m
integer :: i,co,k
logical :: check_closure
real(dp) :: sumrr,sumsr,sumlr,summr
save :: sumrr,sumsr,sumlr,summr
!----------------------------------------------------------------------!
!print*,'------------------'
!print*,'daygpp ',daygpp, ssp%cohort
!print*,'resp_l ',resp_l
!print*,'laiinc ',laiinc
!print*,'nppstore ',ssv(co)%nppstore(1:3)

co = ssp%cohort

daynpp = daygpp - resp_l

if (pft(co)%phen/=0) then

ssv(co)%nppstore(1) = ssv(co)%nppstore(1) + daygpp - resp_l

resp = resp_l

!----------------------------------------------------------------------!
! Check carbon closure.
!----------------------------------------------------------------------!
if (check_closure) then
  old_total_carbon = ssv(co)%lai%tot*12.0/pft(co)%sla/18.0 + ssv(co)%nppstore(1) + &
 ssv(co)%stem%tot + ssv(co)%root%tot + resp + ssv(co)%bio(1) + ssv(co)%bio(2)
endif

!----------------------------------------------------------------------!
! Pay for new leaves.
!----------------------------------------------------------------------!
if (laiinc>0.0) then
  ssv(co)%nppstore(1) = ssv(co)%nppstore(1) - laiinc*msv%mv_leafmol*1.25*12.0
  ssv(co)%nppstore(2) = ssv(co)%nppstore(2) - laiinc*msv%mv_leafmol*1.25*12.0
  resp_m = 0.25*laiinc*msv%mv_leafmol*12.0
!  print*,ssv(co)%lai,laiinc,resp_m
  resp = resp + resp_m
else
  resp_m = 0.0
endif

!----------------------------------------------------------------------!
! Age leaves by one day, kill any which have died of old age,
! adjust by laiinc (+ or -), and then sum to get 'rlai'.
! 'leafls' is an integer variable of leaf lifespan in days.
!----------------------------------------------------------------------!
call LAI_ADD(laiinc,leaflit)

!----------------------------------------------------------------------!
! leaf death not through age mortality 
!----------------------------------------------------------------------!
call LAI_DIST(lmor_sc,leaflit)

!----------------------------------------------------------------------!
! Check carbon closure.
!----------------------------------------------------------------------!
if (check_closure) then
  ans = ssv(co)%lai%tot*12.0/pft(co)%sla/18.0 + ssv(co)%nppstore(1) + &
 ssv(co)%stem%tot + ssv(co)%root%tot + resp + ssv(co)%bio(1) + ssv(co)%bio(2)
  if (abs(ans+leaflit-old_total_carbon) > 1.0e-3) then
    write(*,*) 'Breach of carbon after leaves:',ans+leaflit-old_total_carbon,' g/m^2.'
    write(*,*) ans,leaflit,old_total_carbon
    write(*,*) 'lai ',ssv(co)%lai
    write(*,*) 'sla ',pft(co)%sla
    write(*,*) 'nppstore ',ssv(co)%nppstore(1)
    write(*,*) '%stem ',ssv(co)%stem
    write(*,*) '%root ',ssv(co)%root
    write(*,*) 'resp ',resp
    write(*,*) 'biomass ',ssv(co)%bio(1)
    write(*,*) 'biomass ',ssv(co)%bio(2)
    stop
  endif
endif
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Pay for days roots if veg exists
!----------------------------------------------------------------------!
if ((ssv(co)%nppstore(1)>0.0).and.(daynpp>0.0)) then
  yy = ssv(co)%nppstore(1)*tgp%p_rootfr
else
  yy = 0.0
endif

ssv(co)%nppstore(1) = ssv(co)%nppstore(1) - yy
rootnpp = yy

!----------------------------------------------------------------------!
! Age roots and add todays root npp.
!----------------------------------------------------------------------!
call ROOT_ADD(yy,root_fixed)
ssv(co)%bio(2) = ssv(co)%bio(2) + root_fixed

!----------------------------------------------------------------------!
! Root resperation
!----------------------------------------------------------------------!
call ROOT_DIST(msv%mv_respref,resp_r)
resp = resp + resp_r
rootnpp = rootnpp - resp_r

!----------------------------------------------------------------------!
! Check carbon closure.
!----------------------------------------------------------------------!
if (check_closure) then
  ans = ssv(co)%lai%tot*12.0/pft(co)%sla/18.0 + ssv(co)%nppstore(1) + &
 ssv(co)%stem%tot + ssv(co)%root%tot + resp + ssv(co)%bio(1) + ssv(co)%bio(2)
  if (abs(ans+leaflit-old_total_carbon) > 1.0e-3) then
    write(*,*) 'Breach of carbon after roots:',ans+leaflit-old_total_carbon,' g/m^2.'
    stop
  endif
endif
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Stem resperation, and NPP                                            !
!----------------------------------------------------------------------!
if ((ssv(co)%nppstore(1)>0.0).and.(daynpp>0.0)) then
  yy = ssv(co)%nppstore(1)*tgp%p_stemfr
  ssv(co)%nppstore(1) = ssv(co)%nppstore(1) - yy
else
  yy = 0.0
endif
stemnpp = yy

!----------------------------------------------------------------------!
! Age stems and add todays stem npp.
!----------------------------------------------------------------------!
call STEM_ADD(yy,stem_fixed)
ssv(co)%bio(1) = ssv(co)%bio(1) + stem_fixed

!----------------------------------------------------------------------!
! Stem respiration.
!----------------------------------------------------------------------!
call STEM_DIST(msv%mv_respref,resp_s)
resp = resp + resp_s
stemnpp = stemnpp - resp_s

!----------------------------------------------------------------------!
! Check carbon closure.
!----------------------------------------------------------------------!
if (check_closure) then
  ans = ssv(co)%lai%tot*12.0/pft(co)%sla/18.0 + ssv(co)%nppstore(1) + &
 ssv(co)%stem%tot + ssv(co)%root%tot + resp + ssv(co)%bio(1) + ssv(co)%bio(2)
  if (abs(ans+leaflit-old_total_carbon) > 1.0e-3) then
    write(*,*) 'Breach of carbon after stems:',ans+leaflit-old_total_carbon,' g/m^2.'
    stop
  endif
endif
!----------------------------------------------------------------------!

daynpp = daynpp - resp
daynpp = daygpp - resp
ssv(co)%npp = daynpp

!----------------------------------------------------------------------!
! Check carbon closure.
!----------------------------------------------------------------------!
if (check_closure) then
  total_carbon = ssv(co)%lai%tot*12.0/pft(co)%sla/18.0 + ssv(co)%nppstore(1) + &
 ssv(co)%stem%tot + ssv(co)%root%tot + resp + ssv(co)%bio(1) + ssv(co)%bio(2)
  if (abs(total_carbon+leaflit-old_total_carbon) > 1.0e-3) then
    write(*,*) 'Breach of carbon in phenology:',ans+leaflit-old_total_carbon,' g/m^2.'
    stop
  endif
endif

endif

if (ssp%cohort == -1) then
  print*,'gpp & resp',daygpp,resp_l
  print*,'maint & lai inc ',resp_m,laiinc
  print*,ssp%mnth,ssp%day,resp_s,resp_r
!  print'(2i4,6f10.3)',ssp%mnth,ssp%day,daygpp,daynpp,resp_l,resp_s,resp_r,resp_m
endif

if (ssp%cohort == 1) then
      if ((ssp%day == 1).and.(ssp%mnth == 1)) then
        sumsr=0.0
        sumrr=0.0
        sumlr=0.0
        summr=0.0
      endif
      sumsr = sumsr + resp_s
      sumrr = sumrr + resp_r
      sumlr = sumlr + resp_l
      summr = summr + resp_m
      if ((ssp%day == -30).and.(ssp%mnth == 12)) print'(4f10.3)',sumsr,sumrr,sumlr,summr
endif

end subroutine allocation





!**********************************************************************!
!                                                                      !
!                      suma_add :: phenological_methods                !
!                      --------------------------------                !
!                                                                      !
!  SUBROUTINE suma_add(laiinc)                                         !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Update suma by adding daily increase and removing via mortality.
!! @details suma mortality is hardwired at 360 days.
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine suma_add(laiinc)
!**********************************************************************!
real(dp) :: laiinc
integer :: i,co,k
!----------------------------------------------------------------------!

co = ssp%cohort

!----------------------------------------------------------------------!
! Leaf death through old age.
!----------------------------------------------------------------------!
if (ssv(co)%suma%no > 0) then
  if (ssp%jday-int(ssv(co)%suma%c(1)%age) > 360) then
    k = suma_comp_length-(ssp%jday-int(ssv(co)%suma%c(1)%age)-360)+1
    if (k == 1) then
      do i=1,ssv(co)%suma%no-1
        ssv(co)%suma%c(i)%val = ssv(co)%suma%c(i+1)%val
        ssv(co)%suma%c(i)%age = ssv(co)%suma%c(i+1)%age
      enddo
      ssv(co)%suma%c(ssv(co)%suma%no)%val = 0.0
      ssv(co)%suma%c(ssv(co)%suma%no)%age = 0.0
      ssv(co)%suma%no = ssv(co)%suma%no - 1
    else
      ssv(co)%suma%c(1)%val = ssv(co)%suma%c(1)%val*(real(k-1)/real(k))
    endif
  endif
endif

!----------------------------------------------------------------------!
! Add on lai increase, create new compartment if necessary.
!----------------------------------------------------------------------!
if (laiinc>0.0) then
  if (ssv(co)%suma%no == 0) then
    ssv(co)%suma%no = 1
    ssv(co)%suma%c(ssv(co)%suma%no)%val = 0.0
    ssv(co)%suma%c(ssv(co)%suma%no)%age = ssp%jday
  elseif (ssp%jday-int(ssv(co)%suma%c(ssv(co)%suma%no)%age) >= suma_comp_length) then
    ssv(co)%suma%no = ssv(co)%suma%no + 1
    ssv(co)%suma%c(ssv(co)%suma%no)%val = 0.0
    ssv(co)%suma%c(ssv(co)%suma%no)%age = ssp%jday
  endif
  ssv(co)%suma%c(ssv(co)%suma%no)%val =  ssv(co)%suma%c(ssv(co)%suma%no)%val + laiinc
endif

ssv(co)%suma%tot = 0.0
do i=1,ssv(co)%suma%no
  ssv(co)%suma%tot = ssv(co)%suma%tot + ssv(co)%suma%c(i)%val
enddo

end subroutine suma_add
   




!**********************************************************************!
!                                                                      !
!                   stem_add  :: phenological_methods                  !
!                   ---------------------------------                  !
!                                                                      !
!  SUBROUTINE stem_add(laiinc,leaflit)                                 !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Update stems by adding daily increase and removing via mortality.
!! @details Stem mortality is given in the pft parameters (pft%sls).
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine stem_add(stem_inc,stem_fixed)
!**********************************************************************!
real(dp) :: stem_inc,stem_fixed
integer :: i,co,k
!----------------------------------------------------------------------!

co = ssp%cohort
stem_fixed = 0.0

!----------------------------------------------------------------------!
! Stem death through old age.
!----------------------------------------------------------------------!
if (ssv(co)%stem%no > 0) then
  if (ssp%jday-int(ssv(co)%stem%c(1)%age) > pft(co)%sls) then
    k = stem_comp_length-(ssp%jday-int(ssv(co)%stem%c(1)%age)-pft(co)%sls)+1
    if (k == 1) then
      stem_fixed = ssv(co)%stem%c(1)%val
      do i=1,ssv(co)%stem%no-1
        ssv(co)%stem%c(i)%val = ssv(co)%stem%c(i+1)%val
        ssv(co)%stem%c(i)%age = ssv(co)%stem%c(i+1)%age
      enddo
      ssv(co)%stem%c(ssv(co)%stem%no)%val = 0.0
      ssv(co)%stem%c(ssv(co)%stem%no)%age = 0.0
      ssv(co)%stem%no = ssv(co)%stem%no - 1
    else
      stem_fixed = ssv(co)%stem%c(1)%val/real(k)
      ssv(co)%stem%c(1)%val = ssv(co)%stem%c(1)%val*(real(k-1)/real(k))
    endif
  endif
endif

!----------------------------------------------------------------------!
! Add on stem increase, create new compartment if necessary.
!----------------------------------------------------------------------!
if (stem_inc>0.0) then
  if (ssv(co)%stem%no == 0) then
    ssv(co)%stem%no = 1
    ssv(co)%stem%c(ssv(co)%stem%no)%val = 0.0
    ssv(co)%stem%c(ssv(co)%stem%no)%age = ssp%jday
  elseif (ssp%jday-int(ssv(co)%stem%c(ssv(co)%stem%no)%age) >= stem_comp_length) then
    ssv(co)%stem%no = ssv(co)%stem%no + 1
    ssv(co)%stem%c(ssv(co)%stem%no)%val = 0.0
    ssv(co)%stem%c(ssv(co)%stem%no)%age = ssp%jday
  endif
  ssv(co)%stem%c(ssv(co)%stem%no)%val =  ssv(co)%stem%c(ssv(co)%stem%no)%val + stem_inc
endif

end subroutine stem_add





!**********************************************************************!
!                                                                      !
!                    root_add :: phenology_methods                     !
!                    -----------------------------                     !
!                                                                      !
!  SUBROUTINE root_add(laiinc,leaflit)                                 !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Update roots by adding daily increase and removing via mortality.
!! @details Root mortality is given in the pft parameters (pft%rls).
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine root_add(laiinc,leaflit)
!**********************************************************************!
real(dp) :: laiinc,leaflit
integer i,co,k
!----------------------------------------------------------------------!

co = ssp%cohort
leaflit = 0.0

!----------------------------------------------------------------------!
! Leaf death through old age.
!----------------------------------------------------------------------!
if (ssv(co)%root%no > 0) then
  if (ssp%jday-int(ssv(co)%root%c(1)%age) > pft(co)%rls) then
    k = root_comp_length-(ssp%jday-int(ssv(co)%root%c(1)%age)-pft(co)%rls)+1
    if (k == 1) then
      leaflit = ssv(co)%root%c(1)%val
      do i=1,ssv(co)%root%no-1
        ssv(co)%root%c(i)%val = ssv(co)%root%c(i+1)%val
        ssv(co)%root%c(i)%age = ssv(co)%root%c(i+1)%age
      enddo
      ssv(co)%root%c(ssv(co)%root%no)%val = 0.0
      ssv(co)%root%c(ssv(co)%root%no)%age = 0.0
      ssv(co)%root%no = ssv(co)%root%no - 1
    else
      leaflit = ssv(co)%root%c(1)%val/real(k)
      ssv(co)%root%c(1)%val = ssv(co)%root%c(1)%val*&
 (real(k-1)/real(k))
    endif
  endif
endif

!----------------------------------------------------------------------!
! Add on lai increase, create new compartment if necessary.
!----------------------------------------------------------------------!
if (laiinc>0.0) then
  if (ssv(co)%root%no == 0) then
    ssv(co)%root%no = 1
    ssv(co)%root%c(ssv(co)%root%no)%val = 0.0
    ssv(co)%root%c(ssv(co)%root%no)%age = ssp%jday
  elseif (ssp%jday-int(ssv(co)%root%c(ssv(co)%root%no)%age) &
 >= root_comp_length) then
    ssv(co)%root%no = ssv(co)%root%no + 1
    ssv(co)%root%c(ssv(co)%root%no)%val = 0.0
    ssv(co)%root%c(ssv(co)%root%no)%age = ssp%jday
  endif
  ssv(co)%root%c(ssv(co)%root%no)%val =  &
 ssv(co)%root%c(ssv(co)%root%no)%val + laiinc
endif

end subroutine root_add





!**********************************************************************!
!                                                                      !
!                    lai_add :: phenological_methods                   !
!                    -------------------------------                   !
!                                                                      !
! SUBROUTINE lai_add(laiinc,leaflit)                                   !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine lai_add(laiinc,leaflit)
!**********************************************************************!
real(dp) :: laiinc,leaflit
integer :: i,co,k
!----------------------------------------------------------------------!

co = ssp%cohort

!----------------------------------------------------------------------!
! Leaf death through old age.
!----------------------------------------------------------------------!
if (ssv(co)%lai%no > 0) then
  if (ssp%jday-int(ssv(co)%lai%c(1)%age) > pft(co)%lls) then
    k = lai_comp_length-(ssp%jday-int(ssv(co)%lai%c(1)%age)-pft(co)%lls)+1
    if (k == 1) then
      leaflit = ssv(co)%lai%c(1)%val
      do i=1,ssv(co)%lai%no-1
        ssv(co)%lai%c(i)%val = ssv(co)%lai%c(i+1)%val
        ssv(co)%lai%c(i)%age = ssv(co)%lai%c(i+1)%age
      enddo
      ssv(co)%lai%c(ssv(co)%lai%no)%val = 0.0
      ssv(co)%lai%c(ssv(co)%lai%no)%age = 0.0
      ssv(co)%lai%no = ssv(co)%lai%no - 1
    else
      leaflit = ssv(co)%lai%c(1)%val/real(k)
      ssv(co)%lai%c(1)%val = ssv(co)%lai%c(1)%val*(real(k-1)/real(k))
    endif
  endif
endif

!----------------------------------------------------------------------!
! Add on lai increase, create new compartment if necessary.
!----------------------------------------------------------------------!
if (laiinc>0.0) then
  if (ssv(co)%lai%no == 0) then
    ssv(co)%lai%no = 1
    ssv(co)%lai%c(ssv(co)%lai%no)%val = 0.0
    ssv(co)%lai%c(ssv(co)%lai%no)%age = ssp%jday
  elseif (ssp%jday-int(ssv(co)%lai%c(ssv(co)%lai%no)%age) >= lai_comp_length) then
    ssv(co)%lai%no = ssv(co)%lai%no + 1
    ssv(co)%lai%c(ssv(co)%lai%no)%val = 0.0
    ssv(co)%lai%c(ssv(co)%lai%no)%age = ssp%jday
  endif
  ssv(co)%lai%c(ssv(co)%lai%no)%val =  ssv(co)%lai%c(ssv(co)%lai%no)%val + laiinc

endif

end subroutine lai_add





!**********************************************************************!
!                                                                      !
!                      stem_dist :: phenology_methods                  !
!                      ------------------------------                  !
!                                                                      !
!  SUBROUTINE stem_dist(respref,resp)                                  !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Apply stem respiration to the active stems, and sum the
!! resulting stems.
!! @details Apply stem respiration to each of the stem compartments
!! which span the age, for the current cohort.
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine stem_dist(respref,resp)
!**********************************************************************!
integer :: i,co
real(dp) :: ans,resp,respref
!**********************************************************************!

co = ssp%cohort

resp = 0.0

do i=1,ssv(co)%stem%no
  ans = ssv(co)%stem%c(i)%val*respref
  resp = resp + ans
  ssv(co)%stem%c(i)%val = ssv(co)%stem%c(i)%val - ans
enddo

ssv(co)%stem%tot = 0.0
do i=1,ssv(co)%stem%no
  ssv(co)%stem%tot = ssv(co)%stem%tot + ssv(co)%stem%c(i)%val
enddo

end subroutine stem_dist





!**********************************************************************!
!                                                                      !
!                   root_dist :: phenology_methods                     !
!                   ------------------------------                     !
!                                                                      !
!  SUBROUTINE ROOT_DIST(respref,resp)                                  !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Apply root respiration to the active roots, and sum the
!! resulting roots.
!! @details Apply root respiration to each of the root compartments
!! which span the age, for the current cohort.
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine ROOT_DIST(respref,resp)
!**********************************************************************!
integer :: i,co
real(dp) :: ans,resp,respref
!----------------------------------------------------------------------!

co = ssp%cohort

resp = 0.0

do i=1,ssv(co)%root%no
  ans = ssv(co)%root%c(i)%val*respref
  resp = resp + ans
  ssv(co)%root%c(i)%val = ssv(co)%root%c(i)%val - ans
enddo

ssv(co)%root%tot = 0.0
do i=1,ssv(co)%root%no
  ssv(co)%root%tot = ssv(co)%root%tot + ssv(co)%root%c(i)%val
enddo

end subroutine ROOT_DIST





!**********************************************************************!
!                                                                      !
!                     lai_dist :: phenology_methods                    !
!                     -----------------------------                    !
!                                                                      !
! SUBROUTINE LAI_DIST(lmor_sc,leaflit)                                 !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief (Currently not implemented). Apply leaf motality.
!! @details Leaf mortality is applied to the current cohort.
!! The mortality values are taken from 'lmor_sc' and are applied
!! to the lai compartments which span the age of the leaf.
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine lai_dist(lmor_sc,leaflit)
!**********************************************************************!
integer :: i,co,nlcs
real(dp) :: lmor_sc(3600),ans,leaflit
!----------------------------------------------------------------------!

co = ssp%cohort

do i=1,ssv(ssp%cohort)%lai%no
  ans = lmor_sc(int(ssp%jday-ssv(co)%lai%c(i)%age+1.5))*0.0+1.0
  leaflit = leaflit + ssv(co)%lai%c(i)%val*(1.0 - ans)
  ssv(co)%lai%c(i)%val = ssv(co)%lai%c(i)%val*ans
enddo

! Sum up the lai compartments.
ssv(co)%lai%tot = 0.0
do i=1,ssv(co)%lai%no
  ssv(co)%lai%tot = ssv(co)%lai%tot + ssv(co)%lai%c(i)%val
enddo

leaflit = leaflit*12.0/pft(co)%sla/18.0

end subroutine lai_dist



end module phenology_methods
