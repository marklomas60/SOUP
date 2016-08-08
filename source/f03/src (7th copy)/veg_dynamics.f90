module veg_dynamics

use real_precision
use dims
use system_state
use site_parameters
use state_methods
use func
use tuning_parameters

implicit none

contains

!***********************************************************************
!                                                                      *
!                          SUBROUTINE COVER                            *
!                          ----------------                            *
!                                                                      *
!***********************************************************************
subroutine COVER(nft,tmp,prc,firec,fireres,fprob,ftprop,check_closure)
!***********************************************************************
real(dp) :: npp(max_cohorts),nps(max_cohorts),tmp(12,31),prc(12,31),firec,fprob, &
 ftprop(max_cohorts),fri,norm,ftprop0(max_cohorts),total_carbon,old_total_carbon, &
 mtmp,mprc
integer nft,ft,fireres
logical check_closure

!----------------------------------------------------------------------*
! Check carbon closure.
!----------------------------------------------------------------------*
  if (check_closure) call SUM_CARBON(old_total_carbon,.false.)
!----------------------------------------------------------------------*

do ft=1,ssp%cohorts
  npp(ft) = ssv(ft)%npp
  nps(ft) = ssv(ft)%nps
enddo

!----------------------------------------------------------------------*
! Compute the likelyhood of fire in the current year 'fprob'.          *
! 'find' is the fire index                                             *
!----------------------------------------------------------------------*
call FIRE(prc,tmp,fri,fprob)

!----------------------------------------------------------------------*
! Take off area burnt by fire together with plants past there sell by  *
! date and put this as bare ground ready for new growth 'ngrowth'.     *
! Also shift cover and biomass arrays one to the right.                *
!----------------------------------------------------------------------*
call NEWGROWTH(fprob,npp,nps,fireres,firec)

!----------------------------------------------------------------------*
! Set cover arrays to adjust to ftprop as best they can.               *
! ftprop contains the total proportion of that cover, not the          *
! proportion of bare land to assign. Calculate the new ftprop.         *
!----------------------------------------------------------------------*
norm = 0.0
do ft=1,nft
  ftprop(ft) = ftprop(ft)/100.0
  ftprop0(ft) = 0.0
enddo

do ft=1,ssp%cohorts
  ftprop0(pft(ft)%itag)=ftprop0(pft(ft)%itag) + ssv(ft)%cov
enddo

norm = 0.0d0
do ft=1,nft           
  if (ftprop(ft)>0.0d0) then
    ftprop(ft) = ftprop(ft) - ftprop0(ft)
    if (ftprop(ft)<0.0d0) ftprop(ft)=0.0d0 !can't remove cov
    norm = norm + ftprop(ft)
  endif
enddo

if (ssp%new_cov > 0.0) then
  do ft=1,nft
    ftprop(ft) = ftprop(ft)/norm*ssp%new_cov
  enddo
endif

!----------------------------------------------------------------------*
! Check carbon closure.
!----------------------------------------------------------------------*
if (check_closure) then
  call SUM_CARBON(total_carbon,.false.)
  if (abs(total_carbon+firec-old_total_carbon) > 1.0e-6) then
    write(*,*) 'Breach of carbon closure in COVER:',total_carbon+firec-old_total_carbon,' g/m^2.'
  endif
endif

!***********************************************************************
end subroutine COVER
!***********************************************************************





!***********************************************************************
!                                                                      *
!                          INITIALISE_NEW_COHORTS                      *
!                          **********************                      *
!                                                                      *
!***********************************************************************
subroutine INITIALISE_NEW_COHORTS(nft,ftprop,check_closure)
!***********************************************************************
real(dp) :: ftprop(max_cohorts),sumc,sumftprop,total_carbon,old_total_carbon
integer ft,nft,i,cohort
logical check_closure
!----------------------------------------------------------------------*
! Check carbon closure.
!----------------------------------------------------------------------*
  if (check_closure) call SUM_CARBON(old_total_carbon,.false.)
!----------------------------------------------------------------------*

sumftprop = 0.0
do ft=1,nft
  sumftprop = sumftprop + ftprop(ft)
enddo

!----------------------------------------------------------------------*
! Set cover arrays for this years ft proportions, take carbon from     *
! litter to provide nppstore and canopy.                               *
!----------------------------------------------------------------------*
cohort = ssp%cohorts

do ft=1,nft
  if (ftprop(ft)>0.0) then
    cohort = cohort + 1
    ssp%co2ftmap(ft,1) = ssp%co2ftmap(ft,1) + 1
    ssp%co2ftmap(ft,ssp%co2ftmap(ft,1)+1) = cohort
!----------------------------------------------------------------------*
! Set plant functional type parameterisation for the new cohorts.
!----------------------------------------------------------------------*
    pft(cohort) = pft_tab(ft)

!----------------------------------------------------------------------*
! Set initial values of the system state.
!----------------------------------------------------------------------*
    call INITIALISE_STATE_COHORT(cohort)

    ssv(cohort)%stemfr = ssv(1)%stemfr

    ssv(cohort)%nppstore(1) = pft(cohort)%stemx
    ssv(cohort)%nppstore(2) = pft(cohort)%stemx
    ssv(cohort)%nppstore(3) = pft(cohort)%stemx
    if ((pft(cohort)%mort < 5).and.(pft(cohort)%mort > 0)) then
      if (ssp%co2ftmap(ft,1) > 1) then
        ssv(cohort)%nppstore(1) = ssv(ssp%co2ftmap(ft,2))%nppstore(1)
        ssv(cohort)%nppstore(2) = ssv(ssp%co2ftmap(ft,2))%nppstore(2)
        ssv(cohort)%nppstore(3) = ssv(ssp%co2ftmap(ft,2))%nppstore(3)
      endif
    endif

    ssv(cohort)%cov = ftprop(ft)
    ssv(cohort)%ppm = pft_tab(ft)%ppm0
    ssv(cohort)%hgt = 0.004
    ssv(cohort)%age = 1
    call SET_NEW_SOIL_RES(cohort,ftprop(ft)/sumftprop)

!----------------------------------------------------------------------*
! Take carbon from litter to balance the storage.
!----------------------------------------------------------------------*
    ssv(cohort)%slc = ssv(cohort)%slc - ssv(cohort)%nppstore(1)*ssv(cohort)%cov

    if (ssv(cohort)%slc<0.0) then
!----------------------------------------------------------------------*
! Make up any shortfall with soil carbon.
!----------------------------------------------------------------------*
      sumc = 0.0
      do i=1,8
         sumc = sumc + ssv(cohort)%c(i)
      enddo
      do i=1,8
        ssv(cohort)%c(i) = ssv(cohort)%c(i)*(1.0+ssv(cohort)%slc/sumc/ssv(cohort)%cov)
      enddo
      ssv(cohort)%slc = 0.0
    endif
!----------------------------------------------------------------------*

  endif
enddo
ssp%cohorts = cohort
call RESET_SOIL_RES()

!----------------------------------------------------------------------*
! Check carbon closure.
!----------------------------------------------------------------------*
if (check_closure) then
  call SUM_CARBON(total_carbon,.false.)
  if (abs(total_carbon-old_total_carbon) > 1.0e-3) then
    write(*,*) 'Breach of carbon closure in INITIALISE_NEW_COHORTS:', &
 total_carbon-old_total_carbon,' g/m^2.'
  endif
endif

!***********************************************************************
end subroutine INITIALISE_NEW_COHORTS
!***********************************************************************





!***********************************************************************
!                                                                      *
!                          SUBROUTINE GROWTH                           * 
!                          *****************                           *
!                                                                      *
!***********************************************************************
subroutine GROWTH(nft,lai,stembio,rootbio,check_closure)
!***********************************************************************
real(dp) :: npp(max_cohorts),lai(max_cohorts),nps(max_cohorts),npr(max_cohorts), &
 evp(max_cohorts),rootbio,slc(max_cohorts),rlc(max_cohorts), &
 sln(max_cohorts),rln(max_cohorts),stembio, &
 total_carbon,old_total_carbon,ans
integer nft,ftmor(max_cohorts),ft,i
logical check_closure
!----------------------------------------------------------------------*
! Check carbon closure.
!----------------------------------------------------------------------*
  if (check_closure) call SUM_CARBON(old_total_carbon,.false.)
!----------------------------------------------------------------------*

do ft=1,ssp%cohorts
  npp(ft) = ssv(ft)%npp
  nps(ft) = ssv(ft)%nps
  npr(ft) = ssv(ft)%npr
  evp(ft) = ssv(ft)%evp
enddo

!----------------------------------------------------------------------*
! Initialise litter arrays, and add on leaf litter computed in DOLY.   *
!----------------------------------------------------------------------*
do ft=1,ssp%cohorts
  slc(ft) = 0.0
  rlc(ft) = 0.0
  sln(ft) = 0.0
  rln(ft) = 0.0
enddo

!----------------------------------------------------------------------*
! Thin vegetation where npp is not sufficient to maintain sensible     *
! growth rate.                                                         *
!----------------------------------------------------------------------*
call THIN(nft,npp,lai,nps,evp,slc,check_closure)

!----------------------------------------------------------------------*
!Compute leaf root and stem biomasses.                                 *
!----------------------------------------------------------------------*
stembio = 0.0
rootbio = 0.0
do ft=1,ssp%cohorts
  stembio = stembio + ssv(ft)%bio(1)*ssv(ft)%cov
  rootbio = rootbio + ssv(ft)%bio(2)*ssv(ft)%cov
enddo

do ft=1,ssp%cohorts
  ssv(ft)%npp = npp(ft)
  ssv(ft)%nps = nps(ft)
  ssv(ft)%npr = npr(ft)
  ssv(ft)%evp = evp(ft)
enddo

!----------------------------------------------------------------------*
! Check carbon closure.
!----------------------------------------------------------------------*
if (check_closure) then
  call SUM_CARBON(total_carbon,.false.)
  if (abs(total_carbon-old_total_carbon) > 1.0e-3) then
    write(*,*) 'Breach of carbon closure in GROWTH:',total_carbon-old_total_carbon,' g/m^2.'
  endif
endif

!***********************************************************************
end subroutine GROWTH
!***********************************************************************





!***********************************************************************
!                                                                      *
!                          SUBROUTINE THIN                             *
!                          ***************                             *
!                                                                      *
!*********************************************************************** 
subroutine THIN(nft,npp,lai,nps,evp,slc,check_closure)
!***********************************************************************
real(dp) :: ftmat(max_cohorts),npp(max_cohorts), &
 lai(max_cohorts),nps(max_cohorts),evp(max_cohorts),storelit(max_cohorts), &
 slc(max_cohorts),pbionew,pbioold,no,pbio, &
 ftcov(max_cohorts),covnew(max_age,max_cohorts),ppmnew(max_age,max_cohorts),shv,lmv,nv,pi, &
 hwv(max_age),emv,totno,totcov,scale(max_age),oldbio(max_cohorts,2),dimold, &
 g0,gf,gm,grate(max_age),hgtnew(max_age),dimnew,maxhgt,minhgt,hc1,hc2, &
 tcov1,tcov2,sum,xxx,total_carbon,old_total_carbon
integer nft,ft,i,year,coh,ift
logical check_closure
!----------------------------------------------------------------------*
! Check carbon closure.
!----------------------------------------------------------------------*
  if (check_closure) call SUM_CARBON(old_total_carbon,.false.)
!----------------------------------------------------------------------*

pi = 3.1415926

hc1 = 0.05
hc2 = 2.0

do ft=1,ssp%cohorts
  do i=1,2
    oldbio(ft,i) = ssv(ft)%bio(i)
  enddo
  storelit(ft) = 0.0
enddo

do ft=1,ssp%cohorts
  ftmat(ft) = real(pft(ft)%mort)
enddo

!----------------------------------------------------------------------*
! FT loop for thinning and height competition.                         *
!----------------------------------------------------------------------*
do ft=1,nft
  if ((pft_tab(ft)%gr0>0.0).and.(ssp%co2ftmap(ft,1)>0)) then

!----------------------------------------------------------------------*
! Height competition.                                                  *
!----------------------------------------------------------------------*
    minhgt = 10000.0
    maxhgt =-10000.0
    i = 0
    do coh=2,ssp%co2ftmap(ft,1)+1
      ift = ssp%co2ftmap(ft,coh)
      year = int(ssv(ift)%age+.5)
      if (ssv(ift)%hgt>0) i = i + 1
      scale(year) = 1.0 - ssv(ift)%cov
      if (scale(year)<minhgt)  minhgt = scale(year)
      if (scale(year)>maxhgt)  maxhgt = scale(year)
    enddo

    if (maxhgt-minhgt<-1000.0) then
      do year=1,pft_tab(ft)%mort
        scale(year) = (scale(year) - minhgt)/(maxhgt - minhgt)
        scale(year) = (scale(year)*hc1 + 1.0 - hc1)**hc2
      enddo

      tcov1 = 0.0
      tcov2 = 0.0
      do coh=2,ssp%co2ftmap(ft,1)+1
        ift = ssp%co2ftmap(ft,coh)
        year = int(ssv(ift)%age+0.5)
        covnew(year,ift) = ssv(ift)%cov*scale(year)
        tcov1 = tcov1 + ssv(ift)%cov
        tcov2 = tcov2 + covnew(year,ift)
      enddo

      do coh=2,ssp%co2ftmap(ft,1)+1
        ift = ssp%co2ftmap(ft,coh)
        year = int(ssv(ift)%age+.5)
        if (covnew(year,ift)>0.0) then
          covnew(year,ift) = covnew(year,ift)*tcov1/tcov2
          do i=1,2
            ssv(ift)%bio(i) = ssv(ift)%bio(i)*ssv(ift)%cov/covnew(year,ift)
            oldbio(ift,i) = oldbio(ift,i)*ssv(ift)%cov/covnew(year,ift)
          enddo
          ssv(ift)%ppm = ssv(ift)%ppm*ssv(ift)%cov/covnew(year,ift)
          ssv(ift)%cov = covnew(year,ift)
        endif
      enddo

    endif
!----------------------------------------------------------------------*

    g0 = pft_tab(ft)%gr0
    gf = pft_tab(ft)%grf
    gm = real(pft_tab(ft)%mort)/10.0
    do year=1,pft_tab(ft)%mort
      grate(year) = (gf - g0)/gm*real(year - 1) + g0
      if (real(year)>=gm)  grate(year) = gf
    enddo

    totno = 0.0
    totcov = 0.0
    do coh=2,ssp%co2ftmap(ft,1)+1
      ift = ssp%co2ftmap(ft,coh)
      totno = totno + ssv(ift)%cov*ssv(ift)%ppm
      totcov = totcov + ssv(ift)%cov
    enddo
    totno = totno/totcov

!----------------------------------------------------------------------*
! Compute cover and ppm to sustain a minimum growth rate, put these    *
! values in covnew and ppmnew.                                         *
!----------------------------------------------------------------------*
    do coh=2,ssp%co2ftmap(ft,1)+1
      ift = ssp%co2ftmap(ft,coh)

      year = int(ssv(ift)%age+.5)

      shv = ssv(ift)%stem%tot
      emv = evp(ift)/3600.0/1000.0
      lmv = lai(ift)*1.3/(lai(ift) + 3.0)
      lmv = 1.0
      nv = 1.0

      if (ssv(ift)%ppm*ssv(ift)%cov>0.0) then
!----------------------------------------------------------------------*
! Calculate the increase in diameter produced by stem NPP 'nps'.       *
!----------------------------------------------------------------------*

        if (emv>0.0) then
! hwv = theoretical maximum height (hydrolics)
          hwv(year) = pft(ift)%pdif*lmv/nv*(pft(ift)%xyl*shv/emv/pft(ift)%wden/10000.0)**0.5
        else
          hwv(year) = ssv(ift)%hgt*0.9
        endif

!----------------------------------------------------------------------*
! Calculate new height.
!----------------------------------------------------------------------*
        if (hwv(year)>0.0) then
          hgtnew(year) = (hwv(year)-ssv(ift)%hgt)/hwv(year)*0.5
        else
          hgtnew(year) = 0.0
        endif
        if (hgtnew(year)<0.0)  hgtnew(year) = 0.0
        hgtnew(year) = hgtnew(year) + ssv(ift)%hgt

!----------------------------------------------------------------------*
! pbio = g/individual
!----------------------------------------------------------------------*
        pbioold = ssv(ift)%bio(1)/ssv(ift)%ppm
        pbionew = (ssv(ift)%bio(1) + ssv(ift)%stem%tot*(1.0 - &
 stlit(real(year,dp),ftmat(ift))))/ssv(ift)%ppm

!----------------------------------------------------------------------*
! Old diameter
!----------------------------------------------------------------------*
        if (ssv(ift)%hgt>0.0) then
! dimold = 2.0*(pbioold/1000000.0/ssv(ift)%hgt/
          dimold = 2.0*(pbioold/1000000.0/hgtnew(year)/pi/pft(ift)%wden)**0.5
        else
          dimold = 0.0
        endif

!----------------------------------------------------------------------*
! New diameter
!----------------------------------------------------------------------*
        if (hgtnew(year)>0.0) then
          dimnew = 2.0*(pbionew/1000000.0/hgtnew(year)/pi/pft(ift)%wden)**0.5
        else
          dimnew = 0.0
        endif

        if ((dimnew-dimold)/2.0>grate(year)) then
!----------------------------------------------------------------------*
! No thinning required.                                                *
!----------------------------------------------------------------------*
!          print*,'no thinning'
          ppmnew(year,ift) = ssv(ift)%ppm
          covnew(year,ift) = ssv(ift)%cov
        else
!----------------------------------------------------------------------*
! Thinning required.                                                   *
!----------------------------------------------------------------------*
!          print*,'thinning'
          xxx = (grate(year)+dimold/2.0)**2.0*1000000.0*hgtnew(year)*pi*pft(ift)%wden
          ppmnew(year,ift) = (ssv(ift)%bio(1) + ssv(ift)%stem%tot* &
 (1.0 - stlit(real(year,dp),ftmat(ift)))/100.0)/xxx
          covnew(year,ift) = ssv(ift)%cov

        endif
      else
        pbio = 0.0
        scale(year) = 0.0
        ppmnew(year,ift) = 0.0
        covnew(year,ift) = 0.0
        hgtnew(year) = 0.0
      endif
    enddo

!----------------------------------------------------------------------*
! Correct biomass array, and adjust litter for any thinned trees.      *
!----------------------------------------------------------------------*
    do coh=2,ssp%co2ftmap(ft,1)+1
      ift = ssp%co2ftmap(ft,coh)
      year = int(ssv(ift)%age+0.5)
      if ((ssv(ift)%ppm>0.0).and.(ssv(ift)%cov>0.0)) then
        no = ssv(ift)%ppm*ssv(ift)%cov - ppmnew(year,ift)*covnew(year,ift)

        do i=1,2
          pbio = ssv(ift)%bio(i)/ssv(ift)%ppm
          slc(ift) = slc(ift) + pbio*no

          if (ssv(ift)%cov>0.0) then
            ssv(ift)%bio(i) = (oldbio(ift,i)*ssv(ift)%cov - pbio*no)/covnew(year,ift)
            ssv(ift)%bio(i) = oldbio(ift,i)*ppmnew(year,ift)/ssv(ift)%ppm
          endif
        enddo
        storelit(ift) = storelit(ift) + ssv(ift)%nppstore(1)/ssv(ift)%ppm*no
        ssv(ift)%cov = covnew(year,ift)
        ssv(ift)%ppm = ppmnew(year,ift)
        ssv(ift)%hgt = hgtnew(year)
      else
        ssv(ift)%cov = 0.0
        ssv(ift)%ppm = 0.0
        ssv(ift)%hgt = 0.0
      endif

    enddo
!----------------------------------------------------------------------*
  else
    do coh=2,ssp%co2ftmap(ft,1)+1
      ift = ssp%co2ftmap(ft,coh)
      ssv(ift)%ppm = 0.0
      ssv(ift)%hgt = 0.0
    enddo
  endif
!----------------------------------------------------------------------*
! End of ft loop.
!----------------------------------------------------------------------*
enddo

do ft=1,nft
  ftcov(ft) = 0.0
enddo

do ft=1,ssp%cohorts
  ftcov(pft(ft)%itag) = ftcov(pft(ft)%itag) + ssv(ft)%cov
enddo

!----------------------------------------------------------------------*
! Correct nppstore to account for thinning.                            *
!----------------------------------------------------------------------*
do ft=1,ssp%cohorts
  if ((pft(ft)%gr0>0.0).and.(ssv(ft)%nppstore(1)>0.0).and.(ssv(ft)%cov>0.0)) then
    slc(ft) = slc(ft) +  storelit(ft)
    ssv(ft)%nppstore(1) = (ssv(ft)%nppstore(1)*ssv(ft)%cov - storelit(ft))/ssv(ft)%cov
  endif
enddo

do ft=1,ssp%cohorts
  ssv(ft)%slc = slc(ft)
enddo

!----------------------------------------------------------------------*
! Check carbon closure.
!----------------------------------------------------------------------*
if (check_closure) then
  call SUM_CARBON(total_carbon,.false.)
  if (abs(total_carbon-old_total_carbon) > 1.0e-3) then
    write(*,*) 'Breach of carbon closure in THINNING:',total_carbon-old_total_carbon,' g/m^2.'
    stop
  endif
endif

!***********************************************************************
end subroutine THIN
!***********************************************************************





!***********************************************************************
!                                                                      *
!                             FUNCTION find                            *
!                             *************                            *
!                                                                      *
!***********************************************************************
function find(tmp,prc)
!***********************************************************************
real(dp) :: find,tmp(12),prc(12)
real(dp) :: tmplim,prclim
integer i
!----------------------------------------------------------------------*

tmplim = -5.0
prclim = 50.0

find = 0.0
do i=1,12
  if (tmp(i)>tmplim) then
    if (prc(i)<prclim) then
      find = find + prc(i)/(12.0*prclim)
    else
      find = find + 1/12.0
    endif
  else
    find = find + 1/12.0
  endif
enddo

!***********************************************************************
end function find
!***********************************************************************





!***********************************************************************
!                                                                      *
!                             SUBROUTINE c3c4                          *
!                             ***************                          *
!                                                                      *
!***********************************************************************
subroutine c3c42(ftprop,npp,range)
!***********************************************************************
real(dp) :: ftprop(max_cohorts),npp(max_cohorts),range
real(dp) :: grass,nd

grass = ftprop(2)
nd = npp(2) - npp(3)
ftprop(2) = grass*nd/(2.0*range) + grass/2.0
if (ftprop(2)>grass)  ftprop(2) = grass
if (ftprop(2)<0.0)  ftprop(2) = 0.0
ftprop(3) = grass - ftprop(2)

!***********************************************************************
end subroutine c3c42
!***********************************************************************





!***********************************************************************
!                                                                      *
!                             SUBROUTINE c3c4                          *
!                             ***************                          *
!                                                                      *
!***********************************************************************
subroutine c3c4(ftprop,c3old,c4old,npp,nps)
!***********************************************************************
real(dp) :: ftprop(max_cohorts),npp(max_cohorts),nps(max_cohorts),c3old,c4old
real(dp) :: grass,c3p,c4p,adj

if (c3old+c4old>0.0) then
  c3p = c3old/(c3old + c4old)
  c4p = c4old/(c3old + c4old)
else
  c3p = 0.5
  c4p = 0.5
endif

if (npp(2)*nps(2)+npp(3)*nps(3)>0.0) then
  adj = npp(2)*nps(2)/(npp(2)*nps(2) + npp(3)*nps(3)) - 0.5
else
  adj = 0.0
endif

c3p = c3p + adj
c4p = c4p - adj

if (c3p<0.0) then
  c3p = 0.0
  c4p = 1.0
endif

if (c4p<0.0) then
  c4p = 0.0
  c3p = 1.0
endif

grass = ftprop(2) + ftprop(3)

ftprop(2) = grass*c3p
ftprop(3) = grass*c4p

!***********************************************************************
end subroutine c3c4
!***********************************************************************





!***********************************************************************
!                                                                      *
!                             SUBROUTINE GRASSREC                      *
!                             *******************                      *
!                                                                      *
!***********************************************************************
subroutine GRASSREC(nft,ftprop,gold,x,nat_map)
!***********************************************************************
real(dp) :: ftprop(max_cohorts),gold,x,ntcov,ftt,ftpropo(max_cohorts)
integer nft,ft,nat_map(8)
!----------------------------------------------------------------------*

ntcov = 0.0
ftt = 0.0
do ft=5,8
  ntcov = ntcov + ftprop(nat_map(ft))*ssp%new_cov/100.0
  ftt = ftt + ftprop(nat_map(ft))
enddo

if (ntcov>gold*x) then
  do ft=4,nft
    ftpropo(ft) = ftprop(ft)
    ftprop(ft) = 100.0*gold*x*ftprop(ft)/(ftt*ssp%new_cov)
    ftprop(2) = ftprop(2) + ftpropo(ft) - ftprop(ft)
  enddo

  ntcov = 0.0
  do ft=4,nft
    ntcov = ntcov + ftprop(ft)*ssp%new_cov/100.0
  enddo
  if (abs(ntcov-gold*x)>0.000001) write(11,'(''Treerec subroutine error'')')
endif

!***********************************************************************
end subroutine GRASSREC
!***********************************************************************





!***********************************************************************
!                                                                      *
!                             SUBROUTINE BAREREC                       *
!                             ******************                       *
!                                                                      *
!***********************************************************************
subroutine BAREREC(ftprop,bpaold,x)
!***********************************************************************
real(dp) :: ftprop(max_cohorts),bpaold,x
real(dp) :: nbp,cgcov
!----------------------------------------------------------------------*

! Current growth 'cgcov'
cgcov = 1.0 - ssp%new_cov

nbp = ssp%new_cov*ftprop(1)/100.0

if (nbp<bpaold*(1.0-x)) then
  ssp%new_cov = 1.0 - bpaold*(1.0 - x) - cgcov
  ssv(1)%cov = bpaold*(1.0 - x)
else
  ssv(1)%cov = ssp%new_cov*ftprop(1)/100.0
  ssp%new_cov = ssp%new_cov*(1.0 - ftprop(1)/100.0)
endif

!***********************************************************************
end subroutine BAREREC
!***********************************************************************





!***********************************************************************
!                                                                      *
!                            SUBROUTINE ADDBIO                         *
!                            *****************                         *
! Add on biomass to bio array.                                         *
!***********************************************************************
subroutine ADDBIO(npp,nps,npr)
!***********************************************************************
real(dp) :: npp(max_cohorts),nps(max_cohorts),npr(max_cohorts),npps,nppr
integer ft
!----------------------------------------------------------------------*

do ft=1,ssp%cohorts
  npps = nps(ft)*npp(ft)/100.0
  nppr = npr(ft)*npp(ft)/100.0
  ssv(ft)%bio(1) = ssv(ft)%bio(1) + npps
  ssv(ft)%bio(2) = ssv(ft)%bio(2) + nppr
enddo

!***********************************************************************
end subroutine ADDBIO
!***********************************************************************





!***********************************************************************
!                                                                      *
!                            SUBROUTINE PERC                           *
!                            ***************                           *
!                                                                      *
!***********************************************************************
subroutine PERC(fttags,dof,tmin,ftprop,nft)
!***********************************************************************
real(dp) :: dof(max_cohorts),tmin,ftprop(max_cohorts),dec1
real(dp) :: evbl,evnl,dec,dcnl,drdc,sum
real(dp) :: t1eb,t2eb,t3eb,d1eb,d2eb,d3eb
real(dp) :: t1en,t2en,t3en,d1en,d2en,d3en
real(dp) :: t1d,t2d,t3d,d1d,d2d,d3d
real(dp) :: t1dd,t2dd,t3dd,d1dd,d2dd,d3dd
real(dp) :: t1dn,t2dn,t3dn,d1dn,d2dn,d3dn
real(dp) :: t1d1,t2d1,t3d1,d1d1,d2d1,d3d1,ndof
integer i,ft,nft
character fttags(100)*1000,st1*1000
!----------------------------------------------------------------------*

t1eb = 100.0   ! 0
t2eb = -10.0   ! 1
t3eb = -25.0   ! 0
d1eb =  -1.0   ! 0
d2eb =   0.0   ! 1
d3eb = 150.0   ! 0

t1en = -20.0   ! 0
t2en = -40.0   ! 1
t3en = -60.0   ! 0
d1en = 100.0   ! 0
d2en = 180.0   ! 1
d3en = 260.0   ! 0

t1d =  10.0    ! 0
t2d = -15.0    ! 1  broadleaf
t3d = -50.0    ! 0
d1d =-999.0    ! 0
d2d = 100.0    ! 1
d3d = 140.0    ! 0

t1d1 = 100.0   ! 0
t2d1 =  20.0   ! 1  broadleaf
t3d1 =  10.0   ! 0
d1d1 = 150.0   ! 0
d2d1 = 200.0   ! 1
d3d1 = 400.0   ! 0

t1dd =  12.0   ! 0
t2dd = -15.0   ! 1  broadleaf
t3dd = -30.0   ! 0
d1dd =  75.0   ! 0
d2dd = 200.0   ! 1
d3dd = 400.0   ! 0

t1dn = -60.0   ! 0
t2dn = -70.0   ! 1
t3dn = -80.0   ! 0
d1dn = 180.0   ! 0
d2dn = 240.0   ! 1
d3dn = 300.0   ! 0

nft = 1
st1='Dc_Nl'
ndof = dof(ntags(fttags,st1))
evbl = triscale2(tmin,ndof,t1eb,t2eb,t3eb,d1eb,d2eb,d3eb)
evnl = triscale2(tmin,ndof,t1en,t2en,t3en,d1en,d2en,d3en)
dec  = triscale2(tmin,ndof,t1d ,t2d ,t3d ,d1d ,d2d ,d3d )
drdc = triscale2(tmin,ndof,t1dd,t2dd,t3dd,d1dd,d2dd,d3dd)
dcnl = triscale2(tmin,ndof,t1dn,t2dn,t3dn,d1dn,d2dn,d3dn)
dec1 = triscale2(tmin,ndof,t1d1,t2d1,t3d1,d1d1,d2d1,d3d1)

!IF (tex-1.LT.0.01)  dcnl = 0.0

sum = evbl + evnl + dec + drdc + dcnl + dec1
if (sum>0.001) then
  ftprop(1) = 0.0
  ftprop(2) = 0.0
  ftprop(3) = 0.0
  st1='Ev_Bl'
  ftprop(ntags(fttags,st1)) = 100.0*evbl/sum
  st1='Ev_Nl'
  ftprop(ntags(fttags,st1)) = 100.0*evnl/sum
  st1='Dc_Bl'
  ftprop(ntags(fttags,st1)) = 100.0*(dec + drdc + dec1)/sum
  st1='Dc_Nl'
  ftprop(ntags(fttags,st1)) = 100.0*dcnl/sum
else
  do i=1,nft
    ftprop(i) = 0.0
  enddo
  ftprop(2) = 100.0
endif

ftprop(1) = 100.0
do ft=2,nft
  ftprop(1) = ftprop(1) - ftprop(ft)
enddo

!***********************************************************************
end subroutine PERC
!***********************************************************************





!***********************************************************************
!                                                                      *
!                            SUBROUTINE MKLIT                          *
!                            ****************                          *
! Make litter.                                                         *
!***********************************************************************
subroutine MKLIT(ftmat,slc,rlc,sln,rln,npp,nps,npr)
!***********************************************************************
real(dp) :: ftmat(max_cohorts),slc(max_cohorts),rlc(max_cohorts),sln(max_cohorts), &
 rln(max_cohorts),npp(max_cohorts),nps(max_cohorts),npr(max_cohorts),npps,nppr,sl,rl
integer ft
!----------------------------------------------------------------------*

do ft=1,ssp%cohorts
  if ((pft(ft)%gr0>0.0).and.(ssv(ft)%age>1)) then
    npps = nps(ft)*npp(ft)/100.0
    nppr = npr(ft)*npp(ft)/100.0
    sl = stlit(ssv(ft)%age,ftmat(ft))
    rl = 0.8
    slc(ft) = slc(ft) + sl*npps*ssv(ft)%cov
    rlc(ft) = rlc(ft) + rl*nppr*ssv(ft)%cov
    ssv(ft)%bio(1) = ssv(ft)%bio(1) - sl*npps
    ssv(ft)%bio(2) = ssv(ft)%bio(2) - rl*nppr
  endif
enddo

do ft=1,ssp%cohorts
  sln(ft) = 0.0
  rln(ft) = 0.0
enddo

!***********************************************************************
end subroutine MKLIT
!***********************************************************************





!***********************************************************************
!***********************************************************************
function stlit(age,mat)
!***********************************************************************
real(dp) :: stlit,mat,temp,age
!----------------------------------------------------------------------*

stlit = 0.9*real(age)/mat + 0.1
if (stlit>1.0)  stlit = 1.0

if (age<mat) then
  stlit = 0.0
else
  stlit = 1.0
endif

temp = age/mat
if (temp>1.0) temp = 1.0

stlit = temp**0.5
stlit = stlit*0.6
!stlit = 0.0

!***********************************************************************
end function stlit
!***********************************************************************





!***********************************************************************
!                                                                      *
!                            SUBROUTINE VEGMAT                         *
!                            *****************                         *
! Compute the years to reach maturity for each of the fts.             *
!***********************************************************************
subroutine VEGMAT(nft,ftmor,ftwd,ftxyl,ftpd,evp,lai,npp,nps,ftmat)
!----------------------------------------------------------------------*
real(dp) :: ftwd(max_cohorts),ftxyl(max_cohorts),ftpd(max_cohorts),lai(max_cohorts)
real(dp) :: npp(max_cohorts),nps(max_cohorts),evp(max_cohorts),ftmat(max_cohorts),emxv,wd
real(dp) :: pd,lv,fs,shv,pi,lmvt,nvt,hwvt,dvt,massvt,minvt,ppvt
integer ft,nft,ftmor(max_cohorts)
!----------------------------------------------------------------------*

pi = 3.14159
do ft=2,nft
  if (npp(ft)>1.0e-6) then
    emxv = evp(ft)/3600.0/1000.0
    lv = lai(ft)
    fs = nps(ft)/100.0
    shv = fs*npp(ft)/100.0

    lmvt = lv*1.3/(lv + 3.0)
    nvt = 1.0 + 26.0*exp(-0.9*lv)
    nvt = 1.0

    wd = ftwd(ft)
    pd = ftpd(ft)

    hwvt = pd*lmvt*sqrt(ftxyl(ft)*shv/(emxv*wd*10000.0))/nvt
    dvt = 0.0028*hwvt**1.5
    massvt = pi*(dvt/2.0)**2.0*hwvt*wd*nvt
    minvt = pi*((dvt + 0.001/nvt)/2.0)**2.0*hwvt*wd*nvt - massvt
    ppvt = shv/minvt
    ftmat(ft) = massvt*ppvt/shv
    ftmat(ft) = 1.0
  else
    ftmat(ft) = 0.0
  endif
  if (ftmor(ft)<ftmat(ft)) ftmat(ft) = ftmor(ft)
enddo

!***********************************************************************
end subroutine VEGMAT
!***********************************************************************





!***********************************************************************
!                                                                      *
!                            SUBROUTINE FIRE                           *
!                            ***************                           *
! Compute the fire return interval 'fri' and the probability of a fire *
! in the current year 'fprob'.                                         *
!***********************************************************************
subroutine FIRE(dprc,dtmp,fri,fprob)
!***********************************************************************
real(dp) :: fri,fprob,dtmp(12,31),dprc(12,31),prct(12),prco,lim1,lim2,maxfri
real(dp) :: pow,tlim1,tlim2,totp,indexx,tadj,weight,tmp(12),prc(12)
integer k,no,ind,i,j
!----------------------------------------------------------------------*

do i=1,12
  tmp(i) = 0.0
  prc(i) = 0.0
  do j=1,30
    tmp(i) = tmp(i) + dtmp(i,j)
    prc(i) = prc(i) + dprc(i,j)
  enddo
  tmp(i) = tmp(i)/30
enddo

!----------------------------------------------------------------------*
! Fire model parameters.                                               *
!----------------------------------------------------------------------*
no = 3
lim1 = 150.0
lim2 = 50.0
tlim1 =  0.0
tlim2 = -5.0
weight = 0.5
maxfri = 800.0
pow = 3.0

!----------------------------------------------------------------------*
! Adjust for temperature.                                             *
!----------------------------------------------------------------------*
do i=1,12
  if (tmp(i)>tlim1) then
    tadj = 0.0
  elseif (tmp(i)<tlim2) then
    tadj = lim1
  else
    tadj = (tlim1 - tmp(i))/(tlim1 - tlim2)*lim1
  endif
  prct(i) = min(lim1,prc(i) + tadj)
enddo

!----------------------------------------------------------------------*
! Calculate yearly componet of index.                                  *
!----------------------------------------------------------------------*
totp = 0.0
do k=1,12
  totp = totp + prct(k)/lim1/12.0
enddo

!----------------------------------------------------------------------*
! Calculate monthly component of index.                                *
!----------------------------------------------------------------------*
prco = 0.0
do k=1,no
  ind = minx(prct)
  prco = prco + min(lim2,prct(ind))/real(no)/lim2
  prct(ind) = lim2
enddo

!----------------------------------------------------------------------*
! Compute fire return interval, and convert to probability.            *
!----------------------------------------------------------------------*
indexx = weight*(prco) + (1.0 - weight)*totp

fri = indexx**pow*maxfri

if (fri<2.0)  fri = 2.0
fprob = (1.0 - exp(-1.0/fri))*tgp%p_fprob

!***********************************************************************
end subroutine FIRE
!***********************************************************************





!***********************************************************************
!                                                                      *
!                            SUBROUTINE NEWGROWTH                      *
!                            ********************                      *
! Compute newgrowth and alter cover array accordingly.                 *
!***********************************************************************
subroutine NEWGROWTH(fprob,npp,nps,fireres,firec)
!***********************************************************************
real(dp) ::  fprob,npp(max_cohorts),nps(max_cohorts),tmor,tmor0,npp0,firec,xfprob
integer ft,fireres
!----------------------------------------------------------------------*
xfprob = fprob

firec = 0.0
!----------------------------------------------------------------------*
! Take away veg that has died of old age ie > than pft%mort, and       *
! age the veg by one year.                                             *
!----------------------------------------------------------------------*
do ft=1,max_pftps
  ssp%co2ftmap(ft,1) = 0
enddo

do ft=1,ssp%cohorts
  if (ssv(ft)%age > real(pft(ft)%mort)-0.1) then
    call ACCUMULATE_DIST_SOIL_RES(ft,1.0_dp)
  else
!----------------------------------------------------------------------*
! Age the state structure.                                             *
!----------------------------------------------------------------------*
    ssv(ft)%age = ssv(ft)%age + 1.0
  endif
enddo

!----------------------------------------------------------------------*
! Remove dead cohorts: cover = 0.                                      *
!----------------------------------------------------------------------*
call COMPRESS_STATE()

!----------------------------------------------------------------------*
! Take away veg that is burnt or has died through a hard year ie small *
! LAI.                                                                 *
!----------------------------------------------------------------------*
do ft=1,ssp%cohorts

!----------------------------------------------------------------------*
! 'tmor' is the mortality rate of the forrest based on 'npp'.          *
!----------------------------------------------------------------------*
  npp0 = 0.2
  tmor0 = 6.0
  if (npp(ft)/100.0<0.1) then
    tmor = 1.0
  elseif (npp(ft)/100.0<3.0) then
    tmor = (0.04 - npp0)/3.0*(npp(ft)/100.0) + npp0
  elseif (npp(ft)/100.0<tmor0) then
    tmor = 0.04/(3.0 - tmor0)*(npp(ft)/100.) + &
 3.0*0.04/(tmor0 - 3.0) + 0.04
  else
    tmor = 0.0
  endif
  if (tmor<0.01)  tmor = 0.01

!IF (npp(ft)*nps(ft)/100.0.GT.10.0) THEN
    tmor = 0.002
!ELSE
!  tmor = 1.0 - npp(ft)*nps(ft)/1000.0
!ENDIF
  if (ssv(ft)%age > real(pft(ft)%mort)/2.0) then
    tmor = (1.0 - 2.0*(pft(ft)%mort - ssv(ft)%age)/pft(ft)%mort)/10.0
  else
    tmor = 0.002
  endif
!  tmor = 0.0

!----------------------------------------------------------------------*
! kill off trees with no storage left.                                 *
!----------------------------------------------------------------------*
  if (.not.(ssv(ft)%nppstore(1))>0.1) then
    tmor = 1.0
  endif
!----------------------------------------------------------------------*
!  IF (ssp%iyear == 1)  tmor = 1.0

  if (ssv(ft)%age<real(fireres)) then
    fprob = xfprob
  else
    fprob = 0.0
  endif
  if (fireres<0) fprob = real(-fireres)/1000.0

  call ACCUMULATE_DIST_SOIL_RES(ft,tmor)
  call ACCUMULATE_BURNT_SOIL_RES(ft,fprob,firec)

enddo

call COMPRESS_STATE()

!***********************************************************************
end subroutine NEWGROWTH
!***********************************************************************





!***********************************************************************
!                                                                      *
!                          NATURAL_VEG                                 *
!                          ***********                                 *
!                                                                      *
!***********************************************************************
subroutine NATURAL_VEG(tmp,prc,ftprop,nat_map)
!***********************************************************************
real(dp) :: tmp(12,31),prc(12,31),ftprop(max_cohorts),min_tmp,av_mnth(12)
real(dp) :: x(6,6),evap,tran,bucket,prc_ind
real(dp) :: sum,prop(6)
integer day,mnth,ft,nat_map(8),vdaysoff(5),i,year,daysoff
logical leaves
!----------------------------------------------------------------------*

data (x(1,i),i=1,6)/ -15.0, -10.0, 100.0,   -1.0,   0.0, 150.0/!EB
data (x(2,i),i=1,6)/ -70.0, -45.0, -20.0,  -10.0,   0.0, 350.0/!EN
data (x(3,i),i=1,6)/ -50.0, -15.0,  10.0, -999.0, 100.0, 140.0/!DB
data (x(4,i),i=1,6)/  10.0,  20.0, 100.0,  150.0, 200.0, 400.0/!DB
data (x(5,i),i=1,6)/ -30.0, -15.0,  12.0,  75.0,  200.0, 400.0/!DB
data (x(6,i),i=1,6)/ -80.0, -70.0, -60.0, 200.0,  280.0, 300.0/!DN

!----------------------------------------------------------------------*
! Find average monthly min, and convert to an absolute min.            *
!----------------------------------------------------------------------*
min_tmp = 1000.0
do mnth=1,12
  av_mnth(mnth) = 0.0
  do day=1,30
    if (min_tmp>tmp(mnth,day)) min_tmp = tmp(mnth,day)
    av_mnth(mnth) = av_mnth(mnth) + tmp(mnth,day)     
  enddo
  av_mnth(mnth) = av_mnth(mnth)/30.0     
enddo
min_tmp = minx(av_mnth)*1.29772 - 19.5362

!----------------------------------------------------------------------*
! Find days off index. Based on an imaginary bucket.                   *
!----------------------------------------------------------------------*
evap = 0.5
tran = 0.5

do i=1,5
  vdaysoff(i) = 0
enddo
leaves = .false.

bucket = 0.0

do year=1,10
  daysoff = 0
  do mnth=1,12
    do day=1,30
      bucket = bucket + prc(mnth,day)
      if (leaves) then
        bucket = bucket - evap - tran
      else
        bucket = bucket - evap
      endif
      if (bucket<0.0) then
        leaves = .false.
        bucket = 0.0
      endif 
      if (bucket>ssp%soil_depth/2.0) leaves = .true.
      if (.not.(leaves)) daysoff = daysoff + 1
      if (bucket>ssp%soil_depth) bucket = ssp%soil_depth
    enddo
  enddo
  do i=1,4
    vdaysoff(6-i) = vdaysoff(5-i)
  enddo
  vdaysoff(1) = daysoff
enddo
prc_ind = real(vdaysoff(1)+vdaysoff(2)+vdaysoff(3)+vdaysoff(4))/4.0

!----------------------------------------------------------------------*
! Compute fractions of natural vegetation appropriate to seed new      *
! ground.                                                              *
!----------------------------------------------------------------------*
do i=1,6
  prop(i) =  triscale(min_tmp,x(i,1),x(i,2),x(i,3))*triscale(prc_ind,x(i,4),x(i,5),x(i,6))
enddo

!----------------------------------------------------------------------*
! Map the natural functional types to the list given in the input file.*
!----------------------------------------------------------------------*
ftprop(nat_map(1)) = 0.0
ftprop(nat_map(2)) = 0.0
ftprop(nat_map(3)) = 0.0
ftprop(nat_map(4)) = 0.0
ftprop(nat_map(5)) = prop(1)
ftprop(nat_map(6)) = prop(2)
ftprop(nat_map(7)) = prop(3) + prop(4) + prop(5)
ftprop(nat_map(8)) = prop(6)

!----------------------------------------------------------------------*
! Normalise the fractions.                                             *
!----------------------------------------------------------------------*
sum = 0.0
do ft=1,8
  sum = sum + ftprop(nat_map(ft))
enddo
if (sum>0.0) then
  do ft=1,8
    ftprop(nat_map(ft)) = 100.0*ftprop(nat_map(ft))/sum
  enddo
else
  ftprop(nat_map(3)) = 100.0
endif

!----------------------------------------------------------------------*

!----------------------------------------------------------------------*
!      DO k=-80,100
!        DO j=0,360
!        min_tmp =  (real(k) + 19.5632)/1.29772
!        prc_ind = real(j)
!      DO i=1,6
!        ftprop(i) =  triscale(min_tmp,x(i,1),x(i,2),x(i,3))*
!     &               triscale(prc_ind,x(i,4),x(i,5),x(i,6))
!      ENDDO
!      sum=ftprop(1)+ftprop(2)+ftprop(3)+ftprop(4)+ftprop(5)+ftprop(6)
!        ENDDO
!      ENDDO
!      STOP
!----------------------------------------------------------------------*


!***********************************************************************
end subroutine NATURAL_VEG
!***********************************************************************

end module veg_dynamics

