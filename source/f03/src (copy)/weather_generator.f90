module weather_generator

use real_precision
use data
use metdos
implicit none

contains

!----------------------------------------------------------------------*
!     Interface between ForestETP Weather Generator and SDGVMd
!     Written by Ghislain Picard 30/06/03
!
!----------------------------------------------------------------------*
subroutine EX_CLIM_WEATHER_GENERATOR(stinput,ststats,lat,lon, &
 xlatf,xlatres,xlatresn,xlon0,xlonres,xlonresn,yr0,yrf,tmpv,humv, &
 prcv,mcldv,isite,year0,yearf,du,seed1,seed2,seed3,l_clim,l_stats)
!----------------------------------------------------------------------*
real(dp) :: lat,lon,xlon0,xlatf,xlatres,xlonres
! daily data
integer :: tmpv(500,12,31),humv(500,12,31),prcv(500,12,31)
! monthly data m=monthly
integer :: mtmpvv(4,4,500,12),mhumvv(4,4,500,12),mprcvv(4,4,500,12),mcldvv(4,4,500,12), &
 row,col,cld_default,mtmpv(500,12),mhumv(500,12),mprcv(500,12),mcldv(500,12)
! monthly parameter for the weather generator
real(dp) :: mstdtmpv(12),mstdhumv(12),mraindays(12)
! seed for the random generator of the weather generator
integer :: seed1,seed2,seed3,kode
! local variables
real(dp) :: rain(31),sumprc,correction,ans,tmp_yesterday, hum_yesterday,xx(4,4)
real(dp) :: tmp_autocorr_month,hum_autocorr_month,rrow,rcol,ynorm,xnorm
integer :: year,year0,yearf,recn,ans2(1000),siteno(4,4),i,du,nyears,yr0,mnth, &
 day,isite,yrf,recl1,recl2,xlatresn,xlonresn,fno,nbrain,indx(4,4),ii,jj
character :: fname1*1000,fname2*1000,fname3*1000,stinput*1000,fname4*1000, &
 fname5*1000,fname6*1000,fname7*1000,ststats*1000
logical :: l_clim,l_stats

l_stats = .true.
cld_default = 50

if (seed1==0) seed1 = 1
if (seed2==0) seed2 = 1
if (seed3==0) seed3 = 1

if (du==1) then
!  This works for ftn95
!  recl1 = 4*12
!  recl2 = 6*xlonresn
  recl1 = 4*12 + 2
  recl2 = 6*xlonresn + 2
else
  recl1 = 4*12 + 1
  recl2 = 6*xlonresn + 1
endif

nyears = yearf - year0 + 1

!----------------------------------------------------------------------*
! Find the real(dp) row col corresponding to lat and lon.                  *
!----------------------------------------------------------------------*

rrow = (xlatf - lat)/xlatres
rcol = (lon - xlon0)/xlonres

ynorm = rrow - real(int(rrow))
xnorm = rcol - real(int(rcol))
!----------------------------------------------------------------------*

fno = 90

open(fno+1,file=stinput(1:blank(stinput))//'/maskmap.dat', &
 access='DIRECT',recl=recl2,form='formatted',status='old')

!----------------------------------------------------------------------*
! Find the four sites closest to the desired.                          *
!----------------------------------------------------------------------*
do ii=1,4
  do jj=1,4
    row = int(rrow)+jj-1
    col = int(rcol)+ii-1
    if ((row>=1).and.(row<=xlatresn).and.(col>=1).and. &
 (col<=xlonresn)) then
      recn = row
      read(fno+1,'(720i6)',rec=recn) (ans2(i),i=1,xlonresn)
      siteno(ii,jj) = ans2(col)
      if (siteno(ii,jj)>0) then
        indx(ii,jj) =  1
      else
        indx(ii,jj) =  0
      endif
    else
      indx(ii,jj) = -1
    endif
  enddo
enddo
close(fno+1)
!----------------------------------------------------------------------*

isite = 0
if ((indx(2,2)==1).or.(indx(2,3)==1).or.(indx(3,2)==1).or. &
 (indx(3,3)==1)) then

  l_clim = .true.

  isite = 1

  write(fname1,'(100a)') (stinput(i:i),i=1,blank(stinput)),'/tmp.dat'
  write(fname2,'(100a)') (stinput(i:i),i=1,blank(stinput)),'/hum.dat'
  write(fname3,'(100a)') (stinput(i:i),i=1,blank(stinput)),'/prc.dat'
  write(fname7,'(100a)') (stinput(i:i),i=1,blank(stinput)),'/cld.dat'

  write(fname4,'(100a)') (ststats(i:i),i=1,blank(ststats)),'/stdtmp.dat'
  write(fname5,'(100a)') (ststats(i:i),i=1,blank(ststats)),'/stdhum.dat'
  write(fname6,'(100a)') (ststats(i:i),i=1,blank(ststats)),'/raindays.dat'

  open(fno+1,file=fname1,access='direct',recl=recl1,form='formatted',status='old',iostat=kode)
  if (kode/=0) then
    write(*,*) 'Climate data file does not exist.'
    write(*,'(''"'',A,''/readme.dat"'')') fname1(1:blank(fname1))
    stop
  endif
  open(fno+2,file=fname2,access='direct',recl=recl1,form='formatted',status='old',iostat=kode)
  if (kode/=0) then
    write(*,*) 'Climate data file does not exist.'
    write(*,'(''"'',A,''/readme.dat"'')') fname2(1:blank(fname2))
    stop
  endif
  open(fno+3,file=fname3,access='direct',recl=recl1,form='formatted',status='old',iostat=kode)
  if (kode/=0) then
    write(*,*) 'Climate data file does not exist.'
    write(*,'(''"'',A,''/readme.dat"'')') fname3(1:blank(fname3))
    stop
  endif
  open(fno+7,file=fname7,access='direct',recl=recl1,form='formatted',status='old',iostat=kode)
  if (kode/=0) then
    write(*,*) 'Climate data file does not exist.'
    write(*,'(''"'',A,''/readme.dat"'')') fname7(1:blank(fname7))
    stop
  endif

  do ii=1,4
    do jj=1,4
      if (indx(ii,jj)==1) then
        if (du==1) then
          do year=yr0,yrf
            read(fno+1,'(12i4)',rec=(siteno(ii,jj)-1)*nyears+year-year0+1) &
 (mtmpvv(ii,jj,year-yr0+1,mnth),mnth=1,12)
            read(fno+2,'(12i4)',rec=(siteno(ii,jj)-1)*nyears+year-year0+1) &
 (mhumvv(ii,jj,year-yr0+1,mnth),mnth=1,12)
            read(fno+3,'(12i4)',rec=(siteno(ii,jj)-1)*nyears+year-year0+1) &
 (mprcvv(ii,jj,year-yr0+1,mnth),mnth=1,12)
            if (kode==0) then
            read(fno+7,'(12i4)',rec=(siteno(ii,jj)-1)*nyears+year-year0+1) &
 (mcldvv(ii,jj,year-yr0+1,mnth),mnth=1,12)
            else
              do mnth=1,12
                mcldvv(ii,jj,year-yr0+1,mnth) = cld_default
              enddo
            endif
          enddo
        else
          do year=yr0,yrf
            read(fno+1,'(12i4)',rec=(siteno(ii,jj)-1)*nyears+year-year0+1) &
 (mtmpvv(ii,jj,year-yr0+1,mnth),mnth=1,12)
            read(fno+2,'(12i4)',rec=(siteno(ii,jj)-1)*nyears+year-year0+1) &
 (mhumvv(ii,jj,year-yr0+1,mnth),mnth=1,12)
            read(fno+3,'(12i4)',rec=(siteno(ii,jj)-1)*nyears+year-year0+1) &
 (mprcvv(ii,jj,year-yr0+1,mnth),mnth=1,12)
            if (kode==0) then
            read(fno+7,'(12i4)',rec=(siteno(ii,jj)-1)*nyears+year-year0+1) &
 (mcldvv(ii,jj,year-yr0+1,mnth),mnth=1,12)
            else
              do mnth=1,12
                mcldvv(ii,jj,year-yr0+1,mnth) = cld_default
              enddo
            endif
          enddo
        endif
      else
        do year=yr0,yrf
          do mnth=1,12
            mtmpvv(ii,jj,year-yr0+1,mnth) = 0
            mhumvv(ii,jj,year-yr0+1,mnth) = 0
            mprcvv(ii,jj,year-yr0+1,mnth) = 0
            mcldvv(ii,jj,year-yr0+1,mnth) = 0
          enddo
        enddo
      endif
    enddo
  enddo

!----------------------------------------------------------------------*
! Interpolate climate.                                                 *
!----------------------------------------------------------------------*
  do year=yr0,yrf
    do mnth=1,12
      do ii=1,4
        do jj=1,4
          xx(ii,jj) = real(mtmpvv(ii,jj,year-yr0+1,mnth))
        enddo
      enddo
      call bi_lin(xx,indx,xnorm,ynorm,ans)
      mtmpv(year-yr0+1,mnth) = round(ans)

      do ii=1,4
        do jj=1,4
          xx(ii,jj) = real(mhumvv(ii,jj,year-yr0+1,mnth))
        enddo
      enddo
      call bi_lin(xx,indx,xnorm,ynorm,ans)
      mhumv(year-yr0+1,mnth) = int(ans+0.5)

      do ii=1,4
        do jj=1,4
          xx(ii,jj) = real(mprcvv(ii,jj,year-yr0+1,mnth))
        enddo
      enddo
      call bi_lin(xx,indx,xnorm,ynorm,ans)
      mprcv(year-yr0+1,mnth) = int(ans+0.5)

      do ii=1,4
        do jj=1,4
          xx(ii,jj) = real(mcldvv(ii,jj,year-yr0+1,mnth))
        enddo
      enddo
      call bi_lin(xx,indx,xnorm,ynorm,ans)
      mcldv(year-yr0+1,mnth) = int(ans+0.5)
    enddo
  enddo

  close(fno+1)
  close(fno+2)
  close(fno+3)
  close(fno+7)

! read WG parameters

! read WG calibration in files for stdtmp, stdhum and raindays
  call READ_WG_FILE(du,fname4,ststats,mstdtmpv,lon,lat,l_stats)
  call READ_WG_FILE(du,fname5,ststats,mstdhumv,lon,lat,l_stats)
  call READ_WG_FILE(du,fname6,ststats,mraindays,lon,lat,l_stats)

! this value should be read in a file. It seems to be calibrated for UK.
 tmp_autocorr_month=0.65
 hum_autocorr_month=0.65

! Call the weather generator
! Generate 30 days a month, as read in the EX_CLIM subroutine

! initialise the yesterday temperature to the mean temperature of the 
! first month, first year
 tmp_yesterday=0.1*mtmpv(1,1)
 hum_yesterday=0.1*mhumv(1,1)

 do year=yr0,yrf
    do mnth=1,12
       if (mprcv(year-yr0+1,mnth)<0.0) then
          if (mnth>1) then
             mprcv(year-yr0+1,mnth)=mprcv(year-yr0+1,mnth-1)
          else
             mprcv(year-yr0+1,mnth)=mprcv(year-yr0+1,12)
          endif
       endif
! generate the rain amount 
       call RAIN_DAY_SUB_FULLY_NORMALISED(mraindays(mnth),mprcv(year-yr0+1,mnth),  &
 30,seed1,seed2,seed3,rain)

! copie into xprcv
       sumprc=0.0
       nbrain=0
       do day=1,30
          sumprc=sumprc+rain(day)
          prcv(year-yr0+1,mnth,day)=min(int(10.0*rain(day)+0.5),30000)
          if (rain(day)>1.0e-3) nbrain=nbrain+1
       enddo
       
! generate the temperaure
       correction=0.0

       do day=1,30
          tmp_yesterday=TEMP_DAY_MEAN2_FN(0.1_dp*FLOATINGMEAN(day,mnth,year-yr0+1, &
 yrf-yr0+1,mtmpv),mstdtmpv(mnth),tmp_autocorr_month,tmp_yesterday, &
 seed1,seed2,seed3,1)
          tmpv(year-yr0+1,mnth,day)=int(100.0*tmp_yesterday+0.5)
          correction=correction+real(tmpv(year-yr0+1,mnth,day))
       enddo

       correction=correction/30.0-(100.0*0.1)*real(mtmpv(year-yr0+1,mnth))
! reajust the mean temperature
       ans = 0.0
       do day=1,30
          tmpv(year-yr0+1,mnth,day)=tmpv(year-yr0+1,mnth,day)-int(correction+0.5)
          ans = ans + tmpv(1,mnth,day)
       enddo

! generate the humidity
       correction=0.0
       do day=1,30
          hum_yesterday = HUM_DAY_MEAN2_FN(0.1_dp*FLOATINGMEAN(day,mnth, &
 year-yr0+1,yrf-yr0+1,mhumv),mstdhumv(mnth),hum_autocorr_month, hum_yesterday, &
 seed1,seed2,seed3,1)
          humv(year-yr0+1,mnth,day) =int(100.0*hum_yesterday + 0.5)
          correction=correction+real(humv(year-yr0+1,mnth,day))
       enddo
       correction=correction/30.0-(100.0*0.1)*real(mhumv(year-yr0+1,mnth))
       do day=1,30
          humv(year-yr0+1,mnth,day)=humv(year-yr0+1,mnth,day)-int(correction+0.5)
       enddo
    enddo

 enddo

else
  l_clim = .false.
endif


end subroutine EX_CLIM_WEATHER_GENERATOR




!----------------------------------------------------------------------*
subroutine READ_WG_FILE(du,filename,ststats,values,lon,lat,l_stats)
!----------------------------------------------------------------------*
character :: filename*1000,st1*1000,st2*1000,st3*1000,ststats*1000
real(dp) :: lon,lat,latr,lonr,values(12),statsv(4,4,12),xx(4,4),latf
real(dp) :: lon0,xnorm,ynorm,ans,rrow,rcol
integer :: recl1,du,mnth,ii,jj,i,indx(4,4),latn,lonn
integer :: recn,row,col,kode
logical :: l_stats

if (du==1) then
!  This wordek for ftn95
!  recl1 = 6*12
  recl1 = 6*12 + 2
else
  recl1 = 6*12+1
endif

open(99,file=ststats(1:blank(ststats))//'/readme.dat',status='old',iostat=kode)
if (kode/=0) then
  write(*,*) 'Climate statistics data file does not exist.'
  write(*,'(''"'',A,''/readme.dat"'')') ststats(1:blank(ststats))
  stop
endif

read (99,'(A)') st1
st2='UNIFORM'
st3='LONLAT'

if (stcmp(st1,st2)==1) then
   close(99)
   open(99,file=filename,status='old')
   read(99,*) (values(mnth),mnth=1,12)
   close(99)
elseif (stcmp(st1,st3)==1) then
  read(99,*)
  read(99,*) latf,lon0
  read(99,*)
  read(99,*) latr,lonr
  read(99,*)
  read(99,*) latn,lonn
  close(99)

!----------------------------------------------------------------------*
! Find the real(dp) row col corresponding to lat and lon.                  *
!----------------------------------------------------------------------*
  rrow = (latf - lat)/latr
  rcol = (lon - lon0)/lonr

  ynorm = rrow - real(int(rrow))
  xnorm = rcol - real(int(rcol))

!----------------------------------------------------------------------*
  open(99,file=filename,status='old',form='formatted',access='direct',recl=recl1)

  do ii=1,4
    do jj=1,4
      row = int(rrow)+jj-1
      col = int(rcol)+ii-1
      if ((row>=1).and.(row<=latn).and.(col>=1).and.(col<=lonn)) then
        recn = (row-1)*lonn + col
        read(99,'(12f6.2)',rec=recn) (statsv(ii,jj,i),i=1,12)
        if (statsv(ii,jj,1)<-0.001) then
          indx(ii,jj) = 0
        else
          indx(ii,jj) = 1
        endif
      else
        indx(ii,jj) = -1
      endif
    enddo
  enddo

  if ((indx(2,2)/=1).and.(indx(2,3)/=1).and.(indx(3,2)/=1) &
 .and.(indx(3,3)/=1)) then
    l_stats = .false.
  endif

  do i=1,12
    do ii=1,4
      do jj=1,4
        xx(ii,jj) = statsv(ii,jj,i)
       enddo
    enddo
    call bi_lin(xx,indx,xnorm,ynorm,ans)
    values(i) = ans
  enddo

  close(99)

else
  write(*,*) 'First line of stats files must be UNIFORM or LONLAT'
  write(*,*)  filename(1:blank(filename))
  stop
endif


end subroutine READ_WG_FILE





!----------------------------------------------------------------------*
real(dp) function FLOATINGMEAN(day,mnth,year,lastyear,array)
!----------------------------------------------------------------------*
integer :: day,mnth,year,lastyear,y,m
integer :: array(500,12)
real(dp) :: mean

if (day<=15) then
  if (mnth>1) then
    y=year
    m=mnth-1
  else
    if (year>1) then
      y=year-1
      m=12
    else
      y=year
      m=1
    endif
  endif
  mean=(day*(1.0*array(year,mnth))+(15-day)*(1.0*array(y,m)))/15.0
else
  if (mnth<12) then
    y=year
    m=mnth+1
  else
    if (year<lastyear) then
      y=year+1
      m=1
    else
      y=year
      m=12
    endif
  endif
  mean=((30-day)*(1.0*array(year,mnth))+(day-15)*(1.0*array(y,m)))/15.0

endif
floatingmean=mean

end function FLOATINGMEAN





!----------------------------------------------------------------------*
subroutine RAIN_DAY_SUB_FULLY_NORMALISED(NO_RAIN_DAYS, &
 RAIN_IN_MONTH,DAYS_IN_MONTH,SEED1,SEED2,SEED3,RAIN)
!----------------------------------------------------------------------*
implicit none

!** The routine simulates daily rainfall for all days in a month
!** The typical number of rainy-days and average rainfall for the month
!** Are the only specific Meteorology inputs. Returned is an array
!** (1-dimensional; day) with the simulated rainfall for each day
!** in it.
logical :: WETDAY(31)
real(dp) :: RAIN(31)
real(dp) :: FRACT_WET, RAIN_PER_DAY, NO_RAIN_DAYS
real(dp) :: ZERO
real(dp) :: PROB_WET_WET, PROB_WET_DRY, GAMMA, MU
real(dp) :: GENE_RAIN_AMT, CORRECTION
integer :: DAYS_IN_MONTH, IDAY, rain_in_month
integer :: SEED1, SEED2, SEED3, IP, IP0, CLOOP
real(dp) :: RAND01 , XLO, XHI, gamma2
real(dp) :: xrain_in_month

xrain_in_month = rain_in_month

if (RAIN_IN_MONTH<1.0e-6_dp) then
   do IDAY=1,DAYS_IN_MONTH
      RAIN(IDAY)=0.0_dp
   enddo
   return
endif
if (NO_RAIN_DAYS<1.0) NO_RAIN_DAYS=1.0
!
! **    determine the proportion of the month with wet-days
FRACT_WET = real(NO_RAIN_DAYS)/real(DAYS_IN_MONTH)
!
!**     find the average rainfall per day with rainfall
RAIN_PER_DAY = RAIN_IN_MONTH/real(NO_RAIN_DAYS)
!
!**   find transitional probability of a wet-day following a dry-day
!**   of the month (Geng et al 1986)
PROB_WET_DRY = 0.75*FRACT_WET
!
!**     probability of wet-day following a wet-day (Geng et al 1986)
PROB_WET_WET = 0.25+PROB_WET_DRY

!
!**     Gamma distribution paramaters (Geng et al 1986)
!       GAMMA = -2.16+1.83*RAIN_PER_DAY  
  XLO = 0.0
  XHI=1.0

  RAND01 = RNDF(XLO,XHI, SEED1, SEED2, SEED3)
  GAMMA = 1.0
!
!**     Mu paramaters (Geng et al 1986)
  gamma2 = gamma+(RAND01-0.5)/2.0
  MU = RAIN_PER_DAY/GAMMA2
!
!     first of the month; the chances of it being wet is 'fract_wet'
!     get a random number proportional to the number of days in the 
!     month

GENE_RAIN_AMT=0.0
XLO = 0.0
XHI=1.0
RAND01 = RNDF(XLO,XHI, SEED1, SEED2, SEED3)
do 20 IDAY = 1, DAYS_IN_MONTH
if(IDAY==1) then
!** is the random number<=rain_fract; if so then it rained yesterday
IP0 = 0
if(RAND01<=FRACT_WET) then
  IP0 = 1
endif
endif
RAIN(IDAY) = 0.0
RAND01 = RNDF(XLO, XHI, SEED1, SEED2, SEED3)
if(IP0==1) then
!**       raining on day 0
if(RAND01-PROB_WET_WET <= 0) IP = 1
if(RAND01-PROB_WET_WET > 0) IP = 0
else   ! IP0 = 0 ; its dry on day 0
if(RAND01-PROB_WET_DRY <= 0) IP = 1
if(RAND01-PROB_WET_DRY > 0) IP = 0
endif
IP0 = IP

  if(IP==1) then     ! Its raining today
     WETDAY(IDAY)=.true.
  else
     WETDAY(IDAY)=.false.
  endif
if(WETDAY(IDAY)) then     ! Its raining today
   CLOOP=0
12       continue
ZERO = 0.0

!        TR1 = EXP(-18.42*MU)

!      TR2 = EXP(-18.42*(1-MU))


! Changed GAMMA to GAMMA2 14/12/01 Paul Henshall.
 
  RAIN(IDAY)=GAMMADIST_FN(zero,MU,GAMMA2,SEED1,SEED2,SEED3)
  CLOOP = CLOOP + 1
  if (CLOOP>10) then
     write(*,*) 'Problem in RAIN_DAY_SUB_FULLY_NORMALISED'
     write(*,*) NO_RAIN_DAYS,RAIN_IN_MONTH, DAYS_IN_MONTH
     stop
  endif
  if (rain(iday)>RAIN_PER_DAY*5) goto 12
endif
  GENE_RAIN_AMT=GENE_RAIN_AMT+RAIN(IDAY)
20    continue

! normalised the rain amount
if (GENE_RAIN_AMT>0.0) then
   CORRECTION = RAIN_IN_MONTH/GENE_RAIN_AMT
     
   do IDAY = 1, DAYS_IN_MONTH
      RAIN(IDAY)=RAIN(IDAY)*CORRECTION
   enddo
else
   IDAY=RNDF(XLO,XHI, SEED1, SEED2, SEED3)*DAYS_IN_MONTH+1
   if (IDAY>DAYS_IN_MONTH) IDAY=DAYS_IN_MONTH
   RAIN(IDAY)=RAIN_IN_MONTH
endif

end subroutine RAIN_DAY_SUB_FULLY_NORMALISED



!
!   this subroutine comes from the RAIN_DAY_SUB subroutine in 
!   ForestETP Weather Generator it differs by the fact that the wet
!   days are firstly generated. Then, knowing the number of rain days
!   in the month REALly generated (with respect to average number of
!   wet days), the rain amount is generated. The goal is to reduce the
!   variability in the rain amount for one year to the other, since the
!   measured monthly mean are known and input of the model.





subroutine RAIN_DAY_SUB_NORMALISED(NO_RAIN_DAYS, RAIN_IN_MONTH, &
 DAYS_IN_MONTH,SEED1,SEED2,SEED3,RAIN)
implicit none
!** The routine simulates daily rainfall for all days in a month
!** The typical number of rainy-days and average rainfall for the month 
!** Are the only specific Meteorology inputs. Returned is an array
!** (1-dimensional; day) with the simulated rainfall for each day 
!** in it.
logical WETDAY(31)
real(dp) RAIN(31)
real(dp) FRACT_WET, RAIN_DAY, RAIN_IN_MONTH, NO_RAIN_DAYS
real(dp) ZERO
real(dp) PROB_WET_WET, PROB_WET_DRY, GAMMA, MU
integer DAYS_IN_MONTH, IDAY, GENE_NO_RAIN_DAYS
integer SEED1, SEED2, SEED3, IP, IP0
real(dp)  RAND01 , XLO, XHI, gamma2

! **    determine the proportion of the month with wet-days
FRACT_WET = real(NO_RAIN_DAYS)/real(DAYS_IN_MONTH)
!
!**     find the average rainfall per day with rainfall
!ccccc calculate later      RAIN_DAY = RAIN_IN_MONTH/real(NO_RAIN_DAYS)
!
!**     find transitional probability of a wet-day following a dry-day
!**     of the month (Geng et al 1986)
PROB_WET_DRY = 0.75*FRACT_WET
!
!**     probability of wet-day following a wet-day (Geng et al 1986)
PROB_WET_WET = 0.25+PROB_WET_DRY
!
!**     Gamma distribution paramaters (Geng et al 1986) 
!       GAMMA = -2.16+1.83*RAIN_DAY   
  XLO = 0.0
  XHI=1.0

  RAND01 = RNDF(XLO,XHI, SEED1, SEED2, SEED3)
  GAMMA = 1.0
!
!**     Mu paramaters (Geng et al 1986) 
  gamma2 = gamma+(rand01-0.5)/2.0
!ccccc calculate later      MU = RAIN_DAY/GAMMA2
!
!     first of the month; the chances of it being wet is 'fract_wet'
!     get a random number proportional to the number of days in the 
!     month

GENE_NO_RAIN_DAYS=0
XLO = 0.0
XHI=1.0
RAND01 = RNDF(XLO,XHI, SEED1, SEED2, SEED3)
do 20 IDAY = 1, DAYS_IN_MONTH
  if(IDAY==1) then
!**   is the random number<=rain_fract; if so then it rained yesterday
    IP0 = 0
    if(RAND01<=FRACT_WET) then
      IP0 = 1 
    endif
  endif
  RAIN(IDAY) = 0.0
  RAND01 = RNDF(XLO, XHI, SEED1, SEED2, SEED3)
  if(IP0==1) then
!**       raining on day 0
    if(RAND01-PROB_WET_WET <= 0) IP = 1
    if(RAND01-PROB_WET_WET > 0) IP = 0
  else   ! IP0 = 0 ; its dry on day 0
    if(RAND01-PROB_WET_DRY <= 0) IP = 1
    if(RAND01-PROB_WET_DRY > 0) IP = 0
  endif
  IP0 = IP

  if(IP==1) then     ! Its raining today
     WETDAY(IDAY)=.true.
     GENE_NO_RAIN_DAYS=GENE_NO_RAIN_DAYS+1
  else
     WETDAY(IDAY)=.false.
  endif
20 continue
!
!     Now, the generated number of wet days in the month is know
!     Let's calculate the rain amount
!
RAIN_DAY = RAIN_IN_MONTH/real(GENE_NO_RAIN_DAYS)
MU = RAIN_DAY/GAMMA2

do 120 IDAY = 1, DAYS_IN_MONTH
  if(WETDAY(IDAY)) then     ! Its raining today
112 continue
!112       TR1 = EXP(-18.42*MU)

!  TR2 = EXP(-18.42*(1-MU))

!         RAND01 = RNDF(XLO, XHI, SEED1, SEED2, SEED3)
!         IF(RAND01 - TR1 .LE. 0) THEN
!           S1 = 0
!         ELSE
!           S1 = RAND01*EXP(1/MU)
!         ENDIF
!         RAND01 = RNDF(XLO, XHI, SEED1, SEED2, SEED3)
!         IF(RAND01 - TR2 .LE. 0) THEN
!           S2 = 0
!         ELSE
!           S2 = RAND01*EXP(1/(1-MU))
!         ENDIF
!         IF(S1+S2-1 .LE. 0) THEN
!           IF( (S1+S2).EQ.0 )THEN
!             Z = 0.0
!           ELSE
!             Z = S1/(S1+S2)
!           ENDIF
!         ELSE
!           FLAG = 1
!         ENDIF
!         IF(FLAG.EQ.1) GOTO  30
!
!         RAND01 = RNDF(XLO, XHI, SEED1, SEED2, SEED3)
!         RAIN(IDAY) = -Z*LOG(RAND01)*GAMMA
	  ZERO = 0.0

! Changed GAMMA to GAMMA2 14/12/01 Paul Henshall.
  
     RAIN(IDAY)=GAMMADIST_FN(zero,MU,GAMMA2,SEED1,SEED2,SEED3)
  if(rain(iday)>rain_day*5)goto 112
  endif
120    continue
!
!//! Equation: 32
!//! Description: The routine simulates daily rainfall for all days in \
!//!  a month. The typical number of rainy-days and average rainfall \
!//!  for the month are the only specific Meteorology inputs. Returned \
!//!  is an array (1-dimensional; day) with the simulated \
!//!  rainfall for each day in it. Although the routine uses  \
!//!  Markov-chain principles, the start of each month is INDEPENDENT \
!//!  of what the conditions were on the last day of the previous  \
!//!  month. \\
!//!  The routine calls a random-number gernerator (RNDF) which \
!//!  returns numbers between Xmin and Xmax (in this case 0,1) using \
!//!  the 3 seeds \\
!//!   INPUTS: Number of rainy days in the month; Average rainfall of\
!//!    the month; number of days in the month; 3 seeds \
!//!    for random number generation\\
!//!   OUTPUT: array (dimension:31), with rainfall for each day
!//! Bibliography: S.Evans (SWELTER)
!
!//! Case: Day_{wet,dry}^{i} = 0
!//!  E: Rain_{day}^{i} = 0.0 \\
!//! \\
!//! Case: Day_{wet,dry}^{i} = 1 \\
!//! \\
!//!  E: Rain_{day}^{i} = Z \gamma log_{e} (Rnd_{j1,j2,j3})  \\
!//!   \gamma = -2.16+1.83 Rain_{day,mean}  \\
!//! \\
!//!  Case: S_{1} + S_{2} -1 <= 0    
!//!   E: Z = { {{S_{1}} \over{S_{1}+S_{2}} }  \\
!//! \\
!//!     Case: Rnd_{j1,j2,j3} - Tr_{1} <= 0
!//!      E: S_{1} = 0
!//!     Case: Rnd_{j1,j2,j3} - Tr_{1} > 0
!//!      E: S_{1} = Rnd_{j1,j2,j3} exp(1/\mu) \\
!//! \\
!//!     Case: Rnd_{j1,j2,j3} - Tr_{2} <= 0
!//!      E: S_{2} = 0
!//!     Case: Rnd_{j1,j2,j3} - Tr_{2} > 0
!//!      E: S_{2} = Rnd_{j1,j2,j3} exp(1/(1-\mu))   \\
!//! \\
!//!  Case: S_{1} + S_{2} -1 > 0
!//!     E: re-generate S_{1} and S_{2}  \\
!//!    \\ 
!//!  Tr_{1} = {exp(-18.42)}\over{\mu} \\
!//!  Tr_{2} = {exp(-18.42)}\over{(1-mu)} \\
!//!  \mu = {{Rain_{day,mean}} \over{\gamma}} \\
!//! \\
!//!  Case: Day_{wet,dry}^{i-1} = 1 \\
!//! \\
!//!   Case: Rnd_{j1,j2,j3} - Prob_{wet,wet} <= 0
!//!    E: Day_{wet,dry}^{i} = 1
!//!   Case: Rnd - Prob_{wet,wet} > 0
!//!    E: Day_{wet,dry}^{i} = 0 \\
!//! \\
!//!  Case: Day_{wet,dry}^{i-1} = 0 \\
!//! \\
!//!   Case: Rnd_{j1,j2,j3} - Prob_{wet,dry} <= 0
!//!    E: Day_{wet,dry}^{i} = 1      
!//!   Case: Rnd_{j1,j2,j3} - Prob_{wet,dry} > 0 
!//!    E: Day_{wet,dry}^{i} = 0      \\
!//! \\
!//!  Prob_{wet,wet} = 0.25+Prob{wet,dry} \\
!//!  Prob_{wet,dry} = { {0.75 N_{rain days}} \over{DinM} }
!
!//! Variable: Rain_{day}^{i}
!//! Description: Predicted daily rainfall on day of month, i
!//! Units: mm
!//! Type: real(dp) array (12,31)
!//! SET
!
!//! Variable: Rain_{day,mean}
!//! Description: Mean daily rainfall for the month
!//! Units: mm/day
!//! Type: REAL
!
!//! Variable: Rnd_{j1,j2,j3)
!//! Description: Random number generated between 0 and 1
!//! Units: none
!//! Type: REAL
!
!//! Variable: j1
!//! Description: Seed for random number generator; changed on return
!//! Units: none
!//! Type: integer
!
!//! Variable: j2
!//! Description: Seed for random number generator; changed on return
!//! Units: none
!//! Type: integer
!
!//! Variable: j3
!//! Description: Seed for random number generator; changed on return
!//! Units: none
!//! Type: integer
!
!//! Variable: DinM
!//! Description: Number of days in the month
!//! Units: none
!//! Type: integer
!
end subroutine RAIN_DAY_SUB_NORMALISED




! from ForestETP weather generator.
! the function which generate the temperature
! has been tranformed into humidity generator
! this is the original function:
!       real(dp) FUNCTION TEMP_DAY_MEAN2_FN(TEMP_MEAN_MONTH,
!     1  TEMP_SD_MONTH, TEMP_AUTOCORR_MONTH, TEMP_YESTERDAY, IS1,
!     2  IS2, IS3, NORMSWITCH)
!
!


real(dp) function HUM_DAY_MEAN2_FN(HUM_MEAN_MONTH, &
 HUM_SD_MONTH, HUM_AUTOCORR_MONTH, HUM_YESTERDAY, IS1, &
 IS2, IS3, NORMSWITCH)

!** calculates mean daily temperature  from monthly values   
real(dp) :: HUM_MEAN_MONTH
real(dp) :: HUM_SD_MONTH, HUM_AUTOCORR_MONTH
real(dp) :: HUM_YESTERDAY
real(dp) :: NORM, NORMMN, NORMSD
integer :: IS1, IS2, IS3, NORMSWITCH
! 
!** approximate a function to get a normal distribution, mean=0, sd=1
NORMMN=0.0
NORMSD = 1.0
NORM=NORM_FROM_RAND_FN(NORMMN, NORMSD, IS1, IS2, IS3, NORMSWITCH)
!
HUM_DAY_MEAN2_FN = HUM_MEAN_MONTH + HUM_AUTOCORR_MONTH* &
 (HUM_YESTERDAY - HUM_MEAN_MONTH) + HUM_SD_MONTH* &
 NORM*( (1-(HUM_AUTOCORR_MONTH)**2)**0.5)
!
!//! Equation: 47
!//! Description: Calculates Average daily humidity from monthly \
!//!  means and standard-deviations. Autocorrelation of 0.65 seems \
!//!  to work if one isn't available. \\
!//!   INPUTS: Mean humidity of the month; Standard deviation \
!//!    around the mean on a daily basis; Humidity \
!//!     auto-correlation for each month; Mean humidity of the \
!//!     previous day; 3 random number seeds.
!//! Bibliography: S Evans (SWELTER); Haith, Tubbs & Pickering (1984)\
!//!  Simulation of Pollution by soil erosion & soil nutirnt loss, \
!//!   PUDOC, Wageningen
!
!//! E: Hum_{mean} = Hum_{mean,month} + \rho_{T} (Hum_{day-1} - \
!//!     Hum_{mean,month}) + \sigma_{mT} N (1-\sigma_{mT}^{2})^0.5) \\
!//!  \\
!//!   E: N = { Random, N(0,1) }
!
!//! Variable: Hum_{mean}
!//! Description: Daily mean humidity
!//! Units: Degrees C   (usually)
!//! Type: Double precsion
!//! SET
!
!//! Variable: Hum_{mean,month}
!//! Description: Mean daily humidity (monthly basis)
!//! Units: Degrees C   (usually)
!//! Type: Double precsion
!
!//! Variable: Hum_{day-1}     
!//! Description: Mean daily humidity of the previous day
!//! Units: Degrees C   (usually)
!//! Type: Double precsion
!
!//! Variable: \rho_{T}
!//! Description: Observed first-order auto-correlation of mean \
!//!  daily air humidity for each month
!//! Units: none  (day/day^{-1})
!//! Type: Double precsion
!
!//! Variable: \sigma_{mT}
!//! Description: Standard deviation of daily air humidity about \
!//!  the mean (monthly basis).
!//! Units: Degrees C             
!//! Type: Double precsion
!
!//! Variable: R
!//! Description: Random number between 0,1
!//! Units: none
!//! Type: Double precsion
!
end function HUM_DAY_MEAN2_FN

end module weather_generator

