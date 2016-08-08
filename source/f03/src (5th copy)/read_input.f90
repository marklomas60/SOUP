module read_input

use real_precision
use pft_parameters
use site_parameters
use sdgvm1
use func
use misc_parameters
use tuning_parameters
use open_files

implicit none

contains

!**********************************************************************!
!                                                                      !
!                    read_input_file :: read_input                     !
!                    -----------------------------                     !
!                                                                      !
! SUBROUTINE read_input_file(buff1,stinput,stinit,stoutput,sttxdp,     !
! mask,stlu,ststats,stco2,xlatf,xlon0,xlatres,xlonres,co2const,initise,!
! initiseo,speedc, xspeedc,xseed1,spinl,yr0s,cycle,crand,yr0p,yrfp,    !
! outyears,nyears,yr0,yrf,yearind,idum,yearv,nomdos,otags,omav,ofmt,   !
! outyears1,outyears2,oymd,otagsn,otagsnft,snpshts,snp_no,out_cov,     !
! out_bio,out_bud,out_sen,lutab,grassrc,barerc,fireres,luse,l_b_and_c, !
! soil_chr,topsl,defaulttopsl,sites,latdel,londel,lat_lon,day_mnth,    !
! thty_dys,xlatresn,xlonresn,ilanduse,nft,xyear0,xyearf,lmor_sc,       !
! oymdft,iofnft,sit_grd,du,narg)                                       !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine read_input_file(buff1,stinput,stinit,stoutput,sttxdp,stmask,&
 stlu,ststats,stco2,xlatf,xlon0,xlatres,xlonres,co2const,initise,&
 initiseo,speedc, xspeedc,xseed1,spinl,yr0s,cycle,crand,yr0p,yrfp,&
 outyears,nyears,yr0,yrf,yearind,idum,yearv,nomdos,otags,omav,ofmt,&
 outyears1,outyears2,oymd,otagsn,otagsnft,snpshts,snp_no,out_cov,&
 out_bio,out_bud,out_sen,lutab,grassrc,barerc,fireres,luse,l_b_and_c,&
 soil_chr,topsl,defaulttopsl,sites,latdel,londel,lat_lon,day_mnth,&
 thty_dys,xlatresn,xlonresn,ilanduse,nft,xyear0,xyearf,lmor_sc,&
 oymdft,iofnft,sit_grd,du,narg)
!**********************************************************************!
integer :: fireres,sites,per,site,ibox,jbox,l,recl1, &
 site0,sitef,luse(max_years),ilanduse,persum
real(dp) :: latdel,londel,lat_lon(2000000,2),iadj,jadj, &
 grassrc,barerc,topsl,defaulttopsl,soil_chr(10),lutab(255,100),&
 lmor_sc(3600,max_cohorts)
character :: st1*1000,st2*1000,st3*1000,st4*1000
character :: st5*1000,otags(100)*1000,ofmt(100)*100
character :: stinput*1000,stoutput*1000,stinit*1000,stco2*1000
character :: stmask*1000
character :: sttxdp*1000,stlu*1000,ststats*1000,buff1*80,fttags(100)*1000
integer :: i,du,xyear0,xyearf,xlatresn,xlonresn,kode,ii, &
 xseed1,spinl,yr0s,cycle,yr0p,yrfp,outyears,yr0,yrf,idum, &
 nyears,yearv(max_years),yearind(max_years),nomdos,omav(50),outyears1,&
 oymd,otagsn(50),otagsnft(50),outyears2,snp_no,snpshts(1000),&
 sit_grd,day_mnth,thty_dys,narg,j,iofnft,iofn,ft,nft,k,oymdft
real(dp) :: xlatf,xlon0,lon0,lonf,lat0,latf,xlatres,xlonres,co2const
logical :: initise,initiseo,speedc,xspeedc,crand
logical :: out_cov,out_bio,out_bud,out_sen,l_b_and_c
!----------------------------------------------------------------------!

open(98,file=buff1,status='OLD',iostat=kode)
if (kode/=0) then
  do i=1,80
    st1(i:i) = buff1(i:i)
  enddo
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'SDGVM input file does not exist.'
  write(*,'('' "'',A,''"'')') trim(buff1)
  stop
endif

read(98,'(A)') st1
st2 = 'DOS'
st3 = 'UNIX'

if (stcmp(st1,st2)==1) then
  du = 1
elseif (stcmp(st1,st3)==1) then
  du = 0
else
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'First line of input file must read DOS or UNIX.'
  stop
endif

read(98,'(A)') stinput
call STRIPB(stinput)
st1 = stinput
open(99,file=st1(1:blank(st1))//'/readme.dat',status='OLD',iostat=kode)
if (kode/=0) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'Climate data file does not exist:'
  write(*,'('' "'',A,''/readme.dat"'')') st1(1:blank(st1))
  stop
endif
read(99,'(A)') st1

!----------------------------------------------------------------------!
! Determine whether climate is daily or monthly, from first line of    !
! readme file.                                                         !
!----------------------------------------------------------------------!
st2 = 'DAILY'
st3 = 'MONTHLY'
st4 = 'SITED'
st5 = 'SITEM'
sit_grd = 0
if (stcmp(st1,st2)==1) then
  day_mnth = 1
  thty_dys = 1
  read(99,*)
  read(99,*) xlatf,xlon0
  read(99,*)
  read(99,*) xlatres,xlonres
  read(99,*)
  read(99,*) xlatresn,xlonresn
  read(99,*)
  read(99,*) xyear0,xyearf
elseif (stcmp(st1,st3)==1) then
  day_mnth = 0
  thty_dys = 1
  read(99,*)
  read(99,*) xlatf,xlon0
  read(99,*)
  read(99,*) xlatres,xlonres
  read(99,*)
  read(99,*) xlatresn,xlonresn
  read(99,*)
  read(99,*) xyear0,xyearf
elseif (stcmp(st1,st4)==1) then
  day_mnth = 1
  thty_dys = 0
  sit_grd = 1
  read(99,*)
  read(99,*) xlatf,xlon0
  read(99,*)
  read(99,*) xyear0,xyearf
elseif (stcmp(st1,st5)==1) then
  day_mnth = 0
  thty_dys = 1
  sit_grd = 1
  read(99,*)
  read(99,*) xlatf,xlon0
  read(99,*)
  read(99,*) xyear0,xyearf
else
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) &
 'First line of climat readme.dat file must read DAILY MONTHLY or SITE.'
  stop
endif
close(99)
!----------------------------------------------------------------------!

st1 = stinput

read(98,'(A)') ststats
call STRIPB(ststats)

read(98,'(A)') sttxdp
call STRIPB(sttxdp)

read(98,'(A)') stlu
call STRIPB(stlu)

read(98,'(A)') stco2
call STRIPB(stco2)

read(98,*) co2const

read(98,'(A)') stmask
call STRIPB(stmask)

read(98,*)

narg = 1
read(98,'(A)') stinit
call STRIPB(stinit)
st2 = 'ARGUMENT'
if (stcmp(stinit,st2)==1) then
  narg = narg + 1
  call GETARG(narg,buff1)
  do i=1,80
    stinit(i:i) = buff1(i:i)
  enddo
endif

read(98,'(A)') stoutput
call STRIPB(stoutput)
st2 = 'ARGUMENT'
if (stcmp(stoutput,st2)==1) then
  narg = narg + 1
  call GETARG(narg,buff1)
  do i=1,80
    stoutput(i:i) = buff1(i:i)
  enddo
endif
read(98,*)

!----------------------------------------------------------------------!
! Check if output directory exists.                                    !
!----------------------------------------------------------------------!
st1 = stoutput
open(11,file=st1(1:blank(st1))//'/diag.dat',iostat=kode)
if (kode/=0) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'SDGVM output directory does not exist.'
  write(*,'('' "'',A,''"'')') st1(1:blank(st1))
  stop
endif
write(11,'(''ERRORS'')')
close(11)
!----------------------------------------------------------------------!

read(98,'(1000a)') st1
ii = n_fields(st1)
if ((ii==3).or.(ii==4)) then
  call STRIPBN(st1,i)
  initise = .true.
  if (i==0)  initise = .false.
  call STRIPBN(st1,j)
  initiseo = .true.
  if (j==0)  initiseo = .false.
  call STRIPBN(st1,k)
  speedc = .true.
  if (k==0)  speedc = .false.
  if (ii==4) then
    call STRIPBN(st1,xseed1)
  else
    xseed1 = 1
  endif
else
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'Line 13 must contain 3 or 4 arguments'
  write(*,'('' "'',A,''"'')') st1(1:30)
  stop
endif
xspeedc = speedc

read(98,'(1000a)') st1
call STRIPB(st1)
i = n_fields(st1)
if (i==7) then
  call STRIPBN(st1,spinl)
  call STRIPBN(st1,yr0s)
  call STRIPBN(st1,cycle)
  call STRIPBN(st1,j)
  crand = .true.
  if (j==0)  crand = .false.
  call STRIPBN(st1,yr0p)
  call STRIPBN(st1,yrfp)
  call STRIPBN(st1,outyears)
elseif (i==6) then
  call STRIPBN(st1,spinl)
  call STRIPBN(st1,yr0s)
  call STRIPBN(st1,cycle)
  call STRIPBN(st1,j)
  crand = .true.
  if (j==0)  crand = .false.
  call STRIPBN(st1,yr0p)
  call STRIPBN(st1,yrfp)
  outyears = yrfp - yr0p + 1
elseif (i==5) then
  call STRIPBN(st1,spinl)
  call STRIPBN(st1,yr0s)
  call STRIPBN(st1,cycle)
  call STRIPBN(st1,j)
  crand = .true.
  if (j==0)  crand = .false.
  yr0p = yr0s+1
  yrfp = yr0s
  call STRIPBN(st1,outyears)
elseif (i==4) then
  call STRIPBN(st1,spinl)
  call STRIPBN(st1,yr0s)
  call STRIPBN(st1,cycle)
  call STRIPBN(st1,j)
  crand = .true.
  if (j==0)  crand = .false.
  yr0p = yr0s+1
  yrfp = yr0s
  outyears = cycle + 1
elseif (i==3) then
  spinl = 0
  cycle = max_years
  j = 0
  crand = .false.
  call STRIPBN(st1,yr0p)
  call STRIPBN(st1,yrfp)
  call STRIPBN(st1,outyears)
  yr0s = yr0p
elseif (i==2) then
  spinl = 0
  cycle = max_years
  j = 0
  crand = .false.
  call STRIPBN(st1,yr0p)
  call STRIPBN(st1,yrfp)
  yr0s = yr0p
  outyears = yrfp - yr0p + 1
else
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'Line 14 must contain 2-7 fields'
  stop
endif
nyears = spinl+yrfp-yr0p+1
if (nyears>max_years) then
  write(*,'('' PROGRAM TERMINATED'')')
 write(*,'('' Trying to simulate '',i4,&
 &'' years, maximum allowable is '',i4,''.'')') nyears,max_years
 write(*,*) &
 'Either reduce the length of the simulation, or increase "max_years"'
  write(*,*) &
 'this is set in array_param.txt, you must re-comile after altering'
  write(*,*) 'this file.'
  stop
endif
if ((j/=0).and.(j/=1)) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,'('' Fourth field of line 14 in the input file must be &
 &either 0 or 1.'')')
  write(*,'('' Currently set to '',i5,''.'')') j
  stop
endif
outyears = min(outyears,nyears)

!----------------------------------------------------------------------!
! Set 'yr0' and 'yrf' which are the years of actual climate required   *
! for the run. And check that the climate exists in the climate        !
! database                                                             !
!----------------------------------------------------------------------!
yr0 = min(yr0s,yr0p)
yrf = max(yr0s+min(spinl,cycle)-1,yrfp)
if ((yr0<xyear0).or.(yrf>xyearf)) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,'('' Tyring to use '',i4,''-'',i4,'' climate.'')') yr0,yrf
  write(*,'('' Climate database runs from '',i4,''-'',i4,''.'')') &
 xyear0,xyearf
  stop
endif

!----------------------------------------------------------------------!
! Set up year climate sequence.                                        !
!----------------------------------------------------------------------!
do i=1,max_years
  yearind(i) = i
enddo
idum = 1
do i=1,nyears
  if (i<=spinl) then
    if (crand) then
      if ((mod(i-1,cycle)==0).and.(i>2)) &
 call RANDOMV(yearind,1,cycle,idum)
      yearv(i) = yearind(mod(i-1,cycle)+1) + yr0s - 1
    else
      yearv(i) = mod(i-1,cycle) + yr0s
    endif
  else
    yearv(i) = i - spinl + yr0p - 1
  endif
enddo

!----------------------------------------------------------------------!
! Open monthly/daily PIXEL output files if required.                   !
!----------------------------------------------------------------------!
call OUTPUT_OPTIONS(nomdos,otags,omav,ofmt)

read(98,'(A)') st1
call STRIPB(st1)
st2 = 'PIXEL'
if (stcmp(st1,st2)==0) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) &
 'The first field of line 15 of the input file must read ''PIXEL''.'
  write(*,*) st1(1:30)
  write(*,*) 'Output variable options:'
  write(*,'(1x,20a4)') (otags(i),i=1,15)
  write(*,'(1x,20a4)') (otags(i),i=16,nomdos)
  stop
endif

call SET_PIXEL_OUT(st1,st2,outyears1,nomdos,otagsn,otags,oymd)
outyears1 = min(outyears1,nyears)

st1 = stoutput
iofn = 400
if (outyears1>0) then
if ((oymd==1).or.(oymd==3)) then
  do i=1,50
    if (otagsn(i)==1) then
      iofn = iofn + 1
      open(iofn,&
      file=st1(1:blank(st1))//'/monthly_'//otags(i)(1:3)//'.dat')
    endif
  enddo
endif
if ((oymd==2).or.(oymd==3)) then
  do i=1,50
    if (otagsn(i)==1) then
      iofn = iofn + 1
      open(iofn,&
 file=st1(1:blank(st1))//'/daily_'//otags(i)(1:3)//'.dat')
    endif
  enddo
endif
endif
iofnft = iofn

!----------------------------------------------------------------------!
! Determine whether daily or monthly subpixel outputs are required.    !
!----------------------------------------------------------------------!

read(98,'(A)') st1
call STRIPB(st1)
st2 = 'SUBPIXEL'
if (stcmp(st1,st2)==0) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) &
 'The first field of line 16 of the input file must read ''SUBPIXEL''.'
  write(*,'('' "'',A,''"'')') st1(1:30)
  write(*,*) 'Output variable options:'
  write(*,'(1x,20a4)') (otags(i),i=1,15)
  write(*,'(1x,20a4)') &
 (otags(i),i=16,nomdos),'cov ','bio ','bud ','sen '
  stop
endif

call SET_SUBPIXEL_OUT(st1,st2,outyears2,nomdos,otagsnft,otags,oymdft,&
 out_cov,out_bio,out_bud,out_sen)
outyears2 = min(outyears2,nyears)

!----------------------------------------------------------------------!
! Read in snapshot years.                                              !
!----------------------------------------------------------------------!
read(98,'(A)') st1
call STRIPB(st1)
st2 = 'SNAPSHOTS'
if (stcmp(st1,st2)==0) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'First field of line 17 must read SNAPSHOTS'
  stop
endif
call STRIPBS(st1,st2)
call STRIPB(st1)
call ST2ARR(st1,snpshts,100,snp_no)
read (98,*)

!----------------------------------------------------------------------!
! Read in compulsory functional types.                                 !
!----------------------------------------------------------------------!
pft_tab(1)%tag = 'BARE'
pft_tab(2)%tag = 'CITY'
pft_tab(1)%itag = 1
pft_tab(2)%itag = 2
read(98,'(A)') st1
call STRIPB(st1)
st2 = 'BARE'
if ((stcmp(st2,st1)==0).or.(n_fields(st1)/=2)) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'Line 19 must read "BARE" forllowed by a number (0-1).'
  stop
endif
call STRIPBS(st1,st2)
read(st1,*) pft_tab(1)%mix
read(98,'(A)') st1
call STRIPB(st1)
st2 = 'CITY'
if ((stcmp(st2,st1)==0).or.(n_fields(st1)/=2)) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'Line 19 must read "BARE" forllowed by a number (0-1).'
  stop
endif
call STRIPBS(st1,st2)
read(st1,*) pft_tab(2)%mix
read (98,*)

! Initialise redundant parameterisation
ft = 1
pft_tab(ft)%c3c4 = 0
pft_tab(ft)%phen = 0
pft_tab(ft)%crop = 0.0
pft_tab(ft)%d2h = 0
pft_tab(ft)%mort = 1
pft_tab(ft)%wden = 0.0
pft_tab(ft)%xyl = 0.0
pft_tab(ft)%pdif = 0.0
pft_tab(ft)%sla = 0.0
pft_tab(ft)%lls = 0
pft_tab(ft)%sls = 0
pft_tab(ft)%rls = 0
pft_tab(ft)%lmor = 0.0
pft_tab(ft)%lrat = 0.0
pft_tab(ft)%bbmem = 0
pft_tab(ft)%bb0 = 0.0
pft_tab(ft)%bbmax = 0.0
pft_tab(ft)%bblim = 0.0
pft_tab(ft)%senm = 0
pft_tab(ft)%sens = 0
pft_tab(ft)%senlim = 0.0
pft_tab(ft)%stemx = 0.0
pft_tab(ft)%gr0 = 0.0
pft_tab(ft)%grf = 0.0
pft_tab(ft)%ppm0 = 0.0

ft = 2
pft_tab(ft)%c3c4 = 0
pft_tab(ft)%phen = 0
pft_tab(ft)%crop = 0.0
pft_tab(ft)%d2h = 0
pft_tab(ft)%mort = 1
pft_tab(ft)%wden = 0.0
pft_tab(ft)%xyl = 0.0
pft_tab(ft)%pdif = 0.0
pft_tab(ft)%sla = 0.0
pft_tab(ft)%lls = 0
pft_tab(ft)%sls = 0
pft_tab(ft)%rls = 0
pft_tab(ft)%lmor = 0.0
pft_tab(ft)%lrat = 0.0
pft_tab(ft)%bbmem = 0
pft_tab(ft)%bb0 = 0.0
pft_tab(ft)%bbmax = 0.0
pft_tab(ft)%bblim = 0.0
pft_tab(ft)%senm = 0
pft_tab(ft)%sens = 0
pft_tab(ft)%senlim = 0.0
pft_tab(ft)%stemx = 0.0
pft_tab(ft)%gr0 = 0.0
pft_tab(ft)%grf = 0.0
pft_tab(ft)%ppm0 = 0.0

!----------------------------------------------------------------------!
! Read in functional type parameterisation.                            !
!----------------------------------------------------------------------!

read(98,'(A)') st1
call STRIPB(st1)

if (n_fields(st1)==1) then
  st2 = 'ARGUMENT'

  if (stcmp(st2,st1)==1) then
!----------------------------------------------------------------------!
! From a file.                                                         !
!----------------------------------------------------------------------!
    narg = narg + 1
    call GETARG(narg,buff1)
    do i=1,80
      st1(i:i) = buff1(i:i)
    enddo
  endif

  open(99,file=st1,status='OLD',iostat=kode)
  if (kode/=0) then
    write(*,'('' PROGRAM TERMINATED'')')
    write(*,*) 'Functional type parameterisation file does not exist.'
    write(*,'('' "'',A,''"'')') st1(1:blank(st1))
    stop
  endif

  ft = 2

40      continue
    read(99,'(1000a)',end=30) st1
    ft = ft + 1

      if (ft>max_cohorts) then
        write(*,'('' PROGRAM TERMINATED'')')
        write(*,'('' Maximum number of fts is'',i4,'', &
 &currently more than this are defined.'')') max_cohorts
        write(*,*) &
 'Either increase "max_cohorts" defined in "dims.f90" or decrease'
        write(*,*) &
 'the number of ft parameterisations: this requires re-compilation.'
        stop
      endif

    if (n_fields(st1)/=27) then
      write(*,'('' PROGRAM TERMINATED'')')
      write(*,*) 'The ft parameterisation must contain 27 fields.'
      write(*,'(1x,A,'' has '',i3)') st1(1:blank(st1)),n_fields(st1)
      stop
    endif

    pft_tab(ft)%itag=ft
    read(st1,*) pft_tab(ft)%tag,pft_tab(ft)%mix,pft_tab(ft)%c3c4,&
 pft_tab(ft)%phen, pft_tab(ft)%crop,pft_tab(ft)%d2h,&
 pft_tab(ft)%mort,pft_tab(ft)%wden,pft_tab(ft)%xyl,pft_tab(ft)%pdif,&
 pft_tab(ft)%sla,pft_tab(ft)%lls,pft_tab(ft)%sls,pft_tab(ft)%rls,&
 pft_tab(ft)%lmor,pft_tab(ft)%lrat,pft_tab(ft)%bbmem,&
 pft_tab(ft)%bb0,pft_tab(ft)%bbmax,pft_tab(ft)%bblim,&
 pft_tab(ft)%senm,pft_tab(ft)%sens,pft_tab(ft)%senlim,&
 pft_tab(ft)%stemx,pft_tab(ft)%gr0,pft_tab(ft)%grf,pft_tab(ft)%ppm0

    if (pft_tab(ft)%mort>max_years) then
      write(*,'('' PROGRAM TERMINATED'')')
      write(*,'('' Mortality of '',A,'' is '',i4,'', maximum allowable &
 &is '',i4,''.'')')
      stop
    endif

    if (pft_tab(ft)%sla<0.0) then
      pft_tab(ft)%sla = &
 10.0**(2.35 - 0.39*log10(real(pft_tab(ft)%lls)/30.0))*2.0/10000.0
!      pft_tab(ft)%sla = 10.0**(2.43-0.46*log10(real(pft_tab(ft)%lls)/30.0))
!     &*2.0/10000.0
    endif

    pft_tab(ft)%sla = pft_tab(ft)%sla/tgp%p_sla

    if (pft_tab(ft)%lls<0.0) then
      pft_tab(ft)%lls = &
 int(10.0**((2.35 - log10(pft_tab(ft)%sla*10000.0/2.0))/0.39)*30.0+0.5)
    endif

    if (pft_tab(ft)%mort>max_age) then
      write(*,'('' PROGRAM TERMINATED'')')
      write(*,'('' Maximum age of "'',A,''" is '',i4,'', maximum allowable age is '',i4,''.'')') &
 st1(1:blank(st1)),pft_tab(ft)%mort,max_age
      write(*,*) 'Either increase "max_age" defined in "dims.f90"'
      write(*,*) 'or decrease mortality in the parameterisation.'
      stop
    endif

  goto 40
30 continue
  nft = ft
  read (98,*)

else
!----------------------------------------------------------------------!
! From the input file.                                                 !
!----------------------------------------------------------------------!
  ft = 2
98      continue

    if (ichar(st1(1:1))/=32) then

      ft = ft + 1

      if (ft>max_cohorts) then
        write(*,'('' PROGRAM TERMINATED'')')
        write(*,'('' Maximum number of fts is'',i4,'', currently more &
 &than this are defined.'')') max_cohorts
        write(*,*) 'Either increase "max_cohorts" defined &
 &in "dims.f90" or decrease'
        write(*,*) 'the number of ft parameterisations: this &
 &requires re-compilation.'
        stop
      endif

      if (n_fields(st1)/=27) then
        write(*,'('' PROGRAM TERMINATED'')')
        write(*,*) 'The ft parameterisation must contain 27 fields.'
        write(*,'(1x,A,'' has '',i3)') st1(1:blank(st1)),n_fields(st1)
        stop
      endif

      pft_tab(ft)%itag=ft
      read(st1,*) pft_tab(ft)%tag,pft_tab(ft)%mix,pft_tab(ft)%c3c4,&
 pft_tab(ft)%phen,pft_tab(ft)%crop,pft_tab(ft)%d2h,pft_tab(ft)%mort,&
 pft_tab(ft)%wden,pft_tab(ft)%xyl,pft_tab(ft)%pdif,pft_tab(ft)%sla,&
 pft_tab(ft)%lls,pft_tab(ft)%sls,pft_tab(ft)%rls,pft_tab(ft)%lmor,&
 pft_tab(ft)%lrat,pft_tab(ft)%bbmem,pft_tab(ft)%bb0,pft_tab(ft)%bbmax,&
 pft_tab(ft)%bblim,pft_tab(ft)%senm,pft_tab(ft)%sens,pft_tab(ft)%senlim,&
 pft_tab(ft)%stemx,pft_tab(ft)%gr0,pft_tab(ft)%grf,pft_tab(ft)%ppm0

      if (pft_tab(ft)%sla<0.0) then
        pft_tab(ft)%sla = 10.0**(2.35 - &
 0.39*log10(real(pft_tab(ft)%lls)/30.0))*2.0/10000.0
!      pft_tab(ft)%sla = 10.0**(2.43-&
! 0.46*log10(real(pft_tab(ft)%lls)/30.0))*2.0/10000.0
      endif

      pft_tab(ft)%sla = pft_tab(ft)%sla/tgp%p_sla

      if (pft_tab(ft)%lls<0.0) then
        pft_tab(ft)%lls = int(10.0**((2.35 - log10(pft_tab(ft)%sla*&
 10000.0/2.0))/0.39)*30.0+0.5)
      endif

      if (pft_tab(ft)%mort>max_age) then
        write(*,'('' PROGRAM TERMINATED'')')
        write(*,'('' Maximum age of "'',A,''" is '',i4,'', maximum &
 &allowable age is '',i4,''.'')') &
 st1(1:blank(st1)),pft_tab(ft)%mort,max_age
        write(*,*) 'Either increase "max_age" defined in "dims.f90"'
        write(*,*) 'or decrease mortality in the parameterisation.'
        stop
      endif

      read(98,'(1000a)') st1
      call STRIPB(st1)
    goto 98
  endif
  nft = ft
endif

!----------------------------------------------------------------------!
! Change the units of xylem and water potential difference.            !
!----------------------------------------------------------------------!
do ft=1,nft
  pft_tab(ft)%xyl =  pft_tab(ft)%xyl*1.0e-9
  pft_tab(ft)%pdif  =  pft_tab(ft)%pdif*1.0e+3
enddo

!----------------------------------------------------------------------!
! Create leaf mortality scales values.                                 !
!----------------------------------------------------------------------!
lmor_sc(1,1) = 0.0_dp
lmor_sc(2,1) = 0.0_dp
do ft=3,nft
  do i=1,pft_tab(ft)%lls
    lmor_sc(i,ft)=(real(pft_tab(ft)%lls-i)/&
 real(pft_tab(ft)%lls))**pft_tab(ft)%lmor
  enddo
enddo

!----------------------------------------------------------------------!
! Open site_info file and check output directory exists.               !
!----------------------------------------------------------------------!
st1 = stoutput
open(12,file=st1(1:blank(st1))//'/site_info.dat',iostat=kode)
if (kode/=0) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'Output directory does not exist:'
  write(*,'('' "'',A,''"'')') st1(1:blank(st1))
  stop
endif
write(12, &
'(''COUNTRY           ID     lat      lon     co20  co2f    sand   &
 &silt   bulk     wp     fc    swc    dep       '' &
& ,50A13)') (pft_tab(ft)%tag,ft=1,nft)

!----------------------------------------------------------------------!
! Open monthly/daily ft output files.                                  !
!----------------------------------------------------------------------!
if (outyears2>0) then
  iofn = iofnft
  st1 = stoutput
  if ((oymdft==1).or.(oymdft==3)) then
    do i=1,50
      if (otagsnft(i)==1) then
        do ft=1,nft
          iofn = iofn + 1
          open(iofn,file=st1(1:blank(st1))//'/monthly_'//otags(i) &
 (1:3)//'_'//pft_tab(ft)%tag(1:blank(pft_tab(ft)%tag))//'.dat')
       enddo
      endif
    enddo
  endif
  if ((oymdft==2).or.(oymdft==3)) then
    do i=1,50
      if (otagsnft(i)==1) then
        do ft=1,nft
          iofn = iofn + 1
          open(iofn,file=st1(1:blank(st1))//'/daily_'//otags(i) &
 (1:3)//'_'//pft_tab(ft)%tag(1:blank(pft_tab(ft)%tag))//'.dat')
        enddo
      endif
    enddo
  endif
endif

!----------------------------------------------------------------------!
! Read in landuse index mapping.                                       !
!----------------------------------------------------------------------!
do j=1,255
  do k=1,nft
    lutab(j,k) = 0.0
  enddo
enddo

do ft=1,nft
  fttags(ft) = pft_tab(ft)%tag
enddo

ii = 0
96 continue
  read(98,'(1000a)') st1
  ii = ii + 1
  call STRIPB(st1)
! Read table of conversion from class to ft's proportion
  if (ichar(st1(1:1))/=32) then
    i = n_fields(st1)
    call STRIPBN(st1,j)
    persum = 0
    do k=1,(i-1)/2
      call STRIPBN(st1,per)
      call STRIPBS(st1,st2)
      l = ntags(fttags,st2)
      if (l > 1)  persum = persum + per
      if (l==-9999) then
        write(*,'('' PROGRAM TERMINATED'')')
        write(*,*) 'Error in a tag name in the land use mapping.'
        write(*,&
 '('' Line number'',i3,'', within the mapping list.'')') ii
        write(*,&
 '('' Category number'',i3,'', tag number'',i2,''.'')') j,k
        write(*,'('' "'',A,''"'')') st2(1:blank(st2))
        stop
      endif
      lutab(j,l) = real(per)
    enddo
    lutab(j,1) = 100.0 - real(persum)

  goto 96
endif

!----------------------------------------------------------------------!

read(98,*) grassrc, barerc, fireres
read(98,*)

!----------------------------------------------------------------------!
! Read in sites.                                                       !
!----------------------------------------------------------------------!
latdel = 0.0
londel = 0.0

if ((day_mnth==1).and.(thty_dys==0)) then
!----------------------------------------------------------------------!
! Single site defined in climate file.                                 !
!----------------------------------------------------------------------!
  sites = 1
  lat_lon(sites,1) = xlatf
  lat_lon(sites,2) = xlon0
  read(98,*)
  read(98,*)
else

  read(98,'(1000a)') st1
  call STRIPB(st1)

  if (n_fields(st1)==1) then

    if (du==1) then
      recl1 = 16
    else
      recl1 = 17
    endif

    st2 = 'ARGUMENT'
    if (stcmp(st2,st1)==1) then

      narg = narg + 1
      call GETARG(narg,buff1)
      do i=1,80
        st1(i:i) = buff1(i:i)
      enddo

      narg = narg + 1
      call GETARG(narg,buff1)
      do i=1,80
        st2(i:i) = buff1(i:i)
      enddo
      call STRIPBN(st2,site0)

      narg = narg + 1
      call GETARG(narg,buff1)
      do i=1,80
        st2(i:i) = buff1(i:i)
      enddo
      call STRIPBN(st2,sitef)

      open(99,file=st1,status='OLD',form='FORMATTED',&
 access='DIRECT',recl=recl1,iostat=kode)
      if (kode/=0) then
        write(*,'('' PROGRAM TERMINATED'')')
        write(*,*) 'File containing list of sites does not exist.'
        write(*,'('' "'',A,''"'')') st1(1:blank(st1))
        stop
      endif

      do site=site0,sitef
        read(99,'(F7.3,F9.3)',rec=site) &
 lat_lon(site-site0+1,1),lat_lon(site-site0+1,2)
      enddo
      sites = sitef - site0 + 1

! CLOSE added by Ghislain 15/12/03
      close(99)
    else

      open(99,file=st1,status='OLD',iostat=kode)
      if (kode/=0) then
        write(*,'('' PROGRAM TERMINATED'')')
        write(*,*) 'File containing list of sites does not exist.'
        write(*,'('' "'',A,''"'')') st1(1:blank(st1))
        stop
      endif

      sites = 0
50          read(99,*,end=60) lat_lon(sites+1,1),lat_lon(sites+1,2)
        sites = sites + 1
      goto 50
60          continue

! CLOSE added by Ghislain 15/12/03
      close(99)
    endif

    read(98,*)

  elseif ((n_fields(st1)==2).or.(n_fields(st1)==3)) then

!----------------------------------------------------------------------!
! List of sites in the input file.                                     !
!----------------------------------------------------------------------!
    sites = 1
    read(st1,*) lat_lon(sites,1),lat_lon(sites,2)
97        continue
      read(98,'(1000a)') st1
      call STRIPB(st1)
      if (ichar(st1(1:1))/=32) then
        sites = sites + 1
        read(st1,*) lat_lon(sites,1),lat_lon(sites,2)
      goto 97
    endif

  elseif ((n_fields(st1)==4).or.(n_fields(st1)==5)) then

!----------------------------------------------------------------------!
! Box of sites in input file.                                          !
!----------------------------------------------------------------------!
    read(st1,*) lat0,latf,lon0,lonf
    if ((lat0>latf).or.(lon0>lonf)) then
      write(*,'('' PROGRAM TERMINATED'')')
      write(*,*) 'lat0 and lon0 must be less than latf and lonf.'
    endif
    if (lat0<-90.0)  lat0 = -90.0
    if (latf> 90.0)  latf =  90.0
    if (lon0<-180.0)  lon0 = -180.0
    if (lonf> 180.0)  lonf =  180.0
    sites = 0
    do i=1,int((latf-lat0)/xlatres+0.5)
      do j=1,int((lonf-lon0)/xlonres+0.5)
        sites = sites + 1
        lat_lon(sites,1) = lat0 + (real(i-1)+0.5)*xlatres
        lat_lon(sites,2) = lon0 + (real(j-1)+0.5)*xlonres
      enddo
    enddo
    read(98,*)
    latdel = xlatres
    londel = xlonres

  elseif (n_fields(st1)>=6) then

!----------------------------------------------------------------------!
! Box of sites in input file.                                          !
!----------------------------------------------------------------------!
    read(st1,*) lat0,latf,lon0,lonf,latdel,londel
    if ((lat0>latf).or.(lon0>lonf)) then
      write(*,'('' PROGRAM TERMINATED'')')
      write(*,*) 'lat0 and lon0 must be less than latf and lonf.'
    endif
    if (lat0<-90.0)  lat0 = -90.0
    if (latf> 90.0)  latf =  90.0
    if (lon0<-180.0)  lon0 = -180.0
    if (lonf> 180.0)  lonf =  180.0
    if (latdel<1.0/120.0) latdel = 1.0/120.0
    if (londel<1.0/120.0) londel = 1.0/120.0
    sites = 0
    if (latf>lat0+latdel/2.0) then
      ibox = int((latf-lat0)/latdel+0.5)
      iadj = latdel/2.0
    else
      ibox = 1
      iadj = 0.0
    endif
    if (lonf>lon0+londel/2.0) then
      jbox = int((lonf-lon0)/londel+0.5)
      jadj = londel/2.0
    else
      jbox = 1
      jadj = 0.0
    endif
    do i=1,ibox
      do j=1,jbox
        sites = sites + 1
        lat_lon(sites,1) = lat0 + real(i-1)*latdel + iadj
        lat_lon(sites,2) = lon0 + real(j-1)*londel + jadj
      enddo
    enddo
    read(98,*)
  else
    write(*,'('' PROGRAM TERMINATED'')')
    write(*,*) 'Error in the number of fields for site input'
    stop
  endif

endif

!----------------------------------------------------------------------!
! Read in soil characteristics. Defaults used when zero                !
!----------------------------------------------------------------------!
read(98,'(1000a)') st1
if ((n_fields(st1)<8).or.(n_fields(st1)>9)) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'Eight or nine fields must be present for the soil &
 &characteristics: sand %,'
  write(*,*) 'silt %, bulk density and Soil Depth (cm). Set to &
 &zero if default values'
  write(*,*) 'are required.'
  stop
endif
if (n_fields(st1)==8) then
  read(st1,*) (soil_chr(i),i=1,8)
  topsl = defaulttopsl
else
  read(st1,*) (soil_chr(i),i=1,8),topsl
  if (topsl<1.0e-6) topsl = defaulttopsl
endif
if ((soil_chr(5)<0.0).or.(soil_chr(6)<0.0).or.&
 (soil_chr(7)<0.0)) then
  l_b_and_c = .true.
else
  l_b_and_c = .false.
endif
read(98,'(1000a)') st1

!----------------------------------------------------------------------!
! Open diagnostics file.                                               !
!----------------------------------------------------------------------!
call OPEN_DIAG(stoutput)

!----------------------------------------------------------------------!
! Check sites against land mask and disregard when no land.            !
!----------------------------------------------------------------------!
call LAND_SITE_CHECK(st1,stoutput,sites,lat_lon,latdel,londel,du,stmask)

!----------------------------------------------------------------------!
! Read in type of landuse: 0 = defined by map; 1 = defined explicitly  !
! in the input file; 2 = natural vegetation based on average monthly   !
! temperatures.                                                        !
!----------------------------------------------------------------------!
read(98,*) ilanduse
if (ilanduse==1) then
! Use landuse defined in input file.
  call LANDUSE1(luse,yr0,yrf)
endif
if ((ilanduse<0).or.(ilanduse>2)) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'No landuse defined'
  write(*,*) '0:= Defined from a map.'
  write(*,*) '1:= Defined in the input file.'
  write(*,*) '2:= Natural vegetation.'
  stop
endif

call OPEN_SNAPSHOTS(stoutput,snp_no,snpshts)

!----------------------------------------------------------------------!
! Close the input file.                                                !
!----------------------------------------------------------------------!
close(98)
!----------------------------------------------------------------------!

ssp%nft = nft

end subroutine read_input_file






!**********************************************************************!
!                                                                      !
!                    read_param :: read_input                          !
!                    ------------------------                          !
!                                                                      !
! subroutine read_param(l_regional,site_out,year_out,stver)            !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Read internal parameters from "param.dat" file, and io
!! parameters from "misc_params.dat".
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine read_param(stver)
!**********************************************************************!
integer :: kode
character :: st1*1000,st2*1000,st3*1000,stver*1000
logical :: l_regional

!----------------------------------------------------------------------!
! Read 'misc_params.dat'.                                              !
!----------------------------------------------------------------------!
open(98,file='misc_params.dat',status='OLD',iostat=kode)

if (kode/=0) then
  write(*,*) ' File does not exist: "misc_params.dat"'
  write(*,*) ' Using screen output options: 1 0. '
  write(*,*) ' Using countries, not regions. '
  msp%site_out = 1
  msp%year_out = 0
  msp%l_regional = .false.
  print*,'hello'
else
  read(98,*)
  read(98,*) msp%site_out,msp%year_out
  read(98,*)
  read(98,*) msp%l_regional
endif
close(98)

!----------------------------------------------------------------------!
! Read internal parameters from "param.dat".                           !
!----------------------------------------------------------------------!
open(98,file='inc/param.dat',status='OLD',iostat=kode)
if (kode/=0) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) ' File does not exist: "param.dat"'
  stop
endif

read(98,*)
read(98,*) st1
call STRIPB(st1)
st2 = stver
call STRIPBS(st2,st3)
call STRIPB(st2)
st1=st1(1:blank(st1))
st2=st2(1:blank(st2))
if (stcmp(st1,st2)/=1) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) &
 'Version number mismatch between ''sdgvm0.f'' and ''param.dat''.'
  stop
endif

read(98,*)
read(98,*) tgp%p_sla        ! 4
read(98,*)
read(98,*) tgp%p_stemfr     ! 6
read(98,*)
read(98,*) tgp%p_rootfr     ! 8
read(98,*)
read(98,*) tgp%p_opt        !10
read(98,*)
read(98,*) tgp%p_stmin      !12
read(98,*)
read(98,*) tgp%p_laimem     !14
read(98,*)
read(98,*) tgp%p_resp       !16
read(98,*)
read(98,*) tgp%p_kscale     !18
read(98,*)
read(98,*) tgp%p_nu1,tgp%p_nu2,tgp%p_nu3,tgp%p_nu4
read(98,*)
read(98,*) tgp%p_nleaf
read(98,*)
read(98,*) tgp%p_dresp
read(98,*)
read(98,*) tgp%p_vm
read(98,*)
read(98,*) tgp%p_kgw
read(98,*)
read(98,*) tgp%p_v1,tgp%p_v2,tgp%p_v3
read(98,*)
read(98,*) tgp%p_j1,tgp%p_j2,tgp%p_j3
read(98,*)
read(98,*) tgp%p_pet
read(98,*)
read(98,*) tgp%p_bs
read(98,*)
read(98,*) tgp%p_et
read(98,*)
read(98,*) tgp%p_roff
read(98,*)
read(98,*) tgp%p_roff2
read(98,*)
read(98,*) tgp%p_fprob
read(98,*)
read(98,*) tgp%p_city_dep
close(98)

end subroutine read_param



end module read_input
