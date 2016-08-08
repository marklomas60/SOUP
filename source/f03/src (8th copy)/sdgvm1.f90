module sdgvm1

use real_precision
use dims
use data
use weather_generator
use veg_dynamics
use func

implicit none

contains

!**********************************************************************!
!                                                                      !
!                          outputs :: sdgvm1                           !
!                          -----------------                           !
!                                                                      !
! real(dp) function outputs(daily_out,ind,st1)                         !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Compute outputs from daily_out.
!! @details Four forms of outputs are available and are specified by
!!'st1' which can take any one of the values
!! ['Add', 'Average', 'Max', 'Min'].
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
real(dp) function outputs(daily_out,ind,st1)
!**********************************************************************!
real(dp) ::  daily_out(50,max_cohorts,12,31),omax,ans
integer :: ind,co,mnth,day
character(len=3) :: st1
!----------------------------------------------------------------------!

outputs = 0.0
if ((st1(1:1) == 'A').and.(st1(2:2) == 'd')) then
  do co=1,ssp%cohorts
    ans = 0.0
    do mnth=1,12
      do day=1,30
        ans = ans + daily_out(ind,co,mnth,day)
      enddo
    enddo
    outputs = outputs + ans*ssv(co)%cov
  enddo
elseif ((st1(1:1) == 'A').and.(st1(2:2) == 'v')) then
  do co=1,ssp%cohorts
    ans = 0.0
    do mnth=1,12
      do day=1,30
        ans = ans + daily_out(ind,co,mnth,day)
      enddo
    enddo
    outputs = outputs + ans*ssv(co)%cov/360.0
  enddo
elseif ((st1(1:1) == 'M').or.(st1(2:2) == 'a')) then
  do co=1,ssp%cohorts
    ans = -1000000.0
    do mnth=1,12
      do day=1,30
        if (daily_out(ind,co,mnth,day) > ans) ans = daily_out(ind,co,mnth,day)
      enddo
    enddo
    outputs = outputs + ans*ssv(co)%cov
  enddo
elseif ((st1(1:1) == 'M').or.(st1(2:2) == 'i')) then
  do co=1,ssp%cohorts
    ans = 1000000.0
    do mnth=1,12
      do day=1,30
        if (daily_out(ind,co,mnth,day) < ans) ans = daily_out(ind,co,mnth,day)
      enddo
    enddo
    outputs = outputs + ans*ssv(co)%cov
  enddo
else
  write(*,'(''Error in the third argument to outputs in SDGVM0.'')')
  write(*,'(''Arguments: Add, Average, Max, Min.'')')
  stop
endif

end function outputs





!**********************************************************************!
!                                                                      !
!                          set_landuse :: sdgvm1                       !
!                          ---------------------                       !
!                                                                      !
! subroutine set_landuse(ftprop,ilanduse,tmp,prc,nat_map,nft,cluse,    !
! year,yr0)                                                            !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Set lanuse by allocating the cover of new ground.
!! @details Assign gound available for new growth, through disturbance
!! and mortality, given by ftprop. The values of the desired cover are
!! given in the array 'cluse'. Note that any changes in the desired
!! are only met if new ground is large enough to make up the difference.
!! Over time the cover will reach any desired value and once reached
!! the cover will remain constant.
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine set_landuse(ftprop,ilanduse,tmp,prc,nat_map,nft,cluse,year,yr0)
!**********************************************************************!
integer :: ilanduse,ft,nft,nat_map(8)
integer :: year,yr0
real(dp) :: tmp(12,31),prc(12,31),cluse(max_cohorts,max_years) &
 ,ftprop(max_cohorts)
!----------------------------------------------------------------------!

if (ilanduse==2) then
  call NATURAL_VEG(tmp,prc,ftprop,nat_map)
else ! 0 and 1
  ftprop(1) = 100.0
  do ft=2,nft
    if (check_ft_grow(tmp,ssv(1)%chill,ssv(1)%dschill,ft)==1) then
      ftprop(ft) = cluse(ft,year-yr0+1)
      ftprop(1) = ftprop(1) - ftprop(ft)
    else
      ftprop(ft) = 0.0
    endif
  enddo
  if (ftprop(1)<0.0) then
    do ft=2,nft
      ftprop(ft) = ftprop(ft)*100.0/(100.0 - ftprop(1))
    enddo
    ftprop(1) = 0.0
  endif
endif

end subroutine set_landuse





!**********************************************************************!
!                                                                      !
!                          set_climate :: sdgvm1                       !
!                          ---------------------                       !
!                                                                      !
! subroutine set_climate(xtmpv,xprcv,xhumv,xcldv,withcloudcover,yearv, !
! iyear,tmp,prc,hum,cld,thty_dys,yr0,year)                             !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Set 'tmp' 'hum' 'prc' and 'cld'.
!! @details Assign climate from the raw data which is help in 'x???v'
!! arrays.
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine set_climate(xtmpv,xprcv,xhumv,xcldv,withcloudcover,yearv,&
 iyear,tmp,prc,hum,cld,thty_dys,yr0,year)
!**********************************************************************!
integer,dimension(500,12,31) :: xtmpv,xprcv,xhumv
integer,dimension(500,12) :: xcldv
integer :: yearv(max_years),iyear,mnth,day,thty_dys,yr0,year
logical :: withcloudcover
real(dp) :: tmp(12,31),prc(12,31),hum(12,31),cld(12)
!----------------------------------------------------------------------!

do mnth=1,12
  do day=1,no_days(year,mnth,thty_dys)
    tmp(mnth,day) = real(xtmpv(yearv(iyear)-yr0+1,mnth,day))/100.0
    prc(mnth,day) = real(xprcv(yearv(iyear)-yr0+1,mnth,day))/10.0
    hum(mnth,day) = real(xhumv(yearv(iyear)-yr0+1,mnth,day))/100.0
    if (withcloudcover) then
      cld(mnth) = real(xcldv(yearv(iyear)-yr0+1,mnth))/1000.0
    else
      cld(mnth) = 0.5
    endif
  enddo
enddo

do mnth=1,12
  do day=1,no_days(year,mnth,thty_dys)
    if (hum(mnth,day)<30.0)  hum(mnth,day) = 30.0
    if (hum(mnth,day)>95.0)  hum(mnth,day) = 95.0
  enddo
enddo

end subroutine set_climate





!**********************************************************************!
!                                                                      !
!                          set_co2 :: sdgvm1                           !
!                          -----------------                           !
!                                                                      !
! subroutine set_co2(ca,spinl,iyear,speedc,co2,co2const,year,yr0)      !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Set CO2 value 'ca'.
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine set_co2(ca,spinl,iyear,speedc,co2,co2const,year,yr0)
!**********************************************************************!
real(dp) :: ca,co2const,co2(max_years)
integer :: spinl,iyear,year,yr0
logical :: speedc
!----------------------------------------------------------------------!

if ((spinl>0).and.(iyear>spinl)) then
  speedc = .false.
  ca = co2(year-yr0+1)
else
  if (co2const>0.0) then
    ca = co2const
  else
    ca = co2(year-yr0+1)
  endif
endif

end subroutine set_co2





!**********************************************************************!
!                                                                      !
!                        read_landuse :: sdgvm1                        !
!                        ----------------------                        !
!                                                                      !
! subroutine read_landuse(stlu,ilanduse,yr0,yrf,du,nft,lat,lon,lutab,  !
! luse,cluse,l_lu)                                                     !
!                                                                      !
!> @brief
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine read_landuse(stlu,ilanduse,yr0,yrf,du,nft,lat,lon,lutab,&
 luse,cluse,l_lu)
!**********************************************************************!
integer :: ilanduse,icontinuouslanduse,yr0,yrf,du,year,ft,nft,&
 luse(max_years),nat_map(8)
real(dp) :: lat,lon,cluse(max_cohorts,max_years),lutab(255,100),sum
character :: stlu*1000,fttags(100)*1000,st2*1000
logical :: l_lu
!----------------------------------------------------------------------!

icontinuouslanduse = 1

if (ilanduse==0) then
  if (icontinuouslanduse==0) then
    call EX_LU(stlu,lat,lon,luse,yr0,yrf,du)
! Create the continuous land use (cluse)
    do year=yr0,yrf
      do ft=1,nft
        cluse(ft,year-yr0+1) = lutab(luse(year-yr0+1),ft)
      enddo
    enddo
  else            
    call EX_CLU(stlu,lat,lon,nft,lutab,cluse,yr0,yrf,du,l_lu)
  endif
elseif (ilanduse==1) then
  l_lu = .true.
  do year=yr0,yrf
    do ft=1,nft
      cluse(ft,year-yr0+1) = lutab(luse(year-yr0+1),ft)
    enddo
  enddo
elseif (ilanduse==2) then
  write(*,*) 'Checking natural vegetation types exist:'
  write(*,*) 'BARE CITY C3 C4 Ev_Bl Ev_Nl Dc_Bl Dc_Nl.'
  l_lu = .true.
  st2 = 'BARE'
  nat_map(1) = ntags(fttags,st2)
  st2 = 'CITY'
  nat_map(2) = ntags(fttags,st2)
  st2 = 'C3'
  nat_map(3) = ntags(fttags,st2)
  st2 = 'C4'
  nat_map(4) = ntags(fttags,st2)
  st2 = 'Ev_Bl'
  nat_map(5) = ntags(fttags,st2)
  st2 = 'Ev_Nl'
  nat_map(6) = ntags(fttags,st2)
  st2 = 'Dc_Bl'
  nat_map(7) = ntags(fttags,st2)
  st2 = 'Dc_Nl'
  nat_map(8) = ntags(fttags,st2)
  do year=yr0,yrf
    do ft=1,nft
      cluse(ft,year-yr0+1) = 0.0
    enddo
  enddo
  cluse(1,1) = 100.0
endif

end subroutine read_landuse





!**********************************************************************!
!                                                                      !
!                       read_climate :: sdgvm1                         !
!                       ----------------------                         !
!                                                                      !
! subroutine read_climate(stinput,ststats,lat,lon,xlatf,               !
! xlatres,xlatresn,xlon0,xlonres,xlonresn,yr0,yrf,xtmpv,xhumv,xprcv,   !
! xcldv,isite,xyear0,xyearf,du,seed1,seed2,seed3,l_clim,l_stats,siteno,!
! day_mnth,thty_dys,sit_grd,withcloudcover))                           !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Read internal parameters from "param.dat" file, and io
!! parameters from "misc_params.dat".
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine read_climate(stinput,ststats,lat,lon,xlatf, &
 xlatres,xlatresn,xlon0,xlonres,xlonresn,yr0,yrf,xtmpv,xhumv,xprcv, &
 xcldv,isite,xyear0,xyearf,du,seed1,seed2,seed3,l_clim,l_stats,siteno, &
 day_mnth,thty_dys,sit_grd,withcloudcover)
!**********************************************************************!
character :: stinput*1000,ststats*1000
real(dp) :: lat,lon,xlatf,xlatres,xlon0,xlonres
integer :: xlatresn,xlonresn,yr0,yrf,isite,xyear0,xyearf,du,seed1,seed2,seed3, &
 siteno,day_mnth,thty_dys,sit_grd
integer :: xtmpv(500,12,31),xprcv(500,12,31),xhumv(500,12,31),xcldv(500,12)
logical :: l_clim,l_stats,withcloudcover
!----------------------------------------------------------------------!

if ((day_mnth==1).and.(thty_dys==1)) then
  call EX_CLIM(stinput,lat,lon,xlatf,xlatres,xlatresn,xlon0,xlonres, &
 xlonresn,yr0,yrf,xtmpv,xhumv,xprcv,isite,xyear0,xyearf,siteno,du)
  withcloudcover=.false.
  ssp%latres = xlatres
  ssp%lonres = xlonres
elseif ((day_mnth==0).and.(thty_dys==1).and.(sit_grd==0)) then
  call EX_CLIM_WEATHER_GENERATOR(stinput,ststats,lat,lon,xlatf, &
 xlatres,xlatresn,xlon0,xlonres,xlonresn,yr0,yrf,xtmpv,xhumv,xprcv, &
 xcldv,isite,xyear0,xyearf,du,seed1,seed2,seed3,l_clim,l_stats)
  withcloudcover=.true.
  ssp%latres = xlatres
  ssp%lonres = xlonres
elseif ((day_mnth==1).and.(thty_dys==0)) then
  call EX_CLIM_SITE(stinput,yr0,yrf,xtmpv,xhumv,xprcv,xyear0,xyearf)
  withcloudcover=.false.
  siteno = 1
else
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'Error defining climate to read'
  stop
endif

end subroutine read_climate





!**********************************************************************!
!                                                                      !
!                     get_input_filename :: sdgvm1                     !
!                     ----------------------------                     !
!                                                                      !
! subroutine get_input_filename(st1)                                   !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Read internal parameters from "param.dat" file, and io
!! parameters from "misc_params.dat".
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine get_input_filename(st1)
!**********************************************************************!
character :: st1*80
integer :: IARGC
!----------------------------------------------------------------------!

if (IARGC()>0) then
  call GETARG(1,st1)
else
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) ' Input file must be given as an argument.'
  stop
endif

end subroutine get_input_filename





!**********************************************************************!
!                                                                      !
!                          country :: sdgvm1                           !
!                          -----------------                           !
!                                                                      !
! subroutine country(stmask,lat,lon,country_name,country_id,l_regional)!
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Read internal parameters from "param.dat" file, and io
!! parameters from "misc_params.dat".
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine country(stmask,lat,lon,country_name,country_id,l_regional)
!**********************************************************************!
real(dp) :: lat,lon,adj_lat,adj_lon,xlat,xlon
integer :: country_id,ilat,ilon,ans(360),x,i
character :: country_name*1000,stmask*1000
logical :: l_regional
!----------------------------------------------------------------------!

ilat = int(91.0 - lat)
ilon = int(lon + 181.0)

open(99,file=stmask(1:blank(stmask))//'/country.dat',status='old')
  if (ilat>1) then
    do i=1,ilat-1
      read(99,*)
    enddo
  endif
  read(99,*) ans
  country_id = ans(ilon)
close(99)

!----------------------------------------------------------------------!
! If OCEAN was found try the nearest adjacent squares.                 *
!----------------------------------------------------------------------!
if (country_id==0) then

! Nearest lateral.

  xlat = lat+500-real(int(lat+500.0))
  if (xlat>0.5) then
    adj_lat = 1.0
  else
    adj_lat = -1.0
  endif

  xlon = lon+500-real(int(lon+500.0))
  if (xlon>0.5) then
    adj_lon = 1.0
  else
    adj_lon = -1.0
  endif

  if (abs(xlat-0.5)>abs(xlon-0.5)) then
    adj_lon = 0.0
  else
    adj_lat = 0.0
  endif

  ilat = int(91.0 - (lat + adj_lat))
  ilon = int((lon + adj_lat) + 181.0)

  if (ilon==0) ilon = 360
  if (ilon==361) ilon = 1
  if (ilat==0) ilat = 1
  if (ilat==91) ilat = 90

  open(99,file=stmask(1:blank(stmask))//'/country.dat',status='old')
    if (ilat>1) then
      do i=1,ilat-1
        read(99,*)
      enddo
    endif
    read(99,*) ans
    country_id = ans(ilon)
  close(99)

endif

if (country_id==0) then

! Next nearest lateral.

  if (abs(adj_lat)<0.5) then
    xlat = lat+500-real(int(lat+500.0))
    if (xlat>0.5) then
      adj_lat = 1.0
    else
      adj_lat = -1.0
    endif
  endif

  if (abs(adj_lon)<0.5) then
    xlon = lon+500-real(int(lon+500.0))
    if (xlon>0.5) then
      adj_lon = 1.0
    else
      adj_lon = -1.0
    endif
    adj_lat = 0.0
  endif

  if (abs(adj_lat)>0.5)  adj_lon = 0.0

  ilat = int(91.0 - (lat + adj_lat))
  ilon = int((lon + adj_lat) + 181.0)

  if (ilon==0) ilon = 360
  if (ilon==361) ilon = 1
  if (ilat==0) ilat = 1
  if (ilat==91) ilat = 90

  open(99,file=stmask(1:blank(stmask))//'/country.dat',status='old')
  if (ilat>1) then
    do i=1,ilat-1
      read(99,*)
    enddo
  endif
  read(99,*) ans
  country_id = ans(ilon)
  close(99)

endif

if (country_id==0) then

! Nearest diagonal

  xlat = lat+500-real(int(lat+500.0))
  if (xlat>0.5) then
    adj_lat = 1.0
  else
    adj_lat = -1.0
  endif

  xlon = lon+500-real(int(lon+500.0))
  if (xlon>0.5) then
    adj_lon = 1.0
  else
    adj_lon = -1.0
  endif

  ilat = int(91.0 - (lat + adj_lat))
  ilon = int((lon + adj_lat) + 181.0)

  if (ilon==0) ilon = 360
  if (ilon==361) ilon = 1
  if (ilat==0) ilat = 1
  if (ilat==91) ilat = 90

  open(99,file=stmask(1:blank(stmask))//'/country.dat',status='old')
    if (ilat>1) then
      do i=1,ilat-1
        read(99,*)
      enddo
    endif
    read(99,*) ans
    country_id = ans(ilon)
  close(99)

endif

!----------------------------------------------------------------------!
! Regional or country switch.                                          *
if (.not.(l_regional)) country_id = country_id-mod(country_id,100)
!----------------------------------------------------------------------!

open(99,file=stmask(1:blank(stmask))//'/country_id.dat',status='old')
10    continue
  read(99,'(i6,5x,a15)') x,country_name
  if (x==country_id) goto 20
goto 10
20    continue
close(99)

end subroutine country





!**********************************************************************!
!                                                                      !
!                          outputs :: sdgvm1                           !
!                          -----------------                           !
!                                                                      !
! subroutine co2_0_f(co20,co2f,yearv,yr0,spinl,co2,co2const,nyears)    !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Read internal parameters from "param.dat" file, and io
!! parameters from "misc_params.dat".
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine co2_0_f(co20,co2f,yearv,yr0,spinl,co2,co2const,nyears)
!**********************************************************************!
integer :: iyear,year,spinl,nyears,yr0,yearv(max_years)
real(dp) :: co2const,co20,co2f,co2(max_years)
!----------------------------------------------------------------------!

iyear = 1
year = yearv(iyear)
if ((spinl>0).and.(iyear>spinl)) then
  co20 = co2(year-yr0+1)
else
  if (co2const>0.0) then
    co20 = co2const
  else
    co20 = co2(year-yr0+1)
  endif
endif
iyear = nyears
year = yearv(iyear)
if ((spinl>0).and.(iyear>spinl)) then
  co2f = co2(year-yr0+1)
else
  if (co2const>0.0) then
    co2f = co2const
  else
    co2f = co2(year-yr0+1)
  endif
endif

!**********************************************************************!
end subroutine co2_0_f
!**********************************************************************!





!**********************************************************************!
!                                                                      !
!                          read_soil :: sdgvm1                         !
!                          -------------------                         !
!                                                                      !
! subroutine read_soil(sttxdp,lat,lon,soil_chr,soil_chr2,du,l_soil)    !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Read internal parameters from "param.dat" file, and io
!! parameters from "misc_params.dat".
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine read_soil(sttxdp,lat,lon,soil_chr,soil_chr2,du,l_soil)
!**********************************************************************!
real(dp) :: lat,lon,soil_chr(10),soil_chr2(10)
integer :: du
character :: sttxdp*1000
logical :: l_soil(20)
!----------------------------------------------------------------------!

call EX_SOIL(sttxdp,lat,lon,soil_chr2,du,l_soil)
if (soil_chr(1)>0.01) then
  ssp%sand = soil_chr(1)
  ssp%silt = soil_chr(2)
  l_soil(1) = .true.
  l_soil(2) = .true.
else
  ssp%sand = soil_chr2(1)
  ssp%silt = soil_chr2(2)
endif
ssp%clay = 100.0 - ssp%sand - ssp%silt

if (soil_chr(3)>0.01) then
  ssp%bulk  = soil_chr(3)
  l_soil(3) = .true.
else
  ssp%bulk = soil_chr2(3)
endif

if (soil_chr(4)>0.01) then
  ssp%orgc  = soil_chr(4)
  l_soil(4) = .true.
else
  ssp%orgc = soil_chr2(4)
endif

if (soil_chr(5)>0.01) then
  ssp%wilt  = soil_chr(5)
  ssp%field = soil_chr(6)
  ssp%sat   = soil_chr(7)
  l_soil(5) = .true.
  l_soil(6) = .true.
  l_soil(7) = .true.
else
  ssp%wilt = soil_chr2(5)
  ssp%field = soil_chr2(6)
  ssp%sat = soil_chr2(7)
endif

if (soil_chr(8)>0.01) then
  ssp%soil_depth = soil_chr(8)
  l_soil(8) = .true.
else
  ssp%soil_depth = soil_chr2(8)
endif

end subroutine read_soil





!**********************************************************************!
!                                                                      !
!                       land_site_check :: sdgvm1                      !
!                       -------------------------                      !
!                                                                      !
! subroutine land_site_check(st1,st2,sites,lat_lon,latdel,londel,du,   !
! stmask)                                                              !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Read internal parameters from "param.dat" file, and io
!! parameters from "misc_params.dat".
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine land_site_check(st1,st2,sites,lat_lon,latdel,londel,du,&
 stmask)
!**********************************************************************!
character :: st1*1000,st2*1000,stmask*1000
integer :: sites,du,i,site
real(dp) :: lat_lon(2000000,2),latdel,londel,lat,lon
logical :: lc
!----------------------------------------------------------------------!

open(14,file=st2(1:blank(st2))//'/land_sites.dat')
i = 0
do site=1,sites
  if (mod(site,1)==0) write(*,*) sites,site
  lat = lat_lon(site,1)
  lon = lat_lon(site,2)
  call lorc(du,lat,lon,latdel,londel,stmask,lc)
  if (lc) then
    write(14,'(f7.3,f9.3)') lat,lon
    i = i + 1
    lat_lon(i,1) = lat_lon(site,1)
    lat_lon(i,2) = lat_lon(site,2)
  endif
enddo
close(14)
sites = i

write(*,'(i8,'' potential land sites.'')') sites
if ((n_fields(st1)==3).or.(n_fields(st1)==5).or.(n_fields(st1)==7)) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'The site line of the input file has 3, 5 or 7 fileds.'
  stop
endif

end subroutine land_site_check





!**********************************************************************!
!                                                                      !
!                       output_options :: sdgvm1                       !
!                       ------------------------                       !
!                                                                      !
! subroutine output_options(nomdos,otags,omav,ofmt)                    !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Read internal parameters from "param.dat" file, and io
!! parameters from "misc_params.dat".
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine output_options(nomdos,otags,omav,ofmt)
!**********************************************************************!
integer :: nomdos,omav(50)
character :: otags(100)*1000,ofmt(100)*100

!----------------------------------------------------------------------!
! Read in monthly and/or daily output options.                         !
!----------------------------------------------------------------------!
! file tags.                                                           !
!----------------------------------------------------------------------!
otags(1)  = 'lai'
otags(2)  = 'rof'
otags(3)  = 'evt'
otags(4)  = 'trn'
otags(5)  = 'npp'
otags(6)  = 'gpp'
otags(7)  = 'srp'
otags(8)  = 'nep'
otags(9)  = 'tmp'
otags(10) = 'prc'
otags(11) = 'hum'
otags(12) = 'nps'
otags(13) = 'swf'
otags(14) = 'pet'
otags(15) = 'int'
otags(16) = 'bse'
otags(17) = 'ssm'
otags(18) = 'swc'
otags(19) = 'rsp'
otags(20) = 'qdr'
otags(21) = 'qdf'
otags(22) = 'lfn'
otags(23) = 'lfl'
otags(24) = 'cld'
otags(25) = 'fpr'

nomdos = 25

!----------------------------------------------------------------------!
! Averaged output.                                                     !
!----------------------------------------------------------------------!
omav(1) = 1
omav(2) = 0
omav(3) = 0
omav(4) = 0
omav(5) = 0
omav(6) = 0
omav(7) = 0
omav(8) = 0
omav(9) = 1
omav(10) = 0
omav(11) = 1
omav(12) = 1
omav(13) = 1
omav(14) = 0
omav(15) = 0
omav(16) = 0
omav(17) = 1
omav(18) = 1
omav(19) = 0
omav(20) = 1
omav(21) = 1
omav(22) = 1
omav(23) = 1
omav(24) = 1

!----------------------------------------------------------------------!
! Monthly output format.                                               !
!----------------------------------------------------------------------!
ofmt(1)   ='(i4,1x,12f8.2,f10.1)'
ofmt(2)   ='(i4,1x,12f8.2,f10.1)'
ofmt(3)   ='(i4,1x,12f8.2,f10.1)'
ofmt(4)   ='(i4,1x,12f8.2,f10.1)'
ofmt(5)   ='(i4,1x,12f8.1,f10.3)'
ofmt(6)   ='(i4,1x,12f8.1,f10.1)'
ofmt(7)   ='(i4,1x,12f8.1,f10.1)'
ofmt(8)   ='(i4,1x,12f8.1,f10.1)'
ofmt(9)   ='(i4,1x,12f8.1,f10.1)'
ofmt(10)  ='(i4,1x,12f8.1,f10.1)'
ofmt(11)  ='(i4,1x,12f8.2,f10.1)'
ofmt(12)  ='(i4,1x,12f8.1,f10.1)'
ofmt(13)  ='(i4,1x,12f8.3,f10.3)'
ofmt(14)  ='(i4,1x,12f8.2,f10.1)'
ofmt(15)  ='(i4,1x,12f8.2,f10.1)'
ofmt(16)  ='(i4,1x,12f8.2,f10.1)'
ofmt(17)  ='(i4,1x,12f8.2,f10.5)'
ofmt(18)  ='(i4,1x,12f8.2,f10.3)'
ofmt(19)  ='(i4,1x,12f8.1,f10.1)'
ofmt(20)  ='(i4,1x,12f8.1,f10.1)'
ofmt(21)  ='(i4,1x,12f8.1,f10.1)'
ofmt(22)  ='(i4,1x,12f8.3,f10.1)'
ofmt(23)  ='(i4,1x,12f8.3,f10.1)'
ofmt(24)  ='(i4,1x,12f8.3,f10.1)'
ofmt(25)  ='(i4,1x,12f8.3,f10.1)'

!----------------------------------------------------------------------!
! Daily output format.                                                 !
!----------------------------------------------------------------------!
ofmt(51)  ='(31f7.2)'
ofmt(52)  ='(31f7.2)'
ofmt(53)  ='(31f7.2)'
ofmt(54)  ='(31f7.2)'
ofmt(55)  ='(31f7.1)'
ofmt(56)  ='(31f7.1)'
ofmt(57)  ='(31f7.1)'
ofmt(58)  ='(31f7.1)'
ofmt(59)  ='(31f7.1)'
ofmt(60)  ='(31f7.1)'
ofmt(61)  ='(31f7.2)'
ofmt(62)  ='(31f7.1)'
ofmt(63)  ='(31f7.3)'
ofmt(64)  ='(31f7.2)'
ofmt(65)  ='(31f7.2)'
ofmt(66)  ='(31f7.2)'
ofmt(67)  ='(31f7.4)'
ofmt(68)  ='(31f7.2)'
ofmt(69)  ='(31f7.1)'
ofmt(70)  ='(31f7.2)'
ofmt(71)  ='(31f7.2)'
ofmt(72)  ='(31f7.3)'
ofmt(73)  ='(31f7.3)'
ofmt(74)  ='(31f7.3)'
ofmt(75)  ='(31f7.3)'

end subroutine output_options





!**********************************************************************!
!                                                                      !
!                       set_pixel_out :: sdgvm1                        !
!                       -----------------------                        !
!                                                                      !
! subroutine set_pixel_out(st1,st2,outyears,nomdos,otagsn,otags,oymd)  !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Read internal parameters from "param.dat" file, and io
!! parameters from "misc_params.dat".
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine set_pixel_out(st1,st2,outyears,nomdos,otagsn,otags,oymd)
!**********************************************************************!
integer :: otagsn(50),nomdos,i,j,l,ii,oymd,outyears
character :: otags(100)*1000,st1*1000,st2*1000,st3*1000
character :: st4*1000,st5*1000
!----------------------------------------------------------------------!

call STRIPB(st1)
ii = n_fields(st1)

if (ii>1) then
  call STRIPBS(st1,st2)
  call STRIPB(st1)
  st3 = 'MONTHLY'
  st4 = 'DAILY'
  st5 = 'ALL'
  if (stcmp(st1,st3)==1) then
    oymd = 1
  elseif (stcmp(st1,st4)==1) then
    oymd = 2
  elseif (stcmp(st1,st5)==1) then
    oymd = 3
  else
    write(*,'('' PROGRAM TERMINATED'')')
    write(*,*) & 
 'The second field of line 15 of the input file must read MONTHLY DAILY or ALL.'
    write(*,'('' "'',A,''"'')') st1(1:blank(st1))
    write(*,*) 'Output variable options:'
    write(*,'(1x,20a4)') (otags(j),j=1,15)
    write(*,'(1x,20a4)') (otags(j),j=16,nomdos)
    stop
  endif
  call STRIPBS(st1,st2)

  ii = n_fields(st1)
  if (ii>=1) then
    call STRIPBN(st1,outyears)
    if (outyears==-9999) then
      write(*,'('' PROGRAM TERMINATED'')')
      write(*,*) &
 'The third field of the ''PIXEL'' line must contain an integer, if it exists.'
      stop
    endif
  elseif (ii==0) then
    outyears = 0
  endif
else
  outyears = 0
  oymd = 0
endif

do i=1,50
  otagsn(i) = 0
enddo
if (ii>1) then
  call STRIPB(st1)
  do i=1,ii-1
    call STRIPBS(st1,st2)
    l = ntags(otags,st2)
    if (l/=-9999) then
      otagsn(l) = 1
    else
      write(*,'('' PROGRAM TERMINATED'')')
      write(*,*) 'Error in tag name on the ''PIXEL'' line.'
      write(*,'('' "'',A,''"'')') st2(1:blank(st2))
      write(*,*) 'Available tags:'
      write(*,'(1x,20a4)') (otags(j),j=1,15)
      write(*,'(1x,20a4)') (otags(j),j=16,nomdos)
      stop
    endif 
 enddo
else
  do i=1,nomdos
    otagsn(i) = 1
  enddo
endif

end subroutine set_pixel_out





!**********************************************************************!
!                                                                      !
!                      set_subpixel_out :: sdgvm1                      !
!                      --------------------------                      !
!                                                                      !
! subroutine set_subpixel_out(st1,st2,outyears,nomdos, &               !
! otagsnft,otags,oymdft,out_cov,out_bio,out_bud,out_sen)               !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Read internal parameters from "param.dat" file, and io
!! parameters from "misc_params.dat".
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine set_subpixel_out(st1,st2,outyears,nomdos, &
 otagsnft,otags,oymdft,out_cov,out_bio,out_bud,out_sen)
!**********************************************************************!
integer :: otagsnft(50),nomdos,i,j,l,ii,oymdft,outyears
character :: otags(100)*1000,st1*1000,st2*1000,st3*1000
character :: st4*1000,st5*1000
logical :: out_cov,out_bio,out_bud,out_sen
!----------------------------------------------------------------------!

out_cov = .false.
out_bio = .false.
out_bud = .false.
out_sen = .false.

call STRIPB(st1)
ii = n_fields(st1)
if (ii>1) then
  call STRIPBS(st1,st2)
  call STRIPB(st1)
  st3 = 'MONTHLY'
  st4 = 'DAILY'
  st5 = 'ALL'
  if (stcmp(st1,st3)==1) then
    oymdft = 1
  elseif (stcmp(st1,st4)==1) then
    oymdft = 2
  elseif (stcmp(st1,st5)==1) then
    oymdft = 3
  else
    write(*,'('' PROGRAM TERMINATED'')')
    write(*,*) &
 'The second field of line 16 of the input file must read MONTHLY DAILY or ALL.'
    write(*,'('' "'',A,''"'')') st1(1:blank(st1))
    write(*,*) 'Output variable options:'
    write(*,'(1x,20a4)') (otags(j),j=1,15)
    write(*,'(1x,20a4)') (otags(j),j=16,nomdos),'cov ','bio ','bud ','sen '
    stop
  endif
  call STRIPBS(st1,st2)

  ii = n_fields(st1)
  if (ii>=1) then
    call STRIPBN(st1,outyears)
    if (outyears==-9999) then
      write(*,'('' PROGRAM TERMINATED'')')
      write(*,*) & 
 'The third field of the ''SUB_PIXEL'' line must contain an integer, if it exists.'
      stop
    endif
  elseif (ii==0) then
    outyears = 0
  endif
else
  outyears = 0
  oymdft = 0
endif

do i=1,50
  otagsnft(i) = 0
enddo

if (ii>1) then
  call STRIPB(st1)
  do i=1,ii-1
    call STRIPBS(st1,st2)
    st3 = 'cov'
    if (stcmp(st2,st3)==1) then
      out_cov = .true.
    else
      st3 = 'bio'
      if (stcmp(st2,st3)==1) then
        out_bio = .true.
      else
        st3 = 'bud'
        if (stcmp(st2,st3)==1) then
          out_bud = .true.
        else
          st3 = 'sen'
          if (stcmp(st2,st3)==1) then
            out_sen = .true.
          else
            l = ntags(otags,st2)
            if (l/=-9999) then
              otagsnft(l) = 1
            else
              write(*,'('' PROGRAM TERMINATED'')')
              write(*,*) 'Error in tag name in the ''SUBPIXEL'' line.'
              write(*,'('' "'',A,''"'')') st2(1:blank(st2))
              write(*,*) 'Available tag names:'
              write(*,'(1x,20a4)') (otags(j),j=1,15)
              write(*,'(1x,20a4)') (otags(j),j=16,nomdos),'cov ','bio ','bud ','sen '
              stop
            endif
          endif
        endif
      endif
    endif
  enddo
else
  do i=1,nomdos
    otagsnft(i) = 1
  enddo
  out_cov = .true.
  out_bio = .true.
  out_bud = .true.
  out_sen = .true.
endif

end subroutine set_subpixel_out





!**********************************************************************!
!                                                                      !
!                          mkdlit :: sdgvm1                            !
!                          ----------------                            !
!                                                                      !
! subroutine mkdlit()                                                  !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Read internal parameters from "param.dat" file, and io
!! parameters from "misc_params.dat".
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine mkdlit()
!**********************************************************************!
integer :: ft
!----------------------------------------------------------------------!

do ft=1,ssp%cohorts
 ssv(ft)%dslc = 0.0
 ssv(ft)%drlc = 0.0
 ssv(ft)%dsln = 0.0
 ssv(ft)%drln = 0.0
enddo

do ft=1,ssp%cohorts
  ssv(ft)%dslc = ssv(ft)%dslc + ssv(ft)%slc
  ssv(ft)%drlc = ssv(ft)%drlc + ssv(ft)%rlc
  ssv(ft)%dsln = ssv(ft)%dsln + ssv(ft)%sln
  ssv(ft)%drln = ssv(ft)%drln + ssv(ft)%rln
enddo

end subroutine mkdlit





!**********************************************************************!
!                                                                      !
!                          sum_soilcn :: sdgvm1                        !
!                          --------------------                        !
!                                                                      !
! subroutine sum_soilcn(soilc,soiln,minn)                              !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Read internal parameters from "param.dat" file, and io
!! parameters from "misc_params.dat".
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine sum_soilcn(soilc,soiln,minn)
!**********************************************************************!
real(dp) :: soilc(max_cohorts),soiln(max_cohorts),minn(max_cohorts)
integer :: i,ft
!----------------------------------------------------------------------!

do ft=1,ssp%cohorts
  soilc(ft) = 0.0
  soiln(ft) = 0.0
  minn(ft) = 0.0
  do i=1,8
    soilc(ft) = soilc(ft) + ssv(ft)%c(i)
    soiln(ft) = soiln(ft) + ssv(ft)%n(i)
  enddo
  minn(ft) = ssv(ft)%minn(1) + ssv(ft)%minn(2)
  soiln(ft) = soiln(ft) + minn(ft)
enddo

end subroutine sum_soilcn





!**********************************************************************!
!                                                                      !
!                          swap :: sdgvm1                              !
!                          --------------                              !
!                                                                      !
! subroutine swap(ic0,in0,iminn,c0,n0,minn)                            !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Read internal parameters from "param.dat" file, and io
!! parameters from "misc_params.dat".
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine swap(ic0,in0,iminn,c0,n0,minn)
!**********************************************************************!
real(dp) :: ic0(8),in0(8),iminn(3),c0(8),n0(8),minn(3)
integer :: i
!----------------------------------------------------------------------!

do i=1,8
   ic0(i) = c0(i)
   in0(i) = n0(i)
enddo
do i=1,3
  iminn(i) = minn(i)
enddo

end subroutine swap





!**********************************************************************!
!                                                                      !
!                          mix_water :: sdgvm1                         !
!                          -------------------                         !
!                                                                      !
! subroutine mix_water(ftcov,nft)                                      !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Read internal parameters from "param.dat" file, and io
!! parameters from "misc_params.dat".
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine mix_water(ftcov,nft)
!**********************************************************************!
real(dp) :: ftcov(max_cohorts),mix(max_cohorts)
real(dp) :: share1,share2,share3,share4,share5,share6,sum
integer :: nft,ft,mapft
!----------------------------------------------------------------------!

do ft=1,nft
  ftcov(ft) = 0.0
enddo

do ft=1,ssp%cohorts
  ftcov(pft(ft)%itag) = ftcov(pft(ft)%itag) + ssv(ft)%cov
enddo

do ft=1,ssp%cohorts
  mix(ft) = 1.0 - (1.0 - pft(ft)%mix)**(20.0/360.0)
enddo

share1 = 0.0
share2 = 0.0
share3 = 0.0
share4 = 0.0
share5 = 0.0
share6 = 0.0
sum = 0.0
do ft=1,ssp%cohorts
  mapft = pft(ft)%itag
  share1 = share1 + ssv(ft)%soil_h2o(1)*ssv(ft)%cov*mix(ft)
  share2 = share2 + ssv(ft)%soil_h2o(2)*ssv(ft)%cov*mix(ft)
  share3 = share3 + ssv(ft)%soil_h2o(3)*ssv(ft)%cov*mix(ft)
  share4 = share4 + ssv(ft)%soil_h2o(4)*ssv(ft)%cov*mix(ft)
  share5 = share5 + ssv(ft)%snow       *ssv(ft)%cov*mix(ft)
  share6 = share6 + ssv(ft)%l_snow     *ssv(ft)%cov*mix(ft)
  ssv(ft)%soil_h2o(1) = ssv(ft)%soil_h2o(1)*ssv(ft)%cov*(1.0 - mix(ft))
  ssv(ft)%soil_h2o(2) = ssv(ft)%soil_h2o(2)*ssv(ft)%cov*(1.0 - mix(ft))
  ssv(ft)%soil_h2o(3) = ssv(ft)%soil_h2o(3)*ssv(ft)%cov*(1.0 - mix(ft))
  ssv(ft)%soil_h2o(4) = ssv(ft)%soil_h2o(4)*ssv(ft)%cov*(1.0 - mix(ft))
  ssv(ft)%snow        = ssv(ft)%snow       *ssv(ft)%cov*(1.0 - mix(ft))
  ssv(ft)%l_snow      = ssv(ft)%l_snow     *ssv(ft)%cov*(1.0 - mix(ft))
  sum = sum + ssv(ft)%cov*mix(ft)
enddo

do ft=1,ssp%cohorts
  mapft = pft(ft)%itag
  if (ftcov(mapft)>0.0) then
    if (sum>0.0) then
      ssv(ft)%soil_h2o(1) = (ssv(ft)%soil_h2o(1) + share1*ftcov(mapft)*mix(ft)/sum)/ftcov(mapft)
      ssv(ft)%soil_h2o(2) = (ssv(ft)%soil_h2o(2) + share2*ftcov(mapft)*mix(ft)/sum)/ftcov(mapft)
      ssv(ft)%soil_h2o(3) = (ssv(ft)%soil_h2o(3) + share3*ftcov(mapft)*mix(ft)/sum)/ftcov(mapft)
      ssv(ft)%soil_h2o(4) = (ssv(ft)%soil_h2o(4) + share4*ftcov(mapft)*mix(ft)/sum)/ftcov(mapft)
      ssv(ft)%snow        = (ssv(ft)%snow        + share5*ftcov(mapft)*mix(ft)/sum)/ftcov(mapft)
      ssv(ft)%l_snow      = (ssv(ft)%l_snow      + share6*ftcov(mapft)*mix(ft)/sum)/ftcov(mapft)
    else
      ssv(ft)%soil_h2o(1) = ssv(ft)%soil_h2o(1)/ftcov(mapft)
      ssv(ft)%soil_h2o(2) = ssv(ft)%soil_h2o(2)/ftcov(mapft)
      ssv(ft)%soil_h2o(3) = ssv(ft)%soil_h2o(3)/ftcov(mapft)
      ssv(ft)%soil_h2o(4) = ssv(ft)%soil_h2o(4)/ftcov(mapft)
      ssv(ft)%snow        = ssv(ft)%snow       /ftcov(mapft)
      ssv(ft)%l_snow      = ssv(ft)%l_snow     /ftcov(mapft)
    endif
  endif
enddo

end subroutine mix_water



end module sdgvm1

