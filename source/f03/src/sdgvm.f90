! ============================================================================
! Name        : sdgvm.f90
! Author      : Mark Lomas
! Version     : 0
! Copyright   : Your copyright notice
! Description : sdgvm
! ============================================================================

program sdgvm

use real_precision,     only: dp
use dims,               only: max_cohorts, max_years
use system_state,       only: ssv
use pft_parameters
use state_methods
use sdgvm1
use read_input
use phenology_methods
use doly
use hydrology_methods
use light_methods
use soil_methods
use misc_values

implicit none

character(len=1000), parameter :: stver = 'VERSION current'
real(dp) :: defaulttopsl
parameter(defaulttopsl = 5.0)
logical :: check_closure
parameter(check_closure = .true.)

real(dp), dimension(max_cohorts) :: lai,evt,sresp,rof,gpp,ftprop,nppstoreold, &
 trn,lch,bioo,ht,soilc,soiln,minn,ftcov,covo,flow1,flow2, &
 leafnpp,stemnpp,rootnpp,bioleaf

real(dp) :: lat,lon,ca,resp,soilt,grassrc,tmp(12,31),prc(12,31),hum(12,31), &
 cld(12),latdel,londel,leafper,stemper,rootper,avnpp,avgpp,avlai,co20,co2f, &
 avrof,infix,avnppst,co2const,sum1,oscale,yield,co2(max_years),sumcov, &
 maxcov,maxbio,barerc,avtrn,firec,pet2,f2,f3,avevt,sumbio,kd,kx, &
 stembio,rootbio,sum,lutab(255,100),awl(4),nci(4),gsn,eemm,etmm,rn, &
 xlatf,xlatres,xlon0,xlonres,lat_lon(2000000,2),nleaf,canga,s1in,dp2, &
 h2o,adp(4),sfc(4),sw(4),sswc(4),nupc,swc,swf,ssm,xfprob,xno_fires, &
 nfix,daygpp,evap,tran,roff,pet,daily_out(50,max_cohorts,12,31),xx, &
 ans(12,31),interc,evbs,soil_chr(10),cluse(max_cohorts,max_years), &
 soil_chr2(10),fprob,q,qdirect,qdiff,topsl,hrs,fpr,resp_s,resp_r, &
 lmor_sc(3600,max_cohorts),leaflitter,laiinc,daynpp,resp_l,resp_m, &
 pm_tmp,pm_prc,pm_hum,pm_evbs,pm_swc,pm_swf,pm_ssm

integer :: sites,cycle,yr0,yrf,snp_no,snpshts(1000),snp_year,day,d, &
 isite,du,otagsn(50),otagsnft(50),nft,nat_map(8),ilanduse,siteno,iofn, &
 iofnft,iofngft,nomdos,i,j,k,ft,site,year, &
 covind,bioind,xyearf,mnth,fireres,xyear0,omav(50),year_out, &
 luse(max_years),fno,iyear,oymd,oymdft,iargc,yearind(max_years),idum, &
 outyears,thty_dys,xlatresn,xlonresn,day_mnth,yearv(max_years),nyears,narg,seed1, &
 seed2,seed3,spinl,yr0s,yr0p,yrfp,xseed1,site_dat,site_out, &
 country_id,outyears1,outyears2,budo(max_cohorts),seno(max_cohorts), &
 sit_grd,xtmpv(500,12,31),xprcv(500,12,31),xhumv(500,12,31),xcldv(500,12),co

character(len=1000) :: st1
character(len=1000) :: otags(100)
character(len=100) :: ofmt(100)
character(len=1000) :: stinput
character(len=1000) :: stoutput
character(len=1000) :: stinit
character(len=1000) :: stco2
character(len=1000) :: stmask
character(len=1000) :: country_name
character(len=1000) :: sttxdp
character(len=1000) :: stlu
character(len=1000) :: ststats
character(len=80) :: buff1

character(len=20) :: sttemp

logical :: initise,initiseo,speedc,crand,xspeedc,withcloudcover,l_clim,l_lu, &
 l_soil(20),l_stats,l_regional,out_cov,out_bio,out_bud,out_sen,l_b_and_c

!----------------------------------------------------------------------!
! Get input filename from the command line.                            !
!----------------------------------------------------------------------!
CALL GET_INPUT_FILENAME(buff1)
!buff1='../../SDGVM/output/f03_2.dat'

!----------------------------------------------------------------------!
! Read internal parameters from "param.dat" file, and io               !
! parameters from "misc_params.dat".                                   !
!----------------------------------------------------------------------!
call read_param(stver)

!----------------------------------------------------------------------!
! Read the input file.                                                 !
!----------------------------------------------------------------------!
call read_input_file(buff1,stinput,stinit,stoutput,sttxdp,stmask,stlu,&
 ststats, stco2,xlatf,xlon0,xlatres,xlonres,co2const,initise,initiseo,&
 speedc,xspeedc,xseed1,spinl,yr0s,cycle,crand,yr0p,yrfp,outyears,&
 nyears, yr0,yrf,yearind,idum,yearv,nomdos,otags,omav,ofmt,outyears1,&
 outyears2,oymd,otagsn,otagsnft,snpshts,snp_no,out_cov,out_bio,&
 out_bud,out_sen,lutab,grassrc,barerc,fireres,luse,l_b_and_c,&
 soil_chr,topsl,defaulttopsl,sites,latdel,londel,lat_lon,day_mnth,&
 thty_dys,xlatresn,xlonresn,ilanduse,nft,xyear0,xyearf,lmor_sc,&
 oymdft,iofnft,sit_grd,du,narg)

!----------------------------------------------------------------------!
! Open default output files.                                           !
!----------------------------------------------------------------------!
call open_default(stoutput)

!----------------------------------------------------------------------!
! Open optional output files.                                          !
!----------------------------------------------------------------------!
if (outyears2>0) call OPEN_OPTIONAL(stoutput,nft,out_cov,out_bio,&
 out_bud,out_sen)
iofngft = 200

co2(1) = 0.0

!----------------------------------------------------------------------!
! Open input and/or output files for state vector.                     !
!----------------------------------------------------------------------!
call OPEN_STATE(stinit,stoutput,initise,initiseo)

!----------------------------------------------------------------------!
!  Read co2 file.                                                      !
!----------------------------------------------------------------------!
call READCO2 (stco2,yr0,yrf,co2)
!----------------------------------------------------------------------!

call SET_GOUD_PARAMS()

site_dat = 0

do site=1,sites

  ssp%jday = 5000
  speedc = xspeedc
  lat = lat_lon(site,1)
  lon = lat_lon(site,2)
  seed1 = xseed1
  seed2 = 2*seed1
  seed3 = 3*seed1

!----------------------------------------------------------------------!
! Read in climate.                                                     !
!----------------------------------------------------------------------!
  call READ_CLIMATE(stinput,ststats,lat,lon,xlatf, &
 xlatres,xlatresn,xlon0,xlonres,xlonresn,yr0,yrf,xtmpv,xhumv,xprcv, &
 xcldv,isite,xyear0,xyearf,du,seed1,seed2,seed3,l_clim,l_stats,siteno, &
 day_mnth,thty_dys,sit_grd,withcloudcover)

!----------------------------------------------------------------------!
! Read in soil parameters.                                             !
!----------------------------------------------------------------------!
  call READ_SOIL(sttxdp,lat,lon,soil_chr,soil_chr2,du,l_soil)

!----------------------------------------------------------------------!
! Read in soil parameters.                                             !
!----------------------------------------------------------------------!
  call READ_LANDUSE(stlu,ilanduse,yr0,yrf,du,nft,lat,lon,lutab,luse,&
 cluse,l_lu)

!----------------------------------------------------------------------!
! Driving data exists 'if' statement.                                  !
!----------------------------------------------------------------------!
  if ((l_clim).and.(l_stats).and.(l_soil(1)).and.(l_soil(3)).and. &
 (l_soil(5)).and.(l_soil(8)).and.(l_lu)) then

  site_dat = site_dat + 1
  if (mod(site_dat,max(site_out,1))==min(1,site_out)-1) &
 write(*,'( '' Site no. '',i5,'', Lat ='',f7.3,'', Lon ='',f9.3)') site_dat,lat,lon

!----------------------------------------------------------------------!
! Write lat & lon for output files.                                    !
!----------------------------------------------------------------------!
    call WRITE_LAT_LON(lat,lon,nft,out_cov,out_bio,out_bud,out_sen, &
 oymdft,otagsn,oymd,otagsnft,outyears1,outyears2)

!----------------------------------------------------------------------!
! Computation of hydrological parameters.                              !
!----------------------------------------------------------------------!
    call WSPARAM(l_b_and_c,nupc,awl,kd,kx,nci,infix,adp,topsl,sfc,sw,sswc)

!----------------------------------------------------------------------!
! Extract initial and final C02 values.                                !
!----------------------------------------------------------------------!
    call CO2_0_F(co20,co2f,yearv,yr0,spinl,co2,co2const,nyears)

!----------------------------------------------------------------------!
! Extract the country/state corresponding to the site.                 !
!----------------------------------------------------------------------!
    call COUNTRY(stmask,lat,lon,country_name,country_id,l_regional)

!----------------------------------------------------------------------!
! Write site info to 'site_info.dat'.                                  !
!----------------------------------------------------------------------!
    write(12,'(a15,i6,f9.3,f9.3,1x,2f6.1,1x,2f7.1,f7.3,1x,3f7.3,f7.1,1x,200(2f6.1,1x))') &
 country_name,country_id,lat,lon,co20,co2f,ssp%sand,ssp%silt,ssp%bulk,ssp%wilt, &
 ssp%field,ssp%sat,ssp%soil_depth, &
 (cluse(ft,yearv(1)-yr0+1),cluse(ft,yearv(nyears)-yr0+1),ft=1,nft)

!----------------------------------------------------------------------!
! Initialise the system state.                                         !
!----------------------------------------------------------------------!
    call INITIALISE_STATE(initise,nft,cluse,xtmpv,soilt)

    do iyear=1,nyears

      year = yearv(iyear)

      ssp%iyear = iyear
      ssp%year = year

      nfix = infix

!----------------------------------------------------------------------!
! Set CO2 value 'ca' from 'co2' or 'co2const'.                         !
!----------------------------------------------------------------------!
      call SET_CO2(ca,spinl,iyear,speedc,co2,co2const,year,yr0)

      if (mod(iyear,max(msp%year_out,1))==min(1,msp%year_out)-1) then
        write(*,'('' Year no.'',2i5,'', ca = '',f6.2,'', cohorts = '', i3,''.'')') &
 iyear,year,ca,ssp%cohorts
      endif

!----------------------------------------------------------------------!
! Set 'tmp' 'hum' 'prc' and 'cld', and calc monthly and yearly avs.    *
!----------------------------------------------------------------------!
      call SET_CLIMATE(xtmpv,xprcv,xhumv,xcldv,withcloudcover,yearv,iyear, &
 tmp,prc,hum,cld,thty_dys,yr0,year)

!----------------------------------------------------------------------!
! Set land use through ftprop.                                         !
!----------------------------------------------------------------------!
      call SET_LANDUSE(ftprop,ilanduse,tmp,prc,nat_map,nft,cluse,year,yr0)

      call COVER(nft,tmp,prc,firec,fireres,fprob,ftprop,check_closure)

      call INITIALISE_NEW_COHORTS(nft,ftprop,check_closure)

      call MKDLIT()

      call RESTRICT_COHORT_NUMBERS()

!----------------------------------------------------------------------!
! Initialisations that were in doly at the beginning of the year       !
!----------------------------------------------------------------------!
      do ft=1,ssp%cohorts

        leafnpp(ft) = 0.0
        stemnpp(ft) = 0.0
        rootnpp(ft) = 0.0

        ssv(ft)%evp = 0.0

        budo(ft) = 0
        seno(ft) = 0

!----------------------------------------------------------------------!
! Height (m) ccn  ht(ft) = 0.807*(laimax(ft)**2.13655)                 !
!----------------------------------------------------------------------!
        ht(ft) = 0.807*(5.0**2.13655)
        if (ht(ft)>50.0)  ht(ft)=50.0

      enddo

!----------------------------------------------------------------------!
! START OF MONTH LOOP                                                  !
!----------------------------------------------------------------------!
      do mnth=1,12

        ssp%mnth = mnth

        call SUM_SOILCN(soilc,soiln,minn)

!----------------------------------------------------------------------!
! DAILY LOOP.                                                          !
!----------------------------------------------------------------------!
        do day=1,no_days(year,mnth,thty_dys)

!      CALL HERBIVORY(daily_out)

          ssp%day = day
          ssp%jday = ssp%jday + 1

          fpr=0.0

!----------------------------------------------------------------------!
! Updata daily temperature memory.                                     !
!----------------------------------------------------------------------!
          do i=1,199
            ssp%tmem(201-i) = ssp%tmem(200-i)
          enddo
          ssp%tmem(1) = tmp(mnth,day)

!----------------------------------------------------------------------!
! Radiation calculation.                                               !
!----------------------------------------------------------------------!
          hrs = dayl(lat,no_day(year,mnth,day,thty_dys))
          call PFD(lat,no_day(year,mnth,day,thty_dys),hrs,cld(mnth), &
     qdirect,qdiff,q)

!----------------------------------------------------------------------!
! Mix water resources.                                                 !
!----------------------------------------------------------------------!
          call MIX_WATER(ftcov,nft)

!----------------------------------------------------------------------!

          do ft=1,ssp%cohorts

          fpr=0.0
          ssp%cohort = ft

          if (ssv(ft)%cov>0.0) then

            call SET_MISC_VALUES(pft(ft)%sla,tmp(mnth,day))

!----------------------------------------------------------------------!
! nppstore leafnpp stemnpp rootnpp leaflit stemlit rootlit in mols
!----------------------------------------------------------------------!
            soilt = 0.97*soilt + 0.03*tmp(mnth,day)

            call DOLYDAY(tmp(mnth,day),prc(mnth,day),hum(mnth,day),ca, &
     soilc(ft),soiln(ft),minn(ft),adp,sfc,sw,sswc,awl,kd,kx,daygpp,resp_l,lai(ft), &
     evap,tran,roff,interc,evbs,flow1(ft),flow2(ft),year,mnth,day,pet,ht(ft), &
     thty_dys,ft,lmor_sc(:,pft(ft)%itag), &
     nleaf,leaflitter,hrs,q,qdirect,qdiff,fpr,canga,gsn,rn,check_closure)

            call EVAPOTRANSPIRATION(tmp(mnth,day),hum(mnth,day),rn,canga,gsn,hrs,eemm,etmm)
            pet = eemm
            pet2 = pet

            call HYDROLOGY(adp,sfc,sw,sswc,awl,kd,kx,eemm,etmm,pet2,prc(mnth,day), &
     s1in,tmp(mnth,day),ssv(ft)%lai%tot,evap,tran,roff,interc,evbs,f2,f3,ft)
            flow1(ft) = flow1(ft) + f2/10.0
            flow2(ft) = flow2(ft) + f3/10.0

            call PHENOLOGY(yield,laiinc)

            xx = ssv(ft)%nppstore(1)
            call ALLOCATION(laiinc,daygpp,resp_l,lmor_sc(:,pft(ft)%itag),resp, &
     leaflitter,stemnpp(ft),rootnpp(ft),resp_s,resp_r,resp_m,check_closure)

            ssv(ft)%slc = ssv(ft)%slc + leaflitter*ssv(ft)%cov

            call SOIL_DYNAMICS2(pet,prc(mnth,day),tmp(mnth,day),f2/10.0,f3/10.0,nfix, &
     nci,sresp(ft),lch(ft),ca,site,year,yr0,yrf,speedc,ft,check_closure)

!----------------------------------------------------------------------!
            swc = ssv(ft)%soil_h2o(1)+ssv(ft)%soil_h2o(2)+ssv(ft)%soil_h2o(3)+ssv(ft)%soil_h2o(4)
            swf = (swc-sw(1)-sw(2)-sw(3)-sw(4))/(sfc(1)+sfc(2)+sfc(3)+sfc(4)-sw(1)-sw(2)-sw(3)-sw(4))
            ssm = ssv(ft)%soil_h2o(1)/10.0/topsl
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Compute fire inputs and call fire routine on the last day of the     *
! month.                                                               !
!----------------------------------------------------------------------!
!       CALL BFIRE_INPUTS(pm_tmp,pm_prc,pm_hum,pm_evbs,pm_swc,pm_swf,pm_ssm, &
! tmp(mnth,day),prc(mnth,day),hum(mnth,day),evbs,swc,swf,ssm,day)
!       IF ((day == 30).AND.(ft == ssp%cohorts)) CALL FIRE_B(xfprob,pm_tmp,pm_prc,pm_hum, &
! pm_evbs,pm_swc,pm_swf,pm_ssm,lat,mnth,ssp%cohort)

!----------------------------------------------------------------------!
! Set daily memories for output.                                       !
!----------------------------------------------------------------------!
            daily_out(1,ft,mnth,day) = ssv(ft)%lai%tot
            daily_out(2,ft,mnth,day) = roff
            daily_out(3,ft,mnth,day) = evap+tran
            daily_out(4,ft,mnth,day) = tran
            daily_out(5,ft,mnth,day) = ssv(ft)%npp
!(leafnpp(ft) + stemnpp(ft) +  &
! rootnpp(ft) + ssv(ft)%nppstore(1) - leafold - stemold - rootold - nppsold)*12.0
            daily_out(6,ft,mnth,day) = daygpp
            daily_out(7,ft,mnth,day) = sresp(ft)
            daily_out(8,ft,mnth,day) = daily_out(5,ft,mnth,day) - &
 daily_out(7,ft,mnth,day)
            daily_out(9,ft,mnth,day) = tmp(mnth,day)
            daily_out(10,ft,mnth,day) = prc(mnth,day)
            daily_out(11,ft,mnth,day) = hum(mnth,day)
            daily_out(12,ft,mnth,day) = ssv(ft)%nppstore(1)
            daily_out(13,ft,mnth,day) = swf
            daily_out(14,ft,mnth,day) = pet
            daily_out(15,ft,mnth,day) = interc
            daily_out(16,ft,mnth,day) = evbs
            daily_out(17,ft,mnth,day) = &
 min(1.0,ssv(ft)%soil_h2o(1)/10.0/topsl)
            daily_out(18,ft,mnth,day) = swc
            daily_out(19,ft,mnth,day) = 1.0
!        daily_out(19,ft,mnth,day) = resp
            daily_out(20,ft,mnth,day) = qdirect*hrs*3600.0
            daily_out(21,ft,mnth,day) = qdiff*hrs*3600.0
            daily_out(22,ft,mnth,day) = 1.0
!        daily_out(22,ft,mnth,day) = nleaf
!        daily_out(23,ft,mnth,day) = leaflit(ft) - lflitold
            daily_out(24,ft,mnth,day) = cld(mnth)
            daily_out(25,ft,mnth,day) = fpr
            daily_out(26,ft,mnth,day) = 0.0
            if (pft(ft)%sla > 0.0)  daily_out(26,ft,mnth,day) = &
 ssv(ft)%lai%tot*12.0/pft(ft)%sla/18.0
            daily_out(27,ft,mnth,day) = &
 ssv(ft)%bio(1) + ssv(ft)%stem%tot + ssv(ft)%nppstore(1)
            daily_out(28,ft,mnth,day) = ssv(ft)%bio(2) + ssv(ft)%root%tot
            daily_out(29,ft,mnth,day) = soilc(ft)
            daily_out(30,ft,mnth,day) = soiln(ft)
            daily_out(31,ft,mnth,day) = lch(ft)
            daily_out(32,ft,mnth,day) = firec/360.0
            daily_out(33,ft,mnth,day) = yield

!----------------------------------------------------------------------!
            if ((ssv(ft)%bb==day+(mnth-1)*30).and.(budo(ft)==0)) &
 budo(ft) = ssv(ft)%bb
            if ((ssv(ft)%ss==day+(mnth-1)*30).and.(seno(ft)==0)) &
 seno(ft) = ssv(ft)%ss

          endif ! cov > 0
          enddo
!----------------------------------------------------------------------!
! End of the cohort loop.                                              !
!----------------------------------------------------------------------!
        enddo
!----------------------------------------------------------------------!
! End of the daily loop.                                               !
!----------------------------------------------------------------------!
      enddo
!----------------------------------------------------------------------!
! End of the monthly loop.                                             !
!----------------------------------------------------------------------!

  if (1 < 0) then
  do ft=1,ssp%cohorts

    nppstoreold(ft) = 0.0
    if (abs(leafnpp(ft) + stemnpp(ft) + rootnpp(ft) + &
 ssv(ft)%nppstore(1) - nppstoreold(ft))>1.0e-6) then
      ssv(ft)%npl = leafnpp(ft)/(leafnpp(ft) + stemnpp(ft) + &
 rootnpp(ft) + ssv(ft)%nppstore(1) - nppstoreold(ft))*100.0
      ssv(ft)%nps = stemnpp(ft)/(leafnpp(ft) + stemnpp(ft) + &
 rootnpp(ft) + ssv(ft)%nppstore(1) - nppstoreold(ft))*100.0
      ssv(ft)%npr = rootnpp(ft)/(leafnpp(ft) + stemnpp(ft) + &
 rootnpp(ft) + ssv(ft)%nppstore(1) - nppstoreold(ft))*100.0
    else
      ssv(ft)%npl = 0.0
      ssv(ft)%nps = 0.0
      ssv(ft)%npr = 0.0
      ssv(ft)%npp = 0.0
    endif
    ssv(ft)%npp = (leafnpp(ft) + stemnpp(ft) + rootnpp(ft) + &
 ssv(ft)%nppstore(1) - nppstoreold(ft))

  enddo
  endif

  call GROWTH(nft,lai,stembio,rootbio,check_closure)

!----------------------------------------------------------------------!
! Average outputs by cover proportions.                                !
!----------------------------------------------------------------------!
  avlai   = 0.0
  avnpp   = 0.0
  avnppst = 0.0
  avrof   = 0.0
  avtrn   = 0.0
  avevt   = 0.0
  avgpp   = 0.0
  do ft=1,ssp%cohorts
    avnppst = avnppst + ssv(ft)%cov*ssv(ft)%nppstore(1)
    avrof   = avrof   + ssv(ft)%cov*rof(ft)
    avtrn   = avtrn   + ssv(ft)%cov*trn(ft)
    avevt   = avevt   + ssv(ft)%cov*evt(ft)
    avgpp   = avgpp   + ssv(ft)%cov*gpp(ft)
  enddo

!----------------------------------------------------------------------!
  sumbio = 0.0
  bioind = 0
  covind = 0
  maxbio = 0.0
  maxcov = 0.0
  leafper = 0.0
  stemper = 0.0
  rootper = 0.0

  do ft=1,nft
    bioo(ft) = 0.0
    covo(ft) = 0.0
    if (ssp%co2ftmap(ft,1) > 0) then
      budo(ft) = ssv(ssp%co2ftmap(ft,2))%bb
      seno(ft) = ssv(ssp%co2ftmap(ft,2))%ss
    else
      budo(ft) = 0
      seno(ft) = 0
    endif
    do j=1,ssp%co2ftmap(ft,1)
      co = ssp%co2ftmap(ft,j+1)
      do mnth=1,12
        do day=1,no_days(year,mnth,thty_dys)
          bioo(ft) = bioo(ft) + (daily_out(26,co,mnth,day) + &
 daily_out(27,co,mnth,day) + daily_out(28,co,mnth,day))*ssv(co)%cov
        enddo
      enddo
      covo(ft) = covo(ft) + ssv(co)%cov
      leafper = leafper + ssv(co)%npl*ssv(co)%cov
      stemper = stemper + ssv(co)%nps*ssv(co)%cov
      rootper = rootper + ssv(co)%npr*ssv(co)%cov
    enddo
    bioo(ft) = bioo(ft)/360.0
    sumbio = sumbio + bioo(ft)
    if (covo(ft)>maxcov) then
      maxcov = covo(ft)
      covind = ft
    endif
    if (bioo(ft)>maxbio) then
      maxbio = bioo(ft)
      bioind = ft
    endif
  enddo

!----------------------------------------------------------------------!
! Write outputs.                                                       !
!----------------------------------------------------------------------!
! General.                                                             !
!----------------------------------------------------------------------!
  if (iyear>=nyears-outyears+1) then
     write(21,'('' '',f8.1)',advance='NO') outputs(daily_out,1,'Max')       !lai
    write(22,'('' '',f8.1)',advance='NO') outputs(daily_out,5,'Add')       !npp
    write(23,'('' '',f8.1)',advance='NO') outputs(daily_out,29,'Average')  !soil carbon
    write(24,'('' '',f8.1)',advance='NO') outputs(daily_out,30,'Average')  !soil nitrogen
    write(25,'('' '',f8.1)',advance='NO') outputs(daily_out,5,'Add') - outputs(daily_out,7,'Add') !nep
    write(26,'('' '',f8.1)',advance='NO') min(swc,9999.0)
    write(27,'('' '',f8.1)',advance='NO') outputs(daily_out,26,'Average') + &
 outputs(daily_out,27,'Average') + outputs(daily_out,28,'Average') ! biomass
    write(28,'('' '',i2)',advance='NO')   bioind
    write(29,'('' '',i2)',advance='NO')   covind
    write(30,'('' '',f8.1)',advance='NO') 1.0
    write(31,'('' '',f8.1)',advance='NO') outputs(daily_out,2,'Add')
    write(32,'('' '',f8.2)',advance='NO') outputs(daily_out,32,'Add')
    write(33,'('' '',f8.1)',advance='NO') avnppst
    write(34,'('' '',f8.1)',advance='NO') outputs(daily_out,27,'Average')  !stem biomass
    write(35,'('' '',f8.1)',advance='NO') outputs(daily_out,28,'Average')  !root biomass
    write(36,'('' '',f8.1)',advance='NO') leafper
    write(37,'('' '',f8.1)',advance='NO') stemper
    write(38,'('' '',f8.1)',advance='NO') rootper
    write(39,'('' '',f8.1)',advance='NO') outputs(daily_out,7,'Add')       !soil respiration
    write(40,'('' '',f8.1)',advance='NO') outputs(daily_out,3,'Add')
    write(41,'('' '',f8.1)',advance='NO') outputs(daily_out,6,'Add')       !gpp
    write(42,'('' '',f8.3)',advance='NO') outputs(daily_out,31,'Add')      !soil leaching
    write(43,'('' '',f8.1)',advance='NO') outputs(daily_out,10,'Add') !precipitation
    write(44,'('' '',f9.3)',advance='NO') outputs(daily_out,5,'Add') - outputs(daily_out,7,'Add') - &
 outputs(daily_out,32,'Add') - outputs(daily_out,31,'Add') - outputs(daily_out,33,'Add') !nbp avnpp-avsresp-firec-avlch-avyield
    write(45,'('' '',f8.1)',advance='NO') outputs(daily_out,4,'Add')
    write(46,'('' '',f8.5)',advance='NO') fprob
    write(47,'('' '',f8.2)',advance='NO') outputs(daily_out,9,'Average') !temperature
    write(48,'('' '',f8.2)',advance='NO') outputs(daily_out,11,'Average') !humidity
  endif

!----------------------------------------------------------------------!
! Write optional cov bio bud sen.                                      !
!----------------------------------------------------------------------!
  iofn = iofngft
  if (iyear>=nyears-outyears2+1) then
    do ft=1,nft
      if (out_cov) then
        iofn = iofn + 1
        write(iofn,'('' '',f8.6)',advance='NO') covo(ft)
      endif
      if (out_bio) then
        iofn = iofn + 1
        write(iofn,'('' '',f8.1)',advance='NO') bioo(ft)
      endif
      if (out_bud) then
        iofn = iofn + 1
        write(iofn,'('' '',i8)',advance='NO') budo(ft)
      endif
      if (out_sen) then
        iofn = iofn + 1
        write(iofn,'('' '',i8)',advance='NO') seno(ft)
      endif
    enddo
  endif

!----------------------------------------------------------------------!
! Write monthly PIXEL outputs.                                         !
!----------------------------------------------------------------------!
  iofn = 400
  if (iyear>=nyears-outyears1+1) then
    if ((oymd==1).or.(oymd==3)) then
      do i=1,50
        if (otagsn(i)==1) then
          do mnth=1,12
            ans(mnth,1) = 0.0
          enddo
          do ft=1,nft
            do j=1,ssp%co2ftmap(ft,1)
              co = ssp%co2ftmap(ft,j+1)
              do mnth=1,12
                if (omav(i)==1) then
                  oscale = 1.0/real(no_days(year,mnth,thty_dys))
                else
                  oscale = 1.0
                endif
                do day=1,no_days(year,mnth,thty_dys)
                  ans(mnth,1) = ans(mnth,1) + &
 daily_out(i,co,mnth,day)*ssv(co)%cov*oscale
                enddo
              enddo
            enddo
          enddo
          iofn = iofn + 1
          sum1 = 0.0
          if (omav(i)==1) then
            oscale = 1.0/12.0
          else
            oscale = 1.0
          endif
          do mnth=1,12
            sum1 = sum1 + ans(mnth,1)*oscale
          enddo
          write(iofn,ofmt(i)) year,(ans(mnth,1),mnth=1,12),sum1
        endif
      enddo
    endif

!----------------------------------------------------------------------!
! Write daily PIXEL outputs.                                           !
!----------------------------------------------------------------------!
    if ((oymd==2).or.(oymd==3)) then
      do i=1,50
        if (otagsn(i)==1) then
          iofn = iofn + 1
          write(iofn,'(i4)') year
          do mnth=1,12
            do day=1,no_days(year,mnth,thty_dys)
              ans(mnth,day) = 0.0
            enddo
          enddo
          do ft=1,nft
            do j=1,ssp%co2ftmap(ft,1)
              co = ssp%co2ftmap(ft,j+1)
              do mnth=1,12
                do day=1,no_days(year,mnth,thty_dys)
                  ans(mnth,day) = ans(mnth,day) + &
 daily_out(i,co,mnth,day)*ssv(co)%cov
                enddo
              enddo
            enddo
          enddo
          do mnth=1,12
            write(iofn,ofmt(50+i)) (ans(mnth,day),&
 day=1,no_days(year,mnth,thty_dys))
          enddo
        endif
      enddo
    endif
  endif

  iofn = iofnft
  if (iyear>=nyears-outyears2+1) then
!----------------------------------------------------------------------!
! Write monthly SUBPIXEL outputs.                                      !
!----------------------------------------------------------------------!
    if ((oymdft==1).or.(oymdft==3)) then
      do i=1,50
        if (otagsnft(i)==1) then
          do ft=1,nft
            sumcov = 0.0
            do mnth=1,12
              ans(mnth,1) = 0.0
            enddo
            do j=1,ssp%co2ftmap(ft,1)
              co = ssp%co2ftmap(ft,j+1)
              do mnth=1,12
                if (omav(i)==1) then
                  oscale = 1.0/real(no_days(year,mnth,thty_dys))
                else
                  oscale = 1.0
                endif
                do day=1,no_days(year,mnth,thty_dys)
                  ans(mnth,1) = ans(mnth,1) + &
 daily_out(i,co,mnth,day)*oscale*ssv(co)%cov
                enddo
              enddo
              sumcov = sumcov + ssv(co)%cov
            enddo
            iofn = iofn + 1
            if (sumcov > 0.0) then
              do mnth=1,12
                ans(mnth,1) = ans(mnth,1)/sumcov
              enddo
            endif
            write(iofn,ofmt(i)) year,(ans(mnth,1),mnth=1,12)
          enddo
        endif
      enddo
    endif

!----------------------------------------------------------------------!
! Write daily SUBPIXEL outputs.                                        !
!----------------------------------------------------------------------!
    if ((oymdft==2).or.(oymdft==3)) then
      do i=1,50
        if (otagsnft(i)==1) then
          do ft=1,nft
            iofn = iofn + 1
            write(iofn,'(i4)') year
            do mnth=1,12
              do day=1,no_days(year,mnth,thty_dys)
                ans(mnth,day) = 0.0
              enddo
            enddo
            sumcov = 0.0
            do j=1,ssp%co2ftmap(ft,1)
              co = ssp%co2ftmap(ft,j+1)
              do mnth=1,12
                do day=1,no_days(year,mnth,thty_dys)
                  ans(mnth,day) = ans(mnth,day) + &
 daily_out(i,co,mnth,day)*ssv(co)%cov
                enddo
              enddo
              sumcov = sumcov + ssv(co)%cov
            enddo
            if (sumcov > 0.0) then
              do mnth=1,12
                do day=1,no_days(year,mnth,thty_dys)
                  ans(mnth,day) = ans(mnth,day)/sumcov
                enddo
              enddo
            endif
            do mnth=1,12
              write(iofn,ofmt(50+i)) (ans(mnth,day),&
 day=1,no_days(year,mnth,thty_dys))
            enddo
          enddo
        endif
      enddo
    endif
  endif

!----------------------------------------------------------------------!
! Write snapshots
!----------------------------------------------------------------------!
  snp_year = 1
  if (snp_no>0) then
    if (snp_year<=snp_no) then
      if (year==snpshts(snp_year)) then
        fno = 100 + (snp_year-1)*4
        write(fno+1,'(F7.3,F9.3,4500F9.1)')  &
 lat,lon,(((ssv(ft)%bio(j),j=1,2),ft=1,pft(k)%mort),k=1,nft)
        write(fno+2,'(F7.3,F9.3,1500F12.9)') &
 lat,lon,((ssv(ft)%cov,ft=1,pft(j)%mort),j=1,nft)
        write(fno+3,'(F7.3,F9.3,1500F12.7)') &
 lat,lon,((ssv(ft)%ppm,ft=1,pft(j)%mort),j=1,nft)
        write(fno+4,'(F7.3,F9.3,1500F8.3)')  &
 lat,lon,((ssv(ft)%hgt,ft=1,pft(j)%mort),j=1,nft)
        snp_year = snp_year + 1
      endif
    endif
  endif

!***********************************************************************
enddo ! year loop
!***********************************************************************

!----------------------------------------------------------------------!
!                      End record for default output files             !
!----------------------------------------------------------------------!
  do i=21,48
    write(i,*)
  enddo
  iofn = iofngft
  if (outyears2>0) then
    do ft=1,nft
      if (out_cov) then
        iofn = iofn + 1
        write(iofn,*)
      endif
      if (out_bio) then
        iofn = iofn + 1
        write(iofn,*)
      endif
      if (out_bud) then
        iofn = iofn + 1
        write(iofn,*)
      endif
      if (out_sen) then
        iofn = iofn + 1
        write(iofn,*)
      endif
    enddo
  endif

!----------------------------------------------------------------------!
!                             Write initialisation output              !
!----------------------------------------------------------------------!
!  IF (initiseo) THEN
!  ENDIF

else

  write(11,*) '                  clm stt ssc blk wfs dep lus'
  write(11,'(f7.3,f9.3,1x,20L4)') lat,lon,l_clim,l_stats, &
 l_soil(1),l_soil(3),l_soil(5),l_soil(8),l_lu

!----------------------------------------------------------------------!
! End of climate data exists 'if' statement.                           !
!----------------------------------------------------------------------!
endif

!***********************************************************************
enddo ! site loop
!***********************************************************************

stop
!----------------------------------------------------------------------!
! Open file to record version number, command line, input file and     *
! parameter file.                                                      !
!----------------------------------------------------------------------!
open(13,file=stoutput(1:blank(stoutput))//'/simulation.dat')
call GETARG(0,buff1)
call STRIPB(stver)
write(13,'(A)') stver(1:28)
close(98)
write(13,*)

write(13,'(''************************************************'')')
write(13,'(''* Command line arguments                       !'')')
write(13,'(''************************************************'')')

do i=0,narg
  call GETARG(i,buff1)
  do j=1,80
    st1(j:j) = buff1(j:j)
  enddo
  write(13,'(A,'' '')',advance='NO') st1(1:blank(st1))
enddo
write(13,*)

write(13,'(''************************************************'')')
write(13,'(''* Parameter file (param.dat)                   !'')')
write(13,'(''************************************************'')')

open(98,file='inc/param.dat',status='OLD')
80    read(98,'(A)',end=90) st1
write(13,'(A)') st1(1:last_blank(st1))
goto 80
90    continue
close(98)
write(13,*)

write(13,'(''************************************************'')')
write(13,'(''* Input file                                   !'')')
write(13,'(''************************************************'')')

call GETARG(1,buff1)
open(98,file=buff1,status='OLD')
65    read(98,'(A)',end=70) st1
write(13,'(A)') st1(1:last_blank(st1))
goto 65
70    continue
close(98)

write(13,'(''************************************************'')')
!----------------------------------------------------------------------!

close(11)
close(12)
close(13)
close(14)
do i=21,48
  close(i)
enddo

iofn = iofngft
if (outyears2>0) then
  do ft=1,nft
    if (out_cov) then
      iofn = iofn + 1
      close(iofn)
    endif
    if (out_bio) then
      iofn = iofn + 1
      close(iofn)
    endif
    if (out_bud) then
      iofn = iofn + 1
      close(iofn)
    endif
    if (out_sen) then
      iofn = iofn + 1
      close(iofn)
    endif
  enddo
endif

do i=70,89
  close(i)
enddo
do i=91,93
  close(i)
enddo
do i=98,99
  close(i)
enddo

if (snp_no>0) then
  do i=1,snp_no
    fno = 100+(i-1)*4
    close(fno+1)
    close(fno+1)
    close(fno+3)
    close(fno+4)
  enddo
endif

end program sdgvm


