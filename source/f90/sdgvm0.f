*----------------------------------------------------------------------*
*                       SDGVM0 VERSION 070607                          *
*----------------------------------------------------------------------*
      IMPLICIT NONE

      CHARACTER stver*1000
      PARAMETER(stver = 'VERSION 070607')

      INCLUDE 'array_dims.inc'

      REAL*8 defaulttopsl
      PARAMETER(defaulttopsl = 5.0d0)

      REAL*8 lat,lon,dep,ca,npp(maxnft),lai(maxnft),evp(maxnft),resp
      REAL*8 gpp(maxnft),sresp(maxnft),evt(maxnft),soilt,grassrc
      REAL*8 tmp(12,31),prc(12,31),hum(12,31),cld(12),latdel,londel
      REAL*8 leafper,stemper,rootper,avnpp,avgpp,avlai,avdof,co20,co2f
      REAL*8 avrof,infix,nppsold,avnppst,co2const,ic0(8),in0(8),iminn(3)
      REAL*8 swcold,swcnew,is1,is2,is3,is4,isn,ilsn,tc0(8),tn0(8),sum1
      REAL*8 tminn(3),ts1,ts2,ts3,ts4,tsn,tlsn,oscale,yield(maxnft)
      REAL*8 c0(8,maxnft),n0(8,maxnft),minn(3,maxnft),c0v(8),n0v(8)
      REAL*8 minnv(3),s1(maxnft),s2(maxnft),s3(maxnft),s4(maxnft)
      REAL*8 sn(maxnft),lsn(maxnft),tmin,prctot,npl(maxnft),wtswc
      REAL*8 nps(maxnft),npr(maxnft),dof(maxnft),rof(maxnft),mnthprc(12)
      REAL*8 co2(maxyrs),maxcov,maxbio,ftprop(maxnft),tmem(200),c3old
      REAL*8 cov(maxage,maxnft),bio(maxage,2,maxnft),ppm(maxage,maxnft)
      REAL*8 hgt(maxage,maxnft),wdt(maxage,maxnft),leafdp(3600,maxnft)
      REAL*8 stemdp(360,maxnft),rootdp(360,maxnft),leafv(3600),c4old
      REAL*8 stemv(1000),rootv(1000),yeartmp,yearprc,yearhum,barerc
      REAL*8 laimax(maxnft),laimnth(12,maxnft),ht(maxnft),avtrn,yrpet
      REAL*8 trn(maxnft),slc(maxnft),rlc(maxnft),sln(maxnft),rln(maxnft)
      REAL*8 mnthtmp(12),firec,dslc,drlc,dsln,drln,lch(maxnft),avlch
      REAL*8 ftgr0(maxnft),ftgrf(maxnft),mnthhum(12),ccheck,iadj,jadj
      REAL*8 bioo(maxnft),covo(maxnft),avevt,avsresp,ftppm0(maxnft)
      REAL*8 tsoilc,tsoiln,soilc(maxnft),soiln(maxnft),isoilc,isoiln
      REAL*8 sumbio,ans1,ftstmx(maxnft),leaflit(maxnft),stemlit(maxnft)
      REAL*8 rootlit(maxnft),ftwd(maxnft),ftxyl(maxnft),ftpd(maxnft)
      REAL*8 ftsla(maxnft),ftcov(maxnft),lon0,lonf,ftrat(maxnft),kd,kx
      REAL*8 stembio,rootbio,sum,solcoo,biotoo,lutab(255,100),awl(4)
      REAL*8 nppstore(maxnft),ftlmor(maxnft),bioleaf(maxnft),nci(4)
      REAL*8 nppstoreold(maxnft),xlatf,xlatres,xlon0,xlonres,lat0,latf
      REAL*8 lat_lon(2000000,2),nppstorx(maxnft),avppm,avpet
      REAL*8 fpet(maxnft),SETARANDOM,temp_lat,temp_lon
      REAL*8 nppstor2(maxnft),maxlai(maxnft),ftbb0(maxnft),wtfc,wtwp
      REAL*8 ftbbmax(maxnft),ftbblim(maxnft),ftsslim(maxnft),nleaf
      REAL*8 flow1(maxnft),flow2(maxnft),h2o,adp(4),sfc(4),sw(4),sswc(4)
      REAL*8 leafnpp(maxnft),stemnpp(maxnft),rootnpp(maxnft),srespm,lchm
      REAL*8 evapm(12,maxnft),tranm(12,maxnft),roffm(12,maxnft),avyield
      REAL*8 photm(12,maxnft),avmnpet(maxnft),avmnt,avmnppt,tc,ts,tsi
      REAL*8 bulk,nupc,nfix,daygpp,evap,tran,roff,pet,yrtran,yrevap
      REAL*8 yrroff,daily_out(25,maxnft,12,31),ans(12,31),interc,evbs
      REAL*8 soil_chr(10),leafold,stemold,rootold,ftagh(maxnft)
      REAL*8 cluse(maxnft,maxyrs),soil_chr2(10),fprob,ftmix(maxnft)
      REAL*8 oldlai,qdirect,qdiff,sumadp(360,maxnft),petm(12,maxnft)
      REAL*8 sm_trig(30,maxnft),smtrig(30),topsl,suma(360),dayra
      REAL*8 tsumam(maxnft),stemfr(maxnft),lmor_sc(3600,maxnft)
      REAL*8 wilt,field,sat,orgc,lflitold,fpr,gsm(maxnft)

      INTEGER sites,cycle,yr0,yrf,snp_no,snpshts(1000),dschill(maxnft)
      INTEGER ftbbm(maxnft),ftssm(maxnft),ftsss(maxnft),n_fields
      INTEGER snp_year,ftlls(maxnft),ftsls(maxnft),ftrls(maxnft),day,d
      INTEGER isite,ntags,du,ii,otagsn(50),otagsnft(50),sit_grd
      INTEGER ftmor(maxnft),ftc3(maxnft),nft,site0,sitef,nat_map(8)
      INTEGER ilanduse,siteno,iofn,iofnft,iofngft,recl1,nomdos
      INTEGER icontinuouslanduse,ftphen(maxnft),ftdth(maxnft),kode
      INTEGER i,j,k,l,ft,blank,site,year,age,covind,bioind,xyearf
      INTEGER mnth,no_days,fireres,xyear0,per,omav(50),year_out
      INTEGER dolydo(maxnft),luse(maxyrs),fno,bb(maxnft),bbgs(maxnft)
      INTEGER iyear,oymd,oymdft,dsbb(maxnft),chill(maxnft),persum
      INTEGER stcmp,iargc,yearind(maxyrs),idum,outyears,thty_dys
      INTEGER xlatresn,xlonresn,day_mnth,yearv(maxyrs),nyears,narg
      INTEGER seed1,seed2,seed3,spinl,yr0s,yr0p,yrfp,xseed1,site_dat
      INTEGER ibox,jbox,last_blank,site_out,country_id,outyears1
      INTEGER outyears2,budo(maxnft),seno(maxnft),ss(maxnft),clim_type
      INTEGER check_ft_grow,no_countries,n_param_0,n_param_f,n_param
      REAL*8 xtmpv(500,12,31),xprcv(500,12,31),xhumv(500,12,31)
      REAL*8 xcldv(500,12)

      CHARACTER st1*1000,st2*1000,st3*1000,st4*1000,fttags(100)*1000
      CHARACTER st5*1000,otags(100)*1000,ofmt(100)*100,in2st*10000
      CHARACTER stinput*1000,stoutput*1000,stinit*1000,stco2*1000
      CHARACTER stmask*1000,country_name*1000,countries(100)*20
      CHARACTER sttxdp*1000,stlu*1000,ststats*1000,buff1*80
      CHARACTER param_file*1000,date*8,time*10

      LOGICAL initise,initiseo,speedc,crand,xspeedc,withcloudcover
      LOGICAL l_clim,l_lu,l_soil(20),l_stats,l_regional,l_countries
      LOGICAL out_cov,out_bio,out_bud,out_sen,l_b_and_c,check_c
      LOGICAL land_check,l_parameter

*----------------------------------------------------------------------*
      REAL*8 zs1(maxnft),zs2(maxnft),zs3(maxnft),zs4(maxnft)
      REAL*8 zsn(maxnft),zlsn(maxnft),zic0(8),zin0(8),ziminn(3)
      REAL*8 zdslc,zdrlc,zdsln,zdrln,znpp(maxnft),zlai(maxnft)
      REAL*8 zlat,zlon,zcov(maxage,maxnft),zbio(maxage,2,maxnft)
      REAL*8 zppm(maxage,maxnft),zhgt(maxage,maxnft)
      REAL*8 zwdt(maxage,maxnft),zleafdp(3600,maxnft)
      REAL*8 zstemdp(360,maxnft),zrootdp(360,maxnft)
      REAL*8 znps(maxnft),znpl(maxnft),zevp(maxnft),zdof(maxnft)
      REAL*8 zslc(maxnft),zrlc(maxnft),zsln(maxnft),zrln(maxnft)
      REAL*8 znppstore(maxnft),znppstorx(maxnft),znppstor2(maxnft)
      REAL*8 zmaxlai(maxnft),ztmem(200),zsm_trig(30,maxnft)
      REAL*8 zsumadp(360,maxnft),zstemfr(maxnft)
      INTEGER zbb(maxnft),zbbgs(maxnft),zdsbb(maxnft)
*----------------------------------------------------------------------*

      INCLUDE 'param.inc'

      IF (IARGC().GT.0) THEN
        CALL GETARG(1,buff1)
        narg = 1
      ELSE
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) ' Input file must be given as an argument.'
        STOP
      ENDIF

*      narg = 1
*      buff1='..\..\..\output3\geertje\spin.dat'
*      buff1='..\..\vernon\input.dat'

      no_countries = 0
      zlat = -1000.0d0
      zlon = -1000.0d0

*----------------------------------------------------------------------*
* Read input common parameters.                                        *
*----------------------------------------------------------------------*
      OPEN(98,FILE='param.dat',STATUS='OLD',iostat=kode)
      IF (kode.NE.0) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) ' File does not exist: "param.dat"'
        STOP
      ENDIF

      READ(98,*)
      READ(98,*) st1
      CALL STRIPB(st1)
      st2 = stver
      CALL STRIPBS(st2,st3)
      CALL STRIPB(st2)
      st1=st1(1:blank(st1))
      st2=st2(1:blank(st2))
      IF (stcmp(st1,st2).NE.1) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'Version number mismatch between ''sdgvm0.f'' and ''p
     &aram.dat''.'
        STOP
      ENDIF

      READ(98,*)
      READ(98,*) p_sla        ! 4
      READ(98,*)
      READ(98,*) p_stemfr     ! 6
      READ(98,*)
      READ(98,*) p_rootfr     ! 8
      READ(98,*)
      READ(98,*) p_opt        !10
      READ(98,*)
      READ(98,*) p_stmin      !12
      READ(98,*)
      READ(98,*) p_laimem     !14
      READ(98,*)
      READ(98,*) p_resp       !16
      READ(98,*)
      READ(98,*) p_kscale     !18
      READ(98,*)
      READ(98,*) p_nu1,p_nu2,p_nu3,p_nu4
      READ(98,*)
      READ(98,*) p_nleaf
      READ(98,*)
      READ(98,*) p_dresp
      READ(98,*)
      READ(98,*) p_vm
      READ(98,*)
      READ(98,*) p_kgw
      READ(98,*)
      READ(98,*) p_v1,p_v2,p_v3
      READ(98,*)
      READ(98,*) p_j1,p_j2,p_j3
      READ(98,*)
      READ(98,*) p_pet
      READ(98,*)
      READ(98,*) p_bs
      READ(98,*)
      READ(98,*) p_et
      READ(98,*)
      READ(98,*) p_roff
      READ(98,*)
      READ(98,*) p_roff2
      READ(98,*)
      READ(98,*) p_fprob
      READ(98,*)
      READ(98,*) p_city_dep
      CLOSE(98)

*----------------------------------------------------------------------*
* Read 'misc_params.dat'.                                              *
*----------------------------------------------------------------------*
      OPEN(98,FILE='misc_params.dat',STATUS='OLD',iostat=kode)
      IF (kode.NE.0) THEN
        WRITE(*,*) ' File does not exist: "misc_params.dat"'
        WRITE(*,*) ' Using screen output options: 1 0. '
        WRITE(*,*) ' Using countries, not regions. '
        site_out = 1
        year_out = 0
        l_regional = .FALSE.
      ELSE
        READ(98,*)
        READ(98,*) site_out,year_out
        READ(98,*)
        READ(98,*) l_regional
      ENDIF
      CLOSE(98)

*----------------------------------------------------------------------*
* Read input data.                                                     *
*----------------------------------------------------------------------*
      OPEN(98,FILE=buff1,STATUS='OLD',iostat=kode)
      IF (kode.NE.0) THEN
        DO i=1,80
          st1(i:i) = buff1(i:i)
        ENDDO
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'SDGVM input file does not exist.'
        WRITE(*,'('' "'',A,''"'')') st1(1:blank(st1))
        STOP
      ENDIF

      READ(98,'(A)') st1
      st2 = 'DOS'
      st3 = 'UNIX'

      IF (stcmp(st1,st2).EQ.1) THEN
        du = 1
      ELSEIF (stcmp(st1,st3).EQ.1) THEN
        du = 0
      ELSE
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'First line of input file must read DOS or UNIX.'
        STOP
      ENDIF

      READ(98,'(A)') stinput
      CALL STRIPB(stinput)
      st1 = stinput
      OPEN(99,FILE=st1(1:blank(st1))//'/readme.dat',STATUS='OLD',
     &iostat=kode)
      IF (kode.NE.0) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'Climate data file does not exist:'
        WRITE(*,'('' "'',A,''/readme.dat"'')') st1(1:blank(st1))
        STOP
      ENDIF
      READ(99,'(A)') st1

*----------------------------------------------------------------------*
* Determine whether climate is daily or monthly, from first line of    *
* readme file.                                                         *
*----------------------------------------------------------------------*
      st2 = 'DAILY'
      st3 = 'MONTHLY'
      st4 = 'SITED'
      st5 = 'SITEM'
      sit_grd = 0
      clim_type = 0
      IF (stcmp(st1,st2).EQ.1) THEN
        clim_type = 1
        day_mnth = 1
        thty_dys = 1
        READ(99,*)
        READ(99,*) xlatf,xlon0
        READ(99,*)
        READ(99,*) xlatres,xlonres
        READ(99,*)
        READ(99,*) xlatresn,xlonresn
        READ(99,*)
        READ(99,*) xyear0,xyearf
      ELSEIF (stcmp(st1,st3).EQ.1) THEN
        clim_type = 2
        day_mnth = 0
        thty_dys = 1
        READ(99,*)
        READ(99,*) xlatf,xlon0
        READ(99,*)
        READ(99,*) xlatres,xlonres
        READ(99,*)
        READ(99,*) xlatresn,xlonresn
        READ(99,*)
        READ(99,*) xyear0,xyearf
      ELSEIF (stcmp(st1,st4).EQ.1) THEN
        clim_type = 3
        day_mnth = 1
        thty_dys = 0
        sit_grd = 1
        READ(99,*)
        READ(99,*) xlatf,xlon0
        READ(99,*)
        READ(99,*) xyear0,xyearf
      ELSEIF (stcmp(st1,st5).EQ.1) THEN
        clim_type = 4
        day_mnth = 0
        thty_dys = 1
        sit_grd = 1
        READ(99,*)
        READ(99,*) xlatf,xlon0
        READ(99,*)
        READ(99,*) xyear0,xyearf
      ELSE
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*)
     &'First line of climat readme.dat file must read DAILY MONTHLY or 
     &SITE.'
        STOP
      ENDIF
      CLOSE(99)
*----------------------------------------------------------------------*

      st1 = stinput

      READ(98,'(A)') ststats
      CALL STRIPB(ststats)

      READ(98,'(A)') sttxdp
      CALL STRIPB(sttxdp)

      READ(98,'(A)') stlu
      CALL STRIPB(stlu)

      READ(98,'(A)') stco2
      CALL STRIPB(stco2)

      READ(98,*) co2const

      READ(98,'(A)') stmask
      CALL STRIPB(stmask)

      READ(98,*)

      READ(98,'(A)') stinit
      CALL STRIPB(stinit)
      st2 = 'ARGUMENT'
      IF (stcmp(stinit,st2).EQ.1) THEN
        narg = narg + 1
        CALL GETARG(narg,buff1)
        DO i=1,80
          stinit(i:i) = buff1(i:i)
        ENDDO
      ENDIF

      READ(98,'(A)') stoutput
      CALL STRIPB(stoutput)
      st2 = 'ARGUMENT'
      IF (stcmp(stoutput,st2).EQ.1) THEN
        narg = narg + 1
        CALL GETARG(narg,buff1)
        DO i=1,80
          stoutput(i:i) = buff1(i:i)
        ENDDO
      ENDIF
      READ(98,*)

*----------------------------------------------------------------------*
* Open diagnostics file.                                               *
*----------------------------------------------------------------------*
      st1 = stoutput
      OPEN(11,FILE=st1(1:blank(st1))//'/diag.dat',iostat=kode)
      IF (kode.NE.0) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'SDGVM output directory does not exist.'
        WRITE(*,'('' "'',A,''"'')') st1(1:blank(st1))
        STOP
      ENDIF
      WRITE(11,'(''ERRORS'')')
*----------------------------------------------------------------------*

c change by Ghislain
      READ(98,'(1000a)') st1
      ii = n_fields(st1)
      IF ((ii.EQ.3).OR.(ii.EQ.4)) THEN
        CALL STRIPBN(st1,i)
        initise = .true.
        IF (i.EQ.0)  initise = .false.
        CALL STRIPBN(st1,j)
        initiseo = .true.
        IF (j.EQ.0)  initiseo = .false.
        CALL STRIPBN(st1,k)
        speedc = .true.
        IF (k.EQ.0)  speedc = .false.
        IF (ii.EQ.4) THEN
          CALL STRIPBN(st1,xseed1)
        ELSE
          xseed1 = 1
        ENDIF
      ELSE
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'Line 13 must contain 3 or 4 arguments'
        WRITE(*,'('' "'',A,''"'')') st1(1:30)
        STOP
      ENDIF
      xspeedc = speedc

      READ(98,'(1000a)') st1
      CALL STRIPB(st1)
      i = n_fields(st1)
      IF (i.EQ.7) THEN
        CALL STRIPBN(st1,spinl)
        CALL STRIPBN(st1,yr0s)
        CALL STRIPBN(st1,cycle)
        CALL STRIPBN(st1,j)
        crand = .TRUE.
        IF (j.EQ.0)  crand = .FALSE.
        CALL STRIPBN(st1,yr0p)
        CALL STRIPBN(st1,yrfp)
        CALL STRIPBN(st1,outyears)
      ELSEIF (i.EQ.6) THEN
        CALL STRIPBN(st1,spinl)
        CALL STRIPBN(st1,yr0s)
        CALL STRIPBN(st1,cycle)
        CALL STRIPBN(st1,j)
        crand = .TRUE.
        IF (j.EQ.0)  crand = .FALSE.
        CALL STRIPBN(st1,yr0p)
        CALL STRIPBN(st1,yrfp)
        outyears = yrfp - yr0p + 1
      ELSEIF (i.EQ.5) THEN
        CALL STRIPBN(st1,spinl)
        CALL STRIPBN(st1,yr0s)
        CALL STRIPBN(st1,cycle)
        CALL STRIPBN(st1,j)
        crand = .TRUE.
        IF (j.EQ.0)  crand = .FALSE.
        yr0p = yr0s+1
        yrfp = yr0s
        CALL STRIPBN(st1,outyears)
      ELSEIF (i.EQ.4) THEN
        CALL STRIPBN(st1,spinl)
        CALL STRIPBN(st1,yr0s)
        CALL STRIPBN(st1,cycle)
        CALL STRIPBN(st1,j)
        crand = .TRUE.
        IF (j.EQ.0)  crand = .FALSE.
        yr0p = yr0s+1
        yrfp = yr0s
        outyears = cycle + 1
      ELSEIF (i.EQ.3) THEN
        spinl = 0
        cycle = maxyrs
        j = 0
        crand = .FALSE.
        CALL STRIPBN(st1,yr0p)
        CALL STRIPBN(st1,yrfp)
        CALL STRIPBN(st1,outyears)
        yr0s = yr0p
      ELSEIF (i.EQ.2) THEN
        spinl = 0
        cycle = maxyrs
        j = 0
        crand = .FALSE.
        CALL STRIPBN(st1,yr0p)
        CALL STRIPBN(st1,yrfp)
        yr0s = yr0p
        outyears = yrfp - yr0p + 1
      ELSE
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'Line 14 must contain 2-7 fields'
        STOP
      ENDIF
      nyears = spinl+yrfp-yr0p+1
      IF (nyears.GT.maxyrs) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,'('' Trying to simulate '',i4,'' years, maximum allowabl
     &e is '',i4,''.'')') nyears,maxyrs
        WRITE(*,*) 'Either reduce the length of the simulation, or incr
     &ease "maxyrs"'
        WRITE(*,*) 'this is set in array_param.inc, you must re-comile 
     &after altering'
      WRITE(*,*) 'this file.'
        STOP
      ENDIF
      IF ((j.NE.0).AND.(j.NE.1)) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,'('' Fourth field of line 14 in the input file must be e
     &ither 0 or 1.'')')
        WRITE(*,'('' Currently set to '',i5,''.'')') j
        STOP
      ENDIF
      outyears = min(outyears,nyears)

*----------------------------------------------------------------------*
* Set 'yr0' and 'yrf' which are the years of actual climate required   *
* for the run. And check that the climate exists in the climate        *
* database                                                             *
*----------------------------------------------------------------------*
      yr0 = min(yr0s,yr0p)
      yrf = max(yr0s+min(spinl,cycle)-1,yrfp)
      IF ((yr0.LT.xyear0).OR.(yrf.GT.xyearf)) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,'('' Tyring to use '',i4,''-'',i4,'' climate.'')') 
     &yr0,yrf 
        WRITE(*,'('' Climate database runs from '',i4,''-'',i4,''.'')') 
     &xyear0,xyearf
        STOP
      ENDIF

*----------------------------------------------------------------------*
* Set up year climate sequence.                                        *
*----------------------------------------------------------------------*
      DO i=1,maxyrs
        yearind(i) = i
      ENDDO
      idum = 1
      DO i=1,nyears
        IF (i.LE.spinl) THEN
          IF (crand) THEN
            IF ((mod(i-1,cycle).EQ.0).AND.(i.GT.2))
     & CALL RANDOMV(yearind,1,cycle,idum)
            yearv(i) = yearind(mod(i-1,cycle)+1) + yr0s - 1
          ELSE
            yearv(i) = mod(i-1,cycle) + yr0s
          ENDIF
        ELSE
          yearv(i) = i - spinl + yr0p - 1
        ENDIF
      ENDDO

*----------------------------------------------------------------------*
* Open monthly/daily PIXEL output files if required.                   *
*----------------------------------------------------------------------*
      CALL OUTPUT_OPTIONS(nomdos,otags,omav,ofmt)

      READ(98,'(A)') st1
      CALL STRIPB(st1)
      st2 = 'PIXEL'
      IF (stcmp(st1,st2).EQ.0) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'The first field of line 15 of the input file must
     & read ''PIXEL''.'
        WRITE(*,*) st1(1:30)
        WRITE(*,*) 'Output variable options:'
        WRITE(*,'(1x,20a4)') (otags(i),i=1,15)
        WRITE(*,'(1x,20a4)') (otags(i),i=16,nomdos)
        STOP
      ENDIF

      CALL SET_PIXEL_OUT(st1,st2,outyears1,nomdos,otagsn,otags,
     &oymd)
      outyears1 = min(outyears1,nyears)

      st1 = stoutput
      iofn = 400
      IF (outyears1.GT.0) THEN
      IF ((oymd.EQ.1).OR.(oymd.EQ.3)) THEN
        DO i=1,50
          IF (otagsn(i).EQ.1) THEN
            iofn = iofn + 1
            OPEN(iofn,file=st1(1:blank(st1))//'/monthly_'//otags(i)(1:3)
     &//'.dat')
          ENDIF
        ENDDO
      ENDIF
      IF ((oymd.EQ.2).OR.(oymd.EQ.3)) THEN
        DO i=1,50
          IF (otagsn(i).EQ.1) THEN
            iofn = iofn + 1
            OPEN(iofn,file=st1(1:blank(st1))//'/daily_'//otags(i)(1:3)
     &//'.dat')
          ENDIF
        ENDDO
      ENDIF
      ENDIF
      iofnft = iofn

*----------------------------------------------------------------------*
* Determine whether daily or monthly subpixel outputs are required.    *
*----------------------------------------------------------------------*

      READ(98,'(A)') st1
      CALL STRIPB(st1)
      st2 = 'SUBPIXEL'
      IF (stcmp(st1,st2).EQ.0) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'The first field of line 16 of the input file must
     & read ''SUBPIXEL''.'
        WRITE(*,'('' "'',A,''"'')') st1(1:30)
        WRITE(*,*) 'Output variable options:'
        WRITE(*,'(1x,20a4)') (otags(i),i=1,15)
        WRITE(*,'(1x,20a4)') (otags(i),i=16,nomdos),
     &'cov ','bio ','bud ','sen '
        STOP
      ENDIF

      CALL SET_SUBPIXEL_OUT(st1,st2,outyears2,nomdos,otagsnft,
     &otags,oymdft,out_cov,out_bio,out_bud,out_sen)
      outyears2 = min(outyears2,nyears)

*----------------------------------------------------------------------*
* Read in snapshot years.                                              *
*----------------------------------------------------------------------*
      READ(98,'(A)') st1
      CALL STRIPB(st1)
      st2 = 'SNAPSHOTS'
      IF (stcmp(st1,st2).EQ.0) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'First field of line 17 must read SNAPSHOTS'
        STOP
      ENDIF
      CALL STRIPBS(st1,st2)
      CALL STRIPB(st1)
      CALL ST2ARR(st1,snpshts,100,snp_no)
      READ (98,*)

*----------------------------------------------------------------------*
* Read in compulsory functional types.                                 *
*----------------------------------------------------------------------*
      fttags(1) = 'BARE'
      fttags(2) = 'CITY'
      READ(98,'(A)') st1
      CALL STRIPB(st1)
      st2 = 'BARE'
      IF ((stcmp(st2,st1).EQ.0).OR.(n_fields(st1).NE.2)) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'Line 19 must read "BARE" forllowed by a number (0-1)
     &.'
        STOP
      ENDIF
      CALL STRIPBS(st1,st2)
      READ(st1,*) ftmix(1)
      READ(98,'(A)') st1
      CALL STRIPB(st1)
      st2 = 'CITY'
      IF ((stcmp(st2,st1).EQ.0).OR.(n_fields(st1).NE.2)) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'Line 19 must read "BARE" forllowed by a number (0-1)
     &.'
        STOP
      ENDIF
      CALL STRIPBS(st1,st2)
      READ(st1,*) ftmix(2)
      READ (98,*)

* Initialise redundant parameterisation 
      ft = 1
      ftc3(ft) = 0
      ftphen(ft) = 0
      ftagh(ft) = 0.0d0
      ftdth(ft) = 0
      ftmor(ft) = 0
      ftwd(ft) = 0.0d0
      ftxyl(ft) = 0.0d0
      ftpd(ft) = 0.0d0
      ftsla(ft) = 0.0d0
      ftlls(ft) = 0
      ftsls(ft) = 0
      ftrls(ft) = 0
      ftlmor(ft) = 0.0d0
      ftrat(ft) = 0.0d0
      ftbbm(ft) = 0
      chill(ft) = 0
      dschill(ft) = 0
      ftbb0(ft) = 0.0d0
      ftbbmax(ft) = 0.0d0
      ftbblim(ft) = 0.0d0
      ftssm(ft) = 0
      ftsss(ft) = 0
      ftsslim(ft) = 0.0d0
      ftstmx(ft) = 0.0d0
      ftgr0(ft) = 0.0d0
      ftgrf(ft) = 0.0d0
      ftppm0(ft) = 0.0d0

      ft = 2
      ftc3(ft) = 0
      ftphen(ft) = 0
      ftagh(ft) = 0.0d0
      ftdth(ft) = 0
      ftmor(ft) = 0
      ftwd(ft) = 0.0d0
      ftxyl(ft) = 0.0d0
      ftpd(ft) = 0.0d0
      ftsla(ft) = 0.0d0
      ftlls(ft) = 0
      ftsls(ft) = 0
      ftrls(ft) = 0
      ftlmor(ft) = 0.0d0
      ftrat(ft) = 0.0d0
      ftbbm(ft) = 0
      chill(ft) = 0
      dschill(ft) = 0
      ftbb0(ft) = 0.0d0
      ftbbmax(ft) = 0.0d0
      ftbblim(ft) = 0.0d0
      ftssm(ft) = 0
      ftsss(ft) = 0
      ftsslim(ft) = 0.0d0
      ftstmx(ft) = 0.0d0
      ftgr0(ft) = 0.0d0
      ftgrf(ft) = 0.0d0
      ftppm0(ft) = 0.0d0

*----------------------------------------------------------------------*
* Read in functional type parameterisation.                            *
*----------------------------------------------------------------------*

      READ(98,'(A)') st1
      CALL STRIPB(st1)

      IF (n_fields(st1).EQ.1) THEN
        st2 = 'ARGUMENT'

        IF (stcmp(st2,st1).EQ.1) THEN
*----------------------------------------------------------------------*
* From a file.                                                         *
*----------------------------------------------------------------------*
          narg = narg + 1
          CALL GETARG(narg,buff1)
          DO i=1,80
            st1(i:i) = buff1(i:i)
          ENDDO
        ENDIF

        OPEN(99,FILE=st1,STATUS='OLD',iostat=kode)
        IF (kode.NE.0) THEN
          WRITE(*,'('' PROGRAM TERMINATED'')')
          WRITE(*,*) 'Functional type parameterisation file does not ex
     &ist.'
          WRITE(*,'('' "'',A,''"'')') st1(1:blank(st1))
          STOP
        ENDIF

        ft = 2

40      CONTINUE
          READ(99,'(1000a)',end=30) st1
          ft = ft + 1

            IF (ft.GT.maxnft) THEN
              WRITE(*,'('' PROGRAM TERMINATED'')')
              WRITE(*,'('' Maximum number of fts is'',i4,'', currently m
     &ore than this are defined.'')') maxnft
              WRITE(*,*) 'Either increase "maxnft" defined in "array_dim
     &s.dat" or decrease'
              WRITE(*,*) 'the number of ft parameterisations: this requi
     &res re-compilation.'
              STOP
            ENDIF

          IF (n_fields(st1).NE.27) THEN
            WRITE(*,'('' PROGRAM TERMINATED'')')
            WRITE(*,*) 'The ft parameterisation must contain 27 fields.
     &'
            WRITE(*,'(1x,A,'' has '',i3)') st1(1:blank(st1)),n_fields(st
     &1)
            STOP
          ENDIF

          READ(st1,*) fttags(ft),ftmix(ft),ftc3(ft),ftphen(ft),
     &ftagh(ft),ftdth(ft),ftmor(ft),ftwd(ft),ftxyl(ft),ftpd(ft),
     &ftsla(ft),ftlls(ft),ftsls(ft),ftrls(ft),ftlmor(ft),ftrat(ft),
     &ftbbm(ft),ftbb0(ft),ftbbmax(ft),ftbblim(ft),ftssm(ft),ftsss(ft),
     &ftsslim(ft),ftstmx(ft),ftgr0(ft),ftgrf(ft),ftppm0(ft)

          IF (ftmor(ft).GT.maxyrs) THEN
            WRITE(*,'('' PROGRAM TERMINATED'')')
            WRITE(*,'('' Mortality of '',A,'' is '',i4,'', maximum allow
     &able is '',i4,''.'')')
            STOP
          ENDIF

          IF (ftsla(ft).LT.0.0d0) THEN
            ftsla(ft) = 10.0d0**(2.35d0 -
     &0.39d0*log10(real(ftlls(ft))/30.0d0))*2.0d0/10000.0d0
*      ftsla(ft) = 10.0d0**(2.43d0-0.46d0*log10(real(ftlls(ft))/30.0d0))
*     &*2.0/10000.0
          ENDIF

          ftsla(ft) = ftsla(ft)/p_sla

          IF (ftlls(ft).LT.0.0d0) THEN
            ftlls(ft) = int(10.0d0**((2.35d0 -
     &log10(ftsla(ft)*10000.0d0/2.0d0))/0.39d0)*30.0d0+0.5d0)
          ENDIF

          IF (ftmor(ft).GT.maxage) THEN
            WRITE(*,'('' PROGRAM TERMINATED'')')
            WRITE(*,'('' Maximum age of "'',A,''" is '',i4,'', maximum a
     &llowable age is '',i4,''.'')') st1(1:blank(st1)),ftmor(ft),maxage
            WRITE(*,*) 'Either increase "maxage" defined in "array_dims"
     &'
            WRITE(*,*) 'or decrease mortality in the parameterisation.'
            STOP
          ENDIF

        GOTO 40
30      CONTINUE
        nft = ft
        READ (98,*)

      ELSE
*----------------------------------------------------------------------*
* From the input file.                                                 *
*----------------------------------------------------------------------*
        ft = 2
98      CONTINUE

          IF (ichar(st1(1:1)).NE.32) THEN
            ft = ft + 1

            IF (ft.GT.maxnft) THEN
              WRITE(*,'('' PROGRAM TERMINATED'')')
              WRITE(*,'('' Maximum number of fts is'',i4,'', currently m
     &ore than this are defined.'')') maxnft
              WRITE(*,*) 'Either increase "maxnft" defined in "array_dim
     &s.dat" or decrease'
              WRITE(*,*) 'the number of ft parameterisations: this requi
     &res re-compilation.'
              STOP
            ENDIF

            IF (n_fields(st1).NE.27) THEN
              WRITE(*,'('' PROGRAM TERMINATED'')')
              WRITE(*,*) 'The ft parameterisation must contain 27 field
     &s.'
              WRITE(*,'(1x,A,'' has '',i3)') st1(1:blank(st1)),n_fields(
     &st1)
              STOP
            ENDIF

            READ(st1,*) fttags(ft),ftmix(ft),ftc3(ft),ftphen(ft),
     &ftagh(ft),ftdth(ft),ftmor(ft),ftwd(ft),ftxyl(ft),ftpd(ft),
     &ftsla(ft),ftlls(ft),ftsls(ft),ftrls(ft),ftlmor(ft),ftrat(ft),
     &ftbbm(ft),ftbb0(ft),ftbbmax(ft),ftbblim(ft),ftssm(ft),ftsss(ft),
     &ftsslim(ft),ftstmx(ft),ftgr0(ft),ftgrf(ft),ftppm0(ft)

            IF (ftsla(ft).LT.0.0d0) THEN
              ftsla(ft) = 10.0d0**(2.35d0 -
     &0.39d0*log10(real(ftlls(ft))/30.0d0))*2.0d0/10000.0d0
*      ftsla(ft) = 10.0d0**(2.43d0-0.46d0*log10(real(ftlls(ft))/30.0d0))
*     &*2.0/10000.0
            ENDIF

            ftsla(ft) = ftsla(ft)/p_sla

            IF (ftlls(ft).LT.0.0d0) THEN
              ftlls(ft) = int(10.0d0**((2.35d0 -
     &log10(ftsla(ft)*10000.0d0/2.0d0))/0.39d0)*30.0d0+0.5d0)
            ENDIF

            IF (ftmor(ft).GT.maxage) THEN
              WRITE(*,'('' PROGRAM TERMINATED'')')
              WRITE(*,'('' Maximum age of "'',A,''" is '',i4,'', maximum
     & allowable age is '',i4,''.'')') st1(1:blank(st1)),ftmor(ft),maxag
     &e
              WRITE(*,*) 'Either increase "maxage" defined in "array_dim
     &s.dat"'
              WRITE(*,*) 'or decrease mortality in the parameterisation
     &.'
              STOP
            ENDIF

            READ(98,'(1000a)') st1
            CALL STRIPB(st1)
          GOTO 98        
        ENDIF
        nft = ft
      ENDIF

*----------------------------------------------------------------------*
* Create leaf mortality scales values.                                 *
*----------------------------------------------------------------------*
      lmor_sc(1,1)=0.0d0
      lmor_sc(2,1)=0.0d0
      DO ft=3,nft
        DO i=1,ftlls(ft)
          lmor_sc(i,ft)=(real(ftlls(ft)-i)/real(ftlls(ft)))**ftlmor(ft)
        ENDDO
      ENDDO

*----------------------------------------------------------------------*
* Open site_info file and check output directory exists.               *
*----------------------------------------------------------------------*
      st1 = stoutput
      OPEN(12,FILE=st1(1:blank(st1))//'/site_info.dat',iostat=kode)
      IF (kode.NE.0) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'Output directory does not exist:'
        WRITE(*,'('' "'',A,''"'')') st1(1:blank(st1))
        STOP
      ENDIF
      WRITE(12,'(''COUNTRY           ID     lat      lon     co20  co2f 
     &   sand   silt   bulk   orgc     wp     fc    swc    dep       ''
     &,50A13)') (fttags(i),i=1,nft)

*----------------------------------------------------------------------*
* Open monthly/daily ft output files.                                  *
*----------------------------------------------------------------------*
      IF (outyears2.GT.0) THEN
      iofn = iofnft
      st1 = stoutput
      IF ((oymdft.EQ.1).OR.(oymdft.EQ.3)) THEN
        DO i=1,50
          IF (otagsnft(i).EQ.1) THEN
            DO ft=1,nft
              iofn = iofn + 1
              OPEN(iofn,file=st1(1:blank(st1))//'/monthly_'//otags(i)(1:
     &3)//'_'//fttags(ft)(1:blank(fttags(ft)))//'.dat')
           ENDDO
          ENDIF
        ENDDO
      ENDIF
      IF ((oymdft.EQ.2).OR.(oymdft.EQ.3)) THEN
        DO i=1,50
          IF (otagsnft(i).EQ.1) THEN
            DO ft=1,nft
              iofn = iofn + 1
              OPEN(iofn,file=st1(1:blank(st1))//'/daily_'//otags(i)(1:3)
     &//'_'//fttags(ft)(1:blank(fttags(ft)))//'.dat')
            ENDDO
          ENDIF
        ENDDO
      ENDIF
      ENDIF

      DO i=1,nft
        chill(i) = 0
        dschill(i) = 0
      ENDDO

*----------------------------------------------------------------------*
* Read in landuse index mapping.                                       *
*----------------------------------------------------------------------*
      DO j=1,255
        DO k=1,nft
          lutab(j,k) = 0.0d0
        ENDDO
      ENDDO

      ii = 0
96    CONTINUE
        READ(98,'(1000a)') st1
        ii = ii + 1
        CALL STRIPB(st1)
c     read table of conversion from class to ft's proportion
        IF (ichar(st1(1:1)).NE.32) THEN
          i = n_fields(st1)
          CALL STRIPBN(st1,j)
          persum = 0
          DO k=1,(i-1)/2
            CALL STRIPBN(st1,per)
            persum = persum + per
            CALL STRIPBS(st1,st2)
            l = ntags(fttags,st2)
            IF (l.EQ.-9999) THEN
              WRITE(*,'('' PROGRAM TERMINATED'')')
              WRITE(*,*) 'Error in a tag name in the land use mapping.'
              WRITE(*,'('' Line number'',i3,'', within the mapping list.
     &'')') ii
              WRITE(*,'('' Category number'',i3,'', tag number'',i2,''.'
     &')') j,k
              WRITE(*,'('' "'',A,''"'')') st2(1:blank(st2))
              STOP
            ENDIF
            lutab(j,l) = real(per)
          ENDDO
          lutab(j,1) = 100.0d0 - real(persum)

        GOTO 96        
      ENDIF
*----------------------------------------------------------------------*

      READ(98,*) grassrc, barerc, fireres
      READ(98,*)

*----------------------------------------------------------------------*
* Read in sites.                                                       *
*----------------------------------------------------------------------*
      land_check = .true.

      latdel = 0.0d0
      londel = 0.0d0

      IF (sit_grd.EQ.1) THEN
*----------------------------------------------------------------------*
* Single site defined in climate file.                                 *
*----------------------------------------------------------------------*
        READ(98,'(1000a)') st1
        CALL STRIPB(st1)

        IF (n_fields(st1).EQ.1) THEN
          READ(st1,*) sites
        ELSE
          sites = 1
        ENDIF
        DO site=1,sites
          lat_lon(site,1) = xlatf
          lat_lon(site,2) = xlon0
        ENDDO
        READ(98,*)
      ELSE

        READ(98,'(1000a)') st1
        CALL STRIPB(st1)

        IF (n_fields(st1).EQ.1) THEN
          land_check = .false.

          IF (du.eq.1) THEN
            recl1 = 16
          ELSE
            recl1 = 17
          ENDIF

          st2 = 'ARGUMENT'
          IF (stcmp(st2,st1).EQ.1) THEN

            narg = narg + 1
            CALL GETARG(narg,buff1)
            DO i=1,80
              st1(i:i) = buff1(i:i)
            ENDDO

            narg = narg + 1
            CALL GETARG(narg,buff1)
            DO i=1,80
              st2(i:i) = buff1(i:i)
            ENDDO
            CALL STRIPBN(st2,site0)

            narg = narg + 1
            CALL GETARG(narg,buff1)
            DO i=1,80
              st2(i:i) = buff1(i:i)
            ENDDO
            CALL STRIPBN(st2,sitef)

            OPEN(99,FILE=st1,STATUS='OLD',FORM='FORMATTED',
     &ACCESS='DIRECT',RECL=recl1,iostat=kode)
            IF (kode.NE.0) THEN
              WRITE(*,'('' PROGRAM TERMINATED'')')
              WRITE(*,*) 'File containing list of sites does not exist.
     &'
              WRITE(*,'('' "'',A,''"'')') st1(1:blank(st1))
              STOP
            ENDIF

            DO site=site0,sitef
              READ(99,'(F7.3,F9.3)',REC=site) 
     &lat_lon(site-site0+1,1),lat_lon(site-site0+1,2)
            ENDDO
            sites = sitef - site0 + 1

c CLOSE added by Ghislain 15/12/03
            CLOSE(99)

            READ(98,*)

          ELSE

            OPEN(99,FILE=st1,STATUS='OLD',iostat=kode)
            IF (kode.NE.0) THEN
              WRITE(*,'('' PROGRAM TERMINATED'')')
              WRITE(*,*) 'File containing list of sites does not exist.
     &'
              WRITE(*,'('' "'',A,''"'')') st1(1:blank(st1))
              STOP
            ENDIF

* Check for countries to be run.
            READ(98,'(1000a)') st1
            CALL STRIPB(st1)

            IF (n_fields(st1).GT.0) THEN
              no_countries = n_fields(st1)
              READ(st1,*) (countries(i),i=1,no_countries)
              READ(98,*)
            ENDIF

            sites = 0
50          READ(99,*,END=60) lat_lon(sites+1,1),lat_lon(sites+1,2)
              sites = sites + 1
            GOTO 50
60          CONTINUE

c CLOSE added by Ghislain 15/12/03
            CLOSE(99)
          ENDIF

        ELSEIF ((n_fields(st1).EQ.2).OR.(n_fields(st1).EQ.3)) THEN

*----------------------------------------------------------------------*
* List of sites in the input file.                                     *
*----------------------------------------------------------------------*
          sites = 1
          READ(st1,*) lat_lon(sites,1),lat_lon(sites,2)
97        CONTINUE
            READ(98,'(1000a)') st1
            CALL STRIPB(st1)
            IF (ichar(st1(1:1)).NE.32) THEN
              sites = sites + 1
              READ(st1,*) lat_lon(sites,1),lat_lon(sites,2)
            GOTO 97
          ENDIF

        ELSEIF ((n_fields(st1).EQ.4).OR.(n_fields(st1).EQ.5)) THEN

*----------------------------------------------------------------------*
* Box of sites in input file.                                          *
*----------------------------------------------------------------------*
          READ(st1,*) lat0,latf,lon0,lonf
          IF ((lat0.GT.latf).OR.(lon0.GT.lonf)) THEN
            WRITE(*,'('' PROGRAM TERMINATED'')')
            WRITE(*,*) 'lat0 and lon0 must be less than latf and lonf.'
          ENDIF
          IF (lat0.LT.-90.0d0)  lat0 = -90.0d0
          IF (latf.GT. 90.0d0)  latf =  90.0d0
          IF (lon0.LT.-180.0d0)  lon0 = -180.0d0
          IF (lonf.GT. 180.0d0)  lonf =  180.0d0
          sites = 0
          DO i=1,int((latf-lat0)/xlatres+0.5d0)
            DO j=1,int((lonf-lon0)/xlonres+0.5d0)
              sites = sites + 1
              lat_lon(sites,1) = lat0 + (real(i-1)+0.5d0)*xlatres
              lat_lon(sites,2) = lon0 + (real(j-1)+0.5d0)*xlonres
            ENDDO
          ENDDO
          READ(98,*)
          latdel = xlatres
          londel = xlonres

        ELSEIF (n_fields(st1).GE.6) THEN

*----------------------------------------------------------------------*
* Box of sites in input file.                                          *
*----------------------------------------------------------------------*
          READ(st1,*) lat0,latf,lon0,lonf,latdel,londel
          IF ((lat0.GT.latf).OR.(lon0.GT.lonf)) THEN
            WRITE(*,'('' PROGRAM TERMINATED'')')
            WRITE(*,*) 'lat0 and lon0 must be less than latf and lonf.'
          ENDIF
          IF (lat0.LT.-90.0d0)  lat0 = -90.0d0
          IF (latf.GT. 90.0d0)  latf =  90.0d0
          IF (lon0.LT.-180.0d0)  lon0 = -180.0d0
          IF (lonf.GT. 180.0d0)  lonf =  180.0d0
          IF (latdel.LT.1.0d0/120.0d0) latdel = 1.0d0/120.0d0
          IF (londel.LT.1.0d0/120.0d0) londel = 1.0d0/120.0d0
          sites = 0
          IF (latf.GT.lat0+latdel/2.0d0) THEN
            ibox = int((latf-lat0)/latdel+0.5d0)
            iadj = latdel/2.0d0
          ELSE
            ibox = 1
            iadj = 0.0d0
          ENDIF
          IF (lonf.GT.lon0+londel/2.0d0) THEN
            jbox = int((lonf-lon0)/londel+0.5d0)
            jadj = londel/2.0d0
          ELSE
            jbox = 1
            jadj = 0.0d0
          ENDIF
          DO i=1,ibox
            DO j=1,jbox
              sites = sites + 1
              lat_lon(sites,1) = lat0 + real(i-1)*latdel + iadj
              lat_lon(sites,2) = lon0 + real(j-1)*londel + jadj
            ENDDO
          ENDDO
          READ(98,*)
        ELSE
          WRITE(*,'('' PROGRAM TERMINATED'')')
          WRITE(*,*) 'Error in the number of fields for site input'
          STOP
        ENDIF

      ENDIF

*---------------------------------------------------------------------*
* Check sites against land mask and disregard when no land.           *
*---------------------------------------------------------------------*
      st2 = stoutput
      IF (land_check) THEN
        CALL LAND_SITE_CHECK(st1,st2,sites,lat_lon,latdel,londel,du,
     &stmask,no_countries)
      ENDIF

*----------------------------------------------------------------------*
* Read in soil characteristics. Defaults used when zero                *
*----------------------------------------------------------------------*
      READ(98,'(1000a)') st1
      IF ((n_fields(st1).LT.8).OR.(n_fields(st1).GT.9)) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'Eight or nine fields must be present for the soil 
     &characteristics: sand %,'
        WRITE(*,*) 'silt %, bulk density and Soil Depth (cm). Set to 
     &zero if default values'
        WRITE(*,*) 'are required.'
        STOP
      ENDIF
      IF (n_fields(st1).EQ.8) THEN
        READ(st1,*) (soil_chr(i),i=1,8)
        topsl = defaulttopsl
      ELSE
        READ(st1,*) (soil_chr(i),i=1,8),topsl
        IF (topsl.LT.1.0e-6) topsl = defaulttopsl
      ENDIF
      IF ((soil_chr(5).LT.0.0d0).OR.(soil_chr(6).LT.0.0d0).OR.
     &(soil_chr(7).LT.0.0d0)) THEN
        l_b_and_c = .TRUE.
      ELSE
        l_b_and_c = .FALSE.
      ENDIF
      READ(98,'(1000a)') st1

      READ(98,'(1000a)') st1

*----------------------------------------------------------------------*
* Check for parameter adjustment
*----------------------------------------------------------------------*
      st2 = 'ARGUMENT'
      IF (stcmp(st2,st1).EQ.1) THEN

        l_parameter = .true.

        narg = narg + 1
        CALL GETARG(narg,buff1)
        DO i=1,80
          param_file(i:i) = buff1(i:i)
        ENDDO

        narg = narg + 1
        CALL GETARG(narg,buff1)
        DO i=1,80
          st2(i:i) = buff1(i:i)
        ENDDO
        CALL STRIPBN(st2,n_param_0)

        narg = narg + 1
        CALL GETARG(narg,buff1)
        DO i=1,80
          st2(i:i) = buff1(i:i)
        ENDDO
        CALL STRIPBN(st2,n_param_f)

        sites = n_param_f - n_param_0 + 1
        n_param = n_param_0 - 1
        READ(98,*)
        READ(98,'(1000a)') st1

      ELSE
        l_parameter = .false.
      ENDIF

*----------------------------------------------------------------------*
* Read in type of landuse: 0 = defined by map; 1 = defined explicitly  *
* in the input file; 2 = natural vegetation based on average monthly   *
* temperatures.                                                        *
*----------------------------------------------------------------------*
      READ(st1,*) ilanduse
      IF (ilanduse.EQ.1) THEN
* Use landuse defined in input file.
        CALL LANDUSE1(luse,yr0,yrf)
      ENDIF
      IF ((ilanduse.LT.0).OR.(ilanduse.GT.2)) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'No landuse defined'
        WRITE(*,*) '0:= Defined from a map.'
        WRITE(*,*) '1:= Defined in the input file.'
        WRITE(*,*) '2:= Natural vegetation.'
        STOP
      ENDIF
      CLOSE(98)

*----------------------------------------------------------------------*
* Set up parameters for hard wired funcitonal types 'BARE' and 'CITY'. *
*----------------------------------------------------------------------*
      ftmor(1) = 1
      ftlls(1) = 0
      ftsls(1) = 0
      ftrls(1) = 0
      ftxyl(1) = 0.0d0
      ftpd(1) = 0.0d0
      ftstmx(1) = 0.0d0
      nppstore(1) = 0.0d0
      nppstorx(1) = 0.0d0
      nppstor2(1) = 0.0d0
      tsumam(1) = 0.0d0

      ftmor(2) = 1
      ftlls(2) = 0
      ftsls(2) = 0
      ftrls(2) = 0
      ftxyl(2) = 0.0d0
      ftpd(2) = 0.0d0
      ftstmx(2) = 0.0d0
      nppstore(2) = 0.0d0
      nppstorx(2) = 0.0d0
      nppstor2(2) = 0.0d0
      tsumam(2) = 0.0d0

*----------------------------------------------------------------------*
      srespm = 0.0d0

      nppstore(1) = 0.0d0
      nppstorx(1) = 0.0d0
      nppstor2(1) = 0.0d0

      co2(1) = 0.0d0

      DO ft=1,nft
        ftxyl(ft) =  ftxyl(ft)*1.0e-9
        ftpd(ft)  =  ftpd(ft)*1.0e+3
      ENDDO

      st1 = stoutput
      st2 = stinput
      st3 = stinit

*----------------------------------------------------------------------*
* Open snapshot files.                                                 *
*----------------------------------------------------------------------*
      IF (snp_no.GT.0) THEN
        DO i=1,snp_no
          fno = 100 + (i-1)*4
          st4=in2st(snpshts(i))
          CALL STRIPB(st4)
          OPEN(fno+1,FILE=st1(1:blank(st1))//
     &'/initbio_'//st4(1:blank(st4))//'.dat')
          OPEN(fno+2,FILE=st1(1:blank(st1))//
     &'/initcov_'//st4(1:blank(st4))//'.dat')
          OPEN(fno+3,FILE=st1(1:blank(st1))//
     &'/initppm_'//st4(1:blank(st4))//'.dat')
          OPEN(fno+4,FILE=st1(1:blank(st1))//
     &'/inithgt_'//st4(1:blank(st4))//'.dat')
        ENDDO
      ENDIF

*----------------------------------------------------------------------*
* Open default output files.                                           *
*----------------------------------------------------------------------*
      OPEN(21,FILE=st1(1:blank(st1))//'/lai.dat')
      OPEN(22,FILE=st1(1:blank(st1))//'/npp.dat')
      OPEN(23,FILE=st1(1:blank(st1))//'/scn.dat')
      OPEN(24,FILE=st1(1:blank(st1))//'/snn.dat')
      OPEN(25,FILE=st1(1:blank(st1))//'/nep.dat')
      OPEN(26,FILE=st1(1:blank(st1))//'/swc.dat')
      OPEN(27,FILE=st1(1:blank(st1))//'/biot.dat')
      OPEN(28,FILE=st1(1:blank(st1))//'/bioind.dat')
      OPEN(29,FILE=st1(1:blank(st1))//'/covind.dat')
      OPEN(30,FILE=st1(1:blank(st1))//'/dof.dat')
      OPEN(31,FILE=st1(1:blank(st1))//'/rof.dat')
      OPEN(32,FILE=st1(1:blank(st1))//'/fcn.dat')
      OPEN(33,FILE=st1(1:blank(st1))//'/nppstore.dat')
      OPEN(34,FILE=st1(1:blank(st1))//'/stembio.dat')
      OPEN(35,FILE=st1(1:blank(st1))//'/rootbio.dat')
      OPEN(36,FILE=st1(1:blank(st1))//'/leafper.dat')
      OPEN(37,FILE=st1(1:blank(st1))//'/stemper.dat')
      OPEN(38,FILE=st1(1:blank(st1))//'/rootper.dat')
      OPEN(39,FILE=st1(1:blank(st1))//'/sresp.dat')
      OPEN(40,FILE=st1(1:blank(st1))//'/evt.dat')
      OPEN(41,FILE=st1(1:blank(st1))//'/gpp.dat')
      OPEN(42,FILE=st1(1:blank(st1))//'/lch.dat')
      OPEN(43,FILE=st1(1:blank(st1))//'/prc.dat')
      OPEN(44,FILE=st1(1:blank(st1))//'/nbp.dat')
      OPEN(45,FILE=st1(1:blank(st1))//'/trn.dat')
      OPEN(46,FILE=st1(1:blank(st1))//'/fab.dat')
      OPEN(47,FILE=st1(1:blank(st1))//'/tmp.dat')
      OPEN(48,FILE=st1(1:blank(st1))//'/hum.dat')
      OPEN(49,FILE=st1(1:blank(st1))//'/ppm.dat')
      OPEN(50,FILE=st1(1:blank(st1))//'/pet.dat')

*----------------------------------------------------------------------*
* Open optional yearly cover, biomass, budburst and senescence files.  *
*----------------------------------------------------------------------*
      iofn = 200
      iofngft = iofn
      IF (outyears2.GT.0) THEN
      DO ft=1,nft
        IF (out_cov) THEN
          iofn = iofn + 1
          OPEN(iofn,FILE=st1(1:blank(st1))//'/cov_'//
     &fttags(ft)(1:blank(fttags(ft)))//'.dat')
        ENDIF
        IF (out_bio) THEN
          iofn = iofn + 1
          OPEN(iofn,FILE=st1(1:blank(st1))//'/bio_'//
     &fttags(ft)(1:blank(fttags(ft)))//'.dat')
        ENDIF
        IF (out_bud) THEN
          iofn = iofn + 1
          OPEN(iofn,FILE=st1(1:blank(st1))//'/bud_'//
     &fttags(ft)(1:blank(fttags(ft)))//'.dat')
        ENDIF
        IF (out_sen) THEN
          iofn = iofn + 1
          OPEN(iofn,FILE=st1(1:blank(st1))//'/sen_'//
     &fttags(ft)(1:blank(fttags(ft)))//'.dat')
        ENDIF
      ENDDO
      ENDIF

*----------------------------------------------------------------------*
* Open input files for state vector.                                   *
*----------------------------------------------------------------------*
      IF (.NOT.initise) THEN
        OPEN(70,FILE=st3(1:blank(st3))//'/init.dat',status='old',
     &iostat=kode)
        IF (kode.NE.0) THEN
          WRITE(*,'('' PROGRAM TERMINATED'')')
          WRITE(*,*) 'Initialisation directory/file does not exist.'
          WRITE(*,'('' "'',A,''"'')') st3(1:blank(st3))
          STOP
        ENDIF
        OPEN(71,FILE=st3(1:blank(st3))//'/initbio.dat')
        OPEN(72,FILE=st3(1:blank(st3))//'/initcov.dat')
        OPEN(73,FILE=st3(1:blank(st3))//'/initppm.dat')
        OPEN(74,FILE=st3(1:blank(st3))//'/inithgt.dat')
        OPEN(75,FILE=st3(1:blank(st3))//'/initwdt.dat')
        OPEN(76,FILE=st3(1:blank(st3))//'/initleaf.dat')
        OPEN(77,FILE=st3(1:blank(st3))//'/initstem.dat')
        OPEN(78,FILE=st3(1:blank(st3))//'/initroot.dat')
        OPEN(79,FILE=st3(1:blank(st3))//'/initmisc.dat')
      ENDIF

*----------------------------------------------------------------------*
* Open output files for state vector.                                  *
*----------------------------------------------------------------------*
      OPEN(80,FILE=st1(1:blank(st1))//'/init.dat')
      IF (initiseo) THEN
        OPEN(81,FILE=st1(1:blank(st1))//'/initbio.dat')
        OPEN(82,FILE=st1(1:blank(st1))//'/initcov.dat')
        OPEN(83,FILE=st1(1:blank(st1))//'/initppm.dat')
        OPEN(84,FILE=st1(1:blank(st1))//'/inithgt.dat')
        OPEN(85,FILE=st1(1:blank(st1))//'/initwdt.dat')
        OPEN(86,FILE=st1(1:blank(st1))//'/initleaf.dat')
        OPEN(87,FILE=st1(1:blank(st1))//'/initstem.dat')
        OPEN(88,FILE=st1(1:blank(st1))//'/initroot.dat')
        OPEN(89,FILE=st1(1:blank(st1))//'/initmisc.dat')
      ENDIF

*----------------------------------------------------------------------*
*  Read co2 file.                                                      *
*----------------------------------------------------------------------*
      CALL READCO2 (stco2,yr0,yrf,co2)
*----------------------------------------------------------------------*

      site_dat = 0

*----------------------------------------------------------------------*
*                               Site loop                              *
*----------------------------------------------------------------------*
      DO site=1,sites

        speedc = xspeedc
        swcnew = 0.0d0

        IF (l_parameter) THEN
          lat = lat_lon(1,1)
          lon = lat_lon(1,2)
        ELSE
          lat = lat_lon(site,1)
          lon = lat_lon(site,2)
        ENDIF

        IF (abs(xseed1).EQ.0) THEN
          IF (site.EQ.1) THEN
            seed1 = int(SETARANDOM()*10000.0d0+.5)
            seed2 = 2*seed1
            seed3 = 3*seed1
          ENDIF
        ELSEIF (xseed1.LT.0.0) THEN
          IF (site.EQ.1) THEN
            seed1 = -xseed1
            seed2 = 2*seed1
            seed3 = 3*seed1
          ENDIF
        ELSE
          seed1 = xseed1
          seed2 = 2*seed1
          seed3 = 3*seed1
        ENDIF

      l_countries = .TRUE.
      IF (no_countries.GT.0) THEN
        l_countries = .FALSE.
        CALL COUNTRY(stmask,lat,lon,country_name,country_id,.true.)
        DO i=1,no_countries
          IF (blank(country_name).EQ.blank(countries(i))) THEN
            check_c = .TRUE.
            DO j=1,blank(country_name)
              IF (country_name(j:j).NE.countries(i)(j:j)) 
     &check_c = .FALSE.
            ENDDO
          ELSE
            check_c = .FALSE.
          ENDIF
          IF (check_c) l_countries = .TRUE.
        ENDDO
        CALL COUNTRY(stmask,lat,lon,country_name,country_id,.false.)
        DO i=1,no_countries
          IF (blank(country_name).EQ.blank(countries(i))) THEN
            check_c = .TRUE.
            DO j=1,blank(country_name)
              IF (country_name(j:j).NE.countries(i)(j:j)) 
     &check_c = .FALSE.
            ENDDO
          ELSE
            check_c = .FALSE.
          ENDIF
          IF (check_c) l_countries = .TRUE.
        ENDDO
      ENDIF

      IF (l_countries) THEN
*----------------------------------------------------------------------*
* Read in climate.                                                     *
*----------------------------------------------------------------------*
      IF (clim_type.EQ.1) THEN
*----------------------------------------------------------------------*
* DAILY Gridded.                                                       *
*----------------------------------------------------------------------*
        CALL EX_CLIM(st2,lat,lon,xlatf,xlatres,xlatresn,xlon0,xlonres,
     &xlonresn,yr0,yrf,xtmpv,xhumv,xprcv,isite,xyear0,xyearf,
     &siteno,du)
        withcloudcover=.FALSE.
      ELSEIF (clim_type.EQ.2) THEN
*----------------------------------------------------------------------*
* MONTHLY Gridded.                                                     *
*----------------------------------------------------------------------*
        CALL EX_CLIM_WEATHER_GENERATOR(st2,ststats,lat,lon,xlatf,
     &xlatres,xlatresn,xlon0,xlonres,xlonresn,yr0,yrf,xtmpv,xhumv,xprcv,
     &xcldv,isite,xyear0,xyearf,du,seed1,seed2,seed3,l_clim,l_stats)
        withcloudcover=.TRUE.
      ELSEIF (clim_type.EQ.3) THEN
*----------------------------------------------------------------------*
* DAILY Site                                                           *
*----------------------------------------------------------------------*
        CALL EX_CLIM_SITE(st2,yr0,yrf,xtmpv,xhumv,xprcv,xyear0,xyearf)
        withcloudcover=.FALSE.
        siteno = 1
        l_clim = .TRUE.
        l_stats = .TRUE.
      ELSEIF (clim_type.EQ.4) THEN
*----------------------------------------------------------------------*
* DAILY Site                                                           *
*----------------------------------------------------------------------*
        CALL EX_CLIM_SITE_MONTH(st2,yr0,yrf,xtmpv,xhumv,xprcv,xcldv,
     &xyear0,xyearf)
        withcloudcover=.FALSE.
        siteno = 1
        CALL GENERATE_MONTHLY(ststats,yr0,yrf,xtmpv,xhumv,xprcv,
     &seed1,seed2,seed3)
        l_clim = .TRUE.
        l_stats = .TRUE.
      ELSE
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'Error defining climate to read'
        STOP
      ENDIF

*----------------------------------------------------------------------*
* Read in soil parameters.                                             *
*----------------------------------------------------------------------*
      IF (l_clim) THEN
        CALL SOIL_PARAMETERS(ts,tsi,tc,bulk,orgc,wilt,field,sat,dep,
     &sttxdp,lat,lon,soil_chr,soil_chr2,du,l_soil)
      ENDIF

*----------------------------------------------------------------------*
* IF in parameter adjustment mode read in parameters and perform       *
* necessary adjustments.                                               *
*----------------------------------------------------------------------*
      IF (l_parameter) THEN
        n_param = n_param + 1
        CALL PARAMETER_ADJUSTMENT(param_file,n_param)
        ts = p_sand
        tsi = p_silt
        tc = 100.0 - ts - tsi
        bulk = p_bulk
        orgc = p_orgc
        dep = p_dep
        topsl = p_topsl
      ENDIF

*----------------------------------------------------------------------*
* Read in landuse/cover.                                               *
*----------------------------------------------------------------------*
c     temporaire !
      icontinuouslanduse=1

      IF (ilanduse.EQ.0) THEN
        IF (icontinuouslanduse.EQ.0) THEN
          CALL EX_LU(stlu,lat,lon,luse,yr0,yrf,du)
c     create the continuous land use (cluse)
          DO year=yr0,yrf
            DO ft=1,nft
              cluse(ft,year-yr0+1) = lutab(luse(year-yr0+1),ft)
            ENDDO
          ENDDO
        ELSE            
          CALL EX_CLU(stlu,lat,lon,nft,lutab,cluse,yr0,yrf,du,l_lu)
        ENDIF
      ELSEIF (ilanduse.EQ.1) THEN
        l_lu = .TRUE.
        DO year=yr0,yrf
          DO ft=1,nft
            cluse(ft,year-yr0+1) = lutab(luse(year-yr0+1),ft)
          ENDDO
        ENDDO
      ELSEIF (ilanduse.EQ.2) THEN
        WRITE(*,*) 'Checking natural vegetation types exist:'
        WRITE(*,*) 'BARE CITY C3 C4 Ev_Bl Ev_Nl Dc_Bl Dc_Nl.'
        l_lu = .TRUE.
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
        DO year=yr0,yrf
          DO ft=1,nft
            cluse(ft,year-yr0+1) = 0.0d0
          ENDDO
        ENDDO
        cluse(1,1) = 100.0d0
      ENDIF

*----------------------------------------------------------------------*
      ENDIF
*----------------------------------------------------------------------*
* End of countries check.                                              *
*----------------------------------------------------------------------*

*----------------------------------------------------------------------*
* Start of climate data exists 'if' statement.                         *
*----------------------------------------------------------------------*
      IF ((l_clim).AND.(l_stats).AND.(l_soil(1)).AND.(l_soil(3)).AND.
     &(l_soil(8)).AND.(l_lu).AND.(l_countries)) THEN

        site_dat = site_dat + 1
        IF (mod(site_dat,max(site_out,1)).EQ.min(1,site_out)-1) 
     &WRITE(*,'( '' Site no. '',i5,'', Lat ='',f7.3,'', Lon ='',f9.3,
     &'' CO2 = '',f9.3)') 
     &site_dat,lat,lon,ca

*----------------------------------------------------------------------*
* Write lat lon for general output files.                              *
*----------------------------------------------------------------------*
        DO i=21,50
          WRITE(i,'(f7.3,f9.3,$)') lat,lon
        ENDDO

*----------------------------------------------------------------------*
* Write lat lon for optional yearly ft files                           *
*----------------------------------------------------------------------*
        iofn = iofngft
        IF (outyears2.GT.0) THEN
        DO ft=1,nft
          IF (out_cov) THEN
            iofn = iofn + 1
            WRITE(iofn,'(f7.3,f9.3,$)') lat,lon
          ENDIF
          IF (out_bio) THEN
            iofn = iofn + 1
            WRITE(iofn,'(f7.3,f9.3,$)') lat,lon
          ENDIF
          IF (out_bud) THEN
            iofn = iofn + 1
            WRITE(iofn,'(f7.3,f9.3,$)') lat,lon
          ENDIF
          IF (out_sen) THEN
            iofn = iofn + 1
            WRITE(iofn,'(f7.3,f9.3,$)') lat,lon
          ENDIF
        ENDDO
        ENDIF

*----------------------------------------------------------------------*
* Write lat,lon for monthly/daily output files.                        *
*----------------------------------------------------------------------*
        iofn = 400
        IF (outyears1.GT.0) THEN
        IF ((oymd.EQ.1).OR.(oymd.EQ.3)) THEN
          DO i=1,50
            IF (otagsn(i).EQ.1) THEN
              iofn = iofn + 1
              WRITE(iofn,'(f7.3,f9.3)') lat,lon
            ENDIF
          ENDDO
        ENDIF
        IF ((oymd.EQ.2).OR.(oymd.EQ.3)) THEN
          DO i=1,50
            IF (otagsn(i).EQ.1) THEN
              iofn = iofn + 1
              WRITE(iofn,'(f7.3,f9.3)') lat,lon
            ENDIF
          ENDDO
        ENDIF
        ENDIF

*----------------------------------------------------------------------*
* Write lat,lon for monthly/daily ft output files.                     *
*----------------------------------------------------------------------*
        IF (outyears2.GT.0) THEN
        IF ((oymdft.EQ.1).OR.(oymdft.EQ.3)) THEN
          DO i=1,50
            IF (otagsnft(i).EQ.1) THEN
              DO ft=1,nft
                iofn = iofn + 1
                WRITE(iofn,'(f7.3,f9.3)') lat,lon
             ENDDO
            ENDIF
          ENDDO
        ENDIF
        IF ((oymdft.EQ.2).OR.(oymdft.EQ.3)) THEN
          DO i=1,50
            IF (otagsnft(i).EQ.1) THEN
              DO ft=1,nft
                iofn = iofn + 1
                WRITE(iofn,'(f7.3,f9.3)') lat,lon
              ENDDO
            ENDIF
          ENDDO
        ENDIF
        ENDIF

c other initialisation of water parameter
c npuc is unused
      CALL WSPARAM(l_b_and_c,wilt,field,sat,nupc,awl,kd,kx,nci,infix,
     &adp,sfc,sw,sswc,ts,tc,tsi,bulk,orgc,dep,topsl,wtfc,wtwp,
     &wtswc,l_parameter)

*----------------------------------------------------------------------*
* Write site info to 'site_indo.dat'.                                  *
*----------------------------------------------------------------------*
      CALL CO2_0_F(co20,co2f,yearv,yr0,spinl,co2,co2const,nyears)

      CALL COUNTRY(stmask,lat,lon,country_name,country_id,l_regional)

      WRITE(12,'(a15,i6,f9.3,f9.3,1x,2f6.1,1x,2f7.1,f7.3,f7.2,1x,
     &3f7.3,f7.1,1x,200(2f6.1,1x))')
     & country_name,country_id,lat,lon,co20,co2f,ts,tsi,bulk,orgc,
     &wtwp,wtfc,wtswc,dep,
     &(cluse(ft,yearv(1)-yr0+1),cluse(ft,yearv(nyears)-yr0+1),ft=1,nft)

*----------------------------------------------------------------------*
* Initialise biomass and cover arrays.                                 *
*----------------------------------------------------------------------*
        IF (initise) THEN
          DO ft=1,nft
            IF (cluse(ft,1).GT.0.0d0) THEN
              DO age=1,ftmor(ft)
                cov(age,ft) = cluse(ft,1)/100.0d0/real(ftmor(ft))
                bio(age,1,ft) = 1000.0d0
                bio(age,2,ft) = 100.0d0
              ENDDO
              bio(1,1,1) = 0.0d0
              bio(1,2,1) = 0.0d0
              DO age=1,ftmor(ft)
                ppm(age,ft) = 0.001d0
                hgt(age,ft) = 5.0d0
                wdt(age,ft) = 0.5d0
              ENDDO
              nppstore(ft) = ftstmx(ft)
              nppstorx(ft) = ftstmx(ft)
              nppstor2(ft) = ftstmx(ft)
            ELSE
              DO age=1,ftmor(ft)
                cov(age,ft) = 0.0d0
                bio(age,1,ft) = 0.0d0
                bio(age,2,ft) = 0.0d0
              ENDDO
              bio(1,1,1) = 0.0d0
              bio(1,2,1) = 0.0d0
              DO age=1,ftmor(ft)
                ppm(age,ft) = 0.0d0
                hgt(age,ft) = 0.0d0
                wdt(age,ft) = 0.0d0
              ENDDO
              nppstore(ft) = 0.0d0
              nppstorx(ft) = 0.0d0
              nppstor2(ft) = 0.0d0
            ENDIF
            nppstore(2) = ftstmx(2)
            nppstore(3) = ftstmx(3)
            nppstorx(2) = ftstmx(2)
            nppstorx(3) = ftstmx(3)
            nppstor2(2) = ftstmx(2)
            nppstor2(3) = ftstmx(3)

            slc(ft) = 0.0d0
            rlc(ft) = 0.0d0
            sln(ft) = 0.0d0
            rln(ft) = 0.0d0

            DO day=1,ftlls(ft)
              leafdp(day,ft) = 0.0d0
            ENDDO
            DO day=1,ftsls(ft)
              stemdp(day,ft) = 0.0d0
            ENDDO
            DO day=1,ftrls(ft)
              rootdp(day,ft) = 0.0d0
            ENDDO
            DO day=1,360
              sumadp(day,ft) = 0.0d0
            ENDDO
            lai(ft) = 0.0d0
            npp(ft) = 400.0d0
            nps(ft) = 40.0d0
            npl(ft) = 40.0d0
            evp(ft) = 1.0d0
            dof(ft) = 0.0
            maxlai(ft) = 2.0d0

            bb(ft) = 0
            ss(ft) = 0
            bbgs(ft) = 0
            dsbb(ft) = 0

            DO day=1,30
              sm_trig(day,ft) = 0.0d0
            ENDDO

            tsumam(ft) = 0.0d0
            stemfr(ft) =10.0d0
          ENDDO ! ft

          is1 = 5.0d0
          is2 = 50.0d0
          is3 = 50.0d0
          is4 = 50.0d0
          isn = 0.0d0
          ilsn = 0.0d0
          DO i=1,8
            ic0(i) = 1000.0d0
            in0(i) = 100.0d0
          ENDDO
          ic0(7) = 5000.0d0
          ic0(8) =  5000.0d0
          iminn(1) = 0.0d0
          iminn(2) = 0.0d0
          iminn(3) = 0.0d0
          dslc = 0.0d0
          drlc = 0.0d0
          dsln = 0.0d0
          drln = 0.0d0
          DO day=1,200
            tmem(day) = real(xtmpv(1,(day-1)/30+1,mod(day-1,30)+1))
     &/100.0d0
          ENDDO
          DO ft=1,nft
            s1(ft) = is1
            s2(ft) = is2
            s3(ft) = is3
            s4(ft) = is4
            sn(ft) = isn
            lsn(ft) = ilsn
          ENDDO

        ELSE

          IF ((abs(lat-zlat).GT.0.001).OR.
     &(abs(lon-zlon).GT.0.001)) THEN

          READ(70,*) zlat,zlon
          IF ((abs(lat-zlat).GT.0.001).OR.
     &(abs(lon-zlon).GT.0.001)) THEN
            WRITE(*,'("Error lat and lon dont match.")')
            WRITE(*,*) lat,lon
            WRITE(*,*) zlat,zlon
            STOP
          ENDIF
          READ(70,*) (zs1(ft),ft=1,nft)
          READ(70,*) (zs2(ft),ft=1,nft)
          READ(70,*) (zs3(ft),ft=1,nft)
          READ(70,*) (zs4(ft),ft=1,nft)
          READ(70,*) (zsn(ft),ft=1,nft)
          READ(70,*) (zlsn(ft),ft=1,nft)
          READ(70,*) (zic0(i),i=1,8)
          READ(70,*) (zin0(i),i=1,8)
          READ(70,*) (ziminn(i),i=1,3)
          READ(70,*) zdslc
          READ(70,*) zdrlc
          READ(70,*) zdsln
          READ(70,*) zdrln
          READ(70,*) zlai(4)
          READ(70,*) znpp(4)

          READ(72,*) zlat,zlon,((zcov(i,j),i=1,ftmor(j)),j=2,nft)
          READ(71,*) 
     &zlat,zlon,(((zbio(i,j,k),j=1,2),i=1,ftmor(k)),k=2,nft)
          READ(73,*) zlat,zlon,((zppm(i,j),i=1,ftmor(j)),j=2,nft)
          READ(74,*) zlat,zlon,((zhgt(i,j),i=1,ftmor(j)),j=2,nft)
          READ(75,*) zlat,zlon,((zwdt(i,j),i=1,ftmor(j)),j=2,nft)
          READ(76,*) zlat,zlon,((zleafdp(i,j),i=1,ftlls(j)),j=2,nft)
          READ(77,*) zlat,zlon,((zstemdp(i,j),i=1,ftsls(j)),j=2,nft)
          READ(78,*) zlat,zlon,((zrootdp(i,j),i=1,ftrls(j)),j=2,nft)

          READ(79,*) zlat,zlon
          READ(79,*) (zlai(i),i=1,nft)
          READ(79,*) (znpp(i),i=1,nft)
          READ(79,*) (znps(i),i=1,nft)
          READ(79,*) (znpl(i),i=1,nft)
          READ(79,*) (zevp(i),i=1,nft)
          READ(79,*) (zdof(i),i=1,nft)
          READ(79,*) (zslc(i),i=1,nft)
          READ(79,*) (zrlc(i),i=1,nft)
          READ(79,*) (zsln(i),i=1,nft)
          READ(79,*) (zrln(i),i=1,nft)
          READ(79,*) (znppstore(i),i=1,nft)
          READ(79,*) (znppstorx(i),i=1,nft)
          READ(79,*) (znppstor2(i),i=1,nft)
          READ(79,*) (zbb(i),i=1,nft)
          READ(79,*) (zbbgs(i),i=1,nft)
          READ(79,*) (zdsbb(i),i=1,nft)
          READ(79,*) (zmaxlai(i),i=1,nft)
          READ(79,*) (ztmem(i),i=1,200)
          DO j=1,nft
            READ(79,*) (zsm_trig(i,j),i=1,30)
          ENDDO
          DO ft=3,nft
            READ(79,*) (zsumadp(i,ft),i=1,360)
          ENDDO
          READ(79,*) (zstemfr(ft),ft=1,nft)

          ENDIF

          DO ft=1,nft
            s1(ft) = zs1(ft)
            s2(ft) = zs2(ft)
            s3(ft) = zs3(ft)
            s4(ft) = zs4(ft)
            sn(ft) = zsn(ft)
            lsn(ft) = zlsn(ft)
            chill(ft) = 0
            dschill(ft) = 0
          ENDDO
          DO i=1,8
            ic0(i) = zic0(i)
            in0(i) = zin0(i)
          ENDDO
          DO i=1,3
            iminn(i) = ziminn(i)
          ENDDO
          dslc = zdslc
          drlc = zdrlc
          dsln = zdsln
          drln = zdrln
          lai(4) = zlai(4)
          npp(4) = znpp(4)

          DO j=2,nft
            DO i=1,ftmor(j)
              cov(i,j) = zcov(i,j)
              bio(i,1,j) = zbio(i,1,j)
              bio(i,2,j) = zbio(i,2,j)
              ppm(i,j) = zppm(i,j)
              hgt(i,j) = zhgt(i,j)
              wdt(i,j) = zwdt(i,j)
            ENDDO
          ENDDO
          DO j=2,nft
            DO i=1,ftlls(j)
              leafdp(i,j) = zleafdp(i,j)
            ENDDO
          ENDDO
          DO j=2,nft
            DO i=1,ftsls(j)
              stemdp(i,j) = zstemdp(i,j)
            ENDDO
          ENDDO
          DO j=2,nft
            DO i=1,ftrls(j)
              rootdp(i,j) = zrootdp(i,j)
            ENDDO
          ENDDO

          DO i=1,nft
            lai(i) = zlai(i)
            npp(i) = znpp(i)
            nps(i) = znps(i)
            npl(i) = znpl(i)
            evp(i) = zevp(i)
            dof(i) = zdof(i)
            slc(i) = zslc(i)
            rlc(i) = zrlc(i)
            sln(i) = zsln(i)
            rln(i) = zrln(i)
            nppstore(i) = znppstore(i)
            nppstorx(i) = znppstorx(i)
            nppstor2(i) = znppstor2(i)
            bb(i) = zbb(i)
            bbgs(i) = zbbgs(i)
            dsbb(i) = zdsbb(i)
            maxlai(i) = zmaxlai(i)
          ENDDO

          DO i=1,200
            tmem = ztmem(i)
          ENDDO
          DO j=1,nft
            DO i=1,30
              sm_trig(i,j) = zsm_trig(i,j)
            ENDDO
          ENDDO
          DO ft=3,nft
            DO i=1,360
              sumadp(i,ft) = zsumadp(i,ft)
            ENDDO
          ENDDO
          DO ft=1,nft
            stemfr(ft) = zstemfr(ft)
          ENDDO

        ENDIF

*----------------------------------------------------------------------*
* Ensure cover array sums to 1.                                        *
*----------------------------------------------------------------------*
        bio(1,1,1)=0.0
        bio(1,2,1)=0.0
        hgt(1,1)=0.0
        wdt(1,1)=0.0
        DO ft=2,nft
          leaflit(ft) = 0.0d0
          yield(ft) = 0.0d0
          stemlit(ft) = 0.0d0
          rootlit(ft) = 0.0d0
          ftcov(ft) = 0.0d0
          DO age=1,ftmor(ft)
*            IF (ppm(age,ft).LT.1.0e-16) then
*              cov(age,ft)=0.0
*              bio(age,1,ft)=0.0
*              bio(age,2,ft)=0.0
*              hgt(age,ft)=0.0
*              wdt(age,ft)=0.0
*            ENDIF
            ftcov(ft) = ftcov(ft) + cov(age,ft)
          ENDDO
        ENDDO

        sum = 0.0d0
        DO ft=2,nft
          sum = sum + ftcov(ft)
        ENDDO

        IF (sum.LT.1.0d0) THEN
          cov(1,1) = 1.0d0 - sum
          ftcov(1) = cov(1,1)
        ELSE
          cov(1,1) = 0.0d0
          ftcov(1) = 0.0d0
          DO ft=2,nft
            ftcov(ft) = ftcov(ft)/sum
            DO age=1,ftmor(ft)
              cov(age,ft)=cov(age,ft)/sum
            ENDDO
          ENDDO
        ENDIF
*----------------------------------------------------------------------*

        biotoo = 0.0d0
        DO ft=1,nft
          DO age=1,ftmor(ft)
            biotoo = biotoo + (bio(age,1,ft) + bio(age,2,ft) 
     &)*cov(age,ft) 
          ENDDO
        ENDDO
        solcoo = 0.0d0
        DO i=1,8
          solcoo = solcoo + ic0(i)
        ENDDO

        c3old = cov(1,2)
        c4old = cov(1,3)

        snp_year = 1

        soilt = 10.0d0

*----------------------------------------------------------------------*
*                               Year Loop                              *
*----------------------------------------------------------------------*
      DO iyear=1,nyears

*        if (iyear.eq.3) stop

        year = yearv(iyear)

        DO ft=1,nft
          laimax(ft) = 4.601d0
        ENDDO

        nfix = infix

*----------------------------------------------------------------------*
* Set CO2.                                                             *
*----------------------------------------------------------------------*
        IF ((spinl.gt.0).AND.(iyear.GT.spinl)) THEN
          speedc = .FALSE.
          ca = co2(year-yr0+1)
        ELSE
          IF (co2const.GT.0.0d0) THEN
            ca = co2const
          ELSE
*            ca = co2(year-yr0+1)
            ca = co2(iyear)
          ENDIF
        ENDIF

        IF (mod(iyear,max(year_out,1)).EQ.min(1,year_out)-1) THEN
          WRITE(*,'('' Year no.'',2i5,'', ca = '',f6.2)') iyear,year,ca
        ENDIF

*----------------------------------------------------------------------*
* Set 'tmp' 'hum' 'prc', calculate: minimum temp, total precip.        *
*----------------------------------------------------------------------*
        DO mnth=1,12
          DO day=1,no_days(year,mnth,thty_dys)
            tmp(mnth,day) = 
     &real(xtmpv(yearv(iyear)-yr0+1,mnth,day))/100.0d0
            prc(mnth,day) = 
     &real(xprcv(yearv(iyear)-yr0+1,mnth,day))/10.0d0
            hum(mnth,day) = 
     &real(xhumv(yearv(iyear)-yr0+1,mnth,day))/100.0d0
            IF (withcloudcover) THEN
               cld(mnth) = 
     &real(xcldv(yearv(iyear)-yr0+1,mnth))/1000.0d0
            ELSE
               cld(mnth) = 0.5d0
            ENDIF
          ENDDO
        ENDDO

        tmin = 100.0d0
        yeartmp = 0.0d0
        yearprc = 0.0d0
        yearhum = 0.0d0
        DO mnth=1,12
          mnthprc(mnth) = 0.0d0
          mnthtmp(mnth) = 0.0d0
          mnthhum(mnth) = 0.0d0
          DO day=1,no_days(year,mnth,thty_dys)
            IF (hum(mnth,day).LT.30.0d0)  hum(mnth,day) = 30.0d0
            IF (hum(mnth,day).GT.95.0d0)  hum(mnth,day) = 95.0d0
            mnthprc(mnth) = mnthprc(mnth) + 
     &prc(mnth,day)
            mnthtmp(mnth) = mnthtmp(mnth) + 
     &tmp(mnth,day)/no_days(year,mnth,thty_dys)
            mnthhum(mnth) = mnthhum(mnth) + 
     &hum(mnth,day)/no_days(year,mnth,thty_dys)
          ENDDO ! day
          IF (mnthtmp(mnth).LT.tmin)  tmin = mnthtmp(mnth)
          yeartmp = yeartmp + mnthtmp(mnth)/12.0d0
          yearprc = yearprc + mnthprc(mnth)
          yearhum = yearhum + mnthhum(mnth)/12.0d0
        ENDDO ! month
        tmin = tmin*1.29772d0 - 19.5362d0

*----------------------------------------------------------------------*
* Carbon at the start of the year.                                     *
*----------------------------------------------------------------------*
        DO ft=1,nft
          bioleaf(ft) = 0.0d0
          DO day=1,ftlls(ft)
            bioleaf(ft) = bioleaf(ft) + 
     &leafdp(day,ft)*12.0d0/ftsla(ft)/18.0d0
          ENDDO
        ENDDO

        ans1 = 0.0d0
        DO ft=1,nft
          DO i=1,ftmor(ft)
            ans1 = ans1 + (bio(i,1,ft) + bio(i,2,ft) + bioleaf(ft) +
     &nppstore(ft))*cov(i,ft)
          ENDDO
          ans1 = ans1 + slc(ft) + rlc(ft)
        ENDDO
        ccheck = ans1 + ic0(1) + 
     &ic0(2) + ic0(3) + ic0(4) + ic0(5) + ic0(6) + ic0(7) + ic0(8)

*----------------------------------------------------------------------*
* Set land use through ftprop.                                         *
*----------------------------------------------------------------------*
        IF (ilanduse.EQ.2) THEN
          CALL NATURAL_VEG(tmp,prc,ftprop,nat_map)
        ELSE ! 0 and 1
          ftprop(1) = 100.0d0
          DO ft=2,nft
            IF (check_ft_grow(tmp,ftbbm(ft),ftbb0(ft),ftbbmax(ft),
     &ftbblim(ft),chill(ft),dschill(ft)).EQ.1) THEN
              ftprop(ft) = cluse(ft,year-yr0+1)
              ftprop(1) = ftprop(1) - ftprop(ft)
            ELSE
              ftprop(ft) = 0.0d0
            ENDIF
          ENDDO
          IF (ftprop(1).LT.0.0d0) THEN
            DO ft=2,nft
              ftprop(ft) = ftprop(ft)*100.0d0/(100.0 - ftprop(1))
            ENDDO
            ftprop(1) = 0.0d0
          ENDIF
        ENDIF

*----------------------------------------------------------------------*
        CALL COVER(nft,ftmor,ftppm0,cov,bio,bioleaf,nppstore,
     &npp,nps,mnthtmp,mnthprc,slc,rlc,c3old,c4old,firec,ppm,hgt,fireres,
     &fprob,ftprop,ftstmx,stemdp,rootdp,ftsls,ftrls,ilanduse,nat_map,
     &ic0)

        CALL MKDLIT(nft,ftmor,ftcov,dslc,drlc,dsln,drln,cov,slc,rlc,sln,
     &rln)

        DO i=1,8
          tc0(i) = 0.0d0
          tn0(i) = 0.0d0
        ENDDO

        DO i=1,3
          tminn(i) = 0.0d0
        ENDDO
        ts1 = 0.0d0
        ts2 = 0.0d0
        ts3 = 0.0d0
        ts4 = 0.0d0
        tsn = 0.0d0
        tlsn = 0.0d0

*----------------------------------------------------------------------*
* Find which DOLY runs are required.                                   *
*----------------------------------------------------------------------*
        DO ft=1,nft
          dolydo(ft) = 0
          DO age=1,ftmor(ft)
            IF (cov(age,ft).GT.0.0d0) THEN
              dolydo(ft) = 1
              GOTO 129
            ENDIF
          ENDDO
129       CONTINUE
        ENDDO

*----------------------------------------------------------------------*
* Initialisations that were in doly at the beginning of the year       *
*----------------------------------------------------------------------*

        iminn(3) = iminn(1) + iminn(2)

        isoilc = 0.0d0
        isoiln = 0.0d0
        DO i=1,8
           isoilc = isoilc + ic0(i)
           isoiln = isoiln + in0(i)
        ENDDO
        isoiln = isoiln + iminn(3)

c initialisse for all the ft, even for those doly run is not required
        DO ft=1,nft
          leaflit(ft) = 0.0d0
          yield(ft) = 0.0d0
          stemlit(ft) = 0.0d0
          rootlit(ft) = 0.0d0

          leafnpp(ft) = 0.0d0
          stemnpp(ft) = 0.0d0
          rootnpp(ft) = 0.0d0

          dof(ft) = 0.0d0

          lai(ft) = 0.0d0
          DO day=1,ftlls(ft)
            lai(ft) = lai(ft) + leafdp(day,ft)
          ENDDO

c     water
*           s1(ft) = is1
*           s2(ft) = is2
*           s3(ft) = is3
*           s4(ft) = is4
*           sn(ft) = isn
*           lsn(ft) = ilsn
          DO i=1,8
            c0(i,ft) = ic0(i)
            n0(i,ft) = in0(i)
          ENDDO
          DO i=1,3
            minn(i,ft) = iminn(i)
          ENDDO
          soilc(ft) = isoilc
          soiln(ft) = isoiln

          evp(ft) = 0.0d0

          sresp(ft) = 0.0d0
          lch(ft) = 0.0d0

* Height (m)
cccn           ht(ft) = 0.807d0*(laimax(ft)**2.13655d0)
          ht(ft) = 0.807d0*(5.0d0**2.13655d0)
          IF (ht(ft).GT.50.0d0)  ht(ft)=50.0d0
          laimax(ft)=0.0d0   
*This is to do what doly did... but 
*it's strange. Initialisation to 0.0d0 is ok, but the previous 
*one laimax=4.6d0 above seems not good.


*     Convert nppstore to mols
          nppstore(ft) = nppstore(ft)/12.0d0
          nppstorx(ft) = nppstorx(ft)/12.0d0
          nppstor2(ft) = nppstor2(ft)/12.0d0
          nppstoreold(ft) = nppstore(ft)

          budo(ft) = 0
          seno(ft) = 0

        ENDDO ! ft's

*      print*,'start monthly loop'

*----------------------------------------------------------------------*
* START OF MONTH LOOP                                                  *
*----------------------------------------------------------------------*
      DO mnth=1,12

c     monthly initialisations
        DO ft=1,nft
           flow1(ft) = 0.0d0
           flow2(ft) = 0.0d0

           evapm(mnth,ft) = 0.0d0
           tranm(mnth,ft) = 0.0d0
           roffm(mnth,ft) = 0.0d0
           petm(mnth,ft) = 0.0d0

           photm(mnth,ft) = 0.0d0

           laimnth(mnth,ft) = 0.0d0
           avmnpet(ft) = 0.0d0
        ENDDO

        avmnppt = 0.0d0
        avmnt = 0.0d0
*----------------------------------------------------------------------*
* DAILY LOOP.                                                          *
*----------------------------------------------------------------------*
        DO day=1,no_days(year,mnth,thty_dys)
          fpr=0.0d0

*          if ((year.eq.1906).and.(day.eq.21).and.(mnth.eq.11)) stop

          avmnppt = avmnppt + prc(mnth,day)
          avmnt = avmnt + tmp(mnth,day)/no_days(year,mnth,thty_dys)

*----------------------------------------------------------------------*
* Updata daily climate memory.                                         *
*----------------------------------------------------------------------*
          DO i=1,199
            tmem(201-i) = tmem(200-i)
          ENDDO
          tmem(1) = tmp(mnth,day)

*----------------------------------------------------------------------*
* Mix water resources.                                                 *
*----------------------------------------------------------------------*
          CALL MIX_WATER(s1,s2,s3,s4,sn,lsn,ftcov,ftmix,nft)
*----------------------------------------------------------------------*

          DO ft=1,nft
          IF (dolydo(ft).EQ.1) THEN

            DO d=1,ftlls(ft)
              leafv(d)=leafdp(d,ft)
            ENDDO
            DO d=1,ftsls(ft)
              stemv(d)=stemdp(d,ft)
            ENDDO
            DO d=1,ftrls(ft)
              rootv(d)=rootdp(d,ft)
            ENDDO

            DO d=1,30
              smtrig(d) = sm_trig(d,ft)
            ENDDO

            DO d=1,360
              suma(d) = sumadp(d,ft)
            ENDDO

*This is not really usefull since every ft's have the same minn 
*at the moment... but this may (should) be changed in the future
            DO i=1,3
              minnv(i)=minn(i,ft)
            ENDDO
            daygpp = 0.0d0
            dayra = 0.0d0
            evap = 0.0d0
            tran = 0.0d0
            roff = 0.0d0
            pet =0.0d0
            fpr = 0.0d0

            leafold = leafnpp(ft)
            stemold = stemnpp(ft)
            rootold = rootnpp(ft)
            nppsold = nppstore(ft)

*----------------------------------------------------------------------*
* nppstore mols
* leafnpp  mols
* stemnpp  mols
* rootnpp  mols
* leaflit  mols
* stemlit  mols
* rootlit  mols
*----------------------------------------------------------------------*
            oldlai = lai(ft)

            lflitold = leaflit(ft)
!      print*,'dol ',tmp(mnth,day),prc(mnth,day),hum(mnth,day),cld(mnth), 
!     &ft,soilc(ft),s1(ft),year,mnth
!      stop
            soilt = 0.97d0*soilt + 0.03d0*tmp(mnth,day)

*            print*,'daily soilc 1:3',soilc(1:3)
            CALL DOLYDAY(ftsla(ft),ftc3(ft),ftphen(ft),ftagh(ft),
     &ftdth(ft),ftlls(ft),ftsls(ft),ftrls(ft),ftbbm(ft),
     &ftbb0(ft),ftbbmax(ft),ftbblim(ft),ftssm(ft),ftsss(ft),ftsslim(ft),
     &ftrat(ft),lat,dep,tmp(mnth,day),prc(mnth,day),hum(mnth,day),
     &cld(mnth),ca,soilc(ft),soiln(ft),minnv,s1(ft),s2(ft),s3(ft),
     &s4(ft),sn(ft),lsn(ft),adp,sfc,sw,sswc,awl,kd,kx,daygpp,dayra,
     &lai(ft),nppstore(ft),nppstorx(ft),nppstor2(ft),evp(ft),
     &dof(ft),evap,tran,roff,interc,evbs,flow1(ft),flow2(ft),year,mnth,
     &day,pet,laimax(ft),ht(ft),leafv,stemv,rootv,leaflit(ft),
     &stemlit(ft),rootlit(ft),bb(ft),ss(ft),bbgs(ft),dsbb(ft),
     &tmem,maxlai(ft),thty_dys,wtfc,wtwp,leafnpp(ft),stemnpp(ft),
     &rootnpp(ft),yield(ft),ft,resp,smtrig,qdirect,qdiff,suma,
     &tsumam(ft),stemfr(ft),lmor_sc(:,ft),nleaf,chill(ft),dschill(ft),
     &fpr,gsm(ft))
*----------------------------------------------------------------------*

*----------------------------------------------------------------------*
* Set daily memories for output.                                       *
*----------------------------------------------------------------------*
            daily_out(1,ft,mnth,day) = lai(ft)
            daily_out(2,ft,mnth,day) = roff
            daily_out(3,ft,mnth,day) = evap+tran
            daily_out(4,ft,mnth,day) = tran
            daily_out(5,ft,mnth,day) =(leafnpp(ft) + stemnpp(ft) + 
     &rootnpp(ft) + nppstore(ft) - leafold - stemold - rootold - 
     &nppsold)*12.0d0
            daily_out(6,ft,mnth,day) = daygpp
            daily_out(7,ft,mnth,day) = srespm/30.0d0
            daily_out(8,ft,mnth,day) = daily_out(5,ft,mnth,day) - 
     &daily_out(7,ft,mnth,day)
            daily_out(9,ft,mnth,day) = tmp(mnth,day)
            daily_out(10,ft,mnth,day) = prc(mnth,day)
            daily_out(11,ft,mnth,day) = hum(mnth,day)
            daily_out(12,ft,mnth,day) = nppstore(ft)*12.0d0
            daily_out(13,ft,mnth,day) = (s1(ft) + s2(ft) + s3(ft) + 
     &s4(ft) - sw(1) - sw(2) - sw(3) - sw(4))/(sfc(1) + sfc(2) + 
     &sfc(3) + sfc(4) - sw(1) - sw(2) - sw(3) - sw(4))
*            daily_out(13,ft,mnth,day) = 
*     &(s1(ft) - sw(1))/(sfc(1) - sw(1))
            daily_out(14,ft,mnth,day) = pet
            daily_out(15,ft,mnth,day) = interc
            daily_out(16,ft,mnth,day) = evbs
            daily_out(17,ft,mnth,day) = min(1.0d0,s1(ft)/10.0d0/topsl)
            daily_out(18,ft,mnth,day) = s1(ft) + s2(ft) + s3(ft) + 
     &s4(ft)
            daily_out(19,ft,mnth,day) = resp
            daily_out(20,ft,mnth,day) = qdirect
            daily_out(21,ft,mnth,day) = qdiff
            daily_out(22,ft,mnth,day) = nleaf
            daily_out(23,ft,mnth,day) = leaflit(ft) - lflitold
            daily_out(24,ft,mnth,day) = cld(mnth)
            daily_out(25,ft,mnth,day) = gsm(ft)
*            daily_out(25,ft,mnth,day) = fpr
*            daily_out(25,ft,mnth,day) = dayra

            IF ((bb(ft).EQ.day+(mnth-1)*30).AND.(budo(ft).EQ.0))
     &budo(ft) = bb(ft)
            IF ((ss(ft).EQ.day+(mnth-1)*30).AND.(seno(ft).EQ.0))
     &seno(ft) = ss(ft)

*----------------------------------------------------------------------*
* Sumation of evaporation and transpiration.                           *
*----------------------------------------------------------------------*
            evapm(mnth,ft) = evapm(mnth,ft) + evap
            tranm(mnth,ft) = tranm(mnth,ft) + tran
            roffm(mnth,ft) = roffm(mnth,ft) + roff
            petm(mnth,ft) = petm(mnth,ft) + pet
            avmnpet(ft) = avmnpet(ft) + pet

*----------------------------------------------------------------------*
* Sumation of GPP and NPP.                                             *
*----------------------------------------------------------------------*

            photm(mnth,ft) = photm(mnth,ft) + daygpp
            laimnth(mnth,ft) = laimnth(mnth,ft) + 
     &           lai(ft)/real(no_days(year,mnth,thty_dys))

            DO d=1,ftlls(ft)
              leafdp(d,ft)=leafv(d)
            ENDDO
            DO d=1,ftsls(ft)
              stemdp(d,ft)=stemv(d)
            ENDDO
            DO d=1,ftrls(ft)
              rootdp(d,ft)=rootv(d)
            ENDDO

            DO d=1,30
              sm_trig(d,ft) = smtrig(d)
            ENDDO

            DO d=1,360
              sumadp(d,ft) = suma(d)
            ENDDO

          ELSE ! dolydo ?
            DO i=1,25
              daily_out(i,ft,mnth,day) = 0.0d0
            ENDDO
          ENDIF
          ENDDO ! ft loop
        ENDDO                     ! daily loop

*----------------------------------------------------------------------*
*      Monthly Operation                                               *
*----------------------------------------------------------------------*
        DO ft=1,nft
        IF (dolydo(ft).EQ.1) THEN

          h2o = (s1(ft) + s2(ft) + s3(ft) + s4(ft))
          IF (h2o.LT.0.0d0) h2o = 0.1d0
          DO i=1,8
            c0v(i)=c0(i,ft)
            n0v(i)=n0(i,ft)
          ENDDO
          DO i=1,3
            minnv(i)=minn(i,ft)
          ENDDO

*          print*,avmnpet(ft),avmnppt,avmnt
          CALL DOLYMONTH(ts,tc,avmnpet(ft),avmnppt,avmnt,h2o,flow1(ft),
     &flow2(ft),c0v,n0v,minnv,nfix,nci,dslc,drlc,dsln,srespm,lchm,ca,
     &site,year,yr0,yrf,speedc,soilc(ft),soiln(ft),mnth)

          sresp(ft) = sresp(ft) + srespm
          lch(ft) = lch(ft) + lchm
               
          DO i=1,8
            c0(i,ft)=c0v(i)
            n0(i,ft)=n0v(i)
          ENDDO
          DO i=1,3
            minn(i,ft)=minnv(i)
          ENDDO
               
        ENDIF
        ENDDO

*----------------------------------------------------------------------*
*                             End of year loop                         *
*----------------------------------------------------------------------*
      ENDDO                     ! end of monthly loop

      DO ft=1,nft

        IF (dolydo(ft).EQ.1) THEN

          IF (abs(leafnpp(ft) + stemnpp(ft) + rootnpp(ft) + 
     &nppstore(ft) - nppstoreold(ft)).GT.1.0E-6) THEN
            npl(ft) = leafnpp(ft)/(leafnpp(ft) + stemnpp(ft) + 
     &rootnpp(ft) + nppstore(ft) - nppstoreold(ft))*100.0d0
            nps(ft) = stemnpp(ft)/(leafnpp(ft) + stemnpp(ft) + 
     &rootnpp(ft) + nppstore(ft) - nppstoreold(ft))*100.0d0
            npr(ft) = rootnpp(ft)/(leafnpp(ft) + stemnpp(ft) + 
     &rootnpp(ft) + nppstore(ft) - nppstoreold(ft))*100.0d0
          ELSE
            npl(ft) = 0.0d0
            nps(ft) = 0.0d0
            npr(ft) = 0.0d0
            npp(ft) = 0.0d0
          ENDIF
          npp(ft) = (leafnpp(ft) + stemnpp(ft) + rootnpp(ft) + 
     &nppstore(ft) - nppstoreold(ft))*12.0d0
*          print'(4f8.1)',npp(ft),npl(ft),nps(ft),npr(ft)
*          print*,leafnpp(ft)*12.0d0,stemnpp(ft)*12.0d0,
*     &rootnpp(ft)*12.0d0

c     put c0,n0 and minn into vector
          DO i=1,8
            c0v(i)=c0(i,ft)
            n0v(i)=n0(i,ft)
          ENDDO
          DO i=1,3
            minnv(i)=minn(i,ft)
          ENDDO

          CALL SUMDV(tc0,tn0,tminn,ts1,ts2,ts3,ts4,tsn,tlsn,
     &c0v,n0v,minnv,s1(ft),s2(ft),s3(ft),s4(ft),
     &sn(ft),lsn(ft),ftcov(ft))

c     check water cycle closure
          prctot = 0.0d0
          yrtran = 0.0d0
          yrevap = 0.0d0
          yrroff = 0.0d0
          yrpet = 0.0d0
          gpp(ft) = 0.0d0
          DO mnth=1,12
            DO day=1,no_days(year,mnth,thty_dys)
              prctot = prctot + prc(mnth,day)
            ENDDO
            yrtran = yrtran + tranm(mnth,ft)
            yrevap = yrevap + evapm(mnth,ft)
            yrroff = yrroff + roffm(mnth,ft)
            yrpet = yrpet + petm(mnth,ft)
            gpp(ft) = gpp(ft) + photm(mnth,ft)
          ENDDO

          trn(ft) = yrtran
          evt(ft) = yrtran + yrevap
          rof(ft) = yrroff
          fpet(ft) = yrpet

*           IF (abs(prctot+is1+is2+is3+is4+isn+ilsn
*     &-s1(ft)-s2(ft)-s3(ft)-s4(ft)-sn(ft)-lsn(ft)-
*     &(yrevap+yrtran+yrroff)).GT.0.001d0) then
*             WRITE(*,'('' Water model error in doly'',4f8.1,i3,i5)')
*     & prctot,yrevap,yrtran,yrroff,ft,iyear
*           ENDIF

          lai(ft) = laimax(ft)

*     Convert nppstore back to grams
          nppstore(ft) = nppstore(ft)*12.0d0
          nppstorx(ft) = nppstorx(ft)*12.0d0
          nppstor2(ft) = nppstor2(ft)*12.0d0

        ELSE ! not dolydo
            DO i=1,10
              DO mnth=1,12
                DO day=1,no_days(year,mnth,thty_dys)
                  daily_out(i,ft,mnth,day) = 0.0d0
                ENDDO
              ENDDO
            ENDDO
            npp(ft) = 0.0d0
            lai(ft) = 0.0d0
            laimax(ft) = 0.0d0
            evp(ft) = 0.0d0
            dof(ft) = 0.0d0
            npl(ft) = 0.0d0
            nps(ft) = 0.0d0
            npr(ft) = 0.0d0
            rof(ft) = 0.0d0
            trn(ft) = 0.0d0
            sresp(ft) = 0.0d0
            evt(ft) = 0.0d0
            fpet(ft) = 0.0d0
            gpp(ft) = 0.0d0
            lch(ft) = 0.0d0
            bb(ft) = 0
            ss(ft) = 0
            bbgs(ft) = 0
            dsbb(ft) = 0
            nppstore(ft) = 0.0d0
            nppstorx(ft) = 0.0d0
            nppstor2(ft) = 0.0d0

            DO day=1,ftlls(ft)
               leafdp(day,ft) = 0.0d0
            ENDDO
            DO day=1,ftsls(ft)
               stemdp(day,ft) = 0.0d0
            ENDDO
            DO day=1,ftrls(ft)
               rootdp(day,ft) = 0.0d0
            ENDDO
         ENDIF                  ! doly do
      ENDDO                     ! ft's
c not use
c        do ft=1,nft
c        ans2(ft) = 0.0d0
c        DO day=1,ftlls(ft)
c          ans2(ft) = ans2(ft) + leafdp(day,ft)
c        ENDDO
c        enddo

        npp(1) = 0.0d0
        lai(1) = 0.0d0

        CALL GROWTH(nft,ftmor,ftwd,ftxyl,ftpd,ftgr0,ftgrf,cov,bio,
     &nppstore,npp,lai,nps,npr,evp,slc,rlc,sln,rln,stembio,rootbio,ppm,
     &hgt,leaflit)

        CALL SWAP(ic0,in0,iminn,is1,is2,is3,is4,isn,ilsn,tc0,tn0,
     &tminn,ts1,ts2,ts3,ts4,tsn,tlsn)

*----------------------------------------------------------------------*
* Average outputs by cover proportions.                                *
*----------------------------------------------------------------------*
	avlai   = 0.0d0
	avnpp   = 0.0d0
      avnppst = 0.0d0
	avdof   = 0.0d0
	avrof   = 0.0d0
	avtrn   = 0.0d0
	avsresp = 0.0d0
	avevt   = 0.0d0
	avpet   = 0.0d0
        avgpp   = 0.0d0
        avlch   = 0.0d0
        avyield = 0.0d0
        DO ft=1,nft
          avlai   = avlai   + ftcov(ft)*lai(ft)
          avnpp   = avnpp   + ftcov(ft)*npp(ft)
          avnppst = avnppst + ftcov(ft)*nppstore(ft)
          avdof   = avdof   + ftcov(ft)*dof(ft)
          avrof   = avrof   + ftcov(ft)*rof(ft)
          avtrn   = avtrn   + ftcov(ft)*trn(ft)
          avsresp = avsresp + ftcov(ft)*sresp(ft)
          avevt   = avevt   + ftcov(ft)*evt(ft)
          avpet   = avpet   + ftcov(ft)*fpet(ft)
          avgpp   = avgpp   + ftcov(ft)*gpp(ft)
          avlch   = avlch   + ftcov(ft)*lch(ft)
          avyield = avyield + ftcov(ft)*yield(ft)
        ENDDO
*----------------------------------------------------------------------*

        tsoilc = 0.0d0
        tsoiln = 0.0d0
        DO i=1,8
          tsoilc = tsoilc + ic0(i)
          tsoiln = tsoiln + in0(i)
        ENDDO
        tsoiln = tsoiln + iminn(3)
        DO ft=1,nft
          tsoilc = tsoilc + slc(ft) + rlc(ft)
        ENDDO
        swcold = swcnew
        swcnew = is1 + is2 + is3 + is4 + isn + ilsn
        sumbio = 0.0d0
        bioind = 0
        covind = 0
        maxbio = 0.0d0
        maxcov = 0.0d0
        leafper = 0.0d0
        stemper = 0.0d0
        rootper = 0.0d0

        DO ft=1,nft
          bioleaf(ft) = 0.0d0
          DO day=1,ftlls(ft)
            bioleaf(ft) = bioleaf(ft) + 
     &leafdp(day,ft)*12.0d0/ftsla(ft)/18.0d0
          ENDDO
        ENDDO

        avppm = 0.0
        DO ft=1,nft
          bioo(ft) = 0.0d0
          covo(ft) = 0.0d0
          DO i=1,ftmor(ft)
            bioo(ft) = bioo(ft) + (bio(i,1,ft) +
     &bio(i,2,ft) + bioleaf(ft) + nppstore(ft))*cov(i,ft)
            covo(ft) = covo(ft) + cov(i,ft)
            leafper = leafper + npl(ft)*cov(i,ft)
            stemper = stemper + nps(ft)*cov(i,ft)
            rootper = rootper + npr(ft)*cov(i,ft)
            IF (ftmor(ft).gt.50) THEN
              avppm = avppm + ppm(i,ft)*cov(i,ft)
            ENDIF
          ENDDO
          sumbio = sumbio + bioo(ft)
          IF (covo(ft).GT.maxcov) THEN
            maxcov = covo(ft)
            covind = ft
          ENDIF
          IF (bioo(ft).GT.maxbio) THEN
            maxbio = bioo(ft)
            bioind = ft
          ENDIF
        ENDDO

*----------------------------------------------------------------------*
* Write outputs.                                                       *
*----------------------------------------------------------------------*
* General.                                                             *
*----------------------------------------------------------------------*
        IF (iyear.GE.nyears-outyears+1) THEN
          WRITE(21,'('' '',f8.1,$)') avlai
          WRITE(22,'('' '',f8.1,$)') avnpp
          WRITE(23,'('' '',f8.1,$)') tsoilc
          WRITE(24,'('' '',f8.1,$)') tsoiln
          WRITE(25,'('' '',f8.1,$)') avnpp-avsresp
          WRITE(26,'('' '',f8.1,$)') min(swcnew,9999.0d0)
          WRITE(27,'('' '',f8.1,$)') sumbio
          WRITE(28,'('' '',i2,$)')   bioind
          WRITE(29,'('' '',i2,$)')   covind
          WRITE(30,'('' '',f8.1,$)') avdof
          WRITE(31,'('' '',f8.1,$)') avrof
          WRITE(32,'('' '',f8.2,$)') firec
          WRITE(33,'('' '',f8.1,$)') avnppst
          WRITE(34,'('' '',f8.1,$)') stembio
          WRITE(35,'('' '',f8.1,$)') rootbio
          WRITE(36,'('' '',f8.1,$)') leafper
          WRITE(37,'('' '',f8.1,$)') stemper
          WRITE(38,'('' '',f8.1,$)') rootper
          WRITE(39,'('' '',f8.1,$)') avsresp
          WRITE(40,'('' '',f8.1,$)') avevt
          WRITE(41,'('' '',f8.1,$)') avgpp
          WRITE(42,'('' '',f8.3,$)') avlch
          WRITE(43,'('' '',f8.1,$)') yearprc
          WRITE(44,'('' '',f8.1,$)') avnpp-avsresp-firec-avlch-avyield
          WRITE(45,'('' '',f8.1,$)') avtrn
          WRITE(46,'('' '',f8.5,$)') fprob
          WRITE(47,'('' '',f8.2,$)') yeartmp
          WRITE(48,'('' '',f8.2,$)') yearhum
          WRITE(49,'('' '',f8.5,$)') avppm
          WRITE(50,'('' '',f8.1,$)') avpet
        ENDIF

*----------------------------------------------------------------------*
* Write optional cov bio bud sen.                                      *
*----------------------------------------------------------------------*
        iofn = iofngft
        IF (iyear.GE.nyears-outyears2+1) THEN
          DO ft=1,nft
            IF (out_cov) THEN
              iofn = iofn + 1
              WRITE(iofn,'('' '',f8.6,$)') covo(ft)
            ENDIF
            IF (out_bio) THEN
              iofn = iofn + 1
              WRITE(iofn,'('' '',f8.1,$)') bioo(ft)
            ENDIF
            IF (out_bud) THEN
              iofn = iofn + 1
              WRITE(iofn,'('' '',i8,$)') budo(ft)
            ENDIF
            IF (out_sen) THEN
              iofn = iofn + 1
              WRITE(iofn,'('' '',i8,$)') seno(ft)
            ENDIF
          ENDDO
        ENDIF

*----------------------------------------------------------------------*
* Write monthly PIXEL outputs.                                         *
*----------------------------------------------------------------------*
        iofn = 400
*        print*,'thty_dys'
        IF (iyear.GE.nyears-outyears1+1) THEN
          IF ((oymd.EQ.1).OR.(oymd.EQ.3)) THEN
            DO i=1,50
              IF (otagsn(i).EQ.1) THEN 
                DO mnth=1,12
                  ans(mnth,1) = 0.0d0
                ENDDO
                DO ft=1,nft
                  DO mnth=1,12
                    IF (omav(i).eq.1) THEN
                      oscale = 1.0/real(no_days(year,mnth,thty_dys))
                    ELSE
                      oscale = 1.0
                    ENDIF
                    DO day=1,no_days(year,mnth,thty_dys)
                      ans(mnth,1) = ans(mnth,1) + 
     &daily_out(i,ft,mnth,day)*ftcov(ft)*oscale
                    ENDDO
                  ENDDO
                ENDDO
                iofn = iofn + 1
                sum1 = 0.0d0
                IF (omav(i).eq.1) THEN
                  oscale = 1.0d0/12.0d0
                ELSE
                  oscale = 1.0d0
                ENDIF
                DO mnth=1,12
                  sum1 = sum1 + ans(mnth,1)*oscale
                ENDDO
                WRITE(iofn,ofmt(i)) year,(ans(mnth,1),
     &mnth=1,12),sum1
              ENDIF
            ENDDO
          ENDIF

*----------------------------------------------------------------------*
* Write daily PIXEL outputs.                                           *
*----------------------------------------------------------------------*
          IF ((oymd.EQ.2).OR.(oymd.EQ.3)) THEN
            DO i=1,50
              IF (otagsn(i).EQ.1) THEN
                iofn = iofn + 1
                WRITE(iofn,'(i4)') year
                DO mnth=1,12
                  DO day=1,no_days(year,mnth,thty_dys)
                    ans(mnth,day) = 0.0d0
                  ENDDO
                ENDDO
                DO ft=1,nft
                  DO mnth=1,12
                    DO day=1,no_days(year,mnth,thty_dys)
                      ans(mnth,day) = ans(mnth,day) + 
     &daily_out(i,ft,mnth,day)*ftcov(ft)
                    ENDDO
                  ENDDO
                ENDDO
                DO mnth=1,12
                  WRITE(iofn,ofmt(50+i)) (ans(mnth,day),
     &day=1,no_days(year,mnth,thty_dys))
                ENDDO
              ENDIF
            ENDDO
          ENDIF
        ENDIF

        iofn = iofnft
        IF (iyear.GE.nyears-outyears2+1) THEN
*----------------------------------------------------------------------*
* Write monthly SUBPIXEL outputs.                                      *
*----------------------------------------------------------------------*
          IF ((oymdft.EQ.1).OR.(oymdft.EQ.3)) THEN
            DO i=1,50
              IF (otagsnft(i).EQ.1) THEN 
                DO ft=1,nft
                  DO mnth=1,12
                    ans(mnth,1) = 0.0d0
                  ENDDO
                  DO mnth=1,12
                    IF (omav(i).eq.1) THEN
                      oscale = 1.0/real(no_days(year,mnth,thty_dys))
                    ELSE
                      oscale = 1.0
                    ENDIF
                    DO day=1,no_days(year,mnth,thty_dys)
                      ans(mnth,1) = ans(mnth,1) + 
     &daily_out(i,ft,mnth,day)*oscale
                    ENDDO
                  ENDDO
                  iofn = iofn + 1
                  WRITE(iofn,ofmt(i)) year,(ans(mnth,1),
     &mnth=1,12)
                ENDDO
              ENDIF
            ENDDO
          ENDIF

*----------------------------------------------------------------------*
* Write daily SUBPIXEL outputs.                                        *
*----------------------------------------------------------------------*
          IF ((oymdft.EQ.2).OR.(oymdft.EQ.3)) THEN
            DO i=1,50
              IF (otagsnft(i).EQ.1) THEN
                DO ft=1,nft
                  iofn = iofn + 1
                  WRITE(iofn,'(i4)') year
                  DO mnth=1,12
                    DO day=1,no_days(year,mnth,thty_dys)
                      ans(mnth,day) = 0.0d0
                    ENDDO
                  ENDDO
                  DO mnth=1,12
                    DO day=1,no_days(year,mnth,thty_dys)
                      ans(mnth,day) = ans(mnth,day) + 
     &daily_out(i,ft,mnth,day)
                    ENDDO
                  ENDDO
                  DO mnth=1,12
                    WRITE(iofn,ofmt(50+i)) (ans(mnth,day),
     &day=1,no_days(year,mnth,thty_dys))
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDIF
        ENDIF

*----------------------------------------------------------------------*
* Write snapshots
*----------------------------------------------------------------------*
        IF (snp_no.GT.0) THEN
          IF (snp_year.LE.snp_no) THEN
            IF (year.EQ.snpshts(snp_year)) THEN
              fno = 100 + (snp_year-1)*4
              WRITE(fno+1,'(F7.3,F9.3,4500F9.1)') 
     &lat,lon,(((bio(i,j,k),j=1,2),i=1,ftmor(k)),k=1,nft)
              WRITE(fno+2,'(F7.3,F9.3,1500F12.9)') 
     &lat,lon,((cov(i,j),i=1,ftmor(j)),j=1,nft)
              WRITE(fno+3,'(F7.3,F9.3,1500F12.7)') 
     &lat,lon,((ppm(i,j),i=1,ftmor(j)),j=1,nft)
              WRITE(fno+4,'(F7.3,F9.3,1500F8.3)') 
     &lat,lon,((hgt(i,j),i=1,ftmor(j)),j=1,nft)
              snp_year = snp_year + 1
            ENDIF
          ENDIF
        ENDIF

*----------------------------------------------------------------------*
* Carbon at the end of the year.                                       *
*----------------------------------------------------------------------*
        DO ft=1,nft
          bioleaf(ft) = 0.0d0
          DO day=1,ftlls(ft)
            bioleaf(ft) = bioleaf(ft) + 
     &leafdp(day,ft)*12.0d0/ftsla(ft)/18.0d0
          ENDDO
        ENDDO

        ans1 = 0.0d0
        DO ft=1,nft
          DO i=1,ftmor(ft)
            ans1 = ans1 + (bio(i,1,ft) + bio(i,2,ft) + bioleaf(ft) +
     &nppstore(ft))*cov(i,ft)
          ENDDO
          ans1 = ans1 + slc(ft) + rlc(ft)
        ENDDO
        ccheck = ccheck - (ans1 + tc0(1) + 
     &tc0(2) + tc0(3) + tc0(4) + tc0(5) + tc0(6) + tc0(7) + tc0(8) - 
     &(avnpp-avlch-avsresp-firec)) - avyield

*----------------------------------------------------------------------*
* Check carbon and water balance, write to 'DIAG' if any problems.     *
*----------------------------------------------------------------------*
*        IF (year.GT.yr0) THEN
*          IF (abs(yearprc-avrof-avevt+swcold-swcnew)
*     &.GT.1.0e-3)  WRITE(11,'(f12.6,'' water fault, site '',f7.3,f9.3,
*     &i6,'', year '',2i6)')
*     & yearprc-avrof-avevt+swcold-swcnew,lat,lon,site,iyear,yearv(iyear)
*        ENDIF
        IF (abs(ccheck).GT.1.0e-3) THEN
          IF (.NOT.speedc) THEN
            WRITE(11,'(f12.6,'' C budget not closed, site '',f7.3,f9.3,
     &i6,'', year '',2i6)') ccheck,lat,lon,site,iyear,yearv(iyear)
            WRITE(*,'(''check'',3f12.6)') ccheck
          ENDIF
        ENDIF

      ENDDO
*----------------------------------------------------------------------*
*                             End of year loop                         *
*----------------------------------------------------------------------*

*----------------------------------------------------------------------*
*                      End record for default output files             *
*----------------------------------------------------------------------*
        DO i=21,50
          WRITE(i,*)
        ENDDO
        iofn = iofngft
        IF (outyears2.GT.0) THEN
          DO ft=1,nft
            IF (out_cov) THEN
              iofn = iofn + 1
              WRITE(iofn,*)
            ENDIF
            IF (out_bio) THEN
              iofn = iofn + 1
              WRITE(iofn,*)
            ENDIF
            IF (out_bud) THEN
              iofn = iofn + 1
              WRITE(iofn,*)
            ENDIF
            IF (out_sen) THEN
              iofn = iofn + 1
              WRITE(iofn,*)
            ENDIF
          ENDDO
        ENDIF

*----------------------------------------------------------------------*
*                             Write initialisation output              *
*----------------------------------------------------------------------*
        IF (initiseo) THEN
          WRITE(81,'(F7.3,F9.3)') lat,lon
          WRITE(82,'(F7.3,F9.3)') lat,lon
          WRITE(83,'(F7.3,F9.3)') lat,lon
          WRITE(84,'(F7.3,F9.3)') lat,lon
          WRITE(85,'(F7.3,F9.3)') lat,lon

          DO ft=2,nft
            WRITE(81,'(1000E13.6)') 
     &((max(bio(i,j,ft),0.0),j=1,2),i=1,ftmor(ft))
            WRITE(82,'(1000E13.6)') (max(cov(i,ft),0.0),i=1,ftmor(ft))
            WRITE(83,'(1000E13.6)') (max(ppm(i,ft),0.0),i=1,ftmor(ft))
            WRITE(84,'(1000E13.6)') (max(hgt(i,ft),0.0),i=1,ftmor(ft))
            WRITE(85,'(1000E13.6)') (max(wdt(i,ft),0.0),i=1,ftmor(ft))
          ENDDO

          WRITE(86,'(F7.3,F9.3)') lat,lon
          WRITE(87,'(F7.3,F9.3)') lat,lon
          WRITE(88,'(F7.3,F9.3)') lat,lon
          DO ft=2,nft
            WRITE(86,'(2000E13.6)' ) 
     &(max(leafdp(i,ft),0.0),i=1,ftlls(ft))
            WRITE(87,'(1000E13.6)') 
     &(max(stemdp(i,ft),0.0),i=1,ftsls(ft))
            WRITE(88,'(1000E13.6)') 
     &(max(rootdp(i,ft),0.0),i=1,ftrls(ft))
          ENDDO

          WRITE(89,'(F7.3,F9.3)') lat,lon
          WRITE(89,'(100E16.8)') (lai(ft),ft=1,nft)
          WRITE(89,'(100E16.8)') (npp(ft),ft=1,nft)
          WRITE(89,'(100E16.8)') (nps(ft),ft=1,nft)
          WRITE(89,'(100E16.8)') (npl(ft),ft=1,nft)
          WRITE(89,'(100E16.8)') (evp(ft),ft=1,nft)
          WRITE(89,'(100E16.8)') (dof(ft),ft=1,nft)
          WRITE(89,'(100E16.8)') (slc(ft),ft=1,nft)
          WRITE(89,'(100E16.8)') (rlc(ft),ft=1,nft)
          WRITE(89,'(100E16.8)') (sln(ft),ft=1,nft)
          WRITE(89,'(100E16.8)') (rln(ft),ft=1,nft)
          WRITE(89,'(100E16.8)') (nppstore(ft),ft=1,nft)
          WRITE(89,'(100E16.8)') (nppstorx(ft),ft=1,nft)
          WRITE(89,'(100E16.8)') (nppstor2(ft),ft=1,nft)
          WRITE(89,'(100I4)   ') (bb(ft),ft=1,nft)
          WRITE(89,'(100I4)   ') (bbgs(ft),ft=1,nft)
          WRITE(89,'(100I4)   ') (dsbb(ft),ft=1,nft)
          WRITE(89,'(100E16.8)') (maxlai(ft),ft=1,nft)
          WRITE(89,'(100E16.8)') (tmem(ft),ft=1,200)
          DO ft=1,nft
            WRITE(89,'(100E16.8)') (sm_trig(i,ft),i=1,30)
          ENDDO
          DO ft=3,nft
            WRITE(89,'(360E16.8)') (sumadp(i,ft),i=1,360)
          ENDDO
          WRITE(89,'(100E16.8)') (stemfr(ft),ft=1,nft)
        ENDIF

        WRITE(80,'(F7.3,F9.3)') lat,lon
        WRITE(80,'(100E16.8)') (s1(ft),ft=1,nft)
        WRITE(80,'(100E16.8)') (s2(ft),ft=1,nft)
        WRITE(80,'(100E16.8)') (s3(ft),ft=1,nft)
        WRITE(80,'(100E16.8)') (s4(ft),ft=1,nft)
        WRITE(80,'(100E16.8)') (sn(ft),ft=1,nft)
        WRITE(80,'(100E16.8)') (lsn(ft),ft=1,nft)
        WRITE(80,'(1008E16.8)') ic0
        WRITE(80,'(1008E16.8)') in0
        WRITE(80,'(1003E16.8)') iminn
        WRITE(80,'(100E16.8)') dslc
        WRITE(80,'(100E16.8)') drlc
        WRITE(80,'(100E16.8)') dsln
        WRITE(80,'(100E16.8)') drln
        WRITE(80,'(100E16.8)') avlai
        WRITE(80,'(100E16.8)') avnpp

      ELSE

        WRITE(11,*) '                  clm stt ssc blk wfs dep lus'
        WRITE(11,'(f7.3,f9.3,1x,20L4)') lat,lon,l_clim,l_stats,
     &l_soil(1),l_soil(3),l_soil(5),l_soil(8),l_lu

*----------------------------------------------------------------------*
* End of climate data exists 'if' statement.                           *
*----------------------------------------------------------------------*
      ENDIF

*----------------------------------------------------------------------*
* End of site loop                                                     *
*----------------------------------------------------------------------*
      ENDDO

*----------------------------------------------------------------------*
* Open file to record version number, command line, input file and     *
* parameter file.                                                      *
*----------------------------------------------------------------------*
      st1 = stoutput
      OPEN(13,FILE=st1(1:blank(st1))//'/simulation.dat')
      CALL GETARG(0,buff1)
      DO j=1,80
        st1(j:j) = buff1(j:j)
      ENDDO
      st4 = stver
      CALL STRIPB(st4)
      WRITE(13,'(A)') st4(1:28)
      CLOSE(98)
      WRITE(13,*)

      WRITE(13,'(''************************************************'')')
      WRITE(13,'(''* Command line arguments                       *'')')
      WRITE(13,'(''************************************************'')')

      DO i=0,narg
        CALL GETARG(i,buff1)
        DO j=1,80
          st1(j:j) = buff1(j:j)
        ENDDO
        WRITE(13,'(A)') st1(1:blank(st1))
      ENDDO
      WRITE(13,*)

      WRITE(13,'(''************************************************'')')
      WRITE(13,'(''* Parameter file (param.dat)                   *'')')
      WRITE(13,'(''************************************************'')')

      OPEN(98,file='param.dat',STATUS='OLD')
80    READ(98,'(A)',end=90) st1
      WRITE(13,'(A)') st1(1:last_blank(st1))
      GOTO 80
90    CONTINUE 
      CLOSE(98)
      WRITE(13,*)

      WRITE(13,'(''************************************************'')')
      WRITE(13,'(''* Input file                                   *'')')
      WRITE(13,'(''************************************************'')')

      CALL GETARG(1,buff1)
*      buff1='..\..\lindsay2\lind1.dat'

      OPEN(98,file=buff1,STATUS='OLD')
65    READ(98,'(A)',end=70) st1
      WRITE(13,'(A)') st1(1:last_blank(st1))
      GOTO 65
70    CONTINUE 
      CLOSE(98)

      WRITE(13,'(''************************************************'')')
*----------------------------------------------------------------------*

      CLOSE(11)
      CLOSE(12)
      CLOSE(13)
      CLOSE(14)
      DO i=21,50
        CLOSE(i)
      ENDDO

      iofn = iofngft
      IF (outyears2.GT.0) THEN
        DO ft=1,nft
          IF (out_cov) THEN
            iofn = iofn + 1
            CLOSE(iofn)
          ENDIF
          IF (out_bio) THEN
            iofn = iofn + 1
            CLOSE(iofn)
          ENDIF
          IF (out_bud) THEN
            iofn = iofn + 1
            CLOSE(iofn)
          ENDIF
          IF (out_sen) THEN
            iofn = iofn + 1
            CLOSE(iofn)
          ENDIF
        ENDDO
      ENDIF

      DO i=70,89
        CLOSE(i)
      ENDDO
      DO i=91,93
        CLOSE(i)
      ENDDO
      DO i=98,99
        CLOSE(i)
      ENDDO

      IF (snp_no.GT.0) THEN
        DO i=1,snp_no
          fno = 100+(i-1)*4
          CLOSE(fno+1)
          CLOSE(fno+1)
          CLOSE(fno+3)
          CLOSE(fno+4)
        ENDDO
      ENDIF


      STOP
      END

*----------------------------------------------------------------------*
*                          File I/O table                              *
*----------------------------------------------------------------------*
*  1:
*  2:
*  3:
*  4:
*  5:
*  6:
*  7:
*  8:
*  9:
* 10:
* 11: diag: Diagnostics file
* 12: site: Site information
* 13: simulation.dat: Record of the command, version number, input file
* and parameter file.
* 14: sites: List of land sites from the land sea mask.
* 15:
* 16:
* 17:
* 18:
* 19:
* 20:
* 21: lai: LAI
* 22: npp: Net Primary Productivity (g C /m^2/yr)
* 23: scn: Soil carbon (g/m^2)
* 24: snn: Soil nitrogen (g/m^2)
* 25: nep: Net Ecosystem Productivity (g C /m^2/yr)
* 26: swc: Average soil water content (mm)
* 27: biot: Total biomass (g C /m^2)
* 28: bioind: Dominant ft in terms of biomass
* 29: covind: Dominant ft in terms of cover
* 30: dof: Average number of days that veg is without leaves (days/yr)
* 31: rof: Runoff (mm/yr)
* 32: fcn: Burnt carbon (g C /m^2/yr)
* 33: nppstore: Stored npp (g C /m^2/yr)
* 34: stembio: Stem biomass (g C /m^2)
* 35: rootbio: Root biomass (g C /m^2)
* 36: leafper: Percentage of npp going to leaf production/maintenance
* 37: stemper: Percentage of npp going to stem production/maintenance
* 38: rootper: Percentage of npp going to root production/maintenance
* 39: sresp: Soil respiration (g C /m^2/yr)
* 40: evt: Evapotranspiration (mm/yr)
* 41: gpp: Gross Primary Productivity (g C /m^2/yr)
* 42: lch: Leached soil carbon (g C /m^2/yr)
* 43: prc: Yearly precipitation (mm/year)
* 44: nbp: Net Biosphere production (g C /m^2/yr)
* 45: trn: Transpiration (mm/yr)
* 46: fab: fraction of area burnt (0-1)
* 47: tmp: Average yearly temperature (degrees C)
* 48: hum: Average yearly temperature (degrees C)
* 49:
* 50:
* 51:
* 52:
* 53:
* 54:
* 55:
* 56:
* 57:
* 58:
* 59:
* 60:
* 61:
* 62:
* 63:
* 64:
* 65:
* 66:
* 67:
* 68:
* 69:
* 70: init: Input initialisation file
* 71: initbio: Input biomass state
* 72: initcov: Input cover state
* 73: initppm: Input plants per metre squared
* 74: inithgt: Input height
* 75: initwdt: Input width of stems
* 76: initleaf: Input leaf biomass by age in days
* 77: initstem: Input stem biomass by age in days
* 78: initroot: Input root biomass by age in days
* 79: initmisc: Input miscelaneous state parameters
* 80: init: Output initialisation file
* 81: initbio: Output biomass state
* 82: initcov: Output cover state
* 83: initppm: Output plants per metre squared
* 84: inithgt: Output height
* 85: initwdt: Output width of stems
* 86: initleaf: Output leaf biomass by age in days
* 87: initstem: Output stem biomass by age in days
* 88: initroot: Output root biomass by age in days
* 89: initmisc: Output miscelaneous state parameters
* 90:
* 91: Used to open temporary input files. 
* 92: Used to open temporary input files.
* 93: Used to open temporary input files.
* 94:
* 95:
* 96:
* 97:
* 98: Used to open temporary input files.
* 99: Used to open temporary input files.
*101-200: Snapshots of the system state.
*201-400: cov_'tag': Fractional coverage per ft (0-1)
*         bio_'tag': Biomass per ft (g C /m^2)
*         bud_'tag': Biomass per ft (g C /m^2)
*         sen_'tag': Biomass per ft (g C /m^2)
*401-: Output files for monthly/daily/pixel/subpixel variables specified
*in the input file.
*----------------------------------------------------------------------*


