*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE COUNTRY                          *
*                          ******************                          *
*                                                                      *
* Returns country name and id from a 1 degree * 1 degree map from      *
* NASA GISAS: Global Land Use Datasets: Country Codes                  *
* http:/www.giss.nasa.gov/data/landuse/coutry.html                     *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE COUNTRY(stmask,lat,lon,country_name,country_id,
     &l_regional)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 lat,lon,adj_lat,adj_lon,xlat,xlon
      CHARACTER country_name*1000,stmask*1000
      INTEGER country_id,ilat,ilon,ans(360),x,blank,i
      LOGICAL l_regional

      ilat = int(91.0d0 - lat)
      ilon = int(lon + 181.0d0)

      OPEN(99,file=stmask(1:blank(stmask))//'/country.dat',status='old')
        IF (ilat.GT.1) THEN
          DO i=1,ilat-1
            READ(99,*)
          ENDDO
        ENDIF
        READ(99,*) ans
        country_id = ans(ilon)
      CLOSE(99)

*----------------------------------------------------------------------*
* If OCEAN was found try the nearest adjacent squares.                 *
*----------------------------------------------------------------------*
      IF (country_id.EQ.0) THEN

* Nearest lateral.

        xlat = lat+500-real(int(lat+500.0d0))
        IF (xlat.GT.0.5d0) THEN
          adj_lat = 1.0d0
        ELSE
          adj_lat = -1.0d0
        ENDIF

        xlon = lon+500-real(int(lon+500.0d0))
        IF (xlon.GT.0.5d0) THEN
          adj_lon = 1.0d0
        ELSE
          adj_lon = -1.0d0
        ENDIF

        IF (abs(xlat-0.5d0).GT.abs(xlon-0.5d0)) THEN
          adj_lon = 0.0d0
        ELSE
          adj_lat = 0.0d0
        ENDIF

        ilat = int(91.0d0 - (lat + adj_lat))
        ilon = int((lon + adj_lat) + 181.0d0)

        IF (ilon.EQ.0) ilon = 360
        IF (ilon.EQ.361) ilon = 1
        IF (ilat.EQ.0) ilat = 1
        IF (ilat.EQ.91) ilat = 90

        OPEN(99,file=stmask(1:blank(stmask))//'/country.dat',
     &status='old')
          IF (ilat.GT.1) THEN
            DO i=1,ilat-1
              READ(99,*)
            ENDDO
          ENDIF
          READ(99,*) ans
          country_id = ans(ilon)
        CLOSE(99)

      ENDIF

      IF (country_id.EQ.0) THEN

* Next nearest lateral.

        IF (abs(adj_lat).LT.0.5d0) THEN
          xlat = lat+500-real(int(lat+500.0d0))
          IF (xlat.GT.0.5d0) THEN
            adj_lat = 1.0d0
          ELSE
            adj_lat = -1.0d0
          ENDIF
        ENDIF

        IF (abs(adj_lon).LT.0.5d0) THEN
          xlon = lon+500-real(int(lon+500.0d0))
          IF (xlon.GT.0.5d0) THEN
            adj_lon = 1.0d0
          ELSE
            adj_lon = -1.0d0
          ENDIF
          adj_lat = 0.0d0
        ENDIF

        IF (abs(adj_lat).GT.0.5d0)  adj_lon = 0.0d0

        ilat = int(91.0d0 - (lat + adj_lat))
        ilon = int((lon + adj_lat) + 181.0d0)

        IF (ilon.EQ.0) ilon = 360
        IF (ilon.EQ.361) ilon = 1
        IF (ilat.EQ.0) ilat = 1
        IF (ilat.EQ.91) ilat = 90

        OPEN(99,file=stmask(1:blank(stmask))//'/country.dat',
     &status='old')
        IF (ilat.GT.1) THEN
          DO i=1,ilat-1
            READ(99,*)
          ENDDO
        ENDIF
        READ(99,*) ans
        country_id = ans(ilon)
        CLOSE(99)

      ENDIF

      IF (country_id.EQ.0) THEN

* Nearest diagonal

        xlat = lat+500-real(int(lat+500.0d0))
        IF (xlat.GT.0.5d0) THEN
          adj_lat = 1.0d0
        ELSE
          adj_lat = -1.0d0
        ENDIF

        xlon = lon+500-real(int(lon+500.0d0))
        IF (xlon.GT.0.5d0) THEN
          adj_lon = 1.0d0
        ELSE
          adj_lon = -1.0d0
        ENDIF

        ilat = int(91.0d0 - (lat + adj_lat))
        ilon = int((lon + adj_lat) + 181.0d0)

        IF (ilon.EQ.0) ilon = 360
        IF (ilon.EQ.361) ilon = 1
        IF (ilat.EQ.0) ilat = 1
        IF (ilat.EQ.91) ilat = 90

        OPEN(99,file=stmask(1:blank(stmask))//'/country.dat',
     &status='old')
          IF (ilat.GT.1) THEN
            DO i=1,ilat-1
              READ(99,*)
            ENDDO
          ENDIF
          READ(99,*) ans
          country_id = ans(ilon)
        CLOSE(99)

      ENDIF

*----------------------------------------------------------------------*
* Regional or country switch.                                          *
      IF (.NOT.(l_regional)) country_id = country_id-mod(country_id,100)
*----------------------------------------------------------------------*

      OPEN(99,file=stmask(1:blank(stmask))//'/country_id.dat',
     &status='old')
10    CONTINUE
        read(99,'(i6,5x,a15)') x,country_name
        IF (x.eq.country_id) GOTO 20
      GOTO 10
20    CONTINUE
      CLOSE(99)


      RETURN
      END

*----------------------------------------------------------------------*
*                          SUBROUTINE CO2_0_F                          *
*                          ******************                          *
*----------------------------------------------------------------------*
      SUBROUTINE CO2_0_F(co20,co2f,yearv,yr0,spinl,co2,co2const,nyears)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      INCLUDE 'array_dims.inc'
      INTEGER iyear,year,spinl,nyears,yr0,yearv(maxyrs)
      REAL*8 co2const,co20,co2f,co2(maxyrs)

      iyear = 1
      year = yearv(iyear)
      IF ((spinl.gt.0).AND.(iyear.GT.spinl)) THEN
        co20 = co2(year-yr0+1)
      ELSE
        IF (co2const.GT.0.0d0) THEN
          co20 = co2const
        ELSE
          co20 = co2(year-yr0+1)
        ENDIF
      ENDIF
      iyear = nyears
      year = yearv(iyear)
      IF ((spinl.gt.0).AND.(iyear.GT.spinl)) THEN
        co2f = co2(year-yr0+1)
      ELSE
        IF (co2const.GT.0.0d0) THEN
          co2f = co2const
        ELSE
          co2f = co2(year-yr0+1)
        ENDIF
      ENDIF


      RETURN
      END

*----------------------------------------------------------------------*
*                          SUBROUTINE SOIL_PARAMETERS                  *
*                          **************************                  *
*----------------------------------------------------------------------*
      SUBROUTINE SOIL_PARAMETERS(ts,tsi,tc,bulk,orgc,wilt,field,sat,dep,
     &sttxdp,lat,lon,soil_chr,soil_chr2,du,l_soil)
      IMPLICIT NONE
      REAL*8 lat,lon,soil_chr(10),soil_chr2(10),ts,tsi,tc,bulk,orgc
      REAL*8 wilt,field,sat,dep
      INTEGER du
      CHARACTER sttxdp*1000
      LOGICAL l_soil(20)

      CALL EX_SOIL(sttxdp,lat,lon,soil_chr2,du,l_soil)
      IF (soil_chr(1).GT.0.001d0) THEN
        ts = soil_chr(1)
        tsi = soil_chr(2)
        l_soil(1) = .TRUE.
        l_soil(2) = .TRUE.
      ELSE
        ts = soil_chr2(1)
        tsi = soil_chr2(2)
      ENDIF
      tc = 100.0d0 - ts - tsi

      IF (soil_chr(3).GT.0.001d0) THEN
        bulk = soil_chr(3)
        l_soil(3) = .TRUE.
      ELSE
        bulk = soil_chr2(3)
      ENDIF

      IF (soil_chr(4).GT.0.000000001d0) THEN
        orgc = soil_chr(4)
        l_soil(4) = .TRUE.
      ELSE
        orgc = soil_chr2(4)
      ENDIF

      IF (soil_chr(5).GT.0.001d0) THEN
        wilt = soil_chr(5)
        field = soil_chr(6)
        sat = soil_chr(7)
        l_soil(5) = .TRUE.
        l_soil(6) = .TRUE.
        l_soil(7) = .TRUE.
      ELSE
        wilt = soil_chr2(5)
        field = soil_chr2(6)
        sat = soil_chr2(7)
      ENDIF

      IF (soil_chr(8).GT.0.001d0) THEN
        dep = soil_chr(8)
        l_soil(8) = .TRUE.
      ELSE
        dep = soil_chr2(8)
      ENDIF


      RETURN
      END

*----------------------------------------------------------------------*
*                          SUBROUTINE LAND_SITE_CHECK                  *
*                          **************************                  *
*----------------------------------------------------------------------*
      SUBROUTINE LAND_SITE_CHECK(st1,st2,sites,lat_lon,latdel,londel,
     &du,stmask,no_countries)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      CHARACTER st1*1000,st2*1000,stmask*1000
      INTEGER sites,du,i,blank,site,n_fields,no_countries
      REAL*8 lat_lon(2000000,2),latdel,londel,lat,lon
      LOGICAL lc

      OPEN(14,FILE=st2(1:blank(st2))//'/land_sites.dat')

      i = 0
      DO site=1,sites
        if (mod(site,100).eq.0) print*,sites,site
        lat = lat_lon(site,1)
        lon = lat_lon(site,2)
        CALL lorc(du,lat,lon,latdel,londel,stmask,lc)
        IF (lc) THEN
          WRITE(14,'(f7.3,f9.3)') lat,lon
          i = i + 1
          lat_lon(i,1) = lat_lon(site,1)
          lat_lon(i,2) = lat_lon(site,2)
        ENDIF
      ENDDO
      CLOSE(14)
      sites = i
      WRITE(*,'(i8,'' potential land sites.'')') sites
      IF (no_countries.EQ.0) THEN
      IF ((n_fields(st1).EQ.3).OR.(n_fields(st1).EQ.5).OR.
     &(n_fields(st1).EQ.7)) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'The site line of the input file has 3, 5 or 7 fileds
     &.'
        STOP
      ENDIF
      ENDIF


      RETURN
      END

*----------------------------------------------------------------------*
*                          SUBROUTINE OUTPUT_OPTIONS                   *
*                          *************************                   *
*----------------------------------------------------------------------*
      SUBROUTINE OUTPUT_OPTIONS(nomdos,otags,omav,ofmt)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      INTEGER nomdos,omav(50)
      CHARACTER otags(100)*1000,ofmt(100)*100

*----------------------------------------------------------------------*
* Read in monthly and/or daily output options.                         *
*----------------------------------------------------------------------*
* file tags.                                                           *
*----------------------------------------------------------------------*
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
      otags(25) = 'gsm'

      nomdos = 25

*----------------------------------------------------------------------*
* Averaged output.                                                     *
*----------------------------------------------------------------------*
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
      omav(25) = 1

*----------------------------------------------------------------------*
* Monthly output format.                                               *
*----------------------------------------------------------------------*
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
      ofmt(25)  ='(i4,1x,12f8.3,f10.6)'

*----------------------------------------------------------------------*
* Daily output format.                                                 *
*----------------------------------------------------------------------*
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
      ofmt(75)  ='(31f7.4)'


      RETURN
      END

*----------------------------------------------------------------------*
*                          SUBROUTINE SET_PIXEL_OUT                    *
*                          ************************                    *
*----------------------------------------------------------------------*
      SUBROUTINE SET_PIXEL_OUT(st1,st2,outyears,nomdos,otagsn,otags,
     &oymd)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      INTEGER otagsn(50),nomdos,i,j,ii,n_fields,oymd,stcmp,outyears
      INTEGER blank,l,ntags
      CHARACTER otags(100)*1000,st1*1000,st2*1000,st3*1000
      CHARACTER st4*1000,st5*1000

      CALL STRIPB(st1)
      ii = n_fields(st1)

      IF (ii.GT.1) THEN
        CALL STRIPBS(st1,st2)
        CALL STRIPB(st1)
        st3 = 'MONTHLY'
        st4 = 'DAILY'
        st5 = 'ALL'
        IF (stcmp(st1,st3).EQ.1) THEN
          oymd = 1
        ELSEIF (stcmp(st1,st4).EQ.1) THEN
          oymd = 2
        ELSEIF (stcmp(st1,st5).EQ.1) THEN
          oymd = 3
        ELSE
          WRITE(*,'('' PROGRAM TERMINATED'')')
          WRITE(*,*) 'The second field of line 15 of the input file must
     & read MONTHLY DAILY or ALL.'
          WRITE(*,'('' "'',A,''"'')') st1(1:blank(st1))
          WRITE(*,*) 'Output variable options:'
          WRITE(*,'(1x,20a4)') (otags(j),j=1,15)
          WRITE(*,'(1x,20a4)') (otags(j),j=16,nomdos)
          STOP
        ENDIF
        CALL STRIPBS(st1,st2)

        ii = n_fields(st1)
        IF (ii.GE.1) THEN
          CALL STRIPBN(st1,outyears)
          IF (outyears.EQ.-9999) THEN
            WRITE(*,'('' PROGRAM TERMINATED'')')
            WRITE(*,*) 'The third field of the ''PIXEL'' line must conta
     &in an integer, if it exists.'
            STOP
          ENDIF
        ELSEIF (ii.EQ.0) THEN
          outyears = 0
        ENDIF
      ELSE
        outyears = 0
        oymd = 0
      ENDIF

      DO i=1,50
        otagsn(i) = 0
      ENDDO
      IF (ii.GT.1) THEN
        CALL STRIPB(st1)
        DO i=1,ii-1
          CALL STRIPBS(st1,st2)
          l = ntags(otags,st2)
          IF (l.NE.-9999) THEN
            otagsn(l) = 1
          ELSE
            WRITE(*,'('' PROGRAM TERMINATED'')')
            WRITE(*,*) 'Error in tag name on the ''PIXEL'' line.'
            WRITE(*,'('' "'',A,''"'')') st2(1:blank(st2))
            WRITE(*,*) 'Available tags:'
            WRITE(*,'(1x,20a4)') (otags(j),j=1,15)
            WRITE(*,'(1x,20a4)') (otags(j),j=16,nomdos)
            STOP
          ENDIF 
       ENDDO
      ELSE
        DO i=1,nomdos
          otagsn(i) = 1
        ENDDO
      ENDIF


      RETURN
      END

*----------------------------------------------------------------------*
*                          SUBROUTINE SET_SUBPIXEL_OUT                 *
*                          ***************************                 *
*----------------------------------------------------------------------*
      SUBROUTINE SET_SUBPIXEL_OUT(st1,st2,outyears,nomdos,
     &otagsnft,otags,oymdft,out_cov,out_bio,out_bud,out_sen)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      INTEGER otagsnft(50),nomdos,i,j,ii,n_fields,stcmp,oymdft,outyears
      INTEGER blank,l,ntags
      CHARACTER otags(100)*1000,st1*1000,st2*1000,st3*1000
      CHARACTER st4*1000,st5*1000
      LOGICAL out_cov,out_bio,out_bud,out_sen

      out_cov = .false.
      out_bio = .false.
      out_bud = .false.
      out_sen = .false.

      CALL STRIPB(st1)
      ii = n_fields(st1)
      IF (ii.GT.1) THEN
        CALL STRIPBS(st1,st2)
        CALL STRIPB(st1)
        st3 = 'MONTHLY'
        st4 = 'DAILY'
        st5 = 'ALL'
        IF (stcmp(st1,st3).EQ.1) THEN
          oymdft = 1
        ELSEIF (stcmp(st1,st4).EQ.1) THEN
          oymdft = 2
        ELSEIF (stcmp(st1,st5).EQ.1) THEN
          oymdft = 3
        ELSE
          WRITE(*,'('' PROGRAM TERMINATED'')')
          WRITE(*,*) 'The second field of line 16 of the input file must
     & read MONTHLY DAILY or ALL.'
          WRITE(*,'('' "'',A,''"'')') st1(1:blank(st1))
          WRITE(*,*) 'Output variable options:'
          WRITE(*,'(1x,20a4)') (otags(j),j=1,15)
          WRITE(*,'(1x,20a4)') (otags(j),j=16,nomdos),
     &'cov ','bio ','bud ','sen '
          STOP
        ENDIF
        CALL STRIPBS(st1,st2)

        ii = n_fields(st1)
        IF (ii.GE.1) THEN
          CALL STRIPBN(st1,outyears)
          IF (outyears.EQ.-9999) THEN
            WRITE(*,'('' PROGRAM TERMINATED'')')
            WRITE(*,*) 'The third field of the ''SUB_PIXEL'' line must c
     &ontain an integer, if it exists.'
            STOP
          ENDIF
        ELSEIF (ii.EQ.0) THEN
          outyears = 0
        ENDIF
      ELSE
        outyears = 0
        oymdft = 0
      ENDIF

      DO i=1,50
        otagsnft(i) = 0
      ENDDO

      IF (ii.GT.1) THEN
        CALL STRIPB(st1)
        DO i=1,ii-1
          CALL STRIPBS(st1,st2)
          st3 = 'cov'
          IF (stcmp(st2,st3).EQ.1) THEN
            out_cov = .true.
          ELSE
            st3 = 'bio'
            IF (stcmp(st2,st3).EQ.1) THEN
              out_bio = .true.
            ELSE
              st3 = 'bud'
              IF (stcmp(st2,st3).EQ.1) THEN
                out_bud = .true.
              ELSE
                st3 = 'sen'
                IF (stcmp(st2,st3).EQ.1) THEN
                  out_sen = .true.
                ELSE
                  l = ntags(otags,st2)
                  IF (l.NE.-9999) THEN
                    otagsnft(l) = 1
                  ELSE
                    WRITE(*,'('' PROGRAM TERMINATED'')')
                    WRITE(*,*) 'Error in tag name in the ''SUBPIXEL'' li
     &ne.'
                    WRITE(*,'('' "'',A,''"'')') st2(1:blank(st2))
                    WRITE(*,*) 'Available tag names:'
                    WRITE(*,'(1x,20a4)') (otags(j),j=1,15)
                    WRITE(*,'(1x,20a4)') (otags(j),j=16,nomdos),
     &'cov ','bio ','bud ','sen '
                    STOP
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ELSE
        DO i=1,nomdos
          otagsnft(i) = 1
        ENDDO
        out_cov = .true.
        out_bio = .true.
        out_bud = .true.
        out_sen = .true.
      ENDIF


      RETURN
      END

*----------------------------------------------------------------------*
*                          SUBROUTINE MKDLIT                           *
*                          *****************                           *
*----------------------------------------------------------------------*
      SUBROUTINE MKDLIT(nft,ftmor,ftcov,dslc,drlc,dsln,drln,cov,slc,rlc,
     &sln,rln)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      INCLUDE 'array_dims.inc'
      REAL*8 ftcov(maxnft),dslc,drlc,dsln,drln,cov(maxage,maxnft)
      REAL*8 slc(maxnft),rlc(maxnft),sln(maxnft),rln(maxnft),sum
      INTEGER nft,ft,j,ftmor(maxnft)

      DO ft=1,nft
        ftcov(ft) = 0.0d0
        DO j=1,ftmor(ft)
          ftcov(ft) = ftcov(ft) + cov(j,ft)
        ENDDO
      ENDDO
      sum = 0.0d0
      DO ft=1,nft
        sum = sum + ftcov(ft)
      ENDDO
      DO ft=1,nft
        ftcov(ft) = ftcov(ft)/sum
      ENDDO

      dslc = 0.0d0
      drlc = 0.0d0
      dsln = 0.0d0
      drln = 0.0d0
      DO ft=1,nft
        dslc = dslc + slc(ft)
        drlc = drlc + rlc(ft)
        dsln = dsln + sln(ft)
        drln = drln + rln(ft)
      ENDDO


      RETURN
      END

*----------------------------------------------------------------------*
*                          SUBROUTINE SUMDV                            *
*                          ****************                            *
*----------------------------------------------------------------------*
      SUBROUTINE SUMDV(tc0,tn0,tminn,ts1,ts2,ts3,ts4,tsn,tlsn,c0,n0,
     &minn,s1,s2,s3,s4,sn,lsn,ftcov)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 tc0(8),tn0(8),tminn(3),ts1,ts2,ts3,ts4,tsn,tlsn,c0(8),
     &n0(8),minn(3),s1,s2,s3,s4,sn,lsn,ftcov
      INTEGER i

      DO i=1,8
        tc0(i) = tc0(i) + c0(i)*ftcov
        tn0(i) = tn0(i) + n0(i)*ftcov
      ENDDO
      DO i=1,3
        tminn(i) = tminn(i) + minn(i)*ftcov
      ENDDO
      ts1 = ts1 + s1*ftcov
      ts2 = ts2 + s2*ftcov
      ts3 = ts3 + s3*ftcov
      ts4 = ts4 + s4*ftcov
      tsn = tsn + sn*ftcov
      tlsn = tlsn + lsn*ftcov


      RETURN
      END

*----------------------------------------------------------------------*
*                          SUBROUTINE SWAP                             *
*                          ***************                             *
*----------------------------------------------------------------------*
      SUBROUTINE SWAP(ic0,in0,iminn,is1,is2,is3,is4,isn,ilsn,c0,
     &n0,minn,s1,s2,s3,s4,sn,lsn)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 ic0(8),in0(8),iminn(3),is1,is2,is3,is4,isn,ilsn
      REAL*8 c0(8),n0(8),minn(3),s1,s2,s3,s4,sn,lsn

      INTEGER i

      DO i=1,8
         ic0(i) = c0(i)
         in0(i) = n0(i)
      ENDDO
      DO i=1,3
        iminn(i) = minn(i)
      ENDDO
      is1 = s1
      is2 = s2
      is3 = s3
      is4 = s4
      isn = sn
      ilsn = lsn


      RETURN
      END


*----------------------------------------------------------------------*
*                          SUBROUTINE SWAP                             *
*                          ***************                             *
*----------------------------------------------------------------------*
      SUBROUTINE MIX_WATER(s1,s2,s3,s4,sn,lsn,ftcov,ftmix,nft)
*----------------------------------------------------------------------*
      INCLUDE 'array_dims.inc'
      REAL*8 s1(maxnft),s2(maxnft),s3(maxnft),s4(maxnft),ftcov(maxnft)
      REAL*8 sn(maxnft),lsn(maxnft),ftmix(maxnft),ans,mix(maxnft)
      REAL*8 share1,share2,share3,share4,share5,share6,sum
      INTEGER nft,ft

      ans = 0.0d0
      DO ft=1,nft
        ans=ans+(s1(ft)+s2(ft)+s3(ft)+s4(ft)+sn(ft)+lsn(ft))*ftcov(ft)
      ENDDO

      DO ft=1,nft
        mix(ft) = 1.0d0 - (1.0d0 - ftmix(ft))**(20.0d0/360.0d0)
      ENDDO

      share1 = 0.0d0
      share2 = 0.0d0
      share3 = 0.0d0
      share4 = 0.0d0
      share5 = 0.0d0
      share6 = 0.0d0
      sum = 0.0d0
      DO ft=1,nft
        share1 = share1 + s1(ft)*ftcov(ft)*mix(ft)
        share2 = share2 + s2(ft)*ftcov(ft)*mix(ft)
        share3 = share3 + s3(ft)*ftcov(ft)*mix(ft)
        share4 = share4 + s4(ft)*ftcov(ft)*mix(ft)
        share5 = share5 + sn(ft)*ftcov(ft)*mix(ft)
        share6 = share6 + lsn(ft)*ftcov(ft)*mix(ft)
        s1(ft) = s1(ft)*ftcov(ft)*(1.0d0 - mix(ft))
        s2(ft) = s2(ft)*ftcov(ft)*(1.0d0 - mix(ft))
        s3(ft) = s3(ft)*ftcov(ft)*(1.0d0 - mix(ft))
        s4(ft) = s4(ft)*ftcov(ft)*(1.0d0 - mix(ft))
        sn(ft) = sn(ft)*ftcov(ft)*(1.0d0 - mix(ft))
        lsn(ft) = lsn(ft)*ftcov(ft)*(1.0d0 - mix(ft))
        sum = sum + ftcov(ft)*mix(ft)
      ENDDO

      ans = 0.0d0
      do ft=1,nft
        ans = ans + s1(ft)+ s2(ft)+ s3(ft)+ s4(ft)+ sn(ft)+ lsn(ft)
      enddo

      ans = 0.0d0
      DO ft=1,nft
        IF (ftcov(ft).GT.0.0d0) THEN
          IF (sum.GT.0.0d0) THEN
            s1(ft) = (s1(ft) + share1*ftcov(ft)*mix(ft)/sum)/ftcov(ft)
            s2(ft) = (s2(ft) + share2*ftcov(ft)*mix(ft)/sum)/ftcov(ft)
            s3(ft) = (s3(ft) + share3*ftcov(ft)*mix(ft)/sum)/ftcov(ft)
            s4(ft) = (s4(ft) + share4*ftcov(ft)*mix(ft)/sum)/ftcov(ft)
            sn(ft) = (sn(ft) + share5*ftcov(ft)*mix(ft)/sum)/ftcov(ft)
            lsn(ft) =(lsn(ft)+ share6*ftcov(ft)*mix(ft)/sum)/ftcov(ft)
          ELSE
            s1(ft) = s1(ft)/ftcov(ft)
            s2(ft) = s2(ft)/ftcov(ft)
            s3(ft) = s3(ft)/ftcov(ft)
            s4(ft) = s4(ft)/ftcov(ft)
            sn(ft) = sn(ft)/ftcov(ft)
            lsn(ft) = lsn(ft)/ftcov(ft)
          ENDIF
        ENDIF
      ENDDO

      ans = 0.0d0
      DO ft=1,nft
        ans=ans+(s1(ft)+s2(ft)+s3(ft)+s4(ft)+sn(ft)+lsn(ft))*ftcov(ft)
      ENDDO


      RETURN
      END
