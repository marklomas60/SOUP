      SUBROUTINE GENERATE_MONTHLY(ststats,yr0,yrf,tmpv,humv,prcv,
     &seed1,seed2,seed3)
      IMPLICIT NONE
      REAL*8 tmpv(500,12,31),humv(500,12,31),prcv(500,12,31)
      REAL*8 mtmpv(500,12),mhumv(500,12),mprcv(500,12),lat,lon
      INTEGER yr0,yrf,seed1,seed2,seed3,nbrain,blank,du,i
      CHARACTER ststats*1000
      CHARACTER fname4*1000,fname5*1000,fname6*1000
      LOGICAL l_stats
      REAL*8 mstdtmpv(12),mstdhumv(12),mraindays(12),FLOATINGMEAN
      REAL*8 FLOATINGMEAN2

c seed for the random generator of the weather generator
      INTEGER kode,round,year,mnth,day

c Weather generator functions
      REAL*8 TEMP_DAY_MEAN2_FN
      REAL*8 HUM_DAY_MEAN2_FN

c local variables
      REAL*8 rain(31),sumprc,correction,ans
      REAL*8 tmp_yesterday, hum_yesterday,xx(4,4)
      REAL*8 tmp_autocorr_month,hum_autocorr_month

      lat=0.0
      lon=0.0
      du=1

      DO year=1,yrf-yr0+1
        DO mnth=1,12
          mprcv(year,mnth) = prcv(year,mnth,1)
          mtmpv(year,mnth) = tmpv(year,mnth,1)*10.0
          mhumv(year,mnth) = humv(year,mnth,1)*10.0
        ENDDO
      ENDDO

*----------------------------------------------------------------------*
* Read in stats.
*----------------------------------------------------------------------*
      WRITE(fname4,'(100a)') (ststats(i:i),i=1,blank(ststats)),
     &'/stdtmp.dat'
      WRITE(fname5,'(100a)') (ststats(i:i),i=1,blank(ststats)),
     &'/stdhum.dat'
      WRITE(fname6,'(100a)') (ststats(i:i),i=1,blank(ststats)),
     &'/raindays.dat'

c read WG calibration in files for stdtmp, stdhum and raindays
      CALL READ_WG_FILE(du,fname4,ststats,mstdtmpv,lon,lat,l_stats)
      CALL READ_WG_FILE(du,fname5,ststats,mstdhumv,lon,lat,l_stats)
      CALL READ_WG_FILE(du,fname6,ststats,mraindays,lon,lat,l_stats)
*----------------------------------------------------------------------*

c this value should be read in a file. It seems to be calibrated for UK.
      tmp_autocorr_month=0.65d0
      hum_autocorr_month=0.65d0

c Call the weather generator
c Generate 30 days a month, as read in the EX_CLIM subroutine
c
c initialise the yesterday temperature to the mean temperature of the 
c first month, first year
      tmp_yesterday=0.1d0*mtmpv(1,1)
      hum_yesterday=0.1d0*mhumv(1,1)

      DO year=yr0,yrf

*        sumprc = 0.0
*        do mnth=1,12
*          sumprc = sumprc + real(mprcv(year-yr0+1,mnth))
*        enddo
*        print*,year,sumprc

        DO mnth=1,12
          IF (mprcv(year-yr0+1,mnth).LT.0.0d0) THEN
            IF (mnth.GT.1) THEN
              mprcv(year-yr0+1,mnth)=mprcv(year-yr0+1,mnth-1)
            ELSE
              mprcv(year-yr0+1,mnth)=mprcv(year-yr0+1,12)
            ENDIF
          ENDIF
c generate the rain amount 
          CALL RAIN_DAY_SUB_FULLY_NORMALISED(mraindays(mnth), 
     &1.0d0*mprcv(year-yr0+1,mnth),30, seed1, seed2, seed3, rain)

c copie into xprcv
          sumprc=0.0d0
          nbrain=0
          DO day=1,30
            sumprc=sumprc+rain(day)
            prcv(year-yr0+1,mnth,day)=min(10.0d0*rain(day),30000.0d0)
            IF (rain(day).GT.1.0E-3) nbrain=nbrain+1
          ENDDO
             
c generate the temperaure
          correction=0.0d0
c          print *,'stdtmp=',mstdtmpv(mnth)
          DO day=1,30
            tmp_yesterday=TEMP_DAY_MEAN2_FN(0.1d0*
     &FLOATINGMEAN2(day,mnth,year-yr0+1,yrf-yr0+1,mtmpv),
     &mstdtmpv(mnth),tmp_autocorr_month,tmp_yesterday, 
     &seed1,seed2,seed3,1)

*            print'(3f8.2)',
*     &FLOATINGMEAN(day,mnth,year-yr0+1,yrf-yr0+1,mtmpv),
*     &mtmpv(year-yr0+1,mnth),tmp_yesterday

            tmpv(year-yr0+1,mnth,day)=100.0d0*tmp_yesterday
            correction=correction+real(tmpv(year-yr0+1,mnth,day))
          ENDDO

          correction=correction/30.0d0-(100.0d0*0.1d0)*
     &real(mtmpv(year-yr0+1,mnth))
c reajust the mean temperature
          ans = 0.0
          DO day=1,30
            tmpv(year-yr0+1,mnth,day)=tmpv(year-yr0+1,mnth,day)-
     &correction
            ans = ans + tmpv(1,mnth,day)
          ENDDO

c generate the humidity
          correction=0.0d0
          DO day=1,30
            hum_yesterday=HUM_DAY_MEAN2_FN(0.1d0*
     &FLOATINGMEAN2(day,mnth,year-yr0+1,yrf-yr0+1,mhumv),
     &mstdhumv(mnth),hum_autocorr_month,hum_yesterday, 
     &seed1,seed2,seed3,1)
            humv(year-yr0+1,mnth,day)=100.0d0*hum_yesterday
            correction=correction+real(humv(year-yr0+1,mnth,day))
          ENDDO
          correction=correction/30.0d0-(100.0d0*0.1d0)*
     &real(mhumv(year-yr0+1,mnth))
          DO day=1,30
            humv(year-yr0+1,mnth,day)=humv(year-yr0+1,mnth,day)-
     &correction
          ENDDO
        ENDDO

      ENDDO

      DO mnth=1,12
        do day=1,30
*          print*,tmpv(1,mnth,day)/10.0,mtmpv(1,mnth)
        enddo
      ENDDO

      RETURN
      END



*----------------------------------------------------------------------*
c     Interface between ForestETP Weather Generator and SDGVMd
c     Written by Ghislain Picard 30/06/03
c
*----------------------------------------------------------------------*
      SUBROUTINE EX_CLIM_WEATHER_GENERATOR(stinput,ststats,lat,lon,
     &xlatf,xlatres,xlatresn,xlon0,xlonres,xlonresn,yr0,yrf,tmpv,humv,
     &prcv,mcldv,isite,year0,yearf,du,seed1,seed2,seed3,l_clim,l_stats)

      IMPLICIT NONE

*----------------------------------------------------------------------*
      REAL*8 lat,lon,xlon0,xlatf,xlatres,xlonres
c daily data
      REAL*8 tmpv(500,12,31),humv(500,12,31),prcv(500,12,31)      
c monthly data m=monthly
      INTEGER*2 mtmpvv(4,4,500,12),mhumvv(4,4,500,12),mprcvv(4,4,500,12)
      INTEGER*2 mcldvv(4,4,500,12),row,col,cld_default
      REAL*8 mtmpv(500,12),mhumv(500,12),mprcv(500,12),mcldv(500,12)
c monthly parameter for the weather generator
      REAL*8 mstdtmpv(12),mstdhumv(12),mraindays(12)
c function
      REAL*8 FLOATINGMEAN2

c seed for the random generator of the weather generator
      INTEGER seed1,seed2,seed3,kode,round

c Weather generator functions
      REAL*8 TEMP_DAY_MEAN2_FN
      REAL*8 HUM_DAY_MEAN2_FN

c local variables
      REAL*8 rain(31),sumprc,correction,ans
      REAL*8 tmp_yesterday, hum_yesterday,xx(4,4)
      REAL*8 tmp_autocorr_month,hum_autocorr_month

      INTEGER year,year0,yearf,recn,ans2(1000),siteno(4,4),i,du
      INTEGER nyears,yr0,mnth,day,isite,yrf,blank,recl1,recl2
      INTEGER xlatresn,xlonresn,fno,nbrain,indx(4,4),ii,jj
      CHARACTER fname1*1000,fname2*1000,fname3*1000,stinput*1000
      CHARACTER fname4*1000,fname5*1000,fname6*1000,fname7*1000
      CHARACTER ststats*1000
      LOGICAL l_clim,l_stats

      REAL*8 rrow,rcol,ynorm,xnorm

      l_stats = .TRUE.

      cld_default = 57

      IF (du.eq.1) THEN
        recl1 = 4*12
        recl2 = 6*xlonresn
      ELSE
        recl1 = 4*12 + 1
        recl2 = 6*xlonresn + 1
      ENDIF

      nyears = yearf - year0 + 1

*----------------------------------------------------------------------*
* Find the real row col corresponding to lat and lon.                  *
*----------------------------------------------------------------------*

      rrow = (xlatf - lat)/xlatres
      rcol = (lon - xlon0)/xlonres

      ynorm = rrow - real(int(rrow))
      xnorm = rcol - real(int(rcol))
*----------------------------------------------------------------------*

      fno = 90

      OPEN(fno+1,file=stinput(1:blank(stinput))//'/maskmap.dat',
     &ACCESS='DIRECT',RECL=recl2,FORM='formatted',STATUS='old')

*----------------------------------------------------------------------*
* Find the four sites closest to the desired.                          *
*----------------------------------------------------------------------*
      DO ii=1,4
        DO jj=1,4
          row = int(rrow)+jj-1
          col = int(rcol)+ii-1
          IF ((row.GE.1).AND.(row.LE.xlatresn).AND.(col.GE.1).AND.
     &(col.LE.xlonresn)) THEN
            recn = row
            READ(fno+1,'(720i6)',REC=recn) (ans2(i),i=1,xlonresn)
            siteno(ii,jj) = ans2(col)
            IF (siteno(ii,jj).gt.0) THEN
              indx(ii,jj) =  1
            ELSE
              indx(ii,jj) =  0
            ENDIF
          ELSE
            indx(ii,jj) = -1
          ENDIF
        ENDDO
      ENDDO
      CLOSE(fno+1)
*----------------------------------------------------------------------*

      isite = 0
      IF ((indx(2,2).EQ.1).OR.(indx(2,3).EQ.1).OR.(indx(3,2).EQ.1).OR.
     &(indx(3,3).EQ.1)) THEN

        l_clim = .TRUE.

        isite = 1

        WRITE(fname1,'(100a)') (stinput(i:i),i=1,blank(stinput)),
     &'/tmp.dat'
        WRITE(fname2,'(100a)') (stinput(i:i),i=1,blank(stinput)),
     &'/hum.dat'
        WRITE(fname3,'(100a)') (stinput(i:i),i=1,blank(stinput)),
     &'/prc.dat'
        WRITE(fname7,'(100a)') (stinput(i:i),i=1,blank(stinput)),
     &'/cld.dat'

        WRITE(fname4,'(100a)') (ststats(i:i),i=1,blank(ststats)),
     &'/stdtmp.dat'
        WRITE(fname5,'(100a)') (ststats(i:i),i=1,blank(ststats)),
     &'/stdhum.dat'
        WRITE(fname6,'(100a)') (ststats(i:i),i=1,blank(ststats)),
     &'/raindays.dat'

        OPEN(fno+1,file=fname1,access='direct',recl=recl1,
     &form='formatted',status='old',iostat=kode)
        IF (kode.NE.0) THEN
          WRITE(*,*) 'Climate data file does not exist.'
          WRITE(*,*) 'or record mismatch.'
          WRITE(*,'(''"'',A,''"'')') fname1(1:blank(fname1))
          STOP
        ENDIF

        OPEN(fno+2,file=fname2,access='direct',recl=recl1,
     &form='formatted',status='old',iostat=kode)
        IF (kode.NE.0) THEN
          WRITE(*,*) 'Climate data file does not exist.'
          WRITE(*,*) 'or record mismatch.'
          WRITE(*,'(''"'',A,''"'')') fname2(1:blank(fname2))
          STOP
        ENDIF

        OPEN(fno+3,file=fname3,access='direct',recl=recl1,
     &form='formatted',status='old',iostat=kode)
        IF (kode.NE.0) THEN
          WRITE(*,*) 'Climate data file does not exist.'
          WRITE(*,*) 'or record mismatch.'
          WRITE(*,'(''"'',A,''"'')') fname3(1:blank(fname3))
          STOP
        ENDIF

        OPEN(fno+7,file=fname7,access='direct',recl=recl1,
     &form='formatted',status='old',iostat=kode)
*        IF (kode.NE.0) THEN
*          WRITE(*,*) 'No cloud file found using default value'
*          WRITE(*,*) cld_default
*          WRITE(*,'(''"'',A,''"'')') fname7(1:blank(fname7))
*        ENDIF

        DO ii=1,4
          DO jj=1,4
            IF (indx(ii,jj).EQ.1) THEN
              IF (du.eq.1) THEN
                DO year=yr0,yrf
                  READ(fno+1,'(12i4)',REC=(siteno(ii,jj)-1)*nyears+year
     &-year0+1) (mtmpvv(ii,jj,year-yr0+1,mnth),mnth=1,12)
                  READ(fno+2,'(12i4)',REC=(siteno(ii,jj)-1)*nyears+year
     &-year0+1) (mhumvv(ii,jj,year-yr0+1,mnth),mnth=1,12)
                  READ(fno+3,'(12i4)',REC=(siteno(ii,jj)-1)*nyears+year
     &-year0+1) (mprcvv(ii,jj,year-yr0+1,mnth),mnth=1,12)
                  IF (kode.EQ.0) THEN
                  READ(fno+7,'(12i4)',REC=(siteno(ii,jj)-1)*nyears+year
     &-year0+1) (mcldvv(ii,jj,year-yr0+1,mnth),mnth=1,12)
                  ELSE
                    DO mnth=1,12
                      mcldvv(ii,jj,year-yr0+1,mnth) = cld_default*10
                    ENDDO
                  ENDIF
                ENDDO
              ELSE
                DO year=yr0,yrf
                  READ(fno+1,'(12i4)',REC=(siteno(ii,jj)-1)*nyears+year
     &-year0+1) (mtmpvv(ii,jj,year-yr0+1,mnth),mnth=1,12)
                  READ(fno+2,'(12i4)',REC=(siteno(ii,jj)-1)*nyears+year
     &-year0+1) (mhumvv(ii,jj,year-yr0+1,mnth),mnth=1,12)
                  READ(fno+3,'(12i4)',REC=(siteno(ii,jj)-1)*nyears+year
     &-year0+1) (mprcvv(ii,jj,year-yr0+1,mnth),mnth=1,12)
                  IF (kode.EQ.0) THEN
                  READ(fno+7,'(12i4)',REC=(siteno(ii,jj)-1)*nyears+year
     &-year0+1) (mcldvv(ii,jj,year-yr0+1,mnth),mnth=1,12)
                  ELSE
                    DO mnth=1,12
                      mcldvv(ii,jj,year-yr0+1,mnth) = cld_default*10
                    ENDDO
                  ENDIF
                ENDDO
              ENDIF
            ELSE
              DO year=yr0,yrf
                DO mnth=1,12
                  mtmpvv(ii,jj,year-yr0+1,mnth) = 0
                  mhumvv(ii,jj,year-yr0+1,mnth) = 0
                  mprcvv(ii,jj,year-yr0+1,mnth) = 0
                  mcldvv(ii,jj,year-yr0+1,mnth) = 0
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDDO

*----------------------------------------------------------------------*
* Interpolate climate.                                                 *
*----------------------------------------------------------------------*
        DO year=yr0,yrf
          DO mnth=1,12
            DO ii=1,4
              DO jj=1,4
                xx(ii,jj) = real(mtmpvv(ii,jj,year-yr0+1,mnth))
              ENDDO
            ENDDO
            CALL bi_lin(xx,indx,xnorm,ynorm,ans)
            mtmpv(year-yr0+1,mnth) = ans

            DO ii=1,4
              DO jj=1,4
                xx(ii,jj) = real(mhumvv(ii,jj,year-yr0+1,mnth))
              ENDDO
            ENDDO
            CALL bi_lin(xx,indx,xnorm,ynorm,ans)
            mhumv(year-yr0+1,mnth) = ans

            DO ii=1,4
              DO jj=1,4
                xx(ii,jj) = real(mprcvv(ii,jj,year-yr0+1,mnth))
              ENDDO
            ENDDO
            CALL bi_lin(xx,indx,xnorm,ynorm,ans)
            mprcv(year-yr0+1,mnth) = ans

            DO ii=1,4
              DO jj=1,4
                xx(ii,jj) = real(mcldvv(ii,jj,year-yr0+1,mnth))
              ENDDO
            ENDDO
            CALL bi_lin(xx,indx,xnorm,ynorm,ans)
            mcldv(year-yr0+1,mnth) = ans
          ENDDO
        ENDDO

*----------------------------------------------------------------------*

c     check the humidity bound
c        DO year=yr0,min(yrf,yr0+cycle-1)
c           DO mnth=1,12
c check that the mhumv value are ok. This is a dirty hack to detect the
c bad endianness.
c                 if ((mhumv(year-yr0+1,mnth).NE.-299)
c     &                .AND.(mhumv(year-yr0+1,mnth).NE.-311)
c     &                .AND.(mhumv(year-yr0+1,mnth).NE.-111)
c     &                .AND.(mhumv(year-yr0+1,mnth).NE.-211)
c     &                .AND.(mhumv(year-yr0+1,mnth).NE.-37)
c     &            ) THEN
c                    IF ((mhumv(year-yr0+1,mnth).LT.0)
c     &              .OR.(mhumv(year-yr0+1,mnth).GT.10000)) THEN
c               print *, 'humidity out of range. Check ENDIAN mhumv=',
c     &                      mhumv(year-yr0+1,mnth)
c                    ENDIF
c                 ENDIF
c           ENDDO
c        ENDDO

        CLOSE(fno+1)
        CLOSE(fno+2)
        CLOSE(fno+3)
        CLOSE(fno+7)

c read WG calibration in files for stdtmp, stdhum and raindays
       CALL READ_WG_FILE(du,fname4,ststats,mstdtmpv,lon,lat,l_stats)
       CALL READ_WG_FILE(du,fname5,ststats,mstdhumv,lon,lat,l_stats)
       CALL READ_WG_FILE(du,fname6,ststats,mraindays,lon,lat,l_stats)

c this value should be read in a file. It seems to be calibrated for UK.
       tmp_autocorr_month=0.65d0
       hum_autocorr_month=0.65d0

c Call the weather generator
c Generate 30 days a month, as read in the EX_CLIM subroutine
c
c initialise the yesterday temperature to the mean temperature of the 
c first month, first year
       tmp_yesterday=0.1d0*mtmpv(1,1)
       hum_yesterday=0.1d0*mhumv(1,1)

       DO year=yr0,yrf

*          sumprc = 0.0
*          do mnth=1,12
*            sumprc = sumprc + real(mprcv(year-yr0+1,mnth))
*          enddo
*          print*,year,sumprc

          DO mnth=1,12
             IF (mprcv(year-yr0+1,mnth).LT.0.0d0) THEN
                IF (mnth.GT.1) THEN
                   mprcv(year-yr0+1,mnth)=mprcv(year-yr0+1,mnth-1)
                ELSE
                   mprcv(year-yr0+1,mnth)=mprcv(year-yr0+1,12)
                ENDIF
             ENDIF
c generate the rain amount 
             CALL RAIN_DAY_SUB_FULLY_NORMALISED(mraindays(mnth), 
     1            1.0d0*mprcv(year-yr0+1,mnth),  
     2            30, seed1, seed2, seed3, rain)

c copie into xprcv
             sumprc=0.0d0
             nbrain=0
             DO day=1,30
                sumprc=sumprc+rain(day)
             prcv(year-yr0+1,mnth,day)=min(10.0d0*rain(day),30000.0d0)
                IF (rain(day).GT.1.0E-3) nbrain=nbrain+1
             ENDDO
             
c generate the temperaure
             correction=0.0d0
c             print *,'stdtmp=',mstdtmpv(mnth)
             DO day=1,30
                tmp_yesterday=TEMP_DAY_MEAN2_FN(0.1d0*
     1           FLOATINGMEAN2(day,mnth,year-yr0+1,yrf-yr0+1,mtmpv),
     2           mstdtmpv(mnth),
     3           tmp_autocorr_month, tmp_yesterday, 
     4           seed1, seed2,seed3, 1)
*                print'(3f8.2)',
*     &FLOATINGMEAN2(day,mnth,year-yr0+1,yrf-yr0+1,mtmpv),
*     &mtmpv(year-yr0+1,mnth),tmp_yesterday

                tmpv(year-yr0+1,mnth,day)=100.0d0*tmp_yesterday
                correction=correction+real(tmpv(year-yr0+1,mnth,day))
             ENDDO

             correction=correction/30.0d0-(100.0d0*0.1d0)*
     &            real(mtmpv(year-yr0+1,mnth))
c reajust the mean temperature
             ans = 0.0
             DO day=1,30
                tmpv(year-yr0+1,mnth,day)=tmpv(year-yr0+1,mnth,day)-
     &correction
c                print *,'tmpv=',tmpv(year-yr0+1,mnth,day),
c     &               mtmpv(year-yr0+1,mnth)
                ans = ans + tmpv(1,mnth,day)
             ENDDO

c generate the humidity
             correction=0.0d0
             DO day=1,30
                hum_yesterday=HUM_DAY_MEAN2_FN(0.1d0*
     1           FLOATINGMEAN2(day,mnth,year-yr0+1,yrf-yr0+1,mhumv),
     2           mstdhumv(mnth),
     3           hum_autocorr_month, hum_yesterday, 
     4           seed1, seed2,seed3, 1)
                humv(year-yr0+1,mnth,day)=100.0d0*hum_yesterday
                correction=correction+real(humv(year-yr0+1,mnth,day))
             ENDDO

             correction=correction/30.0d0-(100.0d0*0.1d0)*
     &            real(mhumv(year-yr0+1,mnth))
             DO day=1,30
                humv(year-yr0+1,mnth,day)=humv(year-yr0+1,mnth,day)-
     &correction
                ans = ans + humv(year-yr0+1,mnth,day)
             ENDDO
          ENDDO

       ENDDO

      ELSE
        l_clim = .FALSE.
      ENDIF


      RETURN
      END

*----------------------------------------------------------------------*
      SUBROUTINE READ_WG_FILE(du,filename,ststats,values,lon,lat,
     &l_stats)
*----------------------------------------------------------------------*
      CHARACTER filename*1000,st1*1000,st2*1000,st3*1000,ststats*1000
      REAL*8 lon,lat,latr,lonr,values(12),statsv(4,4,12),xx(4,4),latf
      REAL*8 lon0,xnorm,ynorm,ans,rrow,rcol
      INTEGER stcmp,recl1,du,mnth,blank,ii,jj,i,indx(4,4),latn,lonn
      INTEGER recn,row,col,kode
      LOGICAL l_stats

      IF (du.eq.1) THEN
        recl1 = 6*12
      ELSE
        recl1 = 6*12+1
      ENDIF

      OPEN(99,file=ststats(1:blank(ststats))//'/readme.dat',
     &status='old',iostat=kode)
      IF (kode.NE.0) THEN
        WRITE(*,*) 'Climate statistics data file does not exist.'
        WRITE(*,'(''"'',A,''/readme.dat"'')') ststats(1:blank(ststats))
        STOP
      ENDIF

      READ (99,'(A)') st1
      st2='UNIFORM'
      st3='LONLAT'

      IF (stcmp(st1,st2).EQ.1) THEN
        CLOSE(99)
        OPEN(99,file=filename,status='old',iostat=kode)
        IF (kode.NE.0) THEN
          WRITE(*,*) 'Climate statistics data file does not exist.'
          WRITE(*,'(''"'',A,''"'')') 
     &filename(1:blank(filename))
          STOP
        ENDIF
        READ(99,*) (values(mnth),mnth=1,12)
        CLOSE(99)
      ELSEIF (stcmp(st1,st3).EQ.1) THEN
        READ(99,*)
        READ(99,*) latf,lon0
        READ(99,*)
        READ(99,*) latr,lonr
        READ(99,*)
        READ(99,*) latn,lonn
        CLOSE(99)

*----------------------------------------------------------------------*
* Find the real row col corresponding to lat and lon.                  *
*----------------------------------------------------------------------*
        rrow = (latf - lat)/latr
        rcol = (lon - lon0)/lonr

        ynorm = rrow - real(int(rrow))
        xnorm = rcol - real(int(rcol))

*----------------------------------------------------------------------*
        OPEN(99,FILE=filename,STATUS='old',FORM='formatted',
     &ACCESS='direct',RECL=recl1)

        DO ii=1,4
          DO jj=1,4
            row = int(rrow)+jj-1
            col = int(rcol)+ii-1
            IF ((row.GE.1).AND.(row.LE.latn).AND.(col.GE.1).AND.
     &(col.LE.lonn)) THEN
              recn = (row-1)*lonn + col
              READ(99,'(12f6.2)',rec=recn) (statsv(ii,jj,i),i=1,12)
              IF (statsv(ii,jj,1).LT.-0.001d0) THEN
                indx(ii,jj) = 0
              ELSE
                indx(ii,jj) = 1
              ENDIF
            ELSE
              indx(ii,jj) = -1
            ENDIF
          ENDDO
        ENDDO

        IF ((indx(2,2).NE.1).AND.(indx(2,3).NE.1).AND.(indx(3,2).NE.1)
     &.AND.(indx(3,3).NE.1)) THEN
          l_stats = .FALSE.
        ENDIF

        DO i=1,12
          DO ii=1,4
            DO jj=1,4
              xx(ii,jj) = statsv(ii,jj,i)
             ENDDO
          ENDDO
          CALL bi_lin(xx,indx,xnorm,ynorm,ans)
          values(i) = ans
        ENDDO

        CLOSE(99)

      ELSE
        WRITE(*,*) 'First line of stats files must be UNIFORM or LONLAT'
        WRITE(*,*)  filename(1:blank(filename))
        STOP
      ENDIF


      RETURN
      END

*----------------------------------------------------------------------*
      FUNCTION FLOATINGMEAN (day,mnth,year,lastyear,array)
*----------------------------------------------------------------------*
      INTEGER day,mnth,year,lastyear,y,m
      REAL*8 array(500,12),floatingmean
      REAL*8 mean
      
      IF (day.le.15) THEN
        IF (mnth.GT.1) THEN
          y=year
          m=mnth-1
        ELSE
          IF (year.GT.1) THEN
            y=year-1
            m=12
          ELSE
            y=year
            m=1
          ENDIF
        ENDIF
        mean=(real(day)*(1.0d0*array(year,mnth))+
     &real(15-day)*(1.0d0*array(y,m)))/15.0d0
      ELSE
        IF (mnth.LT.12) THEN
          y=year
          m=mnth+1
        ELSE
          IF (year.LT.lastyear) THEN
            y=year+1
            m=1
          ELSE
            y=year
            m=12
          ENDIF
        ENDIF
        mean=(real(30-day)*(1.0d0*array(year,mnth))+
     &real(day-15)*(1.0d0*array(y,m)))/15.0d0

      ENDIF
      floatingmean=mean


      END

*----------------------------------------------------------------------*
      FUNCTION FLOATINGMEAN2 (day,mnth,year,lastyear,array)
*----------------------------------------------------------------------*
      INTEGER day,mnth,year,lastyear
      REAL*8 array(500,12),floatingmean2
      REAL*8 x0,x,xf,y,xx0,xxf
      INTEGER iyear,imnth

      IF (mnth.GT.1) THEN
        x0 = array(year,mnth-1)
      ELSE
        IF (year.GT.1) THEN
          x0 = array(year-1,12)
        ELSE
          x0 = array(year,12)
        ENDIF
      ENDIF
      x  = array(year,mnth)
      IF (mnth.LT.12) THEN
        xf = array(year,mnth+1)
      ELSE
        IF (year.LT.lastyear) THEN
          xf = array(year+1,1)
        ELSE
          xf = array(year,1)
        ENDIF
      ENDIF

      xx0 = (x0+x)/2.0
      xxf = (xf+x)/2.0
      y = 2.0*x - (xx0 + xxf)/2.0

      IF (day.LE.15) THEN
        floatingmean2 = xx0 + real(day-1)*(y-xx0)/14.0
      ELSE
        floatingmean2 = y + real(day-16)*(xxf-y)/14.0
      ENDIF

*      print'(i3,5f8.2)',day,floatingmean2,x,x0,xf,y
*      if (day.eq.10) print'(3f8.2)',x0,x,xf

      END

*----------------------------------------------------------------------*
      SUBROUTINE RAIN_DAY_SUB_FULLY_NORMALISED(NO_RAIN_DAYS,
     1   RAIN_IN_MONTH, DAYS_IN_MONTH, SEED1, SEED2, SEED3, RAIN)
*----------------------------------------------------------------------*
      IMPLICIT NONE

c** The routine simulates daily rainfall for all days in a month
c** The typical number of rainy-days and average rainfall for the month
c** Are the only specific Meteorology inputs. Returned is an array
c** (1-dimensional; day) with the simulated rainfall for each day
c** in it.
      LOGICAL WETDAY(31)
      REAL*8 RAIN(31)
      REAL*8 FRACT_WET, RAIN_PER_DAY, RAIN_IN_MONTH, NO_RAIN_DAYS
      REAL*8 ZERO
      REAL*8 PROB_WET_WET, PROB_WET_DRY, GAMMA, MU
      REAL*8 GENE_RAIN_AMT, CORRECTION
      INTEGER*4 DAYS_IN_MONTH, IDAY
      INTEGER*4 SEED1, SEED2, SEED3, IP, IP0, CLOOP     
      REAL*8  RNDF, RAND01 , XLO, XHI, gamma2
      REAL*8 GAMMADIST_FN


      IF (RAIN_IN_MONTH.LT.1.0e-6) THEN
         DO IDAY=1,DAYS_IN_MONTH
            RAIN(IDAY)=0.0d0
         ENDDO
         RETURN
      ENDIF
      IF (NO_RAIN_DAYS.LT.1.0d0) NO_RAIN_DAYS=1.0d0
c
c **    determine the proportion of the month with wet-days
      FRACT_WET = DBLE(NO_RAIN_DAYS)/DBLE(DAYS_IN_MONTH)
c
c**     find the average rainfall per day with rainfall
      RAIN_PER_DAY = RAIN_IN_MONTH/DBLE(NO_RAIN_DAYS)
c
c**   find transitional probability of a wet-day following a dry-day
c**   of the month (Geng et al 1986)
      PROB_WET_DRY = 0.75*FRACT_WET
c
c**     probability of wet-day following a wet-day (Geng et al 1986)
      PROB_WET_WET = 0.25+PROB_WET_DRY

c
c**     Gamma distribution paramaters (Geng et al 1986)
c       GAMMA = -2.16+1.83*RAIN_PER_DAY  
        XLO = 0.0
        XHI=1.0

        RAND01 = RNDF(XLO,XHI, SEED1, SEED2, SEED3)
        GAMMA = 1.0
c
c**     Mu paramaters (Geng et al 1986)
        gamma2 = gamma+(rand01-0.5)/2.0
        MU = RAIN_PER_DAY/GAMMA2
c
c     first of the month; the chances of it being wet is 'fract_wet'
c     get a random number proportional to the number of days in the 
c     month

      GENE_RAIN_AMT=0.0
      XLO = 0.0
      XHI=1.0
      RAND01 = RNDF(XLO,XHI, SEED1, SEED2, SEED3)
      DO 20 IDAY = 1, DAYS_IN_MONTH
      IF(IDAY.EQ.1) THEN
c** is the random number<=rain_fract; if so then it rained yesterday
      IP0 = 0
      IF(RAND01.LE.FRACT_WET) THEN
        IP0 = 1
      ENDIF
      ENDIF
      RAIN(IDAY) = 0.0
      RAND01 = RNDF(XLO, XHI, SEED1, SEED2, SEED3)
      IF(IP0.EQ.1) THEN
c**       raining on day 0
      IF(RAND01-PROB_WET_WET .LE. 0) IP = 1
      IF(RAND01-PROB_WET_WET .GT. 0) IP = 0
      ELSE   ! IP0 = 0 ; its dry on day 0
      IF(RAND01-PROB_WET_DRY .LE. 0) IP = 1
      IF(RAND01-PROB_WET_DRY .GT. 0) IP = 0
      ENDIF
      IP0 = IP

        IF(IP.EQ.1) THEN     ! Its raining today
           WETDAY(IDAY)=.TRUE.
        ELSE
           WETDAY(IDAY)=.FALSE.
        ENDIF
      IF(WETDAY(IDAY)) THEN     ! Its raining today
         CLOOP=0
12       CONTINUE
      ZERO = 0.0

c        TR1 = EXP(-18.42*MU)

c      TR2 = EXP(-18.42*(1-MU))


C Changed GAMMA to GAMMA2 14/12/01 Paul Henshall.
       
        RAIN(IDAY)=GAMMADIST_FN(zero,MU,GAMMA2,SEED1,SEED2,SEED3)
        CLOOP = CLOOP + 1
        IF (CLOOP.GT.10) THEN
           print *,'Problem in RAIN_DAY_SUB_FULLY_NORMALISED'
           print *,NO_RAIN_DAYS,RAIN_IN_MONTH, DAYS_IN_MONTH
           stop
        ENDIF
        if(rain(iday).gt.RAIN_PER_DAY*5) goto 12
      ENDIF
        GENE_RAIN_AMT=GENE_RAIN_AMT+RAIN(IDAY)
20    CONTINUE

c normalised the rain amount
c      print *,RAIN_IN_MONTH
      IF (GENE_RAIN_AMT.GT.0.0d0) THEN
         CORRECTION = RAIN_IN_MONTH/GENE_RAIN_AMT
     
         DO IDAY = 1, DAYS_IN_MONTH
            RAIN(IDAY)=RAIN(IDAY)*CORRECTION
         ENDDO
      ELSE
         IDAY=RNDF(XLO,XHI, SEED1, SEED2, SEED3)*DAYS_IN_MONTH+1
         IF (IDAY.GT.DAYS_IN_MONTH) IDAY=DAYS_IN_MONTH
         RAIN(IDAY)=RAIN_IN_MONTH
      ENDIF

      RETURN
      END



c
c   this subroutine comes from the RAIN_DAY_SUB subroutine in 
c   ForestETP Weather Generator it differs by the fact that the wet
c   days are firstly generated. Then, knowing the number of rain days
c   in the month really generated (with respect to average number of
c   wet days), the rain amount is generated. The goal is to reduce the
c   variability in the rain amount for one year to the other, since the
c   measured monthly mean are known and input of the model.


      SUBROUTINE RAIN_DAY_SUB_NORMALISED(NO_RAIN_DAYS, RAIN_IN_MONTH,  
     1  DAYS_IN_MONTH, SEED1, SEED2, SEED3, RAIN)
      IMPLICIT NONE



c** The routine simulates daily rainfall for all days in a month
c** The typical number of rainy-days and average rainfall for the month 
c** Are the only specific Meteorology inputs. Returned is an array
c** (1-dimensional; day) with the simulated rainfall for each day 
c** in it.
      LOGICAL WETDAY(31)
      REAL*8 RAIN(31)
      REAL*8 FRACT_WET, RAIN_DAY, RAIN_IN_MONTH, NO_RAIN_DAYS
      REAL*8 ZERO
      REAL*8 PROB_WET_WET, PROB_WET_DRY, GAMMA, MU
      INTEGER*4 DAYS_IN_MONTH, IDAY, GENE_NO_RAIN_DAYS
      INTEGER*4 SEED1, SEED2, SEED3, IP, IP0
      REAL*8  RNDF, RAND01 , XLO, XHI, gamma2
      REAL*8 GAMMADIST_FN

c
c **    determine the proportion of the month with wet-days
      FRACT_WET = DBLE(NO_RAIN_DAYS)/DBLE(DAYS_IN_MONTH)
c
c**     find the average rainfall per day with rainfall
cccccc calculate later      RAIN_DAY = RAIN_IN_MONTH/DBLE(NO_RAIN_DAYS)
c
c**     find transitional probability of a wet-day following a dry-day
c**     of the month (Geng et al 1986)
      PROB_WET_DRY = 0.75*FRACT_WET
c
c**     probability of wet-day following a wet-day (Geng et al 1986)
      PROB_WET_WET = 0.25+PROB_WET_DRY
c
c**     Gamma distribution paramaters (Geng et al 1986) 
c       GAMMA = -2.16+1.83*RAIN_DAY   
        XLO = 0.0
        XHI=1.0

        RAND01 = RNDF(XLO,XHI, SEED1, SEED2, SEED3)
        GAMMA = 1.0
c
c**     Mu paramaters (Geng et al 1986) 
        gamma2 = gamma+(rand01-0.5)/2.0
cccccc calculate later      MU = RAIN_DAY/GAMMA2
c
c     first of the month; the chances of it being wet is 'fract_wet'
c     get a random number proportional to the number of days in the 
c     month

      GENE_NO_RAIN_DAYS=0
      XLO = 0.0
      XHI=1.0
      RAND01 = RNDF(XLO,XHI, SEED1, SEED2, SEED3)
      DO 20 IDAY = 1, DAYS_IN_MONTH
	IF(IDAY.EQ.1) THEN
c**   is the random number<=rain_fract; if so then it rained yesterday
	  IP0 = 0
	  IF(RAND01.LE.FRACT_WET) THEN
	    IP0 = 1 
	  ENDIF
	ENDIF
	RAIN(IDAY) = 0.0
	RAND01 = RNDF(XLO, XHI, SEED1, SEED2, SEED3)
	IF(IP0.EQ.1) THEN
c**       raining on day 0
	  IF(RAND01-PROB_WET_WET .LE. 0) IP = 1
	  IF(RAND01-PROB_WET_WET .GT. 0) IP = 0
	ELSE   ! IP0 = 0 ; its dry on day 0
	  IF(RAND01-PROB_WET_DRY .LE. 0) IP = 1
	  IF(RAND01-PROB_WET_DRY .GT. 0) IP = 0
	ENDIF
	IP0 = IP

        IF(IP.EQ.1) THEN     ! Its raining today
           WETDAY(IDAY)=.TRUE.
           GENE_NO_RAIN_DAYS=GENE_NO_RAIN_DAYS+1
        ELSE
           WETDAY(IDAY)=.FALSE.
        ENDIF
20    CONTINUE
c
c
c     Now, the generated number of wet days in the month is know
c     Let's calculate the rain amount
c

      RAIN_DAY = RAIN_IN_MONTH/DBLE(GENE_NO_RAIN_DAYS)
      MU = RAIN_DAY/GAMMA2

      DO 120 IDAY = 1, DAYS_IN_MONTH
	IF(WETDAY(IDAY)) THEN     ! Its raining today
112       CONTINUE
*112       TR1 = EXP(-18.42*MU)

*	  TR2 = EXP(-18.42*(1-MU))

c         RAND01 = RNDF(XLO, XHI, SEED1, SEED2, SEED3)
c         IF(RAND01 - TR1 .LE. 0) THEN
c           S1 = 0
c         ELSE
c           S1 = RAND01*EXP(1/MU)
c         ENDIF
c         RAND01 = RNDF(XLO, XHI, SEED1, SEED2, SEED3)
c         IF(RAND01 - TR2 .LE. 0) THEN
c           S2 = 0
c         ELSE
c           S2 = RAND01*EXP(1/(1-MU))
c         ENDIF
c         IF(S1+S2-1 .LE. 0) THEN
c           IF( (S1+S2).EQ.0 )THEN
c             Z = 0.0
c           ELSE
c             Z = S1/(S1+S2)
c           ENDIF
c         ELSE
c           FLAG = 1
c         ENDIF
c         IF(FLAG.EQ.1) GOTO  30
c
c         RAND01 = RNDF(XLO, XHI, SEED1, SEED2, SEED3)
c         RAIN(IDAY) = -Z*LOG(RAND01)*GAMMA
	  ZERO = 0.0

C Changed GAMMA to GAMMA2 14/12/01 Paul Henshall.
        
   	  RAIN(IDAY)=GAMMADIST_FN(zero,MU,GAMMA2,SEED1,SEED2,SEED3)
        if(rain(iday).gt.rain_day*5)goto 112
	ENDIF
120    CONTINUE
c
c//! Equation: 32
c//! Description: The routine simulates daily rainfall for all days in \
c//!  a month. The typical number of rainy-days and average rainfall \
c//!  for the month are the only specific Meteorology inputs. Returned \
c//!  is an array (1-dimensional; day) with the simulated \
c//!  rainfall for each day in it. Although the routine uses  \
c//!  Markov-chain principles, the start of each month is INDEPENDENT \
c//!  of what the conditions were on the last day of the previous  \
c//!  month. \\
c//!  The routine calls a random-number gernerator (RNDF) which \
c//!  returns numbers between Xmin and Xmax (in this case 0,1) using \
c//!  the 3 seeds \\
c//!   INPUTS: Number of rainy days in the month; Average rainfall of\
c//!    the month; number of days in the month; 3 seeds \
c//!    for random number generation\\
c//!   OUTPUT: array (dimension:31), with rainfall for each day
c//! Bibliography: S.Evans (SWELTER)
c
c//! Case: Day_{wet,dry}^{i} = 0
c//!  E: Rain_{day}^{i} = 0.0 \\
c//! \\
c//! Case: Day_{wet,dry}^{i} = 1 \\
c//! \\
c//!  E: Rain_{day}^{i} = Z \gamma log_{e} (Rnd_{j1,j2,j3})  \\
c//!   \gamma = -2.16+1.83 Rain_{day,mean}  \\
c//! \\
c//!  Case: S_{1} + S_{2} -1 <= 0    
c//!   E: Z = { {{S_{1}} \over{S_{1}+S_{2}} }  \\
c//! \\
c//!     Case: Rnd_{j1,j2,j3} - Tr_{1} <= 0
c//!      E: S_{1} = 0
c//!     Case: Rnd_{j1,j2,j3} - Tr_{1} > 0
c//!      E: S_{1} = Rnd_{j1,j2,j3} exp(1/\mu) \\
c//! \\
c//!     Case: Rnd_{j1,j2,j3} - Tr_{2} <= 0
c//!      E: S_{2} = 0
c//!     Case: Rnd_{j1,j2,j3} - Tr_{2} > 0
c//!      E: S_{2} = Rnd_{j1,j2,j3} exp(1/(1-\mu))   \\
c//! \\
c//!  Case: S_{1} + S_{2} -1 > 0
c//!     E: re-generate S_{1} and S_{2}  \\
c//!    \\ 
c//!  Tr_{1} = {exp(-18.42)}\over{\mu} \\
c//!  Tr_{2} = {exp(-18.42)}\over{(1-mu)} \\
c//!  \mu = {{Rain_{day,mean}} \over{\gamma}} \\
c//! \\
c//!  Case: Day_{wet,dry}^{i-1} = 1 \\
c//! \\
c//!   Case: Rnd_{j1,j2,j3} - Prob_{wet,wet} <= 0
c//!    E: Day_{wet,dry}^{i} = 1
c//!   Case: Rnd - Prob_{wet,wet} > 0
c//!    E: Day_{wet,dry}^{i} = 0 \\
c//! \\
c//!  Case: Day_{wet,dry}^{i-1} = 0 \\
c//! \\
c//!   Case: Rnd_{j1,j2,j3} - Prob_{wet,dry} <= 0
c//!    E: Day_{wet,dry}^{i} = 1      
c//!   Case: Rnd_{j1,j2,j3} - Prob_{wet,dry} > 0 
c//!    E: Day_{wet,dry}^{i} = 0      \\
c//! \\
c//!  Prob_{wet,wet} = 0.25+Prob{wet,dry} \\
c//!  Prob_{wet,dry} = { {0.75 N_{rain days}} \over{DinM} }
c
c//! Variable: Rain_{day}^{i}
c//! Description: Predicted daily rainfall on day of month, i
c//! Units: mm
c//! Type: REAL*8 array (12,31)
c//! SET
c
c//! Variable: Rain_{day,mean}
c//! Description: Mean daily rainfall for the month
c//! Units: mm/day
c//! Type: REAL*8
c
c//! Variable: Rnd_{j1,j2,j3)
c//! Description: Random number generated between 0 and 1
c//! Units: none
c//! Type: REAL*8
c
c//! Variable: j1
c//! Description: Seed for random number generator; changed on return
c//! Units: none
c//! Type: Integer
c
c//! Variable: j2
c//! Description: Seed for random number generator; changed on return
c//! Units: none
c//! Type: Integer
c
c//! Variable: j3
c//! Description: Seed for random number generator; changed on return
c//! Units: none
c//! Type: Integer
c
c//! Variable: DinM
c//! Description: Number of days in the month
c//! Units: none
c//! Type: Integer
c
      RETURN
      END





c from ForestETP weather generator.
c the function which generate the temperature
c has been tranformed into humidity generator
c this is the original function:
c       DOUBLE PRECISION FUNCTION TEMP_DAY_MEAN2_FN(TEMP_MEAN_MONTH,
c     1  TEMP_SD_MONTH, TEMP_AUTOCORR_MONTH, TEMP_YESTERDAY, IS1,
c     2  IS2, IS3, NORMSWITCH)
c 
c


      DOUBLE PRECISION FUNCTION HUM_DAY_MEAN2_FN(HUM_MEAN_MONTH,
     1  HUM_SD_MONTH, HUM_AUTOCORR_MONTH, HUM_YESTERDAY, IS1,
     2  IS2, IS3, NORMSWITCH)
      
      IMPLICIT NONE
      EXTERNAL NORM_FROM_RAND_FN
c      
c** calculates mean daily temperature  from monthly values   
      DOUBLE PRECISION HUM_MEAN_MONTH
      DOUBLE PRECISION HUM_SD_MONTH, HUM_AUTOCORR_MONTH
      DOUBLE PRECISION HUM_YESTERDAY
      DOUBLE PRECISION NORM, NORMMN, NORMSD
      DOUBLE PRECISION NORM_FROM_RAND_FN
      INTEGER*4 IS1, IS2, IS3, NORMSWITCH
C 
c** approximate a function to get a normal distribution, mean=0, sd=1
      NORMMN=0.0
      NORMSD = 1.0
      NORM=NORM_FROM_RAND_FN(NORMMN, NORMSD, IS1, IS2, IS3, NORMSWITCH)
c
      HUM_DAY_MEAN2_FN = HUM_MEAN_MONTH + HUM_AUTOCORR_MONTH*
     1   (HUM_YESTERDAY - HUM_MEAN_MONTH) + HUM_SD_MONTH*
     2   NORM*( (1-(HUM_AUTOCORR_MONTH)**2)**0.5)
c
c//! Equation: 47
c//! Description: Calculates Average daily humidity from monthly \
c//!  means and standard-deviations. Autocorrelation of 0.65 seems \
c//!  to work if one isn't available. \\
c//!   INPUTS: Mean humidity of the month; Standard deviation \
c//!    around the mean on a daily basis; Humidity \
c//!     auto-correlation for each month; Mean humidity of the \
c//!     previous day; 3 random number seeds.
c//! Bibliography: S Evans (SWELTER); Haith, Tubbs & Pickering (1984)\
c//!  Simulation of Pollution by soil erosion & soil nutirnt loss, \
c//!   PUDOC, Wageningen
c
c//! E: Hum_{mean} = Hum_{mean,month} + \rho_{T} (Hum_{day-1} - \
c//!     Hum_{mean,month}) + \sigma_{mT} N (1-\sigma_{mT}^{2})^0.5) \\
c//!  \\
c//!   E: N = { Random, N(0,1) }
c
c//! Variable: Hum_{mean}
c//! Description: Daily mean humidity
c//! Units: Degrees C   (usually)
c//! Type: Double precsion
c//! SET
c
c//! Variable: Hum_{mean,month}
c//! Description: Mean daily humidity (monthly basis)
c//! Units: Degrees C   (usually)
c//! Type: Double precsion
c
c//! Variable: Hum_{day-1}     
c//! Description: Mean daily humidity of the previous day
c//! Units: Degrees C   (usually)
c//! Type: Double precsion
c
c//! Variable: \rho_{T}
c//! Description: Observed first-order auto-correlation of mean \
c//!  daily air humidity for each month
c//! Units: none  (day/day^{-1})
c//! Type: Double precsion
c
c//! Variable: \sigma_{mT}
c//! Description: Standard deviation of daily air humidity about \
c//!  the mean (monthly basis).
c//! Units: Degrees C             
c//! Type: Double precsion
c
c//! Variable: R
c//! Description: Random number between 0,1
c//! Units: none
c//! Type: Double precsion
c
      RETURN
      END

