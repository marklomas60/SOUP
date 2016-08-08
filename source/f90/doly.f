*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE DOLYDAY                          *
*                          ******************                          *
*                                                                      *
*           (D)ynamic gl(O)ba(L) phytogeograph(Y) model                *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE DOLYDAY(sla,c3,ftphen,ftagh,ftdth,leafls,stemls,
     &rootls,bbm,bb0,bbmax,bblim,ssm,sss,sslim,lairat,lat,dep,tmp,prc,
     &hum,cld,ca,soilc,soiln,minn,s1,s2,s3,s4,sn,lsn,adp,sfc,sw,sswc,
     &awl,kd,kx,daygpp,dayra,rlai,nppstore,nppstorx,nppstor2,maxevap,
     &daysoff,evap,tran,roff,interc,evbs,flow1,flow2,year,mnth,day,pet,
     &laimax,ht,leafv,stemv,rootv,leaflit,stemlit,rootlit,bb,ss,bbgs,
     &dsbb,tmem,maxlai,thty_dys,wtfc,wtwp,leafnpp,stemnpp,rootnpp,
     &yield,ft,resp,sm_trig,qdirect,qdiff,suma,tsumam,stemfr,
     &lmor_sc,nleaf,chill,dschill,fpr,gsm)

*----------------------------------------------------------------------*
      REAL*8 oi
      PARAMETER(oi = 21000.0d0)
*----------------------------------------------------------------------*
      REAL*8 bb0,bbmax,bblim,sslim,lat,dep,tmp,prc,hum,cld,ca,maxevap
      REAL*8 lairat,interc,evbs,leafv(3600),stemv(1000),rootv(1000),t
      REAL*8 soilc,soiln,rh,tk,rd,rn,pet,petmm,petv,laimax,q,qdiff
      REAL*8 qdirect,hrs,vpd,canga,lam,rho,s,gam,amx,gsn,dp2,pet2,wtfc
      REAL*8 etwt,etmm,ftagh,soilw,maxc,soil2g,can2a,daygpp,daynpp
      REAL*8 tmem(200),gsum,sapresp,sla,nppstorx,leafmol,ht,tf(13)
      REAL*8 eemm,daysoff,amax,respref,suma(360),can2g,cangs,svp,rlai
      REAL*8 wtwp,p,nppstor2,dayl,sum,canres,et,ee,flow1,flow2,tran,rem
      REAL*8 adp(4),evap,f2,sn,lsn,s1,s2,s3,s4,sfc(4),sw(4),awl(4),kd
      REAL*8 kx,f3,minn(3),roff,windspeed,leaflit,lflit,stemlit,smlit
      REAL*8 rootlit,rtlit,npp_eff,eff_dec,eff_age,eff,nppstore,stemnpp
      REAL*8 rootnpp,leafnpp,tfscale,s_gpp,s_rr,s_rn,s_sr,s_sn
      REAL*8 s_lr,s_ln,sswc(4),maxlai,yield,yld,resp,fpr,dayra,gsm
      REAL*8 sm_trig(30),s1in,tsumam,stemfr,lmor_sc(3600),nleaf
      INTEGER leafls,stemls,rootls,bbm,ssm,sss,ftphen,c3,thty_dys,ft
      INTEGER mnth,i,iter,no_day,ndsum(12),lai,day,year,bb,bbgs
      INTEGER ftdth,ss,dsbb,chill,dschill
      LOGICAL veg
      INCLUDE 'param.inc'

      maxlai = 11.9d0
*      print*,rlai,day,mnth

*      if ((day.eq.1).and.(mnth.eq.1)) print*,tmp,prc,hum,cld

*----------------------------------------------------------------------*
* SET parameter for THROUGHFALL.                                       *
*----------------------------------------------------------------------*
       tf(1) = 1.0d0
       tf(2) = 0.95d0
       tf(3) = 0.935d0
       tf(4) = 0.92d0
       tf(5) = 0.905d0
       tf(6) = 0.89d0
       tf(7) = 0.875d0
       tf(8) = 0.86d0
       tf(9) = 0.845d0
       tf(10)= 0.83d0
       tf(11)= 0.815d0
       tf(12)= 0.80d0
       tf(13)= 0.785d0

       tfscale = 2.75d0
       DO i=1,13
         tf(i) = 1.0d0 - tfscale*(1.0d0 - tf(i))
       ENDDO

*----------------------------------------------------------------------*
* Data input for number of days for each month of the year 'ndays',    *
* the mid julian day of each month 'ndsum'.                            *
*----------------------------------------------------------------------*
      DATA ndsum/16,46,75,106,136,167,197,228,259,289,320,350/

      p = 101325.0d0

      s_gpp = 0.0d0
      s_rr = 0.0d0
      s_rn = 0.0d0
      s_lr = 0.0d0
      s_ln = 0.0d0
      s_sr = 0.0d0
      s_sn = 0.0d0

      rem = rlai - int(rlai)
      lai = int(rlai) + 1

*----------------------------------------------------------------------*
* INITIALISATIONS.                                                     *
*----------------------------------------------------------------------*
* Night-time respiration ? CHECK
      rd = 0.82E-6
* maxc, moisture response parameter
      maxc = 690.0d0/622.6d0
* Leaf molecular weight.
      IF (sla.GT.0.0d0) THEN
        leafmol = 1.0d0/(sla*18.0d0)
      ELSE
        leafmol =0.0d0
      ENDIF

c initialisation a faire a chaque fois
      IF ((ft.EQ.1).OR.(ft.EQ.2)) THEN
         veg = .FALSE.
      ELSE
         veg = .TRUE.
      ENDIF

*----------------------------------------------------------------------*
* Water model variables and rain generater.                            *
*----------------------------------------------------------------------*
      iter = 30

c unsused anymore... so not important if set to this value every day. 
c Better: should be remove from nppcalc
      amx = -1000.0d0
      amax = -1000.0d0

*----------------------------------------------------------------------*
      soilw = s1 + s2 + s3 + s4
      soil2g = soilw/(dep*10.0d0)
      IF (soil2g.GT.1.0d0) soil2g = 1.0d0

*----------------------------------------------------------------------*
* Bare ground IF statement                                             *
*----------------------------------------------------------------------*
      IF (.NOT.(veg)) rlai = 0.0d0
*----------------------------------------------------------------------*

      t = tmp
      tk = 273.0d0 + t
      rh = hum
*      abmin = tmp(mnth,day)*1.29772d0 - 19.5362d0
      IF (rh.GT.95.0d0)  rh=95.0d0
      IF (rh.LT.30.0d0)  rh=30.0d0

* Sapwood respiration N.B. 8.3144 J/K/mol= Universal gas constant.
*     sapresp = exp(21.6d0 - 5.367e4/(8.3144d0*tk))
      sapresp = 0.4d0*exp(0.06d0*t)*0.1d0
      sapresp = 0.4d0*exp(0.02d0*t)*0.35d0

* Converted to a monthly figure (N.B. 24 hour day).
*          sapresp = (sapresp*3600.0d0*24.0d0*30.0d0)/1000000.0d0
* Use daily value
      sapresp = (sapresp*3600.0d0*24.0d0)/1000000.0d0

* Totalled for the year.
      IF (soil2g.GT.wtwp) THEN
        respref = sapresp*p_resp*((soil2g - wtwp)/min(1.0d0,(wtfc - wtwp
     &)))**p_kgw
      ELSE
        respref = 0.0d0
      ENDIF

      IF (t.LT.0.0d0) then
        respref = 0.0d0
      ENDIF

*      if (day.eq.1) print'(i4,6f10.4)',mnth,respref,sapresp,t,
*     &((soil2g - wtwp)/min(1.0d0,(wtfc - wtwp)))**p_kgw

*----------------------------------------------------------------------*
* Calculate daily potential evapotraspiration for each month.          *
*----------------------------------------------------------------------*
      hrs = dayl(lat,no_day(year,mnth,day,thty_dys))
      CALL PFD(lat,no_day(year,mnth,day,thty_dys),hrs,cld,qdirect,qdiff,
     &q)
c          print *,'irrad=',q
c     &      ,pfd_without_cloud(lat,no_day(year,mnth,day,thty_dys),hrs)
*      print'(''direct diff '',i3,2f10.5)',mnth,qdirect,qdiff

      t = tmp
      lam = 2500.0d0 - (2.367d0*t)
      s = 48.7d0*exp(0.0532d0*t)
      rn = 0.96d0*(q*1000000.0d0/4.0d0 + 208.0d0 + 6.0d0*t)
      rn = rn*0.52d0
      pet = (1.26d0*s*rn)/(s + 66.0d0)*p_pet
      petmm = (pet*3600.0d0)/(lam*1000.0d0)*10.0d0
      IF (petmm.GT.0.0d0) THEN
        petv = petmm/1.0d0
      ELSE
        petv = 0.0d0
      ENDIF
      pet = petv

*----------------------------------------------------------------------*


*----------------------------------------------------------------------*
c     canga=k^2 u / (log[(z-d)/z0])^2
c     k=von Karman constant. k=0.41
c     z=reference height
c     d=zero plane displacement
c     z0=roughness length

      windspeed= 5.0d0 ! in m/s
      canga = 0.168d0*windspeed/log((200.0d0 - 
     &0.7d0*ht)/(0.1d0*ht))**2

      npp_eff = 0.0d0
      sum = 0.0d0
      eff_dec = 0.75d0
      eff_age = 90.0d0
      DO i=1,leafls
        IF (i.LE.eff_age) THEN
          eff = 1.0d0
        ELSE
          eff = eff_dec + real(leafls - i)*(1.0d0 - eff_dec)/
     &(real(leafls) - eff_age + 1.0d0) 
        ENDIF
        npp_eff = npp_eff +  leafv(i)*eff
        sum = sum + leafv(i)
      ENDDO
      IF (abs(sum-rlai).gt.0.001) THEN
        WRITE(*,*) 'leafv not = rlai ',SUM,RLAI,mnth,day,nppstore
*        STOP
      ENDIF
      IF (sum.GT.0.0d0) THEN
        npp_eff = npp_eff/sum
      ELSE
        npp_eff = 0.0d0
      ENDIF
*      npp_eff = 1.0d0

      DO i=1,359
        suma(361-i) = suma(360-i)
      ENDDO

      IF (veg) THEN
        CALL NPPCALC(npp_eff,c3,maxc,soilc,soiln,minn,soil2g,wtwp,wtfc,
     &rd,rlai,t,rh,ca,oi,rn,qdirect,qdiff,can2a,can2g,canres,suma(1),
     &amx,amax,gsum,hrs,canga/1.3d0,p,mnth,day,nleaf,fpr,gsm)

c      write(*,*) 'o',canga,can2a,can2g,canres,suma,gsum

*----------------------------------------------------------------------*
* Calculate the canopy conductance.                                    *
*----------------------------------------------------------------------*
        cangs = 8.3144d0*tk/p*can2g
        cangs = 1.6d0*cangs             ! convert CO2 gs into H2O gs
        gsn = cangs
*----------------------------------------------------------------------*
      ELSE
        can2a = 0.0d0
        can2g = 0.0d0
        canres = 0.d0
        cangs = 0.0d0
        gsn = 0.0d0
        suma(1) = 0.0d0
      ENDIF

      dp2 = prc
      pet2 = petv

*----------------------------------------------------------------------*
*                      EVAPOTRANSPIRATION                              *
*----------------------------------------------------------------------*
* Penman-Monteith equation for evapotranspiration.                     *
* Units of ET = W/m2  (CHECK) N.B. Conductances in moles.              *
*----------------------------------------------------------------------*
      svp = 6.108d0*exp((17.269d0*t)/(237.3d0 + t))*100.0d0
      vpd = (1.0d0 - rh/100.0d0)*svp
      lam = 2500.0d0 - 2.367d0*t
      rho = 1288.4d0 - 4.103d0*t
      s = 48.7d0*exp(0.0532d0*t)
      gam = 101325.0d0*1.012d0/(0.622d0*lam)

      IF ((rlai.GT.0.1d0).AND.(soil2g.GT.wtwp)) THEN
        et = (s*rn + rho*1.012d0*canga*vpd)/(s + gam*(1.0d0 +
     &canga/gsn))*p_et
c            print *,'gsn=',gsn
c            print *,'vpd=',vpd
c watch dog ajoute par ghislain 
        IF (et.LT.0.0d0) THEN
*           print *,'et is negativ et=',et
           et=0.0d0
        ENDIF

        etwt = (et*3600.0d0*hrs)/lam
        etmm = etwt/1000.0d0
        IF (etmm.GT.0.0d0) THEN
          IF (etmm/hrs.GT.maxevap)  maxevap = etmm/hrs
        ENDIF
      ELSE
        et = 0.0d0
        etwt = 0.0d0
        etmm = 0.0d0
      ENDIF
      etmm = etmm

      ee = (s*rn + rho*1.012d0*canga*vpd)/(s + gam)
      eemm = (ee*3600.0d0*hrs)/(lam*1000.0d0)

c     added by Ghislain 20/10/03
      IF (ee.LT.0.0d0) THEN
*        print *,'ee is negativ ee=',ee,rn
*        print*,canga,s,rho
*        stop
        ee=0.0d0
        eemm=0.0d0
      ENDIF

      pet = eemm
      pet2 = pet
      petv = pet

*----------------------------------------------------------------------*
* Set switch 'stseas' to find when a new growing season has started.   *
* 'stseas = 1' for the first day of a new growing season.              *
*----------------------------------------------------------------------*

      soil2g = soilw/(dep*10.0d0)

      CALL HYDROLOGY(adp,sfc,sw,sswc,awl,kd,kx,tf,s1,s2,s3,s4,sn,lsn,
     &eemm,etmm,pet2,dp2,s1in,t,rlai,evap,tran,roff,interc,evbs,iter,f2,
     &f3,ft,dep)
      soilw = s1 + s2 + s3 + s4
      soil2g = soilw/(dep*10.0d0)
      
      if ((mnth.eq.4).and.(day.eq.-20)) print*,mnth,day,soilw,etmm

      DO i=1,29
        sm_trig(31-i) = sm_trig(30-i)
      ENDDO
      sm_trig(1) = (s1in - 0.0d0*evbs)

      flow1 = flow1 + f2/10.0d0
      flow2 = flow2 + f3/10.0d0

*----------------------------------------------------------------------*
*                         DAILY npp                                    *
*----------------------------------------------------------------------*
      daynpp = can2a*3600.0d0*hrs/1000000.0d0 -
     &canres*3600.0d0*(24.0d0-hrs)/1000000.0d0

*      print*,'gpp & resp ',12.0*can2a*3600.0d0*hrs/1000000.0d0, 
*     &12.0*canres*3600.0d0*(24.0d0-hrs)/1000000.0d0

      if (day.eq.115)  print*,'canres ',mnth,canres

      daygpp= can2a*3600.0d0*hrs*12.0d0/1000000.0d0
      dayra = canres*3600.0d0*(24.0d0-hrs)/1000000.0d0

      s_gpp = s_gpp + can2a*3600.0d0*hrs*12.0d0/1000000.0d0
      s_lr = s_lr + canres*3600.0d0*(24.0d0-hrs)/1000000.0d0

      IF (veg) THEN
        nppstore = nppstore + daynpp
        IF (ftphen.EQ.1) THEN
          CALL PHENOLOGY1(ftagh,ftdth,bbm,bb0,bbmax,bblim,ssm,sss,sslim,
     &nppstore,leafnpp,stemnpp,rootnpp,tran,rlai,lai,lairat,rem,
     &leafls,stemls,rootls,leafmol,respref,soil2g,wtwp,tmem,leafv,stemv,
     &rootv,daysoff,laimax,lflit,smlit,rtlit,mnth,day,s_ln,s_sr,s_sn,
     &s_rr,s_rn,bb,ss,bbgs,dsbb,nppstorx,nppstor2,daynpp,
     &maxlai,wtfc,yld,resp,sm_trig,suma,tsumam,stemfr,lmor_sc,
     &chill,dschill,canres*3600.0d0*(24.0d0-hrs)/1000000.0d0)
        ELSEIF (ftphen.EQ.2) THEN
          CALL PHENOLOGY2(bbm,bb0,bbmax,bblim,ssm,sss,sslim,
     &nppstore,leafnpp,stemnpp,rootnpp,tran,rlai,lai,lairat,rem,leafls,
     &stemls,rootls,leafmol,respref,soil2g,wtwp,tmem,leafv,stemv,rootv,
     &daysoff,laimax,lflit,smlit,rtlit,mnth,day,s_ln,s_sr,s_sn,s_rr,
     &s_rn,bb,ss,bbgs,dsbb,nppstorx,nppstor2,daynpp,maxlai,
     &wtfc,yld,resp,sm_trig,suma,tsumam,stemfr,lmor_sc,chill,
     &dschill,canres*3600.0d0*(24.0d0-hrs)/1000000.0d0)
        ELSE
          WRITE(*,*) 'No phenology defined for ',ftphen
          STOP
        ENDIF
        leaflit = leaflit + lflit*12.0d0*leafmol
        yield = yield + yld*12.0d0*leafmol
        stemlit = stemlit + smlit*12.0d0
        rootlit = rootlit + rtlit*12.0d0
      ENDIF
*----------------------------------------------------------------------*
* End of the DAILY LOOP.                                               *
*----------------------------------------------------------------------*

      qdirect = qdirect*hrs*3600.0d0
      qdiff = qdiff*hrs*3600.0d0


      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE DOLYMONTH                        *
*                          ********************                        *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE DOLYMONTH(ts,tc,avmnpet,avmnppt,avmnt,h2o,
     &flow1,flow2,c0,n0,minn,nfix,nci,slc,rlc,sln,srespm,
     &lchm,ca,site,year,yr0,yrf,speedc,soilc,soiln,mnth)
*----------------------------------------------------------------------*
      IMPLICIT NONE
c     arguments
      REAL*8 ts,tc,avmnpet,avmnppt,avmnt,h2o,flow1,flow2
      REAL*8 slnr,rlnr,scl,rcl,ca
      REAL*8 c0(8),n0(8),minn(3),nfix,soilc,soiln
      LOGICAL speedc
      INTEGER yr0,yrf,year,site,mnth
      REAL*8 srespm,lchm,nci(4)

c     variables
      REAL*8 ar,fms,fmr,tm,ft,cal,cap,cas,csp,csa,kr(8)
      REAL*8 ml,al,c(8),n(8)
      REAL*8 slc,rlc,sln
      REAL*8 scn,t0,t1
      REAL*8 fl(14)
      INTEGER ij

      if (mnth.eq.-4) stop

*----------------------------------------------------------------------*
* Calculate lignin to nitrogen ratios based on the C02 concentration.  *
*----------------------------------------------------------------------*
      slnr = 0.18d0*ca + 2.7d0
      rlnr = 0.18d0*ca + 2.7d0

      scl = 0.12d0   ! Fractions of surface carbon that are lignin.
      rcl = 0.35d0   ! Fractions of root carbon that are lignin.

      IF (slc*sln.GT.1e-6) THEN
        scn = slc/sln ! Ratio of surface carbon to nitrogen.
      ELSE
        scn = 10.0d0
      ENDIF

*----------------------------------------------------------------------*
* Compute monthly accumulation of values required by century.          *
*----------------------------------------------------------------------*
      CALL SETCENPAR(ts,tc,avmnpet,avmnppt,avmnt,h2o,flow2,slnr,
     &rlnr,ar,fms,fmr,scl,rcl,tm,ft,cal,cap,cas,csp,csa,kr)

*----------------------------------------------------------------------*
* Using century model compute the new carbon and nitrogen pools.       *
*----------------------------------------------------------------------*

      t0 = 0.0d0
      t1 = 1.0d0/12.0d0
      srespm = 0.0d0
      lchm = 0.0d0

      IF (soilc.gt.1000.0d0) THEN

        CALL CDYN(c0,c,t0,t1,slc,rlc,scl,rcl,kr,csa,cas,csp,cap,
     &cal,fms,fmr,ft,srespm,lchm,year,yr0,yrf,speedc)
        CALL FLOWS(c0,c,t0,t1,fl,slc,rlc,scl,rcl,fms,fmr,csa,csp,cas,
     &cap)
        CALL NDYN(n0,n,minn,t0,t1,fl,c,flow1,flow2,ts,
     &ml,al,scn,nfix,cal,nci,site)

        IF (srespm.lt.0.0) THEN
          WRITE(*,*) 'Error: soil respiration negative'
          WRITE(*,'(8f8.0)') c0
        ENDIF

      ELSE
        DO ij=1,8
          c(ij) = c0(ij)
          n(ij) = n0(ij)
        ENDDO
        soiln = soiln + minn(3)

        c(1) = c(1) + (1.0d0 - fms)*slc/12.0d0
        c(2) = c(2) + (1.0d0 - fmr)*rlc/12.0d0
        c(5) = c(5) + fms*slc/12.0d0
        c(6) = c(6) + fmr*rlc/12.0d0

        srespm = 0.0d0
        lchm = 0.0d0

      ENDIF

      soilc = 0.0d0
      soiln = 0.0d0
      DO ij=1,8
        soilc = soilc + c(ij)
        soiln = soiln + n(ij)
        c0(ij) = c(ij)
        n0(ij) = n(ij)
      ENDDO
      soiln = soiln + minn(3)

*      CALL SOILCLOSS(c0,soilcl,scl,kr,rcl,ft,cal)


      RETURN
      END
