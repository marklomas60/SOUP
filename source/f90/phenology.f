*----------------------------------------------------------------------*
*                          PHENOLOGY1                                  *
*                          **********                                  *
*                                                                      *
* Used for grasses and crops. Lai being controlled by the available    *
* npp storage.                                                         *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE PHENOLOGY1(ftagh,ftdth,
     &                      bbm,bb0,bbmax,bblim,ssm,sss,sslim,
     &nppstore,leafnpp,stemnpp,rootnpp,tran,rlai,lai,lairat,rem,leafls,
     &stemls,rootls,leafmol,respref,soil2g,half,tmem,leafv,stemv,rootv,
     &daysoff,laimax,leaflit,stemlit,rootlit,mnth,day,s_ln,s_sr,s_sn,
     &s_rr,s_rn,bb,ss,bbgs,dsbb,nppstorx,nppstor2,daynpp,
     &maxlai,wtfc,yield,resp,sm_trig,suma,tsumam,stemfr,lmor_sc,
     &chill,dschill,s_lr)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 ftagh,nppstore,rootnpp,tran,rlai,lairat,leafmol,respref
      REAL*8 soil2g,half,rem,leaflit,leafv(3600),stemv(1000),rootv(1000)
      REAL*8 leafnpp,stemnpp,laiinc,tmem(200),resp,nppstor2,s_ln,s_mr
      REAL*8 daysoff,laimax,stemnppr,rootnppr,nppstorx,wtfc,s_lr
      REAL*8 s_sn,s_rr,s_rn,stemlit,rootlit,daynpp,maxlai,bb0,bbmax,s_sr
      REAL*8 bblim,sslim,bbsum,yield,sm_trig(30),smtrig,xdaynpp
      REAL*8 suma(360),tsuma,tsumam,stemfr,maint,yy,lmor_sc(3600)
      INTEGER lai,leafls,stemls,rootls,mnth,day,ij,i,bb,gs,bbgs,sssum
      INTEGER bb2bbmin,bb2bbmax,bbm,ssm,sss,ss,dsbb
      INTEGER ftdth,harvest,chill,dschill
      INCLUDE 'param.inc'
      real sumrr,sumsr,sumlr,summr
      save sumrr,sumsr,sumlr,summr


      leaflit = 0.0d0
      stemlit = 0.0d0
      rootlit = 0.0d0
      yield = 0.0d0

      bb2bbmin = 285
      bb2bbmax = 390
      gs = 60

      rem = rlai - int(rlai)
      lai = int(rlai) + 1
      resp = 0.0d0

      harvest = 0

      ss = 0

      xdaynpp = daynpp

*----------------------------------------------------------------------*
* Compute active root and stem npp
*----------------------------------------------------------------------*
      stemnppr = 0.0d0
      rootnppr = 0.0d0
      DO ij=1,stemls
        stemnppr = stemnppr + stemv(ij)
      ENDDO
      DO ij=1,rootls
        rootnppr = rootnppr + rootv(ij)
      ENDDO

*----------------------------------------------------------------------*
* Pay for days roots if veg exists
*----------------------------------------------------------------------*
*      IF (tran/13.4/12.0d0*p_rootsc.LT.0.001d0*nppstore) THEN
*        yy = tran/13.4/12.0d0*p_rootsc
*      ELSE
*        yy = 0.001d0*nppstore
*      ENDIF
*      IF (respref.LT.0.00001d0) THEN
*        yy = 0.0d0
*      ENDIF

      IF ((nppstore.GT.0.0d0).AND.(daynpp.GT.0.0d0)) THEN
        yy = nppstore*p_rootfr
      ELSE
        yy = 0.0d0
      ENDIF

      nppstore = nppstore - yy
      xdaynpp = xdaynpp - yy

      rootnpp = rootnpp + yy
      s_rn = s_rn + yy

*----------------------------------------------------------------------*
* Age roots and set daily vector
*----------------------------------------------------------------------*
      DO ij=1,rootls-1
        rootv(rootls+1-ij) = rootv(rootls-ij)
      ENDDO
      rootv(1) = yy
      IF (rootv(1).LT.0.0d0) THEN
        print*,'rootv negativ: rootv=',rootv(1),'tran=',tran
      ENDIF

*----------------------------------------------------------------------*
* Root resperation
*----------------------------------------------------------------------*
      DO ij=1,rootls-1
        resp = resp + respref*rootv(rootls-ij)
        s_rr = s_rr + respref*rootv(rootls-ij)
        rootnpp = rootnpp - respref*rootv(rootls-ij)
        rootv(rootls-ij) = rootv(rootls-ij)*(1.0 - respref)
      ENDDO

*----------------------------------------------------------------------*
* Check for chilling.                                                  *
*----------------------------------------------------------------------*
      IF (chill.EQ.0) THEN
        bbsum = 0.0d0
        DO i=1,20
          IF (tmem(i).LT.-5.0d0)  bbsum = bbsum + max(-10.0d0,tmem(i)+
     &5.0d0)
        ENDDO
        IF (bbsum.LT.-100) THEN
          chill = 1
          dschill = 1
        ENDIF
      ENDIF
      IF (chill.EQ.1) THEN
        dschill = dschill + 1
      ENDIF
      IF (dschill.GT.260) THEN
        chill = 0
        dschill = 0
      ENDIF

*      print*,mnth,day,bbsum,chill,dschill,' chill'

*----------------------------------------------------------------------*
* Bubburst, if no budburst set, and sufficient soil moisture, then     *
* check for  budburst.                                                 *
*----------------------------------------------------------------------*
      IF (((bb.eq.0).AND.(soil2g.GT.half+0.25d0*(wtfc-half))).OR.
     &((dsbb.GT.bb2bbmax).AND.(soil2g.GT.half+0.1d0*(wtfc-half)))) THEN

        smtrig = 0.0d0
        DO i=1,30
          smtrig = smtrig + sm_trig(i)
        ENDDO

*        print*,mnth,day,smtrig
        IF ((smtrig.GT.30.0d0).OR.(dsbb.GT.bb2bbmax)) THEN

*----------------------------------------------------------------------*
* Check for budburst using degree days.                                *
*----------------------------------------------------------------------*
          bbsum = 0.0d0
          DO i=1,bbm
            IF (tmem(i).GT.bb0)  bbsum = bbsum + min(bbmax,tmem(i)-bb0)
          ENDDO
*          print*,mnth,day,bbsum,bblim,dschill

          IF ((real(bbsum).GE.real(bblim)*exp(-0.01d0*real(dschill)))
     &.OR.(dsbb.GT.bb2bbmax)) THEN
*----------------------------------------------------------------------*
* Adjust proportion of gpp going into stem production based on suma.   *
* This is essentially the LAI control.                                 *
*----------------------------------------------------------------------*
            tsuma = 0.0d0
            DO i=1,360
              tsuma = tsuma + suma(i)
            ENDDO
*            print*,'tsuma ',tsuma,stemfr,nppstore
            maint = max(1.0d0,(real(leafls)/360.0d0))
            tsuma = tsuma - leafmol*1.25d0/maint*p_opt
            IF (tsuma.GT.leafmol*1.25d0/maint*p_opt)
     &tsuma = leafmol*1.25d0/maint*p_opt
            IF (tsuma.LT.-leafmol*1.25d0/maint*p_opt)
     &tsuma = -leafmol*1.25d0/maint*p_opt
            tsumam = 0.0d0*tsumam + 1.0d0*tsuma
            stemfr = stemfr + tsumam*p_laimem
            IF (stemfr.LT.10.0d0) stemfr = 10.0d0

*----------------------------------------------------------------------*
* Bud burst occurance
*----------------------------------------------------------------------*
            bb = (mnth-1)*30 + day
            bbgs = 0
            dsbb = 0
            IF (stemfr.LT.0.8d0*nppstore) THEN
              nppstor2 = nppstore - stemfr
            ELSE 
              IF (stemfr.LT.0.75d0*nppstore) THEN
                nppstor2 = nppstore - stemfr
              ELSE
                nppstor2 = nppstore*0.25d0
              ENDIF
              stemfr = stemfr*0.8d0
            ENDIF
            laiinc = (nppstore - 0.0*nppstor2)/leafmol/1.25d0
            nppstorx = nppstore
*            print*,'laiinc ',laiinc,nppstore*12.0
          ENDIF
        ENDIF
      ENDIF

*----------------------------------------------------------------------*
* Compute length of growing season, and set to zero when equal to      *
* growing season.                                                      *
*----------------------------------------------------------------------*
      IF (bb.GT.0)  bbgs = bbgs + 1
      IF (bbgs-gs.gt.bb2bbmin) THEN
        bb = 0
        bbgs = 0
      ENDIF

      IF (dsbb.LT.500) dsbb = dsbb + 1

*----------------------------------------------------------------------*
* Set LAI increase.                                                    *
*----------------------------------------------------------------------*
      IF ((bb.GT.0).AND.(bbgs.LT.gs).AND.(nppstore.GT.0.0d0)) THEN
        laiinc = lairat*(nppstorx - 0.0d0*nppstor2)/leafmol/1.0d0
*        IF (rlai+laiinc.LT.0.5d0)  laiinc = 0.5d0 - rlai
        IF (rlai+laiinc.GT.maxlai)  laiinc = maxlai - rlai
        IF (rlai+laiinc.GT.11.5d0)  laiinc = 11.5d0 - rlai
        IF ((rlai.GT.0).and.(nppstore.LT.0.0d0)) laiinc = 0.0d0
      ELSE
        laiinc = 0.0d0
      ENDIF

*----------------------------------------------------------------------*
* Senescence, if rlai is greater than zero, compute senescence.        *
*----------------------------------------------------------------------*
      IF (rlai.GT.1.0e-6) THEN
        IF (abs(ftagh).LT.1.0E-6) THEN
*----------------------------------------------------------------------*
* Grass senescence
*----------------------------------------------------------------------*
          IF (soil2g.LT.half*0.5d0) THEN
            laiinc = -rlai
            ss = day + (mnth - 1)*30
          ELSEIF (bbgs.GT.100) THEN
*----------------------------------------------------------------------*
* Check for senescence, senescence occurs there are 'sss' days colder  *
* than 'sslim' out of the last 'ssm' days.                             *
*----------------------------------------------------------------------*
            sssum = 0
            DO i=1,ssm
              IF (tmem(i).LT.sslim)  sssum = sssum + 1
            ENDDO
            IF (sssum.GE.sss) THEN
*----------------------------------------------------------------------*
* Senescence event due to temperature.                                 *
*----------------------------------------------------------------------*
              laiinc =-rlai
              ss = day + (mnth - 1)*30
            ENDIF
          ENDIF

        ELSE
*----------------------------------------------------------------------*
* Crop senescence
*----------------------------------------------------------------------*
          IF ((bbgs.GT.leafls).OR.(bbgs.GT.ftdth)) THEN
            harvest = 1
            laiinc =-rlai
            ss = day + (mnth - 1)*30
          ENDIF

        ENDIF

      ENDIF

*----------------------------------------------------------------------*
* Pay for new leaves.
*----------------------------------------------------------------------*
      IF (laiinc.GT.0.0d0) THEN
        nppstore = nppstore - laiinc*leafmol*1.25d0
        nppstorx = nppstorx - laiinc*leafmol*1.25d0
        leafnpp = leafnpp + laiinc*leafmol
        resp = resp + 0.25d0*laiinc*leafmol
        xdaynpp = xdaynpp - laiinc*leafmol*1.25d0

        s_ln = s_ln + laiinc*leafmol*1.25d0
*        print*,rlai,laiinc,s_ln*12.0
      ENDIF
      s_mr = 0.25d0*laiinc*leafmol

*----------------------------------------------------------------------*
* Age leaves by one day, kill any which have died of old age,
* adjust by laiinc (+ or -), and then sum to get 'rlai'.
* 'leafls' is an integer variable of leaf lifespan in days.
*----------------------------------------------------------------------*
      CALL LAIALT(leafv,rlai,laiinc,leafls,leaflit)

      IF (harvest.EQ.1) THEN
        yield = leaflit*ftagh
        leaflit = (1.0d0 - ftagh)*leaflit
      ELSE
        yield = 0.0d0
      ENDIF

      rem = rlai - int(rlai)
      lai = int(rlai) + 1

*----------------------------------------------------------------------*
* leaf death not through age mortality 
*----------------------------------------------------------------------*
*      DO i=1,leafls
*        rlai = rlai - leafv(i)*(1.0d0 - lmor**(1.0d0/dble(leafls)))
*        leaflit = leaflit + leafv(i)*(1.0d0 - 
*     &lmor**(1.0d0/dble(leafls)))
*        leafv(i)=leafv(i)*lmor**(1.0d0/real(leafls))
*      ENDDO
*----------------------------------------------------------------------*
* leaf death not through age mortality 
*----------------------------------------------------------------------*
*      DO i=1,leafls
*        xx = (real(leafls - i)/real(leafls))**lmor
*        rlai = rlai - leafv(i)*(1.0d0 - xx)
*        leaflit = leaflit + leafv(i)*(1.0d0 - xx)
*        leafv(i)=leafv(i)*xx
*      ENDDO
*----------------------------------------------------------------------*
* leaf death not through age mortality 
*----------------------------------------------------------------------*
      rlai = rlai + leaflit
      DO i=1,leafls
        IF (leafv(i).GT.0.0d0) THEN
          IF (leafv(i).GT.1.0e-6) THEN
            leaflit = leaflit + leafv(i)
            leafv(i)=leafv(i)*lmor_sc(i)
            leaflit = leaflit - leafv(i)
          ELSE
            leaflit = leaflit + leafv(i)
            leafv(i) = 0.0d0
          ENDIF
        ENDIF
      ENDDO
      rlai = rlai - leaflit

*----------------------------------------------------------------------*
* Stem resperation, and NPP                                            *
*----------------------------------------------------------------------*
*      nppstore = nppstore - stemnppr*respref
*      resp = resp + stemnppr*respref
*      s_sr = s_sr + stemnppr*respref
*      xdaynpp = xdaynpp - stemnppr*respref

      s_sr = resp
      DO ij=1,stemls-1
        resp = resp + respref*stemv(stemls-ij)
        stemnpp = stemnpp - respref*stemv(stemls-ij)
        stemv(stemls-ij) = stemv(stemls-ij)*(1.0 - respref)
      ENDDO
      s_sr = resp - s_sr

      DO ij=1,stemls-1
        stemv(stemls+1-ij) = stemv(stemls-ij)
      ENDDO

      IF ((nppstore.GT.0.0d0).and.(daynpp.GT.0.0d0)) THEN
        yy = nppstore*p_stemfr
        stemv(1) = yy
        s_sn = s_sn + yy
        stemnpp = stemnpp + yy
        nppstore = nppstore - yy
      ELSE
        stemv(1) = 0.0d0
      ENDIF

      IF (stemnppr.GT.0d0) THEN
        IF ((stemv(1).LT.p_stmin)) then
          stemnpp = stemnpp + p_stmin - stemv(1)
          nppstore = nppstore - p_stmin + stemv(1)
          s_sn = s_sn + p_stmin - stemv(1)
          stemv(1) = p_stmin
        ENDIF
      ENDIF
*----------------------------------------------------------------------*

      IF (rlai.GT.laimax)  laimax = rlai
      IF (.NOT.(rlai.GT.0.0d0))  daysoff = daysoff + 1.0d0

      IF (nppstore.LT.0.0d0) THEN
        DO i=1,leafls
          leaflit = leaflit + leafv(i)
          leafv(i) = 0.0d0
        ENDDO 
c added by Ghislain 07/10/03
        stemlit = 0.0d0
        DO i=1,stemls
          stemlit = stemlit + stemv(i)
          stemv(i) = 0.0d0
        ENDDO 
c added by Ghislain 07/10/03
        rootlit = 0.0d0
        DO i=1,rootls
          rootlit = rootlit + rootv(i)
          rootv(i) = 0.0d0
        ENDDO 

        nppstore = 0.0d0
        rlai = 0.0d0

      ELSE
        stemlit = 0.0d0
        rootlit = 0.0d0
      ENDIF

      if ((day == 1).and.(mnth == 1)) then
        sumsr=0.0
        sumrr=0.0
        sumlr=0.0
        summr=0.0
      endif
      sumsr = sumsr + 12.0*s_sr
      sumrr = sumrr + 12.0*s_rr
      sumlr = sumlr + 12.0*s_lr
      summr = summr + 12.0*s_mr
      if ((day == -30).and.(mnth == 12)) 
     &print'(4f10.3)',sumsr,sumrr,sumlr,summr


      RETURN
      END

*----------------------------------------------------------------------*
*                          PHENOLOGY2                                  *
*                          **********                                  *
*                                                                      *
* Used for trees. Lai being controlled by a weighted memory of npp     *
* storage over the previous 5 years.                                   *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE PHENOLOGY2(bbm,bb0,bbmax,bblim,ssm,sss,sslim,
     &nppstore,leafnpp,stemnpp,rootnpp,tran,rlai,lai,lairat,rem,leafls,
     &stemls,rootls,leafmol,respref,soil2g,half,tmem,leafv,stemv,rootv,
     &daysoff,laimax,leaflit,stemlit,rootlit,mnth,day,s_ln,s_sr,s_sn,
     &s_rr,s_rn,bb,ss,bbgs,dsbb,nppstorx,nppstor2,daynpp,
     &maxlai,wtfc,yield,resp,sm_trig,suma,tsumam,stemfr,lmor_sc,
     &chill,dschill,s_lr)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 nppstore,rootnpp,tran,rlai,lairat,leafmol,respref,soil2g
      REAL*8 half,rem,leaflit,leafv(3600),stemv(1000),rootv(1000)
      REAL*8 leafnpp,stemnpp,laiinc,tmem(200),resp,nppstor2,s_ln
      REAL*8 daysoff,laimax,stemnppr,rootnppr,nppstorx,wtfc,s_lr,s_mr
      REAL*8 s_sn,s_rr,s_rn,stemlit,rootlit,daynpp,maxlai,bb0,bbmax,s_sr
      REAL*8 bblim,sslim,bbsum,yield,sm_trig(30),smtrig,xdaynpp,ans
      REAL*8 suma(360),tsuma,tsumam,stemfr,maint,yy,lmor_sc(3600)
      real sumrr,sumsr,sumlr,summr,resp_r,resp_s,resp_m,resp_l
      SAVE sumrr,sumsr,sumlr,summr
      INTEGER lai,leafls,stemls,rootls,mnth,day,ij,i,bb,gs,bbgs,sssum
      INTEGER bb2bbmin,bb2bbmax,bbm,ssm,sss,ss,dsbb
      INTEGER chill,dschill
      INCLUDE 'param.inc'

      bb2bbmin = 315
      bb2bbmax = 375
      gs = 30

      rem = rlai - int(rlai)
      lai = int(rlai) + 1
      resp = 0.0d0

      yield = 0.0d0

      ss = 0

*      print*
*      print*,'day npp ',mnth,day,daynpp*12
*      print*,nppstore*12,nppstorx*12,nppstor2*12
*      print*,'stemfr lai ',stemfr*12.0,rlai

      xdaynpp = daynpp

*----------------------------------------------------------------------*
* Compute active root and stem npp
*----------------------------------------------------------------------*
      stemnppr = 0.0d0
      rootnppr = 0.0d0
      DO ij=1,stemls
        stemnppr = stemnppr + stemv(ij)
      ENDDO
      DO ij=1,rootls
        rootnppr = rootnppr + rootv(ij)
      ENDDO

*----------------------------------------------------------------------*
* Pay for days roots if veg exists
*----------------------------------------------------------------------*
*      IF (tran/13.4/12.0d0*p_rootsc.LT.0.01d0*nppstore) THEN
*        yy = tran/13.4/12.0d0*p_rootsc
*      ELSE
*        yy = 0.01d0*nppstore
*      ENDIF
*      IF (respref.LT.0.00001d0) THEN
*        yy = 0.0d0
*      ENDIF

      IF ((nppstore.GT.0.0d0).AND.(daynpp.GT.0.0d0)) THEN
        yy = nppstore*p_rootfr
      ELSE
        yy = 0.0d0
      ENDIF

      nppstore = nppstore - yy
      xdaynpp = xdaynpp - yy

      rootnpp = rootnpp + yy
      s_rn = s_rn + yy

*----------------------------------------------------------------------*
* Age roots and set daily vector
*----------------------------------------------------------------------*
      DO ij=1,rootls-1
        rootv(rootls+1-ij) = rootv(rootls-ij)
      ENDDO
*      rootv(1) = (tran/13.4d0)/12.0d0*p_rootsc
      rootv(1) = yy
      IF (rootv(1).LT.0.0d0) THEN
        print*,'rootv negativ: rootv=',rootv(1),'tran=',tran
      ENDIF

*----------------------------------------------------------------------*
* Root resperation
*----------------------------------------------------------------------*
      resp_r = resp
      DO ij=1,rootls-1
        resp = resp + respref*rootv(rootls-ij)
        s_rr = s_rr + respref*rootv(rootls-ij)
        rootnpp = rootnpp - respref*rootv(rootls-ij)
        rootv(rootls-ij) = rootv(rootls-ij)*(1.0 - respref)
      ENDDO
      resp_r = resp - resp_r

*      print*,'root resp ',resp*12.0

*      nppstore = nppstore - rootnppr*respref
*      resp = resp + rootnppr*respref
*      xdaynpp = xdaynpp - rootnppr*respref
*      s_rr = s_rr + rootnppr*respref

*----------------------------------------------------------------------*
* Check for chilling.                                                  *
*----------------------------------------------------------------------*
      IF (chill.EQ.0) THEN
        bbsum = 0.0d0
        DO i=1,20
          IF (tmem(i).LT.-5.0d0)  bbsum = bbsum + max(-10.0d0,tmem(i)+
     &5.0d0)
        ENDDO
        IF (bbsum.LT.-100) THEN
          chill = 1
          dschill = 1
        ENDIF
      ENDIF
      IF (chill.EQ.1) THEN
        dschill = dschill + 1
      ENDIF
      IF (dschill.GT.260) THEN
        chill = 0
        dschill = 0
      ENDIF

*----------------------------------------------------------------------*
* Bubburst, if no budburst set, and sufficient soil moisture, then     *
* check for  budburst.                                                 *
*----------------------------------------------------------------------*
*      IF ((bb.eq.0).AND.(soil2g.GT.half+0.0*(wtfc-half)).AND.
*     &(day+mnth*30.eq.60)) THEN
*      print'(5i4,3f8.2)',mnth,day,bb,dsbb,bb2bbmax,smtrig,soil2g,half
      IF (((bb.eq.0).AND.(soil2g.GT.half+0.5d0*(wtfc-half))).OR.
     &((dsbb.GT.bb2bbmax).AND.(soil2g.GT.half+0.1d0*(wtfc-half)))) THEN
*      IF ((bb.eq.0).AND.(soil2g.GT.half+xxx*(wtfc-half)).AND.
*     &(dsw.GT.dswmin).AND.(dsw.LT.dswmax)) THEN

        smtrig = 0.0d0
        DO i=1,30
          smtrig = smtrig + sm_trig(i)
        ENDDO

        IF ((smtrig.GT.30.0d0).OR.(dsbb.GT.bb2bbmax)) THEN

*----------------------------------------------------------------------*
* Check for budburst using degree days.                                *
*----------------------------------------------------------------------*
          bbsum = 0.0d0
          DO i=1,bbm
            IF (tmem(i).GT.bb0)  bbsum = bbsum + min(bbmax,tmem(i)-bb0)
          ENDDO

          IF ((real(bbsum).GE.real(bblim)*exp(-0.01d0*real(dschill)))
     &.OR.(dsbb.GT.bb2bbmax)) THEN
*----------------------------------------------------------------------*
* Adjust proportion of gpp going into stem production based on suma.   *
* This is essentially the LAI control.                                 *
*----------------------------------------------------------------------*
            tsuma = 0.0d0
            DO i=1,360
              tsuma = tsuma + suma(i)
            ENDDO
!            print*,'tsuma ',tsuma
            maint = max(1.0d0,(real(leafls)/360.0d0)*1.0d0)
            tsuma = tsuma - leafmol*1.25d0/maint*p_opt
            IF (tsuma.GT.leafmol*1.25d0/maint*p_opt)
     &tsuma = leafmol*1.25d0/maint*p_opt
            IF (tsuma.LT.-leafmol*1.25d0/maint*p_opt)
     &tsuma = -leafmol*1.25d0/maint*p_opt
            tsumam = 0.0d0*tsumam + 1.0d0*tsuma
            stemfr = stemfr + tsumam*p_laimem

            IF (stemfr.LT.10.0d0) stemfr = 10.0d0

*----------------------------------------------------------------------*
* Bud burst occurance.                                                 *
*----------------------------------------------------------------------*
            bb = (mnth-1)*30 + day
            bbgs = 0
            dsbb = 0
            IF (stemfr.LT.0.75d0*nppstore) THEN
              nppstor2 = nppstore - stemfr
            ELSE
              IF (stemfr.LT.0.75d0*nppstore) THEN
                nppstor2 = nppstore - stemfr
              ELSE
                nppstor2 = nppstore*0.25d0
              ENDIF
              stemfr = stemfr*0.95
            ENDIF
            laiinc = (nppstore - nppstor2)/leafmol/1.25d0
            nppstorx = nppstore
!            print*,'laiinc ',laiinc,tsuma,stemfr*12.0
          ENDIF
        ENDIF
      ENDIF

*----------------------------------------------------------------------*
* Compute length of growing season, and set to zero when equal to      *
* growing season.                                                      *
*----------------------------------------------------------------------*
      IF (bb.GT.0)  bbgs = bbgs + 1
      IF (bbgs-gs.gt.bb2bbmin) THEN
        bb = 0
        bbgs = 0
      ENDIF

      IF (dsbb.LT.500) dsbb = dsbb + 1

*----------------------------------------------------------------------*
* Set LAI increase.                                                    *
*----------------------------------------------------------------------*
      IF ((bb.GT.0).AND.(bbgs.LT.gs).AND.(nppstore.GT.1.0d0)) THEN
        laiinc = lairat*(nppstorx - nppstor2)/leafmol/1.25d0
*        IF (rlai+laiinc.LT.0.5d0)  laiinc = 0.5d0 - rlai
        IF (rlai+laiinc.GT.maxlai)  laiinc = maxlai - rlai
        IF (rlai+laiinc.GT.11.5d0)  laiinc = 11.5d0 - rlai
        IF ((rlai.GT.0).and.(nppstore.LT.5.0)) laiinc = 0.0d0
      ELSE
        laiinc = 0.0d0
      ENDIF

*----------------------------------------------------------------------*
* Senescence, if rlai is greater than zero, compute senescence.        *
*----------------------------------------------------------------------*
      IF (rlai.GT.1.0e-6) THEN
        IF (soil2g.LT.half*0.0d0) THEN
*----------------------------------------------------------------------*
* Senescence event due to soil moisture.                               *
*----------------------------------------------------------------------*
          laiinc = -rlai
          ss = day + (mnth - 1)*30
*          print*,'senescence soil'
        ELSEIF (bbgs.GT.100) THEN
*----------------------------------------------------------------------*
* Check for senescence, senescence occurs there are 'sss' days colder  *
* than 'sslim' out of the last 'ssm' days.                             *
*----------------------------------------------------------------------*
          sssum = 0
          DO i=1,ssm
            IF (tmem(i).LT.sslim)  sssum = sssum + 1
          ENDDO
          IF (sssum.GE.sss) THEN
*----------------------------------------------------------------------*
* Senescence event due to temperature.                                 *
*----------------------------------------------------------------------*
            laiinc =-rlai
            ss = day + (mnth - 1)*30
*          print*,'senescence temp'
*          print'(100f6.1)',(tmem(i),i=1,ssm)
          ENDIF
        ENDIF
      ENDIF

*----------------------------------------------------------------------*
* Pay for new leaves.
*----------------------------------------------------------------------*
      resp_m = resp - resp_m
      IF (laiinc.GT.0.0d0) THEN
        nppstore = nppstore - laiinc*leafmol*1.25d0
        nppstorx = nppstorx - laiinc*leafmol*1.25d0
        leafnpp = leafnpp + laiinc*leafmol
        resp = resp + 0.25d0*laiinc*leafmol
        xdaynpp = xdaynpp - laiinc*leafmol*1.25d0
 
        s_ln = s_ln + laiinc*leafmol*1.25d0
      ENDIF
*      print*,'maint & lai inc ',12.0*0.25d0*laiinc*leafmol,laiinc
      s_mr = 0.25d0*laiinc*leafmol
      resp_m = resp - resp_m

*      print*,'plus leaf resp ',resp*12.0

*----------------------------------------------------------------------*
* Age leaves by one day, kill any which have died of old age,
* adjust by laiinc (+ or -), and then sum to get 'rlai'.
* 'leafls' is an integer variable of leaf lifespan in days.
*----------------------------------------------------------------------*
      CALL LAIALT(leafv,rlai,laiinc,leafls,leaflit)

      rem = rlai - int(rlai)
      lai = int(rlai) + 1

*----------------------------------------------------------------------*
* leaf death not through age mortality 
*----------------------------------------------------------------------*
      rlai = rlai + leaflit
      DO i=1,leafls
        IF (leafv(i).GT.0.0d0) THEN
          IF (leafv(i).GT.1.0e-6) THEN
            leaflit = leaflit + leafv(i)
            leafv(i)=leafv(i)*lmor_sc(i)
            leaflit = leaflit - leafv(i)
          ELSE
            leaflit = leaflit + leafv(i)
            leafv(i) = 0.0d0
          ENDIF
        ENDIF
      ENDDO
      rlai = rlai - leaflit

*----------------------------------------------------------------------*
* Stem resperation, and NPP                                            *
*----------------------------------------------------------------------*
*      nppstore = nppstore - stemnppr*respref
*      resp = resp + stemnppr*respref
*      s_sr = s_sr + stemnppr*respref
*      xdaynpp = xdaynpp - stemnppr*respref

      resp_s = resp - resp_s
      s_sr = resp
      DO ij=1,stemls-1
        resp = resp + respref*stemv(stemls-ij)
        stemnpp = stemnpp - respref*stemv(stemls-ij)
        stemv(stemls-ij) = stemv(stemls-ij)*(1.0 - respref)
      ENDDO
      s_sr = resp - s_sr
      resp_s = resp - resp_s

      DO ij=1,stemls-1
        stemv(stemls+1-ij) = stemv(stemls-ij)
      ENDDO

      IF ((nppstore.GT.0.0d0).and.(daynpp.GT.0.0d0)) THEN
        yy = nppstore*p_stemfr
        stemv(1) = yy
        s_sn = s_sn + yy
        stemnpp = stemnpp + yy
        nppstore = nppstore - yy
      ELSE
        stemv(1) = 0.0d0
      ENDIF

      IF (stemnppr.GT.0d0) THEN
        IF ((stemv(1).LT.p_stmin)) then
          stemnpp = stemnpp + p_stmin - stemv(1)
          nppstore = nppstore - p_stmin + stemv(1)
          s_sn = s_sn + p_stmin - stemv(1)
          stemv(1) = p_stmin
        ENDIF
      ENDIF

*      print*,mnth,day,stemv(1),stemnppr*respref

*----------------------------------------------------------------------*

      IF (rlai.GT.laimax)  laimax = rlai
      IF (.NOT.(rlai.GT.0.0d0))  daysoff = daysoff + 1.0d0

      IF (nppstore.LT.0.0d0) THEN
        DO i=1,leafls
          leaflit = leaflit + leafv(i)
          leafv(i) = 0.0d0
        ENDDO 
c added by Ghislain 07/10/03
        stemlit = 0.0d0
        DO i=1,stemls
          stemlit = stemlit + stemv(i)
          stemv(i) = 0.0d0
        ENDDO 
c added by Ghislain 07/10/03
        rootlit = 0.0d0
        DO i=1,rootls
          rootlit = rootlit + rootv(i)
          rootv(i) = 0.0d0
        ENDDO 

        nppstore = 0.0d0
        rlai = 0.0d0

      ELSE
        stemlit = 0.0d0
        rootlit = 0.0d0
      ENDIF

      ans = 0.0d0
      DO i=1,stemls
        ans = ans + stemv(i)
      ENDDO 

!      print*,'plus stem resp ',resp*12.0
!      print*,nppstore*12,nppstorx*12,nppstor2*12
!      print*,'stemfr lai ',12.0*stemfr,rlai

!      print*,mnth,day,s_sr*12.0,s_rr*12.0

      if ((day == 1).and.(mnth == 1)) then
        sumsr=0.0
        sumrr=0.0
        sumlr=0.0
        summr=0.0
      endif
      sumsr = sumsr + 12.0*s_sr
      sumrr = sumrr + 12.0*s_rr
      sumlr = sumlr + 12.0*s_lr
      summr = summr + 12.0*s_mr
      if ((day == -30).and.(mnth == 12)) 
     &print'(4f10.3)',sumsr,sumrr,sumlr,summr

      RETURN
      END

*----------------------------------------------------------------------*
*                          LAIALT                                      *
*----------------------------------------------------------------------*
      SUBROUTINE LAIALT(leafv,rlai,laiinc,leafls,leaflit)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 leafv(3600),rlai,laiinc,sum,leaflit
      INTEGER leafls,i

      leaflit = leafv(leafls)
      DO i=1,leafls-1
        leafv(leafls+1-i) = leafv(leafls-i)
      ENDDO
      leafv(1) = 0.0d0

      IF (laiinc.GT.0.0d0) THEN
        leafv(1) = laiinc
      ELSEIF (laiinc+leaflit.LT.0.0d0) THEN
        sum = laiinc+leaflit
        i = leafls
20      CONTINUE
          sum = sum + leafv(i)
          IF (sum.GT.1.0e-6) then
            leaflit = leaflit + leafv(i) - sum
            leafv(i) = sum
          ELSE
            leaflit = leaflit + leafv(i)
            leafv(i) = 0.0d0
            i = i - 1
          ENDIF
        IF ((.NOT.(sum.GT.0.0d0)).AND.(i.GT.0))  GOTO 20
      ENDIF

      rlai = 0.0d0
      DO i=1,leafls
        rlai = rlai +  leafv(i)
      ENDDO


      RETURN
      END

