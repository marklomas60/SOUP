*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE HYDROLOGY                        *
*                          ********************                        *
*                                                                      *
* Inclusion of century water model, using doly evaporation and         *
* interception calculations.                                           *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE HYDROLOGY(adp,sfc,sw,sswc,awl,kd,kx,tf,s1,s2,s3,s4,sn,
     &lsn,eemm,etmm,pet2,dp2,s1in,t,rlai,evap,tran,roff,interc,evbs,
     &iter,f2,f3,ft,dep)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 adp(4),sfc(4),sw(4),sswc(4),awl(4),kd,kx,tf(13),s1,s2,s3,s4
      REAL*8 ladp(4),lsfc(4),lsw(4),lsswc(4),s1in
      REAL*8 sn,lsn,eemm,etmm,pet2,dp2,t,rlai,evap,tran,roff,f2,f3
      REAL*8 bst,kf,interc,ms,fs,rwc(4),w(4),rem,f1,ds,st,ws,sl,sf,sd
      REAL*8 fav,ev,ss,sums,evbs,pet3,dep,ans1
      INTEGER lai,iter,ft,i
      INCLUDE 'param.inc'

      ds = 0.0d0
      st = 0.0d0

      kf = 1.0d0 - kd
      rem = rlai - int(rlai)
      lai = int(rlai) + 1

      pet3 = pet2
      pet2 = 3.0d0*pet2

*----------------------------------------------------------------------*
* Adjustment for the 'CITY' functional type.                           *
*----------------------------------------------------------------------*
      IF (ft.EQ.2) THEN
        DO i=1,4
          ladp(i) = adp(i)*p_city_dep/dep/10.0d0
          lsfc(i) = sfc(i)*p_city_dep/dep/10.0d0
          lsw(i) = sw(i)*p_city_dep/dep/10.0d0
          lsswc(i) = sswc(i)*p_city_dep/dep/10.0d0
        ENDDO
      ELSE
        DO i=1,4
          ladp(i) = adp(i)
          lsfc(i) = sfc(i)
          lsw(i) = sw(i)
          lsswc(i) = sswc(i)
        ENDDO
      ENDIF

      IF (t.LT.0.0d0) THEN
        s1in = 0.0d0
      ELSE
        s1in = dp2
      ENDIF

*----------------------------------------------------------------------*
      IF (t.LT.0.0d0) THEN
*----------------------------------------------------------------------*
* Snow calculation (mm day-1).                                         *
*----------------------------------------------------------------------*
        sn = sn + dp2
        evap = 0.0d0
        interc = 0.0d0
      ELSE
*----------------------------------------------------------------------*
* Interception water loss (evap mm day-1).                             *
*----------------------------------------------------------------------*
        IF (rlai.GT.0) THEN
          interc = dp2*(1.0 - (tf(lai) + rem*(tf(lai+1) - tf(lai))))
          interc = interc*min(1.0d0,(0.5d0 + t/32.0d0))
          evap = eemm
          IF (evap.GT.pet2) evap = pet2
          IF (evap.GT.interc)  evap = interc
          pet2 = pet2 - evap
          dp2 = dp2 - evap
          interc = evap
        ELSE
          evap = 0.0d0
          interc = 0.0d0
        ENDIF

        IF (sn.GT.0.0) THEN
          lsn = lsn + dp2
        ELSE
*----------------------------------------------------------------------*
* Soil water input 0-150mm (s1).                                       *
*----------------------------------------------------------------------*
          s1 = s1 + dp2
        ENDIF
      ENDIF

      IF (sn.GT.0.0d0) THEN
*----------------------------------------------------------------------*
* Snow sublimination.                                                  *
*----------------------------------------------------------------------*
        sl = 0.85d0*pet2/50.0d0
        IF (sl.GT.sn) THEN
          sl = sn
          pet2 = pet2 - sl
          sn = 0.0d0
        ELSE
          sn = sn - sl
          pet2 = pet2 - sl
*----------------------------------------------------------------------*
* Liquid snow input (lsn - mm day-1).                                  *
*----------------------------------------------------------------------*
          ms = 0.0d0
          IF (t.GT.0.0d0) THEN
            ms = t*40.0d0/real(iter)
            IF (ms.GT.sn) THEN
              ms = sn
            ENDIF
            sn = sn - ms
            lsn = lsn + ms
          ENDIF
        ENDIF
      ELSE
        sl = 0.0d0
      ENDIF

*----------------------------------------------------------------------*
* Soil water input for first layer.                                    *
*----------------------------------------------------------------------*
      fs = 0.0d0
      IF ((lsn.GT.0.05d0*(lsn + sn)).AND.(t.GT.0.0d0)) THEN
        fs = lsn - 0.05d0*(lsn + sn)
        lsn = 0.05d0*(lsn + sn)
        s1 = s1 + fs
        s1in = s1in + fs
      ENDIF

*----------------------------------------------------------------------*
* Soil water 150-300mm (s2).                                           *
*----------------------------------------------------------------------*
      f1 = 0.0d0
      IF (s1.GT.lsfc(1)) THEN
        f1 = s1 - lsfc(1)
        s1 = lsfc(1)
      ENDIF
      s2 = s2 + f1

*----------------------------------------------------------------------*
* Soil water 150-300mm (s2).                                           *
*----------------------------------------------------------------------*
      f2 = 0.0d0
      IF (s2.GT.lsfc(2)) THEN
        f2 = s2 - lsfc(2)
        s2 = lsfc(2)
      ENDIF
      s3 = s3 + f2

*----------------------------------------------------------------------*
* Soil water 300-450mm (s3).                                           *
*----------------------------------------------------------------------*
      f3 = 0.0d0
      IF (s3.GT.lsfc(3)) THEN
        f3 = s3 - lsfc(3)
        s3 = lsfc(3)
      ENDIF
      s4 = s4 + f3

*----------------------------------------------------------------------*
* Fill up to saturated water content from bottom up.                   *
*----------------------------------------------------------------------*
      roff = 0.0d0
      ans1 = lsfc(4) + p_roff2*(lsswc(4) - lsfc(4))
      IF (s4.gt.ans1) THEN
        s3 = s3 + s4 - ans1
        s4 = ans1
        ans1 = lsfc(3) + p_roff2*(lsswc(3) - lsfc(3))
        IF (s3.gt.ans1) THEN
          s2 = s2 + s3 - ans1
          s3 = ans1
          ans1 = lsfc(2) + p_roff2*(lsswc(2) - lsfc(2))
          IF (s2.gt.ans1) THEN
            s1 = s1 + s2 - ans1
            s2 = ans1
            ans1 = lsfc(1) + p_roff2*(lsswc(1) - lsfc(1))
            IF (s1.gt.ans1) THEN
              roff = s1 - ans1
              s1 = ans1
            ENDIF
          ENDIF
        ENDIF
      ENDIF

*----------------------------------------------------------------------*
* Deep storage (ds - mm day-1) and stream (st - mm day-1).             *
*----------------------------------------------------------------------*
      sf = 0.0d0
      sd = 0.0d0
      IF (s4.GT.lsfc(4)) THEN
        sf = kf*(s4 - lsfc(4))*p_roff
        sd = kd*(s4 - lsfc(4))*p_roff
        roff = roff + sf + sd
*      ELSE
*        roff = 0.0d0
      ENDIF
      ss = ds*kx
      ds = ds + sd - ss
      st = st + sf + ss
      s4 = s4 - sf - sd

*----------------------------------------------------------------------*
* Calculation of transpiration (tran - mm day-1).                      *
*----------------------------------------------------------------------*
      fav = 0.0d0

      IF (s1.GT.lsw(1)) THEN
        fav = fav + s1 - lsw(1)
        rwc(1) = (s1 - lsw(1))/(lsfc(1)-lsw(1))
        rwc(1) = 0.0d0
      ELSE
        rwc(1) = 0.0d0
      ENDIF

      IF (s2.GT.lsw(2)) THEN
        fav = fav + s2 - lsw(2)
        rwc(2) = (s2 - lsw(2))/(lsfc(2)-lsw(2))
      ELSE
        rwc(2) = 0.0d0
      ENDIF

      IF (s3.GT.lsw(3)) THEN
        fav = fav + s3 - lsw(3)
        rwc(3) = (s3 - lsw(3))/(lsfc(3)-lsw(3))
      ELSE
        rwc(3) = 0.0d0
      ENDIF

      IF (s4.GT.lsw(4)) THEN
        fav = fav + s4 - lsw(4)
        rwc(4) = (s4 - lsw(4))/(lsfc(4) - lsw(4))
      ELSE
        rwc(4) = 0.0d0
      ENDIF

      w(1) = rwc(1)*awl(1)*ladp(1)
      w(2) = rwc(2)*awl(2)*ladp(2)
      w(3) = rwc(3)*awl(3)*ladp(3)
      w(4) = rwc(4)*awl(4)*ladp(4)
      ws = w(1) + w(2) + w(3) + w(4)

      tran = etmm
      IF (tran.GT.pet2)  tran = pet2
      pet2 = pet2 - tran

      IF (ws.GT.1E-6) THEN
        s1 = s1 - tran*w(1)/ws
        s2 = s2 - tran*w(2)/ws
        s3 = s3 - tran*w(3)/ws
        s4 = s4 - tran*w(4)/ws
      ELSE
        sums = s2 + s3 + s4
        IF (sums.GT.1E-6) THEN
           s1 = s1 - s1*tran/sums
           s2 = s2 - s2*tran/sums
           s3 = s3 - s3*tran/sums
           s4 = s4 - s4*tran/sums
c changed by Ghislain 6/10/03
        ELSE
           s2 = 0.0d0
           s3 = 0.0d0
           s4 = 0.0d0
           tran = 0.0d0
        ENDIF
      ENDIF

      IF (s1.LT.0.0d0) THEN
        s2 = s2 + s1
        s1 = 0.0d0
      ENDIF
      IF (s2.LT.0.0d0) THEN
        s3 = s3 + s2
        s2 = 0.0d0
      ENDIF
      IF (s3.LT.0.0d0) THEN
        s4 = s4 + s3
        s3 = 0.0d0
      ENDIF

      IF (s4.LT.0.0d0) THEN
        s3 = s3 + s4
        s4 = 0.0d0
        IF (s3.LT.0.0d0) THEN
          s2 = s2 + s3
          s3 = 0.0d0
          IF (s2.LT.0.0d0) THEN
            s1 = s1 + s2
            s2 = 0.0d0
            IF (s1.LT.0.0d0) THEN
              tran = tran + s1
              IF (tran.LT.0.0d0)  tran = 0.0d0
              s1 = 0.0d0
              WRITE(*,*) 'No water at all.'
            ENDIF
          ENDIF
        ENDIF        
      ENDIF

*----------------------------------------------------------------------*
* Calculation of bare soil evaporation 'evbs' (mm day-1).              *
*----------------------------------------------------------------------*
*        ev = (rwc(1) - 0.25d0)/0.75d0
*        ev = (rwc(1) - 0.01d0)/0.75d0
      ev = (s2 - lsw(2))/(lsfc(2) - lsw(2))
      IF (ev.LT.0.0d0)  ev = 0.0d0
      bst = 0.0d0
      if (t.gt.0) bst = (t/16.0d0)
      evbs = ev*p_bs*0.33d0*pet3*1.3d0*bst

      IF (evbs.GT.pet2)  evbs = pet2
      pet2 = pet2 - evbs

*      IF (evbs.GT.s1) THEN
*        IF (evbs.GT.(s1+s2))  evbs = s1 + s2
*        s2 = s2 - evbs + s1
*        s1 = 0.0d0
*      ELSE
*        s1 = s1 - evbs
*      ENDIF
      IF (evbs.GT.s1) THEN
        evbs = s1
        s1 = 0.0d0
      ELSE
        s1 = s1 - evbs
      ENDIF

      evap = evap + sl + evbs


      RETURN
      END
