*----------------------------------------------------------------------*
*                          SUBROUTINE NPPCALC                          *
*                          ******************                          *
*----------------------------------------------------------------------*
      SUBROUTINE NPPCALC(npp_eff,c3,maxc,soilc,soiln,minn,soil2g,wtwp,
     &wtfc,rd,rlai,
     &t,rh,ca,oi,rn,qdirect,qdiff,can2a,can2g,canres,suma,amx,amax,
     &gsum,hrs,canga,p,mnth,day,nleaf_sum,fpr,gsm)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 suma,sumd,rlai,soilc,soil2g,wtfc,nmult,npp_eff
      REAL*8 soiln,y1,y0,x1,x0,mmult,minn(3),nup,t
      REAL*8 tk,up,kc,ko,tau,sum,nleaf_sum
      REAL*8 rem,can(12),vm(12),jm(12),jmx(12),oi
      REAL*8 canres,dresp(12),drespt(12),upt(12),nupw
      REAL*8 vmx(12),rh,c1,c2,canga,p,ga,nleaf
      REAL*8 q,qdiff,qdirect,cs,coeff
      REAL*8 can2a,can2g,rn,wtwp,kg,maxc,rht,tpav
      REAL*8 tppcv,tpgsv,ca,rd,tpaj,tppcj,tpgsj,hrs
      REAL*8 a(12),gs(12),ci(12),tpac4,tppcc4,tpgsc4,xvmax
      REAl*8 apar,fpr,gsm

      REAL*8 asunlit,ashade,pcsunlit,pcshade,gssunlit,gsshade
* light-limited assimilation rate (j) and irradiance (q) for sunlit
* and shade
      REAL*8 jsunlit(12),jshade(12),qsunlit(12),qshade(12)
* fraction of sunlit and shade
      REAL*8 fsunlit(12),fshade(12)

      REAL*8 ax,amax,amx,gsum,lyr,qt
      INTEGER i,lai,k,c3,mnth,day
      INCLUDE 'param.inc'

      apar=0.0d0

      rem = rlai - int(rlai)
      lai = int(rlai) + 1
      tk = 273.0d0 + t

      suma = 0.0d0
      sumd = 0.0d0
      nleaf_sum = 0.0d0

      kc = exp(35.8d0 - 80.5d0/(0.00831d0*tk))
      ko = exp(9.6d0 - 14.51d0/(0.00831d0*tk))*1000.0d0
      tau = exp(-3.949d0 + 28.99d0/(0.00831d0*tk))

* changed by Ghislain 08/10/03      IF (rlai.GT.0.1d0) THEN
      q = qdiff + qdirect
      IF ((rlai.GT.0.1d0).AND.(q.GT.0.0d0)) THEN

      IF (soil2g.GT.wtwp) THEN
        kg = maxc*((soil2g - wtwp)/(wtfc - wtwp))**p_kgw
      ELSE
        kg = 0.0d0
      ENDIF
      IF (kg.GT.maxc)  kg = maxc

* Nitrogen uptake CHECK.
*      nup = 120.0d0*(exp(-0.00008d0*soilc))
      nupw = p_nu1*(exp(p_nu2*soilc))*kg**p_nu3
      IF (nupw.LT.0.0d0) nupw = 0.0d0

* Nitrogen multiplier.
      nmult = soiln*p_nu4
      IF (nmult.GE.1.0d0)  nmult = 1.0d0

      y1 = 1.3d0
      y0 = 1.0d0
      x1 = 50.0d0
      x0 = 10.0d0
      mmult = minn(3)*(y1 - y0)/(x1 - x0) + (y0*x1 - x0*y1)/(x1 - x0)
      IF (mmult.GT.y1)  mmult = y1
      IF (mmult.LT.y0)  mmult = y0

      nup = nupw*nmult*mmult
      nup = nupw*nmult
      up = nup

      if (up.LT.0.0d0)  up = 0.0d0

*----------------------------------------------------------------------*
* Total the Beer's Law values for each lai up to the current one.      *
*----------------------------------------------------------------------*
      coeff =-0.5d0
      sum = 0.0d0
      DO i=1,lai-1
        k = i - 1
        sum = sum + exp(coeff*real(k))
      ENDDO
      sum = sum + exp(coeff*real(lai - 1))*rem
*----------------------------------------------------------------------*

      DO i=1,lai
        lyr = real(i)
        IF (i.LT.lai) THEN
          can(i) = exp(coeff*(lyr - 1.0d0))/sum
        ELSE
          can(i) = exp(coeff*(lyr - 1.0d0))/sum
        ENDIF

        upt(i) = up*can(i)

        nleaf = (upt(i))*p_nleaf
        IF (i.LT.lai) THEN
          nleaf_sum = nleaf_sum + nleaf
        ELSE
          nleaf_sum = nleaf_sum + nleaf*rem
        ENDIF

        dresp(i) = nleaf*p_dresp
        drespt(i) = dresp(i)

        IF (kg.GT.1.0e-10) THEN

*        vm(i) = nleaf*p_vm
*        jm(i) = (29.1d0 + 1.64*vm(i))

*        IF (jm(i).lt.0.0d0) jm(i) = 0.0d0
*        vm(i) = vm(i)/1000000.0d0
*        jm(i) = jm(i)/1000000.0d0

*        xx = (-(2.0d0*p_j2-150.0d0*p_j3)-(
*     &(2.0d0*p_j2-150.0d0*p_j3)**2.0d0-12.0d0*p_j3*
*     &(p_j1-50.0d0*p_j2+1875.0d0*p_j3))**0.5d0)/6.0d0/p_j1

*        IF (t.GT.xx) THEN
*          jmx(i) = jm(i)*(1.0d0 + p_j1*(t-25.0d0) +
*     &p_j2*(t - 25.0d0)**2 + p_j3*(t - 25.0d0)**3)
*        ELSE
*          jmx(i) = jm(i)*(1.0d0 + p_j1*(xx-25.0d0) +
*     &p_j2*(xx - 25.0d0)**2 + p_j3*(xx - 25.0d0)**3)
*        ENDIF

*        jmx(i) = jm(i)*(0.4d0 + t*0.6/30.0d0)
*        IF (jmx(i).GT.jm(i))  jmx(i) = jm(i)

*        yy = (-(2.0d0*p_v2-150.0d0*p_v3)-(
*     &(2.0d0*p_v2-150.0d0*p_v3)**2.0d0-12.0d0*p_v3*
*     &(p_v1-50.0d0*p_v2+1875.0d0*p_v3))**0.5d0)/6.0d0/p_v
*        IF (t.GT.yy) THEN
*          vmx(i) = vm(i)*(1.0d0 + p_v1*(t - 25.0d0) +
*     &p_v2*(t - 25.0d0)**2 + p_v3*(t - 25.0d0)**3)
*        ELSE
*          vmx(i) = vm(i)*(1.0d0 + p_v1*(yy - 25.0d0) +
*     &p_v2*(yy - 25.0d0)**2 + p_v3*(yy - 25.0d0)**3)
*        ENDIF

*        vmx(i) = vm(i)*(0.4d0 + t*0.6/30.0d0)
*        IF (vmx(i).GT.vm(i))  vmx(i) = vm(i)

*        xx = (0.25 + (1.0d0-0.25d0)/25.0d0*t)**1.5d0

        vm(i) = nleaf*p_vm
        qt = 2.3d0
        vmx(i)=vm(i)*qt**(t/10.0d0)/qt**(2.5d0)
        if (t.gt.30.0d0) vmx(i)=vm(i)*qt**(30.0/10.0d0)/qt**(2.5d0)
        vmx(i) = vmx(i)*npp_eff*kg**p_nu3

*        vmx(i) = vmx(i)*(1.0d0+0.25d0*(1.0d0-ca/35.0d0))

        jmx(i) = (29.1d0 + 1.64d0*vmx(i))
        vmx(i) = vmx(i)*0.000001d0
        jmx(i) = jmx(i)*0.000001d0

        rd = vmx(i)*0.02d0

* change Ghislain
*        irc = q*exp(-0.5d0*lyr)

*        CALL BEERSLAW(lyr-0.5d0,rlai,q,0.0d0,
*     &       fsunlit(i),qsunlit(i), fshade(i),qshade(i))

        CALL GOUDRIAANSLAW(lyr-0.5d0,rlai,qdirect,qdiff,
     &fsunlit(i),qsunlit(i), fshade(i),qshade(i))

        IF (i.LT.lai) THEN
           apar = apar + fsunlit(i)*qsunlit(i)+fshade(i)*qshade(i)
        ELSE
           apar = apar + rem*
     &(fsunlit(i)*qsunlit(i)+fshade(i)*qshade(i))
        ENDIF

*      if (i.eq.1) then
*        print*,jmx(1),qsunlit(1),qshade(1)
*      endif

*        j(i) = 0.24d0*irc/(1.0d0 + (0.24d0**2)*(irc**2)/(jmx(i)**2))
*     &**0.5d0

*     calculate separatly the sunlit and shaded light limited 
*     assimilation rate
        jsunlit(i)=0.24d0*qsunlit(i)/
     &(1.0d0+(0.24d0**2)*(qsunlit(i)**2)/(jmx(i)**2))**0.5d0
        jshade(i)=0.24d0*qshade(i)/
     &(1.0d0+(0.24d0**2)*(qshade(i)**2)/(jmx(i)**2))**0.5d0

        ELSE

          CALL GOUDRIAANSLAW(lyr-0.5d0,rlai,qdirect,qdiff,
     &fsunlit(i),qsunlit(i), fshade(i),qshade(i))

          IF (i.LT.lai) THEN
             apar = apar + fsunlit(i)*qsunlit(i)+fshade(i)*qshade(i)
          ENDIF

          vmx(i) = 0.0d0
          jsunlit(i) = 0.0d0
          jshade(i) = 0.0d0

        ENDIF

      ENDDO

**********************************************************************
*     fAPAR calculation added here from gh nppcalc.f amf251105
**********************************************************************

      fpr=apar/(qdiff+qdirect)
      
*********************************************************************      

      if (mnth.eq.13) print'(6f8.2)',vm(1),
     &vmx(1)*1.0e6,jmx(1)*1.0e6,kg,up,t

*      if ((mnth.eq.8).and.(day.eq.15)) print*,kg**p_nu3,' nppcalc'

      if (day.eq.35) print'(''vmax jsn jsh '',3f6.2)',1.0e6*vmx(1),
     &1.0e6*jsunlit(1),1.0e6*jshade(1)

* removed by Ghislain 24/11/03
* the output xvmx, xjmx and xj are not used anymore
*      CALL XVJMAX(xvmx,xjmx,xj,t,q,sum,nup,oi,xdresp)

      canres = 0.0d0
      DO i=1,lai
        IF (i.LT.lai) THEN
          canres = canres + dresp(i)
        ELSE
          canres = canres + dresp(i)*rem
        ENDIF
      ENDDO

      ELSE
        DO i=1,lai
          can(i) = 0.0d0
          upt(i) = 0.0d0
          dresp(i) = 0.0d0
          drespt(i) = dresp(i)
          vm(i) = 0.0d0
          jm(i) = 0.0d0
          jmx(i) = 0.0d0
          vmx(i) = 0.0d0

          jshade(i) = 0.0d0
          jsunlit(i) = 0.0d0
          fsunlit(i) = 0.0d0
          fshade(i) = 0.0d0
          qsunlit(i) = 0.0d0
          qshade(i) = 0.0d0
        ENDDO
        canres = 0.0d0
      ENDIF

      c1 = 142.4d0 - 4.8d0*t
      IF (c1.GT.80.0d0)  c1 = 80.0d0
      IF (c1.LT.8.0d0)  c1 = 8.0d0
      c2 = 12.7d0 - 0.207d0*t
      IF (c2.GT.10.0d0)  c2 = 10.0d0
      IF (c2.LT.6.9d0)  c2 = 6.9d0

*----------------------------------------------------------------------*
* Assimilation calculations using subroutines ASSVMAX and ASSJ.        *
*----------------------------------------------------------------------*
      can2a = 0.0d0
      can2g = 0.0d0

      DO i=1,lai

        dresp(i) = drespt(i)

        IF (rlai.GT.0.1d0) THEN
          IF ((t.GT.0.0d0).AND.(rn.GT.0.0d0).AND.(soil2g.GT.wtwp).AND.
     &(fshade(i)+fsunlit(i).GT.0.1d0)) THEN
            ga = canga/(8.3144d0*tk/p)/rlai

*----------------------------------------------------------------------*
* Only recompute assimilation every few days to speed things up, but   *
* make sure it is calculated on the first day of any growing season.   *
*----------------------------------------------------------------------*
            rht = rh

            IF (c3.EQ.1) THEN
              CALL ASSVMAX(tpav,tppcv,cs,tpgsv,oi,ca,vmx(i),
     &rh,kg,ko,kc,tau,rd,ga,t,p)
*     make the sum of assimilation, internal concentration and stomatal
*     conductance.
              a(i)=0.0d0
              ci(i)=0.0d0
              gs(i)=0.0d0

              IF (fshade(i).GT.0.01d0) THEN
                CALL ASSJ(tpaj,tppcj,tpgsj,oi,ca,rh,kg,tau,rd,ga,t,p,
     &jshade(i))
                CALL LIMITATIONPROCESS(tpav,tppcv,tpgsv,tpaj,tppcj,
     &tpgsj,ashade,pcshade,gsshade)
                a(i)=a(i)+fshade(i)*ashade
                ci(i)=ci(i)+fshade(i)*pcshade
                gs(i)=gs(i)+fshade(i)*gsshade
              ENDIF

              IF (fsunlit(i).GT.0.01d0) THEN
                CALL ASSJ(tpaj,tppcj,tpgsj,oi,ca,rh,kg,tau,rd,ga,t,p,
     &jsunlit(i))
                CALL LIMITATIONPROCESS(tpav,tppcv,tpgsv,tpaj,tppcj,
     &tpgsj,asunlit,pcsunlit,gssunlit)
                a(i)=a(i)+fsunlit(i)*asunlit
                ci(i)=ci(i)+fsunlit(i)*pcsunlit
                gs(i)=gs(i)+fsunlit(i)*gssunlit
              ENDIF

            ELSE

              tpac4=0.0d0
              tppcc4=0.0d0
              tpgsc4=0.0d0

              IF (fshade(i).GT.0.01d0) THEN
                CALL ASSC4(ashade,pcshade,gsshade,ca,xvmax,rht,kg,t,
     &qshade(i),p,upt(i))
                tpac4=tpac4+fshade(i)*ashade
                tppcc4=tppcc4+fshade(i)*pcshade
                tpgsc4=tpgsc4+fshade(i)*gsshade
              ENDIF
              IF (fsunlit(i).GT.0.01d0) THEN
                CALL ASSC4(asunlit,pcsunlit,gssunlit,ca,xvmax,rht,kg,t,
     &qsunlit(i),p,upt(i))
                tpac4=tpac4+fsunlit(i)*asunlit
                tppcc4=tppcc4+fsunlit(i)*pcsunlit
                tpgsc4=tpgsc4+fsunlit(i)*gssunlit
              ENDIF

              a(i) = tpac4
              gs(i) = tpgsc4
              ci(i) = tppcc4

              IF (gs(i).LT.0.0d0)  gs(i) = 0.0d0
              gs(i) = gs(i)/1000.0d0

            ENDIF
*----------------------------------------------------------------------*
* End of assimilation calculation 'IF' statement.                      *
*----------------------------------------------------------------------*
          ELSE
            a(i) = 0.0d0
            gs(i) = 0.005d0
            ci(i) = ca
          ENDIF
*----------------------------------------------------------------------*
* End of 't & rn > 0' 'IF' statement.                                  *
*----------------------------------------------------------------------*
        ELSE
          a(i) = 0.0d0
          gs(i) = 0.0d0
          dresp(i) = 0.0d0
          ci(i) = ca
        ENDIF
*----------------------------------------------------------------------*
* End of 'leaf' 'IF' statement.                                        *
*----------------------------------------------------------------------*

*----------------------------------------------------------------------*
* At the start of a season recompute dark resperation.                 *
*----------------------------------------------------------------------*
*        IF (stseas.EQ.1)  dresp(i) = drespt(i)

*----------------------------------------------------------------------*
* Cumulate daily assimilation of the last lai layer, 'suma'.           *
*----------------------------------------------------------------------*
        IF (i.EQ.lai) THEN
          IF (lai.GT.1) THEN
            suma = suma + (rem*a(lai) + (1.0d0 - rem)*a(lai - 1))*
     &3600.0d0*hrs - (rem*dresp(lai) + (1.0d0 - rem)*dresp(lai -
     &1))*3600.0d0*(24.0d0 - hrs)
            sumd = sumd + (rem*dresp(lai) + (1.0d0 - rem)*dresp(lai -
     &1))*3600.0d0*(24.0d0 - hrs)
          ELSE
            suma = suma + rem*a(lai)*3600.0d0*hrs - rem*dresp(lai)*
     &3600.0d0*(24.0d0 - hrs)
            sumd = sumd + rem*dresp(lai)*3600.0d0*(24.0d0 - hrs)
            suma = 0.0d0
          ENDIF
        ENDIF
*----------------------------------------------------------------------*

        IF (i.LT.lai) THEN
          can2a = can2a + a(i)
          can2g = can2g + gs(i)
        ELSE
          can2a = can2a + a(i)*rem
          can2g = can2g + gs(i)*rem
        ENDIF

      ENDDO

      gsm = gs(1)

*----------------------------------------------------------------------*
* End of LAI loop.                                                     *
*----------------------------------------------------------------------*

      if (day.eq.35) print'(6f6.2)',(a(i),i=1,6)
      if ((mnth.eq.8).and.(day.eq.35)) print*,can2a

      IF (lai.GT.1) THEN
        gsum = gs(lai - 1)*(1.0d0 - rem) + gs(lai)*rem
      ELSE
        gsum = gs(1)
      ENDIF

      IF ((rlai.GT.0.0d0).AND.(t.GT.0.0d0).AND.(rn.GT.0.0d0)) THEN
        IF (mod(day-1,1).EQ.0) THEN

*          CALL XAMAX(ax,t,oi,xvmx,rd,xj,xdresp,ci(1))
* Temporary change by Ghislain
           ax=0
* end if change
          IF (rlai.LT.-1.0d0) THEN
            IF (ax*rem.GT.amax)  THEN
              amax = a(1)*rem
            ENDIF
            IF (ax*rem.GT.amx)  amx = ax*rem
          ELSE
            IF (ax.GT.amax)  THEN
              amax = a(1)
            ENDIF
            IF (ax.GT.amx)  amx = ax
          ENDIF

        ENDIF
      ENDIF

      suma = suma/1000000.0d0


      RETURN
      END

*----------------------------------------------------------------------*
*     determin which process light, or Rubisco is the limiting         *
*     factor and return the assimilation, partial presure and          *
*     stomatal conductance corresponding to the limiting process       *
*----------------------------------------------------------------------*


      SUBROUTINE LIMITATIONPROCESS(av,pcv,gsv,aj,pcj,gsj,a,pc,gs)
      REAL*8 av,pcv,gsv,aj,pcj,gsj,a,pc,gs

      IF (av.LT.aj) THEN
         a = av
         gs = gsv
         pc = pcv
      ELSE
         a = aj
         gs = gsj
         pc = pcj
      ENDIF


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE ASSVMAX                          *
*                          ******************                          *
*                                                                      *
* ASSVMAX computes the assimilation rate by forming a cubic equation   *
* from the equations defining assimilation rate, stomatal conductance, *
* CO2 partial pressure (assimilation diffusion equaiton) and the       *
* carboxylation rate. The cubic equation is solved and the largest     *
* real root is assumed to be the assimilation rate. The variables here *
* correspond to the Doly paper, in the main program they are such that *
*        po,pa,vmax,g0,g1,rh,ko,kc,kg,tau,rd                           *
*     = oi,ca,vmx(i),c1,c2,rh,ko,kc,kg,tau,rd.                         *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE ASSVMAX(a,ci,cs,gs,oi,ca,vmx,rh,kg,ko,kc,tau,rd,ga,t,p)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 a,ci,gs,oi,ca,vmx,rh,ko,kc,tau,rd,fa0,faf,gam,step
      REAL*8 p,ga,x,y,dv,t,a0,cs,af,fa,kg,ax,bx,cx,gsmin
      INTEGER i

      gsmin = 0.005d0

      x = 3.0d0 + 0.2d0*t
      if (x.LT.0.25d0) x = 0.25d0
      y = 1.5d0

      dv = 0.6108d0*exp(17.269d0*t/(237.3d0 + t))*(1.0d0 - 
     &rh/100.0d0)

      if (5.gt.30) then
      do i=0,200
      a = real(i)/10.0d0
      cs = ca - a*1.0e-6*p/ga
      gs = gsmin + (x*a/(1.0d0 + dv/y)/(cs*10.0d0 - 
     &1.54d0*t))*kg
      ci = cs - a*1.0e-6*p/gs
      gam = 1.0d0 - 0.5d0*oi/tau/ci
      fa0 = a - (gam*vmx*ci/(ci + kc*(1.0d0 + oi/ko)) -
     &rd)*1.0e6
      print'(7f10.3)',a,fa0,cs,gs,ci,gam,(x*a/(1.0d0 + dv/y)/(cs*10.0d0
     & - 1.54d0*t))
      enddo
      endif

      step = 1.0d0
      a = 0.0d0
      cs = ca - a*1.0e-6*p/ga
      gs = gsmin + (x*a/(1.0d0 + dv/y)/(cs*10.0d0 - 
     &1.54d0*t))*kg
      ci = cs - a*1.0e-6*p/gs
      gam = 1.0d0 - 0.5d0*oi/tau/ci
      fa0 = a - (gam*vmx*ci/(ci + kc*(1.0d0 + oi/ko)) -
     &rd)*1.0e6

      a0 = a

50    CONTINUE
        a = a + step
        cs = ca - a*1.0e-6*p/ga
        gs = gsmin + (x*a/(1.0d0 + dv/y)/(cs*10.0d0 - 
     &1.54d0*t))*kg
        ci = cs - a*1.0e-6*p/gs
        gam = 1.0d0 - 0.5d0*oi/tau/ci
        faf = a - (gam*vmx*ci/(ci + kc*(1.0d0 + oi/ko)) -
     &rd)*1.0e6
      IF (faf.LT.0.0d0) THEN
        fa0 = faf
        a0 = a
        GOTO 50
      ENDIF
      af = a

      IF (af.GT.100d0) THEN
        a = 0.0d0
        cs = ca - a*1.0e-6*p/ga
        gs = gsmin + (x*a/(1.0d0 + dv/y)/(cs*10.0d0 - 
     &1.54d0*t))*kg
        ci = cs - a*1.0e-6*p/gs
        gam = 1.0d0 - 0.5d0*oi/tau/ci
      ELSE
        a = (a0 + af)/2.0d0
        cs = ca - a*1.0e-6*p/ga
        gs = gsmin + (x*a/(1.0d0 + dv/y)/(cs*10.0d0 - 
     &1.54d0*t))*kg
        ci = cs - a*1.0e-6*p/gs
        gam = 1.0d0 - 0.5d0*oi/tau/ci
        fa = a - (gam*vmx*ci/(ci + kc*(1.0d0 + oi/ko)) -
     &rd)*1.0e6

        bx = ((a+a0)*(faf-fa)/(af-a)-(af+a)*(fa-fa0)/(a-a0))/(a0-af)
        ax = (faf-fa-bx*(af-a))/(af**2-a**2)
        cx = fa-bx*a-ax*a**2

        IF (abs(ax).GT.0.0d0) THEN
          IF (bx**2-4.0d0*ax*cx.GT.0.0d0) THEN
            a = (-bx+(bx**2-4.0d0*ax*cx)**0.5)/(2.0d0*ax)
        IF (a.GT.af)  a = (-bx-(bx**2-4.0d0*ax*cx)**0.5)/(2.0d0*ax)
          ELSE
            a = 0.0d0
          ENDIF
        ELSE
          IF (abs(bx).GT.0.0d0) THEN
            a =-cx/bx
          ELSE
            a = 0.0d0
          ENDIF
        ENDIF

      ENDIF

        cs = ca - a*1.0e-6*p/ga
        gs = gsmin + (x*a/(1.0d0 + dv/y)/(cs*10.0d0 - 
     &1.54d0*t))*kg
        ci = cs - a*1.0e-6*p/gs
        gam = 1.0d0 - 0.5d0*oi/tau/ci
*      print'(6f10.3)',a,fa0,cs,gs,ci,gam
*      print*

      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE ASSJ                             *
*                          ***************                             *
*                                                                      *
* ASSVMAX computes the assimilation rate by forming a cubic equation   *
* from the equations defining assimilation rate, stomatal conductance, *
* CO2 partial pressure (assimilation diffusion equaiton) and the       *
* carboxylation rate. The cubic equation is solved and the largest     *
* real root is assumed to be the assimilation rate. The variables here *
* correspond to the Doly paper, in the main program they are such that *
*        po,pa,vmax,g0,g1,rh,ko,kc,kg,tau,rd                           *
*     = oi,ca,vmx(i),c1,c2,rh,ko,kc,kg,tau,rd.                         *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE ASSJ(a,ci,gs,oi,ca,rh,kg,tau,rd,ga,t,p,j)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 a,ci,gs,oi,ca,rh,tau,rd,j,gam,fa0,faf,a0,af,step
      REAL*8 p,ga,x,y,dv,t,cs,fa,kg
      REAL*8 ax,bx,cx,gsmin

      gsmin = 0.005d0

      x = 3.0d0 + 0.2d0*t
      if (x.LT.0.25d0) x = 0.25d0
      y = 1.5

      dv = 0.6108d0*exp(17.269d0*t/(237.3d0 + t))*(1.0d0 - 
     &rh/100.0d0)

      step = 1.0d0
      a = 0.0d0
      cs = ca - a*1.0e-6*p/ga
      gs = gsmin + (x*a/(1.0d0 + dv/y)/(cs*10.0d0 - 
     &1.54d0*t))*kg
      ci = cs - a*1.0e-6*p/gs
      gam = 1.0d0 - 0.5d0*oi/tau/ci
      fa0 = a - (gam*j*ci/4.0d0/(ci + oi/tau) -
     &rd)*1.0e6

      a0 = a

50    CONTINUE
        a = a + step
        cs = ca - a*1.0e-6*p/ga
        gs = gsmin + (x*a/(1.0d0 + dv/y)/(cs*10.0d0 - 
     &1.54d0*t))*kg
        ci = cs - a*1.0e-6*p/gs
        gam = 1.0d0 - 0.5d0*oi/tau/ci
        faf = a - (gam*j*ci/4.0d0/(ci + oi/tau) -
     &rd)*1.0e6
      IF (faf.LT.0.0d0) THEN
        fa0 = faf
        a0 = a
        GOTO 50
      ENDIF
      af = a

      IF (af.GT.100d0) THEN
        a = 0.0d0
        cs = ca - a*1.0e-6*p/ga
        gs = gsmin + (x*a/(1.0d0 + dv/y)/(cs*10.0d0 - 
     &1.54d0*t))*kg
        ci = cs - a*1.0e-6*p/gs
        gam = 1.0d0 - 0.5d0*oi/tau/ci
      ELSE

        a = (a0 + af)/2.0d0
        cs = ca - a*1.0e-6*p/ga
        gs = gsmin + (x*a/(1.0d0 + dv/y)/(cs*10.0d0 - 
     &1.54d0*t))*kg
        ci = cs - a*1.0e-6*p/gs
        gam = 1.0d0 - 0.5d0*oi/tau/ci
        fa = a - (gam*j*ci/4.0d0/(ci + oi/tau) -
     &rd)*1.0e6

        bx = ((a+a0)*(faf-fa)/(af-a)-(af+a)*(fa-fa0)/(a-a0))/(a0-af)
        ax = (faf-fa-bx*(af-a))/(af**2-a**2)
        cx = fa-bx*a-ax*a**2

        IF (abs(ax).GT.0.0d0) THEN
          IF (bx**2-4.0d0*ax*cx.GT.0.0d0) THEN
            a = (-bx+(bx**2-4.0d0*ax*cx)**0.5)/(2.0d0*ax)
        IF (a.GT.af)  a = (-bx-(bx**2-4.0d0*ax*cx)**0.5)/(2.0d0*ax)
          ELSE
            a = 0.0d0
          ENDIF
        ELSE
          IF (abs(bx).GT.0.0d0) THEN
            a =-cx/bx
          ELSE
            a = 0.0d0
          ENDIF
        ENDIF

      ENDIF

      cs = ca - a*1.0e-6*p/ga
      gs = gsmin + (x*a/(1.0d0 + dv/y)/(cs*10.0d0 - 
     &1.54d0*t))*kg
      ci = cs - a*1.0e-6*p/gs
      gam = 1.0d0 - 0.5d0*oi/tau/ci

      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE ASSVMAX                          *
*                          ******************                          *
*                                                                      *
* ASSVMAX computes the assimilation rate by forming a cubic equation   *
* from the equations defining assimilation rate, stomatal conductance, *
* CO2 partial pressure (assimilation diffusion equaiton) and the       *
* carboxylation rate. The cubic equation is solved and the largest     *
* real root is assumed to be the assimilation rate. The variables here *
* correspond to the Doly paper, in the main program they are such that *
*        po,pa,vmax,g0,g1,rh,ko,kc,kg,tau,rd                           *
*     = oi,ca,vmx(i),c1,c2,rh,ko,kc,kg,tau,rd.                         *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE ASSVMAX2(a,w,pc,gs,po,pa,vmax,g0,g1,rh,ko,kc,kg,tau,rd
     &,acheck)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 a,w,pc,gs,po,pa,vmax,g0,g1,rh,ko,kc,kg,tau,rd,acheck,r2q3
      REAL*8 ax,bx,cx,dx,ex,fx,gx,hx,p,q,r,qt,rt,th,a1,a2,a3,pi,at,bt,sc
      REAL*8 arg1

      pi = atan(1.0d0)*4.0d0
      sc = 1000000.0d0
*----------------------------------------------------------------------*
* 'a2,b2,c2,d2' are such that 'w=(a2a + b2)/(c2a + d2)'.               *
*----------------------------------------------------------------------*
      ax = ko*pa*vmax*(kg*g1*rh - 160.0d0)
      bx = kg*g0*pa**2.0d0*vmax*ko
      cx =-160.0d0*pa*ko + kc*kg*ko*g1*rh + pa*kg*ko*g1*rh +
     & kc*kg*po*g1*rh
      dx = pa**2.0d0*kg*ko*g0 + kc*kg*po*g0*pa + kc*kg*ko*g0*pa

*----------------------------------------------------------------------*
* 'ex,fx,gx,hx' are such that '0.5po/(taupc)-1=(exa+fx)/(gxa+hx)'.     *
*----------------------------------------------------------------------*
      ex =-tau*pa*kg*g1*rh + 160.0d0*tau*pa + 0.5d0*po*kg*g1*rh
      fx =-tau*pa**2*kg*g0 + 0.5d0*po*kg*g0*pa
      gx = tau*pa*(kg*g1*rh - 160.0d0)
      hx = tau*pa**2.0d0*kg*g0

*----------------------------------------------------------------------*
* Cubic coefficients 'p,q,r'.                                          *
*----------------------------------------------------------------------*
      p = sc*rd + (sc*ax*ex + cx*hx + dx*gx)/(cx*gx)
      q = (sc*rd*cx*hx + sc*rd*dx*gx + sc*ax*fx + sc*bx*ex +
     & dx*hx)/(cx*gx)
      r = (sc*rd*dx*hx + sc*bx*fx)/(cx*gx)

*----------------------------------------------------------------------*
* Sove the cubic equaiton to find 'a'.                                 *
*----------------------------------------------------------------------*
      Qt = (p**2 - 3.0d0*q)/9.0d0
      Rt = (2.0d0*p**3 - 9.0d0*p*q + 27.0d0*r)/54.0d0

      r2q3 = Rt**2 - Qt**3

      IF (r2q3.LT.0.0d0) THEN
        arg1 = Rt/Qt**1.5d0
        IF (arg1.GT.1.0d0)  arg1 = 1.0d0
        IF (arg1.LT.-1.0d0)  arg1 = -1.0d0
        th = acos(arg1)
        a1 =-2.0d0*Qt**0.5d0*cos(th/3.0d0) - p/3.0d0
        a2 =-2.0d0*Qt**0.5d0*cos((th + 2.0d0*pi)/3.0d0) - p/3.0d0
        a3 =-2.0d0*Qt**0.5d0*cos((th + 4.0d0*pi)/3.0d0) - p/3.0d0
        a = a1
        IF (a2.GT.a)  a = a2
        IF (a3.GT.a)  a = a3
      ELSE
        at =-Rt/abs(Rt)*(abs(Rt) + r2q3**0.5d0)**(1.0d0/3.0d0)
        IF (abs(at).GT.1E-6) THEN
          bt = Qt/at
        ELSE
          bt = 0.0d0
        ENDIF
        a1 = at + bt - p/3.0d0
*        a2 =-0.5d0*(at + bt) - p/3.0d0 + i*3.0d0**0.5d0*(at - bt)/2.0d0
*        a3 =-0.5d0*(at + bt) - p/3.0d0 - i*3.0d0**0.5d0*(at - bt)/2.0d0
        a = a1
      ENDIF

*----------------------------------------------------------------------*
* Compute 'gs,pc,w' corresponding to 'a'.                              *
*----------------------------------------------------------------------*
      gs = (g0 + g1*a*rh/pa)*kg
      pc = pa - a*160.0d0/gs
      w = vmax*pc/(pc + kc*(1 + po/ko))
      acheck = (w*(1.0d0 - 0.5*po/(tau*pc)) - rd)*sc


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE ASSJ                             *
*                          ***************                             *
*                                                                      *
* ASSJ computes the assimilation rate by forming a cubic equation      *
* from the equations defining assimilation rate, stomatal conductance, *
* CO2 partial pressure (assimilation diffusion equaiton) and the       *
* carboxylation rate. The cubic equation is solved and the largest     *
* real root is assumed to be the assimilation rate. The variables here *
* correspond to the Doly paper, in the main program they are such that *
*        po,pa,j,g0,g1,rh,kg,tau,rd                                    *
*     = oi,ca,j(i),c1,c2,rh,kg,tau,rd                                  *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE ASSJ2(a,w,pc,gs,po,pa,j,g0,g1,rh,kg,tau,rd,acheck)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 a,w,pc,gs,po,pa,j,g0,g1,rh,kg,tau,rd,acheck,r2q3,arg1
      REAL*8 ax,bx,cx,dx,ex,fx,gx,hx,p,q,r,qt,rt,th,a1,a2,a3,pi,at,bt,sc

      pi = atan(1.0d0)*4.0d0
      sc = 1000000.0d0

*----------------------------------------------------------------------*
* 'a2,b2,c2,d2' are such that 'w=(a2a + b2)/(c2a + d2)'.               *
*----------------------------------------------------------------------*
      ax = j*tau*pa*(kg*g1*rh - 160.0d0)
      bx = j*tau*pa*kg*g0*pa
      cx = 4.0d0*pa*kg*tau*g1*rh - 640.0d0*pa*tau + 4.0d0*po*kg*g1*rh
      dx = 4.0d0*pa**2.0d0*kg*tau*g0 + 4.0d0*po*kg*g0*pa

*----------------------------------------------------------------------*
* 'ex,fx,gx,hx' are such that '0.5po/(taupc)-1=(exa+fx)/(gxa+hx)'.     *
*----------------------------------------------------------------------*
      ex =-tau*pa*kg*g1*rh + 160.0d0*tau*pa + 0.5d0*po*kg*g1*rh
      fx =-tau*pa**2.0d0*kg*g0 + 0.5d0*po*kg*g0*pa
      gx = tau*pa*(kg*g1*rh - 160.0d0)
      hx = tau*pa**2.0d0*kg*g0

*----------------------------------------------------------------------*
* Cubic coefficients 'p,q,r'.                                          *
*----------------------------------------------------------------------*
      p = sc*rd + (sc*ax*ex + cx*hx + dx*gx)/(cx*gx)
      q = (sc*rd*cx*hx + sc*rd*dx*gx + sc*ax*fx + sc*bx*ex +
     & dx*hx)/(cx*gx)
      r = (sc*rd*dx*hx + sc*bx*fx)/(cx*gx)

*----------------------------------------------------------------------*
* Sove the cubic equaiton to find 'a'.                                 *
*----------------------------------------------------------------------*
      Qt = (p**2 - 3.0d0*q)/9.0d0
      Rt = (2.0d0*p**3 - 9.0d0*p*q + 27.0d0*r)/54.0d0

      r2q3 = Rt**2 - Qt**3

      IF (r2q3.LT.0.0d0) THEN
        arg1 = Rt/Qt**1.5d0
        IF (arg1.GT.1.0d0)  arg1 = 1.0d0
        IF (arg1.LT.-1.0d0)  arg1 = -1.0d0
        th = acos(arg1)
        a1 =-2.0d0*Qt**0.5d0*cos(th/3.0d0) - p/3.0d0
        a2 =-2.0d0*Qt**0.5d0*cos((th + 2.0d0*pi)/3.0d0) - p/3.0d0
        a3 =-2.0d0*Qt**0.5d0*cos((th + 4.0d0*pi)/3.0d0) - p/3.0d0
        a = a1
        IF (a2.GT.a)  a = a2
        IF (a3.GT.a)  a = a3
      ELSE
        at =-Rt/abs(Rt)*(abs(Rt) + r2q3**0.5d0)**(1.0d0/3.0d0)
        IF (abs(at).GT.1E-6) THEN
          bt = Qt/at
        ELSE
          bt = 0.0d0
        ENDIF
        a1 = at + bt - p/3.0d0
*        a2 =-0.5d0*(at + bt) - p/3.0d0 + i*3.0d0**0.5d0*(at - bt)/2.0d0
*        a3 =-0.5d0*(at + bt) - p/3.0d0 - i*3.0d0**0.5d0*(at - bt)/2.0d0
        a = a1
      ENDIF

*----------------------------------------------------------------------*
* Compute 'gs,pc,w' corresponding to 'a'.                              *
*----------------------------------------------------------------------*
      gs = (g0 + g1*a*rh/pa)*kg
      pc = pa - a*160.0d0/gs
      w = j*pc/4.0d0/(pc + po/tau)
      acheck = (w*(1.0d0 - 0.5d0*po/(tau*pc)) - rd)*sc


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE XAMAX                            *
*                          ****************                            *
*                                                                      *
* This subroutine calculates the value of 'amax' at the maximum        *
* iradience as opposed to the average irradiance in the main program.  *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE XAMAX(a,xt,oi,vmx,rd,j,dresp,pc)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 a,xt,oi,vmx,rd,j,dresp,pc
      REAL*8 t,tk,c1,c2,kc,ko,tau,tpav,tpwv,tpaj,tpwj,sc
*      REAL*8 tppcj,tppcv,tpgsv,tpgsj,acheck

      t = 2.64d0 + 1.14d0*xt
      tk = 273.0d0 + t
      sc = 1000000.0d0

      c1 = 142.4d0 - 4.8d0*t
      IF (c1.GT.80.0d0)  c1 = 80.0d0
      IF (c1.LT.8.0d0)  c1 = 8.0d0
      c2 = 12.7d0 - 0.207d0*t
      IF (c2.GT.10.0d0)  c2 = 10.0d0
      IF (c2.LT.6.9d0)  c2 = 6.9d0
      kc = exp(35.8d0 - 80.5d0/(0.00831d0*tk))
      ko = exp(9.6d0 - 14.51d0/(0.00831d0*tk))*1000.0d0
      tau = exp(-3.949d0 + 28.99d0/(0.00831d0*tk))

*      rht = rh
*      IF (kg*c2*rh.LT.170.0d0)  rht = 170.0d0/(kg*c2)
*      CALL ASSVMAX(tpav,tpwv,tppcv,tpgsv,oi,ca,vmx,c1,c2,
*     &rht,ko,kc,kg,tau,rd,acheck)
*      CALL ASSJ(tpaj,tpwj,tppcj,tpgsj,oi,ca,j,c1,c2,rht,
*     &kg,tau,rd,acheck)

      tpwv = vmx*pc/(pc + kc*(1.0d0 + oi/ko))
      tpav = (tpwv*(1.0d0 - 0.5d0*oi/(tau*pc)) - rd)*sc

      tpwj = j*pc/(4.0d0*(pc + oi/tau))
      tpaj = (tpwj*(1.0d0 - 0.5d0*oi/(tau*pc)) - rd)*sc

      IF (tpav.LT.tpaj) THEN
        A = tpav
      ELSE
        A = tpaj
      ENDIF
      IF (A.LT.-dresp)  A =-dresp


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE ASSC4                            *
*                          ****************                            *
*                                                                      *
* ASS4 computes the assimilation rate by forming a cubic equation      *
* from the equations defining assimilation rate, stomatal conductance, *
* CO2 partial pressure (assimilation diffusion equaiton) and the       *
* carboxylation rate. The cubic equation is solved and the largest     *
* real root is assumed to be the assimilation rate. This is then       *
* compared with assimilation limited to light dependant and vmax       *
* dependant assimilation and the minimum one is taken. The variables   *
* here correspond to the Doly paper, in the main program in the        *
* folowing mannor.                                                     *
*        po,pa,j,g0,g1,rh,kg,tau,rd                                    *
*     = oi,ca,j(i),c1,c2,rh,kg,tau,rd                                  *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE ASSC4(a,pc,gs,pa,vmax,rh,kg,t,qg,p,up)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 a,pc,gs,pa,vmax,rh,t,qg,up,p,vmq,absx,qrub,f,jqq,dv
      REAL*8 amxt,ac4t,kg

      absx = 0.8d0
      qrub = 0.11d0
      f = 0.6d0
      vmq = 2.0d0

      dv = 0.6108d0*exp(17.269d0*t/(237.3d0 + t))*(1.0d0 - 
     &rh/100.0d0)/101.325d0
* kilo pascals to mol/mol = divide by 101.325

      pc = pa*(0.5d0 - (1.6d0*dv*p/2500.0d0/pa)**0.5d0)

      vmax = up*360.0d0*vmq**((t - 25.0d0)/10.0d0)/(360.0d0 + up)
      vmax = vmax/((1 + exp(0.3d0*(10.0d0 - t)))*(0.8d0 +
     &exp(0.14d0*(t - 36.0d0))))
      vmax = vmax*1.0e-6*1.25d0

      jqq = absx*qrub*f*qg*1.0e+6
      vmq = 2.0d0

      amxt = up*190.0d0/(360.0d0 + up)
      ac4t = amxt*vmq**((t - 25.0d0)/10.0d0)/(1.0d0 + 
     &exp(0.3d0*(10.0d0 - t)))/(0.8d0 + exp(0.14*(t - 36.0d0)))
      ac4t = up/2.0d0*vmq**((t - 25.0d0)/10.0d0)/(1.0d0 + 
     &exp(0.3d0*(15.0d0 - t)))/(0.8d0 + exp(0.14*(t - 33.0d0)))


      IF (ac4t.LT.jqq) THEN
        a = ac4t
      ELSE
        a = jqq
      ENDIF

      gs = a*1.6d0*p/(pa - pc)*1.0e-3


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE ASSC4                            *
*                          ****************                            *
*                                                                      *
* ASS4 computes the assimilation rate by forming a cubic equation      *
* from the equations defining assimilation rate, stomatal conductance, *
* CO2 partial pressure (assimilation diffusion equaiton) and the       *
* carboxylation rate. The cubic equation is solved and the largest     *
* real root is assumed to be the assimilation rate. This is then       *
* compared with assimilation limited to light dependant and vmax       *
* dependant assimilation and the minimum one is taken. The variables   *
* here correspond to the Doly paper, in the main program in the        *
* folowing mannor.                                                     *
*        po,pa,j,g0,g1,rh,kg,tau,rd                                    *
*     = oi,ca,j(i),c1,c2,rh,kg,tau,rd                                  *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE ASSC42(a,w,pc,gs,po,pa,vmax,g0,g1,rh,kg,tau,rd,t,qg,
     &acheck,up)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 a,w,pc,gs,po,pa,vmax,g0,g1,rh,kg,tau,rd,t,qg,acheck,up
      REAL*8 p,q,r,qt,rt,th,a1,a2,a3,pi,at,bt,sc,kp,kmq,pg,l,vmq,absx
      REAL*8 qrub,f,a0,b0,b1,c1,d1,e1,f1,b2,c2,d2,b3,c3,d3,wq,wv
      REAL*8 a4,b4,c4,d4,f2,aq,av,a2t,r2q3,arg1

      pi = atan(1.0d0)*4.0d0
      sc = 1000000.0d0

      kp = 0.7d0
      kmq = 2.0d0
      pg = 101325.0d0
      l = 0.000005d0
      vmq = 2.0d0
      absx = 0.8d0
      qrub = 0.11d0
      f = 0.6d0

      vmax = up*360.0d0*vmq**((t - 25.0d0)/10.0d0)/(360.0d0 + up)
      vmax = vmax/((1 + exp(0.3d0*(10.0d0 - t)))*(0.8d0 +
     &exp(0.14d0*(t - 36.0d0))))
      vmax = vmax*1.0e-6*1.25d0

*----------------------------------------------------------------------*
* 'a' and 'b' are such that 'jc=a pc  - b'.                            *
*----------------------------------------------------------------------*
      a0 = kp*kmq**((t - 25.0d0)/10.0d0)/pg
      b0 = l*kmq**((t - 25.0d0)/10.0d0)/pg

*----------------------------------------------------------------------*
* Reduce equations to a cubic in 'a'.                                  *
*----------------------------------------------------------------------*
      a1 =-b0*sc
      b1 = a0*sc
      c1 =-0.5d0*po
      d1 = tau
      e1 =-tau*rd*sc
      f1 = tau

      a2 = g0*kg*(pa**2.0d0)
      b2 = g0*kg*pa
      c2 =-pa*160.0d0 + g1*rh*kg*pa
      d2 = g1*rh*kg

      a3 = a1*b2 + a2*b1
      b3 = a1*d2 + c2*b1
      c3 = c1*b2 + a2*d1
      d3 = c1*d2 + c2*d1

      f2 =-f1

      a4 = f2*c2*d2
      b4 = f2*(a2*d2 + b2*c2) + d3*b3 + e1*c2*d2
      c4 = f2*a2*b2 + c3*b3 + d3*a3 + e1*(a2*d2 + c2*b2)
      d4 = c3*a3 + e1*a2*b2

*----------------------------------------------------------------------*
* Cubic coefficients 'p,q,r'.                                          *
*----------------------------------------------------------------------*
      p = b4/a4
      q = c4/a4
      r = d4/a4

      a2t = a2
*----------------------------------------------------------------------*
* Sove the cubic equaiton to find 'a'.                                 *
*----------------------------------------------------------------------*
      Qt = (p**2 - 3.0d0*q)/9.0d0
      Rt = (2.0d0*p**3 - 9.0d0*p*q + 27.0d0*r)/54.0d0

      r2q3 = Rt**2 - Qt**3

      IF (r2q3.LT.0.0d0) THEN
        arg1 = Rt/Qt**1.5d0
        IF (arg1.GT.1.0d0)  arg1 = 1.0d0
        IF (arg1.LT.-1.0d0)  arg1 = -1.0d0
        th = acos(arg1)
        a1 =-2.0d0*Qt**0.5d0*cos(th/3.0d0) - p/3.0d0
        a2 =-2.0d0*Qt**0.5d0*cos((th + 2.0d0*pi)/3.0d0) - p/3.0d0
        a3 =-2.0d0*Qt**0.5d0*cos((th + 4.0d0*pi)/3.0d0) - p/3.0d0
        a = a1
        IF (a2.GT.a)  a = a2
        IF (a3.GT.a)  a = a3
      ELSE
        at =-Rt/abs(Rt)*(abs(Rt) + r2q3**0.5d0)**(1.0d0/3.0d0)
        IF (abs(at).GT.1E-6) THEN
          bt = Qt/at
        ELSE
          bt = 0.0d0
        ENDIF
        a1 = at + bt - p/3.0d0
*        a2 =-0.5d0*(at + bt) - p/3.0d0 + i*3.0d0**0.5d0*(at - bt)/2.0d0
*        a3 =-0.5d0*(at + bt) - p/3.0d0 - i*3.0d0**0.5d0*(at - bt)/2.0d0
        a = a1
      ENDIF

      gs = (g0 + g1*a*rh/pa)*kg
      pc = pa - a*160.0d0/gs
      w = a0*pc - b0
      acheck = (w*(1.0d0 - 0.5d0*po/(tau*pc)) - rd)*sc

      a2 = a2t

*----------------------------------------------------------------------*
* Compute light dependent assimilation rate.                           *
*----------------------------------------------------------------------*

      wq = absx*qrub*f*qg*1.0d0
      p = tau*c2
      q = tau*a2 + wq*0.5d0*po*sc*d2 - tau*c2*(wq*sc - rd*sc)
      r = wq*0.5d0*po*sc*b2 - tau*a2*(wq*sc - rd*sc)
      aq = (-q + (abs(q)**2.0d0 - 4.0d0*p*r)**0.5d0)/(2.0d0*p)

*----------------------------------------------------------------------*
* Compute 'vmax' dependent assimilation rate.                          *
*----------------------------------------------------------------------*

      wv = vmax*vmq**((t - 25.0d0)/10.0d0)/((1.0d0 + exp(0.3d0*
     &(10.0d0 - t)))*(0.8d0 + exp(0.14d0*(t - 36.0d0))))

      q = tau*a2 + wv*0.5*po*sc*d2 - tau*c2*(wv*sc - rd*sc)
      r = wv*0.5*po*sc*b2 - tau*a2*(wv*sc - rd*sc)
      av = (-q + (abs(q)**2.0d0 - 4.0d0*p*r)**0.5d0)/(2.0d0*p)

*----------------------------------------------------------------------*
* Find limiting assimilation value.                                    *
*----------------------------------------------------------------------*
      IF (aq.LT.a) THEN
        a = aq
        w = wq
        gs = (g0 + g1*a*rh/pa)*kg
        pc = pa - a*160.0d0/gs
      ENDIF

      IF (av.LT.a) THEN
        a = av
        w = wv
        gs = (g0 + g1*a*rh/pa)*kg
        pc = pa - a*160.0d0/gs
      ENDIF


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE XVJMAX                           *
*                          *****************                           *
*                                                                      *
* This subroutine calculates the value of 'amax' at the maximum        *
* iradience as opposed to the average irradiance in the main program.  *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE XVJMAX(vmx,jmx,j,xt,xq,sum,nup,oi,dresp)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 vmx,jmx,j,xt,xq,sum,nup,oi
      REAL*8 tk,shapex,maxx,conv,cc,ea,num,denom,up,kc,ko,tau,can,upt,am
      REAL*8 vm,jm,irc,q,t,aa,bb,dresp

      t = 2.64d0 + 1.14d0*xt
      tk = 273.0d0 + t
      q = 2.0d0*xq

      nup = 1.5d0*nup

      shapex = 40.8d0 + 0.01d0*t - t**2*0.002d0
      maxx = 0.738d0 - 0.002d0*t
      IF (nup.LE.0.0d0) THEN
         print *,'Problem in XVJMAX nup=',nup
         stop
      ENDIF
      conv = 97.4115d0 - 2.5039d0*log(nup)

      cc = 36.0d0*exp(-t*0.14d0) + 20.0d0
      ea = 81000.0d0*exp(-t*0.14d0) + 50000.0d0
      num = exp(shapex - conv/(0.00831d0*tk))
      denom = 1.0d0 + exp((maxx*tk - 205.9d0)/(0.00831d0*tk))
      up = num/denom

      kc = exp(35.8d0 - 80.5d0/(0.00831d0*tk))
      ko = exp(9.6d0 - 14.51d0/(0.00831d0*tk))*1000.0d0
      tau = exp(-3.949d0 + 28.99d0/(0.00831d0*tk))
      aa = 24.5d0/(24.5d0 + kc + kc*oi/ko)
      bb = 0.5d0/(24.5d0 + kc + kc*oi/ko)


      can = 1.0d0/sum
      upt = up*can
      am = upt*190.0d0/(360.0d0 + upt)

      dresp = exp(cc-(ea/(8.3144d0*tk)))*upt/50.0d0

      vm =-1.0d0*(am + 0.82d0)/(-1.0d0*aa + bb*oi/tau)
      jm = 29.1d0 + 1.64d0*vm
      vm = vm/1000000.0d0
      jm = jm/1000000.0d0

      IF (t.GT.6.34365d0) THEN
        jmx = jm*(1.0d0 + 0.0409d0*(t-25.0d0)
     &- 1.54E-3*(t - 25.0d0)**2 - 9.42E-5*(t - 25.0d0)**3)
      ELSE
        jmx = jm*0.312633d0
      ENDIF

      IF (t.GT.9.51718d0) THEN
        vmx = vm*(1.0d0 + 0.0505d0*(t - 25.0d0)
     &- 0.248E-3*(t - 25.0d0)**2 - 8.09E-5*(t-25.0d0)**3)
      ELSE
        vmx = vm*0.458928d0
      ENDIF

* Mark we have to discuss about this calculation for irc here...
      irc = q*exp(-0.5d0)
      j = 0.24d0*irc/(1.0d0 + (0.24d0**2)*(irc**2)/(jmx**2))
     &**0.5d0


      RETURN
      END

