      REAL*8 a,ci,cs,gs,oi,ca,vmx,rh,kg,ko,kc,tau,rd,ga,t,p,j

      t = 14.0d0
      kg = 1.0d0
      rh = 60.0d0
      vmx = 23.0e-6

      p = 
      ga = 

      tk = 273.0d0 + t
      kc = exp(35.8d0 - 80.5d0/(0.00831d0*tk))
      ko = exp(9.6d0 - 14.51d0/(0.00831d0*tk))*1000.0d0
      tau = exp(-3.949d0 + 28.99d0/(0.00831d0*tk))

      CALL ASSVMAX(a,ci,cs,gs,oi,ca,vmx,rh,kg,ko,kc,tau,rd,ga,t,p)

      CALL ASSJ(a,ci,cs,gs,oi,ca,rh,kg,tau,rd,ga,t,p,j)


      stop
      end

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
      REAL*8 a,ci,gs,oi,ca,vmx,rh,ko,kc,tau,rd,fa0,faf,gam,step
      REAL*8 p,ga,x,y,dv,t,a0,cs,af,fa,kg,ax,bx,cx,gsmin
      INTEGER i

      gsmin = 0.005d0

      x = -1.5d0 + 0.2d0*t
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
      SUBROUTINE ASSJ(a,ci,cs,gs,oi,ca,rh,kg,tau,rd,ga,t,p,j)
*----------------------------------------------------------------------*
      REAL*8 a,ci,gs,oi,ca,rh,tau,rd,j,gam,fa0,faf,a0,af,step
      REAL*8 p,ga,x,y,dv,t,cs,fa,kg
      REAL*8 ax,bx,cx,gsmin

      gsmin = 0.005d0

      x = -1.5d0 + 0.2d0*t
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
