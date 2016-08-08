*----------------------------------------------------------------------*
*                          CENTURY SUBROUTINES                         *
*----------------------------------------------------------------------*

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE WSPARAM2                         *
*                          *******************                         *
*                                                                      *
* WSPARAM2 sets some system parameters for both the water and soils    *
* models.                                                              *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE WSPARAM2(nupc,awl,kd,kx,nci,nfix,adp,sfc,sw,sswc,ts,tc,
     &tsi,bulk,dep,topsl,wtfc,wtwp,swc)
*----------------------------------------------------------------------*   
      IMPLICIT NONE
      REAL*8 nupc,awl(4),kd,kx,nci(4),nfix,adp(4),sfc(4),sw(4),sswc(4)
      REAL*8 ts,tc,tsi,bulk,dep,topsl,b,lsmp,smp,wtfc,wtwp,swc

      awl(1) = 0.05d0 ! Relative root density of the ith layer.
      awl(2) = 0.05d0 !                  "
      awl(3) = 0.03d0 !                  "
      awl(4) = 0.01d0 !                  "
      kd = 0.5d0    ! Fraction of excess water flowing to deep storage.
      kx = 0.2d0    !                   "                 stream.
      nci(1) = 0.04d0       ! NC ratio of litter pools.
      nci(2) = 0.06666667d0 !          "
      nci(3) = 0.04d0       !          "
      nci(4) = 0.06666667d0 !          "
      nupc = 1.0d0   ! Amount of nitrogen uptake ratio to litter.
      nfix = 0.5d0   ! Soil deposition plus N fixation (g/m2/y).

*----------------------------------------------------------------------*
*  Compute the weight based field capacity and wilting point           *
* (g water m^3 / g soil m^3).                                          *
*----------------------------------------------------------------------*
        b = 3.1d0 + (0.157d0*tc) - (0.003d0*ts)
        swc = 50.5d0 - (0.142d0*ts) - (0.037d0*tc)
        lsmp = 1.54d0 - (0.0095d0*ts) + (0.0063d0*tsi)
        smp = 10.0d0**lsmp
        wtfc = ((((306.0d0/smp)**(-1.0d0/b))*swc)/100.0d0)*0.78d0
        wtwp = ((((15300.0d0/smp)**(-1.0d0/b))*swc)/100.0d0)*0.78d0
c        water = (((306.0d0/smp)**(-1.0d0/b))*swc/100.d0)*10000.0d0*dep

*----------------------------------------------------------------------*
* Convert soil layer height into mm of water.                          *
*----------------------------------------------------------------------*
      adp(1) = topsl*10.0d0
      adp(2) = (dep/3.0d0 - topsl)*10.0d0
      adp(3) = dep*10.0d0/3.0d0
      adp(4) = dep*10.0d0/3.0d0
      sfc(1) = wtfc*bulk*adp(1)/1000000.0d0
*      sfc(1) = adp(1)
      sfc(2) = wtfc*bulk*adp(2)/1000000.0d0
      sfc(3) = wtfc*bulk*adp(3)/1000000.0d0
      sfc(4) = wtfc*bulk*adp(4)/1000000.0d0
      sw(1) = wtwp*bulk*adp(1)/1000000.0d0
      sw(2) = wtwp*bulk*adp(2)/1000000.0d0
      sw(3) = wtwp*bulk*adp(3)/1000000.0d0
      sw(4) = wtwp*bulk*adp(4)/1000000.0d0
      sswc(1) = swc*adp(1)/100.0d0
      sswc(2) = swc*adp(2)/100.0d0
      sswc(3) = swc*adp(3)/100.0d0
      sswc(4) = swc*adp(4)/100.0d0

*      CALL FCAP(ts,tc,adp,sfc,sw)


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE WSPARAM                          *
*                          ******************                          *
*                                                                      *
* WSPARAM sets some system parameters for both the water and soils     *
* models.                                                              *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE WSPARAM(l_b_and_c,wilt,field,sat,nupc,awl,kd,kx,nci,
     &nfix,adp,sfc,sw,sswc,sand,clay,silt,bulk,org,dep,topsl,fc,wp,swc,
     &l_parameter)
*----------------------------------------------------------------------*
      REAL*8 nupc,awl(4),kd,kx,nci(4),nfix,adp(4),sfc(4),sw(4),sswc(4)
      REAL*8 sand,clay,silt,bulk,dep,topsl,bc_wres,bc_wqs,ans1,fc,wp,swc
      REAL*8 org,wilt,field,sat,bc_l,bc_a,bc_res
      LOGICAL l_b_and_c,l_parameter
      INCLUDE 'param.inc'

*      org = 0.0d0

      awl(1) = 0.05d0 ! Relative root density of the ith layer.
      awl(2) = 0.05d0 !                  "
      awl(3) = 0.03d0 !                  "
      awl(4) = 0.01d0 !                  "
      kd = 0.5d0    ! Fraction of excess water flowing to deep storage.
      kx = 0.2d0    !                   "                 stream.
      nci(1) = 0.04d0       ! NC ratio of litter pools.
      nci(2) = 0.06666667d0 !          "
      nci(3) = 0.04d0       !          "
      nci(4) = 0.06666667d0 !          "
      nupc = 1.0d0   ! Amount of nitrogen uptake ratio to litter.
      nfix = 0.5d0   ! Soil deposition plus N fixation (g/m2/y).

      bc_wres = 0.0d0
      IF (l_parameter) THEN
        org = p_orgc
        awl(1) = p_awl(1)
        awl(2) = p_awl(2)
        awl(3) = p_awl(3)
        awl(4) = p_awl(4)
        bc_res = p_bc_res
        kx = p_kx
        kd = p_kd
      ENDIF

      IF (l_b_and_c) THEN

*************************************************************************
*                        Brookes Corey
*************************************************************************
        bc_wqs = 1.0-bulk/2.65
        ans1 = -0.475662 +0.005165*sand + 0.002066*silt -
     &  0.023327*clay - 0.040628*bulk +
     &  0.006824*org - 0.000136162*sand**2 -
     &  0.0000954716*silt**2 + 0.000298891*clay**2 -
     &  0.0637*bulk**2 - 0.00031679*org**2 +
     &  0.0000010388*sand**3 + 0.0000010299*silt**3 -
     &  0.0000013251*clay**3

        bc_l = 10**ans1

        ans1 = 0.4104 + 0.002684*sand + 0.006354*silt +
     &0.17766*bulk + 0.00013855*clay**2

        bc_a = 1/10**ans1

        wp = bc_wres + (bc_wqs - bc_wres)*(15300.0*bc_a)**(-bc_l)
        fc = bc_wres + (bc_wqs - bc_wres)*(51.0*bc_a)**(-bc_l)
        swc = 1.0-bulk/2.65

      ELSE

        wp = wilt
        fc = field
        swc = sat

      ENDIF

*----------------------------------------------------------------------*
* Convert soil layer height into mm of water.                          *
*----------------------------------------------------------------------*
      adp(1) = topsl*10.0d0
      adp(2) = (dep/3.0d0 - topsl)*10.0d0
      adp(3) = dep*10.0d0/3.0d0
      adp(4) = dep*10.0d0/3.0d0
      sfc(1) = fc*adp(1)
*      sfc(1) = adp(1)
      sfc(2) = fc*adp(2)
      sfc(3) = fc*adp(3)
      sfc(4) = fc*adp(4)
      sw(1) = wp*adp(1)
      sw(2) = wp*adp(2)
      sw(3) = wp*adp(3)
      sw(4) = wp*adp(4)
      sswc(1) = swc*adp(1)
      sswc(2) = swc*adp(2)
      sswc(3) = swc*adp(3)
      sswc(4) = swc*adp(4)

*      CALL FCAP(ts,tc,adp,sfc,sw)


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE FCAP                             *
*                          ***************                             *
*                                                                      *
* FCAP computes the field capacities and the wilting points of the     *
* soil layers, from the soil sand and clay contents.                   *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE FCAP(ts,tc,adp,sfc,sw)
      IMPLICIT NONE
      REAL*8 ts,tc,adp(4),sfc(4),sw(4)
      REAL*8 t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15
      REAL*8 x,y,class(12)
      INTEGER sumx

      x=ts*100.00d0
      y=tc*100.00d0

      t1 = x - 44.32d0
      t2 = y - 39.33d0
      t3 = y + x - 59.55d0
      IF ((t3.GE.0.0d0).AND.(t2.GE.0.0d0).AND.(t1.LE.0.0d0)) THEN
        sumx = 1
      ELSE
        t4 = y - 34.83d0
        IF ((t1.GE.0.0d0).AND.(t4.GE.0.0d0)) THEN
        sumx = 5
        ELSE
          IF (t2.GE.0.0d0) THEN
            sumx = 2
          ELSE
            t10 = y - 19.78d0
            t15 = y + x - 71.78d0
            IF ((t15.GE.0.0d0).AND.(t10.GE.0.0d0).AND.(t1.GE.0.0d0)) THE
     &N
              sumx = 6
            ELSE
              t5 = y - 26.97d0
              t6 = x - 21.59d0
              IF ((t5.GE.0.0d0).AND.(t6.GE.0.0d0)) THEN
                sumx = 4
              ELSE
	        IF (t5.GE.0.0d0) THEN
	          sumx = 3
                ELSE
	          t7 = y + x - 50.00d0
                  t11 = x - 52.27d0
                  t12 = y - 7.87d0
                  IF ((t7.GE.0.0d0).AND.(t12.GE.0.0d0).AND.(t11.LE.0.0d0
     &))THEN
	            sumx = 7
                  ELSE
	            t13 = y - 0.98d0*x + 69.04d0
	            IF ((t7.GE.0.0d0).AND.(t13.GE.0.0d0)) THEN
	              sumx = 12
                    ELSE
	              t14 = y - 2.31d0*x + 197.92d0
                      IF ((t13.LE.0.0d0).AND.(t14.GE.0.0d0)) THEN
	                sumx = 10
                      ELSE
		        IF (t14.LE.0.0d0) THEN
	                  sumx = 11
                        ELSE
	                  t8 = y + x - 20.00d0
                          t9 = y - 10.00d0
                          IF ((t9.GE.0.0d0).OR.(t8.GE.0.0d0)) THEN
                            sumx = 8
                          ELSE
			    sumx = 9
			  ENDIF
			ENDIF
		      ENDIF
		    ENDIF
		  ENDIF
		ENDIF
	      ENDIF
	    ENDIF
	  ENDIF
	ENDIF
      ENDIF

      class(1) = 0.468d0      ! Clay
      class(2) = 0.468d0      ! Silty Clay
      class(3) = 0.464d0      ! Silty Clay Loam
      class(4) = 0.465d0      ! Clay Loam
      class(5) = 0.406d0      ! Sandy Clay
      class(6) = 0.404d0      ! Sandy Clay Loam
      class(7) = 0.439d0      ! Loam
      class(8) = 0.476d0      ! Silty loam
      class(9) = 0.476d0      ! Silt
      class(10) = 0.421d0     ! Loamy Sand
      class(11) = 0.339d0     ! Sand
      class(12) = 0.434d0     ! Sandy Loam

      sfc(1) = class(sumx)*adp(1)*1.0d0
      sfc(2) = class(sumx)*(adp(2) - adp(1))*1.0d0
      sfc(3) = class(sumx)*(adp(3) - adp(2))*1.0d0
      sfc(4) = class(sumx)*(adp(4) - adp(3))*1.0d0
      sw(1) = 10.0d0
      sw(2) = 10.0d0
      sw(3) = 10.0d0
      sw(4) = 10.0d0


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE SETCENPAR                        *
*                          ********************                        *
*                                                                      *
* SETCENPAR sets the century parameters.                               *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE SETCENPAR(ts,tc,pet,ppt,tmp,h2o,flow2,slnr,rlnr,ar,fms,
     &fmr,scl,rcl,tm,ft,cal,cap,cas,csp,csa,kr)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 ts,tc,pet,ppt,tmp,h2o,flow2,slnr,rlnr,ar,fms,fmr,scl,rcl,
     &tm,ft,cal,cap,cas,csp,csa,kr(8)
      REAL*8 smois,eot,eow,lcs,lcr,t
      INCLUDE 'param.inc'

* ppt and pet are monthly totals

      IF (pet.GT.0.0) THEN
        smois = (h2o + ppt*1.0d0)/pet*0.2d0
      ELSE
	smois = 1000.0d0
      ENDIF
      ar =-eow(smois)*eot(tmp)
      if (tmp.lt.0.0d0) ar = 0.0d0

      fms = 0.85d0 - 0.018d0*slnr
      fmr = 0.85d0 - 0.018d0*rlnr

      IF (fms.LT.0.001d0)  fms = 0.001d0
      IF (fmr.LT.0.001d0)  fmr = 0.001d0

      lcs = exp(-3.0d0*scl)
      lcr = exp(-3.0d0*rcl)

      t = 1.0d0 - ts/100.0d0
      tm = (1.0d0 - 0.75d0*t)

      ft = 0.85d0 - 0.68d0*t
      cal = (flow2/18.0d0)*(0.01d0 + 0.04d0*ts/100.0d0)
      cap = 0.003d0 + 0.032d0*tc/100.0d0
      cas = (1.0d0 - cap - cal - ft)

      csp = 0.003d0 - 0.009d0*tc/100.0d0
      csa = 1.0d0 - csp - 0.55d0

      kr(1) = 3.9d0*ar*lcs*p_kscale
      kr(2) = 4.8d0*ar*lcr*p_kscale
      kr(3) = 7.3d0*ar*tm*p_kscale
      kr(4) = 6.0d0*ar*p_kscale
      kr(5) = 14.8d0*ar*p_kscale
      kr(6) = 18.5d0*ar*p_kscale
      kr(7) = 0.2d0*ar*p_kscale
      kr(8) = 0.0045d0*ar*p_kscale


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE SOILCLOSS                        *
*                          ********************                        *
*                                                                      *
* SOILCLOSS calculates carbon losses in the soil due to resperation    *
* and leaching, given the soil carbon pool values.                     *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE SOILCLOSS(c,tot,scl,kr,rcl,ft,cal)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 c(8),tot,scl,kr(8),rcl,ft,cal

      tot = -0.3d0*scl*kr(1)*c(1) - 0.6d0*(1.0d0 - scl)*kr(1)*c(1) -
     &0.55d0*(1 - rcl)*kr(2)*c(2) - 0.3d0*rcl*kr(2)*c(2) -
     &0.6d0*kr(5)*c(5) - 0.55d0*kr(6)*c(6) - 0.6d0*kr(4)*c(4) -
     &0.55d0*kr(7)*c(7) - 0.55d0*kr(8)*c(8) - ft*kr(3)*c(3) -
     &cal*kr(3)*c(3)


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE CDYN                             *
*                          ***************                             *
*                                                                      *
* CDYN uses the NAG routine D02EBF to solve the differential equaitons *
* given by the carbon dynamics.                                        *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE CDYN(c0,c1,t0,t1,slc,rlc,scl,rcl,kr,csa,cas,csp,cap,
     &cal,fms,fmr,ft,sresp,lch,year,yr0,yrf,speedc)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 c0(8),c1(8),t0,t1,slc,rlc,scl,rcl,kr(8),csa,cas,csp,cap,
     &cal,fms,fmr,sresp,lch
      REAL*8 t01,t11,cp(10),cjac(10,10),ft,sc(8)
      INTEGER i,year,yrf,yr0,j
      REAL*8 cslc,crlc,cscl,crcl,ckr(8),ccsa,ccas,ccsp,ccap,ccal,cfms,
     &cfmr,cft
      LOGICAL speedc
      COMMON /cder/cslc,crlc,cscl,crcl,ckr,ccsa,ccas,ccsp,ccap,ccal,
     &cfms,cfmr,cft

      cslc = slc
      crlc = rlc
      cscl = scl
      crcl = rcl
      ccsa = csa
      ccas = cas
      ccsp = csp
      ccap = cap
      ccal = cal
      cfms = fms
      cfmr = fmr
      cft = ft
      DO i=1,8
	c1(i)=c0(i)
        ckr(i) = kr(i)
      ENDDO

      t01 = t0
      t11 = t1
      CALL DER(c1,cp)
      CALL PDER(cjac)

      IF (speedc) THEN
        sc(7) = 10.0d0*1.0d0
        sc(8) = 500.0d0*1.0d0
        cp(7) = cp(7)*sc(7)
        cp(8) = cp(8)*sc(8)
        DO i=1,10
          DO j=7,8
            cjac(j,i) = cjac(j,i)*(sc(j) - (sc(j) - 1.0d0)/
     &real(yrf - yr0 + 1)*real(year - yr0 + 1))
          ENDDO
        ENDDO
      ENDIF

      CALL DIFFEQ(t01,t11,c1,cp,cjac,sresp,lch)

      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE DIFFEQ                           *
*                          *****************                           *
*                                                                      *
* DIFFEQ gives the Jacobean corresponding to the differential carbon   *
* equations.                                                           *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE DIFFEQ(t01,t11,c1,cp,cjac,sresp,lch)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 t01,t11,c1(8),cp(10),cjac(10,10),secder(10),delt,sresp,lch
      INTEGER i,j
 
      delt = t11 - t01
 
      DO i=1,8
        secder(i) = 0.0d0
        DO j=1,8
          secder(i) = secder(i) + cjac(i,j)*cp(j)
        ENDDO
        c1(i) = c1(i) + cp(i)*delt + secder(i)*delt**2.0d0/2.0d0
      ENDDO
 
      secder(9) = 0.0d0
      secder(10) = 0.0d0
      DO j=1,10
        secder(9) = secder(9) + cjac(9,j)*cp(j)
        secder(10) = secder(10) + cjac(10,j)*cp(j)
      ENDDO
      sresp = sresp + cp(9)*delt + secder(9)*delt**2.0d0/2.0d0
      lch = lch + cp(10)*delt + secder(10)*delt**2.0d0/2.0d0


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE DER                              *
*                          **************                              *
*                                                                      *
* DER returns the carbon derivative calculated from the carbon state.  *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE DER(c,f)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 c(8),f(10)
      REAL*8 cslc,crlc,cscl,crcl,ckr(8),ccsa,ccas,ccsp,ccap,ccal,cfms,
     &cfmr,cft
      COMMON /cder/cslc,crlc,cscl,crcl,ckr,ccsa,ccas,ccsp,ccap,ccal,
     &cfms,cfmr,cft

      f(1) = cslc*(1.0D0 - cfms) + ckr(1)*c(1)
      f(2) = crlc*(1.0D0 - cfmr) + ckr(2)*c(2)
      f(3) =-0.45D0*(1.0D0 - crcl)*ckr(2)*c(2) - 0.45D0*ckr(6)*c(6) -
     & 0.45D0*ckr(8)*c(8) - ckr(7)*c(7)*ccsa + ckr(3)*c(3)
      f(4) =-0.4D0*(1.0D0 - cscl)*ckr(1)*c(1) - 0.4D0*ckr(5)*c(5) +
     & ckr(4)*c(4)
      f(5) = cslc*cfms + ckr(5)*c(5)
      f(6) = crlc*cfmr + ckr(6)*c(6)
      f(7) =-0.7D0*cscl*ckr(1)*c(1) - 0.4D0*ckr(4)*c(4) -
     & 0.7D0*crcl*ckr(2)*c(2) - ckr(3)*c(3)*ccas + ckr(7)*c(7)
      f(8) =-ckr(3)*c(3)*ccap - ckr(7)*c(7)*ccsp + ckr(8)*c(8)
      f(9) =-(0.6-0.3d0*cscl)*ckr(1)*c(1) - (0.55-0.25d0*crcl)*
     &ckr(2)*c(2) - cft*ckr(3)*c(3) -
     & 0.6d0*ckr(4)*c(4) - 0.6d0*ckr(5)*c(5) - 0.55d0*ckr(6)*c(6) -
     & 0.55d0*ckr(7)*c(7) - 0.55d0*ckr(8)*c(8)
      f(10) = -ckr(3)*c(3)*ccal

      
      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE PDER                             *
*                          ***************                             *
*                                                                      *
* PDER gives the Jacobean corresponding to the differential carbon     *
* equations.                                                           *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE PDER(jac)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 jac(10,10)
      INTEGER i,j
      REAL*8 cslc,crlc,cscl,crcl,ckr(8),ccsa,ccas,ccsp,ccap,ccal,cfms,
     &cfmr,cft
      COMMON /cder/cslc,crlc,cscl,crcl,ckr,ccsa,ccas,ccsp,ccap,ccal,
     &cfms,cfmr,cft

      DO i=1,10
	DO j=1,10
	  jac(i,j) = 0.0D0
        ENDDO
      ENDDO
      jac(1,1) = ckr(1)
      jac(2,2) = ckr(2)
      jac(3,2) =-0.45D0*(1.0D0 - crcl)*ckr(2)
      jac(3,3) = ckr(3)
      jac(3,6) =-0.45D0*ckr(6)
      jac(3,7) =-ccsa*ckr(7)
      jac(3,8) =-0.45D0*ckr(8)
      jac(4,1) =-0.4D0*(1.0D0 - cscl)*ckr(1)
      jac(4,4) = ckr(4)
      jac(4,5) =-0.4D0*ckr(5)
      jac(5,5) = ckr(5)
      jac(6,6) = ckr(6)
      jac(7,1) =-0.7D0*cscl*ckr(1)
      jac(7,2) =-0.7D0*crcl*ckr(2)
      jac(7,3) =-ccas*ckr(3)
      jac(7,4) =-0.4D0*ckr(4)
      jac(7,7) = ckr(7)
      jac(8,3) =-ccap*ckr(3)
      jac(8,7) =-ccsp*ckr(7)
      jac(8,8) = ckr(8)
      jac(9,1) =-(0.6 - 0.3d0*cscl)*ckr(1)
      jac(9,2) =-(0.55 - 0.25d0*crcl)*ckr(2)
      jac(9,3) =-cft*ckr(3) 
      jac(9,4) =-0.6d0*ckr(4) 
      jac(9,5) =-0.6d0*ckr(5) 
      jac(9,6) =-0.55d0*ckr(6) 
      jac(9,7) =-0.55d0*ckr(7) 
      jac(9,8) =-0.55d0*ckr(8) 
      jac(9,9) = 0.0d0 
      jac(10,3) = -ckr(3)*ccal

      RETURN
      END


*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE CLSASI                           *
*                          *****************                           *
*                                                                      *
* This subroutine returns the persentages of soil clay sand and silt.  *
* These are estimated from the soil classification on the islscp data  *
* file 'soltex'. This results in the folowing classifications:         *
*                                                                      *
*        ISLSCP Index      Clay  Sand  Silt       Comments             *
*              1             7    80    13        loamy sand           *
*              2            12    62    26        sandy loam           *
*              3            18    42    40        loam                 *
*              4            27    63    10        sandy clay loam      *
*              5            30    35    35        clay loam            *
*              6             0     0     0        ice                  *
*              7            18    42    40        loam                 *
*              0             0     0     0        ocean                *
*----------------------------------------------------------------------*
      SUBROUTINE CLSASI(tex,tc,ts,tsi)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 tc,ts,tsi
      INTEGER tex

      IF ((tex.EQ.0).OR.(tex.EQ.6)) THEN
        tc = 0.0d0
        ts = 0.0d0
        tsi = 0.0d0
      ELSE
        IF (tex.EQ.1) THEN
          tc = 7.0d0
          ts = 80.0d0
        ELSEIF (tex.EQ.2) THEN
          tc = 12.0d0
          ts = 62.0d0
        ELSEIF (tex.EQ.3) THEN
          tc = 18.0d0
          ts = 42.0d0
        ELSEIF (tex.EQ.4) THEN
          tc = 27.0d0
          ts = 63.0d0
        ELSEIF (tex.EQ.5) THEN
          tc = 30.0d0
          ts = 35.0d0
        ELSEIF (tex.EQ.7) THEN
          tc = 18.0d0
          ts = 42.0d0
        ENDIF
        tsi = 100.0d0 - tc -ts
      ENDIF


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE NDYN                             *
*                          ***************                             *
*                                                                      *
* Calculate nitrogen dynamics from the carbon dynamics and the carbon  *
* flows calculated using FLOWS.                                        *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE NDYN(n0,n,minn,t0,t1,fl,c,flow1,flow2,ts,
     &ml,al,scn,nfix,cal,nci,site)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 n0(8),n(8),minn(3),t0,t1,fl(14),c(8),flow1,
     &flow2,ts,ml,al,scn,nfix,cal,nci(4)
      REAL*8 ncsm,nca,ncs,ncp,ncratsm,ncrata,ncrats,ncratp,soiln
      REAL*8 plantn,delt,x,totn,ratminn,f1,f2,nin
      REAL*8 az,bz,cz,dz,ez,fz,gz,hz,iz
      INTEGER i,site

      delt = t1 - t0
      soiln = minn(3)
      plantn = 100.0D0/(scn + 1.0D0)
*...check whether 'scn' is correct, 'scn' = surface carbon/surface
* nitrogen, should this be litter carbon/litter nitrogen.

      nca = 1.0D0/ncrata(soiln)

      DO 10 i=1,8
        n(i)=n0(i)
10    CONTINUE

*----------------------------------------------------------------------*
* Find mineral nitrogen leaching due to the water flows between the    *
* soil layers.                                                         *
*----------------------------------------------------------------------*
*      ml = minn(2)*flow2*delt*(0.2D0 + 0.7D0*ts/100.0d0)/18.0D0
*      minn(2) = minn(2) - ml
*      IF (minn(1).GT.0.0D0) THEN
*	temp = minn(1)*flow1*delt*(0.2D0 + 0.7D0*ts/100.0d0)/18.0D0
*	minn(1) = minn(1) - temp
*       minn(2) = minn(2) + temp
*      ENDIF

*      minn(1) = minn(1) + nfix*delt + (1.0D0 - nup)*(sln + rln)*delt

      IF (flow1.LT.1.0e-10) flow1 = 0.0d0
      IF (flow2.LT.1.0e-10) flow2 = 0.0d0

      f1 = flow1*(0.2D0 + 0.7D0*ts/100.0d0)/18.0D0
      f2 = flow2*(0.2D0 + 0.7D0*ts/100.0d0)/18.0D0
      nin = nfix
      IF (abs(f1-f2).LT.1.0e-6) THEN
        f2 = f2 - 0.001d0
      ENDIF

      IF (flow1.GT.0.0d0) THEN
        bz = nin/f1
        az = minn(1) - bz
        minn(1) = az*exp(-f1*delt) + bz
        IF (flow2.GT.0.0d0) THEN
          dz = f1*az/(f2 - f1)
          ez = f1*bz/f2
          cz = minn(2) - dz - ez
          minn(2) = cz*exp(-f2*delt) + dz*exp(-f1*delt) + ez
          fz = -cz
          gz = -dz*f2/f1
          hz = f2*ez
          iz = -fz - gz
          ml = fz*exp(-f2*delt) + gz*exp(-f1*delt) + hz*delt + iz
        ELSE
          cz = -az
          dz = f1*bz
          ez = minn(2) - cz
          minn(2) = cz*exp(-f1*delt) + dz*delt + ez
          fz = -cz*f2/f1
          gz = f2*dz/2.0d0
          hz = f2*ez
          iz =-fz
          ml = fz*exp(-f1*delt) + gz*delt**2.0d0 + hz*delt + iz
        ENDIF
      ELSE
        IF (flow2.GT.0.0d0) THEN
          minn(1) = nfix*delt + minn(1)
          cz = minn(2)
          minn(2) = cz*exp(-f2*delt)
          ml = -cz*exp(-f2*delt) + cz
        ELSE
          minn(1) = nfix*delt + minn(1)
          minn(2) = minn(2)
        ENDIF
      ENDIF

      IF (.NOT.(minn(2).LT.10000.0d0))  WRITE(11,
     &'(''minn(2) > 10000, '',i5,7f12.4)') site,minn(2),cz,f2,dz,f1,ez,
     &delt

      IF (minn(3).LT.2.0d0) THEN
        al = fl(12)*cal/(15.0d0 - 6.0d0*minn(3))
      ELSE
        al = fl(12)*cal/3.0d0
      ENDIF

*----------------------------------------------------------------------*
* Calculate total nitrogen and then find how much of that is mineral.  *
*----------------------------------------------------------------------*
      totn = minn(1) + minn(2) + n(1) + n(2) + n(3) + n(4) + n(5) +
     & n(6) + n(7) + n(8) - al

      IF (.NOT.(minn(1).LT.10000.0d0))  WRITE(11,
     &'(''minn(1) > 10000, '',i5,4f12.4)') site,flow1,flow2,minn(1),
     &minn(2)

      CALL NBAL(c,totn,x,scn,nci)
      IF (.NOT.(minn(1).LT.10000.0d0)) WRITE(11,
     &'(''minn(1) > 10000, '',i5,4f12.4)')  site,flow1,flow2,minn(1),
     &minn(2)

*----------------------------------------------------------------------*
* Adjust all nitrogen pools accordingly.                               *
*----------------------------------------------------------------------*
      IF (minn(2).GT.0.0d0) THEN
        ratminn = minn(1)/minn(2)
        minn(1) = x/(1.0d0 + 1.0d0/ratminn)
        minn(2) = x/(1.0d0 + ratminn)
        minn(3) = x
      ELSE
        minn(1) = x
        minn(2) = 0.0d0
        minn(3) = x
      ENDIF

      IF (.NOT.(minn(1).LT.10000.0d0))  WRITE(11,
     &'(''minn(1) > 10000, '',i5,6f12.4)') site,flow1,flow2,minn(1),
     &minn(2),x,ratminn

      ncsm = 1.0D0/ncratsm(plantn)
      nca = 1.0D0/ncrata(x)
      ncs = 1.0D0/ncrats(x)
      ncp = 1.0D0/ncratp(x)

      n(1) = c(1)*nci(1)
      n(5) = c(5)*nci(2)
      n(2) = c(2)*nci(3)
      n(6) = c(6)*nci(4)
      n(4) = ncsm*c(4)
      n(8) = ncp*c(8)
      n(7) = ncs*c(7)
      n(3) = nca*c(3)
*----------------------------------------------------------------------*


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE NBAL                             *
*                          ***************                             *
*                                                                      *
* NBAL computes the mineral nitrogen pool, given the value of the      *
* total nitrogen content of the soil. To do this a quartic equation    *
* needs to be solved, this is done using the Newton Raphson method.    *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE NBAL(c,tn,x,scn,nci)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 c(8),tn,x,scn,nci(4)
      REAL*8 const,ncsm,ncratsm,plantn,fx0,fx05,fx1,fx15,fx2,b,at
      INTEGER i

      plantn = 100.0D0/(scn + 1.0D0)
      ncsm = 1.0D0/ncratsm(plantn)

      const = tn - c(1)*nci(1) - c(5)*nci(2) - c(2)*nci(3) - c(6)*nci(4)
     &- c(4)*ncsm

*      cx1 =-7.2
*      cx2 = 102.0 + 7.2*const
*      cx3 = 4.8*c(8) + 1.8*c(7) + 1.2*c(3) - 450.0 - 102.0*const
*      cx4 =-36.0*c(8) - 16.5*c(7) - 14.0*c(3) + 600.0 + 450.0*const
*      cx5 = 60.0*c(8) + 30.0*c(7)  + 40.0*c(3) - 600.0*const

      x=0.0D0
      fx0 = 60.0*c(8) - 36.0*c(8)*x + 4.8*c(8)*x**2 + 30.0*c(7) -
     & 16.5*c(7)*x +1.8*c(7)*x**2 + 40.0*c(3) - 14.0*c(3)*x +
     & 1.2*c(3)*x**2 + 600.0*x - 450.0*x**2 + 102.0*x**3 - 7.2*x**4 -
     & 600.0*const + 450.0*const*x - 102.0*const*x**2 +
     & 7.2*const*x**3
      x=0.5D0
       fx05 = 60.0*c(8) - 36.0*c(8)*x + 4.8*c(8)*x**2 + 30.0*c(7) -
     & 16.5*c(7)*x +1.8*c(7)*x**2 + 40.0*c(3) - 14.0*c(3)*x +
     & 1.2*c(3)*x**2 + 600.0*x - 450.0*x**2 + 102.0*x**3 - 7.2*x**4 -
     & 600.0*const + 450.0*const*x - 102.0*const*x**2 +
     & 7.2*const*x**3
      x=1.0D0
       fx1 = 60.0*c(8) - 36.0*c(8)*x + 4.8*c(8)*x**2 + 30.0*c(7) -
     & 16.5*c(7)*x +1.8*c(7)*x**2 + 40.0*c(3) - 14.0*c(3)*x +
     & 1.2*c(3)*x**2 + 600.0*x - 450.0*x**2 + 102.0*x**3 - 7.2*x**4 -
     & 600.0*const + 450.0*const*x - 102.0*const*x**2 +
     & 7.2*const*x**3
      x=1.5D0
       fx15 = 60.0*c(8) - 36.0*c(8)*x + 4.8*c(8)*x**2 + 30.0*c(7) -
     & 16.5*c(7)*x +1.8*c(7)*x**2 + 40.0*c(3) - 14.0*c(3)*x +
     & 1.2*c(3)*x**2 + 600.0*x - 450.0*x**2 + 102.0*x**3 - 7.2*x**4 -
     & 600.0*const + 450.0*const*x - 102.0*const*x**2 +
     & 7.2*const*x**3
      x=2.0D0
       fx2 = 60.0*c(8) - 36.0*c(8)*x + 4.8*c(8)*x**2 + 30.0*c(7) -
     & 16.5*c(7)*x +1.8*c(7)*x**2 + 40.0*c(3) - 14.0*c(3)*x +
     & 1.2*c(3)*x**2 + 600.0*x - 450.0*x**2 + 102.0*x**3 - 7.2*x**4 -
     & 600.0*const + 450.0*const*x - 102.0*const*x**2 +
     & 7.2*const*x**3

      IF ((fx0.LT.0.0).AND.(fx2.GT.0.0D0)) THEN
	IF (fx15.LT.0.0D0) THEN
	  x = 0.5D0*fx15/(fx15 - fx2) + 1.5D0
	ELSEIF (fx1.LT.0.0D0) THEN
	  x = 0.5D0*fx1/(fx1 - fx15) + 1.0D0
	ELSEIF (fx05.LT.0.0D0) THEN
	  x = 0.5D0*fx05/(fx05 - fx1) + 0.5D0
	ELSE
	  x = 0.5D0*fx0/(fx0 - fx05)
	ENDIF
	DO 10 i=1,4
          at = 60.0D0*c(8) - 36.0D0*c(8)*x + 4.8D0*c(8)*x**2 +
     & 30.0D0*c(7) - 16.5D0*c(7)*x +1.8D0*c(7)*x**2 + 40.0D0*c(3) -
     & 14.0D0*c(3)*x + 1.2D0*c(3)*x**2 + 600.0D0*x - 450.0D0*x**2 +
     & 102.0D0*x**3 - 7.2D0*x**4 - 600.0D0*const + 450.0D0*const*x -
     & 102.0D0*const*x**2 + 7.2D0*const*x**3
	b =-36.0D0*c(8) + 9.6D0*c(8)*x - 16.5D0*c(7) + 3.6D0*c(7)*x -
     & 14.0D0*c(3) + 2.4D0*c(3)*x + 600.0D0 - 900.0D0*x + 306.0D0*x**2 -
     & 28.8D0*x**3 + 450.0D0*const - 204.0D0*const*x + 21.6D0*const*x**2
	  x = x - at/b
10      CONTINUE
      ELSE IF ((fx0.GT.0.0).AND.(fx2.GT.0.0D0)) THEN
	x = 0.0D0
      ELSE IF ((fx0.LT.0.0).AND.(fx2.LT.0.0D0)) THEN
	x = const - c(3)/3.0D0 - c(7)/12.0D0 - c(8)/7.0D0
      ENDIF


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SOUBROUTINE FLOWS                           *
*                          *****************                           *
*                                                                      *
* FLOWS computes the fourteen carbon flows associated with the carbon  *
* dynamics.                                                            *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE FLOWS(c0,c,t0,t1,fl,slc,rlc,scl,rcl,fms,fmr,csa,csp,
     &cas,cap)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 c0(8),c(8),t0,t1,fl(14),slc,rlc,scl,rcl,fms,fmr,csa,csp,
     &cas,cap
      REAL*8 delt

      delt = t1 - t0

      fl(1) = slc*(1.0D0 - fms)*delt
      fl(2) = slc*fms*delt
      fl(7) = rlc*(1.0D0 - fmr)*delt
      fl(8) = rlc*fmr*delt
      fl(3) = (c0(1) - c(1) + fl(1))*scl
      fl(4) = (c0(1) - c(1) + fl(1))*(1.0D0 - scl)
      fl(5) = c0(5) - c(5) + fl(2)
      fl(9) = (c0(2) - c(2) + fl(7))*rcl
      fl(10) = (c0(2) - c(2) + fl(7))*(1.0D0 - rcl)
      fl(11) = c0(6) - c(6) + fl(8)
      fl(6) = c0(4) - c(4) + 0.4D0*fl(4) + 0.4D0*fl(5)

*----------------------------------------------------------------------*
* Solve three simultaneous equations to find fl(12),fl(13) and fl(14). *
*----------------------------------------------------------------------*
      fl(12) = (c0(3) - c(3) + 0.45D0*fl(10) + 0.45D0*fl(11) +
     & csa*(c0(7) - c(7) + 0.7D0*fl(3) + 0.4D0*fl(6) + 0.7D0*fl(9)) +
     & 0.45D0*(c0(8) - c(8) + csp*(c0(7) - c(7) + 0.7D0*fl(3) +
     & 0.4D0*fl(6) + 0.7D0*fl(9))))/(1.0D0 - csa*cas - 0.45D0*csp*cas -
     & 0.45D0*cap)
      fl(13) = c0(7) - c(7) + 0.7D0*fl(3) + 0.4D0*fl(6) + 0.7D0*fl(9) +
     & cas*fl(12)
      fl(14) = c0(8) - c(8) + csp*fl(13) + cap*fl(12)

*     WRITE(*,*) 'c1',c(1)-c0(1),fl(1)-fl(3)-fl(4)
*     WRITE(*,*) 'c2',c(2)-c0(2),fl(7)-fl(9)-fl(10)
*     WRITE(*,*) 'c3',c(3)-c0(3),fl(13)*csa+fl(10)*.45D0+fl(11)*.45D0+
*   & fl(14)*.45D0-fl(12)
*     WRITE(*,*) 'c4',c(4)-c0(4),.4D0*fl(4)+.4D0*fl(5)-fl(6)
*     WRITE(*,*) 'c5',c(5)-c0(5),fl(2)-fl(5)
*     WRITE(*,*) 'c6',c(6)-c0(6),fl(8)-fl(11)
*     WRITE(*,*) 'c7',c(7)-c0(7),.7D0*fl(3)+.4D0*fl(6)+.7D0*fl(9)+
*   & cas*fl(12)-fl(13)
*     WRITE(*,*) 'c8',c(8)-c0(8),csp*fl(13)+cap*fl(12)-fl(14)


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE ICCODE                           *
*                          *****************                           *
*                                                                      *
* This subroutine writes out to the output files with iccode values,   *
* these are as follows:                                                *
*                                                                      *
*                   -1 = snow                                          *
*                   -2 = desert                                        *
*                   -3 = insufficient data                             *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE ICCODE(lat,lon,c)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 lat,lon,c

      WRITE(40,1540) lat, lon, c
      WRITE(41,1540) lat, lon, c
      WRITE(42,1542) lat, lon, c, c

1540  FORMAT(F6.2,1X,F7.2,2X,F6.3)
1542  FORMAT(F6.2,1X,F7.2,2X,F8.1,1X,F7.1)

      RETURN
      END

*----------------------------------------------------------------------*
*                          FUNCTIONS                                   *
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*                                                                      *
*                          FUNCTION NCRATA                             *
*                          ***************                             *
*                                                                      *
* 'ncrata' calculates N:C for the active organic matter pool.          *
*                                                                      *
*----------------------------------------------------------------------*
      FUNCTION ncrata(x)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 ncrata,x

      IF (x.GT.2.0d0) THEN
        ncrata = 3.0d0
      ELSE
	ncrata = 15.0d0 - 6.0d0*x
      ENDIF
      IF (x.LT.0.0d0) THEN
	ncrata = 15.0d0
      ENDIF


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          FUNCTION NCRATS                             *
*                          ***************                             *
*                                                                      *
* 'ncrats' calculates N:C for the slow organic matter pool.            *
*                                                                      *
*----------------------------------------------------------------------*
      FUNCTION ncrats(x)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 ncrats,x

      IF (x.GT.2.0d0) THEN
        ncrats = 12.0d0
      ELSE
	ncrats = 20.0d0 - 4.0d0*x
      ENDIF
      IF (x.LT.0.0d0) THEN
	ncrats = 20.0d0
      ENDIF


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          FUNCTION NCRATP                             *
*                          ***************                             *
*                                                                      *
* 'ncratp' calculates N:C for the passive organic matter pool.         *
*                                                                      *
*----------------------------------------------------------------------*
      FUNCTION ncratp(x)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 ncratp,x

      IF (x.GT.2.0d0) THEN
        ncratp = 7.0d0
      ELSE
	ncratp = 10.0d0 - 1.5d0*x
      ENDIF
      IF (x.LT.0.0d0) THEN
	ncratp = 10.0d0
      ENDIF


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          FUNCTION NCRATSM                            *
*                          ****************                            *
*                                                                      *
* 'ncratsm' calculates N:C for the surface microbe pool.               *
*                                                                      *
*----------------------------------------------------------------------*
      FUNCTION ncratsm(x)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 ncratsm,x

      IF (x.GT.2.0d0) THEN
        ncratsm = 10.0d0
      ELSE
      	ncratsm = 20.0d0 - 5.0d0*x
      ENDIF
      IF (x.LT.0.0d0) THEN
	ncratsm = 20.0d0
      ENDIF


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          FUNCTION EOT                                *
*                          ************                                *
*                                                                      *
* 'eot' calculates the effect of temperature used to calculate the     *
* decay rate multiplier 'a'.                                           *
*                                                                      *
*----------------------------------------------------------------------*
      FUNCTION eot(x)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 x,eot

      IF (x.LT.0.0d0) THEN
	eot = 0.04
      ELSEIF (x.LT.10.0d0) THEN
	eot = (0.16d0 - 0.03d0)*(x - 0.0d0)/(10.0d0 - 0.0d0) + 0.03d0
      ELSEIF (x.LT.20.0d0) THEN
	eot = (0.50d0 - 0.16d0)*(x - 10.0d0)/(20.0d0 - 10.0d0) + 0.16d0
      ELSEIF (x.LT.30.0d0) THEN
	eot = (0.88d0 - 0.50d0)*(x - 20.0d0)/(30.0d0 - 20.0d0) + 0.50d0
      ELSEIF (x.LT.40.0d0) THEN
	eot = (0.96d0 - 0.88d0)*(x - 30.0d0)/(40.0d0 - 30.0d0) + 0.88d0
      ELSEIF (x.LT.47.0d0) THEN
	eot = (0.86d0 - 0.96d0)*(x - 40.0d0)/(47.0d0 - 40.0d0) + 0.96d0
      ELSE
	eot = 0.0d0
      ENDIF


      RETURN
      END


*----------------------------------------------------------------------*
*                                                                      *
*                          FUNCTION EOT1                               *
*                          *************                               *
*                                                                      *
* 'eot1' calculates the effect of temperature used to calculate the    *
* decay rate multiplier 'a'.                                           *
*                                                                      *
*----------------------------------------------------------------------*
      FUNCTION eot1(x)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 x,eot1,a,b

      a = (45.0d0 - x)/10.0d0
      b = exp(0.076d0*(1.0d0 - a**2.63))
      eot1 = (a**0.2d0)*b


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          FUNCTION EOW                                *
*                          ************                                *
*                                                                      *
* 'eow' calculates the effect of water used to calculate the decay     *
* rate multiplier 'a'.                                                 *
*                                                                      *
*----------------------------------------------------------------------*
      FUNCTION eow(x)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 x,eow

      IF (x.LT.0.2d0) THEN
	eow = (0.21d0 - 0.09d0)*(x - 0.0d0)/(0.2d0 - 0.0d0) + 0.09d0
      ELSEIF (x.LT.0.4d0) THEN
	eow = (0.48d0 - 0.21d0)*(x - 0.2d0)/(0.4d0 - 0.2d0) + 0.21d0
      ELSEIF (x.LT.0.6d0) THEN
	eow = (0.82d0 - 0.48d0)*(x - 0.4d0)/(0.6d0 - 0.4d0) + 0.48d0
      ELSEIF (x.LT.0.8d0) THEN
	eow = (0.97d0 - 0.82d0)*(x - 0.6d0)/(0.8d0 - 0.6d0) + 0.82d0
      ELSEIF (x.LT.1.0d0) THEN
	eow = (0.97d0 - 0.97d0)*(x -0.8d0)/(1.0d0 - 0.8d0) + 0.97d0
      ELSEIF (x.LT.2.0d0) THEN
	eow = (0.63d0 - 0.97d0)*(x - 1.0d0)/(2.0d0 - 1.0d0) + 0.97d0
      ELSE
	eow = (0.29d0 - 0.63d0)*(x - 2.0d0)/(3.0d0 - 2.0d0) + 0.63d0
      ENDIF

      IF ((eow.LT.0.0d0).AND.(x.GT.1.0d0)) THEN
	  eow = 0.1d0
      ENDIF

      RETURN
      END


*----------------------------------------------------------------------*
*                                                                      *
*                          FUNCTION EOW1                               *
*                          *************                               *
*                                                                      *
* 'eow1' calculates the effect of water used to calculate the decay    *
* rate multiplier 'a'.                                                 *
*                                                                      *
*----------------------------------------------------------------------*
      FUNCTION eow1(x)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 x,eow1

      IF (x.LT.1.5d0) THEN
        eow1 = 1.0d0/(1.0d0 + 30.0d0*exp(-8.5d0*x))
      ELSE
        eow1 = 1.0d0 - 0.7d0*(x - 1.5d0)/1.5d0
      ENDIF

      IF (eow1.LT.0.0001d0)  eow1 = 0.0001d0


      RETURN
      END

