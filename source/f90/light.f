*----------------------------------------------------------------------*
*                          FUNCTIONS                                   *
*----------------------------------------------------------------------*
*                          DAYLENGTH                                   *
*----------------------------------------------------------------------*
      FUNCTION dayl(lat,day)
*----------------------------------------------------------------------*
      REAL*8 dayl,lat,del,has,tem,conv
      INTEGER day
      PARAMETER (conv = 1.74532925E-2)

      del=-23.4d0*cos(conv*360.0d0*(day + 10.0d0)/365.0d0)
      tem=-tan(lat*conv)*tan(del*conv)

      IF (abs(tem).LT.1.0d0) THEN
	has  = acos(tem)/conv
        dayl = 2.0d0*has/15.0d0
      ELSEIF (tem.GT.0.0d0) THEN
	dayl =  0.0d0
      ELSE
	dayl = 24.0d0
      ENDIF


      RETURN
      END


*----------------------------------------------------------------------*
*                          PHOTON FLUX DENSITY                         *
*----------------------------------------------------------------------*
      SUBROUTINE pfd(lat,day,hrs,cloud,direct,diffuse,total)
*     Ghislain 08/12/03
*----------------------------------------------------------------------*
      REAL*8 lat,hrs,del,rlat,toa,pi,conv,cloud,diffprop,alpha,beta
      REAL*8 sigma,dawn_angle,coscst,direct,diffuse,total,clearness
      INTEGER day
      PARAMETER(conv = 1.74532925E-2, pi = 3.1415927)

      IF (hrs.GT.1E-6) THEN
        del  = -23.4d0*cos(conv*360.0d0*(day+10.0d0)/365.0d0)
        del  = del*conv
        rlat = lat*conv

*     solar radiation in watts/m2 average over the day 
*     (between sun rise en sun set)

        dawn_angle = (hrs/2.0d0)*(2.0d0*pi/24.0d0)

c     top of atmosphere irradiance
        toa = 1370.0d0*   ! solar constant in W/m2
     &       24.0d0/(pi*hrs)*(      ! average over the day (daylength)
     &       sin(rlat)*sin(del)*dawn_angle + 
     &       cos(rlat)*cos(del)*sin(dawn_angle))

c     account for cloud cover

        coscst = 1.0d0-1.3614d0*cos(rlat)
        alpha = 18.0d0 - 64.884d0*coscst
        alpha = alpha*4.1842d0*1.0d4 ! convert from cal/cm to J/m2
        alpha = alpha/(24.0d0*3600.0d0) ! convert into w/m2
        
        beta = 0.682d0 - 0.3183d0*coscst

        sigma = 0.02d0*log(max(cloud,0.001d0))+0.03259d0

c     calculate the total irradiance
        total = toa*(beta-sigma*cloud*10.0d0)-alpha
        IF (total.LT.0.0d0) total=0.0d0

        clearness= total/toa

c     calculate diffuse irradiance (from Forest ETP WG)

        IF(clearness.LT.0.07) THEN
           diffprop = 1.0
        ELSE IF(clearness.LT.0.35) THEN
           diffprop = 1.0-2.3*(clearness-0.07)**2
        ELSE IF(clearness.LT.0.75) THEN
           diffprop = 1.33-1.46*clearness
        ELSE
           diffprop = 0.23d0
        ENDIF

        diffuse=total*diffprop

c PAR is about 48% of the total irradiance. A better formula could be 
c found.

        total = 0.48d0*total
        diffuse = 0.48d0*diffuse ! not sure this is really right
                                 ! ... (Ghislain)

c convert in mol/m2/sec 1W/m2=4.6e-6 mol/s/m2 for PAR radiation

        total = 4.6d-6*total
        diffuse = 4.6d-6*diffuse

        direct= total - diffuse
      ELSE
        total=0.0d0
        diffuse=0.0d0
        direct=0.0d0
      ENDIF


      RETURN
      END







*----------------------------------------------------------------------*
*                          PHOTON FLUX DENSITY                         *
*----------------------------------------------------------------------*
      FUNCTION pfd_without_cloud(lat,day,hrs)
*----------------------------------------------------------------------*
      REAL*8 pfd_without_cloud
      REAL*8 lat,hrs,del,rlat,r,sinbt,pi,conv,cloud,st,sdiff,trans
      INTEGER day
      PARAMETER(conv = 1.74532925E-2, pi = 3.1415927)

      cloud = 0.50d0

      IF (hrs.GT.1E-6) THEN
        del  = -23.4d0*cos(conv*360.0d0*(day+10.0d0)/365.0d0)
        del  = del*conv
        rlat = lat*conv

        trans = 0.3d0 + (0.35d0*(1.0d0 - cloud))

* solar radiation in watts/m2 at midday
        sinbt = sin(rlat)*sin(del) + cos(rlat)*cos(del)

        IF (sinbt.GT.1E-6) THEN
          r = 1370.0d0*trans**(1.0d0/sinbt)*sinbt

          sdiff = 0.3d0*(1.0d0 - trans**(1.0d0/sinbt))*1370.0d0*sinbt

          st = r + sdiff

          pfd_without_cloud = 
     &         2.0d0*st/pi*(0.5d0 + cos(rlat)/2.0d0)*2.1d0/1.0e+6

        ELSE
         pfd_without_cloud = 0.0d0
        ENDIF
      ELSE
        pfd_without_cloud = 0.0d0
      ENDIF


      RETURN
      END


*----------------------------------------------------------------------*
*                          PHOTON FLUX DENSITY                         *
*----------------------------------------------------------------------*
      FUNCTION pfd2(lat,day,hrs)
*----------------------------------------------------------------------*
      REAL*8 pfd2,lat,hrs,del,rlat,r,sinbt,pi,conv,rad
      INTEGER day
      PARAMETER(conv = 1.74532925E-2, pi = 3.1415927)

      IF (hrs.GT.1E-6) THEN
	del  = -23.4d0*cos(conv*360.0d0*(day+10.0d0)/365.0d0)
	del  = del*conv
	rlat = lat*conv

* solar radiation in watts/m2 at midday
        sinbt = sin(rlat)*sin(del) + cos(rlat)*cos(del)
        r = 1360.0d0*0.7d0**1.32d0*sinbt

* convert to Langleys/day
        r = r/0.485d0

* account for cloud cover (assume 0.0)
        rad = r*(-0.21843d0*20.0d0 + 58.94408d0)/100.0d0

* solve for daily maximum (mol/m2/sec)
        pfd2 = (rad*41868.0d0*pi)/(hrs*2.0d0*3600.0d0)
        pfd2 = pfd2*0.5d0*4.255E-6
      ELSE
	pfd2 = 0.0d0
      ENDIF


      RETURN
      END

*----------------------------------------------------------------------*
*                          PHOTON FLUX DENSITY                         *
*                                                                      *
* From Total Downward Shortwave at the Surface (TDSS) - for use with   *
* GCM data.                                                            *
*----------------------------------------------------------------------*
      FUNCTION pfds(tdss,hrs)
*----------------------------------------------------------------------*
      REAL*8 pfds,tdss,hrs,pi
      PARAMETER (pi = 3.1415927)

      IF (hrs.GT.1E-6) THEN
*   solve for daily maximum (mol/m2/sec)
        pfds = (tdss*pi*24.0d0)/(hrs*2.0d0)*4.255E-6
*        pfds = pfds*0.5d0*4.255E-6
      ELSE
	pfds = 0.0
      ENDIF


      RETURN
      END
