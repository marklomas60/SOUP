
      SUBROUTINE GOUDRIAANSLAW(lyr,lai,beamrad,diffrad,
     &     fsunlit,qsunlit, fshade,qshade)
      IMPLICIT NONE
      
      REAL*8 lyr,lai,beamrad,diffrad
      REAL*8 fsunlit,qsunlit,fshade,qshade
      
      REAL*8 soilalbedo,leafalbedo,canopyalbedo
      REAL*8 albedobeam,albedodiff
      REAL*8 m
      REAL*8 kbeam,kdiff,kbeamstar

c     value from Wang
      soilalbedo=0.15d0
c     leaf transmittance and reflectance are equal. Then albedo is the double
      leafalbedo=2.0d0*0.075d0

c     extinction coefficent of the direct beam
      kbeam=0.5d0
c     extinction coefficient of diffuse radiations
      kdiff=0.66d0
c     I found somewhere an approximate relation: kbeam=0.75*kdiff

      m=sqrt(1.0d0-leafalbedo)
      kbeamstar=m*kbeam

      canopyalbedo=2.0d0*kbeam/(kbeam+kdiff)*(1.0d0-m)/(1.0d0+m)
      albedobeam=canopyalbedo+(soilalbedo-canopyalbedo)
     &     *exp(-kbeam*lai)
c     ... hum, I don't know if albedodiff is really different or not... so equal
      albedodiff=albedobeam

      qshade=(1.0d0-albedodiff)*diffrad*kdiff*exp(-kdiff*lyr)
     & +(1.0d0-albedobeam)*beamrad*kbeam*exp(-kbeam*lyr)
     & -(1.0d0-leafalbedo)*beamrad*kbeamstar*exp(-kbeamstar*lyr)
      IF (qshade.LT.0.0d0) qshade=0.0d0

      qsunlit=(1.0d0-leafalbedo)*kbeamstar*beamrad+qshade


      fsunlit=exp(-kbeam*lyr)
      fshade=1-fsunlit
      

      END


      SUBROUTINE BEERSLAW(lyr,beamrad,diffrad,
     &     fsunlit,qsunlit, fshade,qshade)
      IMPLICIT NONE
      REAL*8 lyr,beamrad,diffrad
      REAL*8 fsunlit,qsunlit,fshade,qshade
      REAL*8 kbeam,kdiff

c     extinction coefficent of the direct beam
      kbeam=0.5d0
c     extinction coefficient of diffuse radiations
      kdiff=0.66d0     
      
      fshade=1.0d0
      fsunlit=0

      qsunlit=0.0d0
      qshade=beamrad*kbeam*exp(-kbeam*lyr)+diffrad*kdiff*exp(-kdiff*lyr)

      END


      SUBROUTINE SUNFLECT(lyr,beamrad,diffrad,
     &     fsunlit,qsunlit, fshade,qshade)
      IMPLICIT NONE
      REAL*8 lyr,beamrad,diffrad
      REAL*8 fsunlit,qsunlit,fshade,qshade
      REAL*8 kbeam,kdiff

c     extinction coefficent of the direct beam
      kbeam=0.5d0
c     extinction coefficient of diffuse radiations
      kdiff=0.66d0     
      
      fsunlit=0.5d0
      fshade=0.5d0

      qsunlit=beamrad*kbeam*exp(-kbeam*lyr)+
     &     diffrad*kdiff*exp(-kdiff*lyr)
      qshade=qsunlit


      END
