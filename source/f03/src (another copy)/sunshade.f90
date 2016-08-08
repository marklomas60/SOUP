module sunshade

use real_precision

contains

subroutine SET_GOUD_PARAMS()
real(dp) :: soilalbedo,leafalbedo,kbeam,kdiff,m,kbeamstar,canopyalbedo
common /GOUD_PARAMS/soilalbedo,leafalbedo,kbeam,kdiff,m,kbeamstar,canopyalbedo

soilalbedo = 0.15
leafalbedo = 0.15
kbeam = 0.5
kdiff = 0.66
m = sqrt(1.0-leafalbedo)
kbeamstar = m*kbeam
canopyalbedo=2.0*kbeam/(kbeam+kdiff)*(1.0-m)/(1.0+m)

end subroutine SET_GOUD_PARAMS


subroutine GOUDRIAANSLAW2(lyr,lai,beamrad,diffrad,fsunlit,qsunlit,fshade,qshade, &
 soilalbedo,leafalbedo,kbeam,kdiff,m,kbeamstar,canopyalbedo,albedobeam,albedodiff)
implicit none

real(dp) :: lyr,lai,beamrad,diffrad,fsunlit,qsunlit,fshade,qshade
real(dp) :: soilalbedo,leafalbedo,kbeam,kdiff,m,kbeamstar,canopyalbedo,albedobeam,albedodiff

! soilalbedo: value from Wang
! leafalbedo: leaf transmittance and reflectance are equal. Then albedo is the double
! kbeam: extinction coefficent of the direct beam
! kdiff: extinction coefficient of diffuse radiations
! I found somewhere an approximate relation: kbeam=0.75*kdiff
!... hum, I don't know if albedodiff is REAL*8ly different or not... so equal

qshade=(1.0-albedodiff)*diffrad*kdiff*exp(-kdiff*lyr) &
 +(1.0-albedobeam)*beamrad*kbeam*exp(-kbeam*lyr) &
 -(1.0-leafalbedo)*beamrad*kbeamstar*exp(-kbeamstar*lyr)
if (qshade<0.0) qshade=0.0

qsunlit=(1.0-leafalbedo)*kbeamstar*beamrad+qshade

fsunlit=exp(-kbeam*lyr)
fshade=1-fsunlit

return
end subroutine GOUDRIAANSLAW2


subroutine GOUDRIAANSLAW(lyr,lai,beamrad,diffrad,fsunlit,qsunlit,fshade,qshade)
implicit none

real(dp) :: lyr,lai,beamrad,diffrad
real(dp) :: fsunlit,qsunlit,fshade,qshade

real(dp) :: soilalbedo,leafalbedo,canopyalbedo
real(dp) :: albedobeam,albedodiff
real(dp) :: m
real(dp) :: kbeam,kdiff,kbeamstar

!value from Wang
soilalbedo=0.15
!leaf transmittance and reflectance are equal. Then albedo is the double
leafalbedo=2.0*0.075

!extinction coefficent of the direct beam
kbeam=0.5
!extinction coefficient of diffuse radiations
kdiff=0.66
!I found somewhere an approximate relation: kbeam=0.75*kdiff

m=sqrt(1.0-leafalbedo)
kbeamstar=m*kbeam

canopyalbedo=2.0*kbeam/(kbeam+kdiff)*(1.0-m)/(1.0+m)
albedobeam=canopyalbedo+(soilalbedo-canopyalbedo)*exp(-kbeam*lai)
!... hum, I don't know if albedodiff is REAL*8ly different or not... so equal
albedodiff=albedobeam

qshade=(1.0-albedodiff)*diffrad*kdiff*exp(-kdiff*lyr) &
 +(1.0-albedobeam)*beamrad*kbeam*exp(-kbeam*lyr) &
 -(1.0-leafalbedo)*beamrad*kbeamstar*exp(-kbeamstar*lyr)
if (qshade<0.0) qshade=0.0

qsunlit=(1.0-leafalbedo)*kbeamstar*beamrad+qshade


fsunlit=exp(-kbeam*lyr)
fshade=1-fsunlit


end subroutine GOUDRIAANSLAW


subroutine BEERSLAW(lyr,beamrad,diffrad,fsunlit,qsunlit,fshade,qshade)
implicit none
real(dp) :: lyr,beamrad,diffrad
real(dp) :: fsunlit,qsunlit,fshade,qshade
real(dp) :: kbeam,kdiff

!extinction coefficent of the direct beam
kbeam=0.5
!extinction coefficient of diffuse radiations
kdiff=0.66     

fshade=1.0
fsunlit=0

qsunlit=0.0
qshade=beamrad*kbeam*exp(-kbeam*lyr)+diffrad*kdiff*exp(-kdiff*lyr)

end subroutine BEERSLAW


subroutine SUNFLECT(lyr,beamrad,diffrad,fsunlit,qsunlit,fshade,qshade)
implicit none
real(dp) :: lyr,beamrad,diffrad
real(dp) :: fsunlit,qsunlit,fshade,qshade
real(dp) :: kbeam,kdiff

!extinction coefficent of the direct beam
kbeam=0.5
!extinction coefficient of diffuse radiations
kdiff=0.66     

fsunlit=0.5
fshade=0.5

qsunlit=beamrad*kbeam*exp(-kbeam*lyr)+diffrad*kdiff*exp(-kdiff*lyr)
qshade=qsunlit


end subroutine SUNFLECT

end module sunshade
