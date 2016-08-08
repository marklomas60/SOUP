module metdos

use real_precision

implicit none

contains

!**********************************************************************!
real(dp) function DECLINATION_FN(DAY_OF_YEAR,DAYS_IN_YEAR,PI)
!**********************************************************************!
real(dp) :: MAXDECLINATION, PI
integer ::DAYS_IN_YEAR, DAY_OF_YEAR
	MAXDECLINATION = 23.45
	DECLINATION_FN =MAXDECLINATION*PI/180.0*SIN(PI/180.0*(360*real(DAY_OF_YEAR)/ &
 real(DAYS_IN_YEAR)-80))
!//! Equation: 1
!//! Description: Approximates declination angle for each day: RESULT \
!//!  is in RADIANS. Assume sine wave, with dec=23.45 on 21 June. \\ 
!//!  INPUTS: Day-number of the year; Number of days in the year;\
!//!  mathematical constant, PI.          
!//! Bibliography: No Specific reference to this function.
!
!//! E:  Decl = D_{max} * \pi/180 * SIN(\pi/180*(360*DoY\over{DinY}-80))
!   
!//! Variable: Decl 
!//! Description: returned value of the suns's declination
!//! Units: Radians
!//! SET
!//! Type: REAL
!
!//! Variable: D_{max}
!//! Description: maximum value of the suns's declination
!//! Value: 23.45
!//! Units: Degrees
!//! Type: REAL
!    
!//! Variable: \pi 
!//! Description: Mathematical constant, pi 3.1415...        
!//! Units: Radian
!//! Type: REAL
!
!//! Variable: DoY
!//! Description: Day of the year,(1st jan = 1, 31 dec=365 or 366) 
!//! Units: none          
!//! Type: integer ::
!
!//! Variable: DinY
!//! Description: Total number of Days in year,(leap year=366)
!//! Units: none          
!//! Type: integer ::
!
!**********************************************************************!
end function DECLINATION_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function RADIAN_TO_DEGREE_FN()
!**********************************************************************!
!**   360 degrees is 2*PI radians:

	RADIAN_TO_DEGREE_FN = 180.0/(2.0*ASIN(1.0))
!
!//! Equation: 2
!//! Description: a value for the number of degrees in 1 radian \\
!//!   INPUTS: none
!//! Bibliography: No Specific reference to this function.
!
!//! E:   R2D = 180.0 \over{(2.0*ASIN(1.0)} 
!   
!//! Variable: R2D  
!//! Description: returned value: degrees in 1 radian
!//! Units: Degree/radian
!//! SET
!//! Type: REAL
!
!**********************************************************************!
end function RADIAN_TO_DEGREE_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function DEGREE_TO_RADIAN_FN()
!**********************************************************************!
!**   360 degrees is 2*PI radians
!**   1 degree is approx 0.01745 radians:

	DEGREE_TO_RADIAN_FN = 2.0*ASIN(1.0)/180.0
!
!//! Equation: 3
!//! Description: a value for the number of radians in 1 degree  \\
!//!   INPUTS: none
!//! Bibliography: No Specific reference to this function.
!
!//! E:   D2R = 2.0*ASIN(1.0) \over{180.0} 
!   
!//! Variable: D2R  
!//! Description: returned value: radians in 1 degree
!//! Units: Radian/Degree
!//! SET
!//! Type: REAL
!
!**********************************************************************!
end function DEGREE_TO_RADIAN_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function PIVALUE_FN()
!**********************************************************************!
!**   Returns the value PI (3.14157...)              
	PIVALUE_FN = 2.0*ASIN(1.0)          
!
!//! Equation: 4
!//! Description: Sets the value of Pi \\
!//!   INPUTS: none
!//! Bibliography: No Specific reference to this function.
!
!//! E:   Pi = 2.0*ASIN(1.0) 
!   
!//! Variable: Pi   
!//! Description:  mathematical constant number Pi, 3.14159...
!//! Units: none         
!//! SET
!//! Type: REAL
!
!**********************************************************************!
end function PIVALUE_FN
!**********************************************************************!



!**********************************************************************!
real(dp) function DAYLENGTH_FN(LATITUDE_RAD, DECLINATION,PI)
!**********************************************************************!
!**   calculate the day-length, use declination and latitude            
real(dp) :: DECLINATION, LATITUDE_RAD
real(dp) :: PI, COSHS, HS
real(dp) :: ALTANGLE
	ALTANGLE = -0.833*PI/180.0
	COSHS = -((SIN(LATITUDE_RAD)*SIN(DECLINATION)-SIN(ALTANGLE))/(COS(LATITUDE_RAD)*COS(DECLINATION)))     
	HS = ACOS(COSHS)*180/PI                                        
	DAYLENGTH_FN = HS*2/15.0
!
!//! Equation: 5
!//! Description: calculates the daylength in hours. No Account is  \
!//!   taken of altitude, or longitude.  \\
!//!   INPUTS: Latitude; Solar declination; mathematical constant, PI
!//! Bibliography: No Specific reference to this function.
!
!//! E:   daylen=180*2*ACOS( (SIN(Lat)*SIN(Decl)-SIN(altang)) \over  \
!//!       {(COS(Lat)*COS(Decl))} ) /Pi /15.0  \\
!//!       altang = -0.833 \pi/180.0
!   
!//! Variable: daylen
!//! Description:  appproximation of daylength
!//! Units: hours        
!//! SET
!//! Type: REAL
!
!//! Variable: Lat
!//! Description:  Latitude of site (in radians) 
!//! Units: radians      
!//! Type: REAL
!
!//! Variable: Decl
!//! Description:  declination of sun (radians) 
!//! Units: radians      
!//! Type: REAL
!
!//! Variable: altang
!//! Description: altitude angle (radians); calculated value at sea-level.
!//! Value: -0.14538592....
!//! Units: radians      
!//! Type: REAL
!
!//! Variable: Pi
!//! Description:  mathematical constant number Pi, 3.14159...
!//! Units: none      
!//! Type: REAL
!
!**********************************************************************!
end function DAYLENGTH_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function DAWN_FN(DAYLENGTH)
!**********************************************************************!
!**   calculate the time of dawn: use declination and latitude to
!**   get daylength and go from there value is local solar time.
real(dp) :: DAYLENGTH
	DAWN_FN = 12.0-(DAYLENGTH/2.0)
!
!//! Equation: 6
!//! Description: finds the approximate time of dawn - no account is \
!//!   taken of altitude, longitude etc. \\
!//!   INPUTS: Daylength
!//! Bibliography: No Specific reference to this function.
!
!//! E:   T_{dawn} = 12.0 - (daylen/2.0)
!   
!//! Variable: T_{dawn} 
!//! Description:  Approximate time of dawn.
!//! Units: Local Solar time
!//! SET
!//! Type: REAL
!
!//! Variable: daylen   
!//! Description:  Approximate daylength.   
!//! Units: none         
!//! Type: REAL
!
!**********************************************************************!
end function DAWN_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function DUSK_FN(DAYLENGTH)
!**********************************************************************!
!**   calculate the time of dusk: use declination and latitude            
real(dp) :: DAYLENGTH
	DUSK_FN = 12.0+(DAYLENGTH/2.0)                                         
!
!//! Equation: 7
!//! Description: finds the approximate time of dusk - no account is \
!//!   taken of altitude, longitude etc. \\
!//!   INPUTS: Daylength
!//! Bibliography: No Specific reference to this function.
!
!//! E:   T_{dusk} = 12.0 - (daylen/2.0)
!   
!//! Variable: T_{dusk} 
!//! Description:  Approximate time of dusk.
!//! Units: Local Solar Time
!//! SET
!//! Type: REAL
!
!//! Variable: daylen   
!//! Description:  Approximate daylength.   
!//! Units: none         
!//! Type: REAL
!
!**********************************************************************!
end function DUSK_FN
!**********************************************************************!





!**********************************************************************!
integer function DAY_OF_YEAR_FN(DD,MM,DAYSINYEAR)
!**********************************************************************!
integer ::DAYSINYEAR, LEAPSHIFT
integer ::DD, MM, JUL
!**   returns the day of the year given the month day and year.
!**   checks for leap year using DAYS_IN_YEAR
!**   
!**   These are NOT JULIAN dates (although people often refer mistakenly  
!**   as such).                                                           
	LEAPSHIFT = 0
	if(DAYSINYEAR==366) LEAPSHIFT=1
	if(MM==1) then                                                    
	  JUL = DD                                                        
	else if(MM==2) then    
	  JUL = 31+DD                                                     
	else if(MM==3) then                                               
	  JUL = 31+28+DD + LEAPSHIFT
	else if(MM==4) then                                                
	  JUL = 31+28+31+DD + LEAPSHIFT
	else if(MM==5) then
	  JUL = 31+28+31+30+DD + LEAPSHIFT
	else if(MM==6) then 
	  JUL = 31+28+31+30+31+DD + LEAPSHIFT
	else if(MM==7) then     
	  JUL = 31+28+31+30+31+30+DD + LEAPSHIFT
	else if(MM==8) then 
	  JUL = 31+28+31+30+31+30+31+DD + LEAPSHIFT
	else if(MM==9) then 
	  JUL = 31+28+31+30+31+30+31+31+DD + LEAPSHIFT
	else if(MM==10) then 
	  JUL = 31+28+31+30+31+30+31+31+30+DD + LEAPSHIFT
	else if(MM==11) then 
	  JUL = 31+28+31+30+31+30+31+31+30+31+DD + LEAPSHIFT
	else if(MM==12) then
	  JUL = 31+28+31+30+31+30+31+31+30+31+30+DD + LEAPSHIFT
	endif 
	DAY_OF_YEAR_FN = JUL
!
!//! Equation: 8
!//! Description: finds the day of the year, starting with Jan 1 =1 \\
!//!   INPUTS: Day of Month; Month of year (numeric); Days in Year
!//! Bibliography: No Specific reference to this function.
!
!//! E:   DoY = \sum_{i=0}^{mon-1} Dm_{i} + dd
!   
!//! Variable: DoY      
!//! Description:   Day of the year 
!//! Units: none         
!//! SET
!//! Type: integer ::
!
!//! Variable: mon      
!//! Description:  month number 
!//! Units: none         
!//! Type: integer ::
!
!//! Variable: Dm_{i}
!//! Description:  Days in month (i) 
!//! Units: none         
!//! Type: integer ::
!
!//! Variable: dd       
!//! Description:   day of the current month
!//! Units: none         
!//! Type: integer ::
!
!**********************************************************************!
end function DAY_OF_YEAR_FN
!**********************************************************************!





!**********************************************************************!
integer function DAYS_IN_YEAR_FN(YYYY)
!**********************************************************************!
!**   checks the year and tests if it is a leap-year. Returns the
!**   number of days in the year (ie. 365 or 366). 
integer ::YYYY, MODVAL
	MODVAL=100
	if( MOD(YYYY,MODVAL)==0) then
!**       century - must be divisible by 400
	  MODVAL=400
	  if( MOD(YYYY,MODVAL)==0) then
	    DAYS_IN_YEAR_FN=366
	  else
	    DAYS_IN_YEAR_FN=365
	  endif
	else
!**     not a century - is it divisble by 4
	  MODVAL=4
	  if( MOD(YYYY,MODVAL)==0) then
	    DAYS_IN_YEAR_FN=366
	  else
	    DAYS_IN_YEAR_FN=365
	  endif
	endif
return
!
!//! Equation: 9
!//! Description:  calculates number of days in any year  \\
!//!   INPUTS: Year
!//! Bibliography:  After: L. Halsall, Stats & Computing, Forest \
!//!   Research, Forestry Commission, Alice Holt, Farnham.
!
!//! Case: MOD(year, 100)=0 and MOD(year,400)=0
!//! E:   DinY = 366
!//! Case: MOD(year,100)<>0 and MOD(year,4)=0
!//! E:   DinY = 366
!//! Case: All other conditions
!//! E:   DinY = 365
!   
!//! Variable: DinY     
!//! Description: Number of days in a year:    \\
!//!   If it is a century and divisible by 400, then it is a leap year \\
!//!   If its NOT a century but divisible by 4, then it is a leap year
!//! Units: none         
!//! SET
!//! Type: integer ::
!
!//! Variable: year     
!//! Description:   full year number (eg 1997, 2002 etc).
!//! Units: none         
!//! Type: integer ::
!
!**********************************************************************!
end function DAYS_IN_YEAR_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function DECLINATION2_FN(DAY_OF_YEAR)
!**********************************************************************!
real(dp) :: MAXDECLINATION
integer ::DAY_OF_YEAR
	MAXDECLINATION = 0.39785
	DECLINATION2_FN = MAXDECLINATION*SIN(4.868961+0.017203*DAY_OF_YEAR+0.033446*SIN(6.224111+0.017202*DAY_OF_YEAR))
!
!//! Equation: 15
!//! Description: Approximates declination angle for each day: RESULT \
!//!  is in RADIANS.                                       \\
!//!   INPUTS: Day-number of year
!//! Bibliography:  S. Evans: SWELTER 
!
!//! E:   Decl = D_{max} * sin(4.868961+0.017203 DoY + 0.033446      \
!//!             sin(6.22411 + 0.017202 DoY))
!   
!//! Variable: Decl 
!//! Description: returned value of the suns's declination
!//! Units: Radians
!//! SET
!//! Type: REAL
!
!//! Variable: D_{max}
!//! Description: maximum value of the suns's declination
!//! Value: 0.39785
!//! Units: Radians
!//! Type: REAL
!    
!//! Variable: DoY
!//! Description: Day of the year,(1st jan = 1, 31 dec=365 or 366) 
!//! Units: none          
!//! Type: integer ::
!
!**********************************************************************!
end function DECLINATION2_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function DAYLENGTH2_FN(LATITUDE_RAD, DAY_OF_YEAR,PI)
!**********************************************************************!
!**   calculate the day-length, use declination and latitude            
real(dp) :: LATITUDE_RAD
real(dp) :: PI
real(dp) :: AMPLITUDE
integer ::DAY_OF_YEAR
	AMPLITUDE = exp((7.42+0.045*LATITUDE_RAD*180./PI))/3600.0
	DAYLENGTH2_FN = 12.0 + AMPLITUDE*SIN((DAY_OF_YEAR-79)*0.017121)
!
!//! Equation: 16
!//! Description: calculates the daylength in hours from latitude \
!//! and day-number. no Account is taken of altitude, or longitude. \\ 
!//!   INPUTS: Latitude; Day-number of year, Mathematical constant, PI
!//! Bibliography: S. Evans: SWELTER 
!
!//! E:   daylen= 12 + Amp sin((DoY-79) 0.017121)   \\ 
!//!        Amp={exp(7.42+0.045 lat 180/\pi )} \over{3600}
!   
!//! Variable: daylen
!//! Description:  appproximation of daylength. based on solar time
!//! Units: hours        
!//! SET
!//! Type: REAL
!
!//! Variable: lat
!//! Description:  Latitude of site (in radians) 
!//! Units: radians      
!//! Type: REAL
!
!//! Variable: Amp
!//! Description: amplitude of sine curve at relevant latitude. 
!//! Units: none      
!//! Type: REAL
!
!//! Variable: Pi
!//! Description:  mathematical constant number Pi, 3.14159...
!//! Units: none      
!//! Type: REAL
!
!//! Variable: DoY 
!//! Description: Day of the year,(1st jan = 1, 31 dec=365 or 366) 
!//! Units: none          
!//! Type: integer ::
!
!**********************************************************************!
end function DAYLENGTH2_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function DAWN_HOURANGLE_FN(DECLINATION,LATITUDE_RAD)
!**********************************************************************!
!**   find the sun-rise/sun-set solar hour-angle: ie at DAWN and DUSK
real(dp) :: DECLINATION, LATITUDE_RAD, Temp
  Temp = -TAN(LATITUDE_RAD)*TAN(DECLINATION) 
  if(Temp<-1) Temp = -1.0
  if(Temp>1) Temp =1.0
	DAWN_HOURANGLE_FN = ACOS(Temp)
!
!//! Equation: 17
!//! Description: calculates the sunrise/sunset hour angle for the \
!//!  day (declination) and latitude. Answer is in Radians \\
!//!   INPUTS: Solar declination; Latitude
!//! Bibliography:  S. Evans: SWELTER 
!
!//! E:   hs = acos(-tan(Lat)*tan(Decl))
!   
!//! Variable: hs   
!//! Description: returned value of the solar sunrise/set hour angle \
!//!  uncorrected for longitude.
!//! Units: Radians
!//! SET
!//! Type: REAL
!
!//! Variable: Decl
!//! Description:  declination of sun (radians) 
!//! Units: radians      
!//! Type: REAL
!
!//! Variable: Lat
!//! Description:  Latitude of site (in radians) 
!//! Units: radians      
!//! Type: REAL
!
!**********************************************************************!
end function DAWN_HOURANGLE_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function RAD_ET_DAY2_FN(LATITUDE_RAD,DECLINATION,DAWN_HOURANGLE, SOL_ELLIPSE, SOL_CONST, PI)
!**********************************************************************!
!**   caclulates the extra-terrestrial daily radiation 
real(dp) :: DECLINATION, LATITUDE_RAD
real(dp) :: DAWN_HOURANGLE, SOL_CONST, PI
real(dp) ::  SOL_ELLIPSE

RAD_ET_DAY2_FN = SOL_CONST*SOL_ELLIPSE*(DAWN_HOURANGLE*SIN(LATITUDE_RAD)*SIN(DECLINATION) + &
 COS(LATITUDE_RAD)*COS(DECLINATION)*SIN(DAWN_HOURANGLE))/ PI*86400

!//! Equation: 18
!//! Description: caclulates the daily extra-terrestrial radiation  \
!//!    result in W/m2/day         \\
!//!   INPUTS: Latitude; Solar Declination; Hour-angle of sun-rise/set;\
!//!    Elliptical scalar for suns' orbit; Solar constant, Mathematical\
!//!    constant, Pi
!//! Bibliography:  S. Evans: SWELTER 
!
!//! E:   RAD_{et}= S_c E (hs sin(Lat) sin(Decl) +        \
!//!       cos(lat) cos(decl) sin(hs)) \over{\pi} 86400
!   
!//! Variable: RAD_{et}
!//! Description: daily value of solar radiation outside the atmosphere
!//! Units: W/m2/day
!//! SET
!//! Type: REAL
!
!//! Variable: Decl
!//! Description:  declination of sun (radians) 
!//! Units: radians      
!//! Type: REAL
!
!//! Variable: Lat
!//! Description:  Latitude of site (in radians) 
!//! Units: radians      
!//! Type: REAL
!
!//! Variable: S_c
!//! Description:  Solar constant. Nb there may be slight variations \
!//!    on this depending on the source (approx 1367)
!//! Units: W/m2
!//! Type: REAL
!
!//! Variable: Pi   
!//! Description:  mathematical constant Pi: 3.14159...
!//! Units: none          
!//! Type: real(dp)
!
!//! Variable: hs   
!//! Description: value of the solar sunrise/set hour angle.
!//! Units: Radians
!//! Type: REAL
!
!//! Variable: E
!//! Description: Solar elliptical orbit scalar
!//! Units: none          
!//! Type:  real(dp)
!
!**********************************************************************!
end function RAD_ET_DAY2_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function RAD_ET_DAY_FN(LATITUDE_RAD,DECLINATION,DAYLENGTH,SOL_ELLIPSE,SOL_CONST,PI)
!**********************************************************************!
!**   caclulates the extra-terrestiral daily radiation 
real(dp) :: DECLINATION, LATITUDE_RAD
real(dp) :: SOL_CONST, PI
real(dp) :: DAYLENGTH, NOON
real(dp) :: SOL_ELEV_NOON, SOL_ELLIPSE, SOL_NOON

!  since the routine sol_elev has local time as an input; convert 
!  noon solar time to gmt
NOON = 12.0     ! LST
SOL_ELEV_NOON = SOL_ELEV_FN(LATITUDE_RAD,DECLINATION,NOON,PI)
SOL_NOON = SOL_ELLIPSE*SOL_CONST*SIN(SOL_ELEV_NOON)
RAD_ET_DAY_FN = SOL_NOON*2*DAYLENGTH/(PI*SIN(PI/2.0))*60*60 !  /1000/1000 
!
!//! Equation: 19
!//! Description: caclulates the daily extra-terrestrial radiation  \
!//!    result in W/m2/day                        \\
!//!   INPUTS: Latitude; Longitude; Solar declination; Daylength; \
!//!     Solar elliptical orbit scalar; Day-number of the year;\
!//!     Solar constant; Mathematical constant, pi.
!//! Bibliography: none/FOREST-GROWTH
!
!//! E:   RAD_{et}= S_{12} 2 Day_{l} / (\pi sin(\pi/2)                   \\
!//!       S_{12} = E S_c sin(Sol_el)                                 \\
!//!       Sol_{el} = \pi/2 - Lat + Decl                              \\
!//!       E =  1.033-0.066*\sqrt{(1-((DoY-183)**2)/(183.**2) ) }
!   
!//! Variable: RAD_{et}
!//! Description: daily value of solar radiation outside the atmosphere
!//! Units: W/m2/s
!//! SET
!//! Type: REAL
!
!//! Variable: E   
!//! Description:  Elliptical scalar for the orbit of the sun 
!//! Units: none         
!//! Type: REAL
!
!//! Variable: Decl
!//! Description:  declination of sun (radians) 
!//! Units: radians      
!//! Type: REAL
!
!//! Variable: Lat
!//! Description:  Latitude of site (in radians) 
!//! Units: radians      
!//! Type: REAL
!
!//! Variable: S_c
!//! Description:  Solar constant. Nb there may be slight variations \
!//!    on this depending on the source  (eg 1367.0)
!//! Units: W/m2         
!//! Type: REAL
!
!//! Variable: Pi   
!//! Description:  mathematical constant Pi: 3.14159...
!//! Units: none          
!//! Type: real(dp)
!
!//! Variable: Day_l
!//! Description:  day length
!//! Units: hours         
!//! Type: real(dp)
!
!//! Variable: DoY 
!//! Description: Day of the year,(1st jan = 1, 31 dec=365 or 366) 
!//! Units: none          
!//! Type: integer ::
!
!**********************************************************************!
end function RAD_ET_DAY_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function RAD_ET_DAY3_FN(LATITUDE_RAD,DECLINATION,SOL_ELLIPSE,DAWN,DUSK,SOL_CONST,PI)
!**********************************************************************!
!**   caclulates the extra-terrestrial daily radiation 
real(dp) :: DECLINATION, LATITUDE_RAD
real(dp) :: SOL_CONST, PI, SOL_ELLIPSE
real(dp) :: SOLAR_TIME
real(dp) :: SOL_NOW, SIN_SOL_ELEV, DAWN, DUSK
integer :: i

RAD_ET_DAY3_FN = 0.0
do i=0,12*24
!DO SOLAR_TIME = 0,24, 1/12.0
  solar_time = real(i)*1/12.0
!**     step throughout the day in 5 min intervals
	if((SOLAR_TIME>=DAWN).and.(SOLAR_TIME<DUSK)) then
	  SIN_SOL_ELEV = SIN( SOL_ELEV_FN(LATITUDE_RAD,DECLINATION,SOLAR_TIME,PI))
	  SOL_NOW = MAX(0.0,SOL_ELLIPSE*SOL_CONST*SIN_SOL_ELEV)*60*5
	  RAD_ET_DAY3_FN = RAD_ET_DAY3_FN + SOL_NOW
	endif         
enddo
!
!//! Equation: 20
!//! Description: caclulates the daily extra-terrestrial radiation  \
!//!   by 5 minute interation using solar elevation result in W/m2/day. \
!//!   This functions call a solar elevation fuunction internally.  \\
!//!   INPUTS: Latitude; Longitude; Solar Declination; scalar of\
!//!   sun's elliptical orbit; Day-number of the year; Solar constant;\
!//!   Mathematical constant, pi
!//! Bibliography: C.J.T Spitters et al. (1986); Agric & Forest Met \
!//!   38:217-229
!
!//! E:  RAD_{et}= \sum_{t=0}^{24*12} max(0, S_c E sin(Sol_{el}) 300) \\
!//!     if(t>=dawn and t<= dusk)
!   
!//! Variable: RAD_{et}
!//! Description: daily value of solar radiation outside the atmosphere
!//! Units: W/m2/s
!//! SET
!//! Type: REAL
!
!//! Variable: Sol_{el}
!//! Description:  Solar Elevation; returned through SOL_ELEV_FN
!//! Units: radians      
!//! Type: REAL
!
!//! Variable: Decl
!//! Description:  declination of sun (radians) 
!//! Units: radians      
!//! Type: REAL
!
!//! Variable:  E 
!//! Description:   Solar elliptical scalar     
!//! Units:  none        
!//! Type: REAL
!
!//! Variable: Lat
!//! Description:  Latitude of site (in radians) 
!//! Units: radians      
!//! Type: REAL
!
!//! Variable: Dawn
!//! Description:  Time of Sunrise (Local Solar time)
!//! Units: radians      
!//! Type: REAL
!
!//! Variable: Dusk
!//! Description:  Time of Sunset (Local Solar time)
!//! Units: radians      
!//! Type: REAL
!
!//! Variable: S_c
!//! Description:  Solar constant. Nb there may be slight variations \
!//!    on this depending on the source  (1367.0)
!//! Units: W/m2         
!//! Type: REAL
!
!//! Variable: Pi   
!//! Description:  mathematical constant Pi: 3.14159...
!//! Units: none          
!//! Type: REAL
!
!//! Variable: DoY 
!//! Description: Day of the year,(1st jan = 1, 31 dec=365 or 366) 
!//! Units: none          
!//! Type: integer ::
!
!**********************************************************************!
end function RAD_ET_DAY3_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function VP_SAT_FN(TEMP_C)
!**********************************************************************!
!**   Returns the saturated air pressure for a given temperature
real(dp) :: TEMP_C
!
VP_SAT_FN = 6.1078*exp(17.269*TEMP_C/(TEMP_C+237.3))
!
!//! Equation: 21
!//! Description: Returns saturated vapour pressure for a temperature\\
!//!   INPUTS: Temperature
!//! Bibliography: S Evans: SWELTER, Groff-Gratch equation
!
!//!  E:    VP_{sat} = 6.1078 exp{( {17.269 T}\over{(T+237.3)} )}
!
!//! Variable: VP_{sat}
!//! Description: saturated vapour pressure at a given temperature
!//! Units: mbar
!//! SET
!//! Type: REAL
!
!//! Variable: T
!//! Description: Air temperature 
!//! Units: Degrees  C 
!//! Type: REAL
!
!**********************************************************************!
end function VP_SAT_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function SOL_ELEV_FN(LATITUDE_RAD,DECLINATION,SOLAR_TIME,PI)
!**********************************************************************!
!      function returns the elevation (radians), deviation of the sun from
!      the horizontal to the vertical. Also referred to as the solar 
!      altitude angle.
										
real(dp) :: LATITUDE_RAD
real(dp) :: DECLINATION, H, PI
real(dp) :: H1, SOLAR_TIME
									
H1 = 15.*(SOLAR_TIME-12.0)       ! degrees   
H = H1/180*PI  ! radians                                      
SOL_ELEV_FN = ASIN(COS(LATITUDE_RAD)*COS(DECLINATION)*COS(H) + SIN(LATITUDE_RAD)*SIN(DECLINATION))
!//! Equation: 22
!//! Description: Returns solar elevation angle at time of the day \
!//!    and day of the year.                                      \\
!//!  INPUTS:  Latitude; Longitude; Declination; \\
!//!   Local Solar time; mathematical constant, pi
!//! Bibliography: C GronBeck, 1999; Sun Angle (web pages)
!
!//!  E:   SOL_{el} = asin(sin(Lat) sin(Decl) + cos(Lat) cos(Decl) \
!//!     cos(H)                                                       \\
!//!     H = 15(Lst-12) \pi /180.0 
!
!//! Variable: SOL_{el}
!//! Description:  Solar elevation of the sun
!//! Units: Radians     
!//! SET
!//! Type: REAL
!
!//! Variable: Lat
!//! Description: Latitude of site 
!//! Units:  Radians
!//! Type: REAL
!
!//! Variable: Decl
!//! Description:  declination of sun (radians) 
!//! Units: radians      
!//! Type: REAL
!
!//! Variable: Lst
!//! Description: Local solar time
!//! Units: (hour) 
!//! Type: REAL
!
!//! Variable: \pi
!//! Description: Mathematical constant, pi
!//! Units: (radians)
!//! Type: real(dp)
!
!**********************************************************************!
end function SOL_ELEV_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function AZIMUTH_FN(LATITUDE_RAD,DECLINATION,SOL_ELEV)
!**********************************************************************!
!     function returns the azimuth angle (in radians): Angular distance
!     between SOUTH and the projection of the line of sight to the
!     sun on the ground.
										
real(dp) :: LATITUDE_RAD
real(dp) :: DECLINATION, SOL_ELEV, SIGN, HOURANG

! check the sign of the hour-angle. Firstly convert back to hour-angle
!    using latitude, declination and elevation.
HOURANG = ACOS((SIN(SOL_ELEV)-SIN(LATITUDE_RAD)*SIN(DECLINATION))/(COS(LATITUDE_RAD)*COS(DECLINATION)))
SIGN = 1.0
if(HOURANG<0) SIGN=-1.0
									
AZIMUTH_FN = SIGN*ACOS((SIN(SOL_ELEV)*SIN(LATITUDE_RAD)-SIN(DECLINATION))/(COS(SOL_ELEV)*COS(LATITUDE_RAD) ) )
!//! Equation: 23
!//! Description: Returns azimuth angle at any time of the day and  \
!//!     day of the year. \\
!//!  INPUTS: Latitude; solar declination; Solar elevation; \
!//!  Mathematical constant, pi
!//! Bibliography: C GronBeck, 1999; Sun Angle (web pages)
!
!//!  E:   Az = acos( {sin(Sol_{el}) sin(Lat) - sin(Decl)}   \
!//!     /over{cos(Sol_{el}) cos(Decl)}                      \\
!//!     and: the sign of Az is the same as the sign of the hour-angle.
!
!//! Variable: Az    
!//! Description:  Azimuth angle of the sun. 
!//! Units: Radians     
!//! SET
!//! Type: REAL
!
!//! Variable: Lat
!//! Description: Latitude of site 
!//! Units:  Radians
!//! Type: REAL
!
!//! Variable:  Decl
!//! Description: solar declination
!//! Units:  Radians
!//! Type: REAL
!
!//! Variable: Sol_{el}       
!//! Description: Solar elevation 
!//! Units:  Radians
!//! Type: real(dp)
!
!//! Variable: \pi
!//! Description: Mathematical constant, pi
!//! Units: (radians)
!//! Type: real(dp)
!
!**********************************************************************!
end function AZIMUTH_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function RAD_TERRA_DAY_FN(RAD_ET_DAY,LATITUDE_RAD,CLOUD_COVER,BETA)
!**********************************************************************!
!     converts extra-terrestiral radiation to surface radiation.
real(dp) :: RAD_ET_DAY, LATITUDE_RAD
real(dp) :: CLOUD_COVER
real(dp) :: ALPHA, BETA , SIGMA
!
!**     Angstrom tubidity factor related to aerosol size & properties
ALPHA = 18.0-64.884*(1-1.3614*COS(LATITUDE_RAD))
SIGMA = 0.03259
Alpha = alpha*4.1842*100*100  ! convert from cal/cm to w/m2
!      another fudge to try and transmit less light on heavy cloud days
!      if(cloud_cover.gt.1.0) alpha = alpha*2
!       if(cloud_cover.lt.0.15) then
!         cloud_cover = cloud_cover/1.4
!         alpha = 0.0
!       endif
 sigma = .02*LOG(max(cloud_cover,.001))+Sigma
RAD_TERRA_DAY_FN = RAD_ET_DAY*(BETA-SIGMA*CLOUD_COVER*10)-ALPHA
!
!//! Equation: 24
!//! Description: Converts extra-terrestrial radiation to terrestrial \
!//!   radiation  uncorrected for site altitude.                      \\ 
!//!   INPUTS: Extra-terrestrial radiation; Latitude; cloud-cover; \
!//!     Angstrom (beta) turbidity
!//! Bibliography: S. Evans:SWELTER; Nikolov & Zeller (1992)
!
!//! E:   RAD_{terra} = RAD_{et} (\beta -\sigma 10 Cld) - \alpha   \\
!//!        \alpha = 18.0-64.884*(1-1.3614*cos(Lat))     \\
!//!        \sigma = 0.03259   
!
!//! Variable: RAD_{terra}
!//! Description: Terrestrial radiation on a horizonatal plane \
!//!  for the given latitude and day.
!//! Units: W/m2/d
!//! Type: REAL
!//! SET
!
!//! Variable: RAD_{et}
!//! Description: Extra-Terrestrial radiation above the atmosphere  \
!//!  for the given latitude and day.
!//! Units: W/m2/d
!//! Type: REAL
!
!//! Variable: \alpha  
!//! Description: Angstrom tubidity factor related to aerosol size & \
!//!  properties
!//! Units: none  
!//! Type: REAL
!
!//! Variable: \beta
!//! Description: Angstrom tubidity factor: max clear sky transmittance
!//! Units: none  
!//! Type: REAL
!
!//! Variable: \sigma  
!//! Description: Angstrom tubidity factor: Absorption by cloud cover
!//! Units: none  
!//! Type: REAL
!
!//! Variable: Lat
!//! Description:  Latitude of site (in radians) 
!//! Units: radians      
!//! Type: REAL
!
!//! Variable: Cld     
!//! Description: Cloudiness of the day 
!//! Units: proportion 
!//! Type: REAL
!
!**********************************************************************!
end function RAD_TERRA_DAY_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function DIFFPROP_DAY_FN(TRANS_PROP_DAY)
!**********************************************************************!
!**   determines the proportion of light  that is diffuse on a daily 
!**   basis
real(dp) :: TRANS_PROP_DAY
!
if(TRANS_PROP_DAY<0.07) then
	DIFFPROP_DAY_FN = 1.0
else if(TRANS_PROP_DAY<0.35) then
	DIFFPROP_DAY_FN = 1.0-2.3*(TRANS_PROP_DAY-0.07)**2
else if(TRANS_PROP_DAY<0.75) then
	DIFFPROP_DAY_FN = 1.33-1.46*TRANS_PROP_DAY
else
	DIFFPROP_DAY_FN = 0.23
endif
!
!//! Equation: 25
!//! Description: The proportion of light reaching the ground that \
!//!   is made up of diffuse component                      \\
!//!   INPUTS: proportion of extra-terrestrial radiation to terrestrial \
!//!    radiation transmitted
!//! Bibliography: Spitters Et al (1986), Agric & Forest Meteo  \
!//!  38:217-229; de Jong (1980).
!  
!//! Case: Tran < 0.07
!//!  E: Diff_{prop} = 1.0
!//! Case: 0.07<= Tran < 0.35
!//!  E: Diff_{prop} = 1-0.23(Tran - 0.07)^{2}
!//! Case: 0.35<= Tran < 0.75
!//!  E: Diff_{prop} = 1.33 - 1.46 Tran
!//! Case: 0.75 <= Tran
!//!  E: Diff_{prop} = 0.23
!
!//! Variable: Diff_{prop}
!//! Description: Proportion of daily terrestrial radiation that is \
!//!   diffuse                                              \\
!//!   INPUTS: Extra-terrestrial daily radiation; Terrestrial daily \
!//!   radiation
!//! Units: none (Ratio)
!//! Type: REAL
!//! SET
!
!//! Variable: Tran       
!//! Description: Proportion of Terrestrial radiation to    \
!//!   Extra-terrestrial radiation
!//! Units: none (Ratio)
!//! Type: REAL
!
!**********************************************************************!
end function DIFFPROP_DAY_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function SOLAR_ELLIPSE_FN(DAY_OF_YEAR,DAYS_IN_YEAR,PI)
!**   calculates the scalar of the elliptical orbit of the sun around 
!**   the earth. About +/- 3%; Max in mud-summer
real(dp) :: GAMMA, PI
integer ::DAY_OF_YEAR, DAYS_IN_YEAR
										
GAMMA = 2*PI*(DAY_OF_YEAR-1)/real(DAYS_IN_YEAR)
SOLAR_ELLIPSE_FN = 1.00011+0.034211*COS(GAMMA)+0.00128*SIN(GAMMA)+0.000719*COS(2*GAMMA)+0.000077*SIN(2*GAMMA)
!
!//! Equation: 26
!//! Description: Returns the scalar of the elliptical orbit of the \
!//!   sun around the earth (approx +/- 3%; peak mid-summmer)    \\
!//!   INPUTS: Day-number of the year; Number of days in the year;\\
!//!    mathematical constant, pi
!//! Bibliography: S. Evans: SWELTER
!
!//! E:  Ellipse = 1.00011+0.034211 cos(\Gamma) +              \
!//!     0.00128 sin(\Gamma)+0.000719 cos(2 \Gamma) +         \
!//!     0.000077*sin(2 \Gamma)                           \\
!//!     \Gamma = {2 \pi (DoY-1)} \over{DinY}
!
!//! Variable: Ellipse
!//! Description: Scalar of solar ellipse of the sun around the earth 
!//! Units: none
!//! Type: REAL
!//! SET
!
!//! Variable: \Gamma
!//! Description: Intermediate variable
!//! Units: none
!//! Type: REAL
!
!//! Variable: DoY
!//! Description: Day of the year,(1st jan = 1, 31 dec=365 or 366) 
!//! Units: none          
!//! Type: integer ::
!
!//! Variable: DinY
!//! Description: Total number of Days in year,(leap year=366)
!//! Units: none          
!//! Type: integer ::
!
!//! Variable; \pi
!//! Description: the constant, pi, 3.12415...
!//! Units: (radian)
!//! Type: REAL
!
!**********************************************************************!
end function SOLAR_ELLIPSE_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function SOLAR_ELLIPSE2_FN(DAY_OF_YEAR)
!**********************************************************************!
!**   calculates the scalar of the elliptical orbit of the sun around 
!**   the earth. About +/- 3%; Max in mud-summer
integer ::DAY_OF_YEAR
										
SOLAR_ELLIPSE2_FN = 1.033-0.066*sqrt(1.0-((DAY_OF_YEAR-183)**2)/(183.**2))  
!
!//! Equation: 27
!//! Description: Returns the scalar of the elliptical orbit of the \
!//!   sun around the earth (approx +/- 3%; peak mid-summmer)    \\
!//!   INPUTS: Day-number of the year
!//! Bibliography: 
!
!//! E:  Ellipse = 1.033-0.066 \sqrt{{(1.0-((DoY-183)**2)} \over    \
!//!          {(183.**2)} )}  
!
!//! Variable: Ellipse
!//! Description: Scalar of solar ellipse of the sun around the earth 
!//! Units: none
!//! Type: REAL
!//! SET
!
!//! Variable: DoY
!//! Description: Day of the year,(1st jan = 1, 31 dec=365 or 366) 
!//! Units: none          
!//! Type: integer ::
!
!**********************************************************************!
end function SOLAR_ELLIPSE2_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function SOLAR_TIME_FN(GMT,LONGITUDE_RAD,DAY_OF_YEAR,PI)
!**********************************************************************!
!     Routine calculates the local solar time from GMT. Applies only 
!     to GMT based time.
real(dp) :: GMT, EQNTIME, CORR1, PI
real(dp) :: LONGITUDE_RAD
integer ::DAY_OF_YEAR

!      correct time by equation of time and the longitudinal offset    
 CORR1 = (279.575+0.986*DAY_OF_YEAR)*PI/180.                         
!      Equation of time:                                                
 EQNTIME = (-104.7*SIN(CORR1) + 596.2*SIN(2*CORR1) +    &               
 4.3*SIN(3*CORR1) - 12.7*SIN(4*CORR1) - 429.3*COS(CORR1) -  &     
 2*COS(2*CORR1) + 19.3*COS(3*CORR1))/3600.                       
!      
 SOLAR_TIME_FN = GMT+EQNTIME+LONGITUDE_RAD*180.0/PI/15.0
!
!//! Equation: 28
!//! Description: Converts from GMT (hours) to local solar time;   \
!//!  ie accounts for longitude and the equation of time.     \\
!//!  INPUTS: Grenwich mean-time; longitude; day-number of the year;\
!//!   mathematical constant, pi
!//! Bibliography: Plaz & Grief (1996); Gronbeck 1999 (Web Pages)
!
!//! E:  Sol_{time} = GMT + EQN_{time}+ {{180 Long } \over{\pi 15}}
!
!//! Variable: Sol_{time}
!//! Description: Local solar time  - corrected solar time (accounts \
!//!  for longitude etc.)
!//!  Units: Hours
!//!  Type: REAL
!//!  SET
! 
!//! Variable: GMT
!//! Description: Time of day  (Grenwich mean time)
!//!  Units: Hours
!//!  Type: REAL
!
!//! Variable: EQN_{time}
!//! Description: Equation of time ; corrects for the slightly   \
!//!   elliptical orbit of the earth.
!//!  Units: Hours
!//!  Type: REAL
!
!//! Variable: Long      
!//! Description: Longitude of Site 
!//!  Units: Radians
!//!  Type: REAL
!
!//! Variable: \pi
!//! Description: Mathematical constant, pi: 3.1415...
!//! Units: (radian)
!//! Type: REAL
!
!**********************************************************************!
end function SOLAR_TIME_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function GMT_FN(SOLAR_TIME,LONGITUDE_RAD,DAY_OF_YEAR,PI)
!**********************************************************************!
!     Routine calculates the GMT from solar time. Applies only to GMT
!     based time.
real(dp) :: SOLAR_TIME, EQNTIME, CORR1, PI
real(dp) :: LONGITUDE_RAD
integer ::DAY_OF_YEAR

!      correct time by equation of time and the longitudinal offset    
 CORR1 = (279.575+0.986*DAY_OF_YEAR)*PI/180.                         
!      Equation of time:                                                
 EQNTIME = (-104.7*SIN(CORR1) + 596.2*SIN(2*CORR1) +   &                
 4.3*SIN(3*CORR1) - 12.7*SIN(4*CORR1) - 429.3*COS(CORR1) - &      
 2*COS(2*CORR1) + 19.3*COS(3*CORR1))/3600.                       
!      
 GMT_FN = SOLAR_TIME -(EQNTIME+LONGITUDE_RAD*180.0/PI/15.0)
!
!//! Equation: 29
!//! Description: Converts from local solar time to GMT;   \
!//!  ie accounts for longitude and the equation of time.     \\
!//!  INPUTS: Local Solar time; longitude; day-number of the year; Pi
!//! Bibliography: Plaz & Grief (1996); Gronbeck 1999 (Web Pages)
!
!//! E:  GMT = Sol_{time} -(EQN_{time} +{{180 Long } \over{\pi 15}})
!
!//! Variable: GMT
!//! Description:  Grenwich Mean Time of day (corrected for longitude \
!//!  etc.)
!//!  Units: Hours
!//!  Type: REAL
!//!  SET
! 
!//! Variable: Sol_{time}
!//! Description: Local solar time  - corrected solar time 
!//!  Units: Hours
!//!  Type: REAL
!
!//! Variable: EQN_{time}
!//! Description: Equation of time ; corrects for the slightly   \
!//!   elliptical orbit of the earth.
!//!  Units: Hours
!//!  Type: REAL
!
!//! Variable: Long      
!//! Description: Longitude of Site 
!//!  Units: Radians
!//!  Type: REAL
!
!//! Variable: \pi
!//! Description: Mathematical constant, pi: 3.1415...
!//! Units: (radian)
!//! Type: REAL
!
!**********************************************************************!
end function GMT_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function TEMP_HOUR_FN(HOUR, TEMPMAX, TEMPMIN,PI)
!**********************************************************************!
!**  This function fits a sine curve over the day, ranging from 
!**  TEMPMIN to TEMPMAX, with the maximum when HOUR=14 (ie 2pm)
real(dp) :: HOUR, TEMPMAX, TEMPMIN
real(dp) :: PI, TMEAN
!
TMEAN = ((TEMPMAX-TEMPMIN)/2.0)+TEMPMIN                       
TEMP_HOUR_FN = SIN(((HOUR-2)/24.0*360-90)*PI/180.0)*((TEMPMAX-TEMPMIN)/2.0) + TMEAN
!
!//! Equation: 30
!//! Description: returns the temperature at the time of day from \
!//!    fitting a sine curve with peaks/troughs at the maximum and \
!//!    minimum of the day. Peak assumed to be at time=14:00. Time is \
!//!    assumed to be absolute ie no corrections for latitude are made.\\
!//!  INPUTS: Hour of the day; Maximum temperature of the day; \
!//!    minimum temperature of the day. 
!//! Bibliography: none
!
!//!   E: TEMP_{hour} = ({{T_{max}-T_{min}} \over{2}} + T_{min}) +    \
!//!          sin( { ({{{360(hour-2)}\over{24}} - 90} \over {180} )}  \
!//!          ({{T_{max}-T_{min}} \over{2}})
!//!
!//!  Variable: TEMP_{hour}
!//!  Description: Temperature of day at time 'hour'
!//!  Units: Degrees celcius  (usually)
!//!  Type: real(dp)
!//!  SET      
!
!//!  Variable: T_{max}
!//!  Description: Maximum temperature of the day
!//!  Units: Degrees celcius  (usually)
!//!  Type: REAL
!
!//!  Variable: T_{min}
!//!  Description: Minumum temperature of the day
!//!  Units: Degrees celcius  (usually)
!//!  Type: REAL
!
!//!  Variable: \pi     
!//!  Description: Mathematical value 3.1415..
!//!  Units: none 
!//!  Type: REAL
!
!**********************************************************************!
end function TEMP_HOUR_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function HUMIDITY_INST_FN(VP_UNSAT_BASE,VP_SAT_NOW)
!**********************************************************************!
!**  Function calculates the humidity at any temperature
!**  Values at temperatures below freezing are doubtful.
real(dp) ::  VP_UNSAT_BASE, VP_SAT_NOW
!
 HUMIDITY_INST_FN = 100*VP_UNSAT_BASE/VP_SAT_NOW                                   
!//! Equation: 31
!//! Description: Gets humidity at any time of the day eg 85%. \
!//!   Assumes that vapour pressure is constant through the day. \
!//!   usually we assume that VP is constant through the day and use\
!//!   the '9am' values are used to get the vp; as temperature \
!//!   changes, so will humidity. Beware if temperature < ZERO. \\
!//!   INPUTS: vapour pressure at base time;   \
!//!    Saturated VP at current temp.
!//! Bibliography: LI-Cor Manual ??
!//!  
!//!  E: HUMID_{t} =  {{100 VP_{unsat,base}} \over{VP_{sat,t}} } 
!
!//! Variable: HUMID_{t}
!//! Description: humidity at temperature, t. Based on holding vapour \
!//!  pressure constant through the day. Beware if temperature < ZERO \\
!//! Units: Percentage
!//! Type: REAL
!//! SET
!
!//! Variable: VP_{unsat,base}
!//! Description: vapour pressure at reference time eg 9am
!//! Units: Pascal (Pa)
!//! Type: REAL
!
!//! Variable: VP_{sat,t}
!//! Description: Saturated vapour pressure at temperature, t
!//! Units: Pascal (Pa)
!//! Type: REAL
!
!**********************************************************************!
end function HUMIDITY_INST_FN
!**********************************************************************!





!**********************************************************************!
subroutine RAIN_DAY_SUB(NO_RAIN_DAYS, RAIN_IN_MONTH,DAYS_IN_MONTH, SEED1, SEED2, SEED3, RAIN)
!**********************************************************************!
!** The routine simulates daily rainfall for all days in a month
!** The typical number of rainy-days and average rainfall for the month 
!** Are the only specific Meteorology inputs. Returned is an array
!** (1-dimensional; day) with the simulated rainfall for each day 
!** in it.
real(dp) :: RAIN(31)
real(dp) :: FRACT_WET, RAIN_DAY, RAIN_IN_MONTH, NO_RAIN_DAYS
real(dp) :: ZERO
real(dp) :: PROB_WET_WET, PROB_WET_DRY, GAMMA, MU
integer :: DAYS_IN_MONTH, IDAY
integer :: SEED1, SEED2, SEED3, IP, IP0
real(dp)  :: RAND01 , XLO, XHI, gamma2

! **    determine the proportion of the month with wet-days
FRACT_WET = real(NO_RAIN_DAYS)/real(DAYS_IN_MONTH)
!
!**     find the average rainfall per day with rainfall
RAIN_DAY = RAIN_IN_MONTH/real(NO_RAIN_DAYS)
!
!**     find transitional probability of a wet-day following a dry-day
!**     of the month (Geng et al 1986)
PROB_WET_DRY = 0.75*FRACT_WET
!
!**     probability of wet-day following a wet-day (Geng et al 1986)
PROB_WET_WET = 0.25+PROB_WET_DRY
!
!**     Gamma distribution paramaters (Geng et al 1986) 
!       GAMMA = -2.16+1.83*RAIN_DAY   
  XLO = 0.0
  XHI=1.0

  RAND01 = RNDF(XLO,XHI, SEED1, SEED2, SEED3)
  GAMMA = 1.0
!
!**     Mu paramaters (Geng et al 1986) 
  gamma2 = gamma+(RAND01-0.5)/2.0
MU = RAIN_DAY/GAMMA2
!
!     first of the month; the chances of it being wet is 'fract_wet'
!     get a random number proportional to the number of days in the month
XLO = 0.0
XHI=1.0
RAND01 = RNDF(XLO,XHI, SEED1, SEED2, SEED3)
do 20 IDAY = 1, DAYS_IN_MONTH
	if(IDAY==1) then
!**       is the random number<=rain_fract; if so then it rained yesterday
	  IP0 = 0
	  if(RAND01<=FRACT_WET) then
	    IP0 = 1 
	  endif
	endif
	RAIN(IDAY) = 0.0
	RAND01 = RNDF(XLO, XHI, SEED1, SEED2, SEED3)
	if(IP0==1) then
!**       raining on day 0
	  if(RAND01-PROB_WET_WET <= 0) IP = 1
	  if(RAND01-PROB_WET_WET > 0) IP = 0
	else   ! IP0 = 0 ; its dry on day 0
	  if(RAND01-PROB_WET_DRY <= 0) IP = 1
	  if(RAND01-PROB_WET_DRY > 0) IP = 0
	endif
	IP0 = IP
!
	if(IP==1) then     ! Its raining today
12        continue
!12	  TR1 = EXP(-18.42*MU)

!	  TR2 = EXP(-18.42*(1-MU))


!         RAND01 = RNDF(XLO, XHI, SEED1, SEED2, SEED3)
!         IF(RAND01 - TR1 .LE. 0) THEN
!           S1 = 0
!         ELSE
!           S1 = RAND01*EXP(1/MU)
!         ENDIF
!         RAND01 = RNDF(XLO, XHI, SEED1, SEED2, SEED3)
!         IF(RAND01 - TR2 .LE. 0) THEN
!           S2 = 0
!         ELSE
!           S2 = RAND01*EXP(1/(1-MU))
!         ENDIF
!         IF(S1+S2-1 .LE. 0) THEN
!           IF( (S1+S2).EQ.0 )THEN
!             Z = 0.0
!           ELSE
!             Z = S1/(S1+S2)
!           ENDIF
!         ELSE
!           FLAG = 1
!         ENDIF
!         IF(FLAG.EQ.1) GOTO  30
!
!         RAND01 = RNDF(XLO, XHI, SEED1, SEED2, SEED3)
!         RAIN(IDAY) = -Z*LOG(RAND01)*GAMMA
	  ZERO = 0.0

! Changed GAMMA to GAMMA2 14/12/01 Paul Henshall.
  
  	  RAIN(IDAY)=GAMMADIST_FN(zero,MU,GAMMA2,SEED1,SEED2,SEED3)
  if(rain(iday)>rain_day*5)goto 12
	endif
20    continue
!
!//! Equation: 32
!//! Description: The routine simulates daily rainfall for all days in \
!//!  a month. The typical number of rainy-days and average rainfall \
!//!  for the month are the only specific Meteorology inputs. Returned \
!//!  is an array (1-dimensional; day) with the simulated \
!//!  rainfall for each day in it. Although the routine uses  \
!//!  Markov-chain principles, the start of each month is INDEPENDENT \
!//!  of what the conditions were on the last day of the previous  \
!//!  month. \\
!//!  The routine calls a random-number gernerator (RNDF) which \
!//!  returns numbers between Xmin and Xmax (in this case 0,1) using \
!//!  the 3 seeds \\
!//!   INPUTS: Number of rainy days in the month; Average rainfall of \
!//!    the month; number of days in the month; 3 seeds \
!//!    for random number generation\\
!//!   OUTPUT: array (dimension:31), with rainfall for each day
!//! Bibliography: S.Evans (SWELTER)
!
!//! Case: Day_{wet,dry}^{i} = 0
!//!  E: Rain_{day}^{i} = 0.0 \\
!//! \\
!//! Case: Day_{wet,dry}^{i} = 1 \\
!//! \\
!//!  E: Rain_{day}^{i} = Z \gamma log_{e} (Rnd_{j1,j2,j3})  \\
!//!   \gamma = -2.16+1.83 Rain_{day,mean}  \\
!//! \\
!//!  Case: S_{1} + S_{2} -1 <= 0    
!//!   E: Z = { {{S_{1}} \over{S_{1}+S_{2}} }  \\
!//! \\
!//!     Case: Rnd_{j1,j2,j3} - Tr_{1} <= 0
!//!      E: S_{1} = 0
!//!     Case: Rnd_{j1,j2,j3} - Tr_{1} > 0
!//!      E: S_{1} = Rnd_{j1,j2,j3} exp(1/\mu) \\
!//! \\
!//!     Case: Rnd_{j1,j2,j3} - Tr_{2} <= 0
!//!      E: S_{2} = 0
!//!     Case: Rnd_{j1,j2,j3} - Tr_{2} > 0
!//!      E: S_{2} = Rnd_{j1,j2,j3} exp(1/(1-\mu))   \\
!//! \\
!//!  Case: S_{1} + S_{2} -1 > 0
!//!     E: re-generate S_{1} and S_{2}  \\
!//!    \\ 
!//!  Tr_{1} = {exp(-18.42)}\over{\mu} \\
!//!  Tr_{2} = {exp(-18.42)}\over{(1-mu)} \\
!//!  \mu = {{Rain_{day,mean}} \over{\gamma}} \\
!//! \\
!//!  Case: Day_{wet,dry}^{i-1} = 1 \\
!//! \\
!//!   Case: Rnd_{j1,j2,j3} - Prob_{wet,wet} <= 0
!//!    E: Day_{wet,dry}^{i} = 1
!//!   Case: Rnd - Prob_{wet,wet} > 0
!//!    E: Day_{wet,dry}^{i} = 0 \\
!//! \\
!//!  Case: Day_{wet,dry}^{i-1} = 0 \\
!//! \\
!//!   Case: Rnd_{j1,j2,j3} - Prob_{wet,dry} <= 0
!//!    E: Day_{wet,dry}^{i} = 1      
!//!   Case: Rnd_{j1,j2,j3} - Prob_{wet,dry} > 0 
!//!    E: Day_{wet,dry}^{i} = 0      \\
!//! \\
!//!  Prob_{wet,wet} = 0.25+Prob{wet,dry} \\
!//!  Prob_{wet,dry} = { {0.75 N_{rain days}} \over{DinM} }
!
!//! Variable: Rain_{day}^{i}
!//! Description: Predicted daily rainfall on day of month, i
!//! Units: mm
!//! Type: real(dp) array (12,31)
!//! SET
!
!//! Variable: Rain_{day,mean}
!//! Description: Mean daily rainfall for the month
!//! Units: mm/day
!//! Type: REAL
!
!//! Variable: Rnd_{j1,j2,j3)
!//! Description: Random number generated between 0 and 1
!//! Units: none
!//! Type: REAL
!
!//! Variable: j1
!//! Description: Seed for random number generator; changed on return
!//! Units: none
!//! Type: integer
!
!//! Variable: j2
!//! Description: Seed for random number generator; changed on return
!//! Units: none
!//! Type: integer
!
!//! Variable: j3
!//! Description: Seed for random number generator; changed on return
!//! Units: none
!//! Type: integer
!
!//! Variable: DinM
!//! Description: Number of days in the month
!//! Units: none
!//! Type: integer
!
!**********************************************************************!
end subroutine RAIN_DAY_SUB
!**********************************************************************!




 
!**********************************************************************!
subroutine RAIN_DURATION_SUB(RAIN, RAIN_DURATION, RAIN_INTENSITY,DAYS_IN_MONTH)
!**********************************************************************!
!** This routine estimates the duration of rain on each day for a month
!** Accuracy is suspect.
real(dp) :: RAIN(31) , RAIN_DURATION(31), RAIN_INTENSITY(31)
real(dp) :: RAIN_LIM(10), PRANGE, DUR(10)
integer :: DAYS_IN_MONTH, IDAY
integer :: N(10), RANGE
!
RAIN_LIM(1) = 0
RAIN_LIM(2) = 5
RAIN_LIM(3) = 10
RAIN_LIM(4) = 15
RAIN_LIM(5) = 20
RAIN_LIM(6) = 25
RAIN_LIM(7) = 50
RAIN_LIM(8) = 75
RAIN_LIM(9) = 100
do RANGE = 1,10
  N(RANGE) = 0
enddo
!
do 60 IDAY = 1, DAYS_IN_MONTH
	do 70 RANGE = 1, 8,1
	  if((RAIN(IDAY)>RAIN_LIM(RANGE)).and.(RAIN(IDAY)<=RAIN_LIM(RANGE+1))) then
	    N(RANGE) = N(RANGE)+1
	  endif
70      continue
60    continue
!     
do 80 RANGE=1, 8, 1
	PRANGE = (RAIN_LIM(RANGE+1) - RAIN_LIM(RANGE)) / 25.4
	DUR(RANGE) = N(RANGE)/ (1.39*((PRANGE+0.1)** (-3.55)))
80    continue
!
do 90 IDAY = 1,  DAYS_IN_MONTH
	RAIN_DURATION(IDAY) = 0.0
	RAIN_INTENSITY(IDAY) = 0.0
	if(RAIN(IDAY)>0) then
	  do 100 RANGE =1,8,1
	    if( (RAIN(IDAY)>RAIN_LIM(RANGE)) .and. (RAIN(IDAY)<=RAIN_LIM(RANGE+1)) ) then
	      PRANGE = (RAIN_LIM(RANGE+1) - RAIN_LIM(RANGE)) 
	      RAIN_DURATION(IDAY) = DUR(RANGE)*RAIN(IDAY)/PRANGE*60
	    endif
100       continue
	RAIN_INTENSITY(IDAY) = RAIN(IDAY)/(RAIN_DURATION(IDAY)/60.0)
	endif

90    continue
!
!//! Equation: 33
!//! Description: This routine takes a months-worth of daily rainfall \
!//!  data and calculates the duration of each event and intensity of\
!//!  each event. Only one Event occurs each day. Accuracy seems \\
!//!   dubious. \\
!//!   INPUTS: Daily rainfall (array 31); Duaration of daily \
!//!    rainfall [OUTPUT], (array 31); intensity of rainfall [OUTPUT]\
!//!    ,(array 31); number of days in the month
!//! Bibliography: S. Evans (SWELTER)
!
!//! Case: no rain
!//!  E: Rain_{dur}_{i} = 0.0
!//! Case: Rainfall
!//!  E: Rain_{dur}_{i} = {{60 Rain_{dur:c} Rain_{day}_{i}}     \\
!//!     \over {Range_{c}}}   \\
!//!   Range_{c} = Rainfall category range: x <rain> y  :\\
!//!    [<5, 10, 15, 20, 25, 50, 75, 100] mm   \\
!//!  \\
!//!   Rain_{dur:c}={ {\sum_{i=1}^{DinM} Occurences of Rain_{day}_{i} \\
!//!    within rainfall category c}  \over {1.39 \\
!//!    ((Range_{c}/25.4 + 0.1)^{-3.55})} } \\
!//! \\
!//!   Rain_{inten} = Rain_{day}^{i}/(Rain_{dur}^{i}/60.0)
!//!
!
!//! Variable: Rain_{dur}_{i}
!//! Description: Duration of rainfall on day of month, i
!//! Units: minutes
!//! Type: real(dp) array (31)
!//! SET
!
!//! Variable: Rain_{inten}_{i}
!//! Description: Rainfall intensity on day of month, i
!//! Units: mm /hour
!//! Type: REAL

!//! Variable: Rain_{day}_{i}
!//! Description: Rainfall on day of month, i
!//! Units: mm 
!//! Type: REAL

!//! Variable: DinM
!//! Description: Number of days in the month
!//! Units: none
!//! Type: integer
!
!**********************************************************************!
end subroutine RAIN_DURATION_SUB
!**********************************************************************!





!**********************************************************************!
real(dp) function RNDF(XLO, XHI,J1, J2, J3)
!**********************************************************************!
!** This random number routine w is from Brian Wichman and David Hill
!** (BYTE, March 1987, page 127-128). Apparently the cycle length is
!** about 6.95E12. Original reference: Wichmann & Hill, App stats (1982)        
!** 31(2) pp 188-190 with modification: A Mcleod, App stats (1985)            
!** 34(2) pp 198-200                                             
!
integer :: J1, J2, J3
real(dp) :: R,XLO,XHI
!                                                                   
if(XLO<XHI)then                                            
	J1=171*MOD(J1,177)-2*(J1/177)                                  
	if(J1<0) J1 = J1+30269                                         
!                                                                   
	J2=172*MOD(J2,176)-35*(J2/176)                                 
	if(J2<0)J2=J2+30307                                           
!                                                                     
	J3=170*MOD(J3,178)-63*(J3/178)  
  if(J3<0) J3 = J3+30323                                           
!                                                                     
	R=real(J1)/30269.0+real(J2)/30307.0+real(J3)/30323.0       
!** modification re: A McLeod (1985) app. stats 34(2) 198-200     
!** returns out of loops for efficiency - otherwise it was just   
!**  a question of scaling.                                        
	if(R>0) then                                               
	  RNDF=(R-INT(R))*(XHI-XLO)+XLO                               
	  return                                                      
	endif                                                         
	R=MOD(real(FLOAT(J1))/30269.0 + real(FLOAT(J2))/30307.0 +  real(FLOAT(J3))/30323.0, 1.0)                           
	if(R>=1.0) R=0.999999                                       
	RNDF=(R-INT(R))*(XHI-XLO)+XLO                                 
	return                                                        
endif                                                           
!                                                                     
RNDF=XLO                                                        
!
!//! Equation: 34
!//! Description: Random number generator: uses 3 seeds and returns \
!//!  changed seeds and a random number between Xlo and Xhi\\
!//!   INPUTS: Lowest value of random number; highest value of random \
!//!    number; 3 seeds
!//! Bibliography: Wichmann & Hill (1882), App stats 31(2) pp188-190; \
!//!  A McLeod (1985) app. stats 34(2) pp98-200
!
!//! Case: R > 0
!//!  E: Rnd_{j1,j2,j3}=(R-INT(R))*(X_{hi}-X_{lo})+X_{lo}
!//! Case: R<= 0
!//!  E: RNDF=(R1-INT(R1))*(X_{hi}-X_{lo})+X_{lo}  \\
!//! \\
!//! E: R = j1/30269 + j2/30307 + j3/30323  \\
!//!   \\
!//! E: R1 = MAX(0.999999, MOD(S_{1}/30269 + S_{2}/30307 + \
!//!     S_{3}/30323.0, 1.0)                           \\
!//! \\
!//!   E:  j1 = 171 MOD(j1,177)-2(j}/177)  \\
!//!     IF(j1 < 0  ; j1 = j1 + 30269   \\
!//!  \\
!//!   E:  j2 = 172 MOD(j2,176)-35(j2/176)  \\
!//!     IF(j2 < 0  ; j2 = j2 + 30307   \\
!//!  \\
!//!   E:  j3 = 170 MOD(j3,178)-63(j3/178)  \\
!//!     IF(j3 < 0  ; j3 = j3 + 30323   
!  
!//! Variable: Rnd_{j1,j2,j3}
!//! Description: Returned random number
!//! Units: none
!//! Type: REAL
!//! SET
!
!//! Variable: j1
!//! Description:  Seed for random number generation changed on exit
!//! Units: none
!//! Type: integer
!//! SET
!
!//! Variable: j2
!//! Description:  Seed for random number generation changed on exit
!//! Units: none
!//! Type: integer
!//! SET
!
!//! Variable: j3
!//! Description:  Seed for random number generation changed on exit
!//! Units: none
!//! Type: integer
!//! SET
!
!//! Variable:  X_{lo}
!//! Description:  lowest possible value of random number
!//! Units: none
!//! Type: REAL
!
!//! Variable:  X_{hi}
!//! Description:  highest possible value of random number
!//! Units: none
!//! Type: REAL
!
!**********************************************************************!
end function RNDF
!**********************************************************************!





!**********************************************************************!
real(dp) function TRANS_PROP_DAY_FN(RAD_ET_DAY,RAD_TERRA_SW_DAY)
!**********************************************************************!
!** Evaluates the total trasmission proportion (coefficient) for the day
real(dp) :: RAD_TERRA_SW_DAY
real(dp) :: RAD_ET_DAY
!
TRANS_PROP_DAY_FN = RAD_TERRA_SW_DAY/RAD_ET_DAY
!
!//! Equation: 35 
!//! Description: Works out the proportion of Extra-terrestrial \
!//!  radiation that arrives on the earth's surface as short-wave \
!//!  in a day \\
!//!   INPUTS: Daily short-wave radiation on earth's surface; daily \
!//!    Extra-terrestrial radiation.
!//! Bibliography: S. Evans (SWELTER); Liu & Jordan (1960) Solar Energy \
!//!  4 pp1-19.
!
!//! E:  T_{sw} = {RAD_{terra, sw} /over{RAD_{et}}}
!
!//! Variable: T_{sw}
!//! Description: proportion of Extra-terrestrial radiation that \
!//!  arrives on the earth's surface as short-wave.
!//! Units: none
!//! Type: REAL
!//! SET
!
!//! Variable: RAD_{terra,sw}
!//! Description: Terrestrial short-wave radiation on a horizonatal \
!//!  plane for the given latitude and day.
!//! Units: W/m2/d
!//! Type: REAL
!                                                                    
!//! Variable: RAD_{et}   
!//! Description: Extra-Terrestrial radiation above the atmosphere  \
!//!  for the given latitude and day.
!//! Units: W/m2/d
!//! Type: real(dp)
!
!**********************************************************************!
end function TRANS_PROP_DAY_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function ANGSTROM_BETA_FN(LATITUDE_RAD)
!**********************************************************************!
!** finds the Angstrom turbidity factor (the max clear-sky atmospheric
!** transmittance characteristics at the latitude
real(dp) :: LATITUDE_RAD
!
ANGSTROM_BETA_FN = 0.682 - 0.3183*(1-1.3614*COS(LATITUDE_RAD))
!
!//! Equation: 36
!//! Description: Caculation of the Angstrom beta function value: A \
!//!  turbidity factor for the maximum clear-sky transmittance \
!//!  characteristics at a given latitude. \\
!//!   INPUTS: Latitude
!//! Bibliography: Nikolov & Zeller (1992) Ecological Modelling \
!//!  61:149-168; Gueymard (1989) Agric & Forest Met 45:215-229
!
!//! E: \beta =  0.682 - 0.3183*(1-1.3614*cos(Lat))
!
!//! Variable: \beta
!//! Description: Angstrom tubidity factor: max clear sky transmittance
!//! Units: none  
!//! Type: REAL
!//! SET
!
!//! Variable: Lat
!//! Description:  Latitude of site (in radians) 
!//! Units: radians      
!//! Type: REAL
!
!**********************************************************************!
end function ANGSTROM_BETA_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function DIFFPROP_DAY2_FN(TRANS_PROP_DAY,BETA)
!**********************************************************************!
!**Estimates the proportion of terrestrial radiation that is diffuse
real(dp) :: TRANS_PROP_DAY, BETA
!
DIFFPROP_DAY2_FN = TRANS_PROP_DAY*(1-exp((-0.6*(1-BETA)/TRANS_PROP_DAY) / (BETA-0.4)) )
!
!//! Equation: 37
!//! Description: The proportion of light reaching the ground that \
!//!   is made up of diffuse component                      \\
!//!   INPUTS: Daily extra-terrestrial radiation; Daily terrestrial \
!//!    radiation; Angstrom beta value
!//! Bibliography: S Evans (SWELTER); Bristow & Campbell (1985) Agric \
!//!  & Forest Met 35:123-131; Becker & Weingarten (1991) Agric &  \
!//!  Forest Met 53:347-353
!
!//! E: Diff_{prop} =   Tran(1-exp{({{-0.6(1-\beta)} \over{Tran}} \
!//!     \over{\beta-0.4})}
!
!//! Variable: Diff_{prop}                                           
!//! Description: Proportion of daily terrestrial radiation that is \
!//!   diffuse                                              \\       
!//!   INPUTS: Extra-terrestrial daily radiation; Terrestrial daily \
!//!   radiation                                                     
!//! Units: none (Ratio)                                             
!//! Type: real(dp)
!//! SET                                                             
!                                                                    
!//! Variable: Tran                                                  
!//! Description: Proportion of Terrestrial radiation to    \
!//!   Extra-terrestrial radiation                                   
!//! Units: none (Ratio)                                             
!//! Type: real(dp)
!                                                                    
!//! Variable: \beta
!//! Description: Angstrom tubidity factor: max clear sky transmittance
!//! Units: none  
!//! Type: REAL
!
!**********************************************************************!
end function DIFFPROP_DAY2_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function RAD_TERRA_INST_FN(DAWN, DUSK, DAYLENGTH, HOUR, RAD_TERRA_DAY, PI)
!**********************************************************************!
!** calculate solar radiation for  at any instant in time        
!** Assume sine wave with maximum at noon so ideally use local solar time
real(dp) :: RAD_TERRA_DAY
real(dp) :: DAWN, DUSK, DAYLENGTH, HOUR, SUNTIME
real(dp) :: PI
!
!**Give zero values for radiation at night otherwise set radiation
!
if(HOUR < DAWN .or. (HOUR > DUSK)) then
	  RAD_TERRA_INST_FN = 0.0
else
!** Sinusoidal function for hourly radiation from daily total    
	SUNTIME = HOUR - DAWN
	RAD_TERRA_INST_FN = RAD_TERRA_DAY*PI/(2*DAYLENGTH)*SIN(PI/DAYLENGTH*SUNTIME)/3600.0 
endif
!
!//! Equation: 38
!//! Description: Fit a sine curve over the day - with peak at noon \
!//!  (so local solar time is best) and periodicity over the  \
!//!  daylength \\
!//!   INPUTS:  Latitude; Time of dawn; Time of dusk; Daylength; \
!//!    Time of day; Daily terrestrial radiation ; Pi
!//! Bibliography: Forest-GROWTH/FLUX; Meastro
!
!//! Case: Night
!//!  E: RAD_{now} = 0.0
!//! Case: Dawn <= hour <= Dusk
!//! E: RAD_{now} = {{RAD_{terra, day} \pi sin(\pi/(DAYL T)} \
!//!       \over{(2 DAYL 3600)}}
!
!//! Variable: RAD_{now}
!//! Description: terrestrial radiation at any point in time; attained \
!//!  from fitting a sine curve between dawn and dusk, peak at noon \
!//!  using the daily radiation. Local solar times should be used
!//! Units: W/m^{2} (shortwave radiation)
!//! Type: REAL
!//! SET
!
!//! Variable: Dawn
!//! Description: Time of sunrise 
!//! Units: Hour
!//! Type: REAL
!
!//! Variable: Dusk
!//! Description: Time of sunset  
!//! Units: Hour
!//! Type: REAL
!
!//! Variable: RAD_{terra, day}
!//! Description:  Daily terrestrial radiation 
!//! Units: W/m^{2}/day
!//! Type: REAL
!
!//! Variable: DAYL
!//! Description:  Daylength                   
!//! Units: Hours
!//! Type: REAL
!
!//! Variable: hour
!//! Description: Hour of the day - GMT base
!//! Units: hour
!//! Type: REAL
!
!**********************************************************************!
end function RAD_TERRA_INST_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function WINDSPEED_INST_FN(WIND_RUN_DAY, HOUR)
!**********************************************************************!
!** calculates windspeed at any point in time from a daily run of
!** wind. This is based on a regression to get the daily max and 
!** min windspeed from the run of wind. Winspeed is assumed to be 
!** maximum at 14:00, and changes linearly about that
real(dp) :: WIND_RUN_DAY, HOUR
real(dp) :: WSMIN, WSMINCONST, WSMINCOEFF
real(dp) :: WSMAX, WSMAXCONST, WSMAXCOEFF
!
! ** this is parmaterised from Alice -Holt 1999: Regressing 
!** Run Vs WSMIN gives the coeffs: 0.0075x -0.3181, r2=0.7559
!** Run vs WSMAX gives the coeffs, 0.0142x +0.6976, r2=0.8241
WSMINCONST = -0.3181 ! -0.2253
WSMINCOEFF = 0.0075  !  0.0391
WSMAXCONST = 0.6976  !  0.2822
WSMAXCOEFF = 0.0142  !  0.0444
WSMIN=MAX(0.0,WIND_RUN_DAY*WSMINCOEFF+WSMINCONST)    
WSMAX=MAX(0.0,WIND_RUN_DAY*WSMAXCOEFF+WSMAXCONST)
!
if(HOUR<=13.0) then                  
	WINDSPEED_INST_FN = MAX(0.1,WSMIN+(WSMAX-WSMIN)/13.*(HOUR-1))
else                                  
	WINDSPEED_INST_FN = MAX(0.1,WSMAX-(WSMAX-WSMIN)/11.*(HOUR-13))
endif
!
!//! Equation: 39
!//! Description: Approximates wind-speed at any time of the day \
!//! from the daily run of wind. \\
!//!   INPUTS: Daily run of wind; Time of day
!//! Bibliography: none ( M Broadmeadow: developed from regression at \
!//!  the Straits)
!
!//! Case: Hour <=13
!//!  E: Wind_{now} = max(0.1, Wind_{min}+(Wind_{max}-Wind_{min}) \
!//!     (Hour-1)/13)
!//! Case: Hour >13
!//!  E: Wind_{now} = max(0.1, Wind_{max}+(Wind_{max}-Wind_{min}) \
!//!     (Hour-13)/11) \\
!//!  \\
!//! Wind_{max} = Wind_{run} WS_{max,1} + WS_{max,2} \\
!//!  \\
!//! Wind_{min} = Wind_{run} WS_{min,1} + WS_{min,2} 
!
!//! Variable: Wind_{now}
!//! Description: Windspeeds at current time
!//! Units: m/s
!//! Type: REAL
!//! SET
!
!//! Variable: WS_{max,1}
!//! Description:  Linear regression coefficient of max windspeed Vs \
!//!  run of wind
!//! Units: Km/g / m/s
!//! Type: REAL
!
!//! Variable: WS_{max,2}
!//! Description:  Linear regression constant of max windspeed Vs \
!//!  run of wind
!//! Units: m/s
!//! Type: REAL
!
!//! Variable: WS_{min,1}
!//! Description:  Linear regression coefficient of min windspeed Vs \
!//!  run of wind
!//! Units: Km/g / m/s
!//! Type: REAL
!
!//! Variable: WS_{min,2}
!//! Description:  Linear regression constant of min windspeed Vs \
!//!  run of wind
!//! Units: m/s
!//! Type: REAL
!
!//! Variable: Wind_{run}
!//! Description:  Daily run of wind
!//! Units: Km/d
!//! Type: REAL
!
!//! Variable: Hour
!//! Description: Hour of the day - GMT base
!//! Units: hour
!//! Type: REAL
!
!**********************************************************************!
end function WINDSPEED_INST_FN
!**********************************************************************!




 
!**********************************************************************!
real(dp) function RAD_TERRA_SUNDAY_FN(DAYLENGTH,SUN_HOURS,RAD_ET_DAY)
!**********************************************************************!
!**  Predicts terrestrial radiation from sun-shine hours
!**  Rad_terra= a + b*Sun/daylen + c*(Sun/daylen)**2 + d*(Sun/DAYLEN)**3
!**  Other conversions eg Angstom 1924; Collingbourne 1976 use a linear 
!**  Function; ie 2 params
real(dp) :: CONST, COEFFA, COEFFB, COEFFC, SRATIO
real(dp) :: DAYLENGTH, PRED_GRATIO
real(dp) :: RAD_ET_DAY, SUN_HOURS
!                                                             
CONST   = 0.17262                                       
COEFFA   =0.9901                                        
COEFFB   =-1.043                                        
COEFFC   = 0.585                                        
!                                                             
	  SRATIO = SUN_HOURS/DAYLENGTH
	  if(SRATIO>1) then                                
	    SRATIO = 1                                        
	  endif                                               
!
	 PRED_GRATIO = CONST + COEFFA*SRATIO + COEFFB*SRATIO**2 + COEFFC*SRATIO**3 
	 RAD_TERRA_SUNDAY_FN = RAD_ET_DAY*PRED_GRATIO !J/day
!
!//! Equation: 40
!//! Description: Method of predicting terrestrial radiation from \
!//!  sun-shine hours. Note often only \alpha and \beta are used in  \
!//!  the equation below:  \\
!//!   INPUTS: Daylength; Sunhine-hours; Daily Extra-terrestrial \
!//!   radiation; Mathematical constant, pi
!//! Bibliography: Angstrom (1924) Q.J. R.met. Soc 50,121-5; \
!//!  Collingbourne (1976)  In: The Climate of the British Isles. \
!//!  eds TJ Chandler & S Gregory, Longman, pp74-95; Randle (1997) \
!//!  FC Internal report
!
!//! E: RAD_{terra,day} = RAD_{et,day} (\alpha + \beta*S/S_{0} +  \
!//!  \gamma (S/S_{0})^{2} + \delta (S/S_{0})^{2} )
!
!//! Variable:  RAD_{terra,day}
!//! Description: Daily terrestrial radiation predicted from  \
!//!  sun-shine hours
!//! Units: W/m^{2}/d
!//! Type: REAL
!//! SET
!
!//! Variable: RAD_{et,day}
!//! Description: Daily extra-terrestrial radiation
!//! Units: W/m^{2}/d
!//! Type: REAL
!
!//! Variable: S
!//! Description: Sun-shine hours in the day
!//! Units: hours     
!//! Type: REAL
!
!//! Variable: S_{0}
!//! Description: Daylength
!//! Units: hours     
!//! Type: REAL
!
!//! Variable: /alpha
!//! Description: Constant of regression of  S/S_{0} and Rad_{et,day}. \
!//!  If linear regression this is 'Angstrom \alpha'
!//! Units: W/m^{2}/d 
!//! Type: REAL
!
!//! Variable: /beta 
!//! Description: First order coefficient of regression of  S/S_{0}  \
!//! and Rad_{et,day}. If linear regression this is 'Angstrom \beta'
!//! Units: W/m^{2}/d 
!//! Type: REAL
!
!//! Variable: /gamma
!//! Description:  Second order coefficient of regression of  S/S_{0}. \
!//!  and Rad_{et,day}. If linear regression this is zero
!//! Units: W/m^{2}/d 
!//! Type: REAL
!
!//! Variable: /delta 
!//! Description: Constant of regression of  S/S_{0} and Rad_{et,day}. \
!//!  if linear regression this is zero
!//! Units: W/m^{2}/d 
!//! Type: REAL
!
!//! Variable: \pi
!//! Description: Mathematical constant, pi
!//! UNits: (Radian)
!//! Type: REAL
!
!**********************************************************************!
end function RAD_TERRA_SUNDAY_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function TEMP_AIR_AMPLITUDE_FN(TRANS_PROP_DAY,A_SKY, ANGSTROM_BETA, C_SKY)
!**********************************************************************!
!** determines apriximate air temperature amplitude (tmax- tmin) from
!** transmitted radiation and some parameters
real(dp) :: TRANS_PROP_DAY
real(dp) :: A_SKY, ANGSTROM_BETA, C_SKY, DELTA

DELTA = 1-TRANS_PROP_DAY/ANGSTROM_BETA
if(DELTA<=0) then 
	TEMP_AIR_AMPLITUDE_FN = 0.0
else
	TEMP_AIR_AMPLITUDE_FN = (LOG(DELTA)/(-A_SKY))**(1/C_SKY)
endif
!
!//! Equation: 41
!//! Description: Works out the daily temperature differecne amplitude \
!//!  from the transmittance of radiation; if it is cloudy, then less \
!//!  radiation gets through and there should be less temperature \
!//!  variation \\
!//!   INPUTS: Proportion of Extra-terrestrial light reaching the \
!//!    earth; Coefficient for maximum clear-sky transmittance;  \
!//!   Angstrom \beta value; Coefficient of clear-sky transmittance \
!//!   with amplitude temperature increase. 
!//! Bibliography: S Evans (SWELTER); Bristow & Campbell (1984) Agric \
!//!  & Forest Met 31: 159-166
!
!//! Case: 1-{{T_{sw}} \over{\beta}} <=0
!//!  E: \Delta T = 0.0
!//! Case:  1-{{T_{sw}} \over{\beta}} > 0
!//!  E: \Delta T = (Log_{e}(\delta/(-A_{sky})) )^(1/C_{sky})
!
!//! Variable: \Delta T
!//! Description: Daily amplitude of temperature variation
!//! Units: Degrees C  (usually)
!//! Type: REAL
!//! SET
!
!//! Variable: T_{sw}
!//! Description: proportion of Extra-terrestrial radiation that \
!//!  arrives on the earth's surface as short-wave.
!//! Units: none
!//! Type: REAL
!
!//! Variable: \beta
!//! Description: Angstrom tubidity factor: max clear sky transmittance
!//! Units: none  
!//! Type: REAL
!
!//! Variable: A_{sky}
!//! Description: Coefficient of maximum clear-sky transmittance \
!//!  characteristics
!//! Units: none  
!//! Type: REAL
!
!//! Variable: C_{sky}
!//! Description: Coefficient of maximum clear-sky transmittance \
!//!   with \Delta T increase
!//! Units: none  
!//! Type: REAL
!
!**********************************************************************!
end function TEMP_AIR_AMPLITUDE_FN
!**********************************************************************!



 

!**********************************************************************!
real(dp) function DIFF_PROP_HOUR_FN(DIFF_PROP_DAY,SOL_ELEV)
!**********************************************************************!
!** finds the proportion of idiifuse light at the earth's surface on
!** an hourly basis (ie elevation should be calculated hour by hour
!** and this feeds back the diffuse light proportion for this hour
real(dp) :: DIFF_PROP_DAY, R, K
real(dp) :: SOL_ELEV
!
R = 0.847-1.61*SIN(SOL_ELEV) + 1.04*(SIN(SOL_ELEV)**2)
K = (1.47-R)/1.66         
if(DIFF_PROP_DAY<0.22) then 
	 DIFF_PROP_HOUR_FN = 1.0             
else if(DIFF_PROP_DAY<0.35) then     
	 DIFF_PROP_HOUR_FN = 1-6.4*(DIFF_PROP_DAY-0.22)**2 
else if(DIFF_PROP_DAY<K) then         
	 DIFF_PROP_HOUR_FN = 1.47-1.66*DIFF_PROP_DAY
else                                                          
	 DIFF_PROP_HOUR_FN = R                                             
endif                                                         
!
!//! Equation: 42
!//! Description: Finds the proportion of diffuse/total light at the \
!//!  earth's surface on an hourly basis. Elevation should change \
!//!  hour-by-hour. Usually used when hourly terrestrial radiation \
!//!  levels are known. \\
!//!   INPUTS: Proportion of diffuse light (daily basis); \
!//!    Solar elevation 
!//! Bibliography: Spitters Et al (1986), Agric & Forest Met  \
!//!  38:217-229; de Jong (1980).                                
!
!//! Case: Diff_{prop} <= 0.22
!//!  E: Diff_{prop,h} = 1.0
!//! Case: 0.22 < Diff_{prop} <= 0.35
!//!  E: Diff_{prop,h} = 1 - 6.4(Diff_{prop}-0.22)^{2}
!//! Case: 0.35 < Diff_{prop} <= K  
!//!  E: Diff_{prop,h} = 1.47 - 1.66 Diff_{prop}
!//! Case: K < Diff_{prop} 
!//!  E: Diff_{prop,h} = R     \\
!//!   \\
!//! E: K = (1.47-R)/1.66  \\
!//!  \\
!//! E: R = 0.847-1.61 sin(\beta)
!
!//! Variable: Diff_{prop,h}                                            
!//! Description: Proportion of hourly terrestrial radiation that is \
!//!   diffuse at hour h                                    \\
!//!   INPUTS: Proportion of daily radiation that is diffuse; Solar \
!//!    elevation at hour, h
!//! Units: none (Ratio)                                              
!//! Type: real(dp)
!//! SET
!
!//! Variable: Diff_{prop}                                            
!//! Description: Proportion of daily terrestrial radiation that is \
!//!   diffuse                                              \\
!//! Units: none (Ratio)                                              
!//! Type: real(dp)
!                                                                     
!//! Variable: \beta                                                  
!//! Description:  Solar elevarion at hour, h
!//! Units: Radians
!//! Type: real(dp)
!
!**********************************************************************!
end function DIFF_PROP_HOUR_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function TEMP_MAX_DAY_FN(TEMP_AIR_AMPLITUDE,TEMP_DAY_MEAN)
!**********************************************************************!
!** predicts max air temp from the mean temp and amplitude: based on 
!** tmean = (tmax-tmin)/2
real(dp) :: TEMP_AIR_AMPLITUDE
real(dp) :: TEMP_DAY_MEAN
!
TEMP_MAX_DAY_FN = TEMP_DAY_MEAN + (TEMP_AIR_AMPLITUDE)/2.0
!
!//! Equation: 43
!//! Description: Derives maximum air temperature from the daily mean \
!//!  and daily amplitude. Based on Tmean=(tmax-Tmin)/2  \\
!//!   INPUTS: Daily temperature variation amplitude; Daily mean \
!//!    temperature
!//! Bibliography: none
!
!//! E: T_{max} = T_{mean} + \Delta T /2
!
!//! Variable: T_{max}
!//! Description: Daily maximum temperature
!//! Units: Degrees C   (usually)
!//! Type: Double precsion
!//! SET
!
!//! Variable: T_{mean}
!//! Description: Daily mean temperature
!//! Units: Degrees C   (usually)
!//! Type: Double precsion
!
!//! Variable: \Delta T 
!//! Description: Daily temperature variation amplitude
!//! Units: Degrees C   (usually)
!//! Type: Double precsion
!
!**********************************************************************!
end function TEMP_MAX_DAY_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function TEMP_MIN_DAY_FN(TEMP_AIR_AMPLITUDE,TEMP_DAY_MEAN)
!**********************************************************************!
!** predicts min air temp from the mean temp and amplitude: based on 
!** tmean = (tmax-tmin)/2
real(dp) :: TEMP_AIR_AMPLITUDE
real(dp) :: TEMP_DAY_MEAN
!
TEMP_MIN_DAY_FN = TEMP_DAY_MEAN - (TEMP_AIR_AMPLITUDE)/2.0
!
!//! Equation: 44
!//! Description: Derives minimum air temperature from the daily mean \
!//!  and daily amplitude. Based on Tmean=(tmax-Tmin)/2  \\
!//!   INPUTS: Daily temperature variation amplitude; Daily mean \
!//!    temperature
!//! Bibliography: none
!
!//! E: T_{min} = T_{mean} + \Delta T /2
!
!//! Variable: T_{min}
!//! Description: Daily minimum temperature
!//! Units: Degrees C   (usually)
!//! Type: Double precsion
!//! SET
!
!//! Variable: T_{mean}
!//! Description: Daily mean temperature
!//! Units: Degrees C   (usually)
!//! Type: Double precsion
!
!//! Variable: \Delta T 
!//! Description: Daily temperature variation amplitude
!//! Units: Degrees C   (usually)
!//! Type: Double precsion
!
!**********************************************************************!
end function TEMP_MIN_DAY_FN
!**********************************************************************!




!**********************************************************************!
real(dp) function TEMP_AIR_AMPLITUDE2_FN(TEMP_MAX_DAY,TEMP_MIN_DAY)
!**********************************************************************!
!** gets the daily temperature variation: ie. Tmax-Tmin
real(dp) :: TEMP_MAX_DAY
real(dp) :: TEMP_MIN_DAY
!
TEMP_AIR_AMPLITUDE2_FN = TEMP_MAX_DAY - TEMP_MIN_DAY
!
!//! Equation: 45
!//! Description: Derives the daily temperature variation from the \
!//!  maximum and minimums of the day \\
!//!  INPUTS: Temperature maximum; Temperature minimum
!//! Bibliography: none
!
!//! E: \Delta T = T_{max} - T_{min}
!
!//! Variable: \Delta T
!//! Description: Daily amplitude of temperature variation
!//! Units: Degrees C  (usually)
!//! Type: REAL
!//! SET
!
!//! Variable: T_{min}
!//! Description: Daily minimum temperature
!//! Units: Degrees C   (usually)
!//! Type: Double precsion
!
!//! Variable: T_{max}
!//! Description: Maximum daily temperature
!//! Units: Degrees C   (usually)
!//! Type: Double precsion
!
!**********************************************************************!
end function TEMP_AIR_AMPLITUDE2_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function TEMP_DAY_MEAN_FN(TEMP_MAX_DAY,TEMP_MIN_DAY)
!**********************************************************************!
!** calculates mean daily temperature based on (Tmin+Tmax)/2
real(dp) :: TEMP_MAX_DAY, TEMP_MIN_DAY
!
TEMP_DAY_MEAN_FN = (TEMP_MIN_DAY + TEMP_MAX_DAY)/2.0
!
!//! Equation: 46
!//! Description: calculates Avrage daily temperature from the daily \
!//!  minimum and maximum  \\
!//!   INPUTS: Daily maximum temperature; Daily minimum temperature
!//! Bibliography: none
!
!//! E: Temp_{mean} = (Temp_{min}+Temp_{max})/2.0
!
!//! Variable: T_{mean}
!//! Description: Daily mean temperature
!//! Units: Degrees C   (usually)
!//! Type: Double precsion
!//! SET
!
!//! Variable: T_{min}
!//! Description: Daily minimum temperature
!//! Units: Degrees C   (usually)
!//! Type: Double precsion
!
!//! Variable: T_{max}
!//! Description: Maximum daily temperature
!//! Units: Degrees C   (usually)
!//! Type: Double precsion
!
!**********************************************************************!
end function TEMP_DAY_MEAN_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function TEMP_DAY_MEAN2_FN(TEMP_MEAN_MONTH,TEMP_SD_MONTH,TEMP_AUTOCORR_MONTH,TEMP_YESTERDAY,IS1, &
 IS2, IS3, NORMSWITCH)
!**********************************************************************!
!** calculates mean daily temperature  from monthly values   
real(dp) :: TEMP_MEAN_MONTH
real(dp) :: TEMP_SD_MONTH, TEMP_AUTOCORR_MONTH
real(dp) :: TEMP_YESTERDAY
real(dp) :: NORM, NORMMN, NORMSD
integer ::IS1, IS2, IS3, NORMSWITCH
! 
!** approximate a function to get a normal distribution, mean=0, sd=1
NORMMN=0.0
NORMSD = 1.0
NORM=NORM_FROM_RAND_FN(NORMMN, NORMSD, IS1, IS2, IS3, NORMSWITCH)
!
TEMP_DAY_MEAN2_FN=TEMP_MEAN_MONTH+TEMP_AUTOCORR_MONTH*(TEMP_YESTERDAY-TEMP_MEAN_MONTH)+TEMP_SD_MONTH* &
 NORM*( (1-(TEMP_AUTOCORR_MONTH)**2)**0.5)
!
!//! Equation: 47
!//! Description: Calculates Average daily temperature from monthly \
!//!  means and standard-deviations. Autocorrelation of 0.65 seems \
!//!  to work if one isn't available. \\
!//!   INPUTS: Mean temperature of the month; Standard deviation \
!//!    around the mean on a daily basis; Temperature \
!//!     auto-correlation for each month; Mean temperature of the \
!//!     previous day; 3 random number seeds.
!//! Bibliography: S Evans (SWELTER); Haith, Tubbs & Pickering (1984) \
!//!  Simulation of Pollution by soil erosion & soil nutirnt loss, \
!//!   PUDOC, Wageningen
!
!//! E: Temp_{mean} = Temp_{mean,month} + \rho_{T} (Temp_{day-1} - \
!//!     Temp_{mean,month}) + \sigma_{mT} N (1-\sigma_{mT}^{2})^0.5) \\
!//!  \\
!//!   E: N = { Random, N(0,1) }
!
!//! Variable: Temp_{mean}
!//! Description: Daily mean temperature
!//! Units: Degrees C   (usually)
!//! Type: Double precsion
!//! SET
!
!//! Variable: Temp_{mean,month}
!//! Description: Mean daily temperature (monthly basis)
!//! Units: Degrees C   (usually)
!//! Type: Double precsion
!
!//! Variable: Temp_{day-1}     
!//! Description: Mean daily temperature of the previous day
!//! Units: Degrees C   (usually)
!//! Type: Double precsion
!
!//! Variable: \rho_{T}
!//! Description: Observed first-order auto-correlation of mean \
!//!  daily air temperature for each month
!//! Units: none  (day/day^{-1})
!//! Type: Double precsion
!
!//! Variable: \sigma_{mT}
!//! Description: Standard deviation of daily air temperature about \
!//!  the mean (monthly basis).
!//! Units: Degrees C             
!//! Type: Double precsion
!
!//! Variable: R
!//! Description: Random number between 0,1
!//! Units: none
!//! Type: Double precsion
!
!**********************************************************************!
end function TEMP_DAY_MEAN2_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function WIND_DAY_MEAN_FN(WIND_MONTH_MEAN, &
 WIND_SD_MONTH, WIND_AUTOCORR_MONTH, WIND_YESTERDAY, IS1, IS2, IS3, NORMSWITCH)
!**********************************************************************!

!** Generates windspeed at a daily scale from monthly means, standard 
!** devaitions and autocorrelation (monthly)
real(dp) :: WIND_MONTH_MEAN, WIND_SD_MONTH
real(dp) :: WIND_AUTOCORR_MONTH, WIND_YESTERDAY
real(dp) :: NORM, NORMMN, NORMSD
integer ::IS1, IS2, IS3, NORMSWITCH
!
NORMMN = 0.0
NORMSD =1.0
!** approximate a function to get a normal distribution, mean=0, sd=1
NORM = NORM_FROM_RAND_FN(NORMMN, NORMSD,IS1,IS2,IS3,NORMSWITCH)
!
WIND_DAY_MEAN_FN = WIND_MONTH_MEAN + WIND_AUTOCORR_MONTH* &
 (WIND_YESTERDAY - WIND_MONTH_MEAN) + WIND_SD_MONTH* &
 NORM*( (1-(WIND_AUTOCORR_MONTH)**2)**0.5)
WIND_DAY_MEAN_FN = MAX(WIND_DAY_MEAN_FN, 0.1)
!
!//! Equation: 48
!//! Description: Calculates Average daily windspeed from monthly \
!//!  means and standard-deviations. Autocorrelation of 0.65 seems \
!//!  to work if one isn't available. \\
!//!   INPUTS: Mean daily windspeed of the month; Standard deviation \
!//!    around the mean on a daily basis; windspeed  \
!//!     auto-correlation for each month; Mean windspeed of the \
!//!     previous day; 3 random number seeds.
!//! Bibliography: S Evans (SWELTER); Haith, Tubbs & Pickering (1984) \
!//!  Simulation of Pollution by soil erosion & soil nutirnt loss, \
!//!   PUDOC, Wageningen
!
!//! E: Wind_{mean} = Wind_{mean,month} + \rho_{W} (Wind_{day-1} - \
!//!     Wind_{mean,month}) + \sigma_{mW} N (1-\sigma_{mW}^{2})^0.5) \\
!//!  \\
!//!   E: N = { Random, N(0,1) }
!
!//! Variable: Wind_{mean}
!//! Description: Daily mean winsdpeed   
!//! Units: m/s 
!//! Type: Double precsion
!//! SET
!
!//! Variable: Wind_{mean,month}
!//! Description: Mean daily windspeed (monthly basis)
!//! Units: m/s
!//! Type: Double precsion
!
!//! Variable: Wind_{day-1}     
!//! Description: Mean daily windspeed of the previous day
!//! Units: m/s
!//! Type: Double precsion
!
!//! Variable: \rho_{W}
!//! Description: Observed first-order auto-correlation of mean \
!//!  daily windspeed for each month
!//! Units: none  (day/day^{-1})
!//! Type: Double precsion
!
!//! Variable: \sigma_{mW}
!//! Description: Standard deviation of daily windspeed about \
!//!  the mean (monthly basis).
!//! Units: m/s 
!//! Type: Double precsion
!
!//! Variable: R
!//! Description: Random number between 0,1
!//! Units: none
!//! Type: Double precsion
!
!**********************************************************************!
end function WIND_DAY_MEAN_FN
!**********************************************************************!



     

!c commented by Ghislain 30/06/03 because it does not compile... and I don't need it   
!cc      real(dp) FUNCTION CLOUD_COVER_FN(VP_SAT,RAINFALL_DAY, RH,
!cc     1    is1, is2, is3, pmax) 
!cc      IMPLICIT NONE
!cc
!cc
!cc
!ccc**    Cloudiness for the day
!cc      real(dp) VP_SAT, RAINFALL_DAY, RH, GAMMADIST_FN
!cc      real(dp) RNDF, xlo ,xhi, pmax, mock_rain,mu, gamma, zero
!cc      integer ::is1, is2, is3
!ccc
!cc      CLOUD_COVER_FN = 10-1.155*SQRT(VP_SAT/MAX(RAINFALL_DAY,0.1))
!ccc      CLOUD_COVER_FN = 10-2.5*SQRT(VP_SAT/MAX(RAINFALL_DAY,0.1))
!cc      CLOUD_COVER_FN = MAX(CLOUD_COVER_FN, 0.0)/10.0
!cc      if(RAINFALL_DAY.le.0.001) then
!cc        mu = rh  ! average rain on a wet day
!cc        gamma = 1.0
!cc        zero = 0.0
!cc        mock_rain = GAMMADIST_FN(zero, mu, gamma, is1, is2, is3)
!cc        CLOUD_COVER_FN = 10-2.5*SQRT(VP_SAT/MAX(mock_rain,0.1))
!ccc       CLOUD_COVER_FN = 10-1.155*SQRT(VP_SAT/MAX(mock_rain,0.1))
!cc
!cc        CLOUD_COVER_FN = MAX(CLOUD_COVER_FN, 0.0)/10.0
!cc        CLOUD_COVER_FN=1*CLOUD_COVER_FN
!cc      else  ! when raining adjust
!ccc     increase cloudiness - esp during summer when svp is high and raining
!cc        CLOUD_COVER_FN=1*CLOUD_COVER_FN
!cc        CLOUD_COVER_FN=min(1, (0.007*VP_SAT+0.995)*CLOUD_COVER_FN)
!cc      endif
!cc       
!ccc//! Equation: 49
!ccc//! Description: Proportion of the day that is cloudy \\
!ccc//!  INPUTS: saturated vapour-pressure; daily rainfall
!ccc//! Bibliography: S Evans {SWELTER); Nikolov & Zeller (1992) \
!ccc//!  Ecological Modelling 61:149-168
!ccc
!ccc//! E: Cld = 10-1.155 {(VP_{sat}\over{Rain})^0.5
!ccc
!ccc//! Variable: Cld
!ccc//! Description: Cloudiness of the day
!ccc//! Units: none  (proportion)
!ccc//! Type: REAL
!ccc//! SET
!ccc
!ccc//! Variable: VP_{sat}
!ccc//! Description: Saturated vapour pressure
!ccc//! Units: Pa
!ccc//! Type: REAL
!ccc
!ccc//! Variable: Rain
!ccc//! Description: amount of rainfall on the day
!ccc//! Units: mm
!ccc//! Type: REAL
!ccc
!cc      RETURN
!cc      END





!**********************************************************************!
integer function DAYS_IN_MONTH_FN(MM, DAYSINYEAR)
!**********************************************************************!
integer ::DAYSINYEAR, MM, LEAPSHIFT
!**   returns the day number of days in a month (allows for leap year).

	LEAPSHIFT=DAYSINYEAR-365
	if(MM==1) then                                                    
	  DAYS_IN_MONTH_FN = 31                                           
	else if(MM==2) then    
	  DAYS_IN_MONTH_FN = 28 + LEAPSHIFT                               
	else if(MM==3) then                                               
	  DAYS_IN_MONTH_FN = 31
	else if(MM==4) then                                                
	  DAYS_IN_MONTH_FN = 30
	else if(MM==5) then
	  DAYS_IN_MONTH_FN = 31
	else if(MM==6) then 
	  DAYS_IN_MONTH_FN = 30
	else if(MM==7) then     
	  DAYS_IN_MONTH_FN = 31
	else if(MM==8) then 
	  DAYS_IN_MONTH_FN = 31
	else if(MM==9) then 
	  DAYS_IN_MONTH_FN = 30
	else if(MM==10) then 
	  DAYS_IN_MONTH_FN = 31
	else if(MM==11) then 
	  DAYS_IN_MONTH_FN = 30
	else if(MM==12) then
	  DAYS_IN_MONTH_FN = 31
	endif

!
!//! Equation: 50
!//! Description: finds the number of days in a month \\
!//!   INPUTS: Month of year (numeric); Number of days in the year
!//! Bibliography: No Specific reference to this function.
!
!//! E:   None   
!   
!//! Variable: DinM      
!//! Description:   Number of days in the month 
!//! Units: none         
!//! SET
!//! Type: integer ::
!
!//! Variable: mm      
!//! Description:  month number 
!//! Units: none         
!//! Type: integer ::
!
!//! Variable: DinY     
!//! Description: Days in the year (allows for leap year)
!//! Units: none         
!//! Type: integer ::
!
!**********************************************************************!
end function DAYS_IN_MONTH_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function RADIAN_FROM_DEGREE_FN(DEGREES, PI)
!**********************************************************************!
!**   360 degrees is 2*PI radians:
real(dp) :: DEGREES, PI
	RADIAN_FROM_DEGREE_FN = (2*PI*DEGREES)/360.0
!
!//! Equation: 51
!//! Description: converts degrees to radians \\
!//!   INPUTS: Degrees (decimal); the mathematical constant, PI
!//! Bibliography: No Specific reference to this function.
!
!//! E:   RAD = 2 \pi DEG \over{260.0)
!   
!//! Variable: RAD  
!//! Description: returned value: radian
!//! Units: radian
!//! Type: REAL
!//! SET
!
!//! Variable: \pi  
!//! Description: mathematical constant, pi, 3.14115...
!//! Units: radian
!//! Type: REAL
!
!//! Variable: DEG  
!//! Description: number of degrees to be converted
!//! Units: Degrees
!//! Type: REAL
!
!**********************************************************************!
end function RADIAN_FROM_DEGREE_FN
!**********************************************************************!





!**********************************************************************!
subroutine NORTHEAST_TO_LATLONG_SUB(EASTING, NORTHING, LATITUDE, LONGITUDE, PI)
!**********************************************************************!
!** Uses the National Grid scale factors, true origins etc (after     
!** Airy (1830). These are not the same os the Irish National grid etc
!** though the formaule are the same (just different values)          
 real(dp) :: NORTHING, EASTING, LATITUDE, LONGITUDE, F0
 real(dp) :: P0, L0, N0, E0, A, B, PI
 real(dp) :: P0_RAD, L0_RAD, R, ETASQ, V, PX, PXNEW, ESQ
 real(dp) :: I7, I8, I9, I10, I11, I12, I12A , M,N
!   
!  set the constants for the OS36 UK national grid settings:
 A = 6377563.396 
 B = 6356256.910 
 F0 = 0.9996012717 
 P0 = 49.0     ! 49 deg N  True origin
 L0 = -2.0     !  2 deg W  True origin
 E0 = 400000.0
 N0 = -100000.0
! convert true origins to radians 
P0_RAD = P0*PI/180.0       
L0_RAD = L0*PI/180.0    
!
! Eccentricity squared:                                                
ESQ = (A**2 - B**2)/(A**2)                                        
! Intermediates:                                                      
N = (A-B)/(A+B)          
!                                                                     
! PX new is the estimate for the the Latitude
PX = ((NORTHING-N0)/(A*F0))+P0_RAD                             
! Initialise M for first iteration, then continue iterations until
! desired accuracy is met (1mm)
M = -99999999999999.0
do while ((NORTHING-N0-M)>=1)
	M = B*F0*( (1 +N +5/4.*N**2 +5/4.*N**3)*(PX-P0_RAD)-    &         
 (3*N +3*N**2 +21/8.*N**3) *SIN(PX-P0_RAD) *COS(PX+P0_RAD) +  &
 (15/8.*N**2 +15/8.*N**3)*SIN(2*(PX-P0_RAD))*COS(2*(PX+P0_RAD)) &
 -(35/24.*N**3) *SIN(3*(PX-P0_RAD)) *COS(3*(PX+P0_RAD))  )
!                                                                      
	PXNEW = (NORTHING-N0-M)/(A*F0) + PX                          
	PX = PXNEW                                                   
enddo
!                                                                        
V = A*F0*(1-ESQ*SIN(PX)**2)**(-0.5)                              
R = A*F0*(1-ESQ)*(1-ESQ*SIN(PX)**2)**(-1.5)                      
ETASQ = V/R -1.0
!                                                                        
I7 = TAN(PX)/(2*R*V)                                             
I8 = TAN(PX)/(24*R*V**3)*(5+3*TAN(PX)**2+ETASQ -9*(TAN(PX)**2)*ETASQ**2)                                                  
I9 = TAN(PX)/(720*R*V**5) *(61 +90*TAN(PX)**2 +45*(TAN(PX)**4))
I10 = 1/COS(PX)/V                                                
I11 = 1/COS(PX)/(6*V**3) * (V/R +2*TAN(PX)**2)                   
I12 = 1/COS(PX)/(120*V**5) *(5 +28*TAN(PX)**2 +24*TAN(PX)**4)    
I12A = 1/COS(PX)/(5040*V**7) *(61 +622*TAN(PX)**2 + 1320*TAN(PX)**4 +720*TAN(PX)**6)                          
LATITUDE = PX-I7*(EASTING-E0)**2 +I8*(EASTING-E0)**4-I9*(EASTING-E0)**6                                    
LONGITUDE = L0_RAD +I10*(EASTING-E0) - I11*(EASTING-E0)**3 + I12*(EASTING-E0)**5 - I12A*(EASTING-E0)**7             
!
!//! Equation: 52
!//! Description: Routine to take Easting and Northings to convert\
!//!  to 'conventional' Latitude and longitude (in radians). The \
!//!  contsnts used are specifically set for the UK OS36 system\
!//!  This routine is believed to be correct and is based on the OS\
!//!  manual - it should be noted that accuarcy is Unknown- the best\
!//!  guess is to about 20m. \\
!//!   INPUTS: Easting; Northing ; Latitude; Longitude; Pi\
!//! Bibliography: A Guide to coordinate systems in Great Britain,\
!//!  An introduction to mapping systems and the use of GPS datasets\
!//!  with Ordnance Survey mapping (1999), Ordnance Survey, Southampton.
!
!//! E: See the  reference - far to complicated!!
!
!//! Variable: Latitude
!//! Description: Returned Latitude value 
!//! Units:  Radians (North)
!//! Type: REAL
!//! SET
!
!//! Variable: Longitude
!//! Description: Returned Longitude value 
!//! Units:  Radians (East; for West, negative number returned)
!//! Type: REAL
!//! SET
!
!//! Variable: Easting
!//! Description: Easting coordinate of the OS National Grid
!//! Units: m 
!//! Type: REAL
!
!//! Variable: Northing
!//! Description: Northing coordinate of the OS National Grid
!//! Units: m 
!//! Type: REAL
!
!**********************************************************************!
end subroutine NORTHEAST_TO_LATLONG_SUB
!**********************************************************************!





!**********************************************************************!
real(dp) function DAYDATA_FROM_MONTHDATA_FN(MONTHDATA,DAY_OF_MONTH)
!**********************************************************************!
!** obtains a single element from and array containing a months worth
!** of data.
real(dp) :: MONTHDATA(31)
integer ::DAY_OF_MONTH
!
DAYDATA_FROM_MONTHDATA_FN = MONTHDATA(DAY_OF_MONTH)
!
!//! Equation: 53
!//! Description: Simple 'get' routine, to extract a single bit of\
!//!  data from an array with a months worth of data. \\
!//!  INPUTS: Array with the month's data; Day of the month
!//! Bibliography: none
!
!//! E: DayData = MonthData(dd)
!
!//! Variable: DayData
!//! Description: returned element of the array
!//! Units: Same as contents of MonthData
!//! Type: REAL
!//! SET
!
!//! Variable: MonthData
!//! Description: Array of size 31 with data for a month in it
!//! Units: whatever the data is
!//! Type: REAL
!
!//! Variable: dd
!//! Description: Day of the month - index for referencing the element\
!//!   in the array
!//! Units: None  (Day number of month)
!//! Type: REAL
!
!     RETURN
!**********************************************************************!
end function DAYDATA_FROM_MONTHDATA_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function HUMIDITY_FN(TEMP_AIR_AMPLITUDE)
!**********************************************************************!
!** humidity based on the temperature max-min amplitude
real(dp) :: TEMP_AIR_AMPLITUDE, RHFUN
!
RHFUN = (110.3-(4.68*TEMP_AIR_AMPLITUDE))
!       RHFUN = 77.4-3.289*TEMP_AIR_AMPLITUDE+27.357
if(RHFUN<=5) then
	HUMIDITY_FN = 5.0
else if(RHFUN>=100) then
	HUMIDITY_FN = 100.0
else
	HUMIDITY_FN = RHFUN
endif
!
!//! Equation 54
!//! Description: Returns an estimate of the daily humidity based \
!//!  on the difference between Tmax and Tmin (the amplitude) \\
!//!   INPUTS: amplitude os the temperature difference (Tmax - Tmin)
!//! Bibliography: S Evans (SWELTER)
!
!//! Case: RH_x <= 5
!//!   E: RH = 5.0
!//! Case: RH_x >=100
!//!   E: RH = 100
!//! Case: 5 < RH_x <100
!//!   E: RH = RH_x  \\
!//!  \\
!//! E: RH_x = (110.3-(4.68*\Delta T_{air}))
!
!//! Variable: RH
!//! Description: Relative humidity
!//! Units: %
!//! Type: REAL
!//! SET
!
!//! Variable: \Delta T_{air}
!//! Description: Amplitude of temperature variation during the day
!//! Units: \deg C
!//! Type: REAL
!
!**********************************************************************!
end function HUMIDITY_FN
!**********************************************************************!





!**********************************************************************!
subroutine GET_GRIDSQ_SUB(GRIDSIZE, LATITUDE, LONGITUDE, GRIDLAT, GRIDLONG, PI, LONG_OFFSET)
!**********************************************************************!
!** This routine gets a grid square reference point from latitude and 
!** longitude it assumes that the grid reference point is at the center
!** of the grid square, and therefore any location within the square
!** refers back to the central point. Incomming Latitude and LOngitude
!** is assumed to be in radians, and westerly longitude denoted by
!** use of a negative number. The parameter LONG_OFFSET is used to 
!** to a 360 degree basis if necessary. Ie. 1 deg West is expected as 
!** -1 deg (in radians). If LONG_OFFSET is zero, then the negative 
!** aspect is maintained, if it is 360, then 1 deg W becomes 359 deg (E)
!** At the moment this deals only with the Northern Hemisphere
real(dp) :: GRIDSIZE, LATITUDE, LONGITUDE, GRIDLAT, GRIDLONG
real(dp) :: PI, LONG_OFFSET, SIGN, LATDEG, LONGDEG
!
! convert lat * long to degrees (easier reference)
LATDEG = LATITUDE*180.0/PI
LONGDEG = LONGITUDE*180.0/PI
!
GRIDLAT = INT(LATDEG/GRIDSIZE)*GRIDSIZE+GRIDSIZE/2.0
SIGN = 1
if(LONGDEG<0) SIGN = -1
GRIDLONG = INT(LONGDEG/GRIDSIZE)*GRIDSIZE+GRIDSIZE/2.0*SIGN+LONG_OFFSET
!
!//! Equation: 55
!//! Description: from a given lat-long reference for the Northern  \
!//!  hemisphere, (negative numbers for Westerly longitude), this \
!//!  pinpoints a reference grid point, assuming uniform grid squares\
!//!  centered at the grid reference point. The value returned is in \
!//! degrees. a value offset can be used to deal with E/W longitude\
!//! in slightly different ways. (see description of offset variable)\\
!//!  INPUTS: Grid-square size; latitude of site (radians); \
!//!    longitude of site (radians);\
!//!   Grid-reference latitude [output]; Grid-reference longitude \
!//!   [output]; Constant, Pi, longitude offset
!//! Bibliography: T Randle (2000) - seemed like a good idea at the time
!
!//! E: GR_{lat} = INT(LAT/GR_{size}) GR_{size} + GR_{size)/2.0 \\
!//! \\
!//! E: GR_{long} = INT(LONG/GR_{size}) Gr_{size} + \
!//!      GR_{size) S/2.0 + LAT_{off}  \
!//!    S = 1 if LONG >=0;      Else S=-1 if LONG <0
!
!//! Variable: GR_{lat}
!//! Description: Grid-point center latitude value 
!//! Units: Degrees 
!//! Type: REAL
!//! SET
!
!//! Variable: GR_{long}
!//! Description: Grid-point center longitude value 
!//! Units: Degrees 
!//! Type: REAL
!//! SET
!
!//! Variable: LAT
!//! Description: Latitude of the site (Northern Hemishpere only) - in \
!//!  the above equations LAT has been transformed into degrees, though\
!//!  the input to the routine is radians
!//! Units: Degrees (in calculation); (radian in input)
!//! Type: REAL
!
!//! Variable: LONG
!//! Description: Longitude of the site (Northern Hemishpere only) - in\
!//!  the above equations LONG is transformed into degrees, though\
!//!  the input to the routine is radians
!//! Units: Degrees (in calculation); (radian in input)
!//! Type: REAL
!
!//! Variable: GR_{size}
!//! Description: Size of the grid-square
!//! Units: degrees
!//! Type: REAL
!
!//! Variable:  LAT_{off}
!//! Description: Offset value to deal with different formats of \
!//!  longitude. If Westerly values are needed as negative, then \
!//!  LAT_{off} should be 0; if degrees from east are required, then\
!//!  it should have the value 360
!//! Units: degrees
!//! Type: REAL
!
!**********************************************************************!
end subroutine GET_GRIDSQ_SUB
!**********************************************************************!





!**********************************************************************!
real(dp) function VP_UNSAT2_FN(VP_SAT, HUMIDITY)
!**********************************************************************!
!**   reverses the humidity function and derives the vapour pressure 
!**  (unsaturated) from the saturated VP (at the current temp) and 
!**   humidity
real(dp) :: VP_SAT, HUMIDITY
!
VP_UNSAT2_FN = HUMIDITY*VP_SAT/100.0
!
!//! Equation: 56
!//! Description: Reverses the humidity function and derives the \
!//!  vapour pressure unsaturated) from the saturated VP (at the \
!//!  current temp) and humidity \\
!//!   INPUTS: Saturated vapour pressure; humidity
!//! Bibliogrpahy: none
!
!//! E: VP_{unsat} = {RH VP_{sat}}\over {100.0}
!
!//! Variable: VP_{unsat} 
!//! Description: unsaturated vapour pressure
!//! Units: mbar (or same units as VP_{sat}
!//! Type: REAL
!//! SET
!
!//! Variable: VP_{sat} 
!//! Description: saturated vapour pressure
!//! Units: mbar  (or other eg kPa - will change VP_{unsat} units)
!//! Type: REAL
!
!//! Variable: RH
!//! Description: Relative humidity
!//! Units: %
!//! Type: REAL
!
!**********************************************************************!
end function VP_UNSAT2_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function RAD_TO_RADPAR_FN(RAD, RAD_PERC_PAR)
!**********************************************************************!
!** function to convert 'global' radiation to the radiation falling
!** within the photosynthetically active wavelengths; circa 375-800nm
!** the value for this may vary somewhere between 35% and 60% - a 
!** general approximation seems to be around 45%
real(dp) ::  RAD, RAD_PERC_PAR
!
RAD_TO_RADPAR_FN = RAD*RAD_PERC_PAR/100.0
!
!//! Equation: 57
!//! Description: Converts 'global' radiation to the radiation \
!//!  falling in the photosynthetically active wavelengths; circa \
!//!  375-800nm \\
!//!   INPUTS: Total light; percentage of total light that is \
!//!    photosynthetically active
!//! Bibliography: S Evans (SWELTER); Salisbury & Ross (1992), \
!//!  Plant Physiology, Wadsworth Publishing, California (pp613)
!
!//! E: RAD_{PAR} = RAD RAD_{perc, PAR}/100.0
!
!//! Variable: RAD_{PAR} 
!//! Description: Photosynthetically active radiation (W/m2)
!//! Units: W/m^2
!//! Type: REAL
!//! SET
!
!//! Variable: RAD
!//! Description: Total radiation (W/m2)
!//! Units: W/m^2
!//! Type: REAL
!
!//! Variable: RAD_{perc,PAR}
!//! Description: Percentage of total radiation that is in the PAR range
!//! Units: none (percentage)
!//! Type: REAL
!
!**********************************************************************!
end function RAD_TO_RADPAR_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function RADPAR_TO_RAD_FN(RADPAR, RAD_PERC_PAR)
!**********************************************************************!
!** function to convert radiation falling
!** within the photosynthetically active wavelengths; circa 375-800nm
!** to total radiation 
!** the value for this may vary somewhere between 35% and 60% - a 
!** general approximation seems to be around 45%
real(dp) :: RADPAR, RAD_PERC_PAR
!
RADPAR_TO_RAD_FN = RADPAR/RAD_PERC_PAR*100.0
!
!//! Equation: 58
!//! Description: Converts PAR radiation to the total radiation \\
!//!   INPUTS: photosynthetically active; percentage of total light \
!//!    that is photosynthetically active
!//! Bibliography: S Evans (SWELTER); Salisbury & Ross (1992), \
!//!  Plant Physiology, Wadsworth Publishing, California (pp613)
!
!//! E: RAD = RAD_{PAR} 100.0 / RAD_{perc, PAR}
!
!//! Variable: RAD
!//! Description: Total radiation (W/m2)
!//! Units: W/m^2
!//! Type: REAL
!//! SET
!
!//! Variable: RAD_{PAR}
!//! Description: Photosynthetically active radiation (W/m2)
!//! Units: W/m^2
!//! Type: REAL
!
!//! Variable: RAD_{perc,PAR}
!//! Description: Percentage of total radiation that is in the PAR range
!//! Units: none (percentage)
!//! Type: REAL
!
!**********************************************************************!
end function RADPAR_TO_RAD_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function RADPAR_TO_RADPPFD_FN(RADPAR,PPFDCONV)
!**********************************************************************!
!** function to convert PAR radiation to photosynthtic flux density
!** general approximation seems to be around: 1W (PAR) = 4.5 umol PAR
real(dp) :: RADPAR, PPFDCONV
!
RADPAR_TO_RADPPFD_FN = RADPAR*PPFDCONV
!
!//! Equation: 59
!//! Description: Converts PAR radiation to the PPFD PAR \\
!//!   INPUTS: PAR radiation; conversion of PAR light to PPFD 
!//! Bibliography: S Evans (SWELTER); Salisbury & Ross (1992), \
!//!  Plant Physiology, Wadsworth Publishing, California (pp613)
!
!//! E: RAD_{PPFD} = RAD_{PAR} PPFD_{conv}
!
!//! Variable: RAD_{PPFD} 
!//! Description: Photosynthetic flux density  (umol/m2)
!//! Units: umol (PAR) /m^2
!//! Type: REAL
!//! SET
!
!//! Variable: RAD_{PAR}
!//! Description: Photosynthetically active radiation (W/m2)
!//! Units: W/m^2
!//! Type: REAL
!
!//! Variable: PPFD_{conv}
!//! Description: Conversion factor for convertinq W (PAR) to umol (PAR)
!//! Units: umol(PAR)/m^2 / W(PAR)/m^2
!//! Type: REAL
!
!**********************************************************************!
end function RADPAR_TO_RADPPFD_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function RADPPFD_TO_RADPAR_FN(RADPPFD,PPFDCONV)
!**********************************************************************!
!** function to convert photosynthtic flux density to PAR radiation 
!** general approximation usess around: 1W (PAR) = 4.5 umol PAR
real(dp) :: RADPPFD, PPFDCONV
!
RADPPFD_TO_RADPAR_FN = RADPPFD/PPFDCONV
!
!//! Equation: 60
!//! Description: Converts PPFD PAR to PAR radiation\\
!//!   INPUTS: PPFD PAR ; conversion of PAR light to PPFD 
!//! Bibliography: S Evans (SWELTER); Salisbury & Ross (1992), \
!//!  Plant Physiology, Wadsworth Publishing, California (pp613)
!
!//! E: RAD_{PAR} = RAD_{PPFD} PPFD_{conv}
!
!//! Variable: RAD_{PAR}
!//! Description: Photosynthetically active radiation (W/m2)
!//! Units: W/m^2
!//! Type: REAL
!//! SET
!
!//! Variable: RAD_{PPFD} 
!//! Description: Photosynthetically active flux density
!//! Units: umol (PAR) /m^2
!//! Type: REAL
!
!//! Variable: PPFD_{conv}
!//! Description: Conversion factor for convertinq W (PAR) to umol (PAR)
!//! Units: umol(PAR)/m^2 / W(PAR)/m^2
!//! Type: REAL
!
!**********************************************************************!
end function RADPPFD_TO_RADPAR_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function RAD_DIFF_ABOVE_FN(DIFFPROP_DAY, RAD_TERRA_DAY, ALBEDO_DIFF)
!**********************************************************************!
!** works out the amount of terrestrial radiation that hits a canopy
!** that is diffuse - this is after albedo reflection has occurred.
real(dp) :: RAD_TERRA_DAY, ALBEDO_DIFF
real(dp) :: DIFFPROP_DAY
!
RAD_DIFF_ABOVE_FN = RAD_TERRA_DAY*DIFFPROP_DAY*(1-ALBEDO_DIFF)
!
!//! Equation: 61
!//! Description: The amount of terrestrial radiation above the \
!//!  canopy that is diffuse - after albedo reflection \\
!//!   INPUTS: proportion of terrestrial radiation that is diffuse;\
!//!    terrestrial radiation; Albedo (reflection) for diffuse light
!//! Bibliography: none
!
!//! E: RAD_{abv,D} = RAD_{terra} Diff_{prop} (1-ALB_{D})
!
!//! Variable: RAD_{abv,D}
!//! Description: Diffuse radiation above the canopy after albedo\
!//!  reflection
!//! Units: W/m^2/day
!//! Type: REAL
!//! SET
!
!//! Variable: RAD_{terra}
!//! Description: Terrestrial radiation 
!//! Units: W/m^2/day
!//! Type: REAL
!
!//! Variable: Diff_{prop}
!//! Description: proportion of Terrestrial radiation that is diffuse 
!//! Units: none (proportion)
!//! Type: REAL
!
!//! Variable: ALB_{D}
!//! Description: Albedo (reflection) coefficient for diffuse light 
!//! Units: none (proportion)
!//! Type: REAL
!
!**********************************************************************!
end function RAD_DIFF_ABOVE_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function RAD_BEAM_ABOVE_FN(DIFFPROP_DAY, RAD_TERRA_DAY, ALBEDO_BEAM)
!**********************************************************************!
!** works out the amount of terrestrial radiation that hits a canopy
!** that is direct (beam)- this is after albedo reflection has occurred.
real(dp) :: RAD_TERRA_DAY, ALBEDO_BEAM
real(dp) :: DIFFPROP_DAY
!
RAD_BEAM_ABOVE_FN = RAD_TERRA_DAY*(1-DIFFPROP_DAY)*(1-ALBEDO_BEAM)
!
!//! Equation: 62
!//! Description: The amount of terrestrial radiation above the \
!//!  canopy that is direct (beam) - after albedo reflection \\
!//!   INPUTS: proportion of terrestrial radiation that is diffuse;\
!//!    terrestrial radiation; Albedo (reflection) for beam light
!//! Bibliography: none
!
!//! E: RAD_{abv,B} = RAD_{terra} (1-Diff_{prop}) (1-ALB_{B})
!
!//! Variable: RAD_{abv,B}
!//! Description: Direct (beam) radiation above the canopy after albedo\
!//!  reflection
!//! Units: W/m^2/day
!//! Type: REAL
!//! SET
!
!//! Variable: RAD_{terra}
!//! Description: Terrestrial radiation 
!//! Units: W/m^2/day
!//! Type: REAL
!
!//! Variable: Diff_{prop}
!//! Description: proportion of Terrestrial radiation that is diffuse 
!//! Units: none (proportion)
!//! Type: REAL
!
!//! Variable: ALB_{B}
!//! Description: Albedo (reflection) coefficient for direct (beam)light 
!//! Units: none (proportion)
!//! Type: REAL
!
!**********************************************************************!
end function RAD_BEAM_ABOVE_FN
!**********************************************************************!





!**********************************************************************!
subroutine RAIN_DURATION2_SUB(RAIN, RAIN_DURATION, RAIN_INTENSITY, DAYS_IN_MONTH, IS1, IS2, IS3)
!**********************************************************************!
!** This routine estimates the duration and intensity of rain on each day 
!** for a month Accuracy is suspect.
real(dp) :: RAIN(31), RAIN_DURATION(31), RAIN_INTENSITY(31)
real(dp) :: RAIN_LIM(10), TIMESCALE, FREQUENCY, PROB
real(dp) :: XLO, XHI
integer :: DAYS_IN_MONTH, IDAY, CLASSNO, MAX_RAINTIME, I
integer :: IS1, IS2, IS3
!                        
RAIN_LIM(1) = 0    
RAIN_LIM(2) = 5    
RAIN_LIM(3) = 10   
RAIN_LIM(4) = 15   
RAIN_LIM(5) = 20   
RAIN_LIM(6) = 25   
RAIN_LIM(7) = 50   
RAIN_LIM(8) = 75   
RAIN_LIM(9) = 100  

!** find which rain class the current preciptation is in
do 10 IDAY = 1, DAYS_IN_MONTH
	if(RAIN(IDAY)>0) then
	  if(RAIN(IDAY)>RAIN_LIM(9)) then
	    CLASSNO = 9
	  else
	    do 20 I = 9,2,-1
	      if(RAIN(IDAY)<=RAIN_LIM(I)) CLASSNO = I
20          continue
	  endif
!**       The maximum duration of a rain instance is
	  MAX_RAINTIME=MIN(24,INT(.5*EXP(0.1386*RAIN_LIM(CLASSNO))+.5))
!
!**       The cumulative number of rainy days of a class type per century 
!**       is proportional to the rainfall duration. 
!**       (ie: Rainy days of this magnitude per century is timescale*duration)
!
	  TIMESCALE = 1.39*((RAIN_LIM(CLASSNO)/25.4+0.1)**(-3.55))*10
!
!**       thus the cumulative instances of rain of this magnitude (per century)is

	  FREQUENCY = TIMESCALE*MAX_RAINTIME
!
!**       the probability of rain falling at some stage within time, t, is
!         prob = timescale*t/frequency; thus if we get a random number between 
!**       0.025 and 1 (ignoring the chances of P<0.025; then the rain 
!**       duration (minutes) is
!
	  XLO = 0.025
	  XHI = 1.0
	  PROB = RNDF(XLO, XHI, IS1, IS2, IS3)
	  RAIN_DURATION(IDAY) = PROB/TIMESCALE*FREQUENCY*60
!
!**       the intensity (mm/hour) is
	  RAIN_INTENSITY(IDAY) = 60.0/RAIN_DURATION(IDAY)*RAIN(IDAY)
	else 
	  RAIN_DURATION(IDAY)=0.0
	  RAIN_INTENSITY(IDAY)=0.0
	endif
10    continue

!**********************************************************************!
end subroutine RAIN_DURATION2_SUB
!**********************************************************************!





!**********************************************************************!
real(dp) function VP_SAT2_FN(TEMP_C)
!**********************************************************************!
!**   Returns the saturated air pressure for a given temperature
real(dp) :: TEMP_C
!
if(TEMP_C>=0 ) then
	VP_SAT2_FN = 610.78*exp(17.269*TEMP_C/(TEMP_C+237.3))
else
	     VP_SAT2_FN = exp(-6140.4/(273+TEMP_C)+28.916)
	   endif
!
!//! Equation: 71
!//! Description: Returns saturated vapour pressure for a temperature\\
!//!   INPUTS: Temperature
!//! Bibliography: http://www.natmus.dk/cons/tp/atmcalc/atmoclc1.htm  (26.01.2001)
!
!//! Case: if T>=0
!//!  E:    VP_{sat} = 610.78*exp(T\over{(T+238.3)}*17.2694)   
!//! Case: t<0 
!//!  E:    VP_{sat} = exp(-6140.4\over{(273+T)}+28.916)       
!
!//! Variable: VP_{sat}
!//! Description: saturated vapour pressure at a given temperature\\
!//!  typical ranges are 300-1200 Pa: 0.3-1.2kpa: 3-12mbar
!//! Units: Pa
!//! SET
!//! Type: REAL
!
!//! Variable: T
!//! Description: Air temperature 
!//! Units: Degrees  C 
!//! Type: REAL
!
!**********************************************************************!
end function VP_SAT2_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function TRANS_PROP_DAY2_FN(TEMP_AIR_AMPLITUDE,A_SKY, ANGSTROM_BETA, C_SKY)
!**********************************************************************!
!** determines aproximate transmittance of radiatio to the earth
real(dp) :: TEMP_AIR_AMPLITUDE
real(dp) :: A_SKY, ANGSTROM_BETA, C_SKY

if (TEMP_AIR_AMPLITUDE<=0) then
	      TRANS_PROP_DAY2_FN = ANGSTROM_BETA
else
	      TRANS_PROP_DAY2_FN = ANGSTROM_BETA*(1-EXP(-A_SKY*TEMP_AIR_AMPLITUDE**C_SKY))
endif

!
!//! Equation: 72
!//! Description: Works out the transmittance of radiation through the\
!//!   atmosphere - a reverse engineered approximation of eqn 41. Untested
!//!   and certainly dubious at low temperature amplitudes as we fix it to\
!//!   be Angstrom_beta; ehreas the original allows for transmittance to be\
!//!   >  Angstrom_beta\\
!//!   INPUTS: Air temperature amplitude (max-min);    \
!//!    Coefficient for maximum clear-sky transmittance;  \
!//!   Angstrom \beta value; Coefficient of clear-sky transmittance \
!//!   with amplitude temperature increase. 
!//! Bibliography: S Evans (SWELTER); Bristow & Campbell (1984) Agric \
!//!  & Forest Met 31: 159-166
!
!//! Case: \Delta_T <=0.0
!//!  E: T_{sw} = \beta
!//! Case: \Detla_T > 0
!//!  E: T_{sw} = \beta(1-exp(-A_{sky}*\Delta_T^{C_{sky}})) 
!
!//! Variable: T_{sw}
!//! Description: proportion of Extra-terrestrial radiation that \
!//!  arrives on the earth's surface as short-wave.
!//! Units: none
!//! Type: REAL
!//! SET
!
!//! Variable: \Delta T
!//! Description: Daily amplitude of temperature variation
!//! Units: Degrees C  (usually)
!//! Type: REAL
!
!//! Variable: \beta
!//! Description: Angstrom tubidity factor: max clear sky transmittance
!//! Units: none  
!//! Type: REAL
!
!//! Variable: A_{sky}
!//! Description: Coefficient of maximum clear-sky transmittance \
!//!  characteristics
!//! Units: none  
!//! Type: REAL
!
!//! Variable: C_{sky}
!//! Description: Coefficient of maximum clear-sky transmittance \
!//!   with \Delta T increase
!//! Units: none  
!//! Type: REAL
!
!**********************************************************************!
end function TRANS_PROP_DAY2_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function GAMMADIST_FN(A,B,C, SEED1, SEED2, SEED3)
!**********************************************************************!
real(dp) :: Aconst, Bconst, Cconst, Qconst, Tconst
real(dp) :: Dconst, a, b, c, xlo, xhi, R1, R2, P, Y
real(dp) :: retval, p1, p2, v, z, w, R3
integer :: Flag
integer ::seed1, seed2, seed3
Bconst = c-log(4.0)
Tconst = 4.5
Dconst = 1+log(Tconst)
Cconst = 1+c/exp(1.0)
flag = 1
Xlo = 0.0
xhi = 1.0

if(c<1) then
	  flag =0
100     R1 = RNDF(XLO, XHI, SEED1, SEED2, SEED3)
	  p = Cconst*R1
	  if(P>1) goto 300
  y=p**(1/c)
	  R2 = RNDF(XLO, XHI, SEED1, SEED2, SEED3)
	  if(R2<=exp(-y)) then
	    retval = a+b*y
	    flag = 2
	  endif
	  if(flag==0) goto 100

300     y = -log((Cconst-p)/c)
	  r2 = RNDF(XLO, XHI, SEED1, SEED2, SEED3)
	  if((y>0).and.(R2<=y**(c-1))) then
	    retval = a+b*y
	    flag = 3
  endif        
  if(flag==0) goto 100

elseif(c>1) then
  Aconst = 1/SQRT(2*c-1.0)  
  Qconst = c+ 1/(Aconst)
	  flag=0 
20      p1 = RNDF(XLO, XHI, SEED1, SEED2, SEED3)
	  p2 = RNDF(XLO, XHI, SEED1, SEED2, SEED3)
	  v = Aconst*log(p1/(1-p1))
	  y = c*exp(v)
	  z = p1*p1*p2
	  w = Bconst+Qconst*v - y
	  if(w+Dconst-Tconst*z>=0) then
	    retval =a+b*y
	    flag = 3
	  else
	    if(w>=log(z)) then
	      retval = a+b*y
	      flag = 4
	    endif
	  endif
	  if(flag==0) goto 20
else
  r3 = RNDF(XLO,XHI, SEED1, SEED2, SEED3)
  if(b>0) then
    retval = a-b*log(r3)
  else
    write(*,*) 'problem with gamma dist when c=1'
    stop
  endif
endif

gammadist_fn = retval 
!//! Equation: 73
!//! Description: Derives a random value from the gamma distribution \\
!//!   INPUTS: 3 parameters describing the gamma distribution shape; 3\
!//!    Seeds for random number generation
!//! Bibliography: R Saucier (2000), Computer generation of statistical\ 
!//!  Distributions; (US) Army Research Laboratory, ARL-TR-2168.\
!//!  (http://ftp.arl.mil/random/  (17/5/2001) )
!
!**********************************************************************!
end function GAMMADIST_FN
!**********************************************************************!




     
!**********************************************************************!
real(dp) function NORM_FROM_RAND_FN(MEAN, SD, SEED1,SEED2, SEED3, NORMSWITCH)
!**********************************************************************!
real(dp) :: Mean, SD, NormVal
integer :: Seed1, Seed2, Seed3, NormSwitch
!//! Equation: 86
!//! Description: Generates a value from a pseudo normal distribution \\
!//!   is a shell for calling different methods. \\
!//!   INPUTS: Mean; Sd; 3 Seeds for random number generation; Method id
!//! Bibliography: None
!
if(NormSwitch==1) then
   NormVal = NORMDIST_FN(SEED1,SEED2,SEED3)
!         IF((MEAN.NE.Zero).and.(SD.NE.1)) THEN
   if ((abs(MEAN)>0.0).and.(abs(SD-1.0)>0.0)) then
!          transform the distribution
     NormVal = TRANSNORM_FN(Mean, SD, NormVal)
   endif
else if(NormSwitch==2) then
   NormVal = NORMDIST2_FN(MEAN,SD,SEED1,SEED2,SEED3)
endif
NORM_FROM_RAND_FN = NormVal

!**********************************************************************!
end function NORM_FROM_RAND_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function NORMDIST_FN(SEED1,SEED2,SEED3)
!**********************************************************************!
real(dp) :: xlo, xhi, R1, retval
integer ::seed1, seed2, seed3
 
Xlo = 0.0
xhi = 1.0
R1 = RNDF(XLO, XHI, SEED1, SEED2, SEED3)

Retval = ( R1**0.135 - (1-R1)**0.135 )/0.1975

NORMDIST_FN = Retval
!//! Equation: 74
!//! Description: Derives a value from standard normal distribution\
!//!  (mean 0, sd 1), at random. Reputed to have accuracy of 0.5% in\
!//!  The range 0<=p<=0.9975 (NB random numbers are between 0 and 1) \\
!//!   INPUTS: 3 Seeds for random number generation
!//! Bibliography: S Evans, SWELTER; Haith, DA, Tubbs, LJ and \
!//!  Pickering, NB(1984) Simulation of Pollution by soil erosion and \
!//!  soil nutrient loss; Pudoc, Wageningen.
!
!//! E: N = { R1^{0.135} - (1-R1)^{0.135)} } \over {0.1975}
! 
!//! Variable: SEED1, SEED2, SEED3
!//! Description: integer ::seeds for random number generation
!//! Units: none  
!//! Type: integer
!
!//! Variable: R1
!//! Description: Random number between 0 and 1
!//! Units: none  
!//! Type: REAL
!
!//! Variable: NORMDIST_FN
!//! Description: Simulated value from a normal distribution N(0,1)
!//! Units: none
!//! Type: REAL
!//! SET
!
!**********************************************************************!
end function NORMDIST_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function NORMDIST2_FN(MEAN,SD,SEED1,SEED2,SEED3)
!**********************************************************************!
real(dp) :: MEAN, SD
real(dp) :: xlo, xhi, R1, R2, P, retval
integer ::seed1, seed2, seed3
!
Xlo = -1.0
xhi = 1.0
P = 2.0
!  
10    if (P>=1.0) then 
	R1 = RNDF(XLO, XHI, SEED1, SEED2, SEED3)
	R2 = RNDF(XLO, XHI, SEED1, SEED2, SEED3)
	P = R1*R1 + R2*R2
goto 10
endif

!
Retval = MEAN + SD*R1*SQRT(-2.0*LOG(P)/P)
!     NB . . Natural Log!
!
NORMDIST2_FN = Retval
!
!//! Equation: 75
!//! Description: Derives a random value from the Normal distribution\\
!//!   INPUTS: 2 parameters describing the distribution shape;- Mean \
!//!     And Standard deviation; 3 Seeds for random number generation
!//! Bibliography: R Saucier (2000), Computer generation of statistical\
!//!   Distributions; (US) Army Research Laboratory, ARL-TR-2168; \
!//!   (http://ftp.arl.mil/random/  (17/5/2001) )
!
!//! E: N(\mu, \sigma) = \mu + \sigma R1 \root{(-2 Ln(P)/P)}  \
!//!     Where P = R1^2 + R2^2
! 
!//! Variable: SEED1, SEED2, SEED3
!//! Description: integer ::seeds for random number generation
!//! Units: none  
!//! Type: integer
!
!//! Variable: R1, R2
!//! Description: Random number between -1 and +1
!//! Units: none  
!//! Type: REAL
!
!//! Variable: NORMDIST2_FN
!//! Description: Simulated value from a normal distribution \
!//!    N(\mu,\sigma)
!//! Units: none
!//! Type: REAL
!//! SET
!
!**********************************************************************!
end function NORMDIST2_FN
!**********************************************************************!



	
!      real(dp) FUNCTION LAT_SLOPE_FN(Incl_rad, Lat_rad,
!     1   Aspect_rad, Pi)
!      IMPLICIT NONE

!      real(dp) Lat_rad, Incl_rad, retval, Pi, Aspect_rad
!
!      retval = COS(Incl_rad)*SIN(Lat_rad)+
!     1   SIN(Incl_rad)*COS(Lat_rad)*COS(Aspect_rad) 
!      if(retval.lt.-1) retval = -1.0
!      if(retval.gt.1) retval = 1.0 
!      retval=ASIN(retval)
!
!//! Equation 76
!//! Description: calculates an equivalent latitude for a point, \
!//!    given a slope and inclination\\
!//!   INPUTS: Inclination (from horizontal); Latitude; Aspect; Pi
!//! Bibliography: Swift, L.W. (1976) Algorithm for Solar Radiation on \
!//!   Mountain Slopes Water Resources Research (12) 1: 108-112
!
!//! E: Lat_{slope} = ASIN(COS(Incl_{rad}*SIN(Lat_{rad}+ \
!//!          COS(Lat_{rad})*SIN(Asp_{rad}))
!
!//! Variable:  Incl_{rad}
!//! description: Inclination of the slope (from horizontal)
!//! Units: Raidans 
!//! Type: REAL
!
!//! Variable:  Lat_{rad}
!//! description: latitude of the site
!//! Units: Raidans 
!//! Type: REAL
!
!//! Variable:  Asp_{rad}
!//! description: Aspect (Azimuth) of the slope (from North)
!//! Units: Raidans 
!//! Type: REAL
!
!//! Variable: Pi
!//! Description: mathematical constant, pi
!//! Units: none
!//! Type: REAL
!
!//! Variable: Lat_{slope}
!//! Description: Latitude of an equivalent point, accounting for \
!//!    the slope & aspect
!//! Units: radians
!//! Type: REAL
!//! SET
!
!      LAT_SLOPE_FN = retval
!      RETURN
!      END




!**********************************************************************!
real(dp) function RAD_SLOPE_FN(Declin_rad, Lat_rad, Lat_slope_rad, Inclin_rad, Aspect_rad, Pi)
!**********************************************************************!
real(dp) :: Declin_rad, Lat_rad, Lat_slope_rad, Inclin_rad
real(dp) :: Aspect_rad, Pi
real(dp) :: Time_Offset, T, T6, T7, T1, T0, T3, T2, T8, T9
real(dp) :: Rad0, Rad1, Rad_horiz, Adj

!    calculate an offset for time for the actual and equivalent slopes
Time_Offset = ATAN( SIN(Inclin_rad)*SIN(Aspect_rad) / ( COS(Inclin_rad)*COS(Lat_rad) - &
 SIN(Inclin_rad)*SIN(Lat_rad)*COS(Aspect_rad)) )*Pi/180
!    calculate new 'adjusted' dawn/dusk hour angles
T = DAWN_HOURANGLE_FN(Declin_rad, Lat_slope_rad)
T7 = T-Time_Offset
T6 = -T-Time_Offset
!    claclulate un-adjusted hour-angles
T = DAWN_HOURANGLE_FN(Declin_rad, Lat_rad)
T1 = T
T0 = -T
if(T7<T1) then
	T3 = T1
else
	T3 = T7
endif
if(T6>T0) then
	T2 = T6
else
	T2 = T0
endif
!
!    This bit supposedly guards against double sunrises and 90 degree latitudes 
!   UNTESTED!!!
!    Calculate the potential on horizontal surface
if(T2>=T3) then
	T2 = 0
	T3 = 0
endif
T6 = T6+Pi
Rad0 = 0.0
if(T6>=T1) then
	T8 = T6
	T9 = T1
else
	T7 = T7 - Pi
	if(T7>T0) then
	  T8 = T0
	  T9 = T7
	  Rad0 = SIN(Declin_rad)*SIN(Lat_slope_rad)*(T9-T8)*12/Pi + COS(declin_rad)*COS(Lat_slope_rad)* &
 ( SIN(T9+Time_Offset)-SIN(T8+Time_Offset))*12/Pi
	endif
endif
Rad1 = SIN(Declin_rad)*SIN(Lat_slope_rad)*(T3-T2)*12/Pi + COS(declin_rad)*COS(Lat_slope_rad)* &
 ( SIN(T3+ Time_Offset)-SIN(T2+ Time_Offset))*12/Pi
Rad1= Rad1+Rad0
!
!    Calculate the Potential on a horizontal surface
Adj = 0.0
Rad_Horiz = SIN(Declin_rad)*SIN(Lat_rad)*(T1-T0)*12/Pi + COS(declin_rad)*COS(Lat_rad)* &
 ( SIN(T1+Adj)-SIN(T0+Adj))*12/Pi
RAD_SLOPE_FN = Rad1/Rad_Horiz
 
!//! Equation 77
!//! Description: Calculates the ratio between horizontal surface \
!//!   Direct radiation and on a slope it finds an equivalent latitude \
!//!   for daylength purposes, but uses the correct declinations and \
!//!   elevations. It appears that the day is still spread evenly about \
!//!   noon! - this may account for the non-night light that occurs before \
!//!   the 'new sunrise/set
!//! Bibliography: Swift, L.W. (1976) Algorithm for Solar Radiation on \
!//!   Mountain Slopes. Water Resources Research (12) 1: 108-112
!//! E:  Not even going to try - see the paper!!
!
!**********************************************************************!
end function RAD_SLOPE_FN
!**********************************************************************!





!**********************************************************************!
subroutine RAD_SLOPE_DAY_SUB(RAD_BEAM_SLOPE, RAD_DIFF_SLOPE, &
!**********************************************************************!
 RAD_TERRA_DAY, DIFFPROP_DAY, SLOPE_RAD, &
 ASPECT_RAD, DECLINATION_RAD, LAT_RAD, SUNRISE, PI) 
	   implicit none

real(dp) :: RAD_BEAM_SLOPE, RAD_DIFF_SLOPE
real(dp) :: RAD_TERRA_DAY,DIFFPROP_DAY,SLOPE_RAD,ASPECT_RAD
real(dp) :: DECLINATION_RAD, lAT_RAD, SUNRISE, PI
real(dp) :: CT1, SHF, SMP, W, CT2, CT3, CT, CTS, CTZ, ST
real(dp) :: CTZS, TFR, DIFFR
integer ::SH,i
!
!     calculates the adjustment for beam and diffuse light intercepted on
!     a sloping surface. No account is taken of the difference of length
!     of the surface compared to the flat plane!! It also takes no account
!     of increased shading beyond the slope (if facing away from the sun)-
!     we know nothing about surrounding terrain.
	  
CT1 = SIN(DECLINATION_RAD)*(SIN(LAT_RAD)*COS(SLOPE_RAD)- &
 COS(LAT_RAD)*SIN(SLOPE_RAD)*COS(ASPECT_RAD))
SH = INT(SUNRISE+1)      
SHF = SH - SUNRISE
SMP = SUNRISE+SHF/2.0
W = 15*(SMP-12)*PI/180.0
CT2 = COS(DECLINATION_RAD)*COS(W)*(COS(LAT_RAD)*COS(SLOPE_RAD)+ &
 SIN(LAT_RAD)*SIN(SLOPE_RAD)*COS(ASPECT_RAD))
CT3 = COS(DECLINATION_RAD)*SIN(SLOPE_RAD)*SIN(ASPECT_RAD)*SIN(W)
CT = CT1+CT2+CT3
CTS = CT*SHF
CTZ = COS(LAT_RAD)*COS(DECLINATION_RAD)*COS(W)+ &
 SIN(LAT_RAD)*SIN(DECLINATION_RAD)
CTZS = CTZ*SHF
TFR = CT/CTZ

!      DO 10 ST = SH+0.5, 11.5
do i=0,11
  ST = SH+i+0.5
  W = 15*(ST-12)*PI/180.0
  CT2 = COS(DECLINATION_RAD)*COS(W)*(COS(LAT_RAD)*COS(SLOPE_RAD)+ &
 SIN(LAT_RAD)*SIN(SLOPE_RAD)*COS(ASPECT_RAD))
  CT3 = COS(DECLINATION_RAD)*SIN(SLOPE_RAD)*SIN(ASPECT_RAD)*SIN(W)
  CT = CT1+CT2+CT3
  CTS = CTS+CT
  CTZ = COS(LAT_RAD)*COS(DECLINATION_RAD)*COS(W) &
 +SIN(LAT_RAD)*SIN(DECLINATION_RAD)
  CTZS = CTZS+CTZ
  tfr = ct/ctz
  
enddo
!10    CONTINUE
TFR = CTS/CTZS
!
!     now diffuse
DiffR = COS(SLOPE_RAD/2)*COS(SLOPE_RAD/2) 
     
RAD_BEAM_SLOPE = TFR*RAD_TERRA_DAY*(1-DIFFPROP_DAY)
RAD_DIFF_SLOPE = DiffR*(RAD_TERRA_DAY*DIFFPROP_DAY)

!//! Equation 77
!//! Description: modifies incomming light to account for a slope and aspect\
!//!  effect. Sin in the northern hemispere, north facing slopes may \
!//!  only get diffuse radiation, the light must first be split into \
!//!  direct (beam) and diffuse light. No account is taken of the change in\
!//!  surface length, or of any shading beyond the current slope. Hence \
!//!  total radiation over many surfaces will be greater than over a single,\
!//!  flat plane. The difference is aprox 1/cos(slope).
!//!  Bibliography:Duffie and Beckman (1991). Solar engineering of thermal\
!//!    processes. Wiley. Also See groups.google.com; alt.energy.renewable;\
!//!    Subject: Re: Irradiation; Posted 5/1/1998; accessed 15/6/2001. \
!//!    Diffuse radiation from Montieth J.L. (1973) Principles of \
!//!    Environmantal Physics, Arnold, London.
!//! E: see the references

!**********************************************************************!
end subroutine RAD_SLOPE_DAY_SUB
!**********************************************************************!





!**********************************************************************!
real(dp) function TEMP_ALT_CORR_FN(TEMP, METALT, ALTITUDE)
!**********************************************************************!

!**  Function corrects temperature for altitude - should be applied
!**  directly after the temperaure calculation as will therefore
!**  impact on other functions eg SVP etc
real(dp) ::  TEMP, METALT, ALTITUDE, ALTDIFF_FT
!
 ALTDIFF_FT = (ALTITUDE-METALT)*100/2.54/12.0
 TEMP_ALT_CORR_FN = TEMP-(ALTDIFF_FT)/1000.0*2.0                                 
!//! Equation: 78
!//! Description: corrects temperature for site altitude. correction \
!//!   is approc 2C per 1000 ft \\
!//!   INPUTS: Temperature; altitude of met station (m); \
!//!      Altitude of site (m)
!//! Bibliography:  ??
!//!  
!//!  E: T_{alt} = T-ALT_{diff}*2
!
!//! Variable: T_{alt}
!//! Description: corrected air temperature \\
!//! Units: C
!//! Type: REAL
!//! SET
!
!//! Variable: T
!//! Description: Temperature at met station elevation
!//! Units: C
!//! Type: REAL
!
!//! Variable: ALT_{diff}
!//! Description: Altitude difference between site and met station
!//! Units: Feet
!//! Type: REAL
!
!**********************************************************************!
end function TEMP_ALT_CORR_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function PPT_RAIN_FN(TEMP_DAY_MEAN, RAIN_DAY)
!**********************************************************************!

real(dp) :: ZERO_SNOW_T, ALL_SNOW_T, TEMP_DAY_MEAN, RAIN_DAY
real(dp) :: RAINPART
	ZERO_SNOW_T=2.0
ALL_SNOW_T=-2.0

	if(TEMP_DAY_MEAN>ZERO_SNOW_T) then
   PPT_RAIN_FN = RAIN_DAY*1.0
else if(TEMP_DAY_MEAN<ALL_SNOW_T) then
   PPT_RAIN_FN = 0.0
else
!       some combination of snow and rain
  RAINPART = (TEMP_DAY_MEAN - ALL_SNOW_T)/(ZERO_SNOW_T-ALL_SNOW_T)
  PPT_RAIN_FN = RAINPART*RAIN_DAY
endif
  
!
!//! Equation: 79
!//! Description: Finds the amount of real(dp) rain from Gross rain \
!//!  ie. that which is not snow \\ 
!//!  INPUTS: Mean Temperature of the day; Daily rainfall.
   
!//! Bibliography: S Evans; SWELTER. eq.67
!
!//! Case: T_{mean} > T_{snow}
!//!      E:  PPT_{rain} = PPT_{day}
!//! Case: T_{mean} < T_{all_snow}
!//!      E: PPT_{rain} = 0.0
!//! Case: T_{all_snow} < T_{mean} < T_{snow}
!//!      E: PPT_{rain} = (T_{mean} - T_{all_snow})/    \
!//!           (T_{snow}-T_{all_snow}) * PPT_{day}
!   
!//! Variable: T_{Mean}
!//! Description: Mean daily temperature
!//! Units: Celcius
!//! Type: REAL
!
!//! Variable: PPT_{day}
!//! Description: Gross precipitation
!//! Units: mm/day
!//! Type: REAL
!    
!//! Variable: T_{snow} 
!//! Description: Temperature at which snow may start (above which is all rain)        
!//! Units: Celcius
!//! Type: REAL
!
!//! Variable: T_{all_snow}
!//! Description: Temperature below which is all snow 
!//! Units: Celcius
!//! Type: real(dp)
!
!//! Variable: PPT_{rain}
!//! Description: Amount of precipitation which is real(dp) rain
!//! Units: mm/day  (rain)          
!//! Type: REAL
!//! SET   
!
!**********************************************************************!
end function PPT_RAIN_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function PPT_SNOW_FN(RAIN_DAY, PPT_RAIN)
!**********************************************************************!

real(dp) :: RAIN_DAY, PPT_RAIN
   
	PPT_SNOW_FN = (RAIN_DAY - PPT_RAIN)*1.1811
  
!//! Equation: 80
!//! Description: Finds the amount of snow from Gross rain \\
!//!  INPUTS: Daily rainfall; Precipitation which is rain
   
!//! Bibliography: S Evans; SWELTER. eq.67
!
!//!      E: PPT_{snow} = (Rain{day} - PPT_{Rain})*1.1811
!   
!//! Variable: PPT_{day}
!//! Description: Gross precipitation
!//! Units: mm/day
!//! Type: REAL
!
!//! Variable: PPT_{Rain}
!//! Desciption: Precipitation which is Rain
!//! Units: mm/day  (rain)
!    
!//! Variable: PPT_{snow}
!//! Description: precipitation as snow (snow depth)
!//! Units: mm/day  (snow)          
!//! Type: REAL
!//! SET   
!
!**********************************************************************!
end function PPT_SNOW_FN
!**********************************************************************!





!**********************************************************************!
real(dp) function transnorm_fn(MU, SD, Xval)
!**********************************************************************!

real(dp) :: MU, SD, Xval
!     Transforms a standard normal distribution, N(0,1))
!     to N(Mu, SD)
!
!//! Equation: 81
!//! Description: Transforms a value from a Standard Normal distribution N(0,1) to\
!//!   a 'generic' normal N(mu, sd)\\
!//!   INPUTS: Mean of Distribution; Sd of distribution; Standard Normal value
!//! Bibliography: None (T Houston Pers Comm)
TRANSNORM_FN = Xval*SD+Mu
!
!//! E: N(\mu, \sigma) = N(0,1)*\sigma + \mu
!
!//! variable: \mu
!//! Type: REAL
!//! Description: Mean of the desired normal distribution
!//! Units: Any
!
!//! Variable: \sigma
!//! Type: REAL
!//! Description: Standard deviation of the desired normal distribution
!//! Units: Any
!
!//! Variable: N(0,1)
!//! Type: Double Precsion
!//! Description: Value from the standard normal distribution
!//! Units: Any
!
!//! Variable: N(\mu, \sigma)
!//! Type: REAL
!//! Description: Transformed Normal Value - from the shifted norma distribution
!//! Units: Any
!//! SET
!
!**********************************************************************!
end function transnorm_fn
!**********************************************************************!





!**********************************************************************!
real(dp) function normskew_fn(NormVal, Lambda)
!**********************************************************************!

real(dp) :: NormVal, Lambda
!     Box_cox Normaility transformation of a Normal Distribution
!
!//! Equation: 82
!//! Description: Skews a normal disribution\\
!//!   INPUTS: Value from a normal distribution N(Mu,Sd);\
!//!     And Response variable
!//! Bibliography: C. Coarkin & P Tobias, NIST/SEMATECH Engineering Handbook;\
!//!   Ch. 1.3.3.6\
!//!  http://www.itl.nist.gov/div898/handbook/eda/section3/boxcox.htm (9/8/2001)
!
!//! E: SN(N(\mu,\sigma)) = (N(\mu,\sigma)^\lambda -1)\over{\lambda}
!  
NORMSKEW_FN = ( (Normval**lambda)-1)/lambda
!
!//! Variable: \\lambda
!//! Type: REAL
!//! Description: Response variable - if \lambda= 0; use log of data
!//! Units: none
!
!//! Variable: N(\mu,\sigma)
!//! Type: REAL
!//! Description: value from the normal distribution
!//! Units: Any
!//! Variable: SN(N(\mu,\sigma)) 
!//! Type: REAL
!//! Description: box-cox skewed, (transformed) normal value
!//! Units: Any
!//! SET
!
!**********************************************************************!
end function normskew_fn
!**********************************************************************!





!**********************************************************************!
real(dp) function TAMPFROMRH_FN(RH)
!**********************************************************************!

!** humidity based on the temperature max-min amplitude
real(dp) :: RH

TAMPFROMRH_FN = (RH-110.3_dp)/(-4.68_dp)
!
!//! Equation 83
!//! Description: Returns an estimate of the daily air temperature\
!//!  amplitude based on RH. Inverse of RH equation\\
!//!   INPUTS: Relative humidity (%)
!//! Bibliography: S Evans (SWELTER)
!
!//! E: \Delta T_{air} = (RH-110.3)\over{(-4.68)}
!
!//! Variable: RH
!//! Description: Relative humidity
!//! Units: %
!//! Type: REAL
!
!//! Variable: \Delta T_{air}
!//! Description: Amplitude of temperature variation during the day
!//! Units: \deg C
!//! Type: REAL
!//! SET
!
!**********************************************************************!
end function TAMPFROMRH_FN
!**********************************************************************!

end module metdos
                         
