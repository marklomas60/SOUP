      REAL*8 FUNCTION DECLINATION_FN(DAY_OF_YEAR,DAYS_IN_YEAR,
     1   PI)
      IMPLICIT NONE
c

      REAL*8 MAXDECLINATION, PI
      INTEGER*4 DAYS_IN_YEAR, DAY_OF_YEAR
	MAXDECLINATION = 23.45
	DECLINATION_FN =MAXDECLINATION * PI/180.0 *
     1    SIN(PI/180.0 *(360D0*DBLE(DAY_OF_YEAR)/
     2      DBLE(DAYS_IN_YEAR)-80D0))
c
c//! Equation: 1
c//! Description: Approximates declination angle for each day: RESULT \
c//!  is in RADIANS. Assume sine wave, with dec=23.45 on 21 June. \\ 
c//!  INPUTS: Day-number of the year; Number of days in the year;\
c//!  mathematical constant, PI.          
c//! Bibliography: No Specific reference to this function.
c
c//! E:  Decl = D_{max} * \pi/180 * SIN(\pi/180*(360*DoY\over{DinY}-80))
c   
c//! Variable: Decl 
c//! Description: returned value of the suns's declination
c//! Units: Radians
c//! SET
c//! Type: REAL*8
c
c//! Variable: D_{max}
c//! Description: maximum value of the suns's declination
c//! Value: 23.45
c//! Units: Degrees
c//! Type: REAL*8
c    
c//! Variable: \pi 
c//! Description: Mathematical constant, pi 3.1415...        
c//! Units: Radian
c//! Type: REAL*8
c
c//! Variable: DoY
c//! Description: Day of the year,(1st jan = 1, 31 dec=365 or 366) 
c//! Units: none          
c//! Type: Integer           
c
c//! Variable: DinY
c//! Description: Total number of Days in year,(leap year=366)
c//! Units: none          
c//! Type: Integer           
c
      RETURN                                           
      END                                              


      REAL*8 FUNCTION RADIAN_TO_DEGREE_FN()
      IMPLICIT NONE
c**   360 degrees is 2*PI radians:

	RADIAN_TO_DEGREE_FN = 180.0/(2.0*ASIN(1.0))
c
c//! Equation: 2
c//! Description: a value for the number of degrees in 1 radian \\
c//!   INPUTS: none
c//! Bibliography: No Specific reference to this function.
c
c//! E:   R2D = 180.0 \over{(2.0*ASIN(1.0)} 
c   
c//! Variable: R2D  
c//! Description: returned value: degrees in 1 radian
c//! Units: Degree/radian
c//! SET
c//! Type: REAL*8
c
      RETURN
      END


      REAL*8 FUNCTION DEGREE_TO_RADIAN_FN()
      IMPLICIT NONE
c**   360 degrees is 2*PI radians
c**   1 degree is approx 0.01745 radians:

	DEGREE_TO_RADIAN_FN = 2.0*ASIN(1.0)/180.0
c
c//! Equation: 3
c//! Description: a value for the number of radians in 1 degree  \\
c//!   INPUTS: none
c//! Bibliography: No Specific reference to this function.
c
c//! E:   D2R = 2.0*ASIN(1.0) \over{180.0} 
c   
c//! Variable: D2R  
c//! Description: returned value: radians in 1 degree
c//! Units: Radian/Degree
c//! SET
c//! Type: REAL*8
c
      RETURN
      END


      REAL*8 FUNCTION PIVALUE_FN()
      IMPLICIT NONE
c**   Returns the value PI (3.14157...)              

	PIVALUE_FN = 2.0*ASIN(1.0)          
c
c//! Equation: 4
c//! Description: Sets the value of Pi \\
c//!   INPUTS: none
c//! Bibliography: No Specific reference to this function.
c
c//! E:   Pi = 2.0*ASIN(1.0) 
c   
c//! Variable: Pi   
c//! Description:  mathematical constant number Pi, 3.14159...
c//! Units: none         
c//! SET
c//! Type: REAL*8
c
      RETURN
      END


      REAL*8 FUNCTION DAYLENGTH_FN(LATITUDE_RAD, DECLINATION,
     1  PI)
      IMPLICIT NONE
c**   calculate the day-length, use declination and latitude            

      REAL*8 DECLINATION, LATITUDE_RAD
      REAL*8 PI, COSHS, HS
      REAL*8 ALTANGLE
	ALTANGLE = -0.833*PI/180.0
	COSHS = -((SIN(LATITUDE_RAD)*SIN(DECLINATION)-SIN(ALTANGLE)) /
     1    (COS(LATITUDE_RAD)*COS(DECLINATION)))     
	HS = ACOS(COSHS)*180/PI                                        
	DAYLENGTH_FN = HS*2/15.0
c
c//! Equation: 5
c//! Description: calculates the daylength in hours. No Account is  \
c//!   taken of altitude, or longitude.  \\
c//!   INPUTS: Latitude; Solar declination; mathematical constant, PI
c//! Bibliography: No Specific reference to this function.
c
c//! E:   daylen=180*2*ACOS( (SIN(Lat)*SIN(Decl)-SIN(altang)) \over  \
c//!       {(COS(Lat)*COS(Decl))} ) /Pi /15.0  \\
c//!       altang = -0.833 \pi/180.0
c   
c//! Variable: daylen
c//! Description:  appproximation of daylength
c//! Units: hours        
c//! SET
c//! Type: REAL*8
c
c//! Variable: Lat
c//! Description:  Latitude of site (in radians) 
c//! Units: radians      
c//! Type: REAL*8
c
c//! Variable: Decl
c//! Description:  declination of sun (radians) 
c//! Units: radians      
c//! Type: REAL*8
c
c//! Variable: altang
c//! Description: altitude angle (radians); calculated value at sea-level.
c//! Value: -0.14538592....
c//! Units: radians      
c//! Type: REAL*8
c
c//! Variable: Pi
c//! Description:  mathematical constant number Pi, 3.14159...
c//! Units: none      
c//! Type: REAL*8
c
      RETURN                                                         
      END                                                            


      REAL*8 FUNCTION DAWN_FN(DAYLENGTH)
      IMPLICIT NONE

c**   calculate the time of dawn: use declination and latitude to
c**   get daylength and go from there value is local solar time.
      REAL*8 DAYLENGTH
	DAWN_FN = 12.0-(DAYLENGTH/2.0)
c
c//! Equation: 6
c//! Description: finds the approximate time of dawn - no account is \
c//!   taken of altitude, longitude etc. \\
c//!   INPUTS: Daylength
c//! Bibliography: No Specific reference to this function.
c
c//! E:   T_{dawn} = 12.0 - (daylen/2.0)
c   
c//! Variable: T_{dawn} 
c//! Description:  Approximate time of dawn.
c//! Units: Local Solar time
c//! SET
c//! Type: REAL*8
c
c//! Variable: daylen   
c//! Description:  Approximate daylength.   
c//! Units: none         
c//! Type: REAL*8
c
      RETURN                                                         
      END                                                            


      REAL*8 FUNCTION DUSK_FN(DAYLENGTH)
      IMPLICIT NONE

c**   calculate the time of dusk: use declination and latitude            
      REAL*8 DAYLENGTH
	DUSK_FN = 12.0+(DAYLENGTH/2.0)                                         
c
c//! Equation: 7
c//! Description: finds the approximate time of dusk - no account is \
c//!   taken of altitude, longitude etc. \\
c//!   INPUTS: Daylength
c//! Bibliography: No Specific reference to this function.
c
c//! E:   T_{dusk} = 12.0 - (daylen/2.0)
c   
c//! Variable: T_{dusk} 
c//! Description:  Approximate time of dusk.
c//! Units: Local Solar Time
c//! SET
c//! Type: REAL*8
c
c//! Variable: daylen   
c//! Description:  Approximate daylength.   
c//! Units: none         
c//! Type: REAL*8
c
      RETURN                                                         
      END                                                            


      INTEGER*4 FUNCTION DAY_OF_YEAR_FN(DD,MM,DAYSINYEAR)
      IMPLICIT NONE

      INTEGER*4 DAYSINYEAR, LEAPSHIFT
      INTEGER*4 DD, MM, JUL
c**   returns the day of the year given the month day and year.
c**   checks for leap year using DAYS_IN_YEAR
c**   
c**   These are NOT JULIAN dates (although people often refer mistakenly  
c**   as such).                                                           
	LEAPSHIFT = 0
	IF(DAYSINYEAR.EQ.366) LEAPSHIFT=1
	IF(MM.EQ.1) THEN                                                    
	  JUL = DD                                                        
	ELSE IF(MM.EQ.2) THEN    
	  JUL = 31+DD                                                     
	ELSE IF(MM.EQ.3) THEN                                               
	  JUL = 31+28+DD + LEAPSHIFT
	ELSE IF(MM.EQ.4) THEN                                                
	  JUL = 31+28+31+DD + LEAPSHIFT
	ELSE IF(MM.EQ.5) THEN
	  JUL = 31+28+31+30+DD + LEAPSHIFT
	ELSE IF(MM.EQ.6) THEN 
	  JUL = 31+28+31+30+31+DD + LEAPSHIFT
	ELSE IF(MM.EQ.7) THEN     
	  JUL = 31+28+31+30+31+30+DD + LEAPSHIFT
	ELSE IF(MM.EQ.8) THEN 
	  JUL = 31+28+31+30+31+30+31+DD + LEAPSHIFT
	ELSE IF(MM.EQ.9) THEN 
	  JUL = 31+28+31+30+31+30+31+31+DD + LEAPSHIFT
	ELSE IF(MM.EQ.10) THEN 
	  JUL = 31+28+31+30+31+30+31+31+30+DD + LEAPSHIFT
	ELSE IF(MM.EQ.11) THEN 
	  JUL = 31+28+31+30+31+30+31+31+30+31+DD + LEAPSHIFT
	ELSE IF(MM.EQ.12) THEN
	  JUL = 31+28+31+30+31+30+31+31+30+31+30+DD + LEAPSHIFT
	ENDIF 
	DAY_OF_YEAR_FN = JUL
c
c//! Equation: 8
c//! Description: finds the day of the year, starting with Jan 1 =1 \\
c//!   INPUTS: Day of Month; Month of year (numeric); Days in Year
c//! Bibliography: No Specific reference to this function.
c
c//! E:   DoY = \sum_{i=0}^{mon-1} Dm_{i} + dd
c   
c//! Variable: DoY      
c//! Description:   Day of the year 
c//! Units: none         
c//! SET
c//! Type: Integer          
c
c//! Variable: mon      
c//! Description:  month number 
c//! Units: none         
c//! Type: Integer         
c
c//! Variable: Dm_{i}
c//! Description:  Days in month (i) 
c//! Units: none         
c//! Type: Integer         
c
c//! Variable: dd       
c//! Description:   day of the current month
c//! Units: none         
c//! Type: Integer         
c
      RETURN 
      END 


      INTEGER*4 FUNCTION DAYS_IN_YEAR_FN(YYYY)
      IMPLICIT NONE

c**   checks the year and tests if it is a leap-year. Returns the
c**   number of days in the year (ie. 365 or 366). 
      INTEGER*4 YYYY, MODVAL
	MODVAL=100
	IF( MOD(YYYY,MODVAL).EQ.0) THEN
c**       century - must be divisible by 400
	  MODVAL=400
	  IF( MOD(YYYY,MODVAL).EQ.0) THEN
	    DAYS_IN_YEAR_FN=366
	  ELSE
	    DAYS_IN_YEAR_FN=365
	  ENDIF
	ELSE
c**     not a century - is it divisble by 4
	  MODVAL=4
	  IF( MOD(YYYY,MODVAL).EQ.0) THEN
	    DAYS_IN_YEAR_FN=366
	  ELSE
	    DAYS_IN_YEAR_FN=365
	  ENDIF
	ENDIF
      RETURN
c
c//! Equation: 9
c//! Description:  calculates number of days in any year  \\
c//!   INPUTS: Year
c//! Bibliography:  After: L. Halsall, Stats & Computing, Forest \
c//!   Research, Forestry Commission, Alice Holt, Farnham.
c
c//! Case: MOD(year, 100)=0 and MOD(year,400)=0
c//! E:   DinY = 366
c//! Case: MOD(year,100)<>0 and MOD(year,4)=0
c//! E:   DinY = 366
c//! Case: All other conditions
c//! E:   DinY = 365
c   
c//! Variable: DinY     
c//! Description: Number of days in a year:    \\
c//!   If it is a century and divisible by 400, then it is a leap year \\
c//!   If its NOT a century but divisible by 4, then it is a leap year
c//! Units: none         
c//! SET
c//! Type: Integer          
c
c//! Variable: year     
c//! Description:   full year number (eg 1997, 2002 etc).
c//! Units: none         
c//! Type: Integer         
c
      END



      REAL*8 FUNCTION DECLINATION2_FN(DAY_OF_YEAR)
      IMPLICIT NONE

c
      REAL*8 MAXDECLINATION
      INTEGER*4 DAY_OF_YEAR
	MAXDECLINATION = 0.39785
	DECLINATION2_FN = MAXDECLINATION*SIN(4.868961+0.017203*    
     1    DAY_OF_YEAR + 0.033446*SIN(6.224111 + 0.017202*DAY_OF_YEAR))
c
c//! Equation: 15
c//! Description: Approximates declination angle for each day: RESULT \
c//!  is in RADIANS.                                       \\
c//!   INPUTS: Day-number of year
c//! Bibliography:  S. Evans: SWELTER 
c
c//! E:   Decl = D_{max} * sin(4.868961+0.017203 DoY + 0.033446      \
c//!             sin(6.22411 + 0.017202 DoY))
c   
c//! Variable: Decl 
c//! Description: returned value of the suns's declination
c//! Units: Radians
c//! SET
c//! Type: REAL*8
c
c//! Variable: D_{max}
c//! Description: maximum value of the suns's declination
c//! Value: 0.39785
c//! Units: Radians
c//! Type: REAL*8
c    
c//! Variable: DoY
c//! Description: Day of the year,(1st jan = 1, 31 dec=365 or 366) 
c//! Units: none          
c//! Type: Integer           
c
      RETURN                                           
      END                                              


      REAL*8 FUNCTION DAYLENGTH2_FN(LATITUDE_RAD, DAY_OF_YEAR,
     1  PI)
      IMPLICIT NONE

c**   calculate the day-length, use declination and latitude            
      REAL*8 LATITUDE_RAD
      REAL*8 PI
      REAL*8 AMPLITUDE 
      INTEGER*4 DAY_OF_YEAR
	AMPLITUDE = exp((7.42+0.045*LATITUDE_RAD*180./PI))/3600.0
	DAYLENGTH2_FN = 12.0 + AMPLITUDE*SIN((DAY_OF_YEAR-79)*0.017121)
c
c//! Equation: 16
c//! Description: calculates the daylength in hours from latitude \
c//! and day-number. no Account is taken of altitude, or longitude. \\ 
c//!   INPUTS: Latitude; Day-number of year, Mathematical constant, PI
c//! Bibliography: S. Evans: SWELTER 
c
c//! E:   daylen= 12 + Amp sin((DoY-79) 0.017121)   \\ 
c//!        Amp={exp(7.42+0.045 lat 180/\pi )} \over{3600}
c   
c//! Variable: daylen
c//! Description:  appproximation of daylength. based on solar time
c//! Units: hours        
c//! SET
c//! Type: REAL*8
c
c//! Variable: lat
c//! Description:  Latitude of site (in radians) 
c//! Units: radians      
c//! Type: REAL*8
c
c//! Variable: Amp
c//! Description: amplitude of sine curve at relevant latitude. 
c//! Units: none      
c//! Type: REAL*8
c
c//! Variable: Pi
c//! Description:  mathematical constant number Pi, 3.14159...
c//! Units: none      
c//! Type: REAL*8
c
c//! Variable: DoY 
c//! Description: Day of the year,(1st jan = 1, 31 dec=365 or 366) 
c//! Units: none          
c//! Type: Integer           
c
      RETURN                                                         
      END                                                            


      REAL*8 FUNCTION DAWN_HOURANGLE_FN(DECLINATION, 
     1   LATITUDE_RAD)
      IMPLICIT NONE

c**   find the sun-rise/sun-set solar hour-angle: ie at DAWN and DUSK
      REAL*8 DECLINATION, LATITUDE_RAD, Temp
        Temp = -TAN(LATITUDE_RAD)*TAN(DECLINATION) 
        If(Temp.lt.-1) Temp = -1.0
        if(Temp.gt.1) Temp =1.0
	DAWN_HOURANGLE_FN = ACOS(Temp)
c
c//! Equation: 17
c//! Description: calculates the sunrise/sunset hour angle for the \
c//!  day (declination) and latitude. Answer is in Radians \\
c//!   INPUTS: Solar declination; Latitude
c//! Bibliography:  S. Evans: SWELTER 
c
c//! E:   hs = acos(-tan(Lat)*tan(Decl))
c   
c//! Variable: hs   
c//! Description: returned value of the solar sunrise/set hour angle \
c//!  uncorrected for longitude.
c//! Units: Radians
c//! SET
c//! Type: REAL*8
c
c//! Variable: Decl
c//! Description:  declination of sun (radians) 
c//! Units: radians      
c//! Type: REAL*8
c
c//! Variable: Lat
c//! Description:  Latitude of site (in radians) 
c//! Units: radians      
c//! Type: REAL*8
c
      RETURN
      END

      REAL*8 FUNCTION RAD_ET_DAY2_FN(LATITUDE_RAD,DECLINATION,
     1   DAWN_HOURANGLE, SOL_ELLIPSE, SOL_CONST, PI)
      IMPLICIT NONE

c**   caclulates the extra-terrestrial daily radiation 
      REAL*8 DECLINATION, LATITUDE_RAD
      REAL*8 DAWN_HOURANGLE, SOL_CONST, PI
      REAL*8  SOL_ELLIPSE 
c
      RAD_ET_DAY2_FN = SOL_CONST*SOL_ELLIPSE*(
     1  DAWN_HOURANGLE*SIN(LATITUDE_RAD)*SIN(DECLINATION) + 
     2  COS(LATITUDE_RAD)*COS(DECLINATION)*SIN(DAWN_HOURANGLE) )
     3  / PI*86400
c
c//! Equation: 18
c//! Description: caclulates the daily extra-terrestrial radiation  \
c//!    result in W/m2/day         \\
c//!   INPUTS: Latitude; Solar Declination; Hour-angle of sun-rise/set;\
c//!    Elliptical scalar for suns' orbit; Solar constant, Mathematical\
c//!    constant, Pi
c//! Bibliography:  S. Evans: SWELTER 
c
c//! E:   RAD_{et}= S_c E (hs sin(Lat) sin(Decl) +        \
c//!       cos(lat) cos(decl) sin(hs)) \over{\pi} 86400
c   
c//! Variable: RAD_{et}
c//! Description: daily value of solar radiation outside the atmosphere
c//! Units: W/m2/day
c//! SET
c//! Type: REAL*8
c
c//! Variable: Decl
c//! Description:  declination of sun (radians) 
c//! Units: radians      
c//! Type: REAL*8
c
c//! Variable: Lat
c//! Description:  Latitude of site (in radians) 
c//! Units: radians      
c//! Type: REAL*8
c
c//! Variable: S_c
c//! Description:  Solar constant. Nb there may be slight variations \
c//!    on this depending on the source (approx 1367)
c//! Units: W/m2
c//! Type: REAL*8
c
c//! Variable: Pi   
c//! Description:  mathematical constant Pi: 3.14159...
c//! Units: none          
c//! Type: REAL*8  
c
c//! Variable: hs   
c//! Description: value of the solar sunrise/set hour angle.
c//! Units: Radians
c//! Type: REAL*8
c
c//! Variable: E
c//! Description: Solar elliptical orbit scalar
c//! Units: none          
c//! Type:  REAL*8 
c
      RETURN
      END

      REAL*8 FUNCTION RAD_ET_DAY_FN(LATITUDE_RAD, 
     1  DECLINATION, DAYLENGTH, SOL_ELLIPSE,
     2  SOL_CONST, PI)
      IMPLICIT NONE


c**   caclulates the extra-terrestiral daily radiation 
      REAL*8 DECLINATION, LATITUDE_RAD
      REAL*8 SOL_CONST, PI
      REAL*8 DAYLENGTH, NOON, SOL_ELEV_FN
      REAL*8 SOL_ELEV_NOON, SOL_ELLIPSE, SOL_NOON

c  since the routine sol_elev has local time as an input; convert 
c  noon solar time to gmt
      NOON = 12.0     ! LST
      SOL_ELEV_NOON = SOL_ELEV_FN(LATITUDE_RAD, 
     1   DECLINATION, NOON, PI)
      SOL_NOON = SOL_ELLIPSE*SOL_CONST*SIN(SOL_ELEV_NOON)
      RAD_ET_DAY_FN = SOL_NOON*2*DAYLENGTH/(PI*SIN(PI/2.0) )         
     1     *60*60!  /1000/1000 
c
c//! Equation: 19
c//! Description: caclulates the daily extra-terrestrial radiation  \
c//!    result in W/m2/day                        \\
c//!   INPUTS: Latitude; Longitude; Solar declination; Daylength; \
c//!     Solar elliptical orbit scalar; Day-number of the year;\
c//!     Solar constant; Mathematical constant, pi.
c//! Bibliography: none/FOREST-GROWTH
c
c//! E:   RAD_{et}= S_{12} 2 Day_{l} / (\pi sin(\pi/2)                   \\
c//!       S_{12} = E S_c sin(Sol_el)                                 \\
c//!       Sol_{el} = \pi/2 - Lat + Decl                              \\
c//!       E =  1.033-0.066*\sqrt{(1-((DoY-183)**2)/(183.**2) ) }
c   
c//! Variable: RAD_{et}
c//! Description: daily value of solar radiation outside the atmosphere
c//! Units: W/m2/s
c//! SET
c//! Type: REAL*8
c
c//! Variable: E   
c//! Description:  Elliptical scalar for the orbit of the sun 
c//! Units: none         
c//! Type: REAL*8
c
c//! Variable: Decl
c//! Description:  declination of sun (radians) 
c//! Units: radians      
c//! Type: REAL*8
c
c//! Variable: Lat
c//! Description:  Latitude of site (in radians) 
c//! Units: radians      
c//! Type: REAL*8
c
c//! Variable: S_c
c//! Description:  Solar constant. Nb there may be slight variations \
c//!    on this depending on the source  (eg 1367.0)
c//! Units: W/m2         
c//! Type: REAL*8
c
c//! Variable: Pi   
c//! Description:  mathematical constant Pi: 3.14159...
c//! Units: none          
c//! Type: REAL*8  
c
c//! Variable: Day_l
c//! Description:  day length
c//! Units: hours         
c//! Type: REAL*8  
c
c//! Variable: DoY 
c//! Description: Day of the year,(1st jan = 1, 31 dec=365 or 366) 
c//! Units: none          
c//! Type: Integer           
c
      RETURN
      END

      REAL*8 FUNCTION RAD_ET_DAY3_FN(LATITUDE_RAD, 
     1  DECLINATION, SOL_ELLIPSE, DAWN, DUSK,
     2  SOL_CONST, PI)
      IMPLICIT NONE


c**   caclulates the extra-terrestrial daily radiation 
      REAL*8 DECLINATION, LATITUDE_RAD
      REAL*8 SOL_CONST, PI, SOL_ELLIPSE
      REAL*8 SOLAR_TIME, SOL_ELEV_FN
      REAL*8 SOL_NOW, SIN_SOL_ELEV, DAWN, DUSK
      INTEGER i

      RAD_ET_DAY3_FN = 0.0
      DO i=0,288
        SOLAR_TIME = real(i)/12.0
*      DO SOLAR_TIME = 0,24, 1/12.0
c**     step throughout the day in 5 min intervals
	IF((SOLAR_TIME.GE.DAWN).AND.(SOLAR_TIME.LT.DUSK)) THEN
	  SIN_SOL_ELEV = SIN( SOL_ELEV_FN(LATITUDE_RAD,  
     1     DECLINATION, SOLAR_TIME, PI) )
	  SOL_NOW = MAX(0.0d0,SOL_ELLIPSE*SOL_CONST*SIN_SOL_ELEV)*60*5
	  RAD_ET_DAY3_FN = RAD_ET_DAY3_FN + SOL_NOW
	ENDIF
      ENDDO
c
c//! Equation: 20
c//! Description: caclulates the daily extra-terrestrial radiation  \
c//!   by 5 minute interation using solar elevation result in W/m2/day. \
c//!   This functions call a solar elevation fuunction internally.  \\
c//!   INPUTS: Latitude; Longitude; Solar Declination; scalar of\
c//!   sun's elliptical orbit; Day-number of the year; Solar constant;\
c//!   Mathematical constant, pi
c//! Bibliography: C.J.T Spitters et al. (1986); Agric & Forest Met \
c//!   38:217-229
c
c//! E:  RAD_{et}= \sum_{t=0}^{24*12} max(0, S_c E sin(Sol_{el}) 300) \\
c//!     if(t>=dawn and t<= dusk)
c   
c//! Variable: RAD_{et}
c//! Description: daily value of solar radiation outside the atmosphere
c//! Units: W/m2/s
c//! SET
c//! Type: REAL*8
c
c//! Variable: Sol_{el}
c//! Description:  Solar Elevation; returned through SOL_ELEV_FN
c//! Units: radians      
c//! Type: REAL*8
c
c//! Variable: Decl
c//! Description:  declination of sun (radians) 
c//! Units: radians      
c//! Type: REAL*8
c
c//! Variable:  E 
c//! Description:   Solar elliptical scalar     
c//! Units:  none        
c//! Type: REAL*8
c
c//! Variable: Lat
c//! Description:  Latitude of site (in radians) 
c//! Units: radians      
c//! Type: REAL*8
c
c//! Variable: Dawn
c//! Description:  Time of Sunrise (Local Solar time)
c//! Units: radians      
c//! Type: REAL*8
c
c//! Variable: Dusk
c//! Description:  Time of Sunset (Local Solar time)
c//! Units: radians      
c//! Type: REAL*8
c
c//! Variable: S_c
c//! Description:  Solar constant. Nb there may be slight variations \
c//!    on this depending on the source  (1367.0)
c//! Units: W/m2         
c//! Type: REAL*8
c
c//! Variable: Pi   
c//! Description:  mathematical constant Pi: 3.14159...
c//! Units: none          
c//! Type: REAL*8
c
c//! Variable: DoY 
c//! Description: Day of the year,(1st jan = 1, 31 dec=365 or 366) 
c//! Units: none          
c//! Type: Integer           
c
      RETURN
      END

      REAL*8 FUNCTION VP_SAT_FN(TEMP_C)
      IMPLICIT NONE

c**   Returns the saturated air pressure for a given temperature
      REAL*8 TEMP_C
c
      VP_SAT_FN = 6.1078*exp(17.269*TEMP_C/(TEMP_C+237.3))
c
c//! Equation: 21
c//! Description: Returns saturated vapour pressure for a temperature\\
c//!   INPUTS: Temperature
c//! Bibliography: S Evans: SWELTER, Groff-Gratch equation
c
c//!  E:    VP_{sat} = 6.1078 exp{( {17.269 T}\over{(T+237.3)} )}
c
c//! Variable: VP_{sat}
c//! Description: saturated vapour pressure at a given temperature
c//! Units: mbar
c//! SET
c//! Type: REAL*8
c
c//! Variable: T
c//! Description: Air temperature 
c//! Units: Degrees  C 
c//! Type: REAL*8
c
       RETURN
       END


      REAL*8 FUNCTION SOL_ELEV_FN(LATITUDE_RAD,
     1    DECLINATION, SOLAR_TIME, PI)
      IMPLICIT NONE

c      function returns the elevation (radians), deviation of the sun from
c      the horizontal to the vertical. Also referred to as the solar 
c      altitude angle.
										
      REAL*8 LATITUDE_RAD
      REAL*8 DECLINATION, H, PI
      REAL*8 H1, SOLAR_TIME
									
      H1 = 15.*(SOLAR_TIME-12.0)       ! degrees   
      H = H1/180*PI  ! radians                                      
      SOL_ELEV_FN = ASIN( COS(LATITUDE_RAD)*COS(DECLINATION)*COS(H)
     1   + SIN(LATITUDE_RAD)*SIN(DECLINATION) )
c//! Equation: 22
c//! Description: Returns solar elevation angle at time of the day \
c//!    and day of the year.                                      \\
c//!  INPUTS:  Latitude; Longitude; Declination; \\
c//!   Local Solar time; mathematical constant, pi
c//! Bibliography: C GronBeck, 1999; Sun Angle (web pages)
c
c//!  E:   SOL_{el} = asin(sin(Lat) sin(Decl) + cos(Lat) cos(Decl) \
c//!     cos(H)                                                       \\
c//!     H = 15(Lst-12) \pi /180.0 
c
c//! Variable: SOL_{el}
c//! Description:  Solar elevation of the sun
c//! Units: Radians     
c//! SET
c//! Type: REAL*8
c
c//! Variable: Lat
c//! Description: Latitude of site 
c//! Units:  Radians
c//! Type: REAL*8
c
c//! Variable: Decl
c//! Description:  declination of sun (radians) 
c//! Units: radians      
c//! Type: REAL*8
c
c//! Variable: Lst
c//! Description: Local solar time
c//! Units: (hour) 
c//! Type: REAL*8
c
c//! Variable: \pi
c//! Description: Mathematical constant, pi
c//! Units: (radians)
c//! Type: REAL*8 
c
       RETURN
       END

										
      REAL*8 FUNCTION AZIMUTH_FN(LATITUDE_RAD, DECLINATION,
     1   SOL_ELEV)
      IMPLICIT NONE

c     function returns the azimuth angle (in radians): Angular distance
c     between SOUTH and the projection of the line of sight to the
c     sun on the ground.
										
      REAL*8 LATITUDE_RAD
      REAL*8 DECLINATION, SOL_ELEV, SIGN, HOURANG

c check the sign of the hour-angle. Firstly convert back to hour-angle
c    using latitude, declination and elevation.
      HOURANG = ACOS(                                       
     1   (SIN(SOL_ELEV)-SIN(LATITUDE_RAD)*SIN(DECLINATION))/ 
     2   (COS(LATITUDE_RAD)*COS(DECLINATION) ) )
      SIGN = 1.0
      IF(HOURANG.LT.0) SIGN=-1.0
									
      AZIMUTH_FN = SIGN*ACOS( (SIN(SOL_ELEV)*SIN(LATITUDE_RAD) -
     1    SIN(DECLINATION)) / (COS(SOL_ELEV)*COS(LATITUDE_RAD) ) )
c//! Equation: 23
c//! Description: Returns azimuth angle at any time of the day and  \
c//!     day of the year. \\
c//!  INPUTS: Latitude; solar declination; Solar elevation; \
c//!  Mathematical constant, pi
c//! Bibliography: C GronBeck, 1999; Sun Angle (web pages)
c
c//!  E:   Az = acos( {sin(Sol_{el}) sin(Lat) - sin(Decl)}   \
c//!     /over{cos(Sol_{el}) cos(Decl)}                      \\
c//!     and: the sign of Az is the same as the sign of the hour-angle.
c
c//! Variable: Az    
c//! Description:  Azimuth angle of the sun. 
c//! Units: Radians     
c//! SET
c//! Type: REAL*8
c
c//! Variable: Lat
c//! Description: Latitude of site 
c//! Units:  Radians
c//! Type: REAL*8
c
c//! Variable:  Decl
c//! Description: solar declination
c//! Units:  Radians
c//! Type: REAL*8
c
c//! Variable: Sol_{el}       
c//! Description: Solar elevation 
c//! Units:  Radians
c//! Type: REAL*8 
c
c//! Variable: \pi
c//! Description: Mathematical constant, pi
c//! Units: (radians)
c//! Type: REAL*8 
c
      RETURN
      END


      REAL*8 FUNCTION RAD_TERRA_DAY_FN(RAD_ET_DAY, 
     1   LATITUDE_RAD, CLOUD_COVER, BETA)
      IMPLICIT NONE

c     converts extra-terrestiral radiation to surface radiation.
      REAL*8 RAD_ET_DAY, LATITUDE_RAD 
      REAL*8 CLOUD_COVER
      REAL*8 ALPHA, BETA , SIGMA
c
c**     Angstrom tubidity factor related to aerosol size & properties
      ALPHA = 18.0-64.884*(1-1.3614*COS(LATITUDE_RAD))
      SIGMA = 0.03259
      Alpha = alpha*4.1842*100*100  ! convert from cal/cm to w/m2
c      another fudge to try and transmit less light on heavy cloud days
c      if(cloud_cover.gt.1.0) alpha = alpha*2
c       if(cloud_cover.lt.0.15) then
c         cloud_cover = cloud_cover/1.4
c         alpha = 0.0
c       endif
       sigma = .02*LOG(max(cloud_cover,.001d0))+Sigma
c      print*, 'in rad terra',rad_et_day,beta,sigma,cloud_cover, alpha
      RAD_TERRA_DAY_FN = RAD_ET_DAY*(BETA-SIGMA*CLOUD_COVER*10)-ALPHA
c
c//! Equation: 24
c//! Description: Converts extra-terrestrial radiation to terrestrial \
c//!   radiation  uncorrected for site altitude.                      \\ 
c//!   INPUTS: Extra-terrestrial radiation; Latitude; cloud-cover; \
c//!     Angstrom (beta) turbidity
c//! Bibliography: S. Evans:SWELTER; Nikolov & Zeller (1992)
c
c//! E:   RAD_{terra} = RAD_{et} (\beta -\sigma 10 Cld) - \alpha   \\
c//!        \alpha = 18.0-64.884*(1-1.3614*cos(Lat))     \\
c//!        \sigma = 0.03259   
c
c//! Variable: RAD_{terra}
c//! Description: Terrestrial radiation on a horizonatal plane \
c//!  for the given latitude and day.
c//! Units: W/m2/d
c//! Type: REAL*8
c//! SET
c
c//! Variable: RAD_{et}
c//! Description: Extra-Terrestrial radiation above the atmosphere  \
c//!  for the given latitude and day.
c//! Units: W/m2/d
c//! Type: REAL*8
c
c//! Variable: \alpha  
c//! Description: Angstrom tubidity factor related to aerosol size & \
c//!  properties
c//! Units: none  
c//! Type: REAL*8
c
c//! Variable: \beta
c//! Description: Angstrom tubidity factor: max clear sky transmittance
c//! Units: none  
c//! Type: REAL*8
c
c//! Variable: \sigma  
c//! Description: Angstrom tubidity factor: Absorption by cloud cover
c//! Units: none  
c//! Type: REAL*8
c
c//! Variable: Lat
c//! Description:  Latitude of site (in radians) 
c//! Units: radians      
c//! Type: REAL*8
c
c//! Variable: Cld     
c//! Description: Cloudiness of the day 
c//! Units: proportion 
c//! Type: REAL*8
c
      RETURN
      END

      REAL*8 FUNCTION DIFFPROP_DAY_FN(TRANS_PROP_DAY)
      IMPLICIT NONE

c**   determines the proportion of light  that is diffuse on a daily 
c**   basis
      REAL*8 TRANS_PROP_DAY
c
      IF(TRANS_PROP_DAY.LT.0.07) THEN
	DIFFPROP_DAY_FN = 1.0
      ELSE IF(TRANS_PROP_DAY.LT.0.35) THEN
	DIFFPROP_DAY_FN = 1.0-2.3*(TRANS_PROP_DAY-0.07)**2
      ELSE IF(TRANS_PROP_DAY.LT.0.75) THEN
	DIFFPROP_DAY_FN = 1.33-1.46*TRANS_PROP_DAY
      ELSE
	DIFFPROP_DAY_FN = 0.23
      ENDIF
c
c//! Equation: 25
c//! Description: The proportion of light reaching the ground that \
c//!   is made up of diffuse component                      \\
c//!   INPUTS: proportion of extra-terrestrial radiation to terrestrial \
c//!    radiation transmitted
c//! Bibliography: Spitters Et al (1986), Agric & Forest Meteo  \
c//!  38:217-229; de Jong (1980).
c  
c//! Case: Tran < 0.07
c//!  E: Diff_{prop} = 1.0
c//! Case: 0.07<= Tran < 0.35
c//!  E: Diff_{prop} = 1-0.23(Tran - 0.07)^{2}
c//! Case: 0.35<= Tran < 0.75
c//!  E: Diff_{prop} = 1.33 - 1.46 Tran
c//! Case: 0.75 <= Tran
c//!  E: Diff_{prop} = 0.23
c
c//! Variable: Diff_{prop}
c//! Description: Proportion of daily terrestrial radiation that is \
c//!   diffuse                                              \\
c//!   INPUTS: Extra-terrestrial daily radiation; Terrestrial daily \
c//!   radiation
c//! Units: none (Ratio)
c//! Type: REAL*8
c//! SET
c
c//! Variable: Tran       
c//! Description: Proportion of Terrestrial radiation to    \
c//!   Extra-terrestrial radiation
c//! Units: none (Ratio)
c//! Type: REAL*8
c
      RETURN
      END


      REAL*8 FUNCTION SOLAR_ELLIPSE_FN(DAY_OF_YEAR, 
     1   DAYS_IN_YEAR, PI)
      IMPLICIT NONE

c**   calculates the scalar of the elliptical orbit of the sun around 
c**   the earth. About +/- 3%; Max in mud-summer
      REAL*8 GAMMA, PI
      INTEGER*4 DAY_OF_YEAR, DAYS_IN_YEAR
										
      GAMMA = 2*PI*(DAY_OF_YEAR-1)/DBLE(DAYS_IN_YEAR)
      SOLAR_ELLIPSE_FN = 1.00011+0.034211*COS(GAMMA) +
     1  0.00128*SIN(GAMMA)+0.000719*COS(2*GAMMA) + 0.000077*SIN(2*GAMMA)
c
c//! Equation: 26
c//! Description: Returns the scalar of the elliptical orbit of the \
c//!   sun around the earth (approx +/- 3%; peak mid-summmer)    \\
c//!   INPUTS: Day-number of the year; Number of days in the year;\\
c//!    mathematical constant, pi
c//! Bibliography: S. Evans: SWELTER
c
c//! E:  Ellipse = 1.00011+0.034211 cos(\Gamma) +              \
c//!     0.00128 sin(\Gamma)+0.000719 cos(2 \Gamma) +         \
c//!     0.000077*sin(2 \Gamma)                           \\
c//!     \Gamma = {2 \pi (DoY-1)} \over{DinY}
c
c//! Variable: Ellipse
c//! Description: Scalar of solar ellipse of the sun around the earth 
c//! Units: none
c//! Type: REAL*8
c//! SET
c
c//! Variable: \Gamma
c//! Description: Intermediate variable
c//! Units: none
c//! Type: REAL*8
c
c//! Variable: DoY
c//! Description: Day of the year,(1st jan = 1, 31 dec=365 or 366) 
c//! Units: none          
c//! Type: Integer           
c
c//! Variable: DinY
c//! Description: Total number of Days in year,(leap year=366)
c//! Units: none          
c//! Type: Integer           
c
c//! Variable; \pi
c//! Description: the constant, pi, 3.12415...
c//! Units: (radian)
c//! Type: REAL*8
c
      RETURN
      END


      REAL*8 FUNCTION SOLAR_ELLIPSE2_FN(DAY_OF_YEAR) 
      IMPLICIT NONE

c**   calculates the scalar of the elliptical orbit of the sun around 
c**   the earth. About +/- 3%; Max in mud-summer
      INTEGER*4 DAY_OF_YEAR
										
      SOLAR_ELLIPSE2_FN = 1.033-0.066*sqrt(1.0-((DAY_OF_YEAR-183)**2)
     1   /(183.**2) )  
c
c//! Equation: 27
c//! Description: Returns the scalar of the elliptical orbit of the \
c//!   sun around the earth (approx +/- 3%; peak mid-summmer)    \\
c//!   INPUTS: Day-number of the year
c//! Bibliography: 
c
c//! E:  Ellipse = 1.033-0.066 \sqrt{{(1.0-((DoY-183)**2)} \over    \
c//!          {(183.**2)} )}  
c
c//! Variable: Ellipse
c//! Description: Scalar of solar ellipse of the sun around the earth 
c//! Units: none
c//! Type: REAL*8
c//! SET
c
c//! Variable: DoY
c//! Description: Day of the year,(1st jan = 1, 31 dec=365 or 366) 
c//! Units: none          
c//! Type: Integer           
c
      RETURN
      END


      REAL*8 FUNCTION SOLAR_TIME_FN(GMT, LONGITUDE_RAD, 
     1   DAY_OF_YEAR, PI)
      IMPLICIT NONE

c     Routine calculates the local solar time from GMT. Applies only 
c     to GMT based time.
      REAL*8 GMT, EQNTIME, CORR1, PI
      REAL*8 LONGITUDE_RAD
      INTEGER*4 DAY_OF_YEAR

c      correct time by equation of time and the longitudinal offset    
       CORR1 = (279.575+0.986*DAY_OF_YEAR)*PI/180.                         
c      Equation of time:                                                
       EQNTIME = (-104.7*SIN(CORR1) + 596.2*SIN(2*CORR1) +                   
     1  4.3*SIN(3*CORR1) - 12.7*SIN(4*CORR1) - 429.3*COS(CORR1) -       
     2  2*COS(2*CORR1) + 19.3*COS(3*CORR1))/3600.                       
c      
       SOLAR_TIME_FN = GMT+EQNTIME+LONGITUDE_RAD*180.0/PI/15.0
c
c//! Equation: 28
c//! Description: Converts from GMT (hours) to local solar time;   \
c//!  ie accounts for longitude and the equation of time.     \\
c//!  INPUTS: Grenwich mean-time; longitude; day-number of the year;\
c//!   mathematical constant, pi
c//! Bibliography: Plaz & Grief (1996); Gronbeck 1999 (Web Pages)
c
c//! E:  Sol_{time} = GMT + EQN_{time}+ {{180 Long } \over{\pi 15}}
c
c//! Variable: Sol_{time}
c//! Description: Local solar time  - corrected solar time (accounts \
c//!  for longitude etc.)
c//!  Units: Hours
c//!  Type: REAL*8
c//!  SET
c 
c//! Variable: GMT
c//! Description: Time of day  (Grenwich mean time)
c//!  Units: Hours
c//!  Type: REAL*8
c
c//! Variable: EQN_{time}
c//! Description: Equation of time ; corrects for the slightly   \
c//!   elliptical orbit of the earth.
c//!  Units: Hours
c//!  Type: REAL*8
c
c//! Variable: Long      
c//! Description: Longitude of Site 
c//!  Units: Radians
c//!  Type: REAL*8
c
c//! Variable: \pi
c//! Description: Mathematical constant, pi: 3.1415...
c//! Units: (radian)
c//! Type: REAL*8
c
       RETURN
       END


      REAL*8 FUNCTION GMT_FN(SOLAR_TIME, LONGITUDE_RAD, 
     1   DAY_OF_YEAR, PI)
      IMPLICIT NONE

c     Routine calculates the GMT from solar time. Applies only to GMT
c     based time.
      REAL*8 SOLAR_TIME, EQNTIME, CORR1, PI
      REAL*8 LONGITUDE_RAD
      INTEGER*4 DAY_OF_YEAR

c      correct time by equation of time and the longitudinal offset    
       CORR1 = (279.575+0.986*DAY_OF_YEAR)*PI/180.                         
c      Equation of time:                                                
       EQNTIME = (-104.7*SIN(CORR1) + 596.2*SIN(2*CORR1) +                   
     1  4.3*SIN(3*CORR1) - 12.7*SIN(4*CORR1) - 429.3*COS(CORR1) -       
     2  2*COS(2*CORR1) + 19.3*COS(3*CORR1))/3600.                       
c      
       GMT_FN = SOLAR_TIME -(EQNTIME+LONGITUDE_RAD*180.0/PI/15.0)
c
c//! Equation: 29
c//! Description: Converts from local solar time to GMT;   \
c//!  ie accounts for longitude and the equation of time.     \\
c//!  INPUTS: Local Solar time; longitude; day-number of the year; Pi
c//! Bibliography: Plaz & Grief (1996); Gronbeck 1999 (Web Pages)
c
c//! E:  GMT = Sol_{time} -(EQN_{time} +{{180 Long } \over{\pi 15}})
c
c//! Variable: GMT
c//! Description:  Grenwich Mean Time of day (corrected for longitude \
c//!  etc.)
c//!  Units: Hours
c//!  Type: REAL*8
c//!  SET
c 
c//! Variable: Sol_{time}
c//! Description: Local solar time  - corrected solar time 
c//!  Units: Hours
c//!  Type: REAL*8
c
c//! Variable: EQN_{time}
c//! Description: Equation of time ; corrects for the slightly   \
c//!   elliptical orbit of the earth.
c//!  Units: Hours
c//!  Type: REAL*8
c
c//! Variable: Long      
c//! Description: Longitude of Site 
c//!  Units: Radians
c//!  Type: REAL*8
c
c//! Variable: \pi
c//! Description: Mathematical constant, pi: 3.1415...
c//! Units: (radian)
c//! Type: REAL*8
c
      RETURN
      END


      REAL*8 FUNCTION TEMP_HOUR_FN(HOUR, TEMPMAX, TEMPMIN,PI)
      IMPLICIT NONE

c**  This function fits a sine curve over the day, ranging from 
c**  TEMPMIN to TEMPMAX, with the maximum when HOUR=14 (ie 2pm)
      REAL*8 HOUR, TEMPMAX, TEMPMIN
      REAL*8 PI, TMEAN
c
      TMEAN = ((TEMPMAX-TEMPMIN)/2.0)+TEMPMIN                       
      TEMP_HOUR_FN = SIN(((HOUR-2)/24.0*360-90)*PI/180.0)*
     1   ((TEMPMAX-TEMPMIN)/2.0) + TMEAN
c
c//! Equation: 30
c//! Description: returns the temperature at the time of day from \
c//!    fitting a sine curve with peaks/troughs at the maximum and \
c//!    minimum of the day. Peak assumed to be at time=14:00. Time is \
c//!    assumed to be absolute ie no corrections for latitude are made.\\
c//!  INPUTS: Hour of the day; Maximum temperature of the day; \
c//!    minimum temperature of the day. 
c//! Bibliography: none
c
c//!   E: TEMP_{hour} = ({{T_{max}-T_{min}} \over{2}} + T_{min}) +    \
c//!          sin( { ({{{360(hour-2)}\over{24}} - 90} \over {180} )}  \
c//!          ({{T_{max}-T_{min}} \over{2}})
c//!
c//!  Variable: TEMP_{hour}
c//!  Description: Temperature of day at time 'hour'
c//!  Units: Degrees celcius  (usually)
c//!  Type: REAL*8 
c//!  SET      
c
c//!  Variable: T_{max}
c//!  Description: Maximum temperature of the day
c//!  Units: Degrees celcius  (usually)
c//!  Type: REAL*8
c
c//!  Variable: T_{min}
c//!  Description: Minumum temperature of the day
c//!  Units: Degrees celcius  (usually)
c//!  Type: REAL*8
c
c//!  Variable: \pi     
c//!  Description: Mathematical value 3.1415..
c//!  Units: none 
c//!  Type: REAL*8
c
      RETURN
      END

      REAL*8 FUNCTION HUMIDITY_INST_FN(VP_UNSAT_BASE,
     1   VP_SAT_NOW)
      IMPLICIT NONE

c**  Function calculates the humidity at any temperature
c**  Values at temperatures below freezing are doubtful.
      REAL*8  VP_UNSAT_BASE, VP_SAT_NOW
c
       HUMIDITY_INST_FN = 100*VP_UNSAT_BASE/VP_SAT_NOW                                   
c//! Equation: 31
c//! Description: Gets humidity at any time of the day eg 85%. \
c//!   Assumes that vapour pressure is constant through the day. \
c//!   usually we assume that VP is constant through the day and use\
c//!   the '9am' values are used to get the vp; as temperature \
c//!   changes, so will humidity. Beware if temperature < ZERO. \\
c//!   INPUTS: vapour pressure at base time;   \
c//!    Saturated VP at current temp.
c//! Bibliography: LI-Cor Manual ??
c//!  
c//!  E: HUMID_{t} =  {{100 VP_{unsat,base}} \over{VP_{sat,t}} } 
c
c//! Variable: HUMID_{t}
c//! Description: humidity at temperature, t. Based on holding vapour \
c//!  pressure constant through the day. Beware if temperature < ZERO \\
c//! Units: Percentage
c//! Type: REAL*8
c//! SET
c
c//! Variable: VP_{unsat,base}
c//! Description: vapour pressure at reference time eg 9am
c//! Units: Pascal (Pa)
c//! Type: REAL*8
c
c//! Variable: VP_{sat,t}
c//! Description: Saturated vapour pressure at temperature, t
c//! Units: Pascal (Pa)
c//! Type: REAL*8
c
      RETURN
      END


      SUBROUTINE RAIN_DAY_SUB(NO_RAIN_DAYS, RAIN_IN_MONTH,  
     1  DAYS_IN_MONTH, SEED1, SEED2, SEED3, RAIN)
      IMPLICIT NONE



c** The routine simulates daily rainfall for all days in a month
c** The typical number of rainy-days and average rainfall for the month 
c** Are the only specific Meteorology inputs. Returned is an array
c** (1-dimensional; day) with the simulated rainfall for each day 
c** in it.
      REAL*8 RAIN(31)
      REAL*8 FRACT_WET, RAIN_DAY, RAIN_IN_MONTH, NO_RAIN_DAYS
      REAL*8 ZERO
      REAL*8 PROB_WET_WET, PROB_WET_DRY, GAMMA, MU
      INTEGER*4 DAYS_IN_MONTH, IDAY
      INTEGER*4 SEED1, SEED2, SEED3, IP, IP0
      REAL*8  RNDF, RAND01 , XLO, XHI, gamma2
      REAL*8 GAMMADIST_FN

c
c **    determine the proportion of the month with wet-days
      FRACT_WET = DBLE(NO_RAIN_DAYS)/DBLE(DAYS_IN_MONTH)
c
c**     find the average rainfall per day with rainfall
      RAIN_DAY = RAIN_IN_MONTH/DBLE(NO_RAIN_DAYS)
c
c**     find transitional probability of a wet-day following a dry-day
c**     of the month (Geng et al 1986)
      PROB_WET_DRY = 0.75*FRACT_WET
c
c**     probability of wet-day following a wet-day (Geng et al 1986)
      PROB_WET_WET = 0.25+PROB_WET_DRY
c
c**     Gamma distribution paramaters (Geng et al 1986) 
c       GAMMA = -2.16+1.83*RAIN_DAY   
        XLO = 0.0
        XHI=1.0

        RAND01 = RNDF(XLO,XHI, SEED1, SEED2, SEED3)
        GAMMA = 1.0
c
c**     Mu paramaters (Geng et al 1986) 
        gamma2 = gamma+(rand01-0.5)/2.0
      MU = RAIN_DAY/GAMMA2
c
c     first of the month; the chances of it being wet is 'fract_wet'
c     get a random number proportional to the number of days in the month
      XLO = 0.0
      XHI=1.0
      RAND01 = RNDF(XLO,XHI, SEED1, SEED2, SEED3)
      DO 20 IDAY = 1, DAYS_IN_MONTH
	IF(IDAY.EQ.1) THEN
c**       is the random number<=rain_fract; if so then it rained yesterday
	  IP0 = 0
	  IF(RAND01.LE.FRACT_WET) THEN
	    IP0 = 1 
	  ENDIF
	ENDIF
	RAIN(IDAY) = 0.0
	RAND01 = RNDF(XLO, XHI, SEED1, SEED2, SEED3)
	IF(IP0.EQ.1) THEN
c**       raining on day 0
	  IF(RAND01-PROB_WET_WET .LE. 0) IP = 1
	  IF(RAND01-PROB_WET_WET .GT. 0) IP = 0
	ELSE   ! IP0 = 0 ; its dry on day 0
	  IF(RAND01-PROB_WET_DRY .LE. 0) IP = 1
	  IF(RAND01-PROB_WET_DRY .GT. 0) IP = 0
	ENDIF
	IP0 = IP
c
	IF(IP.EQ.1) THEN     ! Its raining today
12        CONTINUE
*12	  TR1 = EXP(-18.42*MU)

*	  TR2 = EXP(-18.42*(1-MU))


c         RAND01 = RNDF(XLO, XHI, SEED1, SEED2, SEED3)
c         IF(RAND01 - TR1 .LE. 0) THEN
c           S1 = 0
c         ELSE
c           S1 = RAND01*EXP(1/MU)
c         ENDIF
c         RAND01 = RNDF(XLO, XHI, SEED1, SEED2, SEED3)
c         IF(RAND01 - TR2 .LE. 0) THEN
c           S2 = 0
c         ELSE
c           S2 = RAND01*EXP(1/(1-MU))
c         ENDIF
c         IF(S1+S2-1 .LE. 0) THEN
c           IF( (S1+S2).EQ.0 )THEN
c             Z = 0.0
c           ELSE
c             Z = S1/(S1+S2)
c           ENDIF
c         ELSE
c           FLAG = 1
c         ENDIF
c         IF(FLAG.EQ.1) GOTO  30
c
c         RAND01 = RNDF(XLO, XHI, SEED1, SEED2, SEED3)
c         RAIN(IDAY) = -Z*LOG(RAND01)*GAMMA
	  ZERO = 0.0

C Changed GAMMA to GAMMA2 14/12/01 Paul Henshall.
        
  	  RAIN(IDAY)=GAMMADIST_FN(zero,MU,GAMMA2,SEED1,SEED2,SEED3)
        if(rain(iday).gt.rain_day*5)goto 12
	ENDIF
20    CONTINUE
c
c//! Equation: 32
c//! Description: The routine simulates daily rainfall for all days in \
c//!  a month. The typical number of rainy-days and average rainfall \
c//!  for the month are the only specific Meteorology inputs. Returned \
c//!  is an array (1-dimensional; day) with the simulated \
c//!  rainfall for each day in it. Although the routine uses  \
c//!  Markov-chain principles, the start of each month is INDEPENDENT \
c//!  of what the conditions were on the last day of the previous  \
c//!  month. \\
c//!  The routine calls a random-number gernerator (RNDF) which \
c//!  returns numbers between Xmin and Xmax (in this case 0,1) using \
c//!  the 3 seeds \\
c//!   INPUTS: Number of rainy days in the month; Average rainfall of \
c//!    the month; number of days in the month; 3 seeds \
c//!    for random number generation\\
c//!   OUTPUT: array (dimension:31), with rainfall for each day
c//! Bibliography: S.Evans (SWELTER)
c
c//! Case: Day_{wet,dry}^{i} = 0
c//!  E: Rain_{day}^{i} = 0.0 \\
c//! \\
c//! Case: Day_{wet,dry}^{i} = 1 \\
c//! \\
c//!  E: Rain_{day}^{i} = Z \gamma log_{e} (Rnd_{j1,j2,j3})  \\
c//!   \gamma = -2.16+1.83 Rain_{day,mean}  \\
c//! \\
c//!  Case: S_{1} + S_{2} -1 <= 0    
c//!   E: Z = { {{S_{1}} \over{S_{1}+S_{2}} }  \\
c//! \\
c//!     Case: Rnd_{j1,j2,j3} - Tr_{1} <= 0
c//!      E: S_{1} = 0
c//!     Case: Rnd_{j1,j2,j3} - Tr_{1} > 0
c//!      E: S_{1} = Rnd_{j1,j2,j3} exp(1/\mu) \\
c//! \\
c//!     Case: Rnd_{j1,j2,j3} - Tr_{2} <= 0
c//!      E: S_{2} = 0
c//!     Case: Rnd_{j1,j2,j3} - Tr_{2} > 0
c//!      E: S_{2} = Rnd_{j1,j2,j3} exp(1/(1-\mu))   \\
c//! \\
c//!  Case: S_{1} + S_{2} -1 > 0
c//!     E: re-generate S_{1} and S_{2}  \\
c//!    \\ 
c//!  Tr_{1} = {exp(-18.42)}\over{\mu} \\
c//!  Tr_{2} = {exp(-18.42)}\over{(1-mu)} \\
c//!  \mu = {{Rain_{day,mean}} \over{\gamma}} \\
c//! \\
c//!  Case: Day_{wet,dry}^{i-1} = 1 \\
c//! \\
c//!   Case: Rnd_{j1,j2,j3} - Prob_{wet,wet} <= 0
c//!    E: Day_{wet,dry}^{i} = 1
c//!   Case: Rnd - Prob_{wet,wet} > 0
c//!    E: Day_{wet,dry}^{i} = 0 \\
c//! \\
c//!  Case: Day_{wet,dry}^{i-1} = 0 \\
c//! \\
c//!   Case: Rnd_{j1,j2,j3} - Prob_{wet,dry} <= 0
c//!    E: Day_{wet,dry}^{i} = 1      
c//!   Case: Rnd_{j1,j2,j3} - Prob_{wet,dry} > 0 
c//!    E: Day_{wet,dry}^{i} = 0      \\
c//! \\
c//!  Prob_{wet,wet} = 0.25+Prob{wet,dry} \\
c//!  Prob_{wet,dry} = { {0.75 N_{rain days}} \over{DinM} }
c
c//! Variable: Rain_{day}^{i}
c//! Description: Predicted daily rainfall on day of month, i
c//! Units: mm
c//! Type: REAL*8 array (12,31)
c//! SET
c
c//! Variable: Rain_{day,mean}
c//! Description: Mean daily rainfall for the month
c//! Units: mm/day
c//! Type: REAL*8
c
c//! Variable: Rnd_{j1,j2,j3)
c//! Description: Random number generated between 0 and 1
c//! Units: none
c//! Type: REAL*8
c
c//! Variable: j1
c//! Description: Seed for random number generator; changed on return
c//! Units: none
c//! Type: Integer
c
c//! Variable: j2
c//! Description: Seed for random number generator; changed on return
c//! Units: none
c//! Type: Integer
c
c//! Variable: j3
c//! Description: Seed for random number generator; changed on return
c//! Units: none
c//! Type: Integer
c
c//! Variable: DinM
c//! Description: Number of days in the month
c//! Units: none
c//! Type: Integer
c
      RETURN
      END

 
      SUBROUTINE RAIN_DURATION_SUB(RAIN, RAIN_DURATION, RAIN_INTENSITY, 
     1  DAYS_IN_MONTH)
      IMPLICIT NONE

c** This routine estimates the duration of rain on each day for a month
c** Accuracy is suspect.
      REAL*8 RAIN(31) , RAIN_DURATION(31), RAIN_INTENSITY(31)
      REAL*8 RAIN_LIM(10), PRANGE, DUR(10)
      INTEGER*4 DAYS_IN_MONTH, IDAY
      INTEGER*4 N(10), RANGE
c
      RAIN_LIM(1) = 0
      RAIN_LIM(2) = 5
      RAIN_LIM(3) = 10
      RAIN_LIM(4) = 15
      RAIN_LIM(5) = 20
      RAIN_LIM(6) = 25
      RAIN_LIM(7) = 50
      RAIN_LIM(8) = 75
      RAIN_LIM(9) = 100
      DO RANGE = 1,10
        N(RANGE) = 0
      ENDDO
c
      DO 60 IDAY = 1, DAYS_IN_MONTH
	DO 70 RANGE = 1, 8,1
	  IF( (RAIN(IDAY).GT.RAIN_LIM(RANGE)) .AND. 
     1       (RAIN(IDAY).LE.RAIN_LIM(RANGE+1)) ) THEN
	    N(RANGE) = N(RANGE)+1
	  ENDIF
70      CONTINUE
60    CONTINUE
c     
      DO 80 RANGE=1, 8, 1
	PRANGE = (RAIN_LIM(RANGE+1) - RAIN_LIM(RANGE)) / 25.4
	DUR(RANGE) = N(RANGE)/ (1.39*((PRANGE+0.1)** (-3.55)))
80    CONTINUE
c
      DO 90 IDAY = 1,  DAYS_IN_MONTH
	RAIN_DURATION(IDAY) = 0.0
	RAIN_INTENSITY(IDAY) = 0.0
	IF(RAIN(IDAY).GT.0) THEN
	  DO 100 RANGE =1,8,1
	    IF( (RAIN(IDAY).GT.RAIN_LIM(RANGE)) .AND. 
     1         (RAIN(IDAY).LE.RAIN_LIM(RANGE+1)) ) THEN
	      PRANGE = (RAIN_LIM(RANGE+1) - RAIN_LIM(RANGE)) 
	      RAIN_DURATION(IDAY) = DUR(RANGE)*RAIN(IDAY)/
     1         PRANGE*60
	    ENDIF
100       CONTINUE
	RAIN_INTENSITY(IDAY) = RAIN(IDAY)/(RAIN_DURATION(IDAY)/60.0)
	ENDIF

90    CONTINUE
c
c//! Equation: 33
c//! Description: This routine takes a months-worth of daily rainfall \
c//!  data and calculates the duration of each event and intensity of\
c//!  each event. Only one Event occurs each day. Accuracy seems \\
c//!   dubious. \\
c//!   INPUTS: Daily rainfall (array 31); Duaration of daily \
c//!    rainfall [OUTPUT], (array 31); intensity of rainfall [OUTPUT]\
c//!    ,(array 31); number of days in the month
c//! Bibliography: S. Evans (SWELTER)
c
c//! Case: no rain
c//!  E: Rain_{dur}_{i} = 0.0
c//! Case: Rainfall
c//!  E: Rain_{dur}_{i} = {{60 Rain_{dur:c} Rain_{day}_{i}}     \\
c//!     \over {Range_{c}}}   \\
c//!   Range_{c} = Rainfall category range: x <rain> y  :\\
c//!    [<5, 10, 15, 20, 25, 50, 75, 100] mm   \\
c//!  \\
c//!   Rain_{dur:c}={ {\sum_{i=1}^{DinM} Occurences of Rain_{day}_{i} \\
c//!    within rainfall category c}  \over {1.39 \\
c//!    ((Range_{c}/25.4 + 0.1)^{-3.55})} } \\
c//! \\
c//!   Rain_{inten} = Rain_{day}^{i}/(Rain_{dur}^{i}/60.0)
c//!
c
c//! Variable: Rain_{dur}_{i}
c//! Description: Duration of rainfall on day of month, i
c//! Units: minutes
c//! Type: REAL*8 array (31)
c//! SET
c
c//! Variable: Rain_{inten}_{i}
c//! Description: Rainfall intensity on day of month, i
c//! Units: mm /hour
c//! Type: REAL*8

c//! Variable: Rain_{day}_{i}
c//! Description: Rainfall on day of month, i
c//! Units: mm 
c//! Type: REAL*8

c//! Variable: DinM
c//! Description: Number of days in the month
c//! Units: none
c//! Type: Integer
c
      RETURN
      END


      REAL*8 FUNCTION RNDF(XLO, XHI,J1, J2, J3)
      IMPLICIT NONE

c** This random number routine w is from Brian Wichman and David Hill
c** (BYTE, March 1987, page 127-128). Apparently the cycle length is
c** about 6.95E12. Original reference: Wichmann & Hill, App stats (1982)        
c** 31(2) pp 188-190 with modification: A Mcleod, App stats (1985)            
c** 34(2) pp 198-200                                             
c
      INTEGER*4 J1, J2, J3                                          
      REAL*8 R,XLO,XHI
c                                                                   
      IF(XLO.LT.XHI)THEN                                            
	J1=171*MOD(J1,177)-2*(J1/177)                                  
	IF(J1.LT.0) J1 = J1+30269                                         
c                                                                   
	J2=172*MOD(J2,176)-35*(J2/176)                                 
	IF(J2.LT.0)J2=J2+30307                                           
c                                                                     
	J3=170*MOD(J3,178)-63*(J3/178)  
        IF(J3.LT.0) J3 = J3+30323                                           
c                                                                     
	R=DBLE(J1)/30269.0D0+DBLE(J2)/30307.0d0+DBLE(J3)/30323.0D0       
c** modification re: A McLeod (1985) app. stats 34(2) 198-200     
c** returns out of loops for efficiency - otherwise it was just   
C**  a question of scaling.                                        
	IF(R.GT.0) THEN                                               
	  RNDF=(R-INT(R))*(XHI-XLO)+XLO                               
	  RETURN                                                      
	ENDIF                                                         
	R=DMOD(DBLE(FLOAT(J1))/30269.0D0 + DBLE(FLOAT(J2))/30307.0D0 +  
     1     DBLE(FLOAT(J3))/30323.0D0, 1.0D0)                           
	IF(R.GE.1.0) R=0.999999                                       
	RNDF=(R-INT(R))*(XHI-XLO)+XLO                                 
	RETURN                                                        
      ENDIF                                                           
c                                                                     
      RNDF=XLO                                                        
c
c//! Equation: 34
c//! Description: Random number generator: uses 3 seeds and returns \
c//!  changed seeds and a random number between Xlo and Xhi\\
c//!   INPUTS: Lowest value of random number; highest value of random \
c//!    number; 3 seeds
c//! Bibliography: Wichmann & Hill (1882), App stats 31(2) pp188-190; \
c//!  A McLeod (1985) app. stats 34(2) pp98-200
c
c//! Case: R > 0
c//!  E: Rnd_{j1,j2,j3}=(R-INT(R))*(X_{hi}-X_{lo})+X_{lo}
c//! Case: R<= 0
c//!  E: RNDF=(R1-INT(R1))*(X_{hi}-X_{lo})+X_{lo}  \\
c//! \\
c//! E: R = j1/30269 + j2/30307 + j3/30323  \\
c//!   \\
c//! E: R1 = MAX(0.999999, MOD(S_{1}/30269 + S_{2}/30307 + \
c//!     S_{3}/30323.0D0, 1.0)                           \\
c//! \\
c//!   E:  j1 = 171 MOD(j1,177)-2(j}/177)  \\
c//!     IF(j1 < 0  ; j1 = j1 + 30269   \\
c//!  \\
c//!   E:  j2 = 172 MOD(j2,176)-35(j2/176)  \\
c//!     IF(j2 < 0  ; j2 = j2 + 30307   \\
c//!  \\
c//!   E:  j3 = 170 MOD(j3,178)-63(j3/178)  \\
c//!     IF(j3 < 0  ; j3 = j3 + 30323   
c  
c//! Variable: Rnd_{j1,j2,j3}
c//! Description: Returned random number
c//! Units: none
c//! Type: REAL*8
c//! SET
c
c//! Variable: j1
c//! Description:  Seed for random number generation changed on exit
c//! Units: none
c//! Type: Integer
c//! SET
c
c//! Variable: j2
c//! Description:  Seed for random number generation changed on exit
c//! Units: none
c//! Type: Integer
c//! SET
c
c//! Variable: j3
c//! Description:  Seed for random number generation changed on exit
c//! Units: none
c//! Type: Integer
c//! SET
c
c//! Variable:  X_{lo}
c//! Description:  lowest possible value of random number
c//! Units: none
c//! Type: REAL*8
c
c//! Variable:  X_{hi}
c//! Description:  highest possible value of random number
c//! Units: none
c//! Type: REAL*8
c
      RETURN                                                          
      END


      REAL*8 FUNCTION TRANS_PROP_DAY_FN(RAD_ET_DAY, 
     1  RAD_TERRA_SW_DAY)
      IMPLICIT NONE

c** Evaluates the total trasmission proportion (coefficient) for the day
      REAL*8 RAD_TERRA_SW_DAY
      REAL*8 RAD_ET_DAY
c
      TRANS_PROP_DAY_FN = RAD_TERRA_SW_DAY/RAD_ET_DAY
c
c//! Equation: 35 
c//! Description: Works out the proportion of Extra-terrestrial \
c//!  radiation that arrives on the earth's surface as short-wave \
c//!  in a day \\
c//!   INPUTS: Daily short-wave radiation on earth's surface; daily \
c//!    Extra-terrestrial radiation.
c//! Bibliography: S. Evans (SWELTER); Liu & Jordan (1960) Solar Energy \
c//!  4 pp1-19.
c
c//! E:  T_{sw} = {RAD_{terra, sw} /over{RAD_{et}}}
c
c//! Variable: T_{sw}
c//! Description: proportion of Extra-terrestrial radiation that \
c//!  arrives on the earth's surface as short-wave.
c//! Units: none
c//! Type: REAL*8
c//! SET
c
c//! Variable: RAD_{terra,sw}
c//! Description: Terrestrial short-wave radiation on a horizonatal \
c//!  plane for the given latitude and day.
c//! Units: W/m2/d
c//! Type: REAL*8
c                                                                    
c//! Variable: RAD_{et}   
c//! Description: Extra-Terrestrial radiation above the atmosphere  \
c//!  for the given latitude and day.
c//! Units: W/m2/d
c//! Type: REAL*8 
c
      RETURN
      END

      REAL*8 FUNCTION ANGSTROM_BETA_FN(LATITUDE_RAD)
      IMPLICIT NONE

c** finds the Angstrom turbidity factor (the max clear-sky atmospheric
c** transmittance characteristics at the latitude
      REAL*8 LATITUDE_RAD
c
      ANGSTROM_BETA_FN = 0.682 - 0.3183*(1-1.3614*COS(LATITUDE_RAD))
c
c//! Equation: 36
c//! Description: Caculation of the Angstrom beta function value: A \
c//!  turbidity factor for the maximum clear-sky transmittance \
c//!  characteristics at a given latitude. \\
c//!   INPUTS: Latitude
c//! Bibliography: Nikolov & Zeller (1992) Ecological Modelling \
c//!  61:149-168; Gueymard (1989) Agric & Forest Met 45:215-229
c
c//! E: \beta =  0.682 - 0.3183*(1-1.3614*cos(Lat))
c
c//! Variable: \beta
c//! Description: Angstrom tubidity factor: max clear sky transmittance
c//! Units: none  
c//! Type: REAL*8
c//! SET
c
c//! Variable: Lat
c//! Description:  Latitude of site (in radians) 
c//! Units: radians      
c//! Type: REAL*8
c
      RETURN
      END


      REAL*8 FUNCTION DIFFPROP_DAY2_FN(TRANS_PROP_DAY,
     1  BETA)
      IMPLICIT NONE

c**Estimates the proportion of terrestrial radiation that is diffuse
      REAL*8 TRANS_PROP_DAY, BETA
c
      DIFFPROP_DAY2_FN = TRANS_PROP_DAY*(1-exp((-0.6*(1-BETA)/
     1  TRANS_PROP_DAY) / (BETA-0.4)) )
c
c//! Equation: 37
c//! Description: The proportion of light reaching the ground that \
c//!   is made up of diffuse component                      \\
c//!   INPUTS: Daily extra-terrestrial radiation; Daily terrestrial \
c//!    radiation; Angstrom beta value
c//! Bibliography: S Evans (SWELTER); Bristow & Campbell (1985) Agric \
c//!  & Forest Met 35:123-131; Becker & Weingarten (1991) Agric &  \
c//!  Forest Met 53:347-353
c
c//! E: Diff_{prop} =   Tran(1-exp{({{-0.6(1-\beta)} \over{Tran}} \
c//!     \over{\beta-0.4})}
c
c//! Variable: Diff_{prop}                                           
c//! Description: Proportion of daily terrestrial radiation that is \
c//!   diffuse                                              \\       
c//!   INPUTS: Extra-terrestrial daily radiation; Terrestrial daily \
c//!   radiation                                                     
c//! Units: none (Ratio)                                             
c//! Type: REAL*8                                          
c//! SET                                                             
c                                                                    
c//! Variable: Tran                                                  
c//! Description: Proportion of Terrestrial radiation to    \
c//!   Extra-terrestrial radiation                                   
c//! Units: none (Ratio)                                             
c//! Type: REAL*8                                          
c                                                                    
c//! Variable: \beta
c//! Description: Angstrom tubidity factor: max clear sky transmittance
c//! Units: none  
c//! Type: REAL*8
c
      RETURN
      END


      REAL*8 FUNCTION RAD_TERRA_INST_FN(
     1  DAWN, DUSK, DAYLENGTH, HOUR, RAD_TERRA_DAY, PI)             
      IMPLICIT NONE

c** calculate solar radiation for  at any instant in time        
c** Assume sine wave with maximum at noon so ideally use local solar time
      REAL*8 RAD_TERRA_DAY
      REAL*8 DAWN, DUSK, DAYLENGTH, HOUR, SUNTIME
      REAL*8 PI
c
c**Give zero values for radiation at night otherwise set radiation
c
      IF(HOUR .LT. DAWN .OR. (HOUR .GT. DUSK)) THEN
	  RAD_TERRA_INST_FN = 0.0
      ELSE
c** Sinusoidal function for hourly radiation from daily total    
	SUNTIME = HOUR - DAWN
	RAD_TERRA_INST_FN = RAD_TERRA_DAY*PI/(2*DAYLENGTH)*
     1    SIN(PI/DAYLENGTH*SUNTIME)/3600.0 
      ENDIF
c
c//! Equation: 38
c//! Description: Fit a sine curve over the day - with peak at noon \
c//!  (so local solar time is best) and periodicity over the  \
c//!  daylength \\
c//!   INPUTS:  Latitude; Time of dawn; Time of dusk; Daylength; \
c//!    Time of day; Daily terrestrial radiation ; Pi
c//! Bibliography: Forest-GROWTH/FLUX; Meastro
c
c//! Case: Night
c//!  E: RAD_{now} = 0.0
c//! Case: Dawn <= hour <= Dusk
c//! E: RAD_{now} = {{RAD_{terra, day} \pi sin(\pi/(DAYL T)} \
c//!       \over{(2 DAYL 3600)}}
c
c//! Variable: RAD_{now}
c//! Description: terrestrial radiation at any point in time; attained \
c//!  from fitting a sine curve between dawn and dusk, peak at noon \
c//!  using the daily radiation. Local solar times should be used
c//! Units: W/m^{2} (shortwave radiation)
c//! Type: REAL*8
c//! SET
c
c//! Variable: Dawn
c//! Description: Time of sunrise 
c//! Units: Hour
c//! Type: REAL*8
c
c//! Variable: Dusk
c//! Description: Time of sunset  
c//! Units: Hour
c//! Type: REAL*8
c
c//! Variable: RAD_{terra, day}
c//! Description:  Daily terrestrial radiation 
c//! Units: W/m^{2}/day
c//! Type: REAL*8
c
c//! Variable: DAYL
c//! Description:  Daylength                   
c//! Units: Hours
c//! Type: REAL*8
c
c//! Variable: hour
c//! Description: Hour of the day - GMT base
c//! Units: hour
c//! Type: REAL*8
c
      RETURN   
      END     

      
      REAL*8 FUNCTION WINDSPEED_INST_FN(WIND_RUN_DAY, HOUR)
      IMPLICIT NONE

c** calculates windspeed at any point in time from a daily run of
c** wind. This is based on a regression to get the daily max and 
c** min windspeed from the run of wind. Winspeed is assumed to be 
c** maximum at 14:00, and changes linearly about that
      REAL*8 WIND_RUN_DAY, HOUR
      REAL*8 WSMIN, WSMINCONST, WSMINCOEFF
      REAL*8 WSMAX, WSMAXCONST, WSMAXCOEFF
c
c ** this is parmaterised from Alice -Holt 1999: Regressing 
c** Run Vs WSMIN gives the coeffs: 0.0075x -0.3181, r2=0.7559
c** Run vs WSMAX gives the coeffs, 0.0142x +0.6976, r2=0.8241
      WSMINCONST = -0.3181 ! -0.2253
      WSMINCOEFF = 0.0075  !  0.0391
      WSMAXCONST = 0.6976  !  0.2822
      WSMAXCOEFF = 0.0142  !  0.0444
      WSMIN=MAX(0.0d0,WIND_RUN_DAY*WSMINCOEFF+WSMINCONST)    
      WSMAX=MAX(0.0d0,WIND_RUN_DAY*WSMAXCOEFF+WSMAXCONST)
c
      IF(HOUR.LE.13.0) THEN                  
	WINDSPEED_INST_FN = MAX(0.1d0,WSMIN+(WSMAX-WSMIN)/13.*(HOUR-1))
      ELSE                                  
	WINDSPEED_INST_FN = MAX(0.1d0,WSMAX-(WSMAX-WSMIN)/11.*(HOUR-13))
      ENDIF
c
c//! Equation: 39
c//! Description: Approximates wind-speed at any time of the day \
c//! from the daily run of wind. \\
c//!   INPUTS: Daily run of wind; Time of day
c//! Bibliography: none ( M Broadmeadow: developed from regression at \
c//!  the Straits)
c
c//! Case: Hour <=13
c//!  E: Wind_{now} = max(0.1, Wind_{min}+(Wind_{max}-Wind_{min}) \
c//!     (Hour-1)/13)
c//! Case: Hour >13
c//!  E: Wind_{now} = max(0.1, Wind_{max}+(Wind_{max}-Wind_{min}) \
c//!     (Hour-13)/11) \\
c//!  \\
c//! Wind_{max} = Wind_{run} WS_{max,1} + WS_{max,2} \\
c//!  \\
c//! Wind_{min} = Wind_{run} WS_{min,1} + WS_{min,2} 
c
c//! Variable: Wind_{now}
c//! Description: Windspeeds at current time
c//! Units: m/s
c//! Type: REAL*8
c//! SET
c
c//! Variable: WS_{max,1}
c//! Description:  Linear regression coefficient of max windspeed Vs \
c//!  run of wind
c//! Units: Km/g / m/s
c//! Type: REAL*8
c
c//! Variable: WS_{max,2}
c//! Description:  Linear regression constant of max windspeed Vs \
c//!  run of wind
c//! Units: m/s
c//! Type: REAL*8
c
c//! Variable: WS_{min,1}
c//! Description:  Linear regression coefficient of min windspeed Vs \
c//!  run of wind
c//! Units: Km/g / m/s
c//! Type: REAL*8
c
c//! Variable: WS_{min,2}
c//! Description:  Linear regression constant of min windspeed Vs \
c//!  run of wind
c//! Units: m/s
c//! Type: REAL*8
c
c//! Variable: Wind_{run}
c//! Description:  Daily run of wind
c//! Units: Km/d
c//! Type: REAL*8
c
c//! Variable: Hour
c//! Description: Hour of the day - GMT base
c//! Units: hour
c//! Type: REAL*8
c
      RETURN                                
      END                                   

 
      REAL*8 FUNCTION RAD_TERRA_SUNDAY_FN(DAYLENGTH,SUN_HOURS,
     1   RAD_ET_DAY) 
      IMPLICIT NONE

c**  Predicts terrestrial radiation from sun-shine hours
c**  Rad_terra= a + b*Sun/daylen + c*(Sun/daylen)**2 + d*(Sun/DAYLEN)**3
c**  Other conversions eg Angstom 1924; Collingbourne 1976 use a linear 
c**  Function; ie 2 params
      REAL*8 CONST, COEFFA, COEFFB, COEFFC, SRATIO
      REAL*8 DAYLENGTH, PRED_GRATIO
      REAL*8 RAD_ET_DAY, SUN_HOURS
c                                                             
      CONST   = 0.17262                                       
      COEFFA   =0.9901                                        
      COEFFB   =-1.043                                        
      COEFFC   = 0.585                                        
c                                                             
	  SRATIO = SUN_HOURS/DAYLENGTH
	  IF(SRATIO.GT.1) THEN                                
	    SRATIO = 1                                        
	  ENDIF                                               
c
	 PRED_GRATIO = CONST + COEFFA*SRATIO + COEFFB*SRATIO**2 +       
     +   COEFFC*SRATIO**3 
	 RAD_TERRA_SUNDAY_FN = RAD_ET_DAY*PRED_GRATIO !J/day
c
c//! Equation: 40
c//! Description: Method of predicting terrestrial radiation from \
c//!  sun-shine hours. Note often only \alpha and \beta are used in  \
c//!  the equation below:  \\
c//!   INPUTS: Daylength; Sunhine-hours; Daily Extra-terrestrial \
c//!   radiation; Mathematical constant, pi
c//! Bibliography: Angstrom (1924) Q.J. R.met. Soc 50,121-5; \
c//!  Collingbourne (1976)  In: The Climate of the British Isles. \
c//!  eds TJ Chandler & S Gregory, Longman, pp74-95; Randle (1997) \
c//!  FC Internal report
c
c//! E: RAD_{terra,day} = RAD_{et,day} (\alpha + \beta*S/S_{0} +  \
c//!  \gamma (S/S_{0})^{2} + \delta (S/S_{0})^{2} )
c
c//! Variable:  RAD_{terra,day}
c//! Description: Daily terrestrial radiation predicted from  \
c//!  sun-shine hours
c//! Units: W/m^{2}/d
c//! Type: REAL*8
c//! SET
c
c//! Variable: RAD_{et,day}
c//! Description: Daily extra-terrestrial radiation
c//! Units: W/m^{2}/d
c//! Type: REAL*8
c
c//! Variable: S
c//! Description: Sun-shine hours in the day
c//! Units: hours     
c//! Type: REAL*8
c
c//! Variable: S_{0}
c//! Description: Daylength
c//! Units: hours     
c//! Type: REAL*8
c
c//! Variable: /alpha
c//! Description: Constant of regression of  S/S_{0} and Rad_{et,day}. \
c//!  If linear regression this is 'Angstrom \alpha'
c//! Units: W/m^{2}/d 
c//! Type: REAL*8
c
c//! Variable: /beta 
c//! Description: First order coefficient of regression of  S/S_{0}  \
c//! and Rad_{et,day}. If linear regression this is 'Angstrom \beta'
c//! Units: W/m^{2}/d 
c//! Type: REAL*8
c
c//! Variable: /gamma
c//! Description:  Second order coefficient of regression of  S/S_{0}. \
c//!  and Rad_{et,day}. If linear regression this is zero
c//! Units: W/m^{2}/d 
c//! Type: REAL*8
c
c//! Variable: /delta 
c//! Description: Constant of regression of  S/S_{0} and Rad_{et,day}. \
c//!  if linear regression this is zero
c//! Units: W/m^{2}/d 
c//! Type: REAL*8
c
c//! Variable: \pi
c//! Description: Mathematical constant, pi
c//! UNits: (Radian)
c//! Type: REAL*8
c
      RETURN                                                             
      END                                                                


      REAL*8 FUNCTION TEMP_AIR_AMPLITUDE_FN(TRANS_PROP_DAY,
     1  A_SKY, ANGSTROM_BETA, C_SKY) 
      IMPLICIT NONE

c** determines apriximate air temperature amplitude (tmax- tmin) from
c** transmitted radiation and some parameters
      REAL*8 TRANS_PROP_DAY
      REAL*8 A_SKY, ANGSTROM_BETA, C_SKY, DELTA

      DELTA = 1-TRANS_PROP_DAY/ANGSTROM_BETA
      IF(DELTA.LE.0) THEN 
	TEMP_AIR_AMPLITUDE_FN = 0.0
      ELSE
	TEMP_AIR_AMPLITUDE_FN = (LOG(DELTA)/(-A_SKY))**(1/C_SKY)
      ENDIF
c
c//! Equation: 41
c//! Description: Works out the daily temperature differecne amplitude \
c//!  from the transmittance of radiation; if it is cloudy, then less \
c//!  radiation gets through and there should be less temperature \
c//!  variation \\
c//!   INPUTS: Proportion of Extra-terrestrial light reaching the \
c//!    earth; Coefficient for maximum clear-sky transmittance;  \
c//!   Angstrom \beta value; Coefficient of clear-sky transmittance \
c//!   with amplitude temperature increase. 
c//! Bibliography: S Evans (SWELTER); Bristow & Campbell (1984) Agric \
c//!  & Forest Met 31: 159-166
c
c//! Case: 1-{{T_{sw}} \over{\beta}} <=0
c//!  E: \Delta T = 0.0
c//! Case:  1-{{T_{sw}} \over{\beta}} > 0
c//!  E: \Delta T = (Log_{e}(\delta/(-A_{sky})) )^(1/C_{sky})
c
c//! Variable: \Delta T
c//! Description: Daily amplitude of temperature variation
c//! Units: Degrees C  (usually)
c//! Type: REAL*8
c//! SET
c
c//! Variable: T_{sw}
c//! Description: proportion of Extra-terrestrial radiation that \
c//!  arrives on the earth's surface as short-wave.
c//! Units: none
c//! Type: REAL*8
c
c//! Variable: \beta
c//! Description: Angstrom tubidity factor: max clear sky transmittance
c//! Units: none  
c//! Type: REAL*8
c
c//! Variable: A_{sky}
c//! Description: Coefficient of maximum clear-sky transmittance \
c//!  characteristics
c//! Units: none  
c//! Type: REAL*8
c
c//! Variable: C_{sky}
c//! Description: Coefficient of maximum clear-sky transmittance \
c//!   with \Delta T increase
c//! Units: none  
c//! Type: REAL*8
c
      RETURN
      END
 
      
      REAL*8 FUNCTION DIFF_PROP_HOUR_FN(DIFF_PROP_DAY, 
     1  SOL_ELEV)
      IMPLICIT NONE

c** finds the proportion of idiifuse light at the earth's surface on
c** an hourly basis (ie elevation should be calculated hour by hour
c** and this feeds back the diffuse light proportion for this hour
      REAL*8 DIFF_PROP_DAY, R, K
      REAL*8 SOL_ELEV
c
      R = 0.847-1.61*SIN(SOL_ELEV) + 1.04*(SIN(SOL_ELEV)**2)
      K = (1.47-R)/1.66         
      IF(DIFF_PROP_DAY.lt.0.22) then 
	 DIFF_PROP_HOUR_FN = 1.0             
      ELSE IF(DIFF_PROP_DAY.lt.0.35) then     
	 DIFF_PROP_HOUR_FN = 1-6.4*(DIFF_PROP_DAY-0.22)**2 
      ELSE IF(DIFF_PROP_DAY.lt.K) then         
	 DIFF_PROP_HOUR_FN = 1.47-1.66*DIFF_PROP_DAY
      ELSE                                                          
	 DIFF_PROP_HOUR_FN = R                                             
      ENDIF                                                         
c
c//! Equation: 42
c//! Description: Finds the proportion of diffuse/total light at the \
c//!  earth's surface on an hourly basis. Elevation should change \
c//!  hour-by-hour. Usually used when hourly terrestrial radiation \
c//!  levels are known. \\
c//!   INPUTS: Proportion of diffuse light (daily basis); \
c//!    Solar elevation 
c//! Bibliography: Spitters Et al (1986), Agric & Forest Met  \
c//!  38:217-229; de Jong (1980).                                
c
c//! Case: Diff_{prop} <= 0.22
c//!  E: Diff_{prop,h} = 1.0
c//! Case: 0.22 < Diff_{prop} <= 0.35
c//!  E: Diff_{prop,h} = 1 - 6.4(Diff_{prop}-0.22)^{2}
c//! Case: 0.35 < Diff_{prop} <= K  
c//!  E: Diff_{prop,h} = 1.47 - 1.66 Diff_{prop}
c//! Case: K < Diff_{prop} 
c//!  E: Diff_{prop,h} = R     \\
c//!   \\
c//! E: K = (1.47-R)/1.66  \\
c//!  \\
c//! E: R = 0.847-1.61 sin(\beta)
c
c//! Variable: Diff_{prop,h}                                            
c//! Description: Proportion of hourly terrestrial radiation that is \
c//!   diffuse at hour h                                    \\
c//!   INPUTS: Proportion of daily radiation that is diffuse; Solar \
c//!    elevation at hour, h
c//! Units: none (Ratio)                                              
c//! Type: REAL*8                                           
c//! SET
c
c//! Variable: Diff_{prop}                                            
c//! Description: Proportion of daily terrestrial radiation that is \
c//!   diffuse                                              \\
c//! Units: none (Ratio)                                              
c//! Type: REAL*8                                           
c                                                                     
c//! Variable: \beta                                                  
c//! Description:  Solar elevarion at hour, h
c//! Units: Radians
c//! Type: REAL*8                                           
c
      RETURN                                                        
      END                                                           


      REAL*8 FUNCTION TEMP_MAX_DAY_FN(TEMP_AIR_AMPLITUDE, 
     1  TEMP_DAY_MEAN)
      IMPLICIT NONE

c** predicts max air temp from the mean temp and amplitude: based on 
c** tmean = (tmax-tmin)/2
      REAL*8 TEMP_AIR_AMPLITUDE
      REAL*8 TEMP_DAY_MEAN
c
      TEMP_MAX_DAY_FN = TEMP_DAY_MEAN + (TEMP_AIR_AMPLITUDE)/2.0
c
c//! Equation: 43
c//! Description: Derives maximum air temperature from the daily mean \
c//!  and daily amplitude. Based on Tmean=(tmax-Tmin)/2  \\
c//!   INPUTS: Daily temperature variation amplitude; Daily mean \
c//!    temperature
c//! Bibliography: none
c
c//! E: T_{max} = T_{mean} + \Delta T /2
c
c//! Variable: T_{max}
c//! Description: Daily maximum temperature
c//! Units: Degrees C   (usually)
c//! Type: Double precsion
c//! SET
c
c//! Variable: T_{mean}
c//! Description: Daily mean temperature
c//! Units: Degrees C   (usually)
c//! Type: Double precsion
c
c//! Variable: \Delta T 
c//! Description: Daily temperature variation amplitude
c//! Units: Degrees C   (usually)
c//! Type: Double precsion
c
      RETURN
      END


      REAL*8 FUNCTION TEMP_MIN_DAY_FN(TEMP_AIR_AMPLITUDE, 
     1  TEMP_DAY_MEAN)
      IMPLICIT NONE

c** predicts min air temp from the mean temp and amplitude: based on 
c** tmean = (tmax-tmin)/2
      REAL*8 TEMP_AIR_AMPLITUDE
      REAL*8 TEMP_DAY_MEAN
c
      TEMP_MIN_DAY_FN = TEMP_DAY_MEAN - (TEMP_AIR_AMPLITUDE)/2.0
c
c//! Equation: 44
c//! Description: Derives minimum air temperature from the daily mean \
c//!  and daily amplitude. Based on Tmean=(tmax-Tmin)/2  \\
c//!   INPUTS: Daily temperature variation amplitude; Daily mean \
c//!    temperature
c//! Bibliography: none
c
c//! E: T_{min} = T_{mean} + \Delta T /2
c
c//! Variable: T_{min}
c//! Description: Daily minimum temperature
c//! Units: Degrees C   (usually)
c//! Type: Double precsion
c//! SET
c
c//! Variable: T_{mean}
c//! Description: Daily mean temperature
c//! Units: Degrees C   (usually)
c//! Type: Double precsion
c
c//! Variable: \Delta T 
c//! Description: Daily temperature variation amplitude
c//! Units: Degrees C   (usually)
c//! Type: Double precsion
c
      RETURN
      END

      REAL*8 FUNCTION TEMP_AIR_AMPLITUDE2_FN(TEMP_MAX_DAY,
     1    TEMP_MIN_DAY)
      IMPLICIT NONE

c** gets the daily temperature variation: ie. Tmax-Tmin
      REAL*8 TEMP_MAX_DAY
      REAL*8 TEMP_MIN_DAY                                             
c
      TEMP_AIR_AMPLITUDE2_FN = TEMP_MAX_DAY - TEMP_MIN_DAY
c
c//! Equation: 45
c//! Description: Derives the daily temperature variation from the \
c//!  maximum and minimums of the day \\
c//!  INPUTS: Temperature maximum; Temperature minimum
c//! Bibliography: none
c
c//! E: \Delta T = T_{max} - T_{min}
c
c//! Variable: \Delta T
c//! Description: Daily amplitude of temperature variation
c//! Units: Degrees C  (usually)
c//! Type: REAL*8
c//! SET
c
c//! Variable: T_{min}
c//! Description: Daily minimum temperature
c//! Units: Degrees C   (usually)
c//! Type: Double precsion
c
c//! Variable: T_{max}
c//! Description: Maximum daily temperature
c//! Units: Degrees C   (usually)
c//! Type: Double precsion
c
      RETURN
      END


      REAL*8 FUNCTION TEMP_DAY_MEAN_FN(TEMP_MAX_DAY,
     1  TEMP_MIN_DAY)
      IMPLICIT NONE

c** calculates mean daily temperature based on (Tmin+Tmax)/2
      REAL*8 TEMP_MAX_DAY, TEMP_MIN_DAY
c
      TEMP_DAY_MEAN_FN = (TEMP_MIN_DAY + TEMP_MAX_DAY)/2.0
c
c//! Equation: 46
c//! Description: calculates Avrage daily temperature from the daily \
c//!  minimum and maximum  \\
c//!   INPUTS: Daily maximum temperature; Daily minimum temperature
c//! Bibliography: none
c
c//! E: Temp_{mean} = (Temp_{min}+Temp_{max})/2.0
c
c//! Variable: T_{mean}
c//! Description: Daily mean temperature
c//! Units: Degrees C   (usually)
c//! Type: Double precsion
c//! SET
c
c//! Variable: T_{min}
c//! Description: Daily minimum temperature
c//! Units: Degrees C   (usually)
c//! Type: Double precsion
c
c//! Variable: T_{max}
c//! Description: Maximum daily temperature
c//! Units: Degrees C   (usually)
c//! Type: Double precsion
c
      RETURN
      END


      REAL*8 FUNCTION TEMP_DAY_MEAN2_FN(TEMP_MEAN_MONTH,
     1  TEMP_SD_MONTH, TEMP_AUTOCORR_MONTH, TEMP_YESTERDAY, IS1,
     2  IS2, IS3, NORMSWITCH)
      
      IMPLICIT NONE
      EXTERNAL NORM_FROM_RAND_FN


c      
c** calculates mean daily temperature  from monthly values   
      REAL*8 TEMP_MEAN_MONTH
      REAL*8 TEMP_SD_MONTH, TEMP_AUTOCORR_MONTH
      REAL*8 TEMP_YESTERDAY
      REAL*8 NORM, NORMMN, NORMSD
      REAL*8 NORM_FROM_RAND_FN
      INTEGER*4 IS1, IS2, IS3, NORMSWITCH
C 
c** approximate a function to get a normal distribution, mean=0, sd=1
      NORMMN=0.0
      NORMSD = 1.0
      NORM=NORM_FROM_RAND_FN(NORMMN, NORMSD, IS1, IS2, IS3, NORMSWITCH)
c
      TEMP_DAY_MEAN2_FN = TEMP_MEAN_MONTH + TEMP_AUTOCORR_MONTH*
     1   (TEMP_YESTERDAY - TEMP_MEAN_MONTH) + TEMP_SD_MONTH*
     2   NORM*( (1-(TEMP_AUTOCORR_MONTH)**2)**0.5)
c
c//! Equation: 47
c//! Description: Calculates Average daily temperature from monthly \
c//!  means and standard-deviations. Autocorrelation of 0.65 seems \
c//!  to work if one isn't available. \\
c//!   INPUTS: Mean temperature of the month; Standard deviation \
c//!    around the mean on a daily basis; Temperature \
c//!     auto-correlation for each month; Mean temperature of the \
c//!     previous day; 3 random number seeds.
c//! Bibliography: S Evans (SWELTER); Haith, Tubbs & Pickering (1984) \
c//!  Simulation of Pollution by soil erosion & soil nutirnt loss, \
c//!   PUDOC, Wageningen
c
c//! E: Temp_{mean} = Temp_{mean,month} + \rho_{T} (Temp_{day-1} - \
c//!     Temp_{mean,month}) + \sigma_{mT} N (1-\sigma_{mT}^{2})^0.5) \\
c//!  \\
c//!   E: N = { Random, N(0,1) }
c
c//! Variable: Temp_{mean}
c//! Description: Daily mean temperature
c//! Units: Degrees C   (usually)
c//! Type: Double precsion
c//! SET
c
c//! Variable: Temp_{mean,month}
c//! Description: Mean daily temperature (monthly basis)
c//! Units: Degrees C   (usually)
c//! Type: Double precsion
c
c//! Variable: Temp_{day-1}     
c//! Description: Mean daily temperature of the previous day
c//! Units: Degrees C   (usually)
c//! Type: Double precsion
c
c//! Variable: \rho_{T}
c//! Description: Observed first-order auto-correlation of mean \
c//!  daily air temperature for each month
c//! Units: none  (day/day^{-1})
c//! Type: Double precsion
c
c//! Variable: \sigma_{mT}
c//! Description: Standard deviation of daily air temperature about \
c//!  the mean (monthly basis).
c//! Units: Degrees C             
c//! Type: Double precsion
c
c//! Variable: R
c//! Description: Random number between 0,1
c//! Units: none
c//! Type: Double precsion
c
      RETURN
      END

      
      REAL*8 FUNCTION WIND_DAY_MEAN_FN(WIND_MONTH_MEAN, 
     1   WIND_SD_MONTH, WIND_AUTOCORR_MONTH, WIND_YESTERDAY, IS1, IS2, 
     2   IS3, NORMSWITCH)
      IMPLICIT NONE
      EXTERNAL NORM_FROM_RAND_FN



c** Generates windspeed at a daily scale from monthly means, standard 
c** devaitions and autocorrelation (monthly)
      REAL*8 WIND_MONTH_MEAN, WIND_SD_MONTH
      REAL*8 WIND_AUTOCORR_MONTH, WIND_YESTERDAY
      REAL*8 NORM, NORMMN, NORMSD, NORM_FROM_RAND_FN
      INTEGER*4 IS1, IS2, IS3, NORMSWITCH
c
      NORMMN = 0.0
      NORMSD =1.0
c** approximate a function to get a normal distribution, mean=0, sd=1
      NORM = NORM_FROM_RAND_FN(NORMMN, NORMSD,IS1,IS2,IS3,NORMSWITCH)
c
      WIND_DAY_MEAN_FN = WIND_MONTH_MEAN + WIND_AUTOCORR_MONTH*
     1  (WIND_YESTERDAY - WIND_MONTH_MEAN) + WIND_SD_MONTH*
     2   NORM*( (1-(WIND_AUTOCORR_MONTH)**2)**0.5)
      WIND_DAY_MEAN_FN = MAX(WIND_DAY_MEAN_FN, 0.1d0)
c
c//! Equation: 48
c//! Description: Calculates Average daily windspeed from monthly \
c//!  means and standard-deviations. Autocorrelation of 0.65 seems \
c//!  to work if one isn't available. \\
c//!   INPUTS: Mean daily windspeed of the month; Standard deviation \
c//!    around the mean on a daily basis; windspeed  \
c//!     auto-correlation for each month; Mean windspeed of the \
c//!     previous day; 3 random number seeds.
c//! Bibliography: S Evans (SWELTER); Haith, Tubbs & Pickering (1984) \
c//!  Simulation of Pollution by soil erosion & soil nutirnt loss, \
c//!   PUDOC, Wageningen
c
c//! E: Wind_{mean} = Wind_{mean,month} + \rho_{W} (Wind_{day-1} - \
c//!     Wind_{mean,month}) + \sigma_{mW} N (1-\sigma_{mW}^{2})^0.5) \\
c//!  \\
c//!   E: N = { Random, N(0,1) }
c
c//! Variable: Wind_{mean}
c//! Description: Daily mean winsdpeed   
c//! Units: m/s 
c//! Type: Double precsion
c//! SET
c
c//! Variable: Wind_{mean,month}
c//! Description: Mean daily windspeed (monthly basis)
c//! Units: m/s
c//! Type: Double precsion
c
c//! Variable: Wind_{day-1}     
c//! Description: Mean daily windspeed of the previous day
c//! Units: m/s
c//! Type: Double precsion
c
c//! Variable: \rho_{W}
c//! Description: Observed first-order auto-correlation of mean \
c//!  daily windspeed for each month
c//! Units: none  (day/day^{-1})
c//! Type: Double precsion
c
c//! Variable: \sigma_{mW}
c//! Description: Standard deviation of daily windspeed about \
c//!  the mean (monthly basis).
c//! Units: m/s 
c//! Type: Double precsion
c
c//! Variable: R
c//! Description: Random number between 0,1
c//! Units: none
c//! Type: Double precsion
c
      RETURN
      END
     

cc commented by Ghislain 30/06/03 because it does not compile... and I don't need it   
ccc      REAL*8 FUNCTION CLOUD_COVER_FN(VP_SAT,RAINFALL_DAY, RH,
ccc     1    is1, is2, is3, pmax) 
ccc      IMPLICIT NONE
ccc
ccc
ccc
cccc**    Cloudiness for the day
ccc      REAL*8 VP_SAT, RAINFALL_DAY, RH, GAMMADIST_FN
ccc      REAL*8 RNDF, xlo ,xhi, pmax, mock_rain,mu, gamma, zero
ccc      integer*4 is1, is2, is3
cccc
ccc      CLOUD_COVER_FN = 10-1.155*SQRT(VP_SAT/MAX(RAINFALL_DAY,0.1))
cccc      CLOUD_COVER_FN = 10-2.5*SQRT(VP_SAT/MAX(RAINFALL_DAY,0.1))
ccc      CLOUD_COVER_FN = MAX(CLOUD_COVER_FN, 0.0)/10.0
ccc      if(RAINFALL_DAY.le.0.001) then
ccc        mu = rh  ! average rain on a wet day
ccc        gamma = 1.0
ccc        zero = 0.0
ccc        mock_rain = GAMMADIST_FN(zero, mu, gamma, is1, is2, is3)
ccc        CLOUD_COVER_FN = 10-2.5*SQRT(VP_SAT/MAX(mock_rain,0.1))
cccc       CLOUD_COVER_FN = 10-1.155*SQRT(VP_SAT/MAX(mock_rain,0.1))
ccc
ccc        CLOUD_COVER_FN = MAX(CLOUD_COVER_FN, 0.0)/10.0
ccc        CLOUD_COVER_FN=1*CLOUD_COVER_FN
ccc      else  ! when raining adjust
cccc     increase cloudiness - esp during summer when svp is high and raining
ccc        CLOUD_COVER_FN=1*CLOUD_COVER_FN
ccc        CLOUD_COVER_FN=min(1, (0.007*VP_SAT+0.995)*CLOUD_COVER_FN)
ccc      endif
ccc       
cccc//! Equation: 49
cccc//! Description: Proportion of the day that is cloudy \\
cccc//!  INPUTS: saturated vapour-pressure; daily rainfall
cccc//! Bibliography: S Evans {SWELTER); Nikolov & Zeller (1992) \
cccc//!  Ecological Modelling 61:149-168
cccc
cccc//! E: Cld = 10-1.155 {(VP_{sat}\over{Rain})^0.5
cccc
cccc//! Variable: Cld
cccc//! Description: Cloudiness of the day
cccc//! Units: none  (proportion)
cccc//! Type: REAL*8
cccc//! SET
cccc
cccc//! Variable: VP_{sat}
cccc//! Description: Saturated vapour pressure
cccc//! Units: Pa
cccc//! Type: REAL*8
cccc
cccc//! Variable: Rain
cccc//! Description: amount of rainfall on the day
cccc//! Units: mm
cccc//! Type: REAL*8
cccc
ccc      RETURN
ccc      END


      INTEGER*4 FUNCTION DAYS_IN_MONTH_FN(MM, DAYSINYEAR)
      IMPLICIT NONE

      INTEGER*4 DAYSINYEAR, MM, LEAPSHIFT
c**   returns the day number of days in a month (allows for leap year).
c

	LEAPSHIFT=DAYSINYEAR-365
	IF(MM.EQ.1) THEN                                                    
	  DAYS_IN_MONTH_FN = 31                                           
	ELSE IF(MM.EQ.2) THEN    
	  DAYS_IN_MONTH_FN = 28 + LEAPSHIFT                               
	ELSE IF(MM.EQ.3) THEN                                               
	  DAYS_IN_MONTH_FN = 31
	ELSE IF(MM.EQ.4) THEN                                                
	  DAYS_IN_MONTH_FN = 30
	ELSE IF(MM.EQ.5) THEN
	  DAYS_IN_MONTH_FN = 31
	ELSE IF(MM.EQ.6) THEN 
	  DAYS_IN_MONTH_FN = 30
	ELSE IF(MM.EQ.7) THEN     
	  DAYS_IN_MONTH_FN = 31
	ELSE IF(MM.EQ.8) THEN 
	  DAYS_IN_MONTH_FN = 31
	ELSE IF(MM.EQ.9) THEN 
	  DAYS_IN_MONTH_FN = 30
	ELSE IF(MM.EQ.10) THEN 
	  DAYS_IN_MONTH_FN = 31
	ELSE IF(MM.EQ.11) THEN 
	  DAYS_IN_MONTH_FN = 30
	ELSE IF(MM.EQ.12) THEN
	  DAYS_IN_MONTH_FN = 31
	ELSE
          DAYS_IN_MONTH_FN = 0
        ENDIF

c
c//! Equation: 50
c//! Description: finds the number of days in a month \\
c//!   INPUTS: Month of year (numeric); Number of days in the year
c//! Bibliography: No Specific reference to this function.
c
c//! E:   None   
c   
c//! Variable: DinM      
c//! Description:   Number of days in the month 
c//! Units: none         
c//! SET
c//! Type: Integer          
c
c//! Variable: mm      
c//! Description:  month number 
c//! Units: none         
c//! Type: Integer         
c
c//! Variable: DinY     
c//! Description: Days in the year (allows for leap year)
c//! Units: none         
c//! Type: Integer         
c
      RETURN 
      END 


      REAL*8 FUNCTION RADIAN_FROM_DEGREE_FN(DEGREES, PI)
      IMPLICIT NONE

c**   360 degrees is 2*PI radians:
      REAL*8 DEGREES, PI
	RADIAN_FROM_DEGREE_FN = (2*PI*DEGREES)/360.0
c
c//! Equation: 51
c//! Description: converts degrees to radians \\
c//!   INPUTS: Degrees (decimal); the mathematical constant, PI
c//! Bibliography: No Specific reference to this function.
c
c//! E:   RAD = 2 \pi DEG \over{260.0)
c   
c//! Variable: RAD  
c//! Description: returned value: radian
c//! Units: radian
c//! Type: REAL*8
c//! SET
c
c//! Variable: \pi  
c//! Description: mathematical constant, pi, 3.14115...
c//! Units: radian
c//! Type: REAL*8
c
c//! Variable: DEG  
c//! Description: number of degrees to be converted
c//! Units: Degrees
c//! Type: REAL*8
c
      RETURN
      END


      SUBROUTINE NORTHEAST_TO_LATLONG_SUB(EASTING, NORTHING, LATITUDE,
     1   LONGITUDE, PI)
      IMPLICIT NONE   

c** Uses the National Grid scale factors, true origins etc (after     
c** Airy (1830). These are not the same os the Irish National grid etc
c** though the formaule are the same (just different values)          
       REAL*8 NORTHING, EASTING, LATITUDE, LONGITUDE, F0
       REAL*8 P0, L0, N0, E0, A, B, PI
       REAL*8 P0_RAD, L0_RAD, R, ETASQ, V, PX, PXNEW, ESQ
       REAL*8 I7, I8, I9, I10, I11, I12, I12A , M,N
c   
c  set the constants for the OS36 UK national grid settings:
       A = 6377563.396 
       B = 6356256.910 
       F0 = 0.9996012717 
       P0 = 49.0     ! 49 deg N  True origin
       L0 = -2.0     !  2 deg W  True origin
       E0 = 400000.0
       N0 = -100000.0
c convert true origins to radians 
      P0_RAD = P0*PI/180.0       
      L0_RAD = L0*PI/180.0    
c
c Eccentricity squared:                                                
      ESQ = (A**2 - B**2)/(A**2)                                        
c Intermediates:                                                      
      N = (A-B)/(A+B)          
c                                                                     
c PX new is the estimate for the the Latitude
      PX = ((NORTHING-N0)/(A*F0))+P0_RAD                             
c Initialise M for first iteration, then continue iterations until
c desired accuracy is met (1mm)
      M = -99999999999999.0
      DO WHILE ((NORTHING-N0-M).GE.1)
	M = B*F0*( (1 +N +5/4.*N**2 +5/4.*N**3)*(PX-P0_RAD)-             
     1  (3*N +3*N**2 +21/8.*N**3) *SIN(PX-P0_RAD) *COS(PX+P0_RAD) +  
     2  (15/8.*N**2 +15/8.*N**3)*SIN(2*(PX-P0_RAD))*COS(2*(PX+P0_RAD))
     3  -(35/24.*N**3) *SIN(3*(PX-P0_RAD)) *COS(3*(PX+P0_RAD))  )
c                                                                      
	PXNEW = (NORTHING-N0-M)/(A*F0) + PX                          
	PX = PXNEW                                                   
      ENDDO
c                                                                        
      V = A*F0*(1-ESQ*SIN(PX)**2)**(-0.5)                              
      R = A*F0*(1-ESQ)*(1-ESQ*SIN(PX)**2)**(-1.5)                      
      ETASQ = V/R -1.0
c                                                                        
      I7 = DTAN(PX)/(2*R*V)                                             
      I8 = DTAN(PX)/(24*R*V**3)*(5+3*TAN(PX)**2+ETASQ -9*(TAN(PX)**2)*
     1      ETASQ**2)                                                  
      I9 = TAN(PX)/(720*R*V**5) *(61 +90*TAN(PX)**2 +45*(TAN(PX)**4))
      I10 = 1/COS(PX)/V                                                
      I11 = 1/COS(PX)/(6*V**3) * (V/R +2*TAN(PX)**2)                   
      I12 = 1/COS(PX)/(120*V**5) *(5 +28*TAN(PX)**2 +24*TAN(PX)**4)    
      I12A = 1/COS(PX)/(5040*V**7) *(61 +622*TAN(PX)**2 +              
     1       1320*TAN(PX)**4 +720*TAN(PX)**6)                          
      LATITUDE = PX-I7*(EASTING-E0)**2 +I8*(EASTING-E0)**4             
     1          -I9*(EASTING-E0)**6                                    
      LONGITUDE = L0_RAD +I10*(EASTING-E0) - I11*(EASTING-E0)**3 +     
     1          I12*(EASTING-E0)**5 - I12A*(EASTING-E0)**7             
c
c//! Equation: 52
c//! Description: Routine to take Easting and Northings to convert\
c//!  to 'conventional' Latitude and longitude (in radians). The \
c//!  contsnts used are specifically set for the UK OS36 system\
c//!  This routine is believed to be correct and is based on the OS\
c//!  manual - it should be noted that accuarcy is Unknown- the best\
c//!  guess is to about 20m. \\
c//!   INPUTS: Easting; Northing ; Latitude; Longitude; Pi\
c//! Bibliography: A Guide to coordinate systems in Great Britain,\
c//!  An introduction to mapping systems and the use of GPS datasets\
c//!  with Ordnance Survey mapping (1999), Ordnance Survey, Southampton.
c
c//! E: See the  reference - far to complicated!!
c
c//! Variable: Latitude
c//! Description: Returned Latitude value 
c//! Units:  Radians (North)
c//! Type: REAL*8
c//! SET
c
c//! Variable: Longitude
c//! Description: Returned Longitude value 
c//! Units:  Radians (East; for West, negative number returned)
c//! Type: REAL*8
c//! SET
c
c//! Variable: Easting
c//! Description: Easting coordinate of the OS National Grid
c//! Units: m 
c//! Type: REAL*8
c
c//! Variable: Northing
c//! Description: Northing coordinate of the OS National Grid
c//! Units: m 
c//! Type: REAL*8
c
      RETURN                                                           
      END                                                              


      REAL*8 FUNCTION DAYDATA_FROM_MONTHDATA_FN(MONTHDATA,
     1   DAY_OF_MONTH)
      IMPLICIT NONE

c** obtains a single element from and array containing a months worth
c** of data.
      REAL*8 MONTHDATA(31)
      INTEGER*4 DAY_OF_MONTH
c
      DAYDATA_FROM_MONTHDATA_FN = MONTHDATA(DAY_OF_MONTH)
c
c//! Equation: 53
c//! Description: Simple 'get' routine, to extract a single bit of\
c//!  data from an array with a months worth of data. \\
c//!  INPUTS: Array with the month's data; Day of the month
c//! Bibliography: none
c
c//! E: DayData = MonthData(dd)
c
c//! Variable: DayData
c//! Description: returned element of the array
c//! Units: Same as contents of MonthData
c//! Type: REAL*8
c//! SET
c
c//! Variable: MonthData
c//! Description: Array of size 31 with data for a month in it
c//! Units: whatever the data is
c//! Type: REAL*8
c
c//! Variable: dd
c//! Description: Day of the month - index for referencing the element\
c//!   in the array
c//! Units: None  (Day number of month)
c//! Type: REAL*8
c
c     RETURN
      END


      REAL*8 FUNCTION HUMIDITY_FN(TEMP_AIR_AMPLITUDE)
      IMPLICIT NONE

c** humidity based on the temperature max-min amplitude
      REAL*8 TEMP_AIR_AMPLITUDE, RHFUN
c
      RHFUN = (110.3-(4.68*TEMP_AIR_AMPLITUDE))
c       RHFUN = 77.4-3.289*TEMP_AIR_AMPLITUDE+27.357
      IF(RHFUN.LE.5) THEN
	HUMIDITY_FN = 5.0
      ELSE IF(RHFUN.GE.100) THEN
	HUMIDITY_FN = 100.0
      ELSE
	HUMIDITY_FN = RHFUN
      ENDIF
c
c//! Equation 54
c//! Description: Returns an estimate of the daily humidity based \
c//!  on the difference between Tmax and Tmin (the amplitude) \\
c//!   INPUTS: amplitude os the temperature difference (Tmax - Tmin)
c//! Bibliography: S Evans (SWELTER)
c
c//! Case: RH_x <= 5
c//!   E: RH = 5.0
c//! Case: RH_x >=100
c//!   E: RH = 100
c//! Case: 5 < RH_x <100
c//!   E: RH = RH_x  \\
c//!  \\
c//! E: RH_x = (110.3-(4.68*\Delta T_{air}))
c
c//! Variable: RH
c//! Description: Relative humidity
c//! Units: %
c//! Type: REAL*8
c//! SET
c
c//! Variable: \Delta T_{air}
c//! Description: Amplitude of temperature variation during the day
c//! Units: \deg C
c//! Type: REAL*8
c
      RETURN
      END


      SUBROUTINE GET_GRIDSQ_SUB(GRIDSIZE, LATITUDE, LONGITUDE, 
     1  GRIDLAT, GRIDLONG, PI, LONG_OFFSET)
      IMPLICIT NONE

c** This routine gets a grid square reference point from latitude and 
c** longitude it assumes that the grid reference point is at the center
c** of the grid square, and therefore any location within the square
c** refers back to the central point. Incomming Latitude and LOngitude
c** is assumed to be in radians, and westerly longitude denoted by
c** use of a negative number. The parameter LONG_OFFSET is used to 
c** to a 360 degree basis if necessary. Ie. 1 deg West is expected as 
c** -1 deg (in radians). If LONG_OFFSET is zero, then the negative 
c** aspect is maintained, if it is 360, then 1 deg W becomes 359 deg (E)
c** At the moment this deals only with the Northern Hemisphere
      REAL*8 GRIDSIZE, LATITUDE, LONGITUDE, GRIDLAT, GRIDLONG
      REAL*8 PI, LONG_OFFSET, SIGN, LATDEG, LONGDEG
c
c convert lat * long to degrees (easier reference)
      LATDEG = LATITUDE*180.0/PI
      LONGDEG = LONGITUDE*180.0/PI
c
      GRIDLAT = INT(LATDEG/GRIDSIZE)*GRIDSIZE+GRIDSIZE/2.0
      SIGN = 1
      IF(LONGDEG.lt.0) SIGN = -1
      GRIDLONG = INT(LONGDEG/GRIDSIZE)*GRIDSIZE+GRIDSIZE/2.0*SIGN+
     1    LONG_OFFSET
c       print*, latdeg, gridlat, longdeg, gridlong
c
c//! Equation: 55
c//! Description: from a given lat-long reference for the Northern  \
c//!  hemisphere, (negative numbers for Westerly longitude), this \
c//!  pinpoints a reference grid point, assuming uniform grid squares\
c//!  centered at the grid reference point. The value returned is in \
c//! degrees. a value offset can be used to deal with E/W longitude\
c//! in slightly different ways. (see description of offset variable)\\
c//!  INPUTS: Grid-square size; latitude of site (radians); \
c//!    longitude of site (radians);\
c//!   Grid-reference latitude [output]; Grid-reference longitude \
c//!   [output]; Constant, Pi, longitude offset
c//! Bibliography: T Randle (2000) - seemed like a good idea at the time
c
c//! E: GR_{lat} = INT(LAT/GR_{size}) GR_{size} + GR_{size)/2.0 \\
c//! \\
c//! E: GR_{long} = INT(LONG/GR_{size}) Gr_{size} + \
c//!      GR_{size) S/2.0 + LAT_{off}  \
c//!    S = 1 if LONG >=0;      Else S=-1 if LONG <0
c
c//! Variable: GR_{lat}
c//! Description: Grid-point center latitude value 
c//! Units: Degrees 
c//! Type: REAL*8
c//! SET
C
c//! Variable: GR_{long}
c//! Description: Grid-point center longitude value 
c//! Units: Degrees 
c//! Type: REAL*8
c//! SET
c
c//! Variable: LAT
c//! Description: Latitude of the site (Northern Hemishpere only) - in \
c//!  the above equations LAT has been transformed into degrees, though\
c//!  the input to the routine is radians
c//! Units: Degrees (in calculation); (radian in input)
c//! Type: REAL*8
c
c//! Variable: LONG
c//! Description: Longitude of the site (Northern Hemishpere only) - in\
c//!  the above equations LONG is transformed into degrees, though\
c//!  the input to the routine is radians
c//! Units: Degrees (in calculation); (radian in input)
c//! Type: REAL*8
c
c//! Variable: GR_{size}
c//! Description: Size of the grid-square
c//! Units: degrees
c//! Type: REAL*8
c
c//! Variable:  LAT_{off}
c//! Description: Offset value to deal with different formats of \
c//!  longitude. If Westerly values are needed as negative, then \
c//!  LAT_{off} should be 0; if degrees from east are required, then\
c//!  it should have the value 360
c//! Units: degrees
c//! Type: REAL*8
c
      RETURN
      END


      REAL*8 FUNCTION VP_UNSAT2_FN(VP_SAT, HUMIDITY)
      IMPLICIT NONE

c**   reverses the humidity function and derives the vapour pressure 
c**  (unsaturated) from the saturated VP (at the current temp) and 
c**   humidity
      REAL*8 VP_SAT, HUMIDITY
c
      VP_UNSAT2_FN = HUMIDITY*VP_SAT/100.0
c
c//! Equation: 56
c//! Description: Reverses the humidity function and derives the \
c//!  vapour pressure unsaturated) from the saturated VP (at the \
c//!  current temp) and humidity \\
c//!   INPUTS: Saturated vapour pressure; humidity
c//! Bibliogrpahy: none
c
c//! E: VP_{unsat} = {RH VP_{sat}}\over {100.0}
c
c//! Variable: VP_{unsat} 
c//! Description: unsaturated vapour pressure
c//! Units: mbar (or same units as VP_{sat}
c//! Type: REAL*8
c//! SET
c
c//! Variable: VP_{sat} 
c//! Description: saturated vapour pressure
c//! Units: mbar  (or other eg kPa - will change VP_{unsat} units)
c//! Type: REAL*8
c
c//! Variable: RH
c//! Description: Relative humidity
c//! Units: %
c//! Type: REAL*8
c
      RETURN
      END 


      REAL*8 FUNCTION RAD_TO_RADPAR_FN(RAD, RAD_PERC_PAR)
      IMPLICIT NONE

c** function to convert 'global' radiation to the radiation falling
c** within the photosynthetically active wavelengths; circa 375-800nm
c** the value for this may vary somewhere between 35% and 60% - a 
c** general approximation seems to be around 45%
      REAL*8  RAD, RAD_PERC_PAR
c
      RAD_TO_RADPAR_FN = RAD*RAD_PERC_PAR/100.0
c
c//! Equation: 57
c//! Description: Converts 'global' radiation to the radiation \
c//!  falling in the photosynthetically active wavelengths; circa \
c//!  375-800nm \\
c//!   INPUTS: Total light; percentage of total light that is \
c//!    photosynthetically active
c//! Bibliography: S Evans (SWELTER); Salisbury & Ross (1992), \
c//!  Plant Physiology, Wadsworth Publishing, California (pp613)
c
c//! E: RAD_{PAR} = RAD RAD_{perc, PAR}/100.0
c
c//! Variable: RAD_{PAR} 
c//! Description: Photosynthetically active radiation (W/m2)
c//! Units: W/m^2
c//! Type: REAL*8
c//! SET
c
c//! Variable: RAD
c//! Description: Total radiation (W/m2)
c//! Units: W/m^2
c//! Type: REAL*8
c
c//! Variable: RAD_{perc,PAR}
c//! Description: Percentage of total radiation that is in the PAR range
c//! Units: none (percentage)
c//! Type: REAL*8
c
      RETURN
      END


      REAL*8 FUNCTION RADPAR_TO_RAD_FN(RADPAR, RAD_PERC_PAR)
      IMPLICIT NONE

c** function to convert radiation falling
c** within the photosynthetically active wavelengths; circa 375-800nm
c** to total radiation 
c** the value for this may vary somewhere between 35% and 60% - a 
c** general approximation seems to be around 45%
      REAL*8 RADPAR, RAD_PERC_PAR
c
      RADPAR_TO_RAD_FN = RADPAR/RAD_PERC_PAR*100.0
c
c//! Equation: 58
c//! Description: Converts PAR radiation to the total radiation \\
c//!   INPUTS: photosynthetically active; percentage of total light \
c//!    that is photosynthetically active
c//! Bibliography: S Evans (SWELTER); Salisbury & Ross (1992), \
c//!  Plant Physiology, Wadsworth Publishing, California (pp613)
c
c//! E: RAD = RAD_{PAR} 100.0 / RAD_{perc, PAR}
c
c//! Variable: RAD
c//! Description: Total radiation (W/m2)
c//! Units: W/m^2
c//! Type: REAL*8
c//! SET
c
c//! Variable: RAD_{PAR}
c//! Description: Photosynthetically active radiation (W/m2)
c//! Units: W/m^2
c//! Type: REAL*8
c
c//! Variable: RAD_{perc,PAR}
c//! Description: Percentage of total radiation that is in the PAR range
c//! Units: none (percentage)
c//! Type: REAL*8
c
      RETURN
      END


      REAL*8 FUNCTION RADPAR_TO_RADPPFD_FN(RADPAR, 
     1  PPFDCONV)
      IMPLICIT NONE

c** function to convert PAR radiation to photosynthtic flux density
c** general approximation seems to be around: 1W (PAR) = 4.5 umol PAR
      REAL*8 RADPAR, PPFDCONV
c
      RADPAR_TO_RADPPFD_FN = RADPAR*PPFDCONV
c
c//! Equation: 59
c//! Description: Converts PAR radiation to the PPFD PAR \\
c//!   INPUTS: PAR radiation; conversion of PAR light to PPFD 
c//! Bibliography: S Evans (SWELTER); Salisbury & Ross (1992), \
c//!  Plant Physiology, Wadsworth Publishing, California (pp613)
c
c//! E: RAD_{PPFD} = RAD_{PAR} PPFD_{conv}
c
c//! Variable: RAD_{PPFD} 
c//! Description: Photosynthetic flux density  (umol/m2)
c//! Units: umol (PAR) /m^2
c//! Type: REAL*8
c//! SET
c
c//! Variable: RAD_{PAR}
c//! Description: Photosynthetically active radiation (W/m2)
c//! Units: W/m^2
c//! Type: REAL*8
c
c//! Variable: PPFD_{conv}
c//! Description: Conversion factor for convertinq W (PAR) to umol (PAR)
c//! Units: umol(PAR)/m^2 / W(PAR)/m^2
c//! Type: REAL*8
c
      RETURN
      END


      REAL*8 FUNCTION RADPPFD_TO_RADPAR_FN(RADPPFD, 
     1  PPFDCONV)
      IMPLICIT NONE

c** function to convert photosynthtic flux density to PAR radiation 
c** general approximation usess around: 1W (PAR) = 4.5 umol PAR
      REAL*8 RADPPFD, PPFDCONV
c
      RADPPFD_TO_RADPAR_FN = RADPPFD/PPFDCONV
c
c//! Equation: 60
c//! Description: Converts PPFD PAR to PAR radiation\\
c//!   INPUTS: PPFD PAR ; conversion of PAR light to PPFD 
c//! Bibliography: S Evans (SWELTER); Salisbury & Ross (1992), \
c//!  Plant Physiology, Wadsworth Publishing, California (pp613)
c
c//! E: RAD_{PAR} = RAD_{PPFD} PPFD_{conv}
c
c//! Variable: RAD_{PAR}
c//! Description: Photosynthetically active radiation (W/m2)
c//! Units: W/m^2
c//! Type: REAL*8
c//! SET
c
c//! Variable: RAD_{PPFD} 
c//! Description: Photosynthetically active flux density
c//! Units: umol (PAR) /m^2
c//! Type: REAL*8
c
c//! Variable: PPFD_{conv}
c//! Description: Conversion factor for convertinq W (PAR) to umol (PAR)
c//! Units: umol(PAR)/m^2 / W(PAR)/m^2
c//! Type: REAL*8
c
      RETURN
      END


      REAL*8 FUNCTION RAD_DIFF_ABOVE_FN(DIFFPROP_DAY, 
     1  RAD_TERRA_DAY, ALBEDO_DIFF)
      IMPLICIT NONE

c** works out the amount of terrestrial radiation that hits a canopy
c** that is diffuse - this is after albedo reflection has occurred.
      REAL*8 RAD_TERRA_DAY, ALBEDO_DIFF
      REAL*8 DIFFPROP_DAY
c
      RAD_DIFF_ABOVE_FN = RAD_TERRA_DAY*DIFFPROP_DAY*(1-ALBEDO_DIFF)
c
c//! Equation: 61
c//! Description: The amount of terrestrial radiation above the \
c//!  canopy that is diffuse - after albedo reflection \\
c//!   INPUTS: proportion of terrestrial radiation that is diffuse;\
c//!    terrestrial radiation; Albedo (reflection) for diffuse light
c//! Bibliography: none
c
c//! E: RAD_{abv,D} = RAD_{terra} Diff_{prop} (1-ALB_{D})
c
c//! Variable: RAD_{abv,D}
c//! Description: Diffuse radiation above the canopy after albedo\
c//!  reflection
c//! Units: W/m^2/day
c//! Type: REAL*8
c//! SET
c
c//! Variable: RAD_{terra}
c//! Description: Terrestrial radiation 
c//! Units: W/m^2/day
c//! Type: REAL*8
c
c//! Variable: Diff_{prop}
c//! Description: proportion of Terrestrial radiation that is diffuse 
c//! Units: none (proportion)
c//! Type: REAL*8
c
c//! Variable: ALB_{D}
c//! Description: Albedo (reflection) coefficient for diffuse light 
c//! Units: none (proportion)
c//! Type: REAL*8
c
      RETURN
      END

      REAL*8 FUNCTION RAD_BEAM_ABOVE_FN(DIFFPROP_DAY, 
     1  RAD_TERRA_DAY, ALBEDO_BEAM)
      IMPLICIT NONE

c** works out the amount of terrestrial radiation that hits a canopy
c** that is direct (beam)- this is after albedo reflection has occurred.
      REAL*8 RAD_TERRA_DAY, ALBEDO_BEAM
      REAL*8 DIFFPROP_DAY
c
      RAD_BEAM_ABOVE_FN = RAD_TERRA_DAY*(1-DIFFPROP_DAY)*(1-ALBEDO_BEAM)
c
c//! Equation: 62
c//! Description: The amount of terrestrial radiation above the \
c//!  canopy that is direct (beam) - after albedo reflection \\
c//!   INPUTS: proportion of terrestrial radiation that is diffuse;\
c//!    terrestrial radiation; Albedo (reflection) for beam light
c//! Bibliography: none
c
c//! E: RAD_{abv,B} = RAD_{terra} (1-Diff_{prop}) (1-ALB_{B})
c
c//! Variable: RAD_{abv,B}
c//! Description: Direct (beam) radiation above the canopy after albedo\
c//!  reflection
c//! Units: W/m^2/day
c//! Type: REAL*8
c//! SET
c
c//! Variable: RAD_{terra}
c//! Description: Terrestrial radiation 
c//! Units: W/m^2/day
c//! Type: REAL*8
c
c//! Variable: Diff_{prop}
c//! Description: proportion of Terrestrial radiation that is diffuse 
c//! Units: none (proportion)
c//! Type: REAL*8
c
c//! Variable: ALB_{B}
c//! Description: Albedo (reflection) coefficient for direct (beam)light 
c//! Units: none (proportion)
c//! Type: REAL*8
c
      RETURN
      END

      


c----------------------------------------------------------

      SUBROUTINE RAIN_DURATION2_SUB(RAIN, RAIN_DURATION, RAIN_INTENSITY,
     1  DAYS_IN_MONTH, IS1, IS2, IS3)
	   IMPLICIT NONE



c** This routine estimates the duration and intensity of rain on each day 
c** for a month Accuracy is suspect.
      REAL*8 RAIN(31), RAIN_DURATION(31), RAIN_INTENSITY(31)
      REAL*8 RAIN_LIM(10), TIMESCALE, FREQUENCY, PROB, RNDF
      REAL*8 XLO, XHI
      INTEGER*4 DAYS_IN_MONTH, IDAY, CLASSNO, MAX_RAINTIME, I 
      INTEGER*4 IS1, IS2, IS3
c                        
      RAIN_LIM(1) = 0    
      RAIN_LIM(2) = 5    
      RAIN_LIM(3) = 10   
      RAIN_LIM(4) = 15   
      RAIN_LIM(5) = 20   
      RAIN_LIM(6) = 25   
      RAIN_LIM(7) = 50   
      RAIN_LIM(8) = 75   
      RAIN_LIM(9) = 100  

c** find which rain class the current preciptation is in
      DO 10 IDAY = 1, DAYS_IN_MONTH
	IF(RAIN(IDAY).GT.0) THEN
	  IF(RAIN(IDAY).GT.RAIN_LIM(9)) THEN
	    CLASSNO = 9
	  ELSE
	    DO 20 I = 9,2,-1
	      IF(RAIN(IDAY).LE.RAIN_LIM(I)) CLASSNO = I
20          CONTINUE
	  ENDIF
c**       The maximum duration of a rain instance is
	  MAX_RAINTIME=MIN(24,INT(.5*EXP(0.1386*RAIN_LIM(CLASSNO))+.5))
c
c**       The cumulative number of rainy days of a class type per century 
c**       is proportional to the rainfall duration. 
c**       (ie: Rainy days of this magnitude per century is timescale*duration)
c
	  TIMESCALE = 1.39*((RAIN_LIM(CLASSNO)/25.4+0.1)**(-3.55))*10
c
c**       thus the cumulative instances of rain of this magnitude (per century)is

	  FREQUENCY = TIMESCALE*MAX_RAINTIME
c
c**       the probability of rain falling at some stage within time, t, is
c         prob = timescale*t/frequency; thus if we get a random number between 
c**       0.025 and 1 (ignoring the chances of P<0.025; then the rain 
c**       duration (minutes) is
c
	  XLO = 0.025
	  XHI = 1.0
	  PROB = RNDF(XLO, XHI, IS1, IS2, IS3)
	  RAIN_DURATION(IDAY) = PROB/TIMESCALE*FREQUENCY*60
c
c**       the intensity (mm/hour) is
	  RAIN_INTENSITY(IDAY) = 60.0/RAIN_DURATION(IDAY)*RAIN(IDAY)
	ELSE 
	  RAIN_DURATION(IDAY)=0.0
	  RAIN_INTENSITY(IDAY)=0.0
	ENDIF
10    CONTINUE
      RETURN
      END



      REAL*8 FUNCTION VP_SAT2_FN(TEMP_C)
      IMPLICIT NONE

c**   Returns the saturated air pressure for a given temperature
      REAL*8 TEMP_C
c
      IF(TEMP_C.GE.0 ) THEN
	VP_SAT2_FN = 610.78*exp(17.269*TEMP_C/(TEMP_C+237.3))
      ELSE
	     VP_SAT2_FN = exp(-6140.4/(273+TEMP_C)+28.916)
	   Endif
c
c//! Equation: 71
c//! Description: Returns saturated vapour pressure for a temperature\\
c//!   INPUTS: Temperature
c//! Bibliography: http://www.natmus.dk/cons/tp/atmcalc/atmoclc1.htm  (26.01.2001)
c
c//! Case: if T>=0
c//!  E:    VP_{sat} = 610.78*exp(T\over{(T+238.3)}*17.2694)   
c//! Case: t<0 
c//!  E:    VP_{sat} = exp(-6140.4\over{(273+T)}+28.916)       
c
c//! Variable: VP_{sat}
c//! Description: saturated vapour pressure at a given temperature\\
c//!  typical ranges are 300-1200 Pa: 0.3-1.2kpa: 3-12mbar
c//! Units: Pa
c//! SET
c//! Type: REAL*8
c
c//! Variable: T
c//! Description: Air temperature 
c//! Units: Degrees  C 
c//! Type: REAL*8
c
       RETURN
       END


      REAL*8 FUNCTION TRANS_PROP_DAY2_FN(TEMP_AIR_AMPLITUDE,
     1  A_SKY, ANGSTROM_BETA, C_SKY) 
      IMPLICIT NONE

c** determines aproximate transmittance of radiatio to the earth
      REAL*8 TEMP_AIR_AMPLITUDE
      REAL*8 A_SKY, ANGSTROM_BETA, C_SKY

      IF (TEMP_AIR_AMPLITUDE.LE.0) THEN
	      TRANS_PROP_DAY2_FN = ANGSTROM_BETA
      ELSE
	      TRANS_PROP_DAY2_FN = ANGSTROM_BETA*(1-
     1    EXP(-A_SKY*TEMP_AIR_AMPLITUDE**C_SKY))
      ENDIF
      
c
c//! Equation: 72
c//! Description: Works out the transmittance of radiation through the\
c//!   atmosphere - a reverse engineered approximation of eqn 41. Untested
c//!   and certainly dubious at low temperature amplitudes as we fix it to\
c//!   be Angstrom_beta; ehreas the original allows for transmittance to be\
c//!   >  Angstrom_beta\\
c//!   INPUTS: Air temperature amplitude (max-min);    \
c//!    Coefficient for maximum clear-sky transmittance;  \
c//!   Angstrom \beta value; Coefficient of clear-sky transmittance \
c//!   with amplitude temperature increase. 
c//! Bibliography: S Evans (SWELTER); Bristow & Campbell (1984) Agric \
c//!  & Forest Met 31: 159-166
c
c//! Case: \Delta_T <=0.0
c//!  E: T_{sw} = \beta
c//! Case: \Detla_T > 0
c//!  E: T_{sw} = \beta(1-exp(-A_{sky}*\Delta_T^{C_{sky}})) 
c
c//! Variable: T_{sw}
c//! Description: proportion of Extra-terrestrial radiation that \
c//!  arrives on the earth's surface as short-wave.
c//! Units: none
c//! Type: REAL*8
c//! SET
c
c//! Variable: \Delta T
c//! Description: Daily amplitude of temperature variation
c//! Units: Degrees C  (usually)
c//! Type: REAL*8
c
c//! Variable: \beta
c//! Description: Angstrom tubidity factor: max clear sky transmittance
c//! Units: none  
c//! Type: REAL*8
c
c//! Variable: A_{sky}
c//! Description: Coefficient of maximum clear-sky transmittance \
c//!  characteristics
c//! Units: none  
c//! Type: REAL*8
c
c//! Variable: C_{sky}
c//! Description: Coefficient of maximum clear-sky transmittance \
c//!   with \Delta T increase
c//! Units: none  
c//! Type: REAL*8
c
      RETURN
      END


      REAL*8 FUNCTION GAMMADIST_FN(A,B,C, SEED1, SEED2, SEED3)
      IMPLICIT NONE


      REAL*8 Aconst, Bconst, Cconst, Qconst, Tconst
      REAL*8 Dconst, a, b, c, xlo, xhi, R1, R2, P, Y
      REAL*8 retval, p1, p2, v, z, w, rndf, R3
      integer*2 flag
      integer*4 seed1, seed2, seed3
      Bconst = c-log(4.0)
      Tconst = 4.5
      Dconst = 1+log(Tconst)
      Cconst = 1+c/exp(1.0)
      flag = 1
      Xlo = 0.0
      xhi = 1.0
c      print*, 'in gammadist', a,b,c, seed1, seed2, seed3

      if(c.lt.1) then
	  flag =0
100     R1 = RNDF(XLO, XHI, SEED1, SEED2, SEED3)
	  p = Cconst*R1
	  if(P.GT.1) GOTO 300
        y=p**(1/c)
	  R2 = RNDF(XLO, XHI, SEED1, SEED2, SEED3)
	  if(R2.le.exp(-y)) then
	    retval = a+b*y
	    flag = 2
	  endif
	  if(flag.eq.0) goto 100

300     y = -log((Cconst-p)/c)
	  r2 = RNDF(XLO, XHI, SEED1, SEED2, SEED3)
	  if((y.GT.0).AND.(R2.le.y**(c-1))) then
	    retval = a+b*y
	    flag = 3
        endif        
        if(flag.eq.0) goto 100
c        print *,flag,'=',retval

      elseif(c.gt.1) then
        Aconst = 1/SQRT(2*c-1.0)  
        Qconst = c+ 1/(Aconst)
	  flag=0 
20      p1 = RNDF(XLO, XHI, SEED1, SEED2, SEED3)
	  p2 = RNDF(XLO, XHI, SEED1, SEED2, SEED3)
	  v = Aconst*log(p1/(1-p1))
	  y = c*exp(v)
	  z = p1*p1*p2
	  w = Bconst+Qconst*v - y
	  if(w+Dconst-Tconst*z.ge.0) then
	    retval =a+b*y
	    flag = 3
	  else
	    if(w.ge.log(z)) then
	      retval = a+b*y
	      flag = 4
	    endif
	  endif
	  if(flag.eq.0) goto 20
      else
        r3 = RNDF(XLO,XHI, SEED1, SEED2, SEED3)
        if(b.gt.0) then
          retval = a-b*log(r3)
        else
          print*, 'problem with gamma dist when c=1'
          stop
        endif
      endif

      gammadist_fn = retval 
c//! Equation: 73
c//! Description: Derives a random value from the gamma distribution \\
c//!   INPUTS: 3 parameters describing the gamma distribution shape; 3\
c//!    Seeds for random number generation
c//! Bibliography: R Saucier (2000), Computer generation of statistical\ 
c//!  Distributions; (US) Army Research Laboratory, ARL-TR-2168.\
c//!  (http://ftp.arl.mil/random/  (17/5/2001) )
c
      return
      end

     
      REAL*8 FUNCTION NORM_FROM_RAND_FN(MEAN, SD, SEED1,SEED2,
     +     SEED3, NORMSWITCH)

      IMPLICIT NONE

      REAL*8 Mean, SD, NormVal
      REAL*8 NORMDIST_FN, NORMDIST2_FN, TRANSNORM_FN
      Integer*4 Seed1, Seed2, Seed3, NormSwitch
c//! Equation: 86
c//! Description: Generates a value from a pseudo normal distribution \\
c//!   is a shell for calling different methods. \\
c//!   INPUTS: Mean; Sd; 3 Seeds for random number generation; Method id
c//! Bibliography: None
c
      IF(NormSwitch.EQ.1) then
         NormVal = NORMDIST_FN(SEED1,
     +   SEED2,SEED3)
*         IF((MEAN.NE.Zero).and.(SD.NE.1)) THEN
         IF ((abs(MEAN).GT.0.0d0).AND.(abs(SD-1.0d0).GT.0.0d0)) THEN
c          transform the distribution
           NormVal = TRANSNORM_FN(Mean, SD, NormVal)
         ENDIF
      Else if(NormSwitch.eq.2) THEN
         NormVal = NORMDIST2_FN(MEAN,SD,SEED1,SEED2,SEED3)
      EndIF
      NORM_FROM_RAND_FN = NormVal

      Return
      END


      REAL*8 FUNCTION NORMDIST_FN(SEED1,
     +   SEED2,SEED3)
      IMPLICIT NONE

      REAL*8 xlo, xhi, R1, retval, rndf
      integer*4 seed1, seed2, seed3
 
      Xlo = 0.0
      xhi = 1.0
      R1 = RNDF(XLO, XHI, SEED1, SEED2, SEED3)

      Retval = ( R1**0.135 - (1-R1)**0.135 )/0.1975

      NORMDIST_FN = Retval
c//! Equation: 74
c//! Description: Derives a value from standard normal distribution\
c//!  (mean 0, sd 1), at random. Reputed to have accuracy of 0.5% in\
c//!  The range 0<=p<=0.9975 (NB random numbers are between 0 and 1) \\
c//!   INPUTS: 3 Seeds for random number generation
c//! Bibliography: S Evans, SWELTER; Haith, DA, Tubbs, LJ and \
c//!  Pickering, NB(1984) Simulation of Pollution by soil erosion and \
c//!  soil nutrient loss; Pudoc, Wageningen.
C
C//! E: N = { R1^{0.135} - (1-R1)^{0.135)} } \over {0.1975}
C 
c//! Variable: SEED1, SEED2, SEED3
c//! Description: Integer seeds for random number generation
c//! Units: none  
c//! Type: Integer
c
c//! Variable: R1
c//! Description: Random number between 0 and 1
c//! Units: none  
c//! Type: REAL*8
c
c//! Variable: NORMDIST_FN
c//! Description: Simulated value from a normal distribution N(0,1)
c//! Units: none
c//! Type: REAL*8
c//! SET
c

      return
      end



      REAL*8 FUNCTION NORMDIST2_FN(MEAN,SD,SEED1,SEED2,SEED3)
      IMPLICIT NONE


      REAL*8 MEAN, SD
      REAL*8 xlo, xhi, R1, R2, P, retval, rndf 
      integer*4 seed1, seed2, seed3
c
      Xlo = -1.0
      xhi = 1.0
      P = 2.0
c  
10    IF (P.GE.1.0d0) THEN 
	R1 = RNDF(XLO, XHI, SEED1, SEED2, SEED3)
	R2 = RNDF(XLO, XHI, SEED1, SEED2, SEED3)
	P = R1*R1 + R2*R2
      GOTO 10
      ENDIF

c
      Retval = MEAN + SD*R1*SQRT(-2.0*LOG(P)/P)
C     NB . . Natural Log!
c
      NORMDIST2_FN = Retval
c
C//! Equation: 75
c//! Description: Derives a random value from the Normal distribution\\
c//!   INPUTS: 2 parameters describing the distribution shape;- Mean \
c//!     And Standard deviation; 3 Seeds for random number generation
c//! Bibliography: R Saucier (2000), Computer generation of statistical\
c//!   Distributions; (US) Army Research Laboratory, ARL-TR-2168; \
c//!   (http://ftp.arl.mil/random/  (17/5/2001) )
c
c//! E: N(\mu, \sigma) = \mu + \sigma R1 \root{(-2 Ln(P)/P)}  \
c//!     Where P = R1^2 + R2^2
C 
c//! Variable: SEED1, SEED2, SEED3
c//! Description: Integer seeds for random number generation
c//! Units: none  
c//! Type: Integer
c
c//! Variable: R1, R2
c//! Description: Random number between -1 and +1
c//! Units: none  
c//! Type: REAL*8
c
c//! Variable: NORMDIST2_FN
c//! Description: Simulated value from a normal distribution \
c//!    N(\mu,\sigma)
c//! Units: none
c//! Type: REAL*8
c//! SET
c
c      print*, 'done normdist2'
      return
      end

	
*      REAL*8 FUNCTION LAT_SLOPE_FN(Incl_rad, Lat_rad, 
*     1   Aspect_rad, Pi)
*      IMPLICIT NONE

*      REAL*8 Lat_rad, Incl_rad, retval, Pi, Aspect_rad
c
*      retval = COS(Incl_rad)*SIN(Lat_rad)+
*     1   SIN(Incl_rad)*COS(Lat_rad)*COS(Aspect_rad) 
*      if(retval.lt.-1) retval = -1.0
*      if(retval.gt.1) retval = 1.0 
*      retval=ASIN(retval)
c
c//! Equation 76
c//! Description: calculates an equivalent latitude for a point, \
c//!    given a slope and inclination\\
c//!   INPUTS: Inclination (from horizontal); Latitude; Aspect; Pi
c//! Bibliography: Swift, L.W. (1976) Algorithm for Solar Radiation on \
c//!   Mountain Slopes Water Resources Research (12) 1: 108-112
c
c//! E: Lat_{slope} = ASIN(COS(Incl_{rad}*SIN(Lat_{rad}+ \
c//!          COS(Lat_{rad})*SIN(Asp_{rad}))
c
c//! Variable:  Incl_{rad}
c//! description: Inclination of the slope (from horizontal)
c//! Units: Raidans 
c//! Type: REAL*8
c
c//! Variable:  Lat_{rad}
c//! description: latitude of the site
c//! Units: Raidans 
c//! Type: REAL*8
c
c//! Variable:  Asp_{rad}
c//! description: Aspect (Azimuth) of the slope (from North)
c//! Units: Raidans 
c//! Type: REAL*8
c
c//! Variable: Pi
c//! Description: mathematical constant, pi
c//! Units: none
c//! Type: REAL*8
c
c//! Variable: Lat_{slope}
c//! Description: Latitude of an equivalent point, accounting for \
c//!    the slope & aspect
c//! Units: radians
c//! Type: REAL*8
c//! SET
c
*      LAT_SLOPE_FN = retval
*      RETURN
*      END


      REAL*8 FUNCTION RAD_SLOPE_FN(Declin_rad, Lat_rad, 
     1    Lat_slope_rad, Inclin_rad, Aspect_rad, Pi)

      IMPLICIT NONE


      REAL*8 Declin_rad, Lat_rad, Lat_slope_rad, Inclin_rad
      REAL*8 Aspect_rad, Pi
      REAL*8 DAWN_HOURANGLE_FN
      REAL*8 Time_Offset, T, T6, T7, T1, T0, T3, T2, T8, T9
      REAL*8 Rad0, Rad1, Rad_horiz, Adj

c    calculate an offset for time for the actual and equivalent slopes
      Time_Offset = ATAN( SIN(Inclin_rad)*SIN(Aspect_rad) /
     1  ( COS(Inclin_rad)*COS(Lat_rad) -
     2   SIN(Inclin_rad)*SIN(Lat_rad)*COS(Aspect_rad)) )*Pi/180
c    calculate new 'adjusted' dawn/dusk hour angles
c      print*, 'calling dawn:',Declin_rad, Lat_slope_rad, lat_rad
      T = DAWN_HOURANGLE_FN(Declin_rad, Lat_slope_rad)
      T7 = T-Time_Offset
      T6 = -T-Time_Offset
c    claclulate un-adjusted hour-angles
      T = DAWN_HOURANGLE_FN(Declin_rad, Lat_rad)
      T1 = T
      T0 = -T
      If(T7.LT.T1) THEN
	T3 = T1
      Else
	T3 = T7
      Endif
      If(T6.GT.T0) THEN
	T2 = T6
      Else
	T2 = T0
      Endif
c
c    This bit supposedly guards against double sunrises and 90 degree latitudes 
c   UNTESTED!!!
c    Calculate the potential on horizontal surface
      If(T2.GE.T3) THEN
	T2 = 0
	T3 = 0
      ENDIF
      T6 = T6+Pi
      Rad0 = 0.0
      If(T6.Ge.T1) Then
	T8 = T6
	T9 = T1
      Else
	T7 = T7 - Pi
	If(T7.Gt.T0) Then
	  T8 = T0
	  T9 = T7
	  Rad0 = SIN(Declin_rad)*SIN(Lat_slope_rad)*(T9-T8)*12/Pi + 
     1         COS(declin_rad)*COS(Lat_slope_rad)*
     2         ( SIN(T9+Time_Offset)-SIN(T8+Time_Offset))*12/Pi
	Endif
      Endif
      Rad1 = SIN(Declin_rad)*SIN(Lat_slope_rad)*(T3-T2)*12/Pi + 
     1         COS(declin_rad)*COS(Lat_slope_rad)*
     2      ( SIN(T3+ Time_Offset)-SIN(T2+ Time_Offset))*12/Pi
      Rad1= Rad1+Rad0
c     Print*, declin_rad, lat_slope_rad, t2, t3
c
c    Calculate the Potential on a horizontal surface
      Adj = 0.0
      Rad_Horiz = SIN(Declin_rad)*SIN(Lat_rad)*(T1-T0)*12/Pi + 
     1         COS(declin_rad)*COS(Lat_rad)*
     2         ( SIN(T1+Adj)-SIN(T0+Adj))*12/Pi
c      Print*, declin_rad, lat_rad, t0, t1
      RAD_SLOPE_FN = Rad1/Rad_Horiz
c      Print*, Rad1, Rad_Horiz
 

c//! Equation 77
c//! Description: Calculates the ratio between horizontal surface \
c//!   Direct radiation and on a slope it finds an equivalent latitude \
c//!   for daylength purposes, but uses the correct declinations and \
c//!   elevations. It appears that the day is still spread evenly about \
c//!   noon! - this may account for the non-night light that occurs before \
c//!   the 'new sunrise/set
c//! Bibliography: Swift, L.W. (1976) Algorithm for Solar Radiation on \
c//!   Mountain Slopes. Water Resources Research (12) 1: 108-112
c//! E:  Not even going to try - see the paper!!
c
      RETURN
      END


      SUBROUTINE RAD_SLOPE_DAY_SUB(RAD_BEAM_SLOPE, RAD_DIFF_SLOPE,
     1   RAD_TERRA_DAY, DIFFPROP_DAY, SLOPE_RAD, 
     1   ASPECT_RAD, DECLINATION_RAD, LAT_RAD, SUNRISE, PI) 
	   IMPLICIT NONE

      REAL*8 RAD_BEAM_SLOPE, RAD_DIFF_SLOPE
      REAL*8 RAD_TERRA_DAY,DIFFPROP_DAY,SLOPE_RAD,ASPECT_RAD
      REAL*8 DECLINATION_RAD, lAT_RAD, SUNRISE, PI
      REAL*8 CT1, SHF, SMP, W, CT2, CT3, CT, CTS, CTZ, ST
      REAL*8 CTZS, TFR, DIFFR
      integer*4 SH,i
c
c     calculates the adjustment for beam and diffuse light intercepted on
c     a sloping surface. No account is taken of the difference of length
c     of the surface compared to the flat plane!! It also takes no account
c     of increased shading beyond the slope (if facing away from the sun)-
c     we know nothing about surrounding terrain.
	  
      CT1 = SIN(DECLINATION_RAD)*(SIN(LAT_RAD)*COS(SLOPE_RAD)-
     1      COS(LAT_RAD)*SIN(SLOPE_RAD)*COS(ASPECT_RAD))
      SH = INT(SUNRISE+1)      
      SHF = SH - SUNRISE
      SMP = SUNRISE+SHF/2.0
      W = 15*(SMP-12)*PI/180.0
      CT2 = COS(DECLINATION_RAD)*COS(W)*(COS(LAT_RAD)*COS(SLOPE_RAD)+
     1    SIN(LAT_RAD)*SIN(SLOPE_RAD)*COS(ASPECT_RAD))
      CT3 = COS(DECLINATION_RAD)*SIN(SLOPE_RAD)*SIN(ASPECT_RAD)*SIN(W)
      CT = CT1+CT2+CT3
      CTS = CT*SHF
      CTZ = COS(LAT_RAD)*COS(DECLINATION_RAD)*COS(W)+
     1                SIN(LAT_RAD)*SIN(DECLINATION_RAD)
      CTZS = CTZ*SHF
      TFR = CT/CTZ
      
*      DO 10 ST = SH+0.5, 11.5
      DO i=0,11
        ST = SH+i+0.5
        W = 15*(ST-12)*PI/180.0
        CT2 = COS(DECLINATION_RAD)*COS(W)*(COS(LAT_RAD)*COS(SLOPE_RAD)+
     1    SIN(LAT_RAD)*SIN(SLOPE_RAD)*COS(ASPECT_RAD))
        CT3 = COS(DECLINATION_RAD)*SIN(SLOPE_RAD)*SIN(ASPECT_RAD)*SIN(W)
        CT = CT1+CT2+CT3
        CTS = CTS+CT
        CTZ = COS(LAT_RAD)*COS(DECLINATION_RAD)*COS(W)
     1                +SIN(LAT_RAD)*SIN(DECLINATION_RAD)
        CTZS = CTZS+CTZ
        tfr = ct/ctz
        
      ENDDO
*10    CONTINUE
      TFR = CTS/CTZS
c
c     now diffuse
      DiffR = COS(SLOPE_RAD/2)*COS(SLOPE_RAD/2) 
     
      RAD_BEAM_SLOPE = TFR*RAD_TERRA_DAY*(1-DIFFPROP_DAY)
      RAD_DIFF_SLOPE = DiffR*(RAD_TERRA_DAY*DIFFPROP_DAY)

c//! Equation 77
c//! Description: modifies incomming light to account for a slope and aspect\
c//!  effect. Sin in the northern hemispere, north facing slopes may \
c//!  only get diffuse radiation, the light must first be split into \
c//!  direct (beam) and diffuse light. No account is taken of the change in\
c//!  surface length, or of any shading beyond the current slope. Hence \
c//!  total radiation over many surfaces will be greater than over a single,\
c//!  flat plane. The difference is aprox 1/cos(slope).
c//!  Bibliography:Duffie and Beckman (1991). Solar engineering of thermal\
c//!    processes. Wiley. Also See groups.google.com; alt.energy.renewable;\
c//!    Subject: Re: Irradiation; Posted 5/1/1998; accessed 15/6/2001. \
c//!    Diffuse radiation from Montieth J.L. (1973) Principles of \
c//!    Environmantal Physics, Arnold, London.
c//! E: see the references

      RETURN
      END


      REAL*8 FUNCTION TEMP_ALT_CORR_FN(TEMP,
     1   METALT, ALTITUDE)
      IMPLICIT NONE

c**  Function corrects temperature for altitude - should be applied
c**  directly after the temperaure calculation as will therefore
c**  impact on other functions eg SVP etc
      REAL*8  TEMP, METALT, ALTITUDE, ALTDIFF_FT
c
       ALTDIFF_FT = (ALTITUDE-METALT)*100/2.54/12.0
       TEMP_ALT_CORR_FN = TEMP-(ALTDIFF_FT)/1000.0*2.0                                 
c//! Equation: 78
c//! Description: corrects temperature for site altitude. correction \
c//!   is approc 2C per 1000 ft \\
c//!   INPUTS: Temperature; altitude of met station (m); \
c//!      Altitude of site (m)
c//! Bibliography:  ??
c//!  
c//!  E: T_{alt} = T-ALT_{diff}*2
c
c//! Variable: T_{alt}
c//! Description: corrected air temperature \\
c//! Units: C
c//! Type: REAL*8
c//! SET
c
c//! Variable: T
c//! Description: Temperature at met station elevation
c//! Units: C
c//! Type: REAL*8
c
c//! Variable: ALT_{diff}
c//! Description: Altitude difference between site and met station
c//! Units: Feet
c//! Type: REAL*8
c
      RETURN
      END


      REAL*8 FUNCTION PPT_RAIN_FN(TEMP_DAY_MEAN, RAIN_DAY)
      IMPLICIT NONE
c

      REAL*8 ZERO_SNOW_T, ALL_SNOW_T, TEMP_DAY_MEAN, RAIN_DAY
      REAL*8 RAINPART      
	ZERO_SNOW_T=2.0
      ALL_SNOW_T=-2.0

	IF(TEMP_DAY_MEAN.GT.ZERO_SNOW_T) THEN
         PPT_RAIN_FN = RAIN_DAY*1.0
      ELSE IF(TEMP_DAY_MEAN.LT.ALL_SNOW_T) THEN
         PPT_RAIN_FN = 0.0
      ELSE
c       some combination of snow and rain
        RAINPART = (TEMP_DAY_MEAN - ALL_SNOW_T)/
     1         (ZERO_SNOW_T-ALL_SNOW_T)
        PPT_RAIN_FN = RAINPART*RAIN_DAY
      ENDIF
        
c
c//! Equation: 79
c//! Description: Finds the amount of real rain from Gross rain \
c//!  ie. that which is not snow \\ 
c//!  INPUTS: Mean Temperature of the day; Daily rainfall.
         
c//! Bibliography: S Evans; SWELTER. eq.67
c
c//! Case: T_{mean} > T_{snow}
c//!      E:  PPT_{rain} = PPT_{day}
c//! Case: T_{mean} < T_{all_snow}
c//!      E: PPT_{rain} = 0.0
c//! Case: T_{all_snow} < T_{mean} < T_{snow}
c//!      E: PPT_{rain} = (T_{mean} - T_{all_snow})/    \
c//!           (T_{snow}-T_{all_snow}) * PPT_{day}
c   
c//! Variable: T_{Mean}
c//! Description: Mean daily temperature
c//! Units: Celcius
c//! Type: REAL*8
c
c//! Variable: PPT_{day}
c//! Description: Gross precipitation
c//! Units: mm/day
c//! Type: REAL*8
c    
c//! Variable: T_{snow} 
c//! Description: Temperature at which snow may start (above which is all rain)        
c//! Units: Celcius
c//! Type: REAL*8
c
c//! Variable: T_{all_snow}
c//! Description: Temperature below which is all snow 
c//! Units: Celcius
c//! Type: REAL*8          
c
c//! Variable: PPT_{rain}
c//! Description: Amount of precipitation which is real rain
c//! Units: mm/day  (rain)          
c//! Type: REAL*8
c//! SET   
c
      RETURN                                           
      END                                              

      REAL*8 FUNCTION PPT_SNOW_FN(RAIN_DAY, PPT_RAIN)
      IMPLICIT NONE
c

      REAL*8 RAIN_DAY, PPT_RAIN   
   
	PPT_SNOW_FN = (RAIN_DAY - PPT_RAIN)*1.1811
        
c
c//! Equation: 80
c//! Description: Finds the amount of snow from Gross rain \\
c//!  INPUTS: Daily rainfall; Precipitation which is rain
         
c//! Bibliography: S Evans; SWELTER. eq.67
c
c//!      E: PPT_{snow} = (Rain{day} - PPT_{Rain})*1.1811
c   
c//! Variable: PPT_{day}
c//! Description: Gross precipitation
c//! Units: mm/day
c//! Type: REAL*8
c
c//! Variable: PPT_{Rain}
c//! Desciption: Precipitation which is Rain
c//! Units: mm/day  (rain)
c    
c//! Variable: PPT_{snow}
c//! Description: precipitation as snow (snow depth)
c//! Units: mm/day  (snow)          
c//! Type: REAL*8
c//! SET   
c
      RETURN                                           
      END     



      REAL*8 Function TRANSNORM_FN(MU, SD, Xval)
      Implicit none

      REAL*8 MU, SD, Xval
c     Transforms a standard normal distribution, N(0,1))
c     to N(Mu, SD)
c
C//! Equation: 81
c//! Description: Transforms a value from a Standard Normal distribution N(0,1) to\
c//!   a 'generic' normal N(mu, sd)\\
c//!   INPUTS: Mean of Distribution; Sd of distribution; Standard Normal value
c//! Bibliography: None (T Houston Pers Comm)
      TRANSNORM_FN = Xval*SD+Mu
c
c//! E: N(\mu, \sigma) = N(0,1)*\sigma + \mu
c
c//! variable: \mu
c//! Type: REAL*8
c//! Description: Mean of the desired normal distribution
c//! Units: Any
c
c//! Variable: \sigma
c//! Type: REAL*8
c//! Description: Standard deviation of the desired normal distribution
c//! Units: Any
c
c//! Variable: N(0,1)
c//! Type: Double Precsion
c//! Description: Value from the standard normal distribution
c//! Units: Any
c
c//! Variable: N(\mu, \sigma)
c//! Type: REAL*8
c//! Description: Transformed Normal Value - from the shifted norma distribution
c//! Units: Any
c//! SET
c
      Return
      End




      REAL*8 Function NORMSKEW_FN(NormVal, Lambda)
      Implicit none

c
      REAL*8 NormVal, Lambda
c     Box_cox Normaility transformation of a Normal Distribution
c
C//! Equation: 82
c//! Description: Skews a normal disribution\\
c//!   INPUTS: Value from a normal distribution N(Mu,Sd);\
c//!     And Response variable
c//! Bibliography: C. Coarkin & P Tobias, NIST/SEMATECH Engineering Handbook;\
c//!   Ch. 1.3.3.6\
c//!  http://www.itl.nist.gov/div898/handbook/eda/section3/boxcox.htm (9/8/2001)
c
c//! E: SN(N(\mu,\sigma)) = (N(\mu,\sigma)^\lambda -1)\over{\lambda}
c  
      NORMSKEW_FN = ( (Normval**lambda)-1)/lambda
c
c//! Variable: \\lambda
c//! Type: REAL*8
c//! Description: Response variable - if \lambda= 0; use log of data
c//! Units: none
c
c//! Variable: N(\mu,\sigma)
c//! Type: REAL*8
c//! Description: value from the normal distribution
c//! Units: Any
c//! Variable: SN(N(\mu,\sigma)) 
c//! Type: REAL*8
c//! Description: box-cox skewed, (transformed) normal value
c//! Units: Any
c//! SET
c
      Return
      End

      REAL*8 FUNCTION TAMPFROMRH_FN(RH)
      IMPLICIT NONE

c** humidity based on the temperature max-min amplitude
      REAL*8 RH

      TAMPFROMRH_FN = (RH-110.3)/(-4.68)
c
c//! Equation 83
c//! Description: Returns an estimate of the daily air temperature\
c//!  amplitude based on RH. Inverse of RH equation\\
c//!   INPUTS: Relative humidity (%)
c//! Bibliography: S Evans (SWELTER)
c
c//! E: \Delta T_{air} = (RH-110.3)\over{(-4.68)}
c
c//! Variable: RH
c//! Description: Relative humidity
c//! Units: %
c//! Type: REAL*8
c
c//! Variable: \Delta T_{air}
c//! Description: Amplitude of temperature variation during the day
c//! Units: \deg C
c//! Type: REAL*8
c//! SET
c
      RETURN
      END
                               
