      INTEGER max_c,max_y
      PARAMETER (max_c = 50000, max_y = 10000)
      REAL*8 latres,lonres,pi,circ,rad,y,ansv(max_y),area,c_sum(max_c)
      REAL*8 c_area(max_c),ans,lat,lon
      INTEGER years,yr0,yrf,kode,i,count,index
      CHARACTER st1*16,country(max_c)*16

      OPEN(11,file='/data/SDGVM/output/global3_1deg/site_info.dat',
     &STATUS ='old',iostat=kode)
      IF (kode.NE.0) THEN
        WRITE(*,*) 'Directory does not contain ''site_info.dat'' file.'
      ENDIF
      OPEN(12,file='/data/SDGVM/output/global3_1deg/prc.dat',
     &STATUS='old',iostat=kode)
      IF (kode.NE.0) THEN
        WRITE(*,*) 'Directory does not contain input file.'
      ENDIF
      OPEN(21,file='/data/SDGVM/output/global3_1deg/country_sums/sumsprc
     &.dat',iostat=kode)
     
      latres = 1.0d0
      lonres = 1.0d0
      years = 100
      yr0 = 1
      yrf = 79

      pi = 3.1415926d0
      circ = 40008.0d0
      rad = circ/(2.0d0*pi)
      y = circ*latres/360.0d0
      
      DO i=1,max_c
        c_sum(i) = 0.0d0
	c_area(i) = 0.0d0
      ENDDO
      
      count = 0
10    CONTINUE
        count = count + 1

        READ(11,'(A16,i5)',end=20) st1,index
        READ(12,*) lat,lon,(ansv(i),i=1,years)

        IF (index.gt.0) THEN

          DO i=1,16
	    country(index)(i:i) = st1(i:i)
          ENDDO

          ans = 0.0d0
	  DO i=yr0,yrf
	    ans = ans + ansv(i)
	  ENDDO
          ans = ans/real(yrf-yr0+1)

          x = 2.0d0*pi*rad*cos(lat*pi/180.0d0)*lonres/360.0d0
          area = x*y

	  c_sum(index) = c_sum(index) + area*ans
	  c_area(index) = c_area(index) + area

        ENDIF

        IF (count.le.-10000) GOTO 20
      GOTO 10
20    CONTINUE

      CLOSE(11)
      CLOSE(12)

      DO i=1,max_c
        IF (c_area(i).GT.0.0d0) THEN
          c_sum(i) = c_sum(i)/c_area(i)
          c_area(i) = c_area(i)*1.0e-3
	  WRITE(21,'(A16,i6,F8.1,f8.1)') country(i),i,c_sum(i),c_area(i)
	ENDIF
      ENDDO


      STOP
      END
