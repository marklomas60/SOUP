      INTEGER max_c,max_y
      PARAMETER (max_c = 50000, max_y = 100)
      REAL*8 latres,lonres,pi,circ,rad,y,ansv(max_y),area
      REAL*8 c_area(max_c),ans,lat,lon,c_sum(max_c,max_y)
      INTEGER years,yr0,yrf,kode,i,j,count,index
      CHARACTER st1*16,country(max_c)*16

      OPEN(11,file='/home/mrl/networks2/site_info.dat',
     &STATUS ='old',iostat=kode)
      IF (kode.NE.0) THEN
        WRITE(*,*) 'Directory does not contain ''site_info.dat'' file.'
      ENDIF
      OPEN(12,file='/home/mrl/networks2/prc.dat',
     &STATUS='old',iostat=kode)
      IF (kode.NE.0) THEN
        WRITE(*,*) 'Directory does not contain input file.'
      ENDIF
      OPEN(21,file='/home/mrl/networks2/country_prc.dat',iostat=kode)
     
      latres = 0.5d0
      lonres = 0.5d0
      years = 100
      yr0 = 1
      yrf = 100

      pi = 3.1415926d0
      circ = 40008.0d0
      rad = circ/(2.0d0*pi)
      y = circ*latres/360.0d0
      
      DO i=1,max_c
	DO j=yr0,yrf
          c_sum(i,j) = 0.0d0
        ENDDO
	c_area(i) = 0.0d0
      ENDDO

      READ(11,*)
      
      count = 0
10    CONTINUE
        count = count + 1

        READ(11,'(A16,i5)',end=20) st1,index
        READ(12,*) lat,lon,(ansv(i),i=1,years)

        IF (index.gt.0) THEN

          DO i=1,16
	    country(index)(i:i) = st1(i:i)
          ENDDO

          x = 2.0d0*pi*rad*cos(lat*pi/180.0d0)*lonres/360.0d0
          area = x*y

	  DO i=yr0,yrf
	  c_sum(index,i) = c_sum(index,i) + area*ansv(i)
          ENDDO
	  c_area(index) = c_area(index) + area

        ENDIF

*        if (index.EQ.15000) print*,'sp ',ans,lat,lon
*        if (index.EQ.17000) print*,'uk ',ans,lat,lon

        IF (count.le.-10000) GOTO 20
      GOTO 10
20    CONTINUE

      CLOSE(11)
      CLOSE(12)

      DO i=1,max_c
        IF (c_area(i).GT.0.0d0) THEN
	  WRITE(21,'(A16,i6,100f8.1)') country(i),i,
     &(c_sum(i,j)/c_area(i),j=1,yrf-yr0+1)
	ENDIF
      ENDDO


      STOP
      END
