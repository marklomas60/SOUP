      IMPLICIT NONE
      INTEGER nyrs,yr0,no_files
      REAL*8 latres,lonres
      CHARACTER stinput*60,stoutput*60

      PARAMETER(nyrs=100,yr0=1901)
      PARAMETER(stinput  = '/data/SDGVM/output/newnew/')
      PARAMETER(stoutput = '/data/SDGVM/output/newnew/sums_41/')
      PARAMETER(latres=0.5d0,lonres=0.5d0)
      PARAMETER(no_files=18)


      real*8 lat,lon,npp,pi,rad,circ,x,y,area,xx,scale,sum,out2
      real*8 sum1(500),sum2(500),sum3(500),nppt(500),ans,out1
      integer i,j,k,blank,file_no,av,n_poly
      character  stvar(no_files)*60,st*60
      logical pll
      PARAMETER(n_poly=5)
      real*8 plats(n_poly),plons(n_poly)
* 41
      DATA plats/   2.5,   1.0, -15.0, -15.0,   2.5/
      DATA plons/ -59.0, -50.0, -51.5, -59.0, -59.0/

* 42
*      DATA plats/ -13.40, -2.67, -15.00, -20.50, -13.00/
*      DATA plons/ -71.12, -59.0, -59.50, -63.50, -71.12/

* 43
*      DATA plats/   4.00, -4.00, -15.50,  -6.00,   4.00/
*      DATA plons/ -76.50, -60.5, -73.50, -80.00, -76.50/

* 44
*      DATA plats/   2.50,  5.00,  -2.67,  -4.00,   2.50/
*      DATA plons/ -73.40, -59.0, -59.00, -60.50, -73.40/

* 45
*      DATA plats/  -8.00, -1.80,  -4.00, -13.32,  -8.00/
*      DATA plons/ -73.30, -65.0, -60.50, -71.02, -73.30/

      DATA stvar/'gpp.dat','npp.dat','nbp.dat','nbpsum.dat',
     &'lai.dat','evt.dat','trn.dat','rof.dat',
     &'biot.dat','stembio.dat','rootbio.dat','nppstore.dat',
     &'scn.dat','snn.dat','sresp.dat','fcn.dat','prc.dat','tmp.dat'/

      pi = 3.1415926D0
      circ = 40008.0D0
      rad = circ/(2.0D0*pi)
      y = circ*latres/360.0D0

      open(21,file=stoutput(1:blank(stoutput))//'all_sums.dat')
      open(22,file=stoutput(1:blank(stoutput))//'lat_lon.dat')

      do file_no=1,no_files

      scale = 1.0d0

      if (stvar(file_no).eq.'nbpsum.dat') then
        sum = 1.0d0
        k = 3
      else
        sum = 0.0d0
        k = file_no
      endif

      if ((stvar(file_no).eq.'lai.dat').or.
     &(stvar(file_no).eq.'tmp.dat')) then
        av = 1
      else
        av = 0
      endif

      st = stinput
      open(10,file=st(1:blank(st))//stvar(k))

      DO i=1,nyrs
        sum1(i) = 0.0D0
        sum2(i) = 0.0D0
        sum3(i) = 0.0D0
      ENDDO

      j = 0

10    continue

 
        read(10,*,end = 90) lat,lon,(nppt(i),i=1,nyrs)
        CALL POLYGON(lat,lon,plats,plons,n_poly,pll)
        if (file_no.eq.1) then
          if (pll) then
            write(22,'(f7.3,f9.3)') lat,lon
          endif
        endif

        j = j + 1
*        if (mod(j,1000).eq.0.0d0)  print*,j

          if (pll) then

          do i=1,nyrs

            npp = nppt(i)

* Put npp in units of tonnes per kilometer squared.
            npp = npp*scale
            x = 2.0D0*pi*rad*cos(lat*pi/180.0D0)*lonres/360.0D0
            area = x*y
            sum1(i) = sum1(i)+ npp*area
            sum2(i) = sum2(i) + area
            sum3(i) = sum3(i) + area
          enddo

        endif

      goto 10
90    continue

      close(unit=10)

* Areas are in mega KM and npp is in giga tonnes.
      DO i=1,nyrs
        sum1(i) = sum1(i)/1.0e9
        sum2(i) = sum2(i)/1.0e6
        sum3(i) = sum3(i)/1.0e6
      ENDDO

      print'(''Total npp       = '',F9.3,'' giga tonnes'')', 
     &sum1(nyrs)
      print'(''Vegetation area = '',F9.4,'' mega Km^2'')', sum2(nyrs)
      print'(''Land area       = '',F9.4,'' mega Km^2'')', sum3(nyrs)

      open(11,file=stoutput(1:blank(stoutput))//stvar(file_no))

      av = 1
      ans = 0.0d0
      do i=1,nyrs
        xx = 1.0d0
        if (av.EQ.1)  xx = sum2(i)/1000.0d0
        write(11,1000) yr0-1+i,sum1(i)/xx + ans*sum
        if (i.EQ.1) out1 = sum1(i)/xx + ans*sum
        if (i.EQ.nyrs) out2 = sum1(i)/xx + ans*sum
        ans = ans + sum1(i)
1000    format(i6,f12.5)
      enddo

      if (file_no.eq.1) write(21,'(''Area'',f8.3,'' mega Km^2'')') 
     &sum3(nyrs)

      write(21,'(A15,2f12.5)') stvar(file_no),out1,out2


      if (av.eq.1)  write(*,'(''Average         = '',f9.3)') 
     &sum1(nyrs)/xx

      close(unit=11)


      enddo



      stop 
      end



*----------------------------------------------------------------------*
*                                                                      *
*                          FUNCTION blank                              *
*                          **************                              *
*                                                                      *
* 'blank' returns the number of characters before the space in the     *
* string which is its argument.                                        *
*                                                                      *
*----------------------------------------------------------------------*
      FUNCTION blank(st1)
*----------------------------------------------------------------------*
      CHARACTER st1*60
      INTEGER blank

      blank = 0
10    CONTINUE
      blank = blank + 1
      IF ((ichar(st1(blank:blank)).NE.32).AND.(blank.LE.60))  GOTO 10
      blank = blank - 1


      RETURN
      END


      SUBROUTINE POLYGON(lat,lon,plats,plons,n_poly,pll)
      integer n_poly
      real*8 lat,lon,pi,x1,x2,y1,y2,angle,dtheta,theta1,theta2
      real*8 plats(n_poly),plons(n_poly)
      logical pll

      pi = 3.1415926

      angle = 0.0

      do i=1,n_poly-1

        x1 = plats(i) - lat
        y1 = plons(i) - lon
        x2 = plats(i+1) - lat
        y2 = plons(i+1) - lon

        theta1 = atan2(y1,x1)
        theta2 = atan2(y2,x2)
        dtheta = theta2 - theta1
        IF (dtheta.GT. pi) dtheta = dtheta - 2.0*pi
        IF (dtheta.LT.-pi) dtheta = dtheta + 2.0*pi
        angle = angle + dtheta
*        print*,angle,dtheta

      enddo

*      print*,angle*180/pi,angle,angle

      IF (abs(angle).LT.pi) THEN
        pll = .false.
      ELSE
        pll = .true.
      ENDIF
*      print*,pll

      RETURN
      END










