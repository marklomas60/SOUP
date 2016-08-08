      INTEGER nyrs,yr0,no_files
      REAL*8 latres,lonres,lat0,latf,lon0,lonf
      CHARACTER stinput*60

      PARAMETER(nyrs=100,yr0=1901)
      PARAMETER(stinput  = '../../kruger/xxx/')
      PARAMETER(latres=5.0d0,lonres=5.0d0)
      PARAMETER(lat0=-99.0d0,latf =99.0d0,lon0=-400.0d0,lonf=400.0d0)
      PARAMETER(no_files=18)

      real*8 lat,lon,npp,pi,rad,circ,x,y,area,xx,scale,sum,out2
      real*8 sum1(500),sum2(500),sum3(500),nppt(500),ans,out1
      integer i,j,k,blank,file_no,av
      character  stvar(no_files)*60,st*60

      pi = 3.1415926D0
      circ = 40008.0D0
      rad = circ/(2.0D0*pi)
      y = circ*latres/360.0D0

      DATA stvar/'gpp.dat','npp.dat','nbp.dat','nbpsum.dat',
     &'lai.dat',
     &'evt.dat','trn.dat','rof.dat'
     &'biot.dat','stembio.dat','rootbio.dat','nppstore.dat',
     &'scn.dat','snn.dat','sresp.dat',
     &'fcn.dat',
     &'prc.dat','tmp.dat'/

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
      open(21,file=st(1:blank(st))//'sums/'//'all_sums.dat')

      DO i=1,nyrs
        sum1(i) = 0.0D0
        sum2(i) = 0.0D0
        sum3(i) = 0.0D0
      ENDDO

      j = 0

10    continue

 
        read(10,*,end = 90) lat,lon,(nppt(i),i=1,nyrs)

        j = j + 1

          if ((lat.GT.lat0).and.(lat.LT.latf).and.
     &(lon.GT.lon0).and.(lon.LT.lonf)) then

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

      WRITE(*,'(''Total npp       = '',F9.3,'' giga tonnes'')') sum1(nyrs)
      WRITE(*,'(''Vegetation area = '',F9.4,'' mega Km-2'')') sum2(nyrs)
      WRITE(*,'(''Land area       = '',F9.4,'' mega Km-2'')') sum3(nyrs)

      open(11,file=st(1:blank(st))//'sums/'//stvar(file_no))

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



