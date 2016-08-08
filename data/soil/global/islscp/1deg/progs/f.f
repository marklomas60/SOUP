      REAL*8 latf,lon0,lat,lon,dep,ts,tsi,tc,depth
      integer tex,i,j

      open(11,file='soldep.dat')
      open(12,file='soltex.dat')
      
      open(21,file='data.dat')
      open(22,file='lat_lon.dat')

      latf = 89.5d0
      lon0 = -179.5d0

      do i=1,180
        do j=1,360
          read(11,*) lat,lon,depth
          read(12,*) lat,lon,tex
          CALL CLSASI(tex,tc,ts,tsi)
          lat = latf - real(i-1)
          lon = lon0 + real(j-1)
          write(21,'(f7.3,f9.3,4f12.3)') lat,lon,ts,tsi,-999.0d0,depth
          if (ts+tsi.gt.0.0d0) write(22,'(f7.3,f9.3,f8.4)') lat,lon,ts
        enddo
      enddo



      stop
      end


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
      REAL*8 tc,ts,tsi
      INTEGER tex

      IF ((tex.EQ.0).OR.(tex.EQ.6)) THEN
        tc = -999.0d0
        ts = -999.0d0
        tsi =-999.0d0
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

