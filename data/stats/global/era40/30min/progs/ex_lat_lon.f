      real*8 latf,lon0,latdel,londel,ans,lat,lon
      integer latn,lonn,count

      open(11,file='readme.dat')
      read(11,*)
      read(11,*)
      read(11,*) latf,lon0
      read(11,*)
      read(11,*) latdel,londel
      read(11,*)
      read(11,*) latn,lonn
      close(11)

      open(11,file='raindays.dat')
      open(21,file='lat_lon.dat')
      open(22,file='lat_lon2.dat')

      count = 0
10    read(11,'(f6.2)',end=20) ans
        lat = latf - real(count/lonn)*latdel
        lon = lon0 + mod(count,lonn)*londel
        if (ans.gt.0.00001d0) write(21,'(f7.3,f9.3)') lat,lon
        write(22,'(f7.3,f9.3)') lat,lon
        count = count + 1
      goto 10
20    continue

      close(11)
      close(21)

      stop
      end
