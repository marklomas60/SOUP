      real*8 latf,lon0,latdel,londel,ans,lat,lon
      integer latn,lonn,count,ians

      open(11,file='../readme.dat')
      read(11,*)
      read(11,*)
      read(11,*) latf,lon0
      read(11,*)
      read(11,*) latdel,londel
      read(11,*)
      read(11,*) latn,lonn
      close(11)

      open(11,file='../cont_lu-1-1700.dat',status='old')
      open(21,file='../lat_lon.dat')
      open(22,file='../lat_lon2.dat')

      count = 0
10    read(11,*,end=20) ians
        lat = latf - real(count/lonn)*latdel
        lon = lon0 + mod(count,lonn)*londel
        if (ians.lt.255) write(21,'(f7.3,f9.3)') lat,lon
        write(22,'(f7.3,f9.3)') lat,lon
        count = count + 1
      goto 10
20    continue

      close(11)
      close(21)
      close(22)

      stop
      end
