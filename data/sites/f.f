      real lat,lon
      integer i

      open(11,file='gpp.dat')
      open(21,file='uk_10min.dat')

      do i=1,1065
        read(11,*) lat,lon
        write(21,'(f7.3,f9.3)') lat,lon
      enddo

      stop
      end
