      real*8 ans(360,720,12),x(12),latf,lon0,del,lat,lon
      real*8 lat_lon(2,70000)
      integer i,j,k,ians(360,720,12),ix(12)

      latf =  90.0
      lon0 =-180.0d0
      del = 0.5d0

      open(11,file='lat_lon.dat')
        do i=1,67420
          read(11,*) lat_lon(1,i),lat_lon(2,i)
        enddo
      close(11)


      do i=1,360
        do j=1,720
          do k=1,12
            ans(i,j,k) = -1.0d0
            ians(i,j,k) = -1
          enddo
        enddo
      enddo

      open(11,file='stdhum.dat')
      do i=1,67420
        lat = lat_lon(1,i)
        lon = lat_lon(2,i)
        read(11,'(12i4)') ix
        do j=1,12
          ans(int((latf-lat)/del+1.0d0),int((lon-lon0)/del+
     &1.0d0),j) = real(ix(j))/100.0d0
        enddo
      enddo
      close(11)

      open(11,file='../stdhum.dat')
      do i=1,360
        do j=1,720
          write(11,'(12f6.2)') (ans(i,j,k),k=1,12)
        enddo
      enddo
      close(11)


      do i=1,360
        do j=1,720
          do k=1,12
            ans(i,j,k) = -1.0d0
            ians(i,j,k) = -1
          enddo
        enddo
      enddo

      open(11,file='stdtmp.dat')
      do i=1,67420
        lat = lat_lon(1,i)
        lon = lat_lon(2,i)
        read(11,'(12i4)') ix
        do j=1,12
          ans(int((latf-lat)/del+1.0d0),int((lon-lon0)/del+
     &1.0d0),j) = real(ix(j))/100.0d0
        enddo
      enddo
      close(11)

      open(11,file='../stdtmp.dat')
      do i=1,360
        do j=1,720
          write(11,'(12f6.2)') (ans(i,j,k),k=1,12)
        enddo
      enddo
      close(11)


      do i=1,360
        do j=1,720
          do k=1,12
            ans(i,j,k) = -1.0d0
            ians(i,j,k) = -1
          enddo
        enddo
      enddo

      open(11,file='raindays.dat')
      do i=1,67420
        lat = lat_lon(1,i)
        lon = lat_lon(2,i)
        read(11,'(12i4)') ix
        do j=1,12
          ans(int((latf-lat)/del+1.0d0),int((lon-lon0)/del+
     &1.0d0),j) = real(ix(j))/100.0d0
        enddo
      enddo
      close(11)

      open(11,file='../raindays.dat')
      do i=1,360
        do j=1,720
          write(11,'(12f6.2)') (ans(i,j,k),k=1,12)
        enddo
      enddo
      close(11)


      stop
      end
