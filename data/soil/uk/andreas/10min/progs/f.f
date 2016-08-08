      real*8 ans(72,72,12),x(10),latf,lon0,del,lat,lon
      integer i,j,k,ians(72,72,12)

      latf = 60
      lon0 =-10
      del = 1.0d0/6.0d0

      do i=1,72
        do j=1,72
          do k=1,12
            ans(i,j,k) = -1.0d0
            ians(i,j,k) = -1
          enddo
        enddo
      enddo

      open(11,file='data.txt')
      do i=1,1065
        read(11,*) lat,lon,x
        do j=1,10
          if (x(j).gt.998.0d0) x(j) = -999.0d0
          ans(int((latf-lat)/del+1.5d0),int((lon-lon0)/del+
     &1.5d0),j) = x(j)
        enddo
      enddo
      close(11)

      open(11,file='../data.dat')
      do i=1,72
        do j=1,72
          lat = latf - real(i-1)*del
          lon = lon0 + real(j-1)*del
          write(11,'(f7.3,f9.3,4f12.3)') lat,lon,ans(i,j,8),ans(i,j,9),
     &ans(i,j,10),0.0d0
        enddo
      enddo
      close(11)


      stop
      end
