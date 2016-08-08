      real*8 ans(72,72,12),x(12),latf,lon0,del,lat,lon
      integer i,j,k,ians(72,72,12),ix(12)

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


      open(11,file='stdhum.dat')
      read(11,*)
      do i=1,1232
        read(11,*) lon,lat,x
        do j=1,12
          ans(int((latf-lat)/del+0.5d0),int((lon-lon0)/del+
     &0.5d0),j) = x(j)
        enddo
      enddo
      close(11)

      open(11,file='../stdhum.dat')
      do i=1,72
        do j=1,72
          write(11,'(12f6.2)') (ans(i,j,k),k=1,12)
        enddo
      enddo
      close(11)

      open(11,file='stdtmp.dat')
      read(11,*)
      do i=1,1232
        read(11,*) lon,lat,x
        do j=1,12
          ans(int((latf-lat)/del+0.5d0),int((lon-lon0)/del+
     &0.5d0),j) = x(j)
        enddo
      enddo
      close(11)

      open(11,file='../stdtmp.dat')
      do i=1,72
        do j=1,72
          write(11,'(12f6.2)') (ans(i,j,k),k=1,12)
        enddo
      enddo
      close(11)

      open(11,file='raindays.dat')
      read(11,*)
      do i=1,1232
        read(11,*) lon,lat,ix
        do j=1,12
          ans(int((latf-lat)/del+0.5d0),int((lon-lon0)/del+
     &0.5d0),j) = real(ix(j))
        enddo
      enddo
      close(11)

      open(11,file='../raindays.dat')
      do i=1,72
        do j=1,72
          write(11,'(12f6.2)') (ans(i,j,k),k=1,12)
        enddo
      enddo
      close(11)


      stop
      end
