* namal must be a factor of grid *
      integer namal
      parameter (namal = 5)
      integer i,j,k,l,ii,jj
      integer count
      real ans(180,360,7),lat,lon,ans1(7)

      open(11,file='global/islscp2/1deg/0-150/data.dat',status='old')
      open(21,file='global/islscp2/5deg/0-150/data.dat')

      do i=1,180
        print*,i
        do j=1,360
          read(11,*) lat,lon,(ans(i,j,k),k=1,7)
        enddo
      enddo
      close(11)

      print*,'hello'

      do i=1,180/namal
        do j=1,360/namal
          count = 0
          do m=1,7
            ans1(m) = 0.0
          enddo
          do k=1,namal
            do l=1,namal
              ii = (i-1)*namal+k
              jj = (j-1)*namal+l
              if (ans(ii,jj,1).GT.0.0) then
                count = count + 1
                do m=1,7
                  ans1(m) = ans1(m) + ans(ii,jj,m)
                enddo
              endif
            enddo
          enddo
          lat = 90.0 - real(i-1)*namal - real(namal)*0.5
          lon = -180.0 + real(j-1)*namal + real(namal)*0.5
          if (count.GT.0) then
            write(21,'(f7.3,f9.3,8f10.4)') lat,lon,(ans1(m)/real(count),
     &m=1,7),100.0
          else
            write(21,'(f7.3,f9.3,8f10.4)') lat,lon,-99.0,-99.0,
     &-99.0,-99.0,-99.0,-99.0,-99.0,-999.0
          endif
        enddo
      enddo
      close(21)

      stop
      end
