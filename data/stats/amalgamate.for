* namal must be a factor of grid *
      integer namal
      parameter (namal = 10)
      integer i,j,k,l,ii,jj
      integer count
      real ans(360,720,12),lat,lon,ans1(12)

      open(11,file='global/era40/30min/raindays.dat',status='old')
      open(21,file='global/era40/5deg/raindays.dat')

      do i=1,360
        do j=1,720
          read(11,'(12f6.2)') (ans(i,j,k),k=1,12)
        enddo
      enddo
      close(11)

      do i=1,360/namal
        do j=1,720/namal
          count = 0
          do m=1,12
            ans1(m) = 0.0
          enddo
          do k=1,namal
            do l=1,namal
              ii = (i-1)*namal+k
              jj = (j-1)*namal+l
              if (ans(ii,jj,1).GT.-1.0e-6) then
                count = count + 1
                do m=1,12
                  ans1(m) = ans1(m) + ans(ii,jj,m)
                enddo
              endif
            enddo
          enddo
          if (count.GT.0) then
            write(21,'(12f6.2)') ans1/real(count)
          else
            write(21,'(12f6.2)') -1.0,-1.0,
     &-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0
          endif
        enddo
      enddo
      close(21)



      open(11,file='global/era40/30min/stdtmp.dat',status='old')
      open(21,file='global/era40/5deg/stdtmp.dat')

      do i=1,360
        do j=1,720
          read(11,'(12f6.2)') (ans(i,j,k),k=1,12)
        enddo
      enddo
      close(11)

      do i=1,360/namal
        do j=1,720/namal
          count = 0
          do m=1,12
            ans1(m) = 0.0
          enddo
          do k=1,namal
            do l=1,namal
              ii = (i-1)*namal+k
              jj = (j-1)*namal+l
              if (ans(ii,jj,1).GT.-1.0e-6) then
                count = count + 1
                do m=1,12
                  ans1(m) = ans1(m) + ans(ii,jj,m)
                enddo
              endif
            enddo
          enddo
          if (count.GT.0) then
            write(21,'(12f6.2)') ans1/real(count)
          else
            write(21,'(12f6.2)') -1.0,-1.0,
     &-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0
          endif
        enddo
      enddo
      close(21)



      open(11,file='global/era40/30min/stdhum.dat',status='old')
      open(21,file='global/era40/5deg/stdhum.dat')

      do i=1,360
        do j=1,720
          read(11,'(12f6.2)') (ans(i,j,k),k=1,12)
        enddo
      enddo
      close(11)

      do i=1,360/namal
        do j=1,720/namal
          count = 0
          do m=1,12
            ans1(m) = 0.0
          enddo
          do k=1,namal
            do l=1,namal
              ii = (i-1)*namal+k
              jj = (j-1)*namal+l
              if (ans(ii,jj,1).GT.-1.0e-6) then
                count = count + 1
                do m=1,12
                  ans1(m) = ans1(m) + ans(ii,jj,m)
                enddo
              endif
            enddo
          enddo
          if (count.GT.0) then
            write(21,'(12f6.2)') ans1/real(count)
          else
            write(21,'(12f6.2)') -1.0,-1.0,
     &-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0
          endif
        enddo
      enddo
      close(21)


      stop
      end
