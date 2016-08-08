* namal must be a factor of grid *
      integer namal
      parameter (namal = 2)
      integer i,j,k,l,ii,jj
      integer count,fileno
      real ans(360,720),lat,lon,ans1
      character st1*10

      do fileno=1,23

      write(st1,'(i0)') fileno
      print*,st1

      if (fileno.LT.10) then
      open(11,file='30min/cont_lu-'//st1(1:1)//'-2000.dat',
     &status='old')
      open(21,file='1deg/cont_lu-'//st1(1:1)//'-2000.dat')
      else
      open(11,file='30min/cont_lu-'//st1(1:2)//'-2000.dat',
     &status='old')
      open(21,file='1deg/cont_lu-'//st1(1:2)//'-2000.dat')
      endif

      do i=1,292
        do j=1,720
          read(11,*) ans(i,j)
        enddo
      enddo
      close(11)

      do i=293,360
        do j=1,720
          ans(i,j) = 255
        enddo
      enddo

      do i=1,360/namal
        do j=1,720/namal
          count = 0
          ans1 = 0.0
          do k=1,namal
            do l=1,namal
              ii = (i-1)*namal+k
              jj = (j-1)*namal+l
              if ((ans(ii,jj).GE.0).AND.(ans(ii,jj).LT.255)) then
                count = count + 1
                ans1 = ans1 + ans(ii,jj)
              endif
            enddo
          enddo
          if (count.GT.0) then
            write(21,'(i3)') int(real(ans1)/real(count)+0.5)
          else
            write(21,'(i3)') 255
          endif
        enddo
      enddo
      close(21)

      enddo

      stop
      end
