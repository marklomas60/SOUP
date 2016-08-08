* namal must be a factor of grid *
      integer namal,nyears
      parameter (namal = 2, nyears = 106)
      integer map(720,360),i,j,k,l,recl1,ind,mnth,recn
      integer ans1(12),ans2(200,12),year,nsites,count
      integer map2(720,360)

      open(11,file='30min/maskmap.dat',status='old')
      do j=1,360
        read(11,*) (map(i,j),i=1,720)
      enddo
      close(11)

      recl1 = 12*4

      print*,'hum'
      open(11,file='30min/hum.dat',status='old',
     &access='direct',recl=recl1,form='formatted')
      open(21,file='1deg/hum.dat')

      do i=1,360/namal
        do j=1,720/namal
          do year=1,nyears
            do mnth=1,12
              ans2(year,mnth) = 0
            enddo
          enddo
          nsites = 0
          do k=1,namal
            do l=1,namal
              ind = map((j-1)*namal+l,(i-1)*namal+k)
              if (ind.gt.0) then
                nsites = nsites + 1
                do year=1,nyears
                  recn = (ind-1)*nyears + year
                  read(11,'(12i4)',rec=recn) ans1
                  do mnth=1,12
                    ans2(year,mnth) = ans2(year,mnth) + ans1(mnth)
                  enddo
                enddo
              endif
            enddo
          enddo
          if (nsites.gt.0) then
            do year=1,nyears
              write(21,'(12i4)') 
     &(int(real(ans2(year,mnth))/real(nsites)+.5),
     &mnth=1,12)
            enddo
          endif
        enddo
      enddo
      close(11)
      close(21)


      stop
      end
