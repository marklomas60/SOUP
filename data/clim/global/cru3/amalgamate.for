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

      open(22,file='1deg/lat_lon.dat')
      open(23,file='1deg/maskmap.dat')
      open(24,file='1deg/maskmap2.dat')

      print*,'tmp'
      open(11,file='30min/tmp.dat',status='old',
     &access='direct',recl=recl1,form='formatted')
      open(21,file='1deg/tmp.dat')

      do i=1,360
        do j=1,720
          map2(j,i) = 0
        enddo
      enddo

      count = 0
      do i=1,360/namal
*      do i=1,3
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
            count = count + 1
            do year=1,nyears
              write(21,'(12i4)') 
     &(int(real(ans2(year,mnth))/real(nsites)+.5),
     &mnth=1,12)
            enddo
            write(22,'(f7.3,f9.3)') 90.0-(real(i-1)+.5)*real(namal)*.5,
     &(real(j-1)+.5)*real(namal)*.5-180.0
            map2(j,i) = count
          endif
        enddo
      enddo
      close(11)
      close(21)
      close(22)

      do i=1,360/namal
        write(23,'(1000i6)') (map2(j,i),j=1,720/namal)
        write(24,'(1000i1)') (min(map2(j,i),1),j=1,720/namal)
      enddo
      close(23)
      close(24)

      print*,'prc'
      open(11,file='30min/prc.dat',status='old',
     &access='direct',recl=recl1,form='formatted')
      open(21,file='1deg/prc.dat')

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

      print*,'cld'
      open(11,file='30min/cld.dat',status='old',
     &access='direct',recl=recl1,form='formatted')
      open(21,file='1deg/cld.dat')

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

      print*,'dtr'
      open(11,file='30min/dtr.dat',status='old',
     &access='direct',recl=recl1,form='formatted')
      open(21,file='1deg/dtr.dat')

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


      print*,'dtr_stn'
      open(11,file='30min/dtr_stn.dat',status='old',
     &access='direct',recl=recl1,form='formatted')
      open(21,file='1deg/dtr_stn.dat')

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


      print*,'dtr_cstn'
      open(11,file='30min/dtr_cstn.dat',status='old',
     &access='direct',recl=recl1,form='formatted')
      open(21,file='1deg/dtr_cstn.dat')

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


      print*,'frs'
      open(11,file='30min/frs.dat',status='old',
     &access='direct',recl=recl1,form='formatted')
      open(21,file='1deg/frs.dat')

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


      print*,'prc_stn'
      open(11,file='30min/prc_stn.dat',status='old',
     &access='direct',recl=recl1,form='formatted')
      open(21,file='1deg/prc_stn.dat')

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


      print*,'prc_cstn'
      open(11,file='30min/prc_cstn.dat',status='old',
     &access='direct',recl=recl1,form='formatted')
      open(21,file='1deg/prc_cstn.dat')

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


      print*,'tmn'
      open(11,file='30min/tmn.dat',status='old',
     &access='direct',recl=recl1,form='formatted')
      open(21,file='1deg/tmn.dat')

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


      print*,'tmx'
      open(11,file='30min/tmx.dat',status='old',
     &access='direct',recl=recl1,form='formatted')
      open(21,file='1deg/tmx.dat')

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


      print*,'tmp_stn'
      open(11,file='30min/tmp_stn.dat',status='old',
     &access='direct',recl=recl1,form='formatted')
      open(21,file='1deg/tmp_stn.dat')

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


      print*,'tmp_cstn'
      open(11,file='30min/tmp_cstn.dat',status='old',
     &access='direct',recl=recl1,form='formatted')
      open(21,file='1deg/tmp_cstn.dat')

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


      print*,'wet'
      open(11,file='30min/wet.dat',status='old',
     &access='direct',recl=recl1,form='formatted')
      open(21,file='1deg/wet.dat')

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


      print*,'vap'
      open(11,file='30min/vap.dat',status='old',
     &access='direct',recl=recl1,form='formatted')
      open(21,file='1deg/vap.dat')

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
