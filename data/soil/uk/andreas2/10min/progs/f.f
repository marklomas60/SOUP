      integer i,j,k
      real lat,lon,c1,c30,org1,org30,cl30,sa30,si30,bulk30
      real cl1,sa1,si1,bulk1,map(6,60,60),lat0,lon0

      lat0 = 50.0833
      lon0 = -7.0833

      do i=1,60
        do j=1,60
          do k=1,6
            map(k,i,j) = -999
          enddo
          map(3,i,j)= 0.0
          map(6,i,j)= 0.0
        enddo
      enddo

      open(11,file='UKcarbtext_sdgvm10min.txt',status='old')
      open(21,file='../30cm/data.txt')
      open(22,file='../1m/data.txt')
      read(11,*)
      do i=1,1061
        read(11,*) lat,lon,c1,c30,org1,org30,cl30,sa30,si30,
     &bulk30,cl1,sa1,si1,bulk1
        if (sa30.LT.200) then
        ind1 = int((lat - lat0)*6.0 + .5)+1
        ind2 = int((lon - lon0)*6.0 + .5)+1
        map(1,ind1,ind2) = sa30/(sa30+si30+cl30)*100.0
        map(2,ind1,ind2) = si30/(sa30+si30+cl30)*100.0
        map(3,ind1,ind2) = bulk30
        map(4,ind1,ind2) = sa1/(sa1+si1+cl1)*100.0
        map(5,ind1,ind2) = si1/(sa1+si1+cl1)*100.0
        map(6,ind1,ind2) = bulk1
        if ((ind1.gt.60).or.(ind1.lt.1).or.(ind2.gt.60).or.
     &(ind2.lt.1)) stop
        endif
      enddo

      do i=1,60
        do j=1,60
          write(21,'(f7.3,f9.3,4f12.3)') real(60-i)/6.0+lat0,
     &real(j-1)/6.0+lon0,(map(k,61-i,j),k=1,3),0.0
          write(22,'(f7.3,f9.3,4f12.3)') real(60-i)/6.0+lat0,
     &real(j-1)/6.0+lon0,(map(k,61-i,j),k=4,6),0.0
        enddo
      enddo

      stop
      end
