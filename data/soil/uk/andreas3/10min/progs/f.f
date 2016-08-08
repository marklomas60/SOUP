      integer i,j,k,id
      real lat,lon,c1,c30,org1,org30,cl30,sa30,si30,bulk30
      real cl1,sa1,si1,bulk1,map(8,60,60),lat0,lon0,ca,cc

      lat0 = 50.0833
      lon0 = -7.0833

      do i=1,60
        do j=1,60
          do k=1,8
            map(k,i,j) = -999
          enddo
        enddo
      enddo

      open(11,file='CTCD10minUKcarbtext07.txt',status='old')
      open(21,file='../30cm/data.dat')
      open(22,file='../1m/data.dat')
      read(11,*)
20    read(11,*,end=10) lat,lon,id,ca,cc,bulk30,bulk1,org30,
     &org1,cl30,sa30,si30,cl1,sa1,si1
        if (sa30.LT.200) then
        ind1 = int((lat - lat0)*6.0 + .5)+1
        ind2 = int((lon - lon0)*6.0 + .5)+1
        map(1,ind1,ind2) = sa30/(sa30+si30+cl30)*100.0
        map(2,ind1,ind2) = si30/(sa30+si30+cl30)*100.0
        map(3,ind1,ind2) = bulk30
        map(4,ind1,ind2) = org30

        map(5,ind1,ind2) = sa1/(sa1+si1+cl1)*100.0
        map(6,ind1,ind2) = si1/(sa1+si1+cl1)*100.0
        map(7,ind1,ind2) = bulk1
        map(8,ind1,ind2) = org1

        if ((ind1.gt.60).or.(ind1.lt.1).or.(ind2.gt.60).or.
     &(ind2.lt.1)) then
          print*,'error index gone bad'
          stop
        endif
        endif
      goto 20
10    continue

      do i=1,60
        do j=1,60
          write(21,'(f7.3,f9.3,8f10.4)') real(60-i)/6.0+lat0,
     &real(j-1)/6.0+lon0,(map(k,61-i,j),k=1,4),
     &-999.0,-999.0,-999.0,-999.0
          write(22,'(f7.3,f9.3,8f10.4)') real(60-i)/6.0+lat0,
     &real(j-1)/6.0+lon0,(map(k,61-i,j),k=5,8),
     &-999.0,-999.0,-999.0,-999.0
        enddo
      enddo

      stop
      end
