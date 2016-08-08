      integer i,j,ans,xans(259200),xmin,xmax,xmini,xmaxi
      integer count1,count2

      open(11,file='../cont_lu-1-2000.dat',status='old')
      open(12,file='../cont_lu-2-2000.dat',status='old')
      open(13,file='../cont_lu-3-2000.dat',status='old')
      open(14,file='../cont_lu-4-2000.dat',status='old')
      open(15,file='../cont_lu-5-2000.dat',status='old')
      open(16,file='../cont_lu-6-2000.dat',status='old')
      open(17,file='../cont_lu-7-2000.dat',status='old')
      open(18,file='../cont_lu-8-2000.dat',status='old')
      open(19,file='../cont_lu-9-2000.dat',status='old')
      open(20,file='../cont_lu-10-2000.dat',status='old')
      open(21,file='../cont_lu-11-2000.dat',status='old')
      open(22,file='../cont_lu-12-2000.dat',status='old')
      open(23,file='../cont_lu-13-2000.dat',status='old')
      open(24,file='../cont_lu-14-2000.dat',status='old')
      open(25,file='../cont_lu-15-2000.dat',status='old')
      open(26,file='../cont_lu-16-2000.dat',status='old')
      open(27,file='../cont_lu-17-2000.dat',status='old')
      open(28,file='../cont_lu-18-2000.dat',status='old')
      open(29,file='../cont_lu-19-2000.dat',status='old')
      open(30,file='../cont_lu-20-2000.dat',status='old')
      open(31,file='../cont_lu-21-2000.dat',status='old')
      open(32,file='../cont_lu-22-2000.dat',status='old')
      open(33,file='../cont_lu-23-2000.dat',status='old')
      open(34,file='../cont_lu-24-2000.dat',status='old')
      open(35,file='../cont_lu-25-2000.dat',status='old')
      open(36,file='../cont_lu-26-2000.dat',status='old')
      open(37,file='../cont_lu-27-2000.dat',status='old')
      open(38,file='../cont_lu-28-2000.dat',status='old')
      open(39,file='../cont_lu-29-2000.dat',status='old')
      open(40,file='../cont_lu-30-2000.dat',status='old')
      open(41,file='../cont_lu-31-2000.dat',status='old')

      do j=1,259200
        xans(j) = 0
      enddo

      do i=11,41
        print*,i
        do j=1,22270
          read(i,*) ans
          if (ans.LT.255) xans(j) = xans(j) + ans
          if (ans.LT.0) print*,i,j,ans,' xxxxxxxxxxxxxxxxxxxxxxx'
*          if (j.eq.98506) print*,i,ans
        enddo
      enddo

      count1 = 0
      count2 = 0
      xmin = 200
      xmax = 0
      do j=1,22270
        if ((xans(j).NE.100).and.(xans(j).GT.0)) then
          print*,j,xans(j)
          count1 = count1 + 1
          if (xans(j).LT.xmin) then
            xmin = xans(j)
            xmini = j
          endif
          if (xans(j).GT.xmax) then
            xmax = xans(j)
            xmaxi = j
          endif
        else
          count2 = count2 + 1
        endif
      enddo

      print*,xmin,xmini,xmax,xmaxi
      print*,count1,count2

      stop
      end
