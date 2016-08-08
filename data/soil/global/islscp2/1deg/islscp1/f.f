      real*8 lat,lon,sand,silt,bulk,dep,orgc,wilt,field,satw
      integer i

      open(11,file='data_old.dat')
      open(21,file='data.dat')

      wilt = -999.0
      orgc = -999.0
      field = -999.0
      satw = -999.0
      bulk = -999.0

      do i=1,64800
        read(11,*) lat,lon,sand,silt,bulk,dep
        if (dep.lt.0.0000001) dep = -999.0d0
        write(21,'(f7.3,f9.3,14f10.4)') lat,lon,sand,silt,
     &bulk,orgc,wilt,field,satw,dep
      enddo



      stop
      end
