integer id
character st1*16
real lat,lon,ans1(100),ans2(100),ans3(100)

open(11,file='site_info.dat')
open(12,file='tmp.dat')
open(13,file='prc.dat')
open(14,file='hum.dat')
open(22,file='xtmp.dat')
open(23,file='xprc.dat')
open(24,file='xhum.dat')

read(11,*)
10 read(11,*,end=20) st1,id
read(12,*) lat,lon,ans1
read(13,*) lat,lon,ans2
read(14,*) lat,lon,ans3
if (id == 3200) then
  write(22,'(f7.3,f9.3,100f9.2)') lat,lon,ans1
  write(23,'(f7.3,f9.3,100f9.1)') lat,lon,ans2
  write(24,'(f7.3,f9.3,100f9.2)') lat,lon,ans3
endif
goto 10
20 continue

close(11)
close(12)
close(13)
close(14)
close(22)
close(23)
close(24)


stop
end
