integer iargc,blank,i,n_fields,no,sites,j
integer site,kode,k
character buff1*80,st1*1000,st2*1000
real*8 lat,lon,out(1000)

if (IARGC()>0) then
  call GETARG(1,buff1)
else
  write(*,*) 'Directory name must be given as an argument.'
  stop
endif
do i=1,80
  st1(i:i) = buff1(i:i)
enddo

open(21,file=st1(1:blank(st1))//'/lai.dat',status='old',iostat=kode)
if (kode/=0) then
  write(*,*) 'Directory path not found'
  write(*,'(A)') st1(1:blank(st1))
  stop
endif
sites = 0
30    read(21,*,end=20) lat,lon
  sites = sites + 1
  no = -1
10      continue
  read(21,'(A)') st2
  no = no + 1
  if (n_fields(st2)==1) goto 10
goto 30
20    continue

write(*,*) no,sites
close(21)

open(21,file=st1(1:blank(st1))//'/lai.dat',status='old')
open(22,file=st1(1:blank(st1))//'/npp.dat',status='old')
open(23,file=st1(1:blank(st1))//'/scn.dat',status='old')
open(24,file=st1(1:blank(st1))//'/snn.dat',status='old')
open(25,file=st1(1:blank(st1))//'/nep.dat',status='old')
open(26,file=st1(1:blank(st1))//'/swc.dat',status='old')
open(27,file=st1(1:blank(st1))//'/biot.dat',status='old')
open(28,file=st1(1:blank(st1))//'/bioind.dat',status='old')
open(29,file=st1(1:blank(st1))//'/covind.dat',status='old')
open(30,file=st1(1:blank(st1))//'/dof.dat',status='old')
open(31,file=st1(1:blank(st1))//'/rof.dat',status='old')
open(32,file=st1(1:blank(st1))//'/fcn.dat',status='old')
open(33,file=st1(1:blank(st1))//'/nppstore.dat',status='old')
open(34,file=st1(1:blank(st1))//'/stembio.dat',status='old')
open(35,file=st1(1:blank(st1))//'/rootbio.dat',status='old')
open(36,file=st1(1:blank(st1))//'/sresp.dat',status='old')
open(37,file=st1(1:blank(st1))//'/evt.dat',status='old')
open(38,file=st1(1:blank(st1))//'/gpp.dat',status='old')
open(39,file=st1(1:blank(st1))//'/lch.dat',status='old')
open(40,file=st1(1:blank(st1))//'/prc.dat',status='old')
open(41,file=st1(1:blank(st1))//'/nbp.dat',status='old')
open(42,file=st1(1:blank(st1))//'/trn.dat',status='old')
open(43,file=st1(1:blank(st1))//'/fab.dat',status='old')

open(61,file=st1(1:blank(st1))//'/xlai.dat')
open(62,file=st1(1:blank(st1))//'/xnpp.dat')
open(63,file=st1(1:blank(st1))//'/xscn.dat')
open(64,file=st1(1:blank(st1))//'/xsnn.dat')
open(65,file=st1(1:blank(st1))//'/xnep.dat')
open(66,file=st1(1:blank(st1))//'/xswc.dat')
open(67,file=st1(1:blank(st1))//'/xbiot.dat')
open(68,file=st1(1:blank(st1))//'/xbioind.dat')
open(69,file=st1(1:blank(st1))//'/xcovind.dat')
open(70,file=st1(1:blank(st1))//'/xdof.dat')
open(71,file=st1(1:blank(st1))//'/xrof.dat')
open(72,file=st1(1:blank(st1))//'/xfcn.dat')
open(73,file=st1(1:blank(st1))//'/xnppstore.dat')
open(74,file=st1(1:blank(st1))//'/xstembio.dat')
open(75,file=st1(1:blank(st1))//'/xrootbio.dat')
open(76,file=st1(1:blank(st1))//'/xsresp.dat')
open(77,file=st1(1:blank(st1))//'/xevt.dat')
open(78,file=st1(1:blank(st1))//'/xgpp.dat')
open(79,file=st1(1:blank(st1))//'/xlch.dat')
open(80,file=st1(1:blank(st1))//'/xprc.dat')
open(81,file=st1(1:blank(st1))//'/xnbp.dat')
open(82,file=st1(1:blank(st1))//'/xtrn.dat')
open(83,file=st1(1:blank(st1))//'/xfab.dat')

do j=1,23
  do site=1,sites
    read(j+20,*) lat,lon
    k = 0
    do i=1,no
      k=k+1
      write(*,*) j,site,i
      read(j+20,*) out(k)
    enddo
    write(j+60,'(f7.3,f9.3,1000f12.3)') lat,lon,(out(i),i=1,no)
  enddo
  close(j+20)
  close(j+60)
enddo


stop
end
