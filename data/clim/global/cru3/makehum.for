      implicit none
      integer vap(12),tmp(12),hum(12),i,k
      real t,v

      open(11,file='30min\tmp.dat',status='old')
      open(12,file='30min\vap.dat',status='old')
      open(21,file='30min\hum.dat')

10      read(11,'(12i4)',end=20) tmp
        read(12,'(12i4)') vap
        do i=1,12
          t = real(tmp(i))/10.0
          v = real(vap(i))*10000.0
          hum(i) = int(v/(613.75*exp(17.502*t/(240.97+t)))+0.5)
          if (hum(i).gt.1000)  hum(i) = 1000
          if (hum(i).lt.0)  hum(i) = 0
        enddo
        write(21,'(12i4)') hum
      goto 10
20    continue


      stop
      end