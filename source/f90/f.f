      real lat,lon,ans1(200,2,9),ans2(200,9)
      integer i,j,k,ftmor(9),nft
      data ftmor/1, 1, 2, 2, 200, 200, 200, 200, 1/

      nft=9

      open(71,file='/data/SDGVM/andreas/spin/old_initcov.dat',
     &status='old')
      open(81,file='/data/SDGVM/andreas/spin/initcov.dat')
      DO site=1,42000
      READ(71,'(F7.3,F9.3,4000E13.6)') lat,lon,
     &((ans2(i,j),i=1,ftmor(j)),j=2,nft)
      WRITE(81,'(F7.3,F9.3,4000E13.6)') lat,lon,
     &((max(ans2(i,j),0.0),i=1,ftmor(j)),j=2,nft)
      ENDDO
      close(71)
      close(81)

      open(71,file='/data/SDGVM/andreas/spin/old_initbio.dat',
     &status='old')
      open(81,file='/data/SDGVM/andreas/spin/initbio.dat')
      DO site=1,42000
      READ(71,'(F7.3,F9.3,4000E13.6)') lat,lon,
     &(((ans1(i,k,j),k=1,2),i=1,ftmor(j)),j=2,nft)
      WRITE(81,'(F7.3,F9.3,4000E13.6)') lat,lon,
     &(((max(ans1(i,k,j),0.0),k=1,2),i=1,ftmor(j)),j=2,nft)
      ENDDO
      close(71)
      close(81)

      open(71,file='/data/SDGVM/andreas/spin/old_initppm.dat',
     &status='old')
      open(81,file='/data/SDGVM/andreas/spin/initppm.dat')
      DO site=1,42000
      READ(71,'(F7.3,F9.3,4000E13.6)') lat,lon,
     &((ans2(i,j),i=1,ftmor(j)),j=2,nft)
      WRITE(81,'(F7.3,F9.3,4000E13.6)') lat,lon,
     &((max(ans2(i,j),0.0),i=1,ftmor(j)),j=2,nft)
      ENDDO
      close(71)
      close(81)

      open(71,file='/data/SDGVM/andreas/spin/old_inithgt.dat',
     &status='old')
      open(81,file='/data/SDGVM/andreas/spin/inithgt.dat')
      DO site=1,42000
      READ(71,'(F7.3,F9.3,4000E13.6)') lat,lon,
     &((ans2(i,j),i=1,ftmor(j)),j=2,nft)
      WRITE(81,'(F7.3,F9.3,4000E13.6)') lat,lon,
     &((max(ans2(i,j),0.0),i=1,ftmor(j)),j=2,nft)
      ENDDO
      close(71)
      close(81)

      open(71,file='/data/SDGVM/andreas/spin/old_initwdt.dat',
     &status='old')
      open(81,file='/data/SDGVM/andreas/spin/initwdt.dat')
      DO site=1,42000
      READ(71,'(F7.3,F9.3,4000E13.6)') lat,lon,
     &((ans2(i,j),i=1,ftmor(j)),j=2,nft)
      WRITE(81,'(F7.3,F9.3,4000E13.6)') lat,lon,
     &((max(ans2(i,j),0.0),i=1,ftmor(j)),j=2,nft)
      ENDDO
      close(71)
      close(81)

      stop
      end
