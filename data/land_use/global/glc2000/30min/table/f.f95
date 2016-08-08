PROGRAM f

IMPLICIT NONE
INTEGER, DIMENSION(23,720,292) :: c1
INTEGER, DIMENSION(12) :: x
INTEGER :: i,j,k
CHARACTER :: ch*2

!----------------------------------------------------------------------*
!Read in inundated fraction map.
!----------------------------------------------------------------------*
DO k=1,23
  write(ch,'(i0)') k
  OPEN(21,file='..\cont_lu-'//trim(ch)//'-2000.dat',status='old')
  READ(21,*) ((c1(k,i,j),i=1,720),j=1,292)
  CLOSE(21)
ENDDO

OPEN(61,file='classtable.dat')

k=3
x = 0
DO j=1,292
  DO i=1,720
    if (c1(1,i,j).lt.101) then
      write(61,'(f7.3,f9.3,2x,23f6.2)') 90.25-real(j)/2.0,real(i)/2.0-180.25,(real(c1(k,i,j))/100.0,k=1,23)
    endif
  ENDDO
ENDDO

CLOSE(61)

END PROGRAM f