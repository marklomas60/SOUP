program test1

implicit none
integer :: xx(10,5),i,j,y(10,5),z,x(7),row,col,xcount,countAll,xall(360*120),maxAns
integer(kind=1) :: ans
character(len=10) :: str
character :: outc(6172)

ans = 0
print*,ans
ans = 127
print*,ans
print*,ans
ans = -128
print*,ans
print*,ans

!do i=1,5
!  do j=1,10
!    x(j,i) = (i-1)*10 + j
!  enddo
!enddo

!open(11,file='temp.dat',form='unformatted',access='stream')
!write(11) x(:,1)
!write(11) x(:,2)
!write(11) x(:,3)
!write(11) x(:,4)
!write(11) x(:,5)
!write(11) x
!close(11)
!print*,x

!open(11,file='temp.dat',access='stream')
!read(11) y
!close(11)
!print*,y

open(11,file='../../../data/land_mask.dat',status='old',form='formatted')
open(21,file='../../../data/land_mask.bin',form='unformatted',access='stream',status='replace')

xcount = 0
countAll = 0
maxAns = -1000

do row=1,180*120
!do row=1,1
read(11,'(6172a)') (outc(j),j=1,6172)
print*,row,ans
do col=1,6172
  call base72i(outc(col),x)
  if (col==6172) then
    xall(col*7-6:7*col-4) = x(1:3)
  else
    xall(col*7-6:7*col) = x
  endif
enddo

do col=1,5400
  ans = xall(col*8-7)*128+xall(col*8-6)*64+xall(col*8-5)*32+xall(col*8-4)*16+xall(col*8-3)*8+ &
  xall(col*8-2)*4+xall(col*8-1)*2+xall(col*8) - 128
  write(21) ans
  if (ans>maxAns) then
    maxAns = ans
	print*,maxAns
  endif
enddo
enddo

!  do j=1,7
!    countAll = countAll + 1
!    if (x(j)==1) then
!      print*,.true.
!      write(21) .true.
!	  xcount = xcount + 1
!	else
!      write(21) .false.
!	endif
!  enddo
!print*,countAll,xcount

close(11)
close(21)

!print*,real(xcount)/real(countAll)*100.0
!print*,xall

end program test1


subroutine tt()

print*,'hellox'

end subroutine tt



!**********************************************************************!
!                                                                      !
!                          base72I :: data                             !
!                          ---------------                             !
!                                                                      !
! subroutine base72i(c,x)                                              !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine base72i(c,x)
!**********************************************************************!
character :: c
integer :: x(7),i
!----------------------------------------------------------------------!

i=ichar(c)-100
x(1)=i/64
i=i-x(1)*64
x(2)=i/32
i=i-x(2)*32
x(3)=i/16
i=i-x(3)*16
x(4)=i/8
i=i-x(4)*8
x(5)=i/4
i=i-x(5)*4
x(6)=i/2
i=i-x(6)*2
x(7)=i

end subroutine base72i





!**********************************************************************!
!                                                                      !
!                          n7 :: data                                  !
!                          ----------                                  !
!                                                                      !
! Returns the value of a seven digit binary number given as a seven    !
! dimensional array of 0's and 1's                                     !
!                                                                      !
!    integer function n7(x)                                            !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
integer function n7(x)
!----------------------------------------------------------------------!
integer :: x(7)
!----------------------------------------------------------------------!

n7 = 64*x(1)+32*x(2)+16*x(3)+8*x(4)+4*x(5)+2*x(6)+x(7)
 
!**********************************************************************!
end function n7


