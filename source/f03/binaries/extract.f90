program extract
! 21600 rows, 5400 bytes (43200 bits), cols
implicit none
integer, parameter :: bitRows = 21600, bitCols = 43200
real :: lat,lon,del, latf, lon0, londel, latdel
integer :: n, col, row,bit,byte,rByte,bitPos,rr
integer, dimension(8) :: out
integer(kind=1) :: r

latdel = 0.5
londel = 0.5

lat = 53.0
lon = -5.25

del = 1.0/60.0/2.0

latf = 90.0d0 - del/2.0d0
lon0 =-180.0d0 + del/2.0d0

col = int((lon - lon0)/del + 0.5d0)
row = int((latf - lat)/del + 0.5d0)

print*,'row col ',row,col
print*,latf,lon,del

n = min((latdel/del-1.0d0)/2.0d0,(londel/del-1.0d0)/2.0d0)

open(11,file='../../../data/land_mask.bin',form='unformatted',access='stream')

bit = (row-1)*bitCols+col
byte = (bit-1)/8
bitPos =  bit - byte*8

print*,'bit byte bitRem ',bit, byte, bitPos

read(11,pos=byte) r
rr = r + 128
print*,r,rr
close(11)

call byte2bin(rr,out)
print*,out(bitPos)

end program extract


subroutine byte2bin(inByte,bin)
integer, intent(in) :: inByte
integer, intent(out), dimension(8) :: bin

integer :: byte
!----------------------------------------------------------------------!

byte = inByte

bin(1)=byte/128
byte=byte-bin(1)*128
bin(2)=byte/64
byte=byte-bin(2)*64
bin(3)=byte/32
byte=byte-bin(3)*32
bin(4)=byte/16
byte=byte-bin(4)*16
bin(5)=byte/8
byte=byte-bin(5)*8
bin(6)=byte/4
byte=byte-bin(6)*4
bin(7)=byte/2
byte=byte-bin(7)*2
bin(8)=byte

end subroutine byte2bin


