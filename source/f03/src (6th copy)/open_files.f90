module open_files

use real_precision
use func
use pft_parameters

implicit none

contains

!**********************************************************************!
!                                                                      !
!                     open_diag :: open_files                          !
!                     -----------------------                          !
!                                                                      !
! subroutine open_diag(stoutput)                                       !
!----------------------------------------------------------------------!
!> @brief Open diagnostic file.
!! @details This is used to report any potential errors.
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine open_diag(stoutput)
!**********************************************************************!
character :: stoutput*1000
integer :: kode
!----------------------------------------------------------------------!

open(11,file=stoutput(1:blank(stoutput))//'/diag.dat',iostat=kode)
if (kode/=0) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'SDGVM output directory does not exist.'
  write(*,'('' "'',A,''"'')') stoutput(1:blank(stoutput))
  stop
endif
write(11,'(''ERRORS'')')

end subroutine open_diag





!**********************************************************************!
!                                                                      !
!                     open_snapshots :: open_files                     !
!                     ----------------------------                     !
!                                                                      !
! subroutine open_snapshots(stoutput,snp_no,snpshts)                   !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Open snapshot files file.
!! @details These files contain the entire system state for the final
!! day of the years given in the input files. The files can be used
!! the make site 'cartoons' using the utility sorfware.
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine open_snapshots(stoutput,snp_no,snpshts)
!----------------------------------------------------------------------!
character ::stoutput*1000,st1*1000,st4*1000
integer :: i,fno,snp_no,snpshts(1000)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Open snapshot files.                                                 !
!----------------------------------------------------------------------!
if (snp_no>0) then
  do i=1,snp_no
    fno = 100 + (i-1)*4
    st1=in2st(snpshts(i))
    call STRIPB(st4)
    open(fno+1,file=stoutput(1:blank(stoutput))//'/initbio_'&
 //st1(1:blank(stoutput))//'.dat')
    open(fno+2,file=stoutput(1:blank(stoutput))//'/initcov_'&
 //st1(1:blank(stoutput))//'.dat')
    open(fno+3,file=stoutput(1:blank(stoutput))//'/initppm_'&
 //st1(1:blank(stoutput))//'.dat')
    open(fno+4,file=stoutput(1:blank(stoutput))//'/inithgt_'&
 //st1(1:blank(stoutput))//'.dat')
  enddo
endif

end subroutine open_snapshots






!**********************************************************************!
!                                                                      !
!                     open_default :: open_files                       !
!                     --------------------------                       !
!                                                                      !
! subroutine open_default(stoutput)                                    !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Open default output files.
!! @details These files  will contain yearly values.
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine open_default(stoutput)
!**********************************************************************!
character :: stoutput*1000
!----------------------------------------------------------------------!

open(21,file=stoutput(1:blank(stoutput))//'/lai.dat')
open(22,file=stoutput(1:blank(stoutput))//'/npp.dat')
open(23,file=stoutput(1:blank(stoutput))//'/scn.dat')
open(24,file=stoutput(1:blank(stoutput))//'/snn.dat')
open(25,file=stoutput(1:blank(stoutput))//'/nep.dat')
open(26,file=stoutput(1:blank(stoutput))//'/swc.dat')
open(27,file=stoutput(1:blank(stoutput))//'/biomass.dat')
open(28,file=stoutput(1:blank(stoutput))//'/bioind.dat')
open(29,file=stoutput(1:blank(stoutput))//'/covind.dat')
open(30,file=stoutput(1:blank(stoutput))//'/dof.dat')
open(31,file=stoutput(1:blank(stoutput))//'/rof.dat')
open(32,file=stoutput(1:blank(stoutput))//'/fcn.dat')
open(33,file=stoutput(1:blank(stoutput))//'/nppstore.dat')
open(34,file=stoutput(1:blank(stoutput))//'/stembio.dat')
open(35,file=stoutput(1:blank(stoutput))//'/rootbio.dat')
open(36,file=stoutput(1:blank(stoutput))//'/leafper.dat')
open(37,file=stoutput(1:blank(stoutput))//'/stemper.dat')
open(38,file=stoutput(1:blank(stoutput))//'/rootper.dat')
open(39,file=stoutput(1:blank(stoutput))//'/sresp.dat')
open(40,file=stoutput(1:blank(stoutput))//'/evt.dat')
open(41,file=stoutput(1:blank(stoutput))//'/gpp.dat')
open(42,file=stoutput(1:blank(stoutput))//'/lch.dat')
open(43,file=stoutput(1:blank(stoutput))//'/prc.dat')
open(44,file=stoutput(1:blank(stoutput))//'/nbp.dat')
open(45,file=stoutput(1:blank(stoutput))//'/trn.dat')
open(46,file=stoutput(1:blank(stoutput))//'/fab.dat')
open(47,file=stoutput(1:blank(stoutput))//'/tmp.dat')
open(48,file=stoutput(1:blank(stoutput))//'/hum.dat')

end subroutine open_default





!**********************************************************************!
!                                                                      !
!                     open_optional :: open_files                      !
!                     ----------------------------                     !
!                                                                      !
!                                           !
! subroutine open_optional(stoutput,nft,out_cov,out_bio,out_bud,       !
! out_sen)                                                             !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Open optional output files.
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine open_optional(stoutput,nft,out_cov,out_bio,out_bud,out_sen)
!**********************************************************************!
character :: stoutput*1000
integer :: ft,nft,iofn
logical :: out_cov,out_bio,out_bud,out_sen

!----------------------------------------------------------------------!
! Open optional yearly cover, biomass, budburst and senescence files.  !
!----------------------------------------------------------------------!
iofn = 200
do ft=1,nft
  if (out_cov) then
    iofn = iofn + 1
    open(iofn,file=stoutput(1:blank(stoutput))//'/cov_'&
 //pft_tab(ft)%tag(1:blank(pft_tab(ft)%tag))//'.dat')
  endif
  if (out_bio) then
    iofn = iofn + 1
    open(iofn,file=stoutput(1:blank(stoutput))//'/bio_'&
 //pft_tab(ft)%tag(1:blank(pft_tab(ft)%tag))//'.dat')
  endif
  if (out_bud) then
    iofn = iofn + 1
    open(iofn,file=stoutput(1:blank(stoutput))//'/bud_'&
 //pft_tab(ft)%tag(1:blank(pft_tab(ft)%tag))//'.dat')
  endif
  if (out_sen) then
    iofn = iofn + 1
    open(iofn,file=stoutput(1:blank(stoutput))//'/sen_'&
 //pft_tab(ft)%tag(1:blank(pft_tab(ft)%tag))//'.dat')
  endif
enddo

end subroutine open_optional





!**********************************************************************!
!                                                                      !
!                     open_state :: open_files                         !
!                     ------------------------                         !
!                                                                      !
! subroutine open_state(stinit,stoutput,initise,initiseo)              !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Open input and/or output files for state vector.
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine open_state(stinit,stoutput,initise,initiseo)
!**********************************************************************!
character stinit*1000,stoutput*1000
logical initise,initiseo
integer kode

!----------------------------------------------------------------------!
! Open input files for state vector.                                   !
!----------------------------------------------------------------------!
if (.not.initise) then
  open(70,file=stinit(1:blank(stinit))//'/init.dat',status='old',&
 iostat=kode)
  if (kode/=0) then
    write(*,'('' PROGRAM TERMINATED'')')
    write(*,*) 'Initialisation directory/file does not exist.'
    write(*,'('' "'',A,''"'')') stinit(1:blank(stinit))
    stop
  endif
  open(71,file=stinit(1:blank(stinit))//'/initbio.dat')
  open(72,file=stinit(1:blank(stinit))//'/initcov.dat')
  open(73,file=stinit(1:blank(stinit))//'/initppm.dat')
  open(74,file=stinit(1:blank(stinit))//'/inithgt.dat')
  open(76,file=stinit(1:blank(stinit))//'/initleaf.dat')
  open(77,file=stinit(1:blank(stinit))//'/initstem.dat')
  open(78,file=stinit(1:blank(stinit))//'/initroot.dat')
  open(79,file=stinit(1:blank(stinit))//'/initmisc.dat')
endif

!----------------------------------------------------------------------!
! Open output files for state vector.                                  !
!----------------------------------------------------------------------!
open(80,file=stoutput(1:blank(stoutput))//'/init.dat')
if (initiseo) then
  open(81,file=stoutput(1:blank(stoutput))//'/initbio.dat')
  open(82,file=stoutput(1:blank(stoutput))//'/initcov.dat')
  open(83,file=stoutput(1:blank(stoutput))//'/initppm.dat')
  open(84,file=stoutput(1:blank(stoutput))//'/inithgt.dat')
  open(86,file=stoutput(1:blank(stoutput))//'/initleaf.dat')
  open(87,file=stoutput(1:blank(stoutput))//'/initstem.dat')
  open(88,file=stoutput(1:blank(stoutput))//'/initroot.dat')
  open(89,file=stoutput(1:blank(stoutput))//'/initmisc.dat')
endif

end subroutine open_state





!**********************************************************************!
!                                                                      !
!                     write_lat_lon :: open_files                      !
!                     ---------------------------                      !
!                                                                      !
!                                     !
! subroutine WRITE_LAT_LON(lat,lon,nft,out_cov,out_bio,out_bud,        !
! out_sen,oymdft,otagsn,oymd,otagsnft,outyears1,outyears2)             !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Write lat & lon for output files.
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine write_lat_lon(lat,lon,nft,out_cov,out_bio,out_bud,out_sen, &
 oymdft,otagsn,oymd,otagsnft,outyears1,outyears2)
!**********************************************************************!
integer :: i,ft,nft,outyears1,outyears2,iofn,oymdft,otagsn(50),oymd, &
 otagsnft(50)
real(dp):: lat,lon
logical :: out_cov,out_bio,out_bud,out_sen

!----------------------------------------------------------------------!
! Write lat lon for default output files.                              !
!----------------------------------------------------------------------!

do i=21,48
  write(i,'(f7.3,f9.3)',advance='NO') lat,lon
enddo

!----------------------------------------------------------------------!
! Write lat lon for optional yearly ft files                           !
!----------------------------------------------------------------------!
iofn = 200
if (outyears2>0) then
  do ft=1,nft
    if (out_cov) then
      iofn = iofn + 1
      write(iofn,'(f7.3,f9.3)') lat,lon
    endif
    if (out_bio) then
      iofn = iofn + 1
      write(iofn,'(f7.3,f9.3)') lat,lon
    endif
    if (out_bud) then
      iofn = iofn + 1
      write(iofn,'(f7.3,f9.3)') lat,lon
    endif
    if (out_sen) then
      iofn = iofn + 1
      write(iofn,'(f7.3,f9.3)') lat,lon
    endif
  enddo
endif

!----------------------------------------------------------------------!
! Write lat,lon for monthly/daily output files.                        !
!----------------------------------------------------------------------!
iofn = 400
if (outyears1>0) then
  if ((oymd==1).or.(oymd==3)) then
    do i=1,50
      if (otagsn(i)==1) then
        iofn = iofn + 1
        write(iofn,'(f7.3,f9.3)') lat,lon
      endif
    enddo
  endif
  if ((oymd==2).or.(oymd==3)) then
    do i=1,50
      if (otagsn(i)==1) then
        iofn = iofn + 1
        write(iofn,'(f7.3,f9.3)') lat,lon
      endif
    enddo
  endif
endif

!----------------------------------------------------------------------!
! Write lat,lon for monthly/daily ft output files.                     !
!----------------------------------------------------------------------!
if (outyears2>0) then
  if ((oymdft==1).or.(oymdft==3)) then
    do i=1,50
      if (otagsnft(i)==1) then
        do ft=1,nft
          iofn = iofn + 1
          write(iofn,'(f7.3,f9.3)') lat,lon
       enddo
      endif
    enddo
  endif
  if ((oymdft==2).or.(oymdft==3)) then
    do i=1,50
      if (otagsnft(i)==1) then
        do ft=1,nft
          iofn = iofn + 1
          write(iofn,'(f7.3,f9.3)') lat,lon
        enddo
      endif
    enddo
  endif
endif

end subroutine write_lat_lon



end module open_files



