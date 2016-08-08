*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE PARAMETER_ADJUSTMENT             *
*                          *******************************             *
*                                                                      *
* This routine is invoked by the existance of the word 'ARGUMENT' in   *
* the input file after the soil parameter line.                        *
*                                                                      *
* The routine is used to read in parameters values from a file.        *
* The filename is given by the command line argument, two additional   *
* parameters are supplied in the command line which select the first   *
* and last records to be run.                                          *
* All parameters in 'param.inc', can be adjusted.                      *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE PARAMETER_ADJUSTMENT(param_file,n_param)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      CHARACTER param_file*1000
      INTEGER n_param,blank,i,kode
      REAL ans(17)
      INCLUDE 'param.inc'

      OPEN(98,FILE=param_file(1:blank(param_file)),STATUS='OLD',
     &iostat=kode)
      IF (kode.NE.0) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'Parameter adjustment file does not exist.'
        WRITE(*,'('' "'',A,''"'')') param_file(1:blank(param_file))
        STOP
      ENDIF

      IF (n_param.GT.1) THEN
        DO i=1,n_param-1
          READ(98,*)
        ENDDO
      ENDIF
      READ(98,*) ans
      CLOSE(98)

      p_sand = ans(1)
      p_silt = (100.0-p_sand)*ans(2)/100.0
      p_bulk = ans(3)
      p_orgc = ans(4)
      p_awl(1) = ans(5)
      p_awl(2) = ans(6)
      p_awl(3) = ans(7)
      p_awl(4) = ans(8)
      p_bc_res = ans(9)
      p_kx = ans(10)
      p_kd = ans(11)
      p_bs = ans(12)
      p_et = ans(13)
      p_roff = ans(14)
      p_roff2 = ans(15)
      p_dep = ans(16)
      p_topsl = ans(17)


      RETURN
      END
