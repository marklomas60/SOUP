**********************************************************************************
*                                                                                *
* Use backslashes not forward slashes for the input argument directory delimiter.*
*                                                                                *
********************************************************************************** 
      INTEGER nofiles
      PARAMETER(nofiles=27)

      INTEGER iargc,blank,i,n_fields,no,sites,j,year
      INTEGER site,kode
      CHARACTER buff1*80,st1*1000,st2*1000
      REAL lat(100000),lon(100000),out(100000,500)
      character filenames(nofiles)*100
      data filenames/'lai','npp','scn','snn','nep','swc','biot',
     &'bioind','covind','dof','rof','fcn','nppstore','stembio',
     &'rootbio','sresp','evt','gpp','lch','nbp','trn',
     &'fab','prc','hum','tmp','rootper','stemper'/

      IF (IARGC().GT.0) THEN
        CALL GETARG(1,buff1)
      ELSE
        WRITE(*,*) 'Directory name must be given as an argument.'
        STOP
      ENDIF
      DO i=1,80
        st1(i:i) = buff1(i:i)
      ENDDO

      print*,buff1(1:30)

      close(21)
      OPEN(21,FILE=st1(1:blank(st1))//'/lai.dat',status='old',
     &iostat=kode)
      IF (kode.ne.0) THEN
        WRITE(*,*) 'Directory path not found.'
        WRITE(*,'(A)') st1(1:blank(st1))
        STOP
      ENDIF
      sites = 0
30    READ(21,*,end=20) lat(1),lon(1)
        sites = sites + 1
        no = -1
10      CONTINUE
        READ(21,'(A)') st2
        no = no + 1
        IF (n_fields(st2).eq.1) GOTO 10
      GOTO 30
20    CONTINUE

      print*,'Years = ',no,',  Sites = ',sites
      close(21)

      do j=1,nofiles
        print*,j
        print*,st1(1:blank(st1))//
     &'\'//trim(filenames(j))//'.dat'
        open(j+20,file=st1(1:blank(st1))//
     &'\'//trim(filenames(j))//'.dat',status='old')
        do site=1,sites
          read(j+20,*) lat(site),lon(site)
          do year=1,no
            read(j+20,*) out(site,year)
          enddo
        enddo
        close(j+20)

        open(j+20,file=st1(1:blank(st1))//
     &'\'//trim(filenames(j))//'.dat')
        do site=1,sites
          write(j+20,'(f7.3,f9.3,500E14.6)') lat(site),lon(site),
     &(out(site,year),year=1,no)
        enddo
        close(j+20)
      enddo

      
      STOP
      END



*----------------------------------------------------------------------*
*                                                                      *
*                          FUNCTION blank                              *
*                          **************                              *
*                                                                      *
* 'blank' returns the number of characters before the space in the     *
* string which is its argument.                                        *
*                                                                      *
*----------------------------------------------------------------------*
      FUNCTION blank(st1)
*----------------------------------------------------------------------*
      CHARACTER st1*1000
      INTEGER blank

      blank = 0
10    CONTINUE
      blank = blank + 1
      IF ((blank.LE.100).AND.(ichar(st1(blank:blank)).NE.32)
     &     .AND.(ichar(st1(blank:blank)).NE.9))  GOTO 10
c space and tabulation
      blank = blank - 1


      RETURN
      END




*----------------------------------------------------------------------*
*                                                                      *
*                          FUNCTION n_fields                           *
*                          *****************                           *
*                                                                      *
* Reurns the number of space delimeted fields.                         *
*                                                                      *
*----------------------------------------------------------------------*
      FUNCTION n_fields(st1)
*----------------------------------------------------------------------*
      CHARACTER st1*1000,st2*1000,st3*1000
      INTEGER blanks,n_fields

      n_fields = 0
      st2 = st1

10    CONTINUE
      IF (blanks(st2).LT.1000) THEN
        CALL STRIPBS(st2,st3)
        n_fields = n_fields + 1
        GOTO 10
      ENDIF


      RETURN
      END


*----------------------------------------------------------------------*
*                                                                      *
*                          STRIPBS                                     *
*                          *******                                     *
*                                                                      *
* Strips leading blanks and first string from the string which is its  *
* argument.                                                            *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE STRIPBS(st1,st2)
      CHARACTER st1*1000,st2*1000
      INTEGER blanks,nb,i,blank

      nb = blanks(st1)
      DO i=1,1000-nb
        st1(i:i) = st1(i+nb:i+nb)
      ENDDO
      DO i=1000-nb+1,1000
        st1(i:i) = ' '
      ENDDO

      nb = blank(st1)
      DO i=1,nb
        st2(i:i) = st1(i:i)
      ENDDO
      st2(nb+1:nb+1) = ' '

      IF (nb.GT.0) THEN
        DO i=1,1000-nb
          st1(i:i) = st1(i+nb:i+nb)
        ENDDO
        DO i=1000-nb+1,1000
          st1(i:i) = ' '
        ENDDO
      ENDIF


      RETURN
      END



*----------------------------------------------------------------------*
*                                                                      *
*                          FUNCTION blanks                             *
*                          ***************                             *
*                                                                      *
* 'blanks' returns the number of leading blanks.                       *
*                                                                      *
*----------------------------------------------------------------------*
      FUNCTION blanks(st1)
*----------------------------------------------------------------------*
      CHARACTER st1*1000
      INTEGER blanks

      blanks = 0
10    CONTINUE
      blanks = blanks + 1
      IF (blanks.LE.1000) THEN
        IF ((ichar(st1(blanks:blanks)).EQ.32).OR. 
     &(ichar(st1(blanks:blanks)).EQ.9)) GOTO 10
      ENDIF
      blanks = blanks - 1


      RETURN
      END







