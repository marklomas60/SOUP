      real lat,lon
      integer ans(7),i,j,blank,x,y,xans(18)
      character st1(18)*6,st2(7)*6

      st2(1) = '1700'
      st2(2) = '1750'
      st2(3) = '1800'
      st2(4) = '1850'
      st2(5) = '1900'
      st2(6) = '1950'
      st2(7) = '1990'

      st1(1) = '1'
      st1(2) = '2'
      st1(3) = '3'
      st1(4) = '4'
      st1(5) = '5'
      st1(6) = '6'
      st1(7) = '7'
      st1(8) = '8'
      st1(9) = '9'
      st1(10) = '10'
      st1(11) = '11'
      st1(12) = '12'
      st1(13) = '13'
      st1(14) = '14'
      st1(15) = '15'
      st1(16) = '16'
      st1(17) = '17'
      st1(18) = '18'

      open(11,file='landuse.dat')

*----------------------------------------------------------------------*

      do i=1,18
        do j=1,7
          open(21+(i-1)*7+(j-1),file='../cont_lu-'//st1(i)(1:blank(st1(i
     &)))//'-'//st2(j)(1:blank(st2(j)))//'.dat')
        enddo
      enddo

      do i=1,18
        xans(i) = 0
      enddo

      do x=1,360
        print*,x
        do y=1,720
          read(11,*) lat,lon,ans
        IF (ans(1)+ans(2)+ans(3)+ans(4)+ans(5)+ans(6)+ans(7).gt.0) THEN
            do j=1,7
              xans(ans(j)) = 100
              do i=1,18
                write(21+(i-1)*7+(j-1),'(i3)') xans(i)
              enddo
              xans(ans(j)) = 0
            enddo
          ELSE
            do j=1,7
              do i=1,18
                write(21+(i-1)*7+(j-1),'(i3)') 255
              enddo
            enddo
          ENDIF
        enddo
      enddo


      stop
      end

*----------------------------------------------------------------------*
      FUNCTION minx(x)
*----------------------------------------------------------------------*
      real*8 x(12),minxx
      integer minx,i

      minx = -10
      minxx = 100000000.0d0
      DO i=1,12
        IF (x(i).LT.minxx) then
          minxx = x(i)
          minx = i
        ENDIF
      ENDDO


      RETURN
      END

*----------------------------------------------------------------------*
      FUNCTION dfn(age,mature)
*----------------------------------------------------------------------*
      INTEGER age,lsd
      PARAMETER(lsd = 200)
      REAL*8 mature,dfn

      IF (age.GT.mature) THEN
        dfn = 1.0d0 - 0.5d0*(real(age) - mature)/(real(lsd) -
     &mature)
      ELSE
        dfn = (real(age)/mature)**2.0d0
      ENDIF


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          FUNCTION ttg                                *
*                          ************                                *
*                                                                      *
*----------------------------------------------------------------------*
      FUNCTION ttg(i,ttglim)
*----------------------------------------------------------------------*
      INTEGER i,ttglim
      REAL*8 ttg

      ttg = real(ttglim + 1 - i)/real(ttglim)
      IF (ttg.LT.0.0d0)  ttg = 0.0d0

      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          FUNCTION triscale                           *
*                          *****************                           *
*                                                                      *
*----------------------------------------------------------------------*
      FUNCTION triscale(x,y,x1,x2,x3,y1,y2,y3)
*----------------------------------------------------------------------*
      REAL*8 triscale,x,y,x1,x2,x3,y1,y2,y3,s1,s2

      IF ((x.GT.x1).OR.(x.LT.x3).OR.(y.LT.y1).OR.(y.GT.y3)) THEN
        triscale = 0.0d0
      ELSE
        IF (x.GT.x2) THEN
          s1 = (x - x1)/(x2 - x1)
        ELSE
          s1 = (x - x3)/(x2 - x3)
        ENDIF
        IF (y.LT.y2) THEN
          s2 = (y - y1)/(y2 - y1)
        ELSE
          s2 = (y - y3)/(y2 - y3)
        ENDIF
        triscale = s1*s2
      ENDIF


      RETURN
      END


*----------------------------------------------------------------------*
*                                                                      *
*                          FUNCTION no_days                            *
*                          ****************                            *
*                                                                      *
*----------------------------------------------------------------------*
      FUNCTION no_days(year,mnth,thty_dys)
*----------------------------------------------------------------------*
      INTEGER no_days,year,mnth,thty_dys

      IF (thty_dys.EQ.0) THEN
        IF (mnth.EQ.2) THEN
          IF (mod(year,4).EQ.0) THEN
            no_days = 29
          ELSE
            no_days = 28
          ENDIF
          IF (mod(year,100).EQ.0) no_days = 28
          IF (mod(year,400).EQ.0) no_days = 29
        ELSE
          IF ((mnth.EQ.4).OR.(mnth.EQ.6).OR.(mnth.EQ.9).OR.
     &(mnth.EQ.11)) THEN
            no_days = 30
          ELSE
            no_days = 31
          ENDIF
        ENDIF
      ELSE
        no_days = 30
      ENDIF


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          FUNCTION no_day                             *
*                          ***************                             *
*                                                                      *
*----------------------------------------------------------------------*
      FUNCTION no_day(year,mnth,day,thty_dys)
*----------------------------------------------------------------------*
      INTEGER no_day,year,mnth,day,i,no_days,thty_dys

      no_day = day
      IF (mnth.GT.1) THEN
        DO i=1,mnth-1
          no_day = no_day + no_days(year,i,thty_dys)
        ENDDO
      ENDIF


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          ST2ARR                                      *
*                          ******                                      *
*                                                                      *
* Strips leading zeros and first integer from the string which is its  *
* argument.                                                            *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE ST2ARR(st1,x,ix,no)
      INTEGER x(ix),ix,i,no
      CHARACTER st1*1000

      DO i=1,ix
        CALL STRIPBN(st1,x(i))
      ENDDO

      no = 0
      DO i=1,ix
        IF (x(i).GT.0)  no = no + 1
      ENDDO


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          STRIPBN                                     *
*                          *******                                     *
*                                                                      *
* Strips leading zeros and first integer from the string which is its  *
* argument.                                                            *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE STRIPBN(st1,ans)
      CHARACTER st1*1000
      INTEGER blanks,nb,i,ans,blank

      nb = blanks(st1)
      IF (nb.GT.0) THEN
        DO i=1,1000-nb
          st1(i:i) = st1(i+nb:i+nb)
        ENDDO
        DO i=1000-nb+1,1000
          st1(i:i) = ' '
        ENDDO
      ENDIF

      nb = blank(st1)
      ans = 0
      DO i=1,nb
        ans = ans + 10**(nb-i)*(ichar(st1(i:i))-48)
      ENDDO

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
*                          STRIPBS                                     *
*                          *******                                     *
*                                                                      *
* Strips leading zeros and first string from the string which is its   *
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
*                          SUBROUTUINE STRIPB                          *
*                          ******************                          *
*                                                                      *
* Strips leading zeros of its argument.                                *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE STRIPB(st1)
      CHARACTER st1*1000
      INTEGER blanks,nb,i

      nb = blanks(st1)
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
     &        (ichar(st1(blanks:blanks)).EQ.9)) GOTO 10
      ENDIF
      blanks = blanks - 1


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
      CHARACTER st1*1000,st2*1000
      INTEGER ans,blanks,n_fields

      n_fields = 0
      st2 = st1

10    CONTINUE
      IF (blanks(st2).LT.1000) THEN
        CALL STRIPBN(st2,ans)
        n_fields = n_fields + 1
        GOTO 10
      ENDIF


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE IN2ST                            *
*                          ****************                            *
*                                                                      *
* 'in2st' returns the string value of its integer argument.            *
*                                                                      *
*----------------------------------------------------------------------*
      FUNCTION in2st(int1)
*----------------------------------------------------------------------*
      INTEGER length,loc,int,int1
      PARAMETER(length=4)
      CHARACTER in2st*(length)

      int = int1

      DO loc=1,length
       in2st(loc:loc) = ' '
      ENDDO

      loc = length
30    in2st(loc:loc) = char(ichar('0') + mod(int,10))
      int = int/10
      loc = loc - 1
      IF (int.GT.0)  GOTO 30


      RETURN
      END

*----------------------------------------------------------------------*
      FUNCTION ntags(st1,st2)
*----------------------------------------------------------------------*
      CHARACTER st1(100)*1000,st2*1000,st3*1000
      INTEGER ntags,stcmp,i,j

      i = 0
10    CONTINUE
        i = i + 1
        st3 = st1(i)
        j = stcmp(st2,st3)
      IF ((j.EQ.0).AND.(i.LT.50)) GOTO 10

      IF (j.EQ.1) THEN
        ntags = i
      ELSE
        WRITE(*,*) 'error in tags:',st2(1:10)
        STOP
      ENDIF

      RETURN
      END


*----------------------------------------------------------------------*
      FUNCTION stcmp(st1,st2)
*----------------------------------------------------------------------*
      CHARACTER st1*1000,st2*1000
      INTEGER k,stcmp

      stcmp = 1
      k = 0
10    CONTINUE
        k = k + 1
        IF (((ichar(st1(k:k)).EQ.0).OR.(ichar(st1(k:k)).EQ.32)
     &                             .OR.(ichar(st1(k:k)).EQ.9))   
     &.AND. ((ichar(st2(k:k)).EQ.0).OR.(ichar(st2(k:k)).EQ.32)
     &                             .OR.(ichar(st2(k:k)).EQ.9))) THEN
           stcmp = 1
           GOTO 20
        ENDIF
        IF (st1(k:k).NE.st2(k:k)) THEN
          stcmp = 0
          GOTO 20
        ENDIF
      GOTO 10 
20    CONTINUE


      RETURN
      END

*----------------------------------------------------------------------*
      SUBROUTINE RANDOMV(years,yr0,yrf,idum)
*----------------------------------------------------------------------*
      INTEGER yr0,yrf,years(yrf-yr0+1),idum,ranv,ind(500),i,i1,j

      DO i=1,yrf-yr0+1
        ind(i) = i
      ENDDO

      DO i=1,yrf-yr0+1
        i1 = ranv(1,yrf-yr0+2-i,idum)
        years(i) = ind(i1)+yr0-1
        DO j=i1,yrf-yr0+1-i
          ind(j) = ind(j+1)
        ENDDO
      ENDDO


      RETURN
      END

*----------------------------------------------------------------------*
      FUNCTION ranv(i,j,idum)
*----------------------------------------------------------------------*
      REAL*8 ran1
      INTEGER i,j,ranv,idum

      ranv = int(ran1(idum)*(j-i+1))+i


      RETURN
      END


*----------------------------------------------------------------------*
      FUNCTION ran1(idum)
*----------------------------------------------------------------------*
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL*8 ran1,AM,EPS,RNMX
      PARAMETER(IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     & NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2E-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
 
      IF (idum.le.0.or.iy.eq.0) THEN
        idum = max(-idum,1)
        DO j=NTAB+8,1,-1
          k = idum/IQ
          idum = IA*(idum-k*IQ)-IR*k
          IF (idum.lt.0)  idum = idum + IM
          IF (j.le.NTAB)  iv(j) = idum
        ENDDO
        iy = iv(1)
      ENDIF
      k = idum/IQ
      idum = IA*(idum - k*IQ) - IR*k
      IF (idum.lt.0)  idum = idum + IM
      j = 1 + iy/NDIV
      iy = iv(j)
      iv(j) = idum
      ran1 = min(AM*iy,RNMX)


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          FUNCTION round                              *
*                          **************                              *
*                                                                      *
* 'round' calculates the nearest integer value to the argument.        *
*                                                                      *
*----------------------------------------------------------------------*
      FUNCTION round(x)
*----------------------------------------------------------------------*
      REAL*8 x,rem
      INTEGER round

      rem = abs(x - real(int(x)))

      IF (x.GT.0.0d0) THEN
        IF (rem.GT.0.5d0) THEN
          round = int(x) + 1
        ELSE
          round = int(x)
        ENDIF
      ELSE
        IF (rem.LT.0.5d0) THEN
          round = int(x)
        ELSE
          round = int(x) - 1
        ENDIF
      ENDIF


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE MATCHECK                         *
*                          *******************                         *
*                                                                      *
* MATCHECK computes the Euclidean distance between the vector          *
* [a,b,c,d] and the 'n' rows of 'x'. The minimum distance is stored in *
* 'check' and the row index is stored in 'ind'.                        *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE MATCHECK(x,a,b,c,d,check,ind,n)
*----------------------------------------------------------------------*
      REAL*8 x(n,4),a,b,c,d,check,temp
      INTEGER ind,n,i

      check = 100.0
      DO 10 i=1,n
        temp = max(abs(x(i,1) - a),abs(x(i,2) - b),abs(x(i,3) -
     & c),abs(x(i,4) - d))
        IF (temp.LT.check) THEN
          check = temp
          ind = i
        ENDIF
10    CONTINUE


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE PILE4                            *
*                          ****************                            *
*                                                                      *
* PILE4 puts the row vector [a,b,c,d] at the top of the n*4 matrix 'x' *
* and shifts the other rows down one, dumping the last.                *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE PILE4(x,a,b,c,d,n)
*----------------------------------------------------------------------*
      REAL*8 x(n,4),a,b,c,d
      INTEGER n,i,j

      DO i=1,n-1
        DO j=1,4
          x(n-i+1,j) = x(n-i,j)
        ENDDO
      ENDDO
      x(1,1) = a
      x(1,2) = b
      x(1,3) = c
      x(1,4) = d


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE PILE1                            *
*                          ****************                            *
*                                                                      *
* PILE1 puts the scalar 'a' at the head of the 'n' dimensional vector  *
* 'x' and shifts the other elements down one, loosing the last.        *
*----------------------------------------------------------------------*
      SUBROUTINE PILE1(x,a,n)
*----------------------------------------------------------------------*
      REAL*8 x(n),a
      INTEGER n,i

      DO i=1,n-1
        x(n-i+1) = x(n-i)
      ENDDO
      x(1) = a


      RETURN
      END
