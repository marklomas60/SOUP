*----------------------------------------------------------------------*
      FUNCTION SETARANDOM()
*----------------------------------------------------------------------*
      REAL*8 SETARANDOM
      CHARACTER date*8,time*10,rtime*10
      INTEGER i
*----------------------------------------------------------------------*

      CALL DATE_AND_TIME(date,time)
      DO i=1,10
        rtime(i:i) = time(11-i:11-i)
      ENDDO
      READ(rtime,*) SETARANDOM
      SETARANDOM=SETARANDOM/1000.0

      RETURN
      END


      FUNCTION check_ft_grow(tmp,bbm,bb0,bbf,bbl,xchill,xdschill)
      IMPLICIT NONE
      REAL*8 tmp(12,31),bb0,bbf,bbl,bbsum
      INTEGER bbm,jday,day,mnth,jd,check_ft_grow,xchill,xdschill
      INTEGER chill,dschill,i

      chill = xchill
      dschill = xdschill

*----------------------------------------------------------------------*
* Check for budburst using degree days.                                *
*----------------------------------------------------------------------*
      IF (bbm.GT.0) THEN
      jday = 0
10    jday = jday + 1

*----------------------------------------------------------------------*
* Check for chilling.                                                  *
*----------------------------------------------------------------------*
        IF (chill.EQ.0) THEN
          bbsum = 0.0d0
          DO i=1,20
            mnth = (jday+i-1)/30 + 1
            day = jday + i - (mnth - 1)*30
            mnth = mod(mnth-1,12) + 1
            IF (tmp(mnth,day).LT.-5.0d0)  bbsum = bbsum +
     &max(-10.0d0,tmp(mnth,day)+5.0d0)
          ENDDO
          IF (bbsum.LT.-100) THEN
            chill = 1
            dschill = 1
          ENDIF
        ENDIF
        IF (chill.EQ.1) THEN
          dschill = dschill + 1
        ENDIF
        IF (dschill.GT.260) THEN
          chill = 0
          dschill = 0
        ENDIF

*        print*,'in func ',mnth,day,chill,dschill

        bbsum = 0.0d0
        DO i=1,bbm
          mnth = (jday+i-1)/30 + 1
          day = jday + i - (mnth - 1)*30
          mnth = mod(mnth-1,12) + 1
          IF (tmp(mnth,day).GT.bb0)  bbsum = bbsum + 
     &min(bbf,tmp(mnth,day)-bb0)
        ENDDO
*        print*,bbsum,exp(0.01d0*real(dschill))
        IF (real(bbsum).GE.real(bbl)*exp(-0.01d0*real(dschill))) THEN
          check_ft_grow = 1
          GOTO 20
        ENDIF
        IF ((jday.GE.360+bbm)) THEN
          check_ft_grow = 0
          GOTO 20
        ENDIF
      GOTO 10
20    CONTINUE
      ELSE
        check_ft_grow = 1
      ENDIF


      RETURN
      END

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
*                          FUNCTION triscale2                          *
*                          ******************                          *
*                                                                      *
*----------------------------------------------------------------------*
      FUNCTION triscale2(x,y,x1,x2,x3,y1,y2,y3)
*----------------------------------------------------------------------*
      REAL*8 triscale2,x,y,x1,x2,x3,y1,y2,y3,s1,s2

      IF ((x.LT.x1).OR.(x.GT.x3).OR.(y.LT.y1).OR.(y.GT.y3)) THEN
        triscale2 = 0.0d0
      ELSE
        IF (x.LT.x2) THEN
          s1 = (x - x1)/(x2 - x1)
        ELSE
          s1 = (x - x3)/(x2 - x3)
        ENDIF
        IF (y.LT.y2) THEN
          s2 = (y - y1)/(y2 - y1)
        ELSE
          s2 = (y - y3)/(y2 - y3)
        ENDIF
        triscale2 = s1*s2
      ENDIF


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          FUNCTION triscale                           *
*                          *****************                           *
*                                                                      *
*----------------------------------------------------------------------*
      FUNCTION triscale(x,x1,x2,x3)
*----------------------------------------------------------------------*
      REAL*8 triscale,x,x1,x2,x3

      IF ((x.LT.x1).OR.(x.GT.x3)) THEN
        triscale = 0.0d0
      ELSE
        IF (x.LT.x2) THEN
          triscale = (x - x1)/(x2 - x1)
        ELSE
          triscale = (x - x3)/(x2 - x3)
        ENDIF
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
* Strips leading blanks and first integer from the string which is its *
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
* Strips leading blanks and first integer from the string which is its *
* argument.                                                            *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE STRIPBN(st1,ans)
      CHARACTER st1*1000
      INTEGER blanks,nb,i,ans,blank,check
      LOGICAL minus

      nb = blanks(st1)
      IF (nb.lt.1000) THEN

      minus = .false.
      IF (ichar(st1(nb+1:nb+1)).EQ.45) THEN
        minus = .true.
        nb = nb + 1
      ENDIF
      IF (nb.GT.0) THEN
        DO i=1,1000-nb
          st1(i:i) = st1(i+nb:i+nb)
        ENDDO
        DO i=1000-nb+1,1000
          st1(i:i) = ' '
        ENDDO
      ENDIF

      check = 0
      nb = blank(st1)
      ans = 0
      DO i=1,nb
        IF ((ichar(st1(i:i)).LT.48).OR.(ichar(st1(i:i)).GT.57)) 
     &check = 1 
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

      IF (minus)  ans =-ans

      IF (check.EQ.1) ans = -9999

      ELSE
        ans = -9999
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
*                          SUBROUTUINE STRIPB                          *
*                          ******************                          *
*                                                                      *
* Strips leading blnaks of its argument.                               *
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
*                          FUNCTION last_blank                         *
*                          *******************                         *
*                                                                      *
* 'last_blank' returns the number of characters in st1.                *
*                                                                      *
*----------------------------------------------------------------------*
      FUNCTION last_blank(st1)
*----------------------------------------------------------------------*
      CHARACTER st1*1000
      INTEGER last_blank

      last_blank = 1000
10    CONTINUE
      last_blank = last_blank - 1
      IF (last_blank.GE.1) THEN
        IF ((ichar(st1(last_blank:last_blank)).EQ.32).OR.
     &(ichar(st1(last_blank:last_blank)).EQ.9))  GOTO 10
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
*                          SUBROUTINE IN2ST                            *
*                          ****************                            *
*                                                                      *
* 'in2st' returns the string value of its integer argument.            *
*                                                                      *
*----------------------------------------------------------------------*
      FUNCTION in2st(int1)
*----------------------------------------------------------------------*
      INTEGER length,loc,int,int1
      PARAMETER(length=1000)
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
        ntags = -9999
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
