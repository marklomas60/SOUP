*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE COVER                            *
*                          ****************                            *
*----------------------------------------------------------------------*
      SUBROUTINE COVER(nft,ftmor,ftppm0,cov,bio,bioleaf,nppstore,
     &npp,nps,tmp,prc,slc,rlc,c3old,c4old,firec,ppm,hgt,
     &fireres,fprob,ftprop,ftstmx,stemdp,rootdp,ftsls,ftrls,ilanduse,
     &nat_map,ic0)
*----------------------------------------------------------------------*
      INCLUDE 'array_dims.inc'
      REAL*8 cov(maxage,maxnft),bio(maxage,2,maxnft),bioleaf(maxnft)
      REAL*8 nppstore(maxnft),npp(maxnft),nps(maxnft),tmp(12),prc(12)
      REAL*8 slc(maxnft),rlc(maxnft),firec
      REAL*8 ppm(maxage,maxnft),hgt(maxage,maxnft),ftppm0(maxnft),told
      REAL*8 ftstmx(maxnft),stemdp(360,maxnft),rootdp(360,maxnft)
      REAL*8 fprob,ftprop(maxnft),ngcov,gold,c3old,c4old,fri,norm
      REAL*8 grassrc,ic0(8),sumc
      INTEGER ftsls(maxnft),ftrls(maxnft),nft,ftmor(maxnft),year,i,j
      INTEGER ft,fireres,ilanduse,nat_map(8),age

      IF (ilanduse.eq.2) THEN
*----------------------------------------------------------------------*
* Add up grass coverage.                                               *
*----------------------------------------------------------------------*
      gold = 0.0d0
      DO ft=3,4
        DO age=1,ftmor(nat_map(ft))
          gold = gold + cov(age,nat_map(ft))
        ENDDO
      ENDDO

*----------------------------------------------------------------------*
* Add up tree coverage.                                                *
*----------------------------------------------------------------------*
      told = 0.0d0
      DO ft=5,8
        DO age=1,ftmor(nat_map(ft))
          told = told + cov(age,nat_map(ft))
        ENDDO
      ENDDO

*      IF (gold.GT.0.0d0) THEN
*        c3old = cov(1,2) + cov(2,2)
*        c4old = cov(1,3) + cov(2,3)
*      ENDIF

      ENDIF

*----------------------------------------------------------------------*
* Compute the likelyhood of fire in the current year 'fprob'.          *
* 'find' is the fire index                                             *
*----------------------------------------------------------------------*
      CALL FIRE(prc,tmp,fri,fprob)

*----------------------------------------------------------------------*
* Take off area burnt by fire together with plants past there sell by  *
* date and put this as bare ground ready for new growth 'ngrowth'.     *
* Also shift cover and biomass arrays one to the right.                *
*----------------------------------------------------------------------*
      CALL NEWGROWTH(nft,ftmor,cov,ppm,bio,bioleaf,nppstore,hgt,fprob,
     &npp,nps,ngcov,slc,rlc,fireres,firec)
*----------------------------------------------------------------------*

      IF (ilanduse.eq.2) THEN
*----------------------------------------------------------------------*
* Compute percentages of decid, ever and grass for the new growth      *
* ('ngrowth') this year.                                               *
*----------------------------------------------------------------------*
*      CALL PERC(dof,tmin,ftprop,nft)

*----------------------------------------------------------------------*
* Restrict the rate of trees taking over grassland.                    *
*----------------------------------------------------------------------*
      CALL GRASSREC(nft,ftprop,gold,ngcov,grassrc,nat_map)

*----------------------------------------------------------------------*
* Restrict the rate of bare ground reclimation 0.2, means 20% can be   *
* reclaimed every year and set bare ground cover array.                *
*----------------------------------------------------------------------*
*      CALL BAREREC(ftprop,cov,bpaold,ngcov,barerc)

*----------------------------------------------------------------------*
* Compute c4 c3 grass split.                                           *
*----------------------------------------------------------------------*
*      CALL c3c4(ftprop,c3old,c4old,npp,nps)
      ENDIF
*----------------------------------------------------------------------*
* Set cover arrays to adjust to ftprop as best they can.               *
* ftprop contains the total proportion of that cover, not the          *
* proportion of bare land to assign. Calculate the new ftprop.         *
*----------------------------------------------------------------------*
      norm = 0.0d0
      DO ft=1,nft           
        IF (ftprop(ft).GT.0.0d0) THEN
          ftprop(ft) = ftprop(ft)/100.0d0
          DO j=1,ftmor(ft)
            ftprop(ft)=ftprop(ft) - cov(j,ft)
          ENDDO
          IF (ftprop(ft).LT.0.0d0) ftprop(ft)=0.0d0 !can't remove cov
          norm = norm + ftprop(ft)
        ENDIF
      ENDDO
      DO ft=1,nft
        ftprop(ft) = 100.0d0*ftprop(ft)/norm
      ENDDO

      cov(1,1) = ngcov*ftprop(1)/100.0d0

*----------------------------------------------------------------------*
* Set cover arrays for this years ft proportions, take carbon from     *
* litter to provide nppstore and canopy.                               *
*----------------------------------------------------------------------*
      DO ft=2,nft
        IF (ftprop(ft).GT.0.0d0) THEN

*----------------------------------------------------------------------*
* If no veg exists then set nppstore to initial value and reinitialise.*
*----------------------------------------------------------------------*
          year = 1
10        CONTINUE
            IF (cov(year,ft).GT.0.0d0) THEN
              goto 20
            ELSE
              year = year + 1
              IF (year.EQ.ftmor(ft)+1) THEN
                nppstore(ft) = ftstmx(ft)
                DO i=1,ftsls(ft)
                  stemdp(i,ft) = 0.0d0
                ENDDO
                DO i=1,ftrls(ft)
                  rootdp(i,ft) = 0.0d0
                ENDDO
                GOTO 20
              ENDIF
            ENDIF
          GOTO 10
20        CONTINUE
*----------------------------------------------------------------------*

          cov(1,ft) = ftprop(ft)*ngcov/100.0d0
          ppm(1,ft) = ftppm0(ft)
          hgt(1,ft) = 0.004d0
          slc(ft) = slc(ft) - (nppstore(ft) + bioleaf(ft))*cov(1,ft)
          IF (slc(ft).LT.0.0d0) THEN
            sumc = 0.0d0
            DO i=1,8
               sumc = sumc + ic0(i)
            ENDDO
            DO i=1,8
              ic0(i) = ic0(i)*(1.0d0+slc(ft)/sumc)
            ENDDO
            slc(ft) = 0.0d0
          ENDIF

        ENDIF
      ENDDO


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE GROWTH
*                          *****************                           *
*----------------------------------------------------------------------*
      SUBROUTINE GROWTH(nft,ftmor,ftwd,ftxyl,ftpd,ftgr0,ftgrf,cov,bio,
     &nppstore,npp,lai,nps,npr,evp,slc,rlc,sln,rln,stembio,rootbio,ppm,
     &hgt,leaflit)
*----------------------------------------------------------------------*
      INCLUDE 'array_dims.inc'
      REAL*8 ftwd(maxnft),ftxyl(maxnft),ftpd(maxnft),cov(maxage,maxnft)
      REAL*8 bio(maxage,2,maxnft),npp(maxnft),lai(maxnft),nps(maxnft)
      REAL*8 npr(maxnft),evp(maxnft),nppstore(maxnft),ftgr0(maxnft)
      REAL*8 rootbio,slc(maxnft),rlc(maxnft),sln(maxnft),rln(maxnft)
      REAL*8 stembio,ppm(maxage,maxnft),hgt(maxage,maxnft),ftgrf(maxnft)
      REAL*8 leaflit(maxnft),ftmat(maxnft)
      INTEGER nft,ftmor(maxnft),age,ft,i

*----------------------------------------------------------------------*
* Initialise litter arrays, and add on leaf litter computed in DOLY.   *
*----------------------------------------------------------------------*
      DO ft=1,nft
        slc(ft) = 0.0d0
        rlc(ft) = 0.0d0
        sln(ft) = 0.0d0
        rln(ft) = 0.0d0
      ENDDO
      DO ft=1,nft
        DO age=1,ftmor(ft)
          slc(ft) = slc(ft) + leaflit(ft)*cov(age,ft)
        ENDDO
      ENDDO

*----------------------------------------------------------------------*
* Compute maturity of plant fts.                                       *
*----------------------------------------------------------------------*
*      CALL VEGMAT(nft,ftmor,ftwd,ftxyl,ftpd,evp,lai,npp,
*     &     nps,ftmat,maxnft)

      DO ft=1,nft
        ftmat(ft) = real(ftmor(ft))
      ENDDO

*----------------------------------------------------------------------*
* Thin vegetation where npp is not sufficient to maintain sensible     *
* growth rate.                                                         *
*----------------------------------------------------------------------*
      CALL THIN(nft,ftmor,ftmat,ftwd,ftxyl,ftpd,ftgr0,ftgrf,cov,bio,
     &nppstore,npp,lai,nps,evp,slc,ppm,hgt)

*----------------------------------------------------------------------*
* Add on biomasses for the year.                                       *
*----------------------------------------------------------------------*
      CALL ADDBIO(nft,ftmor,npp,nps,npr,bio)

*----------------------------------------------------------------------*
* Compute litter arrays.                                               *
*----------------------------------------------------------------------*
      CALL MKLIT(nft,ftgr0,ftmor,ftmat,cov,bio,slc,rlc,sln,rln,npp,nps,
     &npr)

*----------------------------------------------------------------------*
*Compute leaf root and stem biomasses.                                 *
*----------------------------------------------------------------------*
      stembio = 0.0d0
      rootbio = 0.0d0
      DO ft=1,nft
        DO i=1,ftmor(ft)
          stembio = stembio + bio(i,1,ft)*cov(i,ft)
          rootbio = rootbio + bio(i,2,ft)*cov(i,ft)
        ENDDO
      ENDDO


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE THIN                             *
*                          ***************                             *
*----------------------------------------------------------------------*
      SUBROUTINE THIN(nft,ftmor,ftmat,ftwd,ftxyl,ftpd,ftgr0,ftgrf,cov,
     &bio,nppstore,npp,lai,nps,evp,slc,ppm,hgt)
*----------------------------------------------------------------------*
      INCLUDE 'array_dims.inc'
      REAL*8 ftwd(maxnft),ftmat(maxnft),ftxyl(maxnft),ftpd(maxnft)
      REAL*8 cov(maxage,maxnft),bio(maxage,2,maxnft),npp(maxnft)
      REAL*8 lai(maxnft),nps(maxnft),evp(maxnft),nppstore(maxnft)
      REAL*8 storelit(maxnft),slc(maxnft),ppm(maxage,maxnft)
      REAL*8 hgt(maxage,maxnft),ftgr0(maxnft),ftgrf(maxnft),pbionew
      REAL*8 pbioold,no,pbio,ftcov(maxnft),covnew(maxage,maxnft)
      REAL*8 ppmnew(maxage,maxnft),shv,lmv,nv,pi,hwv(maxage),emv,totno
      REAL*8 totcov,scale(maxage),oldbio(maxage,2,maxnft),dimold,g0,gf
      REAL*8 gm,grate(maxage),stlit,hgtnew(maxage),dimnew,maxhgt,minhgt
      REAL*8 hc1,hc2,tcov1,tcov2,sum,xxx
      INTEGER nft,ftmor(maxnft),ft,i,year,age

      pi = 3.1415926d0

      hc1 = 0.05d0
      hc2 = 2.0d0

      DO ft=1,nft
        DO year=1,ftmor(ft)
          DO i=1,2
            oldbio(year,i,ft) = bio(year,i,ft)
          ENDDO
        ENDDO
        storelit(ft) = 0.0d0
      ENDDO

*----------------------------------------------------------------------*
* FT loop for thinning and height competition.                         *
*----------------------------------------------------------------------*
      DO ft=1,nft
        IF (ftgr0(ft).GT.0.0d0) THEN

        IF (npp(ft)*nps(ft).GT.0.0d0) THEN

*----------------------------------------------------------------------*
* Height competition.                                                  *
*----------------------------------------------------------------------*
          minhgt = 10000.0d0
          maxhgt =-10000.0d0
          i = 0
          DO year=1,ftmor(ft)
            IF (hgt(year,ft).gt.0) i = i + 1
            scale(year) = hgt(year,ft)
            scale(year) = real(i)
            IF (scale(year).LT.minhgt)  minhgt = scale(year)
            IF (scale(year).GT.maxhgt)  maxhgt = scale(year)
          ENDDO

          IF (maxhgt-minhgt.LT.0.0d0) THEN
            DO year=1,ftmor(ft)
              scale(year) = (scale(year) - minhgt)/(maxhgt - minhgt)
              scale(year) = (scale(year)*hc1 + 1.0d0 - hc1)**hc2
            ENDDO

            tcov1 = 0.0d0
            tcov2 = 0.0d0
            DO year=1,ftmor(ft)
              covnew(year,ft) = cov(year,ft)*scale(year)
              tcov1 = tcov1 + cov(year,ft)
              tcov2 = tcov2 + covnew(year,ft)
            ENDDO

            DO year=1,ftmor(ft)
              IF (covnew(year,ft).GT.0.0d0) THEN
              covnew(year,ft) = covnew(year,ft)*tcov1/tcov2
              DO i=1,2
                bio(year,i,ft) = bio(year,i,ft)*cov(year,ft)/
     &covnew(year,ft)
                oldbio(year,i,ft) = oldbio(year,i,ft)*cov(year,ft)/
     &covnew(year,ft)
              ENDDO
              ppm(year,ft) = ppm(year,ft)*cov(year,ft)/covnew(year,ft)
              cov(year,ft) = covnew(year,ft)
              ENDIF
            ENDDO

          ENDIF
*----------------------------------------------------------------------*

*----------------------------------------------------------------------*
* Height competition.                                                  *
*----------------------------------------------------------------------*

          minhgt = 10000.0d0
          maxhgt =-10000.0d0
          i = 0
          DO year=1,ftmor(ft)
            IF (hgt(year,ft).gt.0) i = i + 1
            scale(year) = 1.0d0 - cov(year,ft)
            IF (scale(year).LT.minhgt)  minhgt = scale(year)
            IF (scale(year).GT.maxhgt)  maxhgt = scale(year)
          ENDDO

          IF (maxhgt-minhgt.GT.0.0d0) THEN
            DO year=1,ftmor(ft)
              scale(year) = (scale(year) - minhgt)/(maxhgt - minhgt)
              scale(year) = (scale(year)*hc1 + 1.0d0 - hc1)**hc2
            ENDDO

            tcov1 = 0.0d0
            tcov2 = 0.0d0
            DO year=1,ftmor(ft)
              covnew(year,ft) = cov(year,ft)*scale(year)
              tcov1 = tcov1 + cov(year,ft)
              tcov2 = tcov2 + covnew(year,ft)
            ENDDO

            DO year=1,ftmor(ft)
              IF (covnew(year,ft).GT.0.0d0) THEN
              covnew(year,ft) = covnew(year,ft)*tcov1/tcov2
              DO i=1,2
                bio(year,i,ft) = bio(year,i,ft)*cov(year,ft)/
     &covnew(year,ft)
                oldbio(year,i,ft) = oldbio(year,i,ft)*cov(year,ft)/
     &covnew(year,ft)
              ENDDO
              ppm(year,ft) = ppm(year,ft)*cov(year,ft)/covnew(year,ft)
              cov(year,ft) = covnew(year,ft)
              ENDIF
            ENDDO

          ENDIF
*----------------------------------------------------------------------*

          g0 = ftgr0(ft)
          gf = ftgrf(ft)
          gm = real(ftmor(ft))/10.0
          DO year=1,ftmor(ft)
            grate(year) = (gf - g0)/gm*real(year - 1) + g0
            IF (real(year).GE.gm)  grate(year) = gf
          ENDDO

          totno = 0.0d0
          totcov = 0.0d0
          DO year=1,ftmor(ft)
            totno = totno + cov(year,ft)*ppm(year,ft)
            totcov = totcov + cov(year,ft)
          ENDDO
          totno = totno/totcov

*----------------------------------------------------------------------*
* Compute cover and ppm to sustain a minimum growth rate, put these    *
* values in covnew and ppmnew.                                         *
*----------------------------------------------------------------------*
          DO year=1,ftmor(ft)

            shv = npp(ft)*nps(ft)/100.0d0/100.0d0
            emv = evp(ft)/3600.0d0/1000.0d0
            lmv = lai(ft)*1.3d0/(lai(ft) + 3.0d0)
            lmv = 1.0d0
            nv = 1.0d0

            IF (ppm(year,ft)*cov(year,ft).GT.0.0d0) THEN
*----------------------------------------------------------------------*
* Calculate the increase in diameter produced by stem NPP 'nps'.       *
*----------------------------------------------------------------------*

              IF (emv.GT.0.0d0) THEN
* hwv = theoretical maximum height (hydrolics)
                hwv(year) = ftpd(ft)*lmv/nv*(ftxyl(ft)*shv/emv/ftwd(ft)/
     &10000.0d0)**0.5d0
              ELSE
                hwv(year) = hgt(year,ft)*0.9d0
              ENDIF

* Calculate new height.
              IF (hwv(year).GT.0.0d0) THEN
                hgtnew(year) = (hwv(year)-hgt(year,ft))/hwv(year)*0.5d0
              ELSE
                hgtnew(year) = 0.0d0
              ENDIF
              IF (hgtnew(year).LT.0.0d0)  hgtnew(year) = 0.0d0
              hgtnew(year) = hgtnew(year) + hgt(year,ft)

* pbio = g/individual
              pbioold = bio(year,1,ft)/ppm(year,ft)
              pbionew = (bio(year,1,ft) + npp(ft)*nps(ft)*(1.0 - 
     &stlit(year,ftmat(ft)))/100.0d0)/ppm(year,ft)

*              print'(''pbio'',3f12.1)',pbionew,pbioold,pbionew-pbioold

* Old diameter
              IF (hgt(year,ft).GT.0.0d0) THEN
*                dimold = 2.0d0*(pbioold/1000000.0d0/hgt(year,ft)/
                dimold = 2.0d0*(pbioold/1000000.0d0/hgtnew(year)/
     &pi/ftwd(ft))**0.5d0
              ELSE
                dimold = 0.0d0
              ENDIF

* New diameter
              IF (hgtnew(year).GT.0.0d0) THEN
                dimnew = 2.0d0*(pbionew/1000000.0d0/hgtnew(year)/
     &pi/ftwd(ft))**0.5d0
              ELSE
                dimnew = 0.0d0
              ENDIF

              IF ((dimnew-dimold)/2.0d0.GT.grate(year)) THEN
*----------------------------------------------------------------------*
* No thinning required.                                                *
*----------------------------------------------------------------------*
                ppmnew(year,ft) = ppm(year,ft)
                covnew(year,ft) = cov(year,ft)
              ELSE
*----------------------------------------------------------------------*
* Thinning required.                                                   *
*----------------------------------------------------------------------*
*               print'(i4,3f8.4)',year,dimnew-dimold,dimnew,dimold
                xxx = (grate(year)+dimold/2.0d0)**2.0d0*1000000.0d0*
     &hgtnew(year)*pi*ftwd(ft)
                ppmnew(year,ft) = (bio(year,1,ft) + npp(ft)*nps(ft)*
     &(1.0d0 - stlit(year,ftmat(ft)))/100.0d0)/xxx
                covnew(year,ft) = cov(year,ft)

              ENDIF

            ELSE

              pbio = 0.0d0
              scale(year) = 0.0d0
              ppmnew(year,ft) = 0.0d0
              covnew(year,ft) = 0.0d0
              hgtnew(year) = 0.0d0

            ENDIF

          ENDDO

*----------------------------------------------------------------------*
* Correct biomass array, and adjust litter for any thinned trees.      *
*----------------------------------------------------------------------*
          DO year=1,ftmor(ft)
            IF ((ppm(year,ft).GT.0.0d0).AND.(cov(year,ft).GT.0.0d0))
     &THEN
                no = ppm(year,ft)*cov(year,ft) - 
     &ppmnew(year,ft)*covnew(year,ft)

              DO i=1,2
                pbio = bio(year,i,ft)/ppm(year,ft)
                slc(ft) = slc(ft) + pbio*no

                IF (cov(year,ft).GT.0.0d0) THEN
                    bio(year,i,ft) = (oldbio(year,i,ft)*cov(year,ft) - 
     &pbio*no)/covnew(year,ft)
                    bio(year,i,ft) = oldbio(year,i,ft)*ppmnew(year,ft)/
     &ppm(year,ft)
                ENDIF
              ENDDO
              storelit(ft) = storelit(ft) + nppstore(ft)/ppm(year,ft)*no
              cov(year,ft) = covnew(year,ft)
              ppm(year,ft) = ppmnew(year,ft)
              hgt(year,ft) = hgtnew(year)
            ELSE
              cov(year,ft) = 0.0d0
              ppm(year,ft) = 0.0d0
              hgt(year,ft) = 0.0d0
            ENDIF

          ENDDO
*----------------------------------------------------------------------*

        ELSE
          DO year=1,ftmor(ft)
            ppm(year,ft) = 0.0d0
            hgt(year,ft) = 0.0d0
          ENDDO
        ENDIF

      ENDIF
*----------------------------------------------------------------------*
* End of ft loop.
*----------------------------------------------------------------------*
      ENDDO

*----------------------------------------------------------------------*
* Correct nppstore to account for thinning.                            *
*----------------------------------------------------------------------*
      DO ft=1,nft
        IF (ftgr0(ft).GT.0.0d0) THEN

        IF (storelit(ft).GT.0.0d0) THEN
          ftcov(ft) = 0.0d0
          DO age=1,ftmor(ft)
            ftcov(ft) = ftcov(ft) + cov(age,ft)
          ENDDO
          IF (ftcov(ft).GT.0.0d0) THEN
            slc(ft) = slc(ft) +  storelit(ft)
            nppstore(ft) = (nppstore(ft)*ftcov(ft) - 
     &storelit(ft))/ftcov(ft)
          ENDIF
        ENDIF

        ENDIF
      ENDDO
*----------------------------------------------------------------------*


      RETURN
      END

*----------------------------------------------------------------------*
*                             FUNCTION find                            *
*                             *************                            *
*----------------------------------------------------------------------*
      FUNCTION find(tmp,prc)
*----------------------------------------------------------------------*
      REAL*8 find,tmp(12),prc(12)
      REAL*8 tmplim,prclim
      INTEGER i

      tmplim = -5.0d0
      prclim = 50.0d0

      find = 0.0d0
      DO i=1,12
        IF (tmp(i).GT.tmplim) THEN
          IF (prc(i).LT.prclim) THEN
            find = find + prc(i)/(12.0d0*prclim)
          ELSE
            find = find + 1/12.0d0
          ENDIF
        ELSE
          find = find + 1/12.0d0
        ENDIF
      ENDDO


      RETURN
      END

*----------------------------------------------------------------------*
*                             SUBROUTINE c3c4                          *
*                             ***************                          *
*----------------------------------------------------------------------*
      SUBROUTINE c3c42(ftprop,npp,range)
*----------------------------------------------------------------------*
      INCLUDE 'array_dims.inc'
      REAL*8 ftprop(maxnft),npp(maxnft),range
      REAL*8 grass,nd

      grass = ftprop(2)
      nd = npp(2) - npp(3)
      ftprop(2) = grass*nd/(2.0d0*range) + grass/2.0d0
      IF (ftprop(2).GT.grass)  ftprop(2) = grass
      IF (ftprop(2).LT.0.0d0)  ftprop(2) = 0.0d0
      ftprop(3) = grass - ftprop(2)


      RETURN
      END

*----------------------------------------------------------------------*
*                             SUBROUTINE c3c4                          *
*                             ***************                          *
*----------------------------------------------------------------------*
      SUBROUTINE c3c4(ftprop,c3old,c4old,npp,nps)
*----------------------------------------------------------------------*
      INCLUDE 'array_dims.inc'
      REAL*8 ftprop(maxnft),npp(maxnft),nps(maxnft),c3old,c4old
      REAL*8 grass,c3p,c4p,adj

      IF (c3old+c4old.GT.0.0d0) THEN
        c3p = c3old/(c3old + c4old)
        c4p = c4old/(c3old + c4old)
      ELSE
        c3p = 0.5d0
        c4p = 0.5d0
      ENDIF

      IF (npp(2)*nps(2)+npp(3)*nps(3).GT.0.0d0) THEN
        adj = npp(2)*nps(2)/(npp(2)*nps(2) + npp(3)*nps(3)) - 0.5d0
      ELSE
        adj = 0.0d0
      ENDIF

      c3p = c3p + adj
      c4p = c4p - adj

      IF (c3p.LT.0.0d0) THEN
        c3p = 0.0d0
        c4p = 1.0d0
      ENDIF

      IF (c4p.LT.0.0d0) THEN
        c4p = 0.0d0
        c3p = 1.0d0
      ENDIF

      grass = ftprop(2) + ftprop(3)

      ftprop(2) = grass*c3p
      ftprop(3) = grass*c4p

      RETURN
      END

*----------------------------------------------------------------------*
*                             SUBROUTINE GRASSREC                      *
*                             *******************                      *
*----------------------------------------------------------------------*
      SUBROUTINE GRASSREC(nft,ftprop,gold,ngcov,x,nat_map)
*----------------------------------------------------------------------*
      INCLUDE 'array_dims.inc'
      REAL*8 ftprop(maxnft),gold,ngcov,x,ntcov,ftt,ftpropo(maxnft)
      INTEGER nft,ft,nat_map(8)

      ntcov = 0.0d0
      ftt = 0.0d0
      DO ft=5,8
        ntcov = ntcov + ftprop(nat_map(ft))*ngcov/100.0d0
        ftt = ftt + ftprop(nat_map(ft))
      ENDDO

      IF (ntcov.GT.gold*x) THEN
        DO ft=4,nft
          ftpropo(ft) = ftprop(ft)
          ftprop(ft) = 100.0d0*gold*x*ftprop(ft)/(ftt*ngcov)
          ftprop(2) = ftprop(2) + ftpropo(ft) - ftprop(ft)
        ENDDO

        ntcov = 0.0d0
        DO ft=4,nft
          ntcov = ntcov + ftprop(ft)*ngcov/100.0d0
        ENDDO
        IF (abs(ntcov-gold*x).GT.0.000001d0) WRITE(11,
     &'(''Treerec subroutine error'')')
      ENDIF


      RETURN
      END

*----------------------------------------------------------------------*
*                             SUBROUTINE BAREREC                       *
*                             ******************                       *
*----------------------------------------------------------------------*
      SUBROUTINE BAREREC(ftprop,cov,bpaold,ngcov,x)
*----------------------------------------------------------------------*
      INCLUDE 'array_dims.inc'
      REAL*8 ftprop(maxnft),bpaold,ngcov,x,cov(maxage,maxnft)
      REAL*8 nbp,cgcov

* Current growth 'cgcov'
      cgcov = 1.0d0 - ngcov

      nbp = ngcov*ftprop(1)/100.0d0

      IF (nbp.LT.bpaold*(1.0d0-x)) THEN
	ngcov = 1.0d0 - bpaold*(1.0d0 - x) - cgcov
	cov(1,1) = bpaold*(1.0d0 - x)
      ELSE
	cov(1,1) = ngcov*ftprop(1)/100.0d0
	ngcov = ngcov*(1.0d0 - ftprop(1)/100.0d0)
      ENDIF


      RETURN
      END

*----------------------------------------------------------------------*
*                            SUBROUTINE ADDBIO                         *
*                            *****************                         *
* Add on biomass to bio array.                                         *
*----------------------------------------------------------------------*
      SUBROUTINE ADDBIO(nft,ftmor,npp,nps,npr,bio)
*----------------------------------------------------------------------*
      INCLUDE 'array_dims.inc'
      REAL*8 npp(maxnft),nps(maxnft),npr(maxnft),bio(maxage,2,maxnft)
      REAL*8 npps,nppr
      INTEGER nft,ftmor(maxnft),ft,age

      DO ft=1,nft
        npps = nps(ft)*npp(ft)/100.0d0
        nppr = npr(ft)*npp(ft)/100.0d0
        DO age=1,ftmor(ft)
          bio(age,1,ft) = bio(age,1,ft) + npps
          bio(age,2,ft) = bio(age,2,ft) + nppr
        ENDDO
      ENDDO


      RETURN
      END

*----------------------------------------------------------------------*
*                            SUBROUTINE PERC                           *
*                            ***************                           *
*----------------------------------------------------------------------*
      SUBROUTINE PERC(fttags,dof,tmin,ftprop,nft)
*----------------------------------------------------------------------*
      INCLUDE 'array_dims.inc'
      REAL*8 dof(maxnft),tmin,ftprop(maxnft),dec1
      REAL*8 evbl,evnl,dec,dcnl,drdc,triscale2,sum
      REAL*8 t1eb,t2eb,t3eb,d1eb,d2eb,d3eb
      REAL*8 t1en,t2en,t3en,d1en,d2en,d3en
      REAL*8 t1d,t2d,t3d,d1d,d2d,d3d
      REAL*8 t1dd,t2dd,t3dd,d1dd,d2dd,d3dd
      REAL*8 t1dn,t2dn,t3dn,d1dn,d2dn,d3dn
      REAL*8 t1d1,t2d1,t3d1,d1d1,d2d1,d3d1,ndof
      INTEGER i,ft,ntags,nft
      CHARACTER fttags(100)*1000,st1*1000


      t1eb = 100.0d0   ! 0
      t2eb = -10.0d0   ! 1
      t3eb = -25.0d0   ! 0
      d1eb =  -1.0d0   ! 0
      d2eb =   0.0d0   ! 1
      d3eb = 150.0d0   ! 0

      t1en = -20.0d0   ! 0
      t2en = -40.0d0   ! 1
      t3en = -60.0d0   ! 0
      d1en = 100.0d0   ! 0
      d2en = 180.0d0   ! 1
      d3en = 260.0d0   ! 0

      t1d =  10.0d0    ! 0
      t2d = -15.0d0    ! 1  broadleaf
      t3d = -50.0d0    ! 0
      d1d =-999.0d0    ! 0
      d2d = 100.0d0    ! 1
      d3d = 140.0d0    ! 0

      t1d1 = 100.0d0   ! 0
      t2d1 =  20.0d0   ! 1  broadleaf
      t3d1 =  10.0d0   ! 0
      d1d1 = 150.0d0   ! 0
      d2d1 = 200.0d0   ! 1
      d3d1 = 400.0d0   ! 0

      t1dd =  12.0d0   ! 0
      t2dd = -15.0d0   ! 1  broadleaf
      t3dd = -30.0d0   ! 0
      d1dd =  75.0d0   ! 0
      d2dd = 200.0d0   ! 1
      d3dd = 400.0d0   ! 0

      t1dn = -60.0d0   ! 0
      t2dn = -70.0d0   ! 1
      t3dn = -80.0d0   ! 0
      d1dn = 180.0d0   ! 0
      d2dn = 240.0d0   ! 1
      d3dn = 300.0d0   ! 0

      nft = 1
      st1='Dc_Nl'
*      ndof = dof(ntags(fttags,st1,nft))
      ndof = 100.0d0
      evbl = triscale2(tmin,ndof,t1eb,t2eb,t3eb,d1eb,d2eb,d3eb)
      evnl = triscale2(tmin,ndof,t1en,t2en,t3en,d1en,d2en,d3en)
      dec  = triscale2(tmin,ndof,t1d ,t2d ,t3d ,d1d ,d2d ,d3d )
      drdc = triscale2(tmin,ndof,t1dd,t2dd,t3dd,d1dd,d2dd,d3dd)
      dcnl = triscale2(tmin,ndof,t1dn,t2dn,t3dn,d1dn,d2dn,d3dn)
      dec1 = triscale2(tmin,ndof,t1d1,t2d1,t3d1,d1d1,d2d1,d3d1)

*      IF (tex-1.LT.0.01d0)  dcnl = 0.0d0

      sum = evbl + evnl + dec + drdc + dcnl + dec1
      IF (sum.GT.0.001d0) THEN
        ftprop(1) = 0.0d0
        ftprop(2) = 0.0d0
        ftprop(3) = 0.0d0
        st1='Ev_Bl'
        ftprop(ntags(fttags,st1,nft)) = 100.0d0*evbl/sum
        st1='Ev_Nl'
        ftprop(ntags(fttags,st1,nft)) = 100.0d0*evnl/sum
        st1='Dc_Bl'
        ftprop(ntags(fttags,st1,nft)) = 100.0d0*(dec + drdc + dec1)/sum
        st1='Dc_Nl'
        ftprop(ntags(fttags,st1,nft)) = 100.0d0*dcnl/sum
      ELSE
        DO i=1,nft
          ftprop(i) = 0.0d0
        ENDDO
        ftprop(2) = 100.0d0
      ENDIF

      ftprop(1) = 100.0d0
      DO ft=2,nft
	  ftprop(1) = ftprop(1) - ftprop(ft)
      ENDDO


      RETURN
      END

*----------------------------------------------------------------------*
*                            SUBROUTINE MKLIT                          *
*                            ****************                          *
* Make litter.                                                         *
*----------------------------------------------------------------------*
      SUBROUTINE MKLIT(nft,ftgr0,ftmor,ftmat,cov,bio,slc,rlc,sln,rln,
     &npp,nps,npr)
*----------------------------------------------------------------------*
      INCLUDE 'array_dims.inc'
      REAL*8 ftmat(maxnft),cov(maxage,maxnft),bio(maxage,2,maxnft)
      REAL*8 slc(maxnft),rlc(maxnft),sln(maxnft),ftgr0(maxnft)
      REAL*8 rln(maxnft),npp(maxnft),nps(maxnft),npr(maxnft),npps,nppr
      REAL*8 sl,rl,stlit
      INTEGER nft,ftmor(maxnft),ft,age

      DO ft=1,nft
        IF (ftgr0(ft).GT.0.0d0) THEN
          npps = nps(ft)*npp(ft)/100.0d0
          nppr = npr(ft)*npp(ft)/100.0d0
          DO age=2,ftmor(ft)
            sl = stlit(age,ftmat(ft))
            rl = 0.5d0
            slc(ft) = slc(ft) + sl*npps*cov(age,ft)
            rlc(ft) = rlc(ft) + rl*nppr*cov(age,ft)
            bio(age,1,ft) = bio(age,1,ft) - sl*npps
            bio(age,2,ft) = bio(age,2,ft) - rl*nppr
          ENDDO
        ENDIF
      ENDDO

      DO ft=1,nft
        sln(ft) = 0.0d0
        rln(ft) = 0.0d0
      ENDDO


      RETURN
      END

*----------------------------------------------------------------------*
      FUNCTION stlit(age,mat)
*----------------------------------------------------------------------*
      REAL*8 stlit,mat,temp
      INTEGER age

      stlit = 0.9d0*real(age)/mat + 0.1d0
      IF (stlit.GT.1.0)  stlit = 1.0d0

      IF (age.LT.mat) then
        stlit = 0.0d0
      ELSE
        stlit = 1.0d0
      ENDIF

      temp = age/mat
      IF (temp.GT.1.0d0) temp = 1.0d0

      stlit = temp**0.5d0
      stlit = stlit*0.1d0
*      stlit = 0.0d0

      RETURN
      END


*----------------------------------------------------------------------*
*                            SUBROUTINE VEGMAT                         *
*                            *****************                         *
* Compute the years to reach maturity for each of the fts.             *
*----------------------------------------------------------------------*
      SUBROUTINE VEGMAT(nft,ftmor,ftwd,ftxyl,ftpd,evp,lai,npp,nps,ftmat)
*----------------------------------------------------------------------*
      INCLUDE 'array_dims.inc'
      REAL*8 ftwd(maxnft),ftxyl(maxnft),ftpd(maxnft),lai(maxnft)
      REAL*8 npp(maxnft),nps(maxnft),evp(maxnft),ftmat(maxnft),emxv,wd
      REAL*8 xyl,pd,lv,fs,shv,pi,lmvt,nvt,hwvt,dvt,massvt,minvt,ppvt
      INTEGER ft,nft,ftmor(maxnft)

      pi = 3.14159d0
      DO ft=2,nft
	  IF (npp(ft).GT.1.0e-6) THEN
          emxv = evp(ft)/3600.0d0/1000.0d0
          lv = lai(ft)
          fs = nps(ft)/100.0d0
          shv = fs*npp(ft)/100.0d0

          lmvt = lv*1.3d0/(lv + 3.0d0)
          nvt = 1.0d0 + 26.0d0*exp(-0.9d0*lv)
          nvt = 1.0d0

          wd = ftwd(ft)
          xyl = ftxyl(ft)
          pd = ftpd(ft)

          hwvt = pd*lmvt*sqrt(xyl*shv/(emxv*wd*10000.0d0))
     &/nvt
          dvt = 0.0028d0*hwvt**1.5d0
          massvt = pi*(dvt/2.0d0)**2.0d0*hwvt*wd*nvt
          minvt = pi*((dvt + 0.001d0/nvt)/2.0d0)**2.0d0*hwvt*wd*nvt
     &- massvt
          ppvt = shv/minvt
          ftmat(ft) = massvt*ppvt/shv
          ftmat(ft) = 1.0d0
	  ELSE
	    ftmat(ft) = 0.0d0
        ENDIF

        IF (ftmor(ft).LT.ftmat(ft)) ftmat(ft) = ftmor(ft)
      ENDDO


      RETURN
      END

*----------------------------------------------------------------------*
*                            SUBROUTINE FIRE                           *
*                            ***************                           *
* Compute the fire return interval 'fri' and the probability of a fire *
* in the current year 'fprob'.                                         *
*----------------------------------------------------------------------*
      SUBROUTINE FIRE(prc,tmp,fri,fprob)
*----------------------------------------------------------------------*
      REAL*8 fri,fprob,tmp(12),prc(12),prct(12),prco,lim1,lim2,maxfri
      REAL*8 pow,tlim1,tlim2,totp,indexx,tadj,weight
      INTEGER k,no,ind,minx,i
      INCLUDE 'param.inc'

*----------------------------------------------------------------------*
* Fire model parameters.                                               *
*----------------------------------------------------------------------*
      no = 3
      lim1 = 150.0d0
      lim2 = 50.0d0
      tlim1 =  0.0d0
      tlim2 = -5.0d0
      weight = 0.5d0
      maxfri = 800.0d0
      pow = 3.0d0

*----------------------------------------------------------------------*
* Adjust for temperature.                                             *
*----------------------------------------------------------------------*
      DO i=1,12
        IF (tmp(i).GT.tlim1) THEN
          tadj = 0.0d0
        ELSEIF (tmp(i).LT.tlim2) THEN
          tadj = lim1
        ELSE
          tadj = (tlim1 - tmp(i))/(tlim1 - tlim2)*lim1
        ENDIF
        prct(i) = min(lim1,prc(i) + tadj)
      ENDDO

*----------------------------------------------------------------------*
* Calculate yearly componet of index.                                  *
*----------------------------------------------------------------------*
      totp = 0.0d0
      DO k=1,12
        totp = totp + prct(k)/lim1/12.0d0
      ENDDO

*----------------------------------------------------------------------*
* Calculate monthly component of index.                                *
*----------------------------------------------------------------------*
      prco = 0.0d0
      DO k=1,no
        ind = minx(prct)
        prco = prco + min(lim2,prct(ind))/real(no)/lim2
        prct(ind) = lim2
      ENDDO

*----------------------------------------------------------------------*
* Compute fire return interval, and convert to probability.            *
*----------------------------------------------------------------------*
      indexx = weight*(prco) + (1.0d0 - weight)*totp

      fri = indexx**pow*maxfri

      IF (fri.LT.2.0d0)  fri = 2.0d0
      fprob = (1.0d0 - exp(-1.0d0/fri))*p_fprob


      RETURN
      END

*----------------------------------------------------------------------*
*                            SUBROUTINE NEWGROWTH                      *
*                            ********************                      *
* Compute newgrowth and alter cover array accordingly.                 *
*----------------------------------------------------------------------*
      SUBROUTINE NEWGROWTH(nft,ftmor,cov,ppm,bio,bioleaf,nppstore,hgt,
     &fprob,npp,nps,ngcov,slc,rlc,fireres,firec)
*----------------------------------------------------------------------*
      INCLUDE 'array_dims.inc'
      REAL*8 bio(maxage,2,maxnft),cov(maxage,maxnft),ppm(maxage,maxnft)
      REAL*8 hgt(maxage,maxnft),fprob,npp(maxnft),nppstore(maxnft)
      REAL*8 nps(maxnft),ngcov,slc(maxnft),rlc(maxnft),bioleaf(maxnft)
      REAL*8 tmor,tmor0,npp0,firec,xfprob
      INTEGER nft,ftmor(maxnft),ft,age,fireres

      xfprob = fprob

      firec = 0.0d0
*----------------------------------------------------------------------*
* Take away veg that has died of old age ie > than ftmor(ft), and      *
* shift cover array on one year.                                       *
*----------------------------------------------------------------------*
      ngcov = 0.0d0
      DO ft=1,nft
        ngcov = ngcov + cov(ftmor(ft),ft)
        slc(ft) = slc(ft)  + (bio(ftmor(ft),1,ft) + bioleaf(ft) +
     &nppstore(ft))*cov(ftmor(ft),ft)
        rlc(ft) = rlc(ft)  + bio(ftmor(ft),2,ft)*cov(ftmor(ft),ft)
      ENDDO
      CALL SHIFT(ftmor,cov,ppm,bio,hgt,1,nft)

*----------------------------------------------------------------------*
* Take away veg that is burnt or has died through a hard year ie small *
* LAI.                                                                 *
*----------------------------------------------------------------------*
      DO ft=2,nft

*----------------------------------------------------------------------*
* 'tmor' is the mortality rate of the forrest based on 'npp'.          *
*----------------------------------------------------------------------*
        npp0 = 0.2d0
        tmor0 = 6.0d0
        IF (npp(ft)/100.0d0.LT.0.1d0) THEN
          tmor = 1.0d0
        ELSEIF (npp(ft)/100.0d0.LT.3.0d0) THEN
          tmor = (0.04d0 - npp0)/3.0d0*(npp(ft)/100.0d0) + npp0
        ELSEIF (npp(ft)/100.0d0.LT.tmor0) THEN
          tmor = 0.04d0/(3.0d0 - tmor0)*(npp(ft)/100.d0) + 
     &3.0d0*0.04d0/(tmor0 - 3.0d0) + 0.04d0
        ELSE
          tmor = 0.0d0
        ENDIF
        IF (tmor.LT.0.01d0)  tmor = 0.01d0

*        IF (npp(ft)*nps(ft)/100.0d0.GT.10.0d0) THEN
          tmor = 0.002d0
*        ELSE
*          tmor = 1.0d0 - npp(ft)*nps(ft)/1000.0d0
*          if (ft.eq.6) print*,'hello ',tmor
*        ENDIF

*----------------------------------------------------------------------*
* kill off trees not able to sustain any growth rate                   *
*----------------------------------------------------------------------*
        IF ((.NOT.(npp(ft)*nps(ft).GT.0.0d0)).OR.
     &(.NOT.(nppstore(ft)).GT.0.0d0)) THEN
          tmor = 1.0d0
*          WRITE(*,*) 'Death ft npp nps store: ',ft,npp(ft),nps(ft),
*     &nppstore(ft)
        ENDIF
*----------------------------------------------------------------------*

        DO age=2,ftmor(ft)
          IF (age.LT.fireres) THEN
            fprob = xfprob
          ELSE
            fprob = 0.0d0
          ENDIF
          IF (fireres.LT.0) fprob = real(-fireres)/1000.0d0

          ngcov = ngcov + cov(age,ft)*(fprob - fprob*tmor + tmor)
          slc(ft) = slc(ft) + (bio(age,1,ft) + bioleaf(ft)+ 
     &nppstore(ft))*(tmor - 0.2d0*fprob*tmor + 0.2d0*fprob)*
     &cov(age,ft)
          rlc(ft) = rlc(ft) + bio(age,2,ft)*(tmor - fprob*tmor + fprob)
     &*cov(age,ft)
          firec = firec + (bio(age,1,ft) + bioleaf(ft) + nppstore(ft))*
     &(0.8d0*fprob - 0.8d0*fprob*tmor)*cov(age,ft)
          cov(age,ft) = cov(age,ft)*(1.0d0 - fprob)*(1.0d0 - tmor)
        ENDDO

      ENDDO

      fprob = xfprob


      RETURN
      END

*----------------------------------------------------------------------*
*                            SUBROUTINE SHIFT                          *
*                            ****************                          *
* Shift shifts elements of an array to the right and leavs the first   *
* element equal to zero.                                               *
*----------------------------------------------------------------------*
      SUBROUTINE SHIFT(ftmor,cov,ppm,bio,hgt,ft1,ft2)
*----------------------------------------------------------------------*
      INCLUDE 'array_dims.inc'
      REAL*8 cov(maxage,maxnft),ppm(maxage,maxnft),bio(maxage,2,maxnft)
      REAL*8 hgt(maxage,maxnft)
      INTEGER ftmor(maxnft),ft1,ft2,ft,age,i

      DO ft=ft1,ft2
        DO age=2,ftmor(ft)
          cov(ftmor(ft)-age+2,ft) = cov(ftmor(ft)-age+1,ft)
          ppm(ftmor(ft)-age+2,ft) = ppm(ftmor(ft)-age+1,ft)
          hgt(ftmor(ft)-age+2,ft) = hgt(ftmor(ft)-age+1,ft)
          DO i=1,2
            bio(ftmor(ft)-age+2,i,ft) = bio(ftmor(ft)-age+1,i,ft)
          ENDDO
        ENDDO
        cov(1,ft) = 0.0d0
        ppm(1,ft) = 0.0d0
        hgt(1,ft) = 0.0d0
        bio(1,1,ft) = 0.0d0
        bio(1,2,ft) = 0.0d0
      ENDDO


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          NATURAL_VEG                                 *
*                          ***********                                 *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE NATURAL_VEG(tmp,prc,ftprop,nat_map)
*----------------------------------------------------------------------*
      INCLUDE 'array_dims.inc'
      REAL*8 tmp(12,31),prc(12,31),ftprop(maxnft),min_tmp,av_mnth(12)
      REAL*8 x(6,6),evap,tran,bucket,dep,prc_ind
      REAL*8 triscale,sum,prop(6)
      INTEGER day,mnth,ft,nat_map(8),vdaysoff(5),i,year,daysoff,minx
      LOGICAL leaves
      DATA (x(1,i),i=1,6)/ -15.0, -10.0, 100.0,   -1.0,   0.0, 150.0/!EB
      DATA (x(2,i),i=1,6)/ -70.0, -45.0, -20.0,  -10.0,   0.0, 350.0/!EN
      DATA (x(3,i),i=1,6)/ -50.0, -15.0,  10.0, -999.0, 100.0, 140.0/!DB
      DATA (x(4,i),i=1,6)/  10.0,  20.0, 100.0,  150.0, 200.0, 400.0/!DB
      DATA (x(5,i),i=1,6)/ -30.0, -15.0,  12.0,  75.0,  200.0, 400.0/!DB
      DATA (x(6,i),i=1,6)/ -80.0, -70.0, -60.0, 200.0,  280.0, 300.0/!DN

*----------------------------------------------------------------------*
* Find average monthly min, and convert to an absolute min.            *
*----------------------------------------------------------------------*
      min_tmp = 1000.0d0
      DO mnth=1,12
        av_mnth(mnth) = 0.0d0
        DO day=1,30
          IF (min_tmp.GT.tmp(mnth,day)) min_tmp = tmp(mnth,day)
          av_mnth(mnth) = av_mnth(mnth) + tmp(mnth,day)     
        ENDDO
        av_mnth(mnth) = av_mnth(mnth)/30.0d0     
      ENDDO
      min_tmp = av_mnth(minx(av_mnth))*1.29772d0 - 19.5362d0

*----------------------------------------------------------------------*
* Find days off index. Based on an imaginary bucket.                   *
*----------------------------------------------------------------------*
      dep = 100.0d0
      evap = 0.5d0
      tran = 0.5d0

      DO i=1,5
        vdaysoff(i) = 0
      ENDDO
      leaves = .false.

      bucket = 0.0d0

      DO year=1,10
        daysoff = 0
        DO mnth=1,12
          DO day=1,30
            bucket = bucket + prc(mnth,day)
            IF (leaves) THEN
              bucket = bucket - evap - tran
            ELSE
              bucket = bucket - evap
            ENDIF
            IF (bucket.LT.0.0) THEN
              leaves = .false.
              bucket = 0.0
            ENDIF 
            IF (bucket.gt.dep/2.0) leaves = .true.
            IF (.NOT.(leaves)) daysoff = daysoff + 1
            IF (bucket.GT.dep) bucket = dep
          ENDDO
        ENDDO
        DO i=1,4
          vdaysoff(6-i) = vdaysoff(5-i)
        ENDDO
        vdaysoff(1) = daysoff
      ENDDO
      prc_ind = real(vdaysoff(1)+vdaysoff(2)+vdaysoff(3)+vdaysoff(4))
     &/4.0d0

*----------------------------------------------------------------------*
* Compute fractions of natural vegetation appropriate to seed new      *
* ground.                                                              *
*----------------------------------------------------------------------*
      DO i=1,6
        prop(i) =  triscale(min_tmp,x(i,1),x(i,2),x(i,3))*
     &             triscale(prc_ind,x(i,4),x(i,5),x(i,6))
      ENDDO

*----------------------------------------------------------------------*
* Map the natural functional types to the list given in the input file.*
*----------------------------------------------------------------------*
      ftprop(nat_map(1)) = 0.0d0
      ftprop(nat_map(2)) = 0.0d0
      ftprop(nat_map(3)) = 0.0d0
      ftprop(nat_map(4)) = 0.0d0
      ftprop(nat_map(5)) = prop(1)
      ftprop(nat_map(6)) = prop(2)
      ftprop(nat_map(7)) = prop(3) + prop(4) + prop(5)
      ftprop(nat_map(8)) = prop(6)

*----------------------------------------------------------------------*
* Normalise the fractions.                                             *
*----------------------------------------------------------------------*
      sum = 0.0d0
      DO ft=1,8
        sum = sum + ftprop(nat_map(ft))
      ENDDO
      IF (sum.GT.0.0d0) THEN
        DO ft=1,8
          ftprop(nat_map(ft)) = 100.0d0*ftprop(nat_map(ft))/sum
        ENDDO
      ELSE
        ftprop(nat_map(3)) = 100.0d0
      ENDIF

*----------------------------------------------------------------------*

*----------------------------------------------------------------------*
*      DO k=-80,100
*        DO j=0,360
*        min_tmp =  (real(k) + 19.5632d0)/1.29772d0
*        prc_ind = real(j)
*      DO i=1,6
*        ftprop(i) =  triscale(min_tmp,x(i,1),x(i,2),x(i,3))*
*     &               triscale(prc_ind,x(i,4),x(i,5),x(i,6))
*      ENDDO
*      sum=ftprop(1)+ftprop(2)+ftprop(3)+ftprop(4)+ftprop(5)+ftprop(6)
*        IF (sum.LT.1.0e-6) print'(2i4)',k,j
*        ENDDO
*      ENDDO
*      STOP
*----------------------------------------------------------------------*


      RETURN
      END
