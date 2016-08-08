      real*8 sand(360),silt(360),clay(360),bulk(360),satw(360),xxx(10)
      real*8 lat,lon,dep,orgc(360),wilt(360),field(360),wp,fc,swc
      real*8 sdgvm_wilt(360),sdgvm_field(360),sdgvm_satw(360),topsl
      REAL*8 nupc,awl(4),kd,kx,nci(4),nfix,adp(4),sfc(4),sw(4),sswc(4)

      integer i,j,ii(360)

      open(11,file='islscp2_soils_1deg/soil_bulk_dens0-150_1d.asc')
      open(12,file='islscp2_soils_1deg/soil_sand_perc0-150_1d.asc')
      open(13,file='islscp2_soils_1deg/soil_silt_perc0-150_1d.asc')
      open(14,file='islscp2_soils_1deg/soil_clay_perc0-150_1d.asc')
      open(15,file='islscp2_soils_1deg/soil_org_carb_perc0-150_1d.asc')
      open(16,file='islscp2_soils_1deg/soil_wilting_point0-150_1d.asc')
      open(17,file='islscp2_soils_1deg/soil_field_cap0-150_1d.asc')
      open(18,file='islscp2_soils_1deg/soil_sat_wat_cont0-150_1d.asc')
      open(19,file='islscp2_soils_1deg/soil_pawc0-150_1d.asc')
      open(20,file='islscp2_soils_1deg/soil_res_wat_cont0-150_1d.asc')
      open(21,file='islscp2_soils_1deg/soil_ksat0-150_1d.asc')
      open(22,file='islscp2_soils_1deg/soil_org_carb_dens0-150_1d.asc')
      open(23,file='islscp2_soils_1deg/soil_nitrogen_dens0-150_1d.asc')
      open(24,file='islscp2_soils_1deg/soil_therm_cap0_0-150_1d.asc')
      open(25,file='islscp2_soils_1deg/soil_therm_cap10_0-150_1d.asc')
      open(26,file='islscp2_soils_1deg/soil_therm_cap50_0-150_1d.asc')
      open(27,file='islscp2_soils_1deg/soil_therm_cap100_0-150_1d.asc')
      open(28,file='islscp2_soils_1deg/soil_texture0-150_1d.asc')

      open(29,file='../../../islscp/1deg_new/data.dat',status='old')

      open(31,file='../0-150/data.dat')
      open(32,file='../0-150/data_sdgvm.dat')

      dep =-99.0d0

      DO i=1,180
        read(11,*) bulk
        read(12,*) sand
        read(13,*) silt
        read(14,*) clay
        read(15,*) orgc
        read(16,*) wilt
        read(17,*) field
        read(18,*) satw
        DO j=1,360
          IF (wilt(j).GT.0.0d0) THEN
            wilt(j) = wilt(j)/1500.0d0
            field(j) = field(j)/1500.0d0
          ENDIF
        ENDDO

        lat =  90.5d0 - real(i)

        DO j=1,360
          lon = -180.5d0 + real(j)
          read(29,*) xxx
          write(31,'(f7.3,f9.3,14f10.4)') lat,lon,sand(j),silt(j),
     &bulk(j),orgc(j),wilt(j),field(j),satw(j),xxx(10)
        ENDDO

        DO j=1,360
          if (sand(j).gt.0.0d0) then
      CALL WSPARAM(nupc,awl,kd,kx,nci,nfix,adp,sfc,sw,sswc,sand(j),
     &clay(j),silt(j),bulk(j)*1000000.0d0,100,topsl,fc,wp,swc,orgc(j))
          sdgvm_wilt(j) = wp
          sdgvm_field(j) = fc
          sdgvm_satw(j) = swc
*          print*,wp,fc,swc
*          stop
          endif
        ENDDO

        DO j=1,360
          lon = -180.5d0 + real(j)
          if (sand(j).gt.0.0d0) 
     &write(32,'(f7.3,f9.3,5f8.2,1x,3f8.4,1x,3f8.4)') lat,lon,
     &sand(j),silt(j),bulk(j),orgc(j),wilt(j),
     &field(j),satw(j),sdgvm_wilt(j),sdgvm_field(j),sdgvm_satw(j)
        ENDDO
      ENDDO


      stop
      end



*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE WSPARAM                          *
*                          ******************                          *
*                                                                      *
* WSPARAM sets some system parameters for both the water and soils     *
* models.                                                              *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE WSPARAM(nupc,awl,kd,kx,nci,nfix,adp,sfc,sw,sswc,sand,
     &clay,silt,bulk,dep,topsl,fc,wp,swc,org)
*----------------------------------------------------------------------*
      REAL*8 nupc,awl(4),kd,kx,nci(4),nfix,adp(4),sfc(4),sw(4),sswc(4)
      REAL*8 sand,clay,silt,bulk,dep,topsl,bc_wres,bc_wqs,ans1,fc,wp,swc
      REAL*8 bulkx,org

      bulkx = bulk*1.0e-6
      org = 0.0d0

      awl(1) = 0.05d0 ! Relative root density of the ith layer.
      awl(2) = 0.05d0 !                  "
      awl(3) = 0.03d0 !                  "
      awl(4) = 0.01d0 !                  "
      kd = 0.5d0    ! Fraction of excess water flowing to deep storage.
      kx = 0.2d0    !                   "                 stream.
      nci(1) = 0.04d0       ! NC ratio of litter pools.
      nci(2) = 0.06666667d0 !          "
      nci(3) = 0.04d0       !          "
      nci(4) = 0.06666667d0 !          "
      nupc = 1.0d0   ! Amount of nitrogen uptake ratio to litter.
      nfix = 0.5d0   ! Soil deposition plus N fixation (g/m2/y).

*************************************************************************
*                        Brookes Corey
*************************************************************************
      bc_wres = 0.0d0
      bc_wqs = 1.0-bulkx/2.65
      ans1 = -0.475662 +0.005165*sand + 0.002066*silt -
     &  0.023327*clay - 0.040628*bulkx +
     &  0.006824*org - 0.000136162*sand**2 -
     &  0.0000954716*silt**2 + 0.000298891*clay**2 -
     &  0.0637*bulkx**2 - 0.00031679*org**2 +
     &  0.0000010388*sand**3 + 0.0000010299*silt**3 -
     &  0.0000013251*clay**3

      bc_l = 10**ans1

      ans1 = 0.4104 + 0.002684*sand + 0.006354*silt +
     &0.17766*bulkx + 0.00013855*clay**2

      bc_a = 1/10**ans1

      wp = bc_wres + (bc_wqs - bc_wres)*(15300.0*bc_a)**(-bc_l)
      fc = bc_wres + (bc_wqs - bc_wres)*(51.0*bc_a)**(-bc_l)
      swc = 1.0-bulkx/2.65

*----------------------------------------------------------------------*
* Convert soil layer height into mm of water.                          *
*----------------------------------------------------------------------*
      adp(1) = topsl*10.0d0
      adp(2) = (dep/3.0d0 - topsl)*10.0d0
      adp(3) = dep*10.0d0/3.0d0
      adp(4) = dep*10.0d0/3.0d0
      sfc(1) = fc*adp(1)
*      sfc(1) = adp(1)
      sfc(2) = fc*adp(2)
      sfc(3) = fc*adp(3)
      sfc(4) = fc*adp(4)
      sw(1) = wp*adp(1)
      sw(2) = wp*adp(2)
      sw(3) = wp*adp(3)
      sw(4) = wp*adp(4)
      sswc(1) = swc*adp(1)
      sswc(2) = swc*adp(2)
      sswc(3) = swc*adp(3)
      sswc(4) = swc*adp(4)

*      CALL FCAP(ts,tc,adp,sfc,sw)


      RETURN
      END


