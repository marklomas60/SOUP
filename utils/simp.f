      real*8 sand,silt,clay,bulk,ans1
      real*8 a_saxton,b_saxton,wp_saxton,fc_saxton,sw_saxton
      real*8 bc_wres,bc_wqs,bs_a,bc_l,wp_bc,fc_bc,sw_bc
      real*8 wp_ch,fc_ch,sw_ch,b,lsmp,smp

      sand = 50.0d0
      silt = 30.0d0
      clay = 100.0d0 - sand - silt
      bulk = 1.33d0
      org = 10.0d0

*************************************************************************
*                            Saxton
*************************************************************************
      b_saxton=-3.14-0.00222*clay**2-3.484/10**5*sand**2*clay
      a_saxton=100.0*exp(-4.396-0.0715*clay-4.88/10**4*sand**2-4.28/
     &(10**5)*(sand**2)*clay)
      wp_saxton = (1500.0/a_saxton)**(1/b_saxton)
      fc_saxton = exp((2.302-log(a_saxton))/b_saxton)
      sw_saxton = 0.332 - 0.0007251*sand + 0.1276*log10(clay)

*************************************************************************
*                        Brookes Corey
*************************************************************************
      bc_wres = 0.0d0
      bc_wqs = 1.0-bulk/2.65
      ans1 = -0.475662 +0.005165*sand + 0.002066*silt -
     &  0.023327*clay -0.040628*bulk +
     &  0.006824*org -0.000136162*sand**2 -
     &  0.0000954716*silt**2 +0.000298891*clay**2 -
     &  0.0637*bulk**2 -0.00031679*org**2 +
     &  0.0000010388*sand**3 +0.0000010299*silt**3 -
     &  0.0000013251*clay**3

      bc_l = 10**ans1

      ans1 = 0.4104+0.002684*sand+0.006354*silt+
     &0.17766*bulk+0.00013855*clay**2

      bc_a = 1/10**ans1

      wp_bc = bc_wres + (bc_wqs - bc_wres)*(15300.0*bc_a)**(-bc_l)
      fc_bc = bc_wres + (bc_wqs - bc_wres)*(51.0*bc_a)**(-bc_l)
      sw_bc = 1.0-bulk/2.65

**************************************************************************
*                Cosby Hornberger Clapp Ginn                      
**************************************************************************  
      b = 3.1d0 + (0.157d0*clay) - (0.003d0*sand)
      swc = 50.5d0 - (0.142d0*sand) - (0.037d0*clay)
      lsmp = 1.54d0 - (0.0095d0*sand) + (0.0063d0*silt)
      smp = 10.0d0**lsmp
      fc_ch = ((((306.0d0/smp)**(-1.0d0/b))*swc)/100.0d0)*0.78d0*bulk
      wp_ch = ((((15300.0d0/smp)**(-1.0d0/b))*swc)/100.0d0)*0.78d0*bulk
      sw_ch = (((306.0d0/smp)**(-1.0d0/b))*swc/100.d0)*bulk

      print*,wp_saxton,fc_saxton,sw_saxton
      print*,wp_bc,fc_bc,sw_bc
      print*,wp_ch,fc_ch,sw_ch

      print*,sw_saxton-wp_saxton,fc_saxton-wp_saxton
      print*,sw_bc-wp_bc,fc_bc-wp_bc
      print*,sw_ch-wp_ch,fc_ch-wp_ch


      stop
      end

