module photo
  
  !Module defining functions related with CO2 assimilation
  implicit none
  private

  ! functions defined here
  public ::                    &
       gross_ph               ,& ! gross photosynthesis (kgC m-2 y-1)
       laia                   ,& ! leaf area index(m2 m-2) 
       f_four                 ,& ! auxiliar function (calculates f4sun or f4shade or sunlai) 
       spec_leaf_area         ,& ! specific leaf area (m2 g-1)
       water_stress_modifier  ,& ! F5 - water stress modifier (dimensionless)
       photosynthesis_rate    ,& ! leaf level CO2 assimilation rate (molCO2 m-2 s-1)
       canopy_resistence      ,& ! Canopy resistence (from Medlyn et al. 2011a) (s/m) == m s-1
       vapor_p_defcit         ,& ! Vapor pressure defcit  (kPa)
       tetens                    ! Maximum vapor pressure (hPa)

contains
  
  !=================================================================
  !=================================================================

  function gross_ph(f1,cleaf,sla) result(ph)
    ! Returns gross photosynthesis rate (kgC m-2 y-1) 
    use global_pars, only: r4
    implicit none

    real(kind=r4),intent(in) :: f1    !mol(CO2) m-2 m-1 ≃ 1 cm s-1 == 100 mm s-1    
    real(kind=r4),intent(in) :: cleaf !kgC m-2 
    real(kind=r4),intent(in) :: sla   !m2 g-1
    real(kind=r4) :: ph
    
    real(kind=r4) :: f4sun
    real(kind=r4) :: f4shade
    
    f4sun = f_four(1,cleaf,sla)
    f4shade = f_four(2,cleaf,sla)
    ph = 0.012*31557600.0*f1*f4sun*f4shade
    
  end function gross_ph
  
  !=================================================================
  !=================================================================
  
  function laia(cleaf,sla) result(lai)
    ! Returns Leaf Area Index m2 m-2
    
    use global_pars, only: r4
    implicit none
    
    real(kind=r4),intent(in) :: cleaf !kgC m-2 
    real(kind=r4),intent(in) :: sla   !m2 g-1
    real(kind=r4) :: lai
    
    lai  = ((cleaf * 1000.)  * sla) ! * 1000 transform kg to g - laia64 in m2 m-2
    
  end function laia
  
  !=================================================================
  !=================================================================
  
  function f_four(fs,cleaf,sla) result(lai_ss)
    use global_pars, only: i4, r4, r8
    use photo_par, only: p26, p27
    implicit none
    
    integer(kind=i4),intent(in) :: fs !function mode 1 == f4sun; 2 == f4shade; >2 == f4
    
    real(kind=r4),intent(in) :: cleaf ! carbon in leaf (kg m-2)
    real(kind=r4),intent(in) :: sla   ! specific leaf area (m2 g-1)
    real(kind=r4) :: lai_ss           ! leaf area index (m2 m-2)
    
    real(kind=r4) :: lai
    real(kind=r8) :: sunlai
    real(kind=r8) :: shadelai
    
    lai = laia(cleaf,sla) 
    sunlai = (1.0-(exp(-p26*lai)))/p26
    lai_ss = real(sunlai,r4)  ! to be used for heterothrophic respiration
    shadelai = lai - sunlai
    
    !Scaling-up to canopy level (dimensionless)
    !------------------------------------------
    lai_ss = real(sunlai,r4)  ! to be used for heterothrophic respiration
    
    if(fs .gt. 2) then 
       return
    endif
   
    !Sun/Shade approach to canopy scaling !Based in de Pury & Farquhar (1997)
    !------------------------------------------------------------------------
    if(fs .eq. 1) then
       ! f4sun
       lai_ss = real((1.0-(exp(-p26*sunlai)))/p26,r4) !sun 90 degrees
       return
    endif
    
    if(fs .eq. 2) then
       !f4shade
       lai_ss = real((1.0-(exp(-p27*shadelai)))/p27,r4) !sun ~20 degrees
       return
    endif
  end function f_four
 
  !=================================================================
  !================================================================= 
  
  function spec_leaf_area(tau_leaf) result(sla)
    ! based on JeDi DGVM 
    use global_pars, only : r4
    implicit none
    
    real(kind=r4),intent(in) :: tau_leaf
    real(kind=r4) :: sla
    
    real(kind=r4) :: leaf_t_months
    real(kind=r4) :: leaf_t_coeff
    real(kind=r4) :: leaf_turnover
    
    !Leaf area index (m2/m2)
    leaf_t_months = tau_leaf*12. ! turnover time in months
    leaf_t_coeff = leaf_t_months/100. !1 - 100 months == ~ 1/12 to 8.3 years (TRY-kattge et al. 2011; Jedi-Pavlick 2012)
    if (leaf_t_coeff .gt. 1.) leaf_t_coeff = 1. 
    leaf_turnover =  (365.0/12.0) * (10. **(2.0*leaf_t_coeff))
    sla = (3e-2 * (365.0/leaf_turnover)**(-0.46))     
    
  end function spec_leaf_area
  
  !=================================================================
  !=================================================================
  
  function water_stress_modifier(w, cfroot, rc, ep) result(f5)
    use global_pars, only: r4, r8, csru, wmax, alfm, gm, rcmin
    implicit none
    
    real(kind=r4),intent(in) :: w      !soil water mm
    real(kind=r4),intent(in) :: cfroot !carbon in fine roots kg m-2
    real(kind=r4),intent(in) :: rc     !Canopy resistence 1/(micromol(CO2) m-2 s-1)
    real(kind=r4),intent(in) :: ep     !potential evapotranspiration
    real(kind=r4) :: f5
    
    
    real(kind=r8) :: pt
    real(kind=r8) :: gc
    real(kind=r8) :: wa
    real(kind=r8) :: d
    real(kind=r8) :: f5_64
    
    wa = w/wmax
    
    pt = csru*(cfroot*1000.)*wa  !(based in Pavlick et al. 2013; *1000. converts kgC/m2 to gC/m2)
    if(rc .gt. 0.0) then
       gc = (1./rc)  ! s/m
    else
       gc =  1.0/rcmin ! BIANCA E HELENA - Mudei este esquema..   
    endif                     ! tentem entender o algoritmo
    
    d =(ep * alfm) / (1. + gm/gc) !(based in Gerten et al. 2004)
    
    if(d .gt. 0.0) then
       f5_64 = pt/d
       f5_64 = exp(-1 * f5_64)
       f5_64 = 1.0 - f5_64
    else
       f5_64 = wa
    endif
    
    f5 = real(f5_64,4)      
  end function water_stress_modifier
  
  !=================================================================
  !=================================================================
  
  function canopy_resistence(vpd_in,f1_in,g1) result(rc2_in)
    ! return stomatal resistence based on Medlyn et al. 2011a
    use global_pars, only:r4 ,r8, ca, rcmax, rcmin
    
    implicit none
    
    real(kind=r4),intent(in) :: f1_in    !Photosynthesis (molCO2/m2/s)
    real(kind=r4),intent(in) :: vpd_in   !hPa
    real(kind=r4),intent(in) :: g1      
    real(kind=r4) :: rc2_in              !Canopy resistence (s/m)
    
    !     Internal
    !     --------
    real(kind=r8) :: f1b      !Photosynthesis (micromolCO2/m2/s)
    real(kind=r8) :: gs2      !Canopy conductance (S/M)
    real(kind=r8) :: gs       !Canopy RESISTENCE (molCO2/m2/s)
    real(kind=r8) :: g0       !MINIMUM stomatal conductance FITTED PARAMETER (DIMENSIONLESS)
    real(kind=r8) :: D1       !kPA
    real(kind=r8) :: aa
    
    f1b = f1_in * 1e6        ! mol 2 µmol : conversion factor = 1e6
    aa = f1b / ca              
    g0 = 0.00003           
    
    if(f1_in .le. 0.0) then 
       rc2_in = rcmax
       return
    endif
    if(vpd_in .gt. 0.1) then
       goto 10
    else
       rc2_in = rcmin
       return
    endif
10  continue
    if (vpd_in .gt. 1.5) then
       rc2_in = rcmax
       return
    else
       D1 = sqrt(vpd_in) 
       gs = g0 + 1.6 * (1. + (g1/D1)) * (aa) !Based on Medlyn et al. 2011
       if(gs .le. 0.0) then
          rc2_in = rcmax
          return
       endif
    endif
    gs2 = (gs/41.)
    if(gs2 .le. 0.0) rc2_in = rcmax
    if(gs2 .gt. 0.0) rc2_in = real((gs2**(-1)),4)
    if(rc2_in .lt. rcmin) rc2_in = rcmin
  end function canopy_resistence
  
  !=================================================================
  !=================================================================

  function vapor_p_defcit(t,rh) result(vpd_0)
    use global_pars, only: r4
    implicit none
    
    real(kind=r4),intent(in) :: t
    real(kind=r4),intent(in) :: rh
    
    real(kind=r4) :: es
    real(kind=r4) :: vpd_ac
    real(kind=r4) :: vpd_0
    ! ext func
    !real(kind=r4) :: tetens
    
    es = tetens(t)  
    !VPD-REAL = Actual vapor pressure
    vpd_ac = es * rh       ! RESULT in hPa == mbar! we want kPa (DIVIDE by 10.)
    !Vapor Pressure Deficit
    vpd_0 = (es - vpd_ac) / 10.
  end function vapor_p_defcit

  !=================================================================
  !=================================================================
  
  function photosynthesis_rate(vm,temp,p0,ipar,ll) result(f1a)
    !returns instantaneous photosynthesis rate at leaf level (molCO2/m2/s)
    use global_pars
    use photo_par
    implicit none
    
    real(kind=r4),intent(in) :: vm
    real(kind=r4),intent(in) :: temp
    real(kind=r4),intent(in) :: p0
    real(kind=r4),intent(in) :: ipar
    logical(kind=l1),intent(in) :: ll
    
    real(kind=r8) :: f2,f3            !Michaelis-Menten CO2/O2 constant (Pa)
    real(kind=r8) :: mgama           !Photo-respiration compensation point (Pa)
    real(kind=r8) :: rmax, r
    real(kind=r8) :: ci
    real(kind=r8) :: jp1
    real(kind=r8) :: jp2
    real(kind=r8) :: jp
    real(kind=r8) :: jc
    real(kind=r8) :: jl
    real(kind=r8) :: je
    real(kind=r8) :: b,c,c2,b2,es,j1,j2
    real(kind=r8) :: delta, delta2,aux_ipar
    real(kind=r8) :: f1a
    ! ext func
    !real(kind=r4) :: tetens
    
    !============================================================
    !Photo-respiration compensation point (Pa)
    mgama = p3/(p8*(p9**(p10*(temp-p11))))
    
    !Michaelis-Menten CO2 constant (Pa)     
    f2 = p12*(p13**(p10*(temp-p11)))
    
    !Michaelis-Menten O2 constant (Pa)
    f3 = p14*(p15**(p10*(temp-p11)))

    !Saturation vapour pressure (hPa)
    es = real(tetens(temp), kind=r8)
    
    !Saturated mixing ratio (kg/kg)     
    rmax = 0.622*(es/(p0-es))
    
    !Moisture deficit at leaf level (kg/kg)    
    r = -0.315*rmax
    
    !Internal leaf CO2 partial pressure (Pa)
    ci = p19* (1.-(r/p20)) * (ca-mgama) + mgama
    !==============================================================
    
    !Rubisco carboxilation limited photosynthesis rate (molCO2/m2/s)
    jc = vm*((ci-mgama)/(ci+(f2*(1.+(p3/f3)))))
    
    !Light limited photosynthesis rate (molCO2/m2/s)  
    if (ll) then
       aux_ipar= ipar * 1.0_r8  
    else
       aux_ipar = ipar - (ipar * 0.20_r8)
    endif
    jl = p4*(1.0-p5)*aux_ipar*((ci-mgama)/(ci+(p6*mgama)))
    
    ! Transport limited photosynthesis rate (molCO2/m2/s) (RuBP) (re)generation
    ! ---------------------------------------------------
    je = p7*vm
    
    !Jp (minimum between jc and jl)
    !------------------------------
    b = (-1.)*(jc+jl)
    c = jc*jl
    delta = (b**2)-4.0*a*c
    if(delta .eq. 0.0_r8)then
       jp = (-b) / (2 * a)
    else if(delta .gt. 0.0_r8) then
       jp1 = (-b-(sqrt(delta)))/(2.0*a)
       jp2 = (-b+(sqrt(delta)))/(2.0*a)
       jp = amin1(jp1,jp2)
    else
       jp = 0.0_r8
    endif
    
    
    !Leaf level gross photosynthesis (minimum between jc, jl and je)
    !---------------------------------------------------------------
    b2 = (-1.)*(jp+je)
    c2 = jp*je
    delta2 = (b2**2)-4.0*a2*c2
    if(delta2 .eq. 0.0_r8)then
       f1a = (-b2) / (2.0 * a2)
    else if(delta2 .gt. 0.0_r8) then
       j1 = (-b2-(sqrt(delta2)))/(2.0*a2)
       j2 = (-b2+(sqrt(delta2)))/(2.0*a2)
       f1a = amin1(j1,j2)
    else
       f1a = 0.0_r8
    endif
    
  end function photosynthesis_rate

  !=================================================================
  !=================================================================
  
  function tetens(t) result(es)
    ! returns Saturation Vapor Pressure (hPa), using Buck equation
    
    ! buck equation...references:
    ! http://www.hygrometers.com/wp-content/uploads/CR-1A-users-manual-2009-12.pdf
    ! Hartmann 1994 - Global Physical Climatology p.351
    ! https://en.wikipedia.org/wiki/Arden_Buck_equation#CITEREFBuck1996
    
    use global_pars
    implicit none
    
    real(kind=r4),intent( in) :: t
    real(kind=r4) :: es
    
    if (t .ge. 0.) then
       es = 6.1121 * exp((18.678-(t/234.5))*(t/(257.14+t))) ! mbar == hPa
       return
    else
       es = 6.1115 * exp((23.036-(t/333.7))*(t/(279.82+t))) ! mbar == hPa
       return
    endif
    
  end function tetens
  
end module photo
