module photo
  
  !Module defining functions related with CO2 assimilation
  implicit none
  private

  ! functions defined here
  public ::                    &
       gross_ph               ,& ! gross photosynthesis (kgC m-2 y-1)
       leaf_area_index        ,& ! leaf area index(m2 m-2) 
       f_four                 ,& ! auxiliar function (calculates f4sun or f4shade or sunlai) 
       spec_leaf_area         ,& ! specific leaf area (m2 g-1)
       water_stress_modifier  ,& ! F5 - water stress modifier (dimensionless)
       photosynthesis_rate    ,& ! leaf level CO2 assimilation rate (molCO2 m-2 s-1)
       canopy_resistence      ,& ! Canopy resistence (from Medlyn et al. 2011a) (s/m) == m s-1
       vapor_p_defcit         ,& ! Vapor pressure defcit  (kPa)
       tetens                 ,& ! Maximum vapor pressure (hPa)
       m_resp                 ,& ! maintenance respiration (plants)
       g_resp                 ,& ! growth Respiration (kg m-2 yr-1)
       carbon2                ,& ! soil + litter + heterothrophic respiration
       pft_area_frac          ,& ! area fraction by biomass
       pft_par                ,& ! aux subroutine to read pls data
       spinup                 ,&
       ascii2bin              ,&
       allocation
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
  
  function leaf_area_index(cleaf,sla) result(lai)
    ! Returns Leaf Area Index m2 m-2
    
    use global_pars, only: r4
    implicit none
    
    real(kind=r4),intent(in) :: cleaf !kgC m-2 
    real(kind=r4),intent(in) :: sla   !m2 g-1
    real(kind=r4) :: lai
    
    lai  = ((cleaf * 1000.)  * sla) ! * 1000 transform kg to g - laia64 in m2 m-2
    
  end function leaf_area_index
  
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
    
    lai = leaf_area_index(cleaf,sla) 
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
       gc = (1.0_r8/rc)  ! s/m
    else
       gc =  1.0_r8/rcmin ! BIANCA E HELENA - Mudei este esquema..   
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
  
  function photosynthesis_rate(vm,temp,p0,ipar,ll) result(f1ab)
    !returns instantaneous photosynthesis rate at leaf level (molCO2/m2/s)
    use global_pars
    use photo_par
    implicit none
    
    real(kind=r4),intent(in) :: vm
    real(kind=r4),intent(in) :: temp
    real(kind=r4),intent(in) :: p0
    real(kind=r4),intent(in) :: ipar
    logical(kind=l1),intent(in) :: ll
    real(kind=r4) :: f1ab
        

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
    ci = p19* (1.-(r/p20)) * ((ca/9.901)-mgama) + mgama
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

   f1ab = real(f1a,kind=r4)
    
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


  !====================================================================
  !====================================================================

  
  function m_resp(temp,cl1,cf1,ca1,ocp_coeff) result(rm)
    use global_pars, only: r4,r8,ncl,ncf,ncs
    implicit none

    real(kind=r4), intent(in) :: temp
    real(kind=r4), intent(in) :: cl1
    real(kind=r4), intent(in) :: cf1
    real(kind=r4), intent(in) :: ca1
    real(kind=r4), intent(in) :: ocp_coeff
    real(kind=r4) :: rm
    
    real(kind=r8) :: csa, rm64, rml64, rmf64, rms64
    
    
    !   Autothrophic respiration
    !   ========================
    !   Maintenance respiration (kgC/m2/yr) (based in Ryan 1991)
    
    csa= 0.05_r8 * ca1           !sapwood carbon content (kgC/m2). 5% of woody tissues (Pavlick, 2013)
    
    rml64 = (ncl * cl1) * 27. * exp(0.03*temp)
    
    rmf64 = (ncf * cf1) * 27. * exp(0.03*temp)
    
    rms64 = (ncs * csa) * 27. * exp(0.03*temp)
    
    rm64 = (rml64 + rmf64 + rms64)
    
    rm = real(rm64,kind=r4) * ocp_coeff

    if (rm.lt.0) then
       rm = 0.0
    endif

  end function m_resp

  !====================================================================
  !====================================================================
  
  function g_resp(beta_leaf,beta_awood, beta_froot, ocp_coeff) result(rg)
    use global_pars, only: r4,r8
    implicit none

    real(kind=r4), intent(in) :: beta_leaf
    real(kind=r4), intent(in) :: beta_froot
    real(kind=r4), intent(in) :: beta_awood
    real(kind=r4), intent(in) :: ocp_coeff
    real(kind=r4) :: rg
    
    real(kind=r8) :: csai, rg64, rgl64, rgf64, rgs64
    
    !     Autothrophic respiration
    !     Growth respiration (KgC/m2/yr)(based in Ryan 1991; Sitch et al.
    !     2003; Levis et al. 2004)         
    
    csai =  (beta_awood * 0.05)
    
    rgl64 = 0.25 * beta_leaf * 365.
    
    rgf64 =  0.25* beta_froot * 365.
    
    rgs64 = (0.25 * csai * 365.)
    
    rg64 = (rgl64 + rgf64 + rgs64)

    rg = real(rg64,kind=r4) * ocp_coeff
    
    if (rg.lt.0) then
       rg = 0.0
    endif
    
  end function g_resp
  
  !====================================================================
  !====================================================================
  
  subroutine carbon2 (tsoil,f5c,evap,laia,d_litter,cl,cs,hr)
    use global_pars
    use photo_par
    implicit none
    !     Variables
    !     =========
    !     Inputs
    !     ------
    real(kind=r4),intent(in) :: tsoil                !Mean monthly soil temperature (oC)
    real(kind=r4),intent(in) :: f5c                  !Stress response to soil moisture (dimensionless)
    real(kind=r4),intent(in) :: evap                 !Actual evapotranspiration (mm/day)
    real(kind=r4),intent(in) :: laia
    real(kind=r4),intent(in) :: d_litter
    !     Outputs 
    !     -------
    real(kind=r4),intent(out) :: cl                   !Litter carbon (kgC/m2)
    real(kind=r4),intent(out) :: cs                   !Soil carbon (kgC/m2)
    real(kind=r4),intent(out) :: hr                   !Heterotrophic (microbial) respiration (kgC/m2)
    
    !     Internal
    !     --------
    real(kind=r8) :: lf                   !Litterfall (kgC/m2)
    real(kind=r8) :: f6                   !Litter decayment function
    real(kind=r8) :: f7                   !Soil carbon storage function
    !     
    !     Initialize
    !     ----------
    !     
    lf  = 0.0
    f6  = 0.0
    f7  = 0.0
    cl  = 0.0
    cs  = 0.0
    
    !     Litter decayment function                                             !Controlled by annual evapotranspiration
    !     -------------------------
    f6 = 1.16*10.**(-1.4553+0.0014175*(evap*365.0))
    
    !     Soil carbon storage function                                          !Controlled by temperature
    !     ----------------------------    
    f7 = p32**(p10*(tsoil-p11))
    
    !     Litterfall (kgC/m2)
    !     ------------------
    lf = p33 * (laia + d_litter)
    
    !     Litter carbon (kgC/m2)
    !     ----------------------  
    cl = real(lf/f6, kind=r4)
    
    !     Soil carbon(kgC/m2)
    !     -------------------
    cs = real(((p34*cl)/(p35*f7))*f5c, kind=r4)
    
    !     Respiration minimum and maximum temperature
    !     -------------------------------------------
    !     
    if ((tsoil.ge.-10.0).and.(tsoil.le.50.0)) then
       hr = real(p36*(cl*(f6**2)+(cs*f5c*evap*(f7**2))),kind=r4) !Litter and Soil!respectively
    else
       hr = 0.0_r4               !Temperature above/below respiration windown
    endif
  
  end subroutine carbon2
  
  !====================================================================
  !====================================================================
  
  SUBROUTINE PFT_AREA_FRAC(CLEAF, CFROOT, CAWOOD, OCP_COEFFS, OCP_WOOD)
    use global_pars
    implicit none
    
    integer(kind=i4),parameter :: npft = npls ! plss futuramente serao
    
    REAL(kind=r4),dimension(npft),intent( in) :: CLEAF, CFROOT, CAWOOD
    REAL(kind=r4),dimension(npft),intent(out) :: OCP_COEFFS
    logical,dimension(npft),intent(out) :: OCP_WOOD
    REAL(kind=r4),dimension(npft) :: TOTAL_BIOMASS_PFT,TOTAL_W_PFT
    INTEGER(kind=i4) :: P,I
    INTEGER(kind=i4),dimension(1) :: MAX_INDEX
    REAL(kind=r4) :: TOTAL_BIOMASS, TOTAL_WOOD
    
    TOTAL_BIOMASS = 0.0
    TOTAL_WOOD = 0.0
    
    do p = 1,npft
       TOTAL_W_PFT(P) = 0.0
       TOTAL_BIOMASS_PFT(P) = 0.0
       OCP_COEFFS(P) = 0.0
       OCP_WOOD(P) = .FALSE.
    enddo
    
    DO P = 1,NPFT
       TOTAL_BIOMASS_PFT(P) = CLEAF(P) + CFROOT(P) + CAWOOD(P) ! only sapwood
       TOTAL_BIOMASS = TOTAL_BIOMASS + TOTAL_BIOMASS_PFT(P)
       TOTAL_WOOD = TOTAL_WOOD + CAWOOD(P)
       TOTAL_W_PFT(P) = CAWOOD(P)
    ENDDO
    
    !     GRID CELL OCCUPATION COEFFICIENTS
    IF(TOTAL_BIOMASS .GT. 0.0) THEN
       DO P = 1,NPFT   
          OCP_COEFFS(P) = TOTAL_BIOMASS_PFT(P) / TOTAL_BIOMASS
          IF(OCP_COEFFS(P) .LT. 0.0) OCP_COEFFS(P) = 0.0
       ENDDO
    ELSE
       DO P = 1,NPFT
          OCP_COEFFS(P) = 0.0
       ENDDO
    ENDIF
    
    !     GRIDCELL PFT LIGTH LIMITATION BY WOOD CONTENT 
    IF(TOTAL_WOOD .GT. 0.0) THEN
       MAX_INDEX = MAXLOC(TOTAL_W_PFT)
       I = MAX_INDEX(1)
       OCP_WOOD(I) = .TRUE.
    ENDIF
    
  END SUBROUTINE PFT_AREA_FRAC
  
  !====================================================================
  !====================================================================

  subroutine pft_par(par, dt)
    use global_pars
    implicit none
    
    
    integer(kind=i4),intent(in) :: par            ! parameter number 
    real(kind=r4), dimension(npls),intent(out) :: dt
    
    ! ['g1','vcmax','tleaf','twood','troot','aleaf','awood','aroot']
    !     dt1 = g1
    !     dt2 = vcmax
    !     dt3 = tleaf
    !     dt4 = twood
    !     dt5 = tfroot
    !     dt6 = aleaf
    !     dt7 = awood
    !     dt8 = aroot
    
    open(45,file='./pls.bin',status='old',&
         &form='unformatted',access='direct',recl=4*npls)
    
    if(par .gt. 0 .and. par .lt. 10) then
       read(45,rec=par) dt
    else
       print*, 'search failed'
    endif
    close(45)
    return
  end subroutine pft_par
  
  !====================================================================
  !====================================================================
  
  subroutine spinup(nppot,cleafini,cfrootini,cawoodini)
    use global_pars
    implicit none
    
    integer(kind=i4),parameter :: npfts = npls
    integer(kind=i4),parameter :: ntl=30000
    
    !   c     inputs
    integer(kind=i4) :: i6, kk, k
    
    real(kind=r4),intent(in) :: nppot
    real(kind=r4) :: sensitivity
    
    !c     outputs
    real(kind=r4),intent(out) :: cleafini(npfts)
    real(kind=r4),intent(out) :: cawoodini(npfts)
    real(kind=r4),intent(out) :: cfrootini(npfts)
    real(kind=r8) :: cleafi_aux(ntl)
    real(kind=r8) :: cfrooti_aux(ntl)
    real(kind=r8) :: cawoodi_aux(ntl)
    
    
    real(kind=r4) :: aleaf(npfts)             !npp percentage alocated to leaf compartment
    real(kind=r4) :: aawood (npfts)           !npp percentage alocated to aboveground woody biomass compartment
    real(kind=r4) :: afroot(npfts)            !npp percentage alocated to fine roots compartmentc 
    real(kind=r4) :: tleaf(npfts)             !turnover time of the leaf compartment (yr)
    real(kind=r4) :: tawood (npfts)           !turnover time of the aboveground woody biomass compartment (yr)
    real(kind=r4) :: tfroot(npfts)            !turnover time of the fine roots compartment
    
    
    call pft_par(6, aleaf)
    call pft_par(7, aawood)
    call pft_par(8, afroot)
    call pft_par(3, tleaf)
    call pft_par(4, tawood)
    call pft_par(5, tfroot)
    
    
    sensitivity = 1.01
    if(nppot .lt. 0.0) goto 200
    ! nppot = nppot/real(npfts,kind=r4)
    do i6=1,npfts
       do k=1,ntl
          if (k.eq.1) then
             cleafi_aux (k) =  aleaf(i6)*(nppot)
             cawoodi_aux(k) = aawood(i6)*(nppot)
             cfrooti_aux(k) = afroot(i6)*(nppot)
             
          else
             if(aawood(i6) .gt. 0.0) then
                cleafi_aux(k) = ((aleaf(i6)*(nppot))-(cleafi_aux(k-1)&
                     &/(tleaf(i6)))) + cleafi_aux(k-1)
                cawoodi_aux(k) = ((aawood(i6)*(nppot))-(cawoodi_aux(k&
                     &-1)/(tawood(i6)))) + cawoodi_aux(k-1)
                cfrooti_aux(k) = ((afroot(i6)*(nppot))-(cfrooti_aux(k&
                     &-1)/(tfroot(i6)))) + cfrooti_aux(k-1)
             else
                cleafi_aux(k) = ((aleaf(i6)*(nppot))-(cleafi_aux(k-1)&
                     &/(tleaf(i6)))) + cleafi_aux(k-1)
                cawoodi_aux(k) = 0.0
                cfrooti_aux(k) = ((afroot(i6)*(nppot))-(cfrooti_aux(k&
                     &-1)/(tfroot(i6)))) + cfrooti_aux(k-1)
             endif
             
             kk =  nint(k*0.66)
             if(cawoodi_aux(kk) .gt. 0.0) then
                if((cfrooti_aux(k)/cfrooti_aux(kk).lt.sensitivity).and.&
                     &(cleafi_aux(k)/cleafi_aux(kk).lt.sensitivity).and.&
                     &(cawoodi_aux(k)/cawoodi_aux(kk).lt.sensitivity)) then
                   
                   cleafini(i6) = real(cleafi_aux(k),4) ! carbon content (kg m-2)
                   cfrootini(i6) = real(cfrooti_aux(k),4)
                   cawoodini(i6) = real(cawoodi_aux(k),4)
                   exit
                ENDIF
             else
                if((cfrooti_aux(k)&
                     &/cfrooti_aux(kk).lt.sensitivity).and.&
                     &(cleafi_aux(k)/cleafi_aux(kk).lt.sensitivity)) then
                   
                   cleafini(i6) = real(cleafi_aux(k),4) ! carbon content (kg m-2)
                   cfrootini(i6) = real(cfrooti_aux(k),4)
                   cawoodini(i6) = 0.0
                   exit
                endif
             endif
          endif
       enddo                  !nt
    enddo                     ! npfts 
200 continue
  end subroutine spinup
  
 !=================================================================
 !=================================================================
 
  subroutine ascii2bin(file_in, file_out, nx1, ny1)
    use global_pars
    implicit none
    
    character*30, intent(in) :: file_in, file_out
    integer(kind=i4),intent(in) :: nx1, ny1
    
    integer(kind=i4) :: i, j
    
    real(kind=r4),allocatable,dimension(:,:) :: arr_in
    
    allocate(arr_in(nx1,ny1))
    
    open (unit=11,file=file_in,status='old',form='formatted',access='sequential',&
         action='read')
    
    
    open (unit=21,file=file_out,status='unknown',&
         form='unformatted',access='direct',recl=nx1*ny1*4)
    
    
    do j = 1, ny1 ! for each line do
       read(11,*) (arr_in(i,j), i=1,nx1) ! read all elements in line j (implicit looping)
       !write(*,*) arr_in(:,j) 
    end do
    
    write(21,rec=1) arr_in   
    close(11)
    close(21)
    
  end subroutine ascii2bin
  !====================================================================
  !====================================================================

   !=====================================================================
 !c     subroutine allocation calculates the daily carbon content of each
 !c     compartment
 !c     
 !c     code written by Bianca Rius & David Lapola (27.Ago.2015)
 !c     
 !c=====================================================================
      
  subroutine allocation (pft,npp,scl1,sca1,scf1,scl2,sca2,scf2,bio_litter)
    use global_pars
    implicit none
    
    integer(kind=i4),parameter :: npfts = npls
    
    !     variables
    integer(kind=i4),intent(in) :: pft   
    real(kind=r4),intent(in) :: npp   !potential npp (KgC/m2/yr)
    real(kind=rbig) :: npp_aux        !auxiliary variable to calculate potential npp in KgC/m2/day
    real(kind=r4),intent( in) :: scl1 !previous day carbon content on leaf compartment (KgC/m2)
    real(kind=r4),intent(out) :: scl2 !final carbon content on leaf compartment (KgC/m2)
    real(kind=r4),intent( in) :: sca1 !previous day carbon content on aboveground woody biomass compartment(KgC/m2)
    real(kind=r4),intent(out) :: sca2 !final carbon content on aboveground woody biomass compartment (KgC/m2)
    real(kind=r4),intent( in) :: scf1 !previous day carbon content on fine roots compartment (KgC/m2)
    real(kind=r4),intent(out) :: scf2 !final carbon content on fine roots compartment (KgC/m2)      
    real(kind=r4),intent(out) :: bio_litter
    real(kind=rbig) :: scf2_128 = 0.0, sca2_128 = 0.0, scl2_128 = 0.0
    
    real(kind=r4), dimension(npfts) :: aleaf             !npp percentage allocated compartment
    real(kind=r4), dimension(npfts) :: aawood
    real(kind=r4), dimension(npfts) :: afroot
    real(kind=r4), dimension(npfts) :: tleaf             !turnover time (yr)
    real(kind=r4), dimension(npfts) :: tawood
    real(kind=r4), dimension(npfts) :: tfroot            
    
    call pft_par(6, aleaf)
    call pft_par(7, aawood)
    call pft_par(8, afroot)
    call pft_par(3, tleaf)
    call pft_par(4, tawood)
    call pft_par(5, tfroot)
    
    !     Carbon content of each compartment(KgC/m2)
    !c     
    !c     
    !c     initialization
    if(((scl1 .lt. 1e-12) .or. (scf1 .lt. 1e-12)) .and. (sca1 .gt. 0.0)) then
       bio_litter = scl1 + scf1 + sca1
       scl2 = 0.0
       scf2 = 0.0
       sca2 = 0.0
      goto 10
   endif
   npp_aux = npp/365.0       !transform (KgC/m2/yr) in (KgC/m2/day)
   scl2_128 = scl1 + (aleaf(pft) * npp_aux) -(scl1 /(tleaf(pft)*365.0))
   scf2_128 = scf1 +(afroot(pft) * npp_aux)-(scf1 /(tfroot(pft)*365.0))
   if(aawood(pft) .gt. 0.0) then
      sca2_128 = sca1 +(aawood(pft)*npp_aux)-(sca1/(tawood(pft)*365.0))
      sca2 = real(sca2_128,r4)
   else
      sca2 = 0.0
   endif
   
   scf2 = real(scf2_128,r4)
   scl2 = real(scl2_128,r4)
   
   
   if(scl2_128 .lt. 1e-12) scl2 = 0.0
   if(scf2_128 .lt. 1e-12) scf2 = 0.0
   if(sca2_128 .lt. 1e-12) sca2 = 0.0
   
10 continue
   return
   
 end subroutine allocation
 
end module photo

! =================-----------------==================----------------

module water
  
  ! this module defines functions related to surface water balance
  implicit none
  private

  ! functions defined here:
  
  public ::              &
       soil_temp        ,&
       soil_temp_sub    ,&
       penman           ,&
       evpot2           ,&
       available_energy ,&
       runoff          

  
contains
  
  !=================================================================
  !=================================================================

  subroutine soil_temp_sub(temp, tsoil)
  ! Calcula a temperatura do solo. Aqui vamos mudar no futuro!
  ! a tsoil deve ter relacao com a et realizada...
  ! a profundidade do solo (H) e o coef de difusao (DIFFU) devem ser
  ! variaveis (MAPA DE SOLO?; agua no solo?)
  use global_pars
  implicit none
  integer(kind=i4),parameter :: m = ntimes
  
  real(kind=r4),dimension(m), intent( in) :: temp ! future __ make temps an allocatable array
  real(kind=r4),dimension(m), intent(out) :: tsoil
   
  ! internal vars
  
  integer(kind=i4) :: n, k
  real(kind=r4) :: t0 = 0.0
  real(kind=r4) :: t1 = 0.0

  tsoil = -9999.0

  do n=1,1200 !run to attain equilibrium
     k = mod(n,m)
     if (k.eq.0) k = 12
     t1 = (t0*exp(-1.0/TAU) + (1.0 - exp(-1.0/TAU)))*temp(k)
     tsoil(k) = (t0 + t1)/2.0
     t0 = t1
  enddo
end subroutine soil_temp_sub
  
  !=================================================================
  !=================================================================
  
  function soil_temp(t0,temp) result(tsoil)
    use global_pars, only: i4, r4, H, TAU, DIFFU
    implicit none
    
    real(kind=r4),intent( in) :: temp
    real(kind=r4),intent( in) :: t0! future __ make temps an allocatable array
    real(kind=r4) :: tsoil
 
    real(kind=r4) :: t1 = 0.0
    
    t1 = (t0*exp(-1.0/TAU) + (1.0 - exp(-1.0/TAU)))*temp
    tsoil = (t0 + t1)/2.0
  end function soil_temp

  !=================================================================
  !=================================================================
  
  function penman (spre,temp,ur,rn,rc2) result(evap)
    use global_pars, only: r4, rcmin, rcmax
    use photo, only: tetens
    implicit none
    
    
    real(kind=r4),intent(in) :: spre                 !Surface pressure (mbar)
    real(kind=r4),intent(in) :: temp                 !Temperature (oC)
    real(kind=r4),intent(in) :: ur                   !Relative humidity (0-1,dimensionless)
    real(kind=r4),intent(in) :: rn                   !Radiation balance (W/m2)
    real(kind=r4),intent(in) :: rc2                  !Canopy resistence (s/m)
    
    real(kind=r4) :: evap                             !Evapotranspiration (mm/day)
    !     Parameters
    !     ----------
    real(kind=r4) :: ra, h5, t1, t2, es, es1, es2, delta_e, delta
    real(kind=r4) :: gama, gama2
    
    
    ra = rcmin
    h5 = 0.0275               !mb-1
    
    !     Delta
    !     -----
    t1 = temp + 1.
    t2 = temp - 1.
    es1 = tetens(t1)       !Saturation partial pressure of water vapour at temperature T
    es2 = tetens(t2)
    
    delta = (es1-es2)/(t1-t2) !mbar/oC
    !     
    !     Delta_e
    !     -------
    es = tetens (temp)
    delta_e = es*(1. - ur)    !mbar
    
    if ((delta_e.ge.(1./h5)-0.5).or.(rc2.ge.rcmax)) evap = 0.
    if ((delta_e.lt.(1./h5)-0.5).or.(rc2.lt.rcmax)) then
       !     Gama and gama2
       !     --------------
       gama  = spre*(1004.)/(2.45e6*0.622)
       gama2 = gama*(ra + rc2)/ra
       
       !     Real evapotranspiration
       !     -----------------------     
       evap = (delta* rn + (1.20*1004./ra)*delta_e)/(delta+gama2) !W/m2
       evap = evap*(86400./2.45e6) !mm/day
       evap = amax1(evap,0.)  !Eliminates condensation
    endif
  end function penman
  
  !=================================================================
  !=================================================================

  function available_energy(temp) result(ae)
    use global_pars, only: r4
    implicit none
    
    real(kind=r4),intent(in) :: temp
    real(kind=r4) :: ae
    
    ae = 2.895 * temp + 52.326 !from NCEP-NCAR Reanalysis data  
  end function  available_energy

  !=================================================================
  !=================================================================

  
  function runoff(wa) result(roff)
    use global_pars, only: r4
    implicit none
    
    real(kind=r4),intent(in) :: wa
    real(kind=r4):: roff
    
    !  roff = 38.*((w/wmax)**11.) ! [Eq. 10]
    roff = 11.5*((wa)**6.6) !from NCEP-NCAR Reanalysis data  
  end function  runoff
  
  !=================================================================
  !=================================================================
  
  function evpot2 (spre,temp,ur,rn) result(evap) 
    use global_pars, only: r4, rcmin, rcmax
    use photo, only: tetens
    implicit none
    
    !     Inputs
    
    real(kind=r4),intent(in) :: spre                 !Surface pressure (mb)
    real(kind=r4),intent(in) :: temp                 !Temperature (oC)
    real(kind=r4),intent(in) :: ur                   !Relative humidity (0-1,dimensionless)
    real(kind=r4),intent(in) :: rn                   !Radiation balance (W/m2)
    !     Output
    !     ------
    !     
    real(kind=r4) :: evap                 !Evapotranspiration (mm/day)
    !     Parameters
    !     ----------
    real(kind=r4) :: ra, t1, t2, es, es1, es2, delta_e, delta
    real(kind=r4) :: gama, gama2, rc
    
    ra = rcmin            !s/m
    
    !     Delta
    !     -----     
    t1 = temp + 1.
    t2 = temp - 1.
    es1 = tetens(t1)
    es2 = tetens(t2)
    delta = (es1-es2)/(t1-t2) !mb/oC
    !     
    !     Delta_e
    !     -------
    !     
    es = tetens (temp)
    delta_e = es*(1. - ur)    !mb
    !     
    !     Stomatal Conductance
    !     --------------------
    !     
    rc = rcmin
    !     
    !     Gama and gama2
  !     --------------
    !     
    gama  = spre*(1004.)/(2.45e6*0.622)
    gama2 = gama*(ra + rc)/ra
    !     
    !     Potencial evapotranspiration (without stress)
    !     ---------------------------------------------
    !     
    evap =(delta*rn + (1.20*1004./ra)*delta_e)/(delta+gama2) !W/m2
    evap = evap*(86400./2.45e6) !mm/day
    evap = amax1(evap,0.)     !Eliminates condensation
  end function evpot2
  
  !=================================================================
  !=================================================================
  
end module water
