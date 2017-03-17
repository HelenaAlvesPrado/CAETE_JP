module productivity
  implicit none
  private

  public :: prod


contains
  
  subroutine prod(pft,ocp_pft,light_limit,temp,p0,w,&
       ipar,rh,emax,cl1,ca1,cf1,beta_leaf,beta_awood,& ! inputs
       beta_froot,ph,ar,nppa,laia,f5,f1,vpd,rm,rg,rc) ! outputs
    
    use global_pars
    use photo
    use water
    
    implicit none  
    !=Variables/Parameters
    !=====================
    integer(kind=i4),parameter :: q=npls
    !Input
    !-----     
    
    integer(kind=i4), intent(in) :: pft
    
    real(kind=r4), intent(in) :: ocp_pft              !PFT area occupation (%)
    real(kind=r4), intent(in) :: temp                 !Mean monthly temperature (oC)
    real(kind=r4), intent(in) :: p0                   !Mean surface pressure (hPa)
    real(kind=r4), intent(in) :: w                    !Soil moisture (dimensionless)
    real(kind=r4), intent(in) :: ipar                 !Incident photosynthetic active radiation (w/m2)'
    real(kind=r4), intent(in) :: rh,emax              !Relative humidity/MAXIMUM EVAPOTRANSPIRATION
    real(kind=r4), intent(in) :: cl1, cf1, ca1        !Carbon in plant tissues (kg/m2)
    real(kind=r4), intent(in) :: beta_leaf            !npp allocation to carbon pools (kg/m2/day)
    real(kind=r4), intent(in) :: beta_awood
    real(kind=r4), intent(in) :: beta_froot
    
    logical, intent(in) :: light_limit                !True for no ligth limitation
    
    !     Output
    !     ------
    real(kind=r4), intent(out) :: ph                   !Canopy gross photosynthesis (kgC/m2/yr)
    real(kind=r4), intent(out) :: rc                   !Stomatal resistence (not scaled to canopy!) (s/m)
    real(kind=r4), intent(out) :: laia                 !Autotrophic respiration (kgC/m2/yr)
    real(kind=r4), intent(out) :: ar                   !Leaf area index (m2 leaf/m2 area)
    real(kind=r4), intent(out) :: nppa                 !Net primary productivity (kgC/m2/yr) 
    real(kind=r4), intent(out) :: vpd            
    real(kind=r4), intent(out) :: f5                   !Water stress response modifier (unitless) 
    real(kind=r4), intent(out) :: rm                   !autothrophic respiration (kgC/m2/day)
    real(kind=r4), intent(out) :: rg 
    
    !     Internal
    !     --------
    
    real(kind=r4) :: f1       !Leaf level gross photosynthesis (molCO2/m2/s)
    real(kind=r4) :: f1a      !auxiliar_f1
    
    real(kind=r4),dimension(q) :: tleaf             !leaf turnover time (yr)
    real(kind=r4),dimension(q) :: p21               !Maximum carboxilation rate (micromolC m-2 d-1)
    real(kind=r4),dimension(q) :: g1
    
    real(kind=r4) :: sla          !specific leaf area (m2/kg)
        
    !getting pft parameters
    

    call pft_par(2, p21)
    call pft_par(3, tleaf)
    call pft_par(1, g1)
    
    !     ==============
    !     Photosynthesis 
    !     ==============
    
    !Rubisco maximum carboxylaton rate (molCO2/m2/s)
    !-----------------------------------------------
    
    f1a = photosynthesis_rate(p21(pft),temp,p0,ipar,light_limit)
    
    ! VPD
    !========
    vpd = vapor_p_defcit(temp,rh)
    
    !Stomatal resistence
    !===================
    rc = canopy_resistence(vpd, f1a, g1(pft))
    
    !     Water stress response modifier (dimensionless)
    !     ----------------------------------------------
    f5 =  water_stress_modifier(w, cf1, rc, emax)
    
    !     Photosysthesis minimum and maximum temperature
    !     ----------------------------------------------
    
    if ((temp.ge.-10.0).and.(temp.le.50.0)) then
       f1 = f1a * f5 ! :water stress factoF
    else
       f1 = 0.0               !Temperature above/below photosynthesis windown
    endif
    
    
    !     Leaf area index (m2/m2)
    laia = leaf_area_index(cl1,spec_leaf_area(tleaf(pft)))
    
    !     Canopy gross photosynthesis (kgC/m2/yr)
    !     =======================================x
    ph =  gross_ph(f1,cl1,sla) * ocp_pft       ! kg m-2 year-1

    
    !     Autothrophic respiration
    !     ========================
    !     Maintenance respiration (kgC/m2/yr) (based in Ryan 1991)
    rm = m_resp(temp,cl1,cf1,ca1,ocp_pft)
  
    ! c     Growth respiration (KgC/m2/yr)(based in Ryan 1991; Sitch et al.
    ! c     2003; Levis et al. 2004)         
    rg = g_resp(beta_leaf,beta_awood, beta_froot, ocp_pft) 
    
    if (rg.lt.0) then
       rg = 0.0
    endif
    
    !     c Autotrophic (plant) respiration -ar- (kgC/m2/yr)
    !     Respiration minimum and maximum temperature
    !     -------------------------------------------     
    if ((temp.ge.-10.0).and.(temp.le.50.0)) then
       ar = rm + rg
    else
       ar = 0.0               !Temperature above/below respiration windown
    endif
    
    !c     -----------------------------------------------------------------
    !     NPP
    !     ============
    !     Productivity
    !     ============
    !     Net primary productivity(kgC/m2/yr)
    !     ====================================
    nppa = ph - ar
    
    if(nppa .lt. 0.0) nppa = 0.0
    
  end subroutine prod
  
end module productivity
