module budget_caete
  implicit none
  private
  
  public :: budget
  
contains
  
  subroutine budget (month,w1,g1,s1,ts,temp,prec,p0,ipar,rh&
       &,cl1_pft,ca1_pft,cf1_pft,w2,g2,s2,smavg,ruavg,evavg,epavg&
       &,phavg,aravg,nppavg,laiavg,clavg,csavg,hravg,rcavg,rmavg,rgavg&
       &,cleafavg_pft,cawoodavg_pft,cfrootavg_pft,ocpavg)
    
    use global_pars
    use photo
    use water
    use productivity
    
    implicit none
    integer(kind=i4),parameter :: npft = npls
    
    
    !     ----------------------------INPUTS-------------------------------
    !        
    integer(kind=i4),intent(in) :: month             !Actual month (1-12)
    
    real(kind=r4),intent(in) :: w1(npft)             !Initial (previous month last day) soil moisture storage (mm)
    real(kind=r4),intent(in) :: g1(npft)             !Initial soil ice storage (mm)
    real(kind=r4),intent(in) :: s1(npft)             !Initial overland snow storage (mm)
    real(kind=r4),intent(inout) :: cl1_pft(npft)        ! initial BIOMASS cleaf compartment
    real(kind=r4),intent(inout) :: cf1_pft(npft)        !                 froot
    real(kind=r4),intent(inout) :: ca1_pft(npft)        !                 cawood
    
    real(kind=r4),intent(in) :: ts                   !Soil temperature (oC)
    real(kind=r4),intent(in) :: temp                 !Surface air temperature (oC)
    real(kind=r4),intent(in) :: prec                 !Precipitation (mm/day)
    real(kind=r4),intent(in) :: p0                   !Surface pressure (mb)
    real(kind=r4),intent(in) :: ae                   !Available energy (W/m2)
    real(kind=r4),intent(in) :: ca                   !Atmospheric carbon
    real(kind=r4),intent(in) :: ipar                 !Incident photosynthetic active radiation
    real(kind=r4),intent(in) :: rh                   !Relative humidity
    
    !     ----------------------------OUTPUTS------------------------------
    real(kind=r4),intent(out) :: epavg                !Maximum evapotranspiration monthly average (mm/day)
    real(kind=r4),intent(out) :: w2(npft)             !Final (last day) soil moisture storage (mm)
    real(kind=r4),intent(out) :: g2(npft)             !Final soil ice storage (mm)
    real(kind=r4),intent(out) :: s2(npft)             !Final overland snow storage (mm)
    real(kind=r4),intent(out) :: smavg(npft)          !Snowmelt monthly average (mm/day)
    real(kind=r4),intent(out) :: ruavg(npft)          !Runoff monthly average (mm/day)
    real(kind=r4),intent(out) :: evavg(npft)          !Actual evapotranspiration monthly average (mm/day)
    real(kind=r4),intent(out) :: phavg(npft)          !Monthly photosynthesis
    real(kind=r4),intent(out) :: aravg(npft)          !Monthly autotrophic respiration
    real(kind=r4),intent(out) :: nppavg(npft)         !Monthly NPP (average between PFTs)
    real(kind=r4),intent(out) :: laiavg(npft)         !Monthly leaf area Index
    real(kind=r4),intent(out) :: clavg(npft)          !Monthly carbon litter
    real(kind=r4),intent(out) :: csavg(npft)          !Monthly carbon soil
    real(kind=r4),intent(out) :: hravg(npft)          !Monthly heterotrophic respiration
    real(kind=r4),intent(out) :: rcavg(npft)          !Monthly canopy resistence
    real(kind=r4),intent(out) :: rmavg(npft),rgavg(npft) ! maintenance/growth respiration
    real(kind=r4),intent(out) :: cleafavg_pft(npft) ! Carbon in plant tissues
    real(kind=r4),intent(out) :: cawoodavg_pft(npft) 
    real(kind=r4),intent(out) :: cfrootavg_pft(npft)
    real(kind=r4),intent(out) :: ocpavg(npft)
    
    !     -----------------------Internal Variables------------------------
    integer(kind=i4) :: p ,i
    
    real(kind=r4) :: alfa_leaf(npft), alfa_awood(npft), alfa_froot(npft)
    real(kind=r4) :: beta_leaf(npft), beta_awood(npft), beta_froot(npft)
    
    !     RELATED WITH GRIDCELL OCUPATION
    
    real(kind=r4 ) :: ocp_coeffs(npft),ocp_mm(npft)
    logical :: ocp_wood(npft)
    
    !     WBM COMMUNICATION (water balance)
    real(kind=r4 ) :: tsnow                !Temperature threshold for snowfall (oC)
    real(kind=r4 ) :: tice                 !Temperature threshold for soil freezing (oC)
    real(kind=r4 ) :: psnow                !Snowfall (mm/day)
    real(kind=r4 ) :: prain                !Rainfall (mm/day)
    real(kind=r4 ) :: rimelt(npft)               !Runoff due to soil ice melting
    real(kind=r4 ) :: smelt(npft)                !Snowmelt (mm/day)
    real(kind=r4 ) :: w   (npft)           !Daily soil moisture storage (mm)
    real(kind=r4 ) :: g   (npft)           !Daily soil ice storage (mm)
    real(kind=r4 ) :: s   (npft)           !Daily overland snow storage (mm)
    real(kind=r4 ) :: ds  (npft)
    real(kind=r4 ) :: dw  (npft)
    real(kind=r4 ) :: roff(npft)           !Total runoff
    real(kind=r4 ) :: evap(npft)           !Actual evapotranspiration (mm/day)
    real(kind=r4 ) :: emax
        
    !c     Carbon Cycle
    real(kind=r4 ) ::  ph  (npft)           !Canopy gross photosynthesis (kgC/m2/yr)
    real(kind=r4 ) ::  ar  (npft)           !Autotrophic respiration (kgC/m2/yr)
    real(kind=r4 ) ::  nppa(npft)           !Net primary productivity / auxiliar
    real(kind=r4 ) ::  laia(npft)           !Leaf area index (m2 leaf/m2 area)
    real(kind=r4 ) ::  cl  (npft)           !Litter carbon (kgC/m2)
    real(kind=r4 ) ::  cs  (npft)           !Soil carbon (kgC/m2) 
    real(kind=r4 ) ::  hr  (npft)           !Heterotrophic (microbial) respiration (kgC/m2/yr)
    real(kind=r4 ) ::  rc2 (npft)           !Canopy resistence (s/m)
    real(kind=r8 ) ::  f1  (npft)           !
    real(kind=r4 ) ::  f5  (npft)           !Photosynthesis (mol/m2/s)
    real(kind=r4 ) ::  vpd (npft)           !Vapor Pressure deficit
    real(kind=r4 ) ::  rm(npft)             ! maintenance & growth a.resp
    real(kind=r4 ) ::  rg(npft)
    real(kind=r4 ) ::  cl1(npft),cf1(npft),ca1(npft) ! carbon pre-allocation 
    real(kind=r4 ) ::  cl2(npft),cf2(npft),ca2(npft) ! carbon pos-allocation
    real(kind=r4 ) ::  dl(npft)      ! litter resulting from pft exclusion
    


    do p = 1,npft
       rc2(p) = 0.0
       f1(p)  = 0.0
       f5(p) = 0.0
    enddo

    tsnow = -1.0
    tice  = -2.5
    
    !     Precipitation
    !     =============     
    psnow = 0.0
    prain = 0.0
    if (temp.lt.tsnow) then
       psnow = prec/real(ndmonth(month))
    else
       prain = prec/real(ndmonth(month))
    endif
    
    !     Initialization
    !     --------------
    epavg = 0.0
    do p = 1,npft
       w(p)       = w1(p)     ! hidrological pools state vars  
       g(p)       = g1(p)
       s(p)       = s1(p)
       smavg(p)   = 0.0       ! accumulators for i days of month m for pft p
       ruavg(p)   = 0.0
       evavg(p)   = 0.0
       rcavg(p)   = 0.0
       laiavg(p)  = 0.0
       phavg(p)   = 0.0
       aravg(p)   = 0.0
       nppavg(p)  = 0.0
       clavg(p)   = 0.0
       csavg(p)   = 0.0
       hravg(p)   = 0.0
       rmavg(p)   = 0.0
       rgavg(p)   = 0.0
       ocpavg(p)  = 0.0
       ocp_mm(p)  = 0.0
       
       alfa_leaf(p) = 0.0
       alfa_awood(p) = 0.0
       alfa_froot(p) = 0.0
       
    enddo
    
    !     Numerical integration
    !     ---------------------
    do i=1,ndmonth(month)
       emax  = 0.0
       
       cl1 = cl1_pft
       ca1 = ca1_pft
       cf1 = cf1_pft
       
       beta_leaf = alfa_leaf
       beta_awood = alfa_awood
       beta_froot = alfa_froot
       
       nppa  = 0.0
       ph    = 0.0
       ar    = 0.0
       laia  = 0.0
       f5    = 0.0
       f1    = 0.0
       vpd   = 0.0
       rc2   = 0.0
       rm    = 0.0
       rg    = 0.0
       dl    = 0.0        
       if ((i.eq.1).and.(month.eq.1)) then    
          beta_leaf = 0.00000001
          beta_awood = 0.00000001
          beta_froot = 0.00000001
       endif
       
       !     Grid cell area fraction (%) ocp_coeffs(pft(1), pft(2), ...,pft(p))
       !     =================================================================     
       CALL PFT_AREA_FRAC(CL1, CF1, CA1, OCP_COEFFS, OCP_WOOD) ! def in productivity1.f
       
       !     Maximum evapotranspiration   (emax)
       !     =================================
       emax = evpot2(p0,temp,rh,available_energy(temp))
       
       !     Productivity (ph, aresp, vpd, rc2 & etc.) for each PFT
       !     =================================
       do p = 1,npft
          
          call prod(p,ocp_coeffs(p),OCP_WOOD(P),temp,p0,w(p)&
               &,ipar,rh,emax,cl1(p),ca1(p),cf1(p),beta_leaf(p)&
               &,beta_awood(p),beta_froot(p),ph(p),ar(p),nppa(p)&
               &,laia(p),f5(p),f1(p),vpd(p),rm(p),rg(p),rc2(p))
          
          
          !c     Carbon allocation (carbon content on each compartment)
          !     =====================================================
          call allocation (p, nppa(p), cl1(p), ca1(p),cf1(p),cl2(p),&
               & ca2(p), cf2(p), dl(p)) 
          
          
          alfa_leaf(p)  = cl2(p) - cl1(p)
          alfa_awood(p) = ca2(p) - ca1(p)
          alfa_froot(p) = cf2(p) - cf1(p)
          
          !     Snow budget
          !     ===========     
          smelt(p) = 2.63 + 2.55*temp + 0.0912*temp*prain !Snowmelt (mm/day)
          smelt(p) = amax1(smelt(p),0.)
          smelt(p) = amin1(smelt(p),s(p)+psnow)
          ds(p) = psnow - smelt(p)
          s(p) = s(p) + ds(p)
          
          !     Water budget
          !     ============
          if (ts.le.tice) then !Frozen soil
             g(p) = g(p) + w(p) !Soil moisture freezes
             w(p) = 0.0
             roff(p) = smelt(p) + prain !mm/day
             evap(p) = 0.0
             ph(p) = 0.0
             ar(p) = 0.0
             nppa(p) = 0.0
             laia(p) = 0.0
             cl(p) = 0.0
             cs(p) = 0.0
             hr(p) = 0.0
          else                !Non-frozen soil
             w(p) = w(p) + g(p)
             g(p) = 0.0
             rimelt(p) = 0.0
             if (w(p).gt.wmax) then
                rimelt(p) = w(p) - wmax !Runoff due to soil ice melting
                w(p) = wmax
             endif
             
             roff(p) = runoff(w(p)/wmax)       !Soil moisture runoff (roff, mm/day)
             
             evap(p) = penman(p0,temp,rh,ae,rc2(p)) !Actual evapotranspiration (evap, mm/day)
             dw(p) = prain + smelt(p) - evap(p) - roff(p)
             w(p) = w(p) + dw(p)
             if (w(p).gt.wmax) then
                roff(p) = roff(p) + (w(p) - wmax)
                w(p) = wmax
             endif
             if (w(p).lt.0.) w(p) = 0.
             roff(p) = roff(p) + rimelt(p) !Total runoff
             
             
             !     Carbon cycle (Microbial respiration, litter and soil carbon)
             !     ============================================================     
             call carbon2 (ts,f5(p),evap(p),laia(p),dl(p),cl(p),cs(p),hr(p))   
          endif
          
          !     Accumulate daily budgets weighted by occupation coefficients
          
          ocp_mm(p) = ocp_mm(p) + ocp_coeffs(p)
          
          if(p .eq. 1) epavg = epavg + emax !mm/day
          smavg(p) = smavg(p) + smelt(p)
          ruavg(p) = ruavg(p) + roff(p)  ! mm day-1
          evavg(p) = evavg(p) + evap(p) * ocp_coeffs(p)  ! mm day-1
          rcavg(p) = rcavg(p) + rc2(p) * ocp_coeffs(p) ! s m -1
          
          phavg(p) = phavg(p) +   ph(p) !kgC/m2/day
          aravg(p) = aravg(p) +   ar(p)  !kgC/m2/year
          nppavg(p) = nppavg(p) + nppa(p) !kgC/m2/day
          
          laiavg(p) = laiavg(p) + laia(p)
          clavg(p) = clavg(p) + cl(p) * ocp_coeffs(p) !kgC/m2/day
          csavg(p) = csavg(p) + cs(p) * ocp_coeffs(p) !kgC/m2/day
          hravg(p) = hravg(p) + hr(p) * ocp_coeffs(p) !kgC/m2/day
          rmavg(p) = rmavg(p) + rm(p) * ocp_coeffs(p)
          rgavg(p) = rgavg(p) + rg(p) * ocp_coeffs(p)
          cleafavg_pft(p)  =  cl2(p)
          cawoodavg_pft(p) =  ca2(p)
          cfrootavg_pft(p) =  cf2(p)
          cl1_pft(p) = cl2(p)     ! Adicionado ------ para fazer transforcaoes diarias
          ca1_pft(p) = ca2(p)     ! Adicionado ------ para fazer transforcaoes diarias
          cf1_pft(p) = cf2(p)     ! Adicionado ------ para fazer transforcaoes diarias
       enddo
    enddo                     ! end ndmonth loop
    
    !     Final calculations
    !     ------------------
    !     monthly values
    do p=1,NPFT
       if (p .eq. 1) epavg = epavg/real(ndmonth(month))
       w2(p) = w(p) * (ocp_mm(p)/real(ndmonth(month)))
       g2(p) = g(p) * (ocp_mm(p)/real(ndmonth(month)))
       s2(p) = s(p) * (ocp_mm(p)/real(ndmonth(month)))
       smavg(p) = smavg(p)/real(ndmonth(month))
       ruavg(p) = ruavg(p)/real(ndmonth(month))
       evavg(p) = evavg(p)/real(ndmonth(month))
       rcavg(p) = rcavg(p)/real(ndmonth(month))
       phavg(p) = (phavg(p)/365.0) * 12.0 !kgC/m2/yr
       aravg(p) = (aravg(p)/365.0) * 12.0   !kgC/m2/yr
       nppavg(p) = (nppavg(p)/365.0) * 12.0 !kgC/m2/yr
       laiavg(p) = (laiavg(p)/365.0) * 12.0
       clavg(p) = (clavg(p)/365.0) * 12.0 !kgC/m2
       csavg(p) = (csavg(p)/365.0) * 12.0 !kgC/m2
       hravg(p) = (hravg(p)/365.0) * 12.0 !kgC/m2/yr
       rmavg(p) = (rmavg(p)/365.0) * 12.0 
       rgavg(p) = (rgavg(p)/365.0) * 12.0
       ocpavg(p) = ocp_coeffs(p) * 100.
    enddo
    return
  end subroutine budget
  
end module budget_caete
