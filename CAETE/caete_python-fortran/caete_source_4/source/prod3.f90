!=================================================================
!code based in CPTEC-PVM2
!LAST UPDATE: 28 Jan 2017 00:37:49 BRST JP 
! =================================================================

subroutine wbm (prec,temp,p0,ca,par,rhs,cleaf_ini,cawood_ini&
     &,cfroot_ini,emaxm, tsoil, photo_pft,aresp_pft,npp_pft,lai_pft&
     &,clit_pft,csoil_pft, hresp_pft,rcm_pft,runom_pft,evapm_pft&
     &,wsoil_pft,rm_pft,rg_pft,cleaf_pft,cawood_pft,cfroot_pft,grid_area)
  
  use global_pars
  implicit none
  integer(kind=i4),parameter :: q = npls ! plss futuramente serao
  integer(kind=i4),parameter :: nt = ntimes ! (meses) ! modelo estac. --

  
  !c     --------------------------I N P U T S----------------------------
  real(kind=r4),intent(in) :: ca                   !CO2 atmospheric concentration (ppmv)
  real(kind=r4),intent(in) :: p0(nt)         !Atmospheric pressure (mb)
  real(kind=r4),intent(in) :: prec(nt)       !Precipitation (mm/month)
  real(kind=r4),intent(in) :: temp(nt)       !Temperature (oC)
  real(kind=r4),intent(in) :: par(nt)        !IPAR (Ein/m2/s)
  real(kind=r4),intent(in) :: rhs(nt)        !Relative humidity
  
  real(kind=r4),intent(in) ::  cleaf_ini(q)  ! Initial carbon content in leaves (kg m-2)
  real(kind=r4),intent(in) :: cawood_ini(q) ! Initial carbon content in aboveground wood (kg m-2)
  real(kind=r4),intent(in) :: cfroot_ini(q) ! Initial carbon content in fineroots (kg m-2)
  !C     -----------------------------E N D-------------------------------
      
  !c     -------------------------O U T P U T S---------------------------
       
  real(kind=r4),intent(out) :: tsoil(nt)     !soil temperature
  
  real(kind=r4),intent(out) :: photo_pft(nt,q) !Monthly photosynthesis   (kgC m-2)
  real(kind=r4),intent(out) :: aresp_pft(nt,q) !Monthly autotrophic res  (kgC m-2)
  real(kind=r4),intent(out) :: npp_pft  (nt,q) !Monthly net primary produ (kgC m-2)
      
  real(kind=r4),intent(out) :: lai_pft  (nt,q) !Monthly leaf area index
  real(kind=r4),intent(out) :: clit_pft (nt,q) !Monthly litter carbon
  real(kind=r4),intent(out) :: csoil_pft(nt,q) !Monthly soil carbon
  real(kind=r4),intent(out) :: hresp_pft(nt,q) !Monthly het resp  (kgC/m2)
  real(kind=r4),intent(out) :: rcm_pft  (nt,q) 

  real(kind=r4),intent(out) :: emaxm    (nt) !Max.evapotranspiration (kg m-2 day-1)
  real(kind=r4),intent(out) :: runom_pft(nt,q) !Runoff 
  real(kind=r4),intent(out) :: evapm_pft(nt,q) !Actual evapotranspiration        
  real(kind=r4),intent(out) :: wsoil_pft(nt,q) !Soil moisture (mm)
  
  !c     NOVOS OUTPUTS DA BIANCA
  real(kind=r4),intent(out) :: rm_pft   (nt,q)
  real(kind=r4),intent(out) :: rg_pft   (nt,q)
  real(kind=r4),intent(out) :: cleaf_pft (q) ! leaf biomass (KgC/m2)
  real(kind=r4),intent(out) :: cawood_pft(q) ! aboveground wood biomass (KgC/m2)
  real(kind=r4),intent(out) :: cfroot_pft(q)  ! fine root biomass      
  real(kind=r4),intent(out) :: grid_area(q)   ! gridcell area fraction of pfts!
  
  !  c     --------------------------------E N D----------------------------
  
  !  c     ------------------------- internal variables---------------------
  
  integer(kind=i4) :: k, kk, n, mes, nerro, p
      
  real(kind=r4) :: wg0  (nt)      !Previous year soil moisture
  real(kind=r4) :: wsoilt(nt)     !soil water (check wbm equilibrium)
  real(kind=r4) :: gsoilt(nt)     !soil ice   (check wbm equilibrium)
  real(kind=r4) :: gsoil(nt,q)    !Soil ice 
  real(kind=r4) :: ssoil(nt,q)    !Soil snow
  real(kind=r4) :: snowm(nt,q)    !Snowmelt

  real(kind=r4) :: sini(q), sfim(q)
  real(kind=r4) :: wini(q), wfim(q)
  real(kind=r4) :: gini(q), gfim(q)

  real(kind=r4) :: cleaf1_pft (q)
  real(kind=r4) :: cawood1_pft(q)
  real(kind=r4) :: cfroot1_pft(q)
  real(kind=r4) :: wood(q)
!     outputs for budget (these variables store monthly values from budget)
  real(kind=r4) :: epmes ! equal for all pfts - potential evapotranspiration(PET)(mm day-1)
  real(kind=r4) :: rmes(q),phmes(q),smes(q),rcmes(q),hrmes(q)
  real(kind=r4) :: nppmes(q),laimes(q), armes(q),clmes(q),csmes(q)
  real(kind=r4) :: emes(q),rmmes(q),rgmes(q)
  real(kind=r4) :: cleafmes(q),cawoodmes(q),cfrootmes(q), gridocpmes(q) 
  real(kind=r4) :: wsaux1,dwww,wmax     !auxiliar to equilibrium check
  real(kind=r4) :: ae                   !Available energy     
  real(kind=r4) :: pr,spre,ta,td,ipar,ru

  real(kind=r4) :: leaf0(q), froot0(q),awood0(q)
  real(kind=r4) :: biomass, biomass0
  real(kind=r4) :: sensi
  logical :: check
    
  
  !     ================      
  !     Soil temperature
  !     ================
  call soil_temp(temp,tsoil)
!     ============
!     Water budget
!     ============
!     
!     For all grid points
!     -------------------      

  do p = 1,q   
     cleaf1_pft (p) =  cleaf_ini(p)
     cawood1_pft(p) = cawood_ini(p)
     cfroot1_pft(p) = cfroot_ini(p)
     
     leaf0(p) = cleaf_ini(p)
     froot0(p) = cfroot_ini(p)
     awood0(p) = cawood_ini(p)
     
     cleaf_pft(p)  = 0.0 ! leaf biomass (KgC/m2)
     cawood_pft(p) = 0.0 ! aboveground biomass (KgC/m2)
     cfroot_pft(p) = 0.0 ! fine root biomass (KgC/m2)
     grid_area(p)  = 0.0 ! gridcell area fraction of pfts(%)
  enddo
  
  ! c     Write to track program execution     
  !            if ((mod(j,ny).eq.0).and.(mod(i,10).eq.0))
  !     &          write(*,*) 'working: ', (real(i)/real(nx))*100.0, '%' 
          
!     Initialize variables
  do k=1,nt
     wg0  (k) = -1.0 !Soil water content in preceeding year integration
     emaxm(k) =  0.0 !Maximum evapotranspiration
     do p=1,q
        photo_pft(k,p) = 0.0 !Monthly photosynthesis (kgC/m2)
        aresp_pft(k,p) = 0.0 !Monthly autotrophic respiration (kgC/m2)
        npp_pft  (k,p) = 0.0 !Monthly net primary productivity (average between PFTs) (kgC/m2)
        lai_pft  (k,p) = 0.0 !Monthly leaf area index
        clit_pft (k,p) = 0.0 !Monthly litter carbon
        csoil_pft(k,p) = 0.0 !Monthly soil carbon
        hresp_pft(k,p) = 0.0 !Monthly heterotrophic respiration (kgC/m2)
        rcm_pft  (k,p) = 0.0
        gsoil    (k,p) = 0.0 !Soil ice
        ssoil    (k,p) = 0.0 !Soil snow
        runom_pft(k,p) = 0.0 !Runoff
        evapm_pft(k,p) = 0.0 !Actual evapotranspiration        
        wsoil_pft(k,p) = 0.0 !Soil moisture (mm)
        rm_pft   (k,p) = 0.0
        rg_pft   (k,p) = 0.0 
     enddo
  enddo
  !     Set some variables
  !     ------------------
  do p = 1,q
     wini(p)  = 0.01  !Soil moisture_initial condition (mm)
     gini(p)  = 0.0   !Soil ice_initial condition (mm)
     sini(p)  = 0.0   !Overland snow_initial condition (mm)
  enddo

  !     =================
  !     START INTEGRATION
  !     =================
  n = 0
10 continue
  n = n + 1
  
  k = mod(n,12)
  if (k.eq.0) k = 12
  mes = k
  spre = p0(k) * 0.01 ! transforamando de Pascal pra mbar (hPa)
  td = tsoil(k)
  ta = temp(k)
  pr = prec(k)
  ipar = par(k)
  ru = rhs(k)
  ae = 2.895*ta+52.326 !Available energy (W/m2) - From NCEP-NCAR reanalysis data
!     
!     Monthly water budget
!     ====================
  
  call budget (mes,wini,gini,sini,td,ta,pr,spre,ae,ca,ipar,ru&
       &,cleaf1_pft,cawood1_pft,cfroot1_pft ,wfim,gfim, sfim,smes&
       &,rmes,emes,epmes,phmes,armes,nppmes,laimes,clmes,csmes,hrmes&
       &,rcmes,rmmes,rgmes,cleafmes,cawoodmes,cfrootmes, gridocpmes)

  do p=1,q
     if(p .eq. 1) emaxm(k) = epmes
     gsoil    (k,p) = gfim(p)
     ssoil    (k,p) = sfim(p)
     wsoil_pft(k,p) = wfim(p)
     snowm    (k,p) = smes(p)
     runom_pft(k,p) = rmes(p)
     evapm_pft(k,p) = emes(p)
     rcm_pft  (k,p) = rcmes(p)
     lai_pft  (k,p) = laimes(p)
     photo_pft(k,p) = phmes(p)
     aresp_pft(k,p) = armes(p)
     npp_pft  (k,p) = nppmes(p)
     clit_pft (k,p) = clmes(p)
     csoil_pft(k,p) = csmes(p)
     hresp_pft(k,p) = hrmes(p)
     rm_pft   (k,p) = rmmes(p)
     rg_pft   (k,p) = rgmes(p)
                  
     wini(p) = wfim(p)
     gini(p) = gfim(p)
     sini(p) = sfim(p)
     
     cleaf1_pft(p)  = cleafmes(p) 
     cawood1_pft(p) = cawoodmes(p)
     cfroot1_pft(p) = cfrootmes(p)
     
     if(k .eq. 12) then
        cleaf_pft (p) = cleafmes(p)
        cawood_pft(p) = cawoodmes(p)
        cfroot_pft(p) = cfrootmes(p)
        grid_area (p) = gridocpmes(p)
     endif
  enddo
               
!     Check if equilibrium is attained
!     --------------------------------
  if (k.eq.12) then
     wsoilt(k) = 0.0
     gsoilt(k) = 0.0

     do p = 1,q
        wsoilt(k) = wsoilt(k) + wsoil_pft(k,p)
        gsoilt(k) = gsoilt(k) + gsoil(k,p)
     enddo
                  
     wmax = 500.
     nerro = 0
     
     do kk=1,nt
        wsaux1 = wsoilt(kk) + gsoilt(kk)   
        dwww = (wsaux1 - wg0(kk)) / wmax
        if (abs(dwww).gt.0.01) nerro = nerro + 1
     enddo
                  
     if (nerro.ne.0) then
        
        do kk=1,12
           wg0(kk) = wsoilt(kk) + gsoilt(kk)
        enddo
     else
        goto 100
     endif
  endif
  goto 10
               
  !     PFTs equilibrium check
  !     ==================================================
  !     tentativa 1 - usando variacao no pool de C vegetal
  !     --------------------------------------------------
  if (k.eq.12) then
     nerro = 0
     biomass = 0.0
     biomass0 = 0.0
     check = .false.
     sensi = 1.0! (kg/m2/y) if biomas change .le. sensi: equilibrium
     ! brienen et al. 2015 mean biomass change in Amazon forest - 1995 
     call pft_par(4,wood)
     
     do p = 1,q
        if(cleaf_pft(p) .gt. 0.0 .and.cfroot_pft(p).gt. 0.0)&
             & then
           ! GRASS
           if(wood(p) .le. 0.0) then
              check = .true.
              biomass = biomass + cleaf1_pft(p) + cfroot1_pft(p)  
              biomass0 = biomass0 + leaf0(p) + froot0(p)
              ! WOOD
           else
              check = .true.
              biomass = biomass + cleaf1_pft(p) + cfroot1_pft(p) +&
                   & cawood1_pft(p)
              biomass0 = biomass0 + leaf0(p) + froot0(p) + awood0(p)
           endif
        endif
     enddo
     if(check) then
        if (abs(biomass0-biomass) .gt. sensi) then
           do p=1,q
              leaf0(p) = cleaf1_pft(p)
              froot0(p) = cfroot1_pft(p)
              awood0(p) = cawood1_pft(p)
              nerro = nerro + 1
           enddo
        endif
     endif
     if(nerro .gt. 0) then
        !print*, abs(biomass-biomass0), n
        goto 10
     else
        !print*, 'eq attained',n
        continue
     endif
  endif
100 continue !mudar para linha 245 para incluir pft check
  return
  
end subroutine wbm




subroutine budget (month,w1,g1,s1,ts,temp,prec,p0,ae,ca,ipar,rh&
     &,cl1_pft,ca1_pft,cf1_pft,w2,g2,s2,smavg,ruavg,evavg,epavg&
     &,phavg,aravg,nppavg,laiavg,clavg,csavg,hravg,rcavg,rmavg,rgavg&
     &,cleafavg_pft,cawoodavg_pft,cfrootavg_pft,ocpavg)

  use global_pars
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
  real(kind=r4 ) :: wmax                 !Soil moisture availability (mm)
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
  real(kind=r8 ) :: wapft                     !Maximum evapotranspiration
      
  integer(kind=i4) :: ndmonth(12)       !Number of months
  data ndmonth /31,28,31,30,31,30,31,31,30,31,30,31/ !Number of days for each month
  
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
      

  !     Initialize Parameters
  !     ---------------------------------------
  do p = 1,npft
     rc2(p) = 0.0
     f1(p)  = 0.0
     f5(p) = 0.0
  enddo
  
  !    more Parameters
  !     ----------
  wmax  = 500.0
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
     do p=1,npft
        cl1(p) = cl1_pft(p)
        ca1(p) = ca1_pft(p)
        cf1(p) = cf1_pft(p)
        
        beta_leaf(p) = alfa_leaf(p)
        beta_awood(p) = alfa_awood(p)
        beta_froot(p) = alfa_froot(p)
        
        nppa(p)  = 0.0
        ph(p)    = 0.0
        ar(p)    = 0.0
        laia(p)  = 0.0
        f5(p)    = 0.0
        f1(p)    = 0.0
        vpd(p)   = 0.0
        rc2(p)   = 0.0
        rm(p)    = 0.0
        rg(p)    = 0.0
        
        if ((i.eq.1).and.(month.eq.1)) then    
           beta_leaf(p) = 0.00000001
           beta_awood(p) = 0.00000001
           beta_froot(p)= 0.00000001
        endif
     enddo
         
     !     Grid cell area fraction (%) ocp_coeffs(pft(1), pft(2), ...,pft(p))
     !     =================================================================     
     CALL PFT_AREA_FRAC(CL1, CF1, CA1, OCP_COEFFS, OCP_WOOD) ! def in productivity1.f
         
     !     Maximum evapotranspiration   (emax)
     !     ================================= 
     call evpot2 (p0,temp,rh,ae,emax)
             
     !     Productivity (ph, aresp, vpd, rc2 & etc.) for each PFT
     !     =================================
     do p = 1,npft
        
        call productivity1 (p,ocp_coeffs(p),OCP_WOOD(P),temp,p0,w(p)&
             &,wmax,ca,ipar,rh,emax,cl1(p),ca1(p),cf1(p),beta_leaf(p)&
             &,beta_awood(p),beta_froot(p),ph(p),ar(p),nppa(p)&
             &,laia(p),f5(p),f1(p),vpd(p),rm(p),rg(p),rc2(p))

            
        !c     Carbon allocation (carbon content on each compartment)
        !     =====================================================
        call allocation (p, nppa(p), cl1(p), ca1(p),cf1(p),cl2(p),&
             & ca2(p), cf2(p)) 

            
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
           
           wapft = (w(p)/wmax)
           call runoff (wapft,roff(p))       !Soil moisture runoff (roff, mm/day)
           
           call penman (p0,temp,rh,ae,rc2(p),evap(p)) !Actual evapotranspiration (evap, mm/day)
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
           call carbon2 (ts,f5(p),evap(p),laia(p),cl(p),cs(p),hr(p))   
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




subroutine soil_temp(temp, tsoil)
  ! Calcula a temperatura do solo. Aqui vamos mudar no futuro!
  ! a tsoil deve ter relacao com a et realizada...
  ! a profundidade do solo (H) e o coef de difusao (DIFFU) devem ser
  ! variaveis (MAPA DE SOLO?; agua no solo?)
  use global_pars
  implicit none
  integer(kind=i4),parameter :: m = ntimes ! (meses) ! modelo estac. --

  !parameters
  real(kind=r4), parameter :: H = 1.0                         ! soil layer thickness (meters)
  real(kind=r4), parameter :: DIFFU = 4.e7 * (30.0 * 86400.0) ! soil thermal diffusivity (m2/mes)
  real(kind=r4), parameter :: TAU = (H ** 2) / (2.0 * DIFFU)  ! e-folding times (months) 
  ! i/o

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
end subroutine soil_temp




subroutine productivity1 (pft,ocp_pft,ligth_limit,temp,p0,w,&
     wmax,ca,ipar,rh,emax,cl1,ca1,cf1,beta_leaf,beta_awood,& ! inputs
     beta_froot,ph,ar,nppa,laia,f5,f1,vpd,rm,rg,rc) ! outputs

  use global_pars
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
  real(kind=r4), intent(in) :: w,wmax            !Soil moisture (dimensionless)
  real(kind=r4), intent(in) :: ca                   !Atmospheric CO2 concentration (Pa)
  real(kind=r4), intent(in) :: ipar                 !Incident photosynthetic active radiation (w/m2)'
  real(kind=r4), intent(in) :: rh,emax                   !Relative humidity
  real(kind=r4), intent(in) :: cl1, cf1, ca1        !Carbon in plant tissues (kg/m2)
  real(kind=r4), intent(in) :: beta_leaf            !npp allocation to carbon pools (kg/m2/day)
  real(kind=r4), intent(in) :: beta_awood
  real(kind=r4), intent(in) :: beta_froot

  logical, intent(in) :: ligth_limit                !True for no ligth limitation
  
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
  real(kind=rbig) :: f5_64                           !f5 auxiliar(more precision)
  
  real(kind=r4) :: es, es2,wa                        !Saturation partial pressure (hPa) == mbar
  real(kind=r4) :: aux_ipar
  
  real(kind=r8) :: vm                                !Rubisco maximum carboxylaton rate (molCO2/m2/s)
  real(kind=r8) :: mgama                             !Photo-respiration compensation point (Pa)
  real(kind=r8) :: rmax                              !Saturated mixing ratio (kg/kg)
  real(kind=r8) :: r                                 !Moisture deficit at leaf level (kg/kg)
  real(kind=r8) :: ci                                !Internal leaf CO2 partial pressure (Pa)
  real(kind=r8) :: a,b,c,a2,b2,c2                    !Auxiliars
  real(kind=r8) :: delta,delta2                      !Auxiliars
  real(kind=r8) :: sunlai,shadelai                   !Sunlai/Shadelai
  real(kind=r8) :: vpd64, vpd_real                 !Vapor pressure deficit (kPa)
  real(kind=r8) :: laia64, ph64, ar64
  real(kind=r8) :: rm64, rms64, rml64
  real(kind=r8) :: rmf64, rg64, rgf64
  real(kind=r8) :: rgs64, rgl64
  real(kind=r8) :: leaf_t_months
  real(kind=r8) :: leaf_turnover
  real(kind=r8) :: leaf_t_coeff
  
  ! Rubisco, light and transport limited photosynthesis rate
  ! --------------------------------------------------------
  real(kind=r8) :: jc,jl,je,jp,jp1,jp2,j1,j2 !Auxiliars    
  real(kind=r8) :: f1       !Leaf level gross photosynthesis (molCO2/m2/s)
  real(kind=r8) :: f1a      !auxiliar_f1
  real(kind=r8) :: f2       !Michaelis-Menton CO2 constant (Pa)
  real(kind=r8) :: f3       !Michaelis-Menton O2 constant (Pa)
  real(kind=r8) :: f4,f4sun,f4shade !Scaling-up LAI to canopy level (dimensionless)
      
  real(kind=r4),dimension(q) :: tleaf             !leaf turnover time (yr)
  real(kind=r4),dimension(q) :: p21               !Maximum carboxilation rate (micromolC m-2 d-1)

  real(kind=r8) :: sla          !specific leaf area (m2/kg)
  real(kind=r8) :: csa      !sapwood compartment´s carbon content (5% of woody tissues) (kgC/m2)
  real(kind=r8) :: ncl      !leaf N:C ratio (gN/gC)
  real(kind=r8) :: ncf      !fine roots N:C ratio (gN/gC)
  real(kind=r8) :: ncs      !sapwood N:C ratio(gN/gC)
  real(kind=r8) :: csai
  
  real(kind=r8) :: pt,d     !taxa potencial de fornecimento para transpiração (mm/dia)
  real(kind=r8) :: csru     !Specific root water uptake (0.5 mm/gC/dia; based in jedi
  real(kind=r8) :: alfm     !maximum Priestley-Taylor coefficient (based in Gerten et al. 2004
  real(kind=r8) :: gm       !scaling conductance (dia/mm)(based in Gerten et al. 2004;
  real(kind=r8) :: gc       !Canopy conductance (dia/mm)(based in Gerten et al. 2004;

  ! SOME MODEL PARAMETERS
  real(kind=r8) ::  p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14
  real(kind=r8) ::  p15,p19,p20,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31


  ! setting some parameter values
  a   = 0.8300              !Photosynthesis co-limitation coefficient
  a2  = 0.930               !Photosynthesis co-limitation coefficient
  p3  = 21200.0             !Atmospheric oxygen concentration (Pa)
  p4  = 0.080               !Quantum efficiency (mol electrons/Ein)
  p5  = 0.150               !Light scattering rate
  p6  = 2.0                 !Parameter for jl
  p7  = 0.50                !Ratio of light limited photosynthesis to Rubisco carboxylation
  p8  = 5200.0              !Photo-respiration compensation point
  p9  = 0.570               !Photosynthesis co-limitation coefficient
  p10 = 0.100               !Q10 function
  p11 = 25.0                !Q10 function reference temperature (oC)
  p12 = 30.0                !Michaelis-Menten constant for CO2 (Pa)
  p13 = 2.100               !Michaelis-Menten constant for CO2
  p14 = 30000.0             !Michaelis-Menten constant for O2 (Pa)
  p15 = 1.20                !Michaelis-Menten constant for O2
  p19 = 0.90                !Maximum ratio of internal to external CO2
  p20 = 0.10                !Critical humidity deficit (kg/kg)
  p22 = 2.0                 !Rubisco carboxylation rate
  p23 = 0.30                !Rubisco carboxylation rate
  p24 = 36.0                !Rubisco carboxylation rate (oC)
  p25 = 0.000008            !Maximum gross photosynthesis rate (molCO2/m2/s)
  p26 = 0.50                !light extinction coefficient for IPAR/sun (0.5/sen90)
  p27 = 1.50                !light extinction coefficient for IPAR/shade (0.5/sen20)
  p28 = 0.500               !Soil moisture at field capacity
  p29 = 0.205               !Soil moisture at wilting point
  p30 = 0.015               !Ratio of respiration to Rubisco carboxylation rates
  p31 = 3.850               !Whole plant to leaf respiration ratio
  
  !     getting pft parameters
  call pft_par(2, p21)
  call pft_par(6, tleaf)

  !     ==============
  !     Photosynthesis 
  !     ==============
  
  !     Rubisco maximum carboxylaton rate (molCO2/m2/s)
  !     -----------------------------------------------
  vm = (p21(pft))
  !c      vm = (p21(pft)*(p22**(p10*(temp-p11))))/ !Free-range parameter --> 0.0358>vm>840 (micromol)
  !c     &    (1.0+exp(p23*(temp-p24)))

  
  !     Photo-respiration compensation point (Pa)
  !     -----------------------------------------
  mgama = p3/(p8*(p9**(p10*(temp-p11))))

  
  !     Michaelis-Menton CO2 constant (Pa)
  !     ----------------------------------     
  f2 = p12*(p13**(p10*(temp-p11)))


  !     Michaelis-Menton O2 constant (Pa)
  !     ---------------------------------
  f3 = p14*(p15**(p10*(temp-p11)))

      
  !     Saturation partial pressure of water vapour (Pa)
  call tetens(temp,es)

  !     Saturated mixing ratio (kg/kg)
  !     -----------------------------     
  rmax = 0.622*(es/(p0-es))
     
  !     Moisture deficit at leaf level (kg/kg)
  !     --------------------------------------    
  r = -0.315*rmax

  !     Internal leaf CO2 partial pressure (Pa)
  !     ---------------------------------------     
  ci = p19* (1.-(r/p20)) * (ca-mgama) + mgama

  !     Rubisco carboxilation limited photosynthesis rate (molCO2/m2/s)
  !     ---------------------------------------------------------------     
  jc = vm*((ci-mgama)/(ci+(f2*(1.+(p3/f3)))))

  !     Light limited photosynthesis rate (molCO2/m2/s)
  !     -----------------------------------------------
  if (ligth_limit) then
     aux_ipar= ipar   
  else
     aux_ipar = ipar - (ipar * 0.20)
  endif
  jl = p4*(1.0-p5)*aux_ipar*((ci-mgama)/(ci+(p6*mgama)))
  
  !     Transport limited photosynthesis rate (molCO2/m2/s)
  !     ---------------------------------------------------
  je = p7*vm
  
  !     Jp (minimum between jc and jl)
  !     ------------------------------   
  a = 0.83
  b = (-1.)*(jc+jl)
  c = jc*jl
  delta = (b**2)-4.0*a*c
  
  jp1=(-b-(sqrt(delta)))/(2.0*a)
  jp2=(-b+(sqrt(delta)))/(2.0*a)
  jp= amin1(jp1,jp2)
  
  !     Leaf level gross photosynthesis (minimum between jc, jl and je)
  !     ---------------------------------------------------------------
  a2 = 0.93
  b2 = (-1.)*(jp+je)
  c2 = jp*je
  delta2 = (b2**2)-4.0*a2*c2
  
  j1=(-b2-(sqrt(delta2)))/(2.0*a2)
  j2=(-b2+(sqrt(delta2)))/(2.0*a2)
  f1a = amin1(j1,j2)
  
  !     VPD
  !     ===
  !     buck equation...references:
  !     http://www.hygrometers.com/wp-content/uploads/CR-1A-users-manual-2009-12.pdf
  !     Hartmann 1994 - Global Physical Climatology p.351
  !     https://en.wikipedia.org/wiki/Arden_Buck_equation#CITEREFBuck1996
  
  !     ES2 = VPD-POTENTIAL - Saturation Vapor Pressure (millibar)
  call tetens(temp, es2)
  
  !     VPD-REAL = Actual vapor pressure
  vpd_real = es2 * rh       ! RESULTS are IN hPa == mbar! we want kPa (DIVIDE by 10.)
  
  !     Vapor Pressure Deficit
  vpd64 = (es2 - VPD_REAL) / 10. 
  vpd = real(vpd64, 4)
  
  !     Soil water
  !     ==========      
  wa = (w/wmax)
  
  !     Stomatal resistence
  !     ===================
  call canopy_resistence(pft , vpd64, f1a, rc)
  
  !     Water stress response modifier (dimensionless)
  !     ----------------------------------------------
  
!  c     Water stress response modifier (dimensionless)
!  c     [f5 ; Eq. 21]
!  ! vamos deixar o F5 como antigamente ate este problema ser resolvido
!  if (wa.gt.0.5) then
!     f5 = 1.0               !Not too lower in e.g. Amazonian dry season
!  else if((wa.ge.0.205).and.(wa.le.0.5)) then
!     f5 = (wa-0.205)/(0.5-0.205)
!  else if(wa.lt.0.205) then
!     f5 = wa !Below wilting point f5 accompains wa (then Sahara is well represented)
!  endif
  
  
   csru = 0.5 
   pt = csru*(cf1*1000.)*wa  !(based in Pavlick et al. 2013; *1000. converts kgC/m2 to gC/m2)
   alfm = 1.391
   gm = 3.26 * 86400.           !(*86400 transform s/mm to dia/mm)    
   if(rc .gt. 0.001) then
      gc = rc !* 1.15741e-08 ! transfor s/m  to dia/mm)  !testamos, nao muda nada! Bia vai rever
      gc = (1./gc)  ! molCO2/mm2/dia
   else
      gc =  1.0/0.001 ! BIANCA E HELENA - Mudei este esquema..   
   endif                     ! tentem entender o algoritmo
!!$c                                ! e tenham certeza que faz sentido ecologico
   d =(emax*alfm)/(1. + gm/gc) !(based in Gerten et al. 2004)
!!$!     BIanca- Eu alterei a estrutura desta equacao para evitar erros
!!$!     Isso faz a mesma coisa que o calculo que vc implementou - jp
   if(d .gt. 0.0) then
      f5_64 = pt/d
      f5_64 = exp(-1 * f5_64)
      f5_64 = 1.0 - f5_64
   else
      f5_64 = wa  ! eu acrescentei esta parte caso d seja igual a zero
  endif          ! nao sei se faz sentido!
!!$
      f5 = real(f5_64,4) !esta funcao transforma o f5 (double precision)
!     em real (32 bits) 
      
!     Photosysthesis minimum and maximum temperature
!     ----------------------------------------------
      
  if ((temp.ge.-10.0).and.(temp.le.50.0)) then
     f1 = f1a * (f5 + 0.0d0) !f5:water stress factor-- Notem que aqui a tranformacao eh de 128 pra 64 bits
  else
     f1 = 0.0               !Temperature above/below photosynthesis windown
  endif
  
  !     Leaf area index (m2/m2)
  leaf_t_months = tleaf(pft)*12. ! turnover time in months
  leaf_t_coeff = leaf_t_months/100. !1 - 100 months == ~ 1/12 to 8.3 years (TRY-kattge et al. 2011; Jedi-Pavlick 2012)
  !  c      if(leaf_t_months .gt. 0) print*, leaf_t_months, leaf_t_coeff, pft
  if (leaf_t_coeff .gt. 1.) leaf_t_coeff = 1. 
  leaf_turnover =  (365.0/12.0) * (10. **(2.0*leaf_t_coeff))
  sla = (3e-2 * (365.0/leaf_turnover)**(-0.46))     
  laia64 = ((cl1*1000.)  * sla) ! * 1000 transform kg to g - laia64 in m2 m-2
  !  c      if(laia64 .gt. 0.0) print*, laia64, 'jp'
  !     ---------------------------------
  !     Specifc Leaf Area----------------
  !      sla=((0.0300*1000.)*((365./(((tleaf(pft))/365.)/12.))**(-0.46)))    
  !      laia64 = (cl1 * 365.0 * sla)
  
!!$
!!$c      sla=(0.03*(365/(tleaf(pft)/365))**(-0.46))
!!$c      laia64 = (cl1 * 1000 * 365.0 * sla)
!!$c      if(laia64 .gt. 0.0) print*, laia64, 'bia'
  !     LAI
  !     ------
  sunlai = (1.0-(exp(-p26*laia64)))/p26
  !     --------
  shadelai = laia64 - sunlai
  
  !     Scaling-up to canopy level (dimensionless)
  !     ------------------------------------------
  f4 = (1.0-(exp(-p26*laia64)))/p26 !Sun 90 degrees in the whole canopy, to be used for respiration
  
  !     Sun/Shade approach to canopy scaling !Based in de Pury & Farquhar (1997)
  !     ------------------------------------------------------------------------
  f4sun = ((1.0-(exp(-p26*sunlai)))/p26) !sun 90 degrees
  f4shade = ((1.0-(exp(-p27*shadelai)))/p27) !sun ~20 degrees
  
  laia  = real(laia64,4) ! real((f4sun + f4shade), 4) ! pra mim faz sentido que a laia final seja
  ! a soma da lai em nivel de dossel (sun + shade) - jp
  !     Canopy gross photosynthesis (kgC/m2/yr)
  !     =======================================
  !     (0.012 converts molCO2 to kgC)
  !     (31557600 converts seconds to year [with 365.25 days])
  ph64 = 0.012*31557600.0*f1*f4sun*f4shade
  ph = real(ph64, 4) * ocp_pft       ! kg m-2 year-1
!!$  c      write(*,*) 'ph', ph
!!$  c     ============================================================
!!$  c     Autothrophic respiration
  !     ========================
  !     Maintenance respiration (kgC/m2/yr) (based in Ryan 1991)
  csa= 0.05 * ca1           !sapwood carbon content (kgC/m2). 5% of woody tissues (Pavlick, 2013)
  
  ncl = 1./29.               !(gN/gC) from lpj3 
  ncf = 1./29.               !(gN/gC)
  ncs = 1./330.              !(gN/gC)
  
  rml64 = (ncl * cl1) * 27. * exp(0.03*temp)
  
  rmf64 = (ncf * cf1) * 27. * exp(0.03*temp)
  
  rms64 = (ncs * csa) * 27. * exp(0.03*temp)
  
  rm64 = (rml64 + rmf64 + rms64)
  rm = real(rm64, 4)
  
  ! c     Growth respiration (KgC/m2/yr)(based in Ryan 1991; Sitch et al.
  ! c     2003; Levis et al. 2004)         
  
  csai =  (beta_awood * 0.05)
  
  rgl64 = 0.25 * beta_leaf * 365.
  
  rgf64 =  0.25* beta_froot * 365.
  
  rgs64 = (0.25 * csai * 365.)
  
  rg64 = (rgl64 + rgf64 + rgs64)
  rg = real(rg64,4) 
  
  if (rg.lt.0) then
     rg = 0.0
  endif
  
  !     c Autotrophic (plant) respiration -ar- (kgC/m2/yr)
  !     Respiration minimum and maximum temperature
  !     -------------------------------------------     
  if ((temp.ge.-10.0).and.(temp.le.50.0)) then
     ar64 = rm64 + rg64
     ar = real(ar64,4) * ocp_pft
     if(ar .lt. 0.00001) ar = 0.0
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
  !     No futuro proximo poderiamos usar uma npp negativa para indicar uma perda de carbono
  !     dos tecidos vegetais... e uma possivel extincao do pft/pls da celula de grid.
  return
end subroutine productivity1



!     ==================================================================
subroutine canopy_resistence(pft,vpd_in,f1_in,rc2_in)
  use global_pars
  implicit none
  integer(kind=i4),parameter :: q = npls

  !     Variables/Parameters
  !     =========
  
  !     Inputs
  !     ------      
  integer(kind=i4),intent(in) ::  pft
  real(kind=r8),intent(in) :: f1_in    !Photosynthesis (molCO2/m2/s)
  real(kind=r8),intent(in) :: vpd_in   !kPa
      
  !     Outputs
  !     -------
  real(kind=r4),intent(out) :: rc2_in               !Canopy resistence (s/m)
      
!     Internal
!     --------
  real(kind=r4),dimension(q) :: g1
  real(kind=r4) :: rcmax, rcmin
  
  real(kind=r8) :: f1b      !Photosynthesis (micromolCO2/m2/s)
  real(kind=r8) :: gs2      !Canopy conductance (m/s)
  real(kind=r8) :: gs       !Canopy conductance (molCO2/m2/s)
  real(kind=r8) :: g0       !Residual stomatance conductance
  real(kind=r8) :: D1       !kPA
  real(kind=r8) :: aa

  call pft_par(1, g1)
  
  f1b = (f1_in*10e5)        ! Helena - Mudei algumas coisas aqui
  aa = (f1b/363.)           ! Entenda o algoritmo e tenha certeza de que  
  g0 = 0.01                 ! condiz com a realidade esperada =)
  rcmax = 553.000
  rcmin = 100.000
  
  if(f1_in .le. 0.0) then 
     rc2_in = rcmax
     goto 110
  endif
  if(vpd_in .gt. 0.1) then
     goto 10
  else
     rc2_in = rcmin
     goto 110
  endif
10 continue
  if (vpd_in .gt. 0.95) then
     rc2_in = rcmax
     goto 110
  else
     D1 = sqrt(vpd_in)
     gs = g0 + 1.6 * (1. + (g1(pft)/D1)) * (aa) !Based on Medlyn et al. 2011
     if(gs .le. 0.0) then
        rc2_in = rcmax
        goto 110
     endif
  endif
  
  gs2 = (gs/41.)
  if(gs2 .le. 0.0) rc2_in = rcmax
  if(gs2 .gt. 0.0) rc2_in = real((gs2**(-1)),4)
110 continue
  return
end subroutine canopy_resistence



!     =====================================
!     Microbial (heterotrophic) respiration
!     =====================================
subroutine carbon2 (tsoil,f5c,evap,laia,cl,cs,hr)
  implicit none
  !     Variables
  !     =========
  integer,parameter :: r4 = kind(0.0)
     
  !     Inputs
  !     ------
  real(kind=r4),intent(in) :: tsoil                !Mean monthly soil temperature (oC)
  real(kind=r4),intent(in) :: f5c                  !Stress response to soil moisture (dimensionless)
  real(kind=r4),intent(in) :: evap                 !Actual evapotranspiration (mm/day)
  real(kind=r4),intent(in) :: laia

  !     Outputs 
  !     -------
  real(kind=r4),intent(out) :: cl                   !Litter carbon (kgC/m2)
  real(kind=r4),intent(out) :: cs                   !Soil carbon (kgC/m2)
  real(kind=r4),intent(out) :: hr                   !Heterotrophic (microbial) respiration (kgC/m2)

  !     Internal
  !     --------
  real(kind=r4) :: lf                   !Litterfall (kgC/m2)
  real(kind=r4) :: f6                   !Litter decayment function
  real(kind=r4) :: f7                   !Soil carbon storage function
  real(kind=r4) :: p10                  !Q10 function
  real(kind=r4) :: p11                  !Q10 function reference temperature (oC)
  real(kind=r4) :: p32                  !Q10 parameter of soil respiration sensibility to temperature
  real(kind=r4) :: p33                  !Litterfall (kg/C/yr)
  real(kind=r4) :: p34                  !Average fraction of litter carbon lost to atmosphere
  real(kind=r4) :: p35                  !Carbon soil turnover (1/20yr)
  real(kind=r4) :: p36                  !Specific rate heterotrophic respiration
  !     
  !     Initialize
  !     ----------
  !     
  lf  = 0.0
  f6  = 0.0
  f7  = 0.0
  cl  = 0.0
  cs  = 0.0
  p10 = 0.10
  p11 = 25.0
  p32 = 2.00
  p33 = 0.10
  p34 = 0.30
  p35 = 0.05
  p36 = 0.25
  
  !     Litter decayment function                                             !Controlled by annual evapotranspiration
  !     -------------------------
  f6 = 1.16*10.**(-1.4553+0.0014175*(evap*365.0))
  
  !     Soil carbon storage function                                          !Controlled by temperature
  !     ----------------------------    
  f7 = p32**(p10*(tsoil-p11))

  !     Litterfall (kgC/m2)
  !     ------------------
  lf = p33 * laia

  !     Litter carbon (kgC/m2)
  !     ----------------------  
  cl = lf/f6
  
  !     Soil carbon(kgC/m2)
  !     -------------------
  cs = ((p34*cl)/(p35*f7))*f5c
  
  !     Respiration minimum and maximum temperature
  !     -------------------------------------------
  !     
  if ((tsoil.ge.-10.0).and.(tsoil.le.50.0)) then
     hr = p36*(cl*(f6**2)+(cs*f5c*evap*(f7**2))) !Litter and Soil!respectively
  else
     hr = 0.0               !Temperature above/below respiration windown
  endif
!     
  return
end subroutine carbon2
!     ===================================================
!!$
!!$subroutine critical_value(var)
!!$  implicit none
!!$  real,intent(inout) :: var
!!$      
!!$  if(abs(var) .lt. 0.000001) var = 0.0
!!$  
!!$  return
!!$end subroutine critical_value
!!$
!!$
!!$      
!!$subroutine critical_value2(var1)
!!$  implicit none
!!$  double precision var1
!!$  
!!$  if(abs(var1) .lt. 0.00000001) var1 = 0.0
!!$  
!!$  return
!!$end subroutine critical_value2
!!$
!!$
!!$      
!!$subroutine critical_value3(var2)
!!$  implicit none
!!$  real*16 var2
!!$  
!!$  if(abs(var2) .lt. 0.000000000001) var2 = 0.0
!!$  
!!$  return
!!$end subroutine critical_value3
!     =================================================================

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
     TOTAL_BIOMASS_PFT(P) = CLEAF(P) + CFROOT(P) + (CAWOOD(P)*0.05) ! only sapwood
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
  
  RETURN
END SUBROUTINE PFT_AREA_FRAC


!     =========================================================

subroutine penman (spre,temp,ur,rn,rc2,evap)
  use global_pars
  implicit none
  !     Inputs
  !     ------
  !     
  
  real(kind=r4),intent(in) :: spre                 !Surface pressure (mb)
  real(kind=r4),intent(in) :: temp                 !Temperature (oC)
  real(kind=r4),intent(in) :: ur                   !Relative humidity (0-1,dimensionless)
  real(kind=r4),intent(in) :: rn                   !Radiation balance (W/m2)
  real(kind=r4),intent(in) :: rc2                  !Canopy resistence (s/m)
  !     
  !     
  !     Output
  !     ------
  !     
  real(kind=r4),intent(out):: evap                 !Evapotranspiration (mm/day)
  !     
  !     Parameters
  !     ----------
  real(kind=r4) :: ra, h5, t1, t2, es, es1, es2, delta_e, delta
  real(kind=r4) :: gama, gama2
  
  
  ra = 100.                  !s/m
  h5 = 0.0275               !mb-1
  !     
  !     Delta
  !     -----
  t1 = temp + 1.
  t2 = temp - 1.
  call tetens(t1,es1)       !Saturation partial pressure of water vapour at temperature T
  call tetens(t2,es2)
  delta = (es1-es2)/(t1-t2) !mbar/oC
  !     
  !     Delta_e
  !     -------
  call tetens (temp,es)
  delta_e = es*(1. - ur)    !mbar
  
  if ((delta_e.ge.(1./h5)-0.5).or.(rc2.ge.550.0)) evap = 0.
  if ((delta_e.lt.(1./h5)-0.5).or.(rc2.lt.550.0)) then
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
  
  return
end subroutine penman
!     ============================================
!     
subroutine evpot2 (spre,temp,ur,rn,evap) 
  use global_pars
  implicit none
  
  !     Inputs
  
  real(kind=r4),intent(in) :: spre                 !Surface pressure (mb)
  real(kind=r4),intent(in) :: temp                 !Temperature (oC)
  real(kind=r4),intent(in) :: ur                   !Relative humidity (0-1,dimensionless)
  real(kind=r4),intent(in) :: rn                   !Radiation balance (W/m2)
  !     Output
  !     ------
  !     
  real(kind=r4),intent(out):: evap                 !Evapotranspiration (mm/day)
  !     
  !     Parameters
  !     ----------
  real(kind=r4) :: ra, t1, t2, es, es1, es2, delta_e, delta
  real(kind=r4) :: gama, gama2, rc, rcmin

  ra      = 100.            !s/m
  rcmin   = 100.            !s/m
  !     
  !     Delta
  !     -----
  !     
  t1 = temp + 1.
  t2 = temp - 1.
  call tetens(t1,es1)
  call tetens(t2,es2)
  delta = (es1-es2)/(t1-t2) !mb/oC
  !     
  !     Delta_e
  !     -------
  !     
  call tetens (temp,es)
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
  evap = (delta*rn + (1.20*1004./ra)*delta_e)/(delta+gama2) !W/m2
  evap = evap*(86400./2.45e6) !mm/day
  evap = amax1(evap,0.)     !Eliminates condensation
  !     
  return
end subroutine evpot2
!     
!     =================================================================

!     
subroutine runoff (wa,roff)
  use global_pars
  implicit none
  
  real(kind=r8), intent( in) :: wa
  real(kind=r4), intent(out) :: roff
  real(kind=rbig) :: roff64
  
  roff64 = 38.*(wa**11.)
!  c      roff64 = 11.5*(wa**6.6) * 1000. !From NCEP-NCAR Reanalysis data
   roff = real(roff64, 4)
   return
 end subroutine runoff
!     
!     =================================================================    
 subroutine tetens (t,es)  !Saturation Vapor Pressure (hPa)!
   use global_pars
   implicit none
   
   real(KIND=R4),intent( in) :: t
   real(KIND=R4),intent(out) :: es
   
   if (t.ge.0.) then
      es = 6.1121 * exp((18.678-(t/234.5))*(t/(257.14+t))) ! mbar == hPa
   else
      es = 6.1115 * exp((23.036-(t/333.7))*(t/(279.82+t))) ! mbar == hPa
   endif
   return
 end subroutine tetens
 


 
 !=====================================================================
 !c     subroutine allocation calculates the daily carbon content of each
 !c     compartment
 !c     
 !c     code written by Bianca Rius & David Lapola (27.Ago.2015)
 !c     
 !c=====================================================================
      
 subroutine allocation (pft, npp ,scl1,sca1,scf1,scl2,sca2,scf2)
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

   real(kind=rbig) :: scf2_128, sca2_128, scl2_128
   
   real(kind=r4), dimension(npfts) :: aleaf             !npp percentage allocated compartment
   real(kind=r4), dimension(npfts) :: aawood
   real(kind=r4), dimension(npfts) :: afroot
   real(kind=r4), dimension(npfts) :: tleaf             !turnover time (yr)
   real(kind=r4), dimension(npfts) :: tawood
   real(kind=r4), dimension(npfts) :: tfroot            
   
   call pft_par(3, aleaf)
   call pft_par(4, aawood)
   call pft_par(5, afroot)
   call pft_par(6, tleaf)
   call pft_par(7, tawood)
   call pft_par(8, tfroot)
    
   !     Carbon content of each compartment(KgC/m2)
   !c     
   !c     
   !c     initialization
   if((scl1 .lt. 1e-12) .or. (scf1 .lt. 1e-12)) then
      IF(NPP .lt. 1e-12) THEN
         scl2 = 0.0
         scf2 = 0.0
         sca2 = 0.0 
         goto 10
      ENDIF
   endif
   npp_aux = npp/365.0       !transform (KgC/m2/yr) in (KgC/m2/day)
   scl2_128 = scl1 + (aleaf(pft) * npp_aux) -(scl1 /(tleaf(pft)*365.0))
   scf2_128 = scf1 +(afroot(pft) * npp_aux)-(scf1 /(tfroot(pft)*365.0))
   if(aawood(pft) .gt. 0.0) then
      sca2_128 = sca1 +(aawood(pft)*npp_aux)-(sca1/(tawood(pft)*365.0))
      sca2 = real(sca2_128,4)
   else
      sca2 = 0.0
   endif
   
   scf2 = real(scf2_128,4)
   scl2 = real(scl2_128,4)
   
   
   if(scl2 .lt. 0.0) scl2 = 0.0
   if(scf2 .lt. 0.0) scf2 = 0.0
   if(sca2 .lt. 0.0) sca2 = 0.0
   
10 continue
   return
 end subroutine allocation
 

 
 !     ==============================================================     
 subroutine pft_par(par, dt) !!!!!!!!!  mudamos os valores de dt
   implicit none

   integer,parameter :: i4 = kind(0)
   integer,parameter :: r4 = kind(0.0)
   
   integer(kind=i4),parameter :: vars = 12
   
   !     input
   integer(kind=i4),intent(in) :: par            ! parameter number 
   real(kind=r4),dimension(vars),intent(out) :: dt
   
   real(kind=r4),dimension(vars) :: dt1,dt2,dt3,dt4,dt5,dt6,dt7,dt8
   
   
   !     dt1 = g1
   !     dt2 = p21 
   !     dt3 = aleaf
   !     dt4 = aawood
   !     dt5 = afroot
   !     dt6 = tleaf
   !     dt7 = tawood
   !     dt8 = tfroot
     
   !     PFTS (Leaf Gas Exchange Database
   
   !     1 = 
   !     2 = Tropical Deciduous Tree        
   !     3 = Tropical Deciduous Savanna 
   !     4 = Tropical Evergreen Savanna 
   !     5 = Tropical Deciduous Grass
   !     6 = Temperate Evergreen Tree 
   !     7 = Temperate Deciduous Tree
   !     8 = Temperate Deciduous Grass
   !     9 = Temperate Deciduous Shrub 
   !    10 = Temperate Evergreen Shrub 
   !    11 = Boreal Evergreen Tree  
   !    12 = Boreal Deciduous Tree
   
   !PFT       1       2       3       4       5       6       7      8       9       10      11      12      
   data dt1/3.77,   4.15,   2.98,   7.18,   4.5,    3.37,   4.64,   4.4,    4.6,    3.92,   1.5,    2.72/   !g1
   data dt2/5.9E-5, 3.4E-5, 3.2E-5, 6.8E-5, 3.1E-5, 5.1E-5, 3.3E-5, 3.1E-5, 3.1E-5, 4.4E-5, 4.2E-5, 4.0E-5/ !p21
   data dt3/0.30,   0.35,   0.35,   0.30,   0.45,   0.30,   0.35,   0.45,   0.35,   0.40,   0.30,   0.35/   !aleaf
   data dt4/0.35,   0.35,   0.20,   0.25,   0.0,    0.35,   0.40,   0.00,   0.20,   0.20,   0.40,   0.35/   !aawood
   data dt5/0.35,   0.30,   0.45,   0.45,   0.55,   0.35,   0.25,   0.55,   0.45,   0.40,   0.30,   0.30/   !afroot
   data dt6/2.0,    1.0,    1.0,    2.0,    2.0,    3.0,    1.0,    2.0,    1.0,    3.0,    3.0,    1.0/    !tleaf
   data dt7/30.0,   30.0,   20.0,   20.0,   0.0,    50.0,   50.0,   0.0,    35.0,   35.0,   50.0,   50.0/   !tawood
   data dt8/3.0,    3.0,    3.5,    3.0,    3.0,    3.5,    3.5,    3.5,    3.5,    3.5,    3.5,    3.5/    !tfroot

   if(par .eq. 1 ) then      ! g1
      dt(:) = dt1(:)
   else if(par .eq. 2) then  ! p21
      dt(:) = dt2(:)
   else if(par .eq. 3) then  ! aleaf
      dt(:) = dt3(:)
   else if(par .eq. 4) then  ! awood
      dt(:) = dt4(:)
   else if(par .eq. 5) then  ! afroot
      dt(:) = dt5(:)
   else if(par .eq. 6) then  ! tleaf
      dt(:) = dt6(:)
   else if(par .eq. 7) then  ! tawood
      dt(:) = dt7(:)
   else if(par .eq. 8) then  ! tfroot
      dt(:) = dt8(:)
   else
      print*, "your search failed"
   endif
   
   return
 end subroutine pft_par
 
 !c     ==================================================
 subroutine spinup(nppot,cleafini,cfrootini,cawoodini)
   !c     &     cbwoodini,cstoini,cotherini,crepini) 
   IMPLICIT NONE
   integer,parameter :: i4 = kind(0)
   integer,parameter :: r4 = kind(0.0)
   integer,parameter :: r8 = kind(0.0D0)
   
   integer(kind=i4),parameter :: npfts = 12
   integer(kind=i4),parameter :: nt=30000
   
   !   c     inputs
   integer(kind=i4) :: i6, kk, k
   
   real(kind=r4),intent(in) :: nppot
   real(kind=r4) :: sensitivity
   
   !c     outputs
   real(kind=r4),intent(out) :: cleafini(npfts)
   real(kind=r4),intent(out) :: cawoodini(npfts)
   real(kind=r4),intent(out) :: cfrootini(npfts)
   real(kind=r8) :: cleafi_aux(nt)
   real(kind=r8) :: cfrooti_aux(nt)
   real(kind=r8) :: cawoodi_aux(nt)
   
   
   real(kind=r4) :: aleaf(npfts)             !npp percentage alocated to leaf compartment
   real(kind=r4) :: aawood (npfts)           !npp percentage alocated to aboveground woody biomass compartment
   real(kind=r4) :: afroot(npfts)            !npp percentage alocated to fine roots compartmentc 
   real(kind=r4) :: tleaf(npfts)             !turnover time of the leaf compartment (yr)
   real(kind=r4) :: tawood (npfts)           !turnover time of the aboveground woody biomass compartment (yr)
   real(kind=r4) :: tfroot(npfts)            !turnover time of the fine roots compartment
   
   
   call pft_par(3, aleaf)
   call pft_par(4, aawood)
   call pft_par(5, afroot)
   call pft_par(6, tleaf)
   call pft_par(7, tawood)
   call pft_par(6, tfroot)
   
   
   sensitivity = 1.001
   if(nppot .lt. 0.0) goto 200
   do i6=1,npfts
      do k=1,nt
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
   return
 end subroutine spinup
 !     ================================
 
