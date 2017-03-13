module water_balance
  implicit none
  private
  
  public :: wbm
  
  
contains
  subroutine wbm (prec,temp,p0,par,rhs,cleaf_ini,cawood_ini&
       &,cfroot_ini,emaxm, tsoil, photo_pft,aresp_pft,npp_pft,lai_pft&
       &,clit_pft,csoil_pft, hresp_pft,rcm_pft,runom_pft,evapm_pft&
       &,wsoil_pft,rm_pft,rg_pft,cleaf_pft,cawood_pft,cfroot_pft,grid_area)
    
    use global_pars
    use budget
    use photo, only:soil_temp_sub
    
    implicit none
    
    integer(kind=i4),parameter :: q = npls ! plss futuramente serao
    integer(kind=i4),parameter :: nt = ntimes ! (meses) ! modelo estac. --
  
    
    !     --------------------------I N P U T S----------------------------
    real(kind=r4),intent(in) :: p0(nt)         !Atmospheric pressure (mb)
    real(kind=r4),intent(in) :: prec(nt)       !Precipitation (mm/month)
    real(kind=r4),intent(in) :: temp(nt)       !Temperature (oC)
    real(kind=r4),intent(in) :: par(nt)        !IPAR (Ein/m2/s)
    real(kind=r4),intent(in) :: rhs(nt)        !Relative humidity
    
    real(kind=r4),intent(in) ::  cleaf_ini(q)  ! Initial carbon content in leaves (kg m-2)
    real(kind=r4),intent(in) :: cawood_ini(q)  ! Initial carbon content in aboveground wood (kg m-2)
    real(kind=r4),intent(in) :: cfroot_ini(q)  ! Initial carbon content in fineroots (kg m-2)
    !     -----------------------------E N D-------------------------------
    
    !     -------------------------O U T P U T S---------------------------
    
    real(kind=r4),intent(out) :: tsoil(nt)       !soil temperature
    
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
    call soil_temp_sub(temp,tsoil)
    !     ============
    !     Water budget
    !     ============
    !     
    !     For all grid points
    !     -------------------     
    cleaf1_pft =  cleaf_ini
    cawood1_pft = cawood_ini
    cfroot1_pft = cfroot_ini
    
    leaf0 = cleaf_ini
    froot0 = cfroot_ini
    awood0 = cawood_ini
    
    cleaf_pft  = 0.0 ! leaf biomass (KgC/m2)
    cawood_pft = 0.0 ! aboveground biomass (KgC/m2)
    cfroot_pft = 0.0 ! fine root biomass (KgC/m2)
    grid_area  = 0.0 ! gridcell area fraction of pfts(%)
    
    wg0 = -1.0 !Soil water content in preceeding year integration
    emaxm =  0.0 !Maximum evapotranspiration
    photo_pft = 0.0 !Monthly photosynthesis (kgC/m2)
    aresp_pft = 0.0 !Monthly autotrophic respiration (kgC/m2)
    npp_pft   = 0.0 !Monthly net primary productivity (average between PFTs) (kgC/m2)
    lai_pft   = 0.0 !Monthly leaf area index
    clit_pft  = 0.0 !Monthly litter carbon
    csoil_pft = 0.0 !Monthly soil carbon
    hresp_pft = 0.0 !Monthly heterotrophic respiration (kgC/m2)
    rcm_pft   = 0.0
    gsoil     = 0.0 !Soil ice
    ssoil     = 0.0 !Soil snow
    runom_pft = 0.0 !Runoff
    evapm_pft = 0.0 !Actual evapotranspiration        
    wsoil_pft = 0.0 !Soil moisture (mm)
    rm_pft    = 0.0
    rg_pft    = 0.0
    
    wini  = 0.01  !Soil moisture_initial condition (mm)
    gini  = 0.0   !Soil ice_initial condition (mm)
    sini  = 0.0   !Overland snow_initial condition (mm)
    
    !     =================
    !     START INTEGRATION
    !     =================
    
    n = 0
10  continue
    n = n + 1  
    k = mod(n,12)
    if (k.eq.0) k = 12
    mes = k
    spre = p0(k) * 0.01 ! transforamando de Pascal pra mbar (hPa)
    td = tsoil(k)
    ta = temp(k)
    pr = prec(k)
    ipar = par(k)/2.18e5
    ru = rhs(k) / 100.
    ae = 2.895*ta+52.326 !Available energy (W/m2) - From NCEP-NCAR reanalysis data
    !     
    !     Monthly water budget
    !     ====================
    epmes = 0.0  
    wfim = 0.0
    gfim = 0.0
    sfim = 0.0
    smes = 0.0
    rmes = 0.0
    emes = 0.0
    phmes = 0.0
    armes = 0.0
    nppmes = 0.0
    laimes = 0.0
    clmes = 0.0
    csmes = 0.0
    hrmes = 0.0
    rcmes = 0.0
    rmmes = 0.0
    rgmes = 0.0
    cleafmes = 0.0
    cawoodmes = 0.0
    cfrootmes = 0.0
    gridocpmes = 0.0
    
    
    call budget (mes,wini,gini,sini,td,ta,pr,spre,ipar,ru&
         &,cleaf1_pft,cawood1_pft,cfroot1_pft ,wfim,gfim, sfim,smes&
         &,rmes,emes,epmes,phmes,armes,nppmes,laimes,clmes,csmes,hrmes&
         &,rcmes,rmmes,rgmes,cleafmes,cawoodmes,cfrootmes, gridocpmes)
    
    emaxm(k) = epmes
    do p=1,q
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
    enddo
    
    wini = wfim
    gini = gfim
    sini = sfim
    
    cleaf1_pft  = cleafmes 
    cawood1_pft = cawoodmes
    cfroot1_pft = cfrootmes
    
    if(k .eq. 12) then
       cleaf_pft  = cleafmes
       cawood_pft = cawoodmes
       cfroot_pft = cfrootmes
       grid_area  = gridocpmes
    endif
    
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
  ! continue 100 goes here for pft check
  
  
  !     PFTs equilibrium check
  !     ==================================================
  !     tentativa 1 - usando variacao no pool de C vegetal
  !     --------------------------------------------------
  if (k.eq.12) then
     nerro = 0
     biomass = 0.0
     biomass0 = 0.0
     check = .false.
     sensi = 0.5! (kg/m2/y) if biomas change .le. sensi: equilibrium
     ! brienen et al. 2015 mean biomass change in Amazon forest - 1995 - 1.0 Mg/ha/y
     ! here I'm using a larger value (5 Mg/ha/y) 
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
  ! continue 100 goes here to exclude pft_check
100 continue  
  
end subroutine wbm

end module water_balance
