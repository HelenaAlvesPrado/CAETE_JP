c234567
      subroutine wbm (prec,temp,lsmk,p0,ca,par,
     &     cleafini,cfrootini,cawoodini,
     &     cbwoodini,cstoini,cotherini,crepini, ! inputs
     &     emaxm, tsoil, photo_pft,aresp_pft,npp_pft,lai_pft, ! out
     &     clit_pft,csoil_pft, hresp_pft,rcm_pft,runom_pft,
     &     evapm_pft,wsoil_pft,
     &     monrml,monrmf,monrms,monrm,
     &     monrgl,monrgf,monrgs,monrg,
     &     cleaf_pft,cfroot_pft,cawood_pft,cbwood_pft, !outputs
     &     csto_pft, crep_pft, cother_pft)
      
c     
c=======================================================================
c     
c     Water balance model (WBM). From monthly climatologies of
c     precipitation and surface temperature, the WBM calculates the
c     environmental variables.
c     
c     05Jul2005, MDO: ae, rh & runoff are changed.
c     11Jul2005, MDO: wsoil2 is written (for testing purpose).
c     31Ago2006, DML: carbon cycle is included
c     
c=======================================================================
c     
c     i/o variables
      integer, parameter :: npft=3
      integer, parameter ::nx=192,ny=96
      
      real lsmk(nx,ny)
      real   p0(nx,ny)
      real prec(nx,ny,12)
      real temp(nx,ny,12)
      real  par(nx,ny,12)
      
      real  cleafini(nx,ny,npft)
      real cfrootini(nx,ny,npft)
      real cawoodini(nx,ny,npft)
      real cbwoodini(nx,ny,npft)
      real   cstoini(nx,ny,npft)
      real cotherini(nx,ny,npft)
      real   crepini(nx,ny,npft)
      
C     OUTPUTS PARA WBM
      real emaxm(nx,ny,12)
      real tsoil(nx,ny,12)
      
c     agora as variaveis para pfts
      real photo_pft(nx,ny,12,npft) !Monthly photosynthesis   (kgC/m2)
      real aresp_pft(nx,ny,12,npft) !Monthly autotrophic res  (kgC/m2)
      real   npp_pft(nx,ny,12,npft) !Monthly net primary produ (kgC/m2)
      
      real   lai_pft(nx,ny,12,npft) !Monthly leaf area index
      real  clit_pft(nx,ny,12,npft) !Monthly litter carbon
      real csoil_pft(nx,ny,12,npft) !Monthly soil carbon
      real hresp_pft(nx,ny,12,npft) !Monthly het resp          (kgC/m2)
      real   rcm_pft(nx,ny,12,npft) !stomatal (or canopy?) resistence
      
c     VARIAVEIS HIDROLOGICAS IMPORTANTES   
      real runom_pft(nx,ny,12,npft) !Runoff
      real evapm_pft(nx,ny,12,npft) !Actual evapotranspiration        
      real wsoil_pft(nx,ny,12,npft) !Soil moisture (mm)
      
      
c     carbon in plant tissues
      real  cleaf_pft(nx,ny,12,npft)
      real cawood_pft(nx,ny,12,npft)
      real cbwood_pft(nx,ny,12,npft)
      real cfroot_pft(nx,ny,12,npft)
      real   csto_pft(nx,ny,12,npft)
      real   crep_pft(nx,ny,12,npft)
      real cother_pft(nx,ny,12,npft)
c     maintenance and growth respiration
      
      real monrml(nx,ny,12,npft)
      real monrmf(nx,ny,12,npft)
      real monrms(nx,ny,12,npft)
      real  monrm(nx,ny,12,npft)
      real monrgl(nx,ny,12,npft)
      real monrgf(nx,ny,12,npft)
      real monrgs(nx,ny,12,npft)
      real  monrg(nx,ny,12,npft)

c     --------------------------------E N D-----------------------------


!     VARIAVEIS INTERNAS
 
 
      integer i, j, k, p, kk, mes, nerro
      real, parameter :: no_data = -9999.0
      
      real gsoil(nx,ny,12,npft)    !Soil ice 
      real ssoil(nx,ny,12,npft)    !Soil snow
      real snowm(nx,ny,12,npft)    !Snowmelt
      real   wg0(nx,ny,12,npft)    !Moisture of the previous year
      
      

C     soil temperature parameters and vars
      real, parameter :: H = 1.0 !Soil layer(m) 
      real, parameter :: diffu = 4.e-7*(30.*86400.0) !Soil thermal diffusivity (m2/month)
      real, parameter :: tau = (H**2)/(2.0*diffu) !E-folding time (months)
C      real, parameter :: auxs = -100.0 !Auxiliar for calculation of Snpp
      real t0,t1

      
      real wsaux1
      real ae,dwww,wmax
      real sini,gfim,sfim,gini,wfim,wini
      real pr,ice,spre,ta,td,ipar
      real rmes,phmes,smes,rcmes,hrmes,nppmes,laimes
      real armes,clmes,csmes,emes,epmes

c     Soil temperature
c     ---------------
c     for all grid points
      do i=1,nx
         do j=1,ny
c     
c     initialize soil temperature
            do k=1,12
               tsoil(i,j,k) = no_data
            enddo
c     
c     only for land grid points
            if (int(lsmk(i,j)).ne.0) then
               t0 = 0.          !initialization
               do n=1,1200      !100 yr (1200 months) run to attain equilibrium
                  k = mod(n,12)
                  if (k.eq.0) k = 12
                  t1 = t0*exp(-1.0/tau) + (1.0 - exp(-1.0/tau))*temp(i,j
     $                 ,k)
                  tsoil(i,j,k) = (t0 + t1)/2.0
                  t0 = t1
               enddo
            endif
c     
         enddo
      enddo
      
      
c     Water budget
c     ------------
c     
c     for all grid points
      do i=1,nx
         do j=1,ny
c     
c     write to track program execution
            if ((mod(j,ny).eq.0).and.(mod(i,10).eq.0))
     &           write(*,*) 'water balance:',i

            do k=1,12
               emaxm(i,j,k) = no_data !average maximum evapotranspiration (mm/day)
               do p=1,npft
                  photo_pft(i,j,k,p) = no_data
                  aresp_pft(i,j,k,p) = no_data
                  npp_pft(i,j,k,p) = no_data
                  lai_pft(i,j,k,p) = no_data
                  wsoil_pft(i,j,k,p) = no_data !soil moisture (mm)
                  runom_pft(i,j,k,p) = no_data !average runoff (mm/day)
                  evapm_pft(i,j,k,p) = no_data !average actual evapotranspiration (mm/day)
                  
                  rcm_pft(i,j,k,p) = no_data !average canopy resistance (s/m)
                  
                  gsoil(i,j,k,p) = no_data !soil ice (mm)
                  ssoil(i,j,k,p) = no_data !soil snow (mm)
                  snowm(i,j,k,p) = no_data !average snowmelt (mm/day)
                  wg0(i,j,k,p) = no_data !soil moisture of the previous year (mm)
                  
                  monrml(i,j,k,p) = no_data !mean monthly leaf maintenance respiration (KgC/m2/day)
                  monrmf(i,j,k,p) = no_data !mean monthly fine roots maintenance respiration (KgC/m2/day)
                  monrms(i,j,k,p) = no_data !mean monthly sapwood maintenance respiration (KgC/m2/day)
                  monrm(i,j,k,p) = no_data !mean monthly total maintenance respiration (KgC/m2/day)
                  monrgl(i,j,k,p) = no_data !mean monthly leaf growth respiration (KgC/m2/day)
                  monrgf(i,j,k,p) = no_data !mean monthly fine roots growth respiration (KgC/m2/day)
                  monrgs(i,j,k,p) = no_data !mean monthly sapwood growth respiration (KgC/m2/day)
                  monrg(i,j,k,p) = no_data !mean monthly total growth respiration (KgC/m2/day)
                  
                  
                  cleaf_pft(i,j,k,p)= no_data !mean monthly leaf biomass for all the PFTs (KgC/m2)
                  cawood_pft(i,j,k,p)= no_data !mean monthly aboveground biomass 
                  cfroot_pft(i,j,k,p)= no_data !mean monthly leaf biomass
                  cbwood_pft(i,j,k,p)= no_data !mean monthly aboveground biomass 
                  csto_pft(i,j,k,p)= no_data !mean monthly leaf biomass 
                  crep_pft(i,j,k,p)= no_data !mean monthly aboveground biomass 
                  cother_pft(i,j,k,p)= no_data !mean monthly aboveground biomass 
                  
               enddo
            enddo
            
c     
c     only for land grid points
            if (int(lsmk(i,j)).ne.0) then
               
               spre = p0(i,j)   !surface pressure (mb)
               
               do p=1,npft
!     initialize variables
                  wini  = 0.01  !soil moisture initial condition (mm)
                  gini  = 0.0   !soil ice initial condition (mm)
                  sini  = 0.0   !overland snow initial condition (mm)
                  
                  
                  do k=1,12   
                     wg0(i,j,k,p) = -1.0
                  enddo
                  
!     initialize here variables containing initial values
!     of npp alloc. for veg_pools 
                  
                  cleaf_ini  = cleafini(i,j,p)
                  cfroot_ini = cfrootini(i,j,p)
                  cawood_ini = cawoodini(i,j,p)
                  cbwood_ini = cbwoodini(i,j,p)
                  csto_ini   = cstoini(i,j,p)
                  cother_ini = cotherini(i,j,p)
                  crep_ini   = crepini(i,j,p)
                  
                  
                  
                  
c     start integration
                  n = 0
 10               continue
                  n = n + 1
                  
c     
c     pre-processing
                  k = mod(n,12)
                  if (k.eq.0) k = 12
                  mes = k
                  td = tsoil(i,j,k)
                  ta = temp(i,j,k)
                  pr = prec(i,j,k)
                  ipar = par(i,j,k)
                  
!     ae = 2.26457*ta + 67.5876 !available energy (W/m2) [Eq. 8]
                  ae = 2.895*ta + 52.326 !from NCEP-NCAR Reanalysis data
                  
c     monthly water budget
                  call budget (p,mes,wini,gini,sini,td,ta,pr,spre,ae,ca,
     &                 ipar,cleaf_ini,cawood_ini,cfroot_ini
     $                 ,cbwood_ini,csto_ini,cother_ini,crep_ini,
     $                 wfim,gfim,sfim,smes,rmes,emes,epmes,phmes,armes
     $                 ,nppmes,laimes,rmlmes,rmfmes,rmsmes,rmmes,rglmes
     $                 ,rgfmes,rgsmes,rgmes,clmes,csmes,hrmes,rcmes
     $                 ,cleafavg,cawoodavg,cfrootavg,cbwoodavg,cstoavg
     $                 ,crepavg,cotheravg)  
                  
                  
c     update variables
                  
                  
                  
                  if(pft .eq. 1) emaxm(i,j,k) = epmes
                  
                  gsoil(i,j,k,p) = gfim
                  ssoil(i,j,k,p)= sfim
                  snowm(i,j,k,p) = smes
                  
                  runom_pft(i,j,k,p) = rmes
                  evapm_pft(i,j,k,p) = emes
                  wsoil_pft(i,j,k,p) = wfim
                  
                  npp_pft(i,j,k,p) = nppmes
                  photo_pft(i,j,k,p) = phmes
                  aresp_pft(i,j,k,p) = armes
                  rcm_pft(i,j,k,p) = rcmes
                  lai_pft(i,j,k,p) = laimes
                  
                  clit_pft(i,j,k,p) = clmes
                  csoil_pft(i,j,k,p) = csmes
                  hresp_pft(i,j,k,p) = hrmes
                  
                  monrml(i,j,k,p) = rmlmes
                  monrmf(i,j,k,p) = rmfmes
                  monrms(i,j,k,p) = rmsmes
                  monrm(i,j,k,p) = rmmes
                  monrgl(i,j,k,p) = rglmes
                  monrgf(i,j,k,p) = rgfmes
                  monrgs(i,j,k,p) = rgsmes
                  monrg(i,j,k,p) = rgmes
                  
                  cleaf_pft(i,j,k,p) = cleafavg
                  cawood_pft(i,j,k,p) = cawoodavg
                  cfroot_pft(i,j,k,p) = cfrootavg
                  cbwood_pft(i,j,k,p) = cbwoodavg
                  csto_pft(i,j,k,p) = cstoavg
                  crep_pft(i,j,k,p) = crepavg
                  cother_pft(i,j,k,p) = cotheravg
                  
                  wini = wfim
                  gini = gfim
                  sini = sfim 
c     
c     check if equilibrium is attained (k=12)
                  if (k.eq.12) then
                     wmax = 500.
                     nerro = 0
                     do kk=1,12
                        dwww = (wsoil_pft(i,j,kk,p)+gsoil(i,j,kk,p)
     $                       -wg0(i,j,kk,p))/wmax
                        if (abs(dwww).gt.0.001) nerro = nerro + 1
                     enddo
                     if (nerro.ne.0) then
                        do kk=1,12
                           wg0(i,j,kk,p) = wsoil_pft(i,j,kk,p) + gsoil(i
     $                          ,j,kk,p)
                        enddo
                     else
                        goto 100
                     endif
                  endif
c     
                  goto 10
 100              continue
                  
               enddo            ! pfts loop
            endif               ! lsmk branch
         enddo                  ! ny loop
      enddo                     ! nx loop
      END SUBROUTINE WBM
      
c=======================================================================
      
      
c$$$c     monthly water budget
c$$$                  call budget (p,mes,wini,gini,sini,td,ta,pr,spre,ae,ca,
c$$$     &                 ipar,cleaf_ini,cawood_ini,cfroo_ini
c$$$     $                 ,cbwood_ini,csto_ini,cother_ini,crep_ini,
c$$$     $                 wfim,gfim,sfim,smes,rmes,emes,epmes,phmes,armes
c$$$     $                 ,nppmes,laimes,rmlmes,rmfmes,rmsmes,rmmes,rglmes
c$$$     $                 ,rgfmes,rgsmes,rgmes,clmes,csmes,hrmes,rcmes
c$$$     $                 ,cleafavg,cawoodavg,cfrootavg,cbwoodavg,cstoavg
c$$$     $                 ,crepavg,cotheravg)  
                  
!     os valores iniciais para alloc e C em cada veg_pool precisam ser
!     inputs
      
      subroutine budget(pft,month,w1,g1,s1,ts,temp,prec,p0,ae,ca
     $     ,ipar,cleaf_ini,cawood_ini,cfroot_ini,cbwood_ini,csto_ini
     $     ,cother_ini,crep_ini,w2,g2,s2,smavg,ruavg,evavg,epavg,phavg
     $     ,aravg,nppavg,laiavg,rmlavg,rmfavg,rmsavg,rmavg,rglavg,rgfavg
     $     ,rgsavg,rgavg,clavg,csavg,hravg,rcavg,cleafavg,cawoodavg
     $     ,cfrootavg,cbwoodavg,cstoavg,crepavg,cotheravg)
c=======================================================================
c     
c     Surface water (soil moisture, snow and ice) budget for a single
c     month.
c     
c     I/O variables
c     -------------
c     input  month : actual month (1-12)
c     w1    : initial (previous month last day) soil moisture storage
c     (mm)
c     g1    : initial soil ice storage (mm)
c     s1    : initial overland snow storage (mm)
c     tsoil : soil temperature (oC)
c     temp  : surface air temperature (oC)
c     prec  : precipitation (mm/day)
c     p0    : surface pressure (mb)
c     ae    : available energy (W/m2)
c     cleaf1: initial (previous month last day) carbon content on leaf
c     compartment (kgC/m2)
c     output w2    : final (last day) soil moisture storage (mm)
c     g2    : final soil ice storage (mm)
c     s2    : final overland snow storage (mm)
c     smavg : snowmelt monthly average (mm/day)
c     ruavg : runoff monthly average (mm/day)
c     evavg : actual evapotranspiration monthly average (mm/day)
c     epavg : maximum evapotranspiration monthly average (mm/day)
c     cleaf2: month final carbon content on leaf compartment (kgC/m2)
c     cleafavg: average leaf carbon content (kgC/m2)
c     
c=======================================================================
c
c       helena
c            subroutine budget (pft,month,w1,g1,s1,ts,temp,prec,p0,ae,
c     &     ca,ipar,w2,g2,s2,smavg,ruavg,evavg,epavg,phavg,
c     &     aravg,nppavg,laiavg,clavg,csavg,hravg,rcavg)
c-------------
c
c     i/o variablesa

      integer month, pft
      real w1,g1,s1,ts,temp,prec,p0,ae,ca,ipar,
     &     w2,g2,s2,smavg,ruavg,evavg,epavg,rcavg,
     &     phavg,aravg,nppavg,laiavg,
     &     clavg,csavg,hravg     
c     internal variables
      real rh,wmax,tsnow,tice
      real psnow,prain
      real w,g,s
      real rimelt,smelt,roff,evap,emax
      integer ndmonth(12)       !number of days for each month
      data ndmonth /31,28,31,30,31,30,31,31,30,31,30,31/
c     carbon cycle

      real ph,ar,nppa,nppb,laia,cl,cs,hr,
     &     rm,rml,rmf,rms,rg,rgl,rgf,rgs,
     &     rmlavg,rmfavg,rmsavg,rmavg,
     &     rglavg,rgfavg,rgsavg,rgavg,
     &     sla

c     carbon allocation
      real cleaf_ini,cawood_ini,cfroot_ini,cbwood_ini,csto_ini,
     &     crep_ini, cl1,cl2,ca1,ca2,cf1,cf2,cb1,cb2,cs1,cs2,
     &     cr1,cr2,co1,co2, cl2_pft,ca2_pft,cf2_pft,cb2_pft,cs2_pft,
     &     cr2_pft,co2_pft
  
c     
c     parameters
c     rh    = 0.6   !relative humidity (adimensional)
      rh    = 0.685             !from NCEP-NCAR Reanalysis data
      wmax  = 500.0             !soil moisture availability (mm)
      tsnow = -1.0              !temperature threshold for snowfall (oC)
      tice  = -2.5              !temperature threshold for soil freezing (oC)
c     
c     precipitation [Eq. 3]
      psnow = 0.0
      prain = 0.0
      
      if (temp.lt.tsnow) then
         psnow = prec/real(ndmonth(month)) !snowfall (mm/day)
      else
         prain = prec/real(ndmonth(month)) !rainfall (mm/day)
      endif
c     
c     initialization
      w = w1                    !w = daily soil moisture storage (mm)
      g = g1                    !g = daily soil ice storage (mm)
      s = s1                    !s = daily overland snow storage (mm)
      
      smavg = 0.
      ruavg = 0.
      evavg = 0.
      epavg = 0.

      rcavg = 0.
      phavg = 0.
      aravg = 0.
      nppavg = 0.
      clavg = 0.
      csavg = 0.
      hravg = 0.
      laiavg = 0.

      rmlavg = 0.
      rmfavg = 0.
      rmsavg = 0.
      rmavg = 0.
      rglavg = 0.
      rgfavg = 0.
      rgsavg = 0.
      rgavg = 0.
      
      cleafavg = 0.
      cawoodavg = 0.
      cfrootavg = 0.
      cbwoodavg = 0.
      cstoavg = 0.
      crepavg = 0.
      cotheravg = 0.

      
c     numerical integration
      do i=1,ndmonth(month)

         if ((i .eq. 1) .and. (month.eq.1)) then  ! se for mes 1, pegue os valores vindos de spinup (inputs de wbm)
            cl1 = cleaf_pft  ! variavel que recebe o tamanho do C pool (Kg m-2)
            ca1 = cawood_pft    ! membros do lado direito sao inputs para a budget
            cf1 = cfroot_pft    ! dentro da budget apenas valores escalares! se n complica d+
            cb1 = cbwood_pft
            cs1 = csto_pft  
            cr1 = crep_pft
            co1 = cother_pft
            
            beta_leaf=0. 
            beta_awood=0.
            beta_bwood=0.
            beta_froot=0.
            
            
         else                   ! PEGUE OS VALORES DO DIA ANTERIOR
            cl1 = cl2_pft       !transforma o valor do dia anterior no valor atual
            ca1 = ca2_pft
            cf1 = cf2_pft
            cb1 = cb2_pft
            cs1 = cs2_pft       !transforma o valor do dia anterior no valor atual
            cr1 = cr2_pft
            co1 = co2_pft
            
            beta_leaf = alfa_leaf
            beta_awood = alfa_awood
            beta_froot = alfa_froot
            beta_bwood = alfa_bwood
            
         endif
         
c     carbon cycle (photosynthesis, plant respiration and NPP)
         call carbon1 (pft, temp,p0,w,wmax,ca,ipar,ts,emax !output !input !input
     $        ,rc2,cl1,ca1,cf1,cb1,beta_leaf, beta_awood, beta_bwood
     $        ,beta_froot,ph,ar,nppa,laia,f5,rm,rml,rmf,rms
     $        ,rg,rgl,rgf,rgs)        
         
         
         
!     carbon allocation (carbon content on each compartment)     
         call allocation (pft, nppa,cl1,ca1,cf1,cb1,cs1,cr1,co1 !output !input
     $        ,cl2,ca2,cf2,cb2,cs2,cr2,co2)  
         
         alfa_leaf  = cl2 - cl1 !it calculates the difference btween the biomass of the last and actual day (for growth respiration)
         alfa_awood = ca2 - ca1 !it calculates the difference btween the biomass of the last and actual day (for growth respiration)  
         alfa_bwood = cb2 - cb1 !it calculates the difference btween the biomass of the last and actual day (for growth respiration) 
         alfa_froot = cf2 - cf1 !it calculates the difference btween the biomass of the last and actual day (for growth respiration) 
         
         
c     maximum evapotranspiration (emax)
         call evpot2 (p0,temp,rh,ae,emax)
c     
c     snow budget
         smelt = 2.63 + 2.55*temp + 0.0912*temp*prain !snowmelt (mm/day) [Eq. 4]
         smelt = amax1(smelt,0.)
         smelt = amin1(smelt,s+psnow)
         ds = psnow - smelt     ![Eq. 2]
         s = s + ds
c     
c     water budget
         if (tsoil.le.tice) then !frozen soil
            g = g + w           !soil moisture freezes
            w = 0.0
            roff = smelt + prain
            evap = 0.0
            ph = 0.0
            ar = 0.0
            nppa = 0.0
            laia = 0.0
            cl = 0.0
            cs = 0.0
            hr = 0.0
            rc2 = 100.0         !default value, equal to aerodynamic resistance (below)
         else                   !non-frozen soil
            w = w + g           !soil ice melts
            g = 0.0
            rimelt = 0.0
            if (w.gt.wmax) then
               rimelt = w - wmax !runoff due to soil ice melting
               w = wmax
            endif
c     
c     Canopy resistance (based in Sellers et al. 1996; SiB2)
c     (rc2 ; s/m) [Eq. 32]
c     [NPP*2.64e-6 converts kgC/m2/yr to molCO2/m2/s]
c     [p0*100 convertes hPa (mb) to Pa]
c     
            nppb = amax1(nppa,0.05)
            rc2 = (ca/(9.0*(nppb*2.64e-5)*0.685*(p0*100)))
            call runoff (w,wmax,roff) !soil moisture runoff (roff, mm/day) [Eq. 10]
            call penman (p0,temp,w,wmax,rh,ae,rc2,evap) !actual evapotranspiration (evap, mm/day)
            dw = prain + smelt - evap - roff ![Eq. 1]
            w = w + dw
            if (w.gt.wmax) then
               roff = roff + (w - wmax)
               w = wmax
            endif
            if (w.lt.0.) w = 0.
            roff = roff + rimelt !total runoff
c     carbon cycle (Microbial respiration, litter and soil carbon)
            call carbon2 (ts,f5,evap,laia, !input
     &           cl,cs,hr)      !output
         endif
c     
c     updating monthly values
         if (pft.eq.1) epavg = epavg + emax 
         smavg = smavg + smelt
         ruavg = ruavg + roff
         evavg = evavg + evap
         rcavg = rcavg + rc2

         phavg = phavg + ph/365.0 !kgC/m2
         aravg = aravg + ar/365.0 !kgC/m2
         nppavg = nppavg + nppa/365.0 !kgC/m2
         laiavg = laiavg + laia/365.0
         clavg = clavg + cl/365.0
         csavg = csavg + cs/365.0
         hravg = hravg + hr/365.0 !kgC/m2
         
         rmlavg = rmlavg + rml/365.0
         rmfavg = rmfavg + rmf/365.0
         rmsavg = rmsavg + rms/365.0
         rmavg = rmavg + rm/365.0
         rglavg = rglavg + rgl/365.0
         rgfavg = rgfavg + rgf/365.0
         rgsavg = rgsavg + rgs/365.0
         rgavg = rgavg + rg/365.0  	  
         
         cleafavg = cleafavg + cl2 !KgC/m2
         cawoodavg = cawoodavg + ca2
         cfrootavg = cfrootavg + cf2
         cbwoodavg = cbwoodavg + cb2
         cstoavg = cstoavg + cs2
         crepavg = crepavg + cr2
         cotheravg = cotheravg + co2 
      enddo


      
!     atualizacao de variaveis

      
      w2 = w
      g2 = g
      s2 = s
      cl2_pft = cl2 !guarda o valor final para o C pool (kg m-2)
      ca2_pft = ca2
      cf2_pft = cf2
      cb2_pft = cb2
      cs2_pft = cs2
      cr2_pft = cr2
      co2_pft = co2

      
      cleafavg = cleafavg/real(ndmonth(month)) !monthly carbon content on leaf compart
      cawoodavg = cawoodavg/real(ndmonth(month))
      cfrootavg = cfrootavg/real(ndmonth(month))
      cbwoodavg = cbwoodavg/real(ndmonth(month))
      cstoavg = cstoavg/real(ndmonth(month))
      crepavg = crepavg/real(ndmonth(month))
      cotheravg = cotheravg/real(ndmonth(month))

!     atualizando outputs de budget
      if(pft .eq. 1) epavg = epavg/real(ndmonth(month))
      smavg = smavg/real(ndmonth(month))
      ruavg = ruavg/real(ndmonth(month))
      evavg = evavg/real(ndmonth(month))
      rcavg = rcavg/real(ndmonth(month))
      phavg = phavg*12.0        !kgC/m2/yr
      aravg = aravg*12.0        !kgC/m2/yr
      rmlavg = rmlavg*12.0
      rmfavg = rmfavg*12.0
      rmsavg = rmsavg*12.0
      rmavg = rmavg*12.0
      rglavg = rglavg*12.0
      rgfavg = rgfavg*12.0
      rgsavg = rgsavg*12.0
      rgavg = rgavg*12.0
      nppavg = nppavg*12.0      !kgC/m2/yr
      laiavg = laiavg*12.0
      clavg = clavg*12.0        !kgC/m2
      csavg = csavg*12.0        !kgC/m2
      hravg = hravg*12.0        !kgC/m2/yr

      return
      
      end subroutine budget
      
      
c     
c======================================================================
c     
      subroutine penman (spre,temp,w,wmax,ur,rn,rc2,evap)
c     
c     Entradas
c     --------
c     spre   = pressao aa supeficie (mb)
c     temp   = temperatura (oC)
c     w      = grau de saturacao (0-1,adimensional)
c     ur     = umidade relativa  (0-1,adimensional)
c     rn     = saldo de radiacao (W m-2)
c     rc2    = resistencia do dossel (s/m)
c     
c     Saida
c     -----
c     evap  = evapotranspiracao (mm/dia)
c     
      real spre,temp,w,wmax,ur,rn,rc2,evap
c     
c     parametros
      ra =    100.              !s/m
      h5    = 0.0275            !mb-1
c     
c     delta
      t1 = temp + 1.
      t2 = temp - 1.
      call tetens(t1,es1)
      call tetens(t2,es2)
      delta = (es1-es2)/(t1-t2) !mb/oC
c     
c     delta_e
      call tetens (temp,es)
      delta_e = es*(1. - ur)    !mb
c     
      if ((delta_e.ge.(1./h5)-0.5).or.(rc2.ge.4500)) evap = 0.
      if ((delta_e.lt.(1./h5)-0.5).or.(rc2.lt.4500)) then
c     gama e gama2
         gama  = spre*(1004.)/(2.45e6*0.622)
         gama2 = gama*(ra + rc2)/ra
c     evapotranspiracao real
         evap = (delta*rn + (1.20*1004./ra)*delta_e)/(delta+gama2) ! W/m2
         evap = evap*(86400./2.45e6) ! mm/dia
         evap = amax1(evap,0.)  !elimina condensacao
      endif
c     
      return
      end
c     
c======================================================================
c     
      subroutine evpot2 (spre,temp,ur,rn,evap)
c     
c     Entradas
c     --------
c     spre   = pressao aa supeficie (mb)
c     temp   = temperatura (oC)
c     ur     = umidade relativa  (0-1,adimensional)
c     rn     = saldo de radiacao (W m-2)
c     
c     Saida
c     -----
c     evap  = evapotranspiracao potencial sem estresse (mm/dia)
c     
      real spre,temp,ur,rn,evap
c     
c     parametros
      ra =    100.              !s/m
      rcmin = 100.              !s/m
c     
c     delta
      t1 = temp + 1.
      t2 = temp - 1.
      call tetens(t1,es1)
      call tetens(t2,es2)
      delta = (es1-es2)/(t1-t2) !mb/oC
c     
c     delta_e
      call tetens (temp,es)
      delta_e = es*(1. - ur)    !mb
c     
c     resistencia estomatica
      rc = rcmin
c     
c     gama e gama2
      gama  = spre*(1004.)/(2.45e6*0.622)
      gama2 = gama*(ra + rc)/ra
c     
c     evapotranspiracao potencial sem estresse
      evap = (delta*rn + (1.20*1004./ra)*delta_e)/(delta+gama2) ! W/m2
      evap = evap*(86400./2.45e6) ! mm/dia
      evap = amax1(evap,0.)     !elimina condensacao
c     
      
      return
      end
c     
c=====================================================================
c     
      subroutine runoff (w,wmax,roff)
      real w,roff
c     roff = 38.*((w/wmax)**11.) ! [Eq. 10]
      roff = 11.5*((w/wmax)**6.6) !from NCEP-NCAR Reanalysis data 
      return
      end
c     
c=====================================================================
c     
      subroutine tetens (t,es)
      real t,es
      if (t.ge.0.) then
         es = 6.1078*exp((7.5*t/(237.3+t))*log(10.))
      else
         es = 6.1078*exp((9.5*t/(265.5+t))*log(10.))
      endif
      return
      end
!     NEW SUBROUTINE SPINUP
      

c======================================================================
      subroutine spinup(nppot,
     &     cleafini,cfrootini,cawoodini,
     &     cbwoodini,cstoini,cotherini,crepini) 

      IMPLICIT NONE

      integer, parameter :: nt=5000
      integer, parameter :: npft=3
c     inputs
      integer i6, kk, k
      
      real :: nppot
      real :: sensitivity,sensitivity2

      
c     outputs
      real :: cleafini(npft)
      real :: cawoodini(npft)
      real :: cbwoodini(npft)
      real :: cfrootini(npft)
      real :: cstoini(npft)
      real :: cotherini(npft)
      real :: crepini(npft)

c     internal vars

      real cleafi_aux(nt)
      real cfrooti_aux(nt)
      real cawoodi_aux(nt)
      real cbwoodi_aux(nt)
      real cstoi_aux(nt)
      real cotheri_aux(nt)
      real crepi_aux(nt)

      
    
      real aleaf(3)             !npp percentage alocated to leaf compartment
      data aleaf /0.35,0.35,0.45/
      real aawood (3)           !npp percentage alocated to aboveground woody biomass compartment
      data aawood /0.40,0.40,0.001/
      real afroot(3)            !npp percentage alocated to fine roots compartment
      data afroot /0.25,0.25,0.55/ 
      real abwood(3)            !npp percentage alocated to belowground woody biomass compartment
      data abwood /0.10,0.10,0.001/
      real asto(3)              !npp percentage alocated to storage compartment
      data asto /0.10,0.10,0.10/
      real arep(3)              !npp percentage alocated to reproduction compartment
      data arep /0.15,0.15,0.10/
      real aother(3)            !npp percentage alocated to other compartment
      data aother /0.05,0.05,0.06/ 
 
      real tleaf(3)             !turnover time of the leaf compartment (yr)
      data tleaf /1.0,0.5,1.0/ 
      real tawood (3)           !turnover time of the aboveground woody biomass compartment (yr)
      data tawood /30.0,30.0,30.0/
      real tfroot(3)            !turnover time of the fine roots compartment
      data tfroot /1.0,1.0,1.0/
      real tbwood (3)           !turnover time of the belowground woody biomass compartment
      data tbwood /40.0,40.0,40.0/
      real tsto  (3)            !turnover time of the storage compartmentturn
      data tsto /5.0,5.0,5.0/ 
      real trep (3)             !turnover time of the reproduction compartment
      data trep /0.25,0.25,0.25/ 
      real tother (3)           !turnover time of the other compartment
      data tother /0.12,0.12,0.12/

      sensitivity = 1.10
      sensitivity2 = 1.40

      do i6=1,npft
         do k=1,nt
            if (k.eq.1) then
               cleafi_aux(k) = aleaf(i6)*(nppot)
               cawoodi_aux(k) = aawood(i6)*(nppot)
               cfrooti_aux(k) = afroot(i6)*(nppot)
               cbwoodi_aux(k) = abwood(i6)*(nppot)
               cstoi_aux(k) = asto(i6)*(nppot)
               cotheri_aux(k) = aother(i6)*(nppot/365)
               crepi_aux(k) = arep(i6)*(nppot/365)
            else
               cleafi_aux(k) = ((aleaf(i6)*(nppot))-
     &              (cleafi_aux(k-1)/(tleaf(i6)))) + cleafi_aux(k-1)
               cawoodi_aux(k) = ((aawood(i6)*(nppot))-
     &              (cawoodi_aux(k-1)/(tawood(i6)))) + cawoodi_aux(k-1)
               cfrooti_aux(k) = ((afroot(i6)*(nppot))-
     &              (cfrooti_aux(k-1)/(tfroot(i6)))) + cfrooti_aux(k-1)
               cbwoodi_aux(k) = ((abwood(i6)*(nppot))-
     &              (cbwoodi_aux(k-1)/(tbwood(i6)))) + cbwoodi_aux(k-1)
               cstoi_aux(k) = ((asto(i6)*(nppot))-
     &              (cstoi_aux(k-1)/(tsto(i6)))) + cstoi_aux(k-1)
               cotheri_aux(k) = ((aother(i6)*(nppot))-
     &              (cotheri_aux(k-1)/(tother(i6)*365))) + cotheri_aux(k
     $              -1)
               crepi_aux(k) = ((arep(i6)*(nppot))-
     &              (crepi_aux(k-1)/(trep(i6)*365))) + crepi_aux(k-1)
               
               kk =  int(k*0.66)
               if((cfrooti_aux(k)/cfrooti_aux(kk).lt.sensitivity)
     $              .and.(cleafi_aux(k)/cleafi_aux(kk).lt.sensitivity)
     $              .and.(cawoodi_aux(k)/cawoodi_aux(kk).lt.
     $              sensitivity2).and.(cbwoodi_aux(k)/cbwoodi_aux(kk)
     $              .lt.sensitivity2).and.(cstoi_aux(k)
     $              /cstoi_aux(kk).lt.sensitivity).and.(cotheri_aux(k)
     $              /cotheri_aux(kk).lt.sensitivity).and.(crepi_aux(k)
     $              /crepi_aux(kk).lt.sensitivity))   then
                  
                  
                  cleafini(i6) = cleafi_aux(k) ! carbon content (kg m-2) 
                  cawoodini(i6) = cawoodi_aux(k)
                  cfrootini(i6) = cfrooti_aux(k)
                  cbwoodini(i6) = cbwoodi_aux(k)
                  cstoini(i6) = cstoi_aux(k)
                  cotherini(i6) = cotheri_aux(k)
                  crepini(i6) = crepi_aux(k)
                  exit
                  
               endif
            endif
         enddo
      enddo
      
      
      return
      end subroutine spinup
