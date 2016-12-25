C23456
      subroutine wbm (prec,temp,lsmk,p0,ca,par, 
     &     cleaf_ini, cawood_ini, cfroot_ini   
     &     emaxm, tsoil, photo_pft,aresp_pft,npp_pft,lai_pft,
     &     clit_pft,csoil_pft, hresp_pft,rcm_pft,runom_pft,
     &     evapm_pft,wsoil_pft,rml_pft,rmf_pft,rms_pft,rm_pft,rgl_pft
     $     ,rgf_pft,rgs_pft,rg_pft,cleaf_pft,
                  cawood_pft,cfroot_pft)
!     implicit none


!     =================================================================
!     Water balance model (WBM).
!     From monthly climatologies of precipitation and surface
!     temperature, 
!     the WBM calculates the environmental variables.
!     
!     05 Jul 2005, MDO: ae, rh & runoff are changed
!     11 Jul 2005, MDO: wsoil2 is written (for testing purpose)
!     31 Ago 2006, DML: carbon cycle is included
C     MODIFIED BY HAP AND JPDF 15-12-2016 -- ADAPTING TO CAETE TROIA 
!     =================================================================
!     
!     Variables
!     =========
!     
      integer,parameter :: nx=720,ny=360,q=3
      real,parameter :: no_data = -9999.0
      
!     
c     --------------------------I N P U T S----------------------------
      real ca
      real p0(nx,ny,12)         !Atmospheric pressure (mb)
      real lsmk(nx,ny)          !Land=1/Ocean=0
      real prec(nx,ny,12)       !Precipitation (mm/month)
      real temp(nx,ny,12)       !Temperature (oC)
      real par(nx,ny,12)        !IPAR (Ein/m2/s)
C     -----------------------------E N D-------------------------------
      
c     -------------------------O U T P U T S---------------------------
      
c     primeiro algumas variaveis que nao sao dependentes de pfts
      real emaxm(nx,ny,12)
      real tsoil(nx,ny,12)
      
c     agora as variaveis para pfts
      real photo_pft(nx,ny,12,q) !Monthly photosynthesis   (kgC/m2)
      real aresp_pft(nx,ny,12,q) !Monthly autotrophic res  (kgC/m2)
      real npp_pft(nx,ny,12,q)  !Monthly net primary produ (kgC/m2)
      
      real lai_pft(nx,ny,12,q)  !Monthly leaf area index
      real clit_pft(nx,ny,12,q) !Monthly litter carbon
      real csoil_pft(nx,ny,12,q) !Monthly soil carbon
      real hresp_pft(nx,ny,12,q) !Monthly het resp          (kgC/m2)
      real rcm_pft(nx,ny,12,q) 
      
c     VARIAVEIS HIDROLOGICAS IMPORTANTES   
      real runom_pft(nx,ny,12,q) !Runoff
      real evapm_pft(nx,ny,12,q) !Actual evapotranspiration        
      real wsoil_pft(nx,ny,12,q) !Soil moisture (mm)


      real rml_pft(nx,ny,12,q)
      real rmf_pft(nx,ny,12,q)
      real rms_pft(nx,ny,12,q)
      real rm_pft(nx,ny,12,q)
      
      real rgl_pft(nx,ny,12,q)
      real rgf_pft(nx,ny,12,q)
      real rgs_pft(nx,ny,12,q)
      real rg_pft(nx,ny,12,q)
      
      real cleaf_pft(nx,ny,12,q) !monthly npp alloc to leaf biomass (KgC/m2)
      real cawood_pft(nx,ny,12,q) !monthly npp alloc to aboveground biomass (KgC/m2)
      real cfroot_pft(nx,ny,12,q) !monthly npp alloc to fine root

      
c --------------------------------E N D--------------------------------

c     VARIAVEIS HIDROLOGICAS IMPORTANTES   
      real gsoil(nx,ny,12,q)      !Soil ice 
      real ssoil(nx,ny,12,q)      !Soil snow
      real snowm(nx,ny,12,q)      !Snowmelt
      real wg0(nx,ny,12,q)        !Moisture of the previous year

      integer i, j, k, kk, mes, nerro, p

      real, parameter :: H = 1.0 !Soil layer(m) 
      real, parameter :: diffu = 4.e-7*(30.*86400.0) !Soil thermal diffusivity (m2/month)
      real, parameter :: tau = (H**2)/(2.0*diffu) !E-folding time (months)
      real, parameter :: auxs = -100.0 !Auxiliar for calculation of Snpp
      real t0,t1

      
      real wsaux1
      real ae,dwww,wmax
      real sini,gfim,sfim,gini,wfim,wini
      real pr,ice,spre,ta,td,ipar
      real rmes,phmes,smes,rcmes,hrmes,nppmes,laimes
      real armes,clmes,csmes,emes,epmes,rmlmes,rmfmes, rmsmes, rmmes
      real rglmes, rgfmes,rgsmes,rgmes, cleafmes, cawoodmes,cfrootmes
      
      
!     Soil temperature
!     ================
!     For all grid points
!     -------------------
c234567
!     
      do i=1,nx
         do j=1,ny
!     
!     Initialize soil temperature
!     ---------------------------
!     
            do k=1,12
               tsoil(i,j,k) = no_data
            enddo
!     
!     Only for land grid points
!     -------------------------
!     
            if (int(lsmk(i,j)).ne.0) then
               t0 = 0.          !Initialization
               do n = 1,1200    !1200 months run to attain equilibrium
                  k = mod(n,12)
                  if (k.eq.0) k = 12
                  t1 = t0*exp(-1.0/tau)+(1.0-exp(-1.0/tau))*temp(i,j,k) 
                  tsoil(i,j,k) = (t0 + t1)/2.0
                  t0 = t1
               enddo
            endif
!     
         enddo
      enddo

c     agora nos ja temos a tsoil para todas as celulas de grid para todos os meses 
!     
c23456
!     ============
!     Water budget
!     ============
!     
!     For all grid points
!     -------------------
!     
      do i=1,nx
         do j=1,ny     
!     Write to track program execution     
            if ((mod(j,ny).eq.0).and.(mod(i,10).eq.0))
     &           write(*,*) 'working:',i

!     Initialize variables

            do k=1,12
               
!     atencao com emaxm... nao precisa calcular para totos os pfts ..
               
               emaxm(i,j,k)  = no_data !Maximum evapotranspiration
               
               do p=1,q
                  
                  photo_pft(i,j,k,p)  = no_data !Monthly photosynthesis (kgC/m2)
                  aresp_pft(i,j,k,p)  = no_data !Monthly autotrophic respiration (kgC/m2)
                  npp_pft(i,j,k,p)  = no_data !Monthly net primary productivity (average between PFTs) (kgC/m2)
                  
                  lai_pft(i,j,k,p)  = no_data !Monthly leaf area index
                  clit_pft(i,j,k,p)  = no_data !Monthly litter carbon
                  csoil_pft(i,j,k,p)  = no_data !Monthly soil carbon
                  hresp_pft(i,j,k,p)  = no_data !Monthly heterotrophic respiration (kgC/m2)
                  
                  rcm_pft(i,j,k,p)  = no_data
                  wg0(i,j,k,p)    = no_data
                  gsoil(i,j,k,p)  = no_data !Soil ice
                  ssoil(i,j,k,p)  = no_data !Soil snow
                  
c     VARIAVEIS HIDROLOGICAS IMPORTANTES   
                  runom_pft(i,j,k,p)  = no_data !Runoff
                  evapm_pft(i,j,k,p)  = no_data !Actual evapotranspiration        
                  wsoil_pft(i,j,k,p)  = no_data !Soil moisture (mm)
                  
C     NOVAS VARIAVEIS VEG POOLS  tem que declarar
                  rml_pft(i,j,k,p) = no_data
                  rmf_pft(i,j,k,p) = no_data
                  rms_pft(i,j,k,p) = no_data
                  rm_pft(i,j,k,p)  = no_data
                  
                  rgl_pft(i,j,k,p) = no_data
                  rgf_pft(i,j,k,p) = no_data
                  rgs_pft(i,j,k,p) = no_data
                  rg_pft(i,j,k,p)  = no_data
                  
                  cleaf_pft(i,j,k,p)= no_data !monthly npp alloc to leaf biomass (KgC/m2)
                  cawood_pft(i,j,k,p)= no_data !monthly npp alloc to aboveground biomass (KgC/m2)
                  cfroot_pft(i,j,k,p)= no_data !monthly npp alloc to fine root biomass (KgC/m2)

               enddo              
            enddo
!     
!     Only for land grid points
!     -------------------------
!     
            if (nint(lsmk(i,j)).ne.0) then
               do k=1,12
                  spre = p0(i,j,k) !Surface pressure (mbar)==kPa
               enddo
               
!     agora nos vamos aplicar esta integracao q vezes dentro da wbm ... uma para cada pft
               do p = 1,q       ! PARA CADA PFT
                  
                  
!     
!     Set some variables
!     ------------------
!     
                  wini  = 0.01  !Soil moisture_initial condition (mm)
                  gini  = 0.0   !Soil ice_initial condition (mm)
                  sini  = 0.0   !Overland snow_initial condition (mm)
                  
                  cleaf1_pft  =  cleafini_pft(i,j,p) !inital leaf biomass for each PFT from spinup (KgC/m2) 
                  cawood1_pft = cawoodini_pft(i,j,p)
                  cfroot1_pft = cfrootini_pft(i,j,p)
!     Initialization
!     --------------
!     
                  do k=1,12
                     wg0(i,j,k,p) = -1.0               
                  enddo
!     
!     Start integration
!     -----------------
!     
                  
                  
                  
                  n = 0
 10               continue
                  n = n + 1
!     
!     Pre-processing
!     --------------
!     
                  k = mod(n,12)
                  if (k.eq.0) k = 12
                  mes = k
                  td = tsoil(i,j,k)
                  ta = temp(i,j,k)
                  pr = prec(i,j,k)
                  ipar = par(i,j,k)
                  ae = 2.895*ta+52.326 !Available energy (W/m2) - From NCEP-NCAR reanalysis data
!     
!     Monthly water budget
!     ====================

                  call budget (p,mes,wini,gini,sini,td,ta,pr,spre,ae,ca,
     &                 ipar,cleaf1_pft, cawood1_pft, cfroot1_pft,wfim
     $                 ,gfim, sfim, clfim, cafim,crfim,smes,rmes,emes
     $                 ,epmes,phmes,armes,nppmes,laimes,clmes,csmes
     $                 ,hrmes,rcmes,rmlmes,rmfmes, rmsmes, rmmes, rglmes
     $                 , rgfmes,rgsmes,rgmes, cleafmes, cawoodmes,
     $                 cfrootmes) 
!     
!     Update variables
!     ----------------
!     
                  
                  if (p.eq.1) emaxm(i,j,k) = epmes
                  gsoil(i,j,k,p) = gfim
                  ssoil(i,j,k,p) = sfim
                  wsoil_pft(i,j,k,p) = wfim
                  snowm(i,j,k,p) = smes
                  runom_pft(i,j,k,p) = rmes
                  evapm_pft(i,j,k,p) = emes
                  
                  rcm_pft(i,j,k,p)   = rcmes
                  lai_pft(i,j,k,p)   = laimes
                  photo_pft(i,j,k,p) = phmes
                  aresp_pft(i,j,k,p) = armes
                  npp_pft(i,j,k,p)   = nppmes
                  clit_pft(i,j,k,p)  = clmes
                  csoil_pft(i,j,k,p) = csmes
                  hresp_pft(i,j,k,p) = hrmes

                  rml_pft(i,j,k,p) = rmlmes
                  rmf_pft(i,j,k,p) = rmfmes
                  rms_pft(i,j,k,p) = rmsmes
                  rm_pft(i,j,k,p)  = rmmes
                  
                  rgl_pft(i,j,k,p) = rglmes
                  rgf_pft(i,j,k,p) = rgfmes
                  rgs_pft(i,j,k,p) = rgsmes
                  rg_pft(i,j,k,p)  = rgmes
                  
                  cleaf_pft_pft(i,j,k,p)  = cleafmes
                  cawood_pft_pf (i,j,k,p) = cawoodmes
                  cfroot_pft_pft(i,j,k,p) = cfrootmes
                  
                  
                  wini         = wfim
                  gini         = gfim
                  sini         = sfim
                  cleaf1_pft  =  clfim 
                  cawood1_pft = cafim
                  cfroot1_pft = cfim
                        
!     Check if equilibrium is attained
!     --------------------------------
!     
                  if (k.eq.12) then
                     wmax = 500.
                     nerro = 0
                     do kk=1,12
                        wsaux1 = wsoil_pft(i,j,kk,p) + gsoil(i,j,kk,p)
                        dwww = (wsaux1 - wg0(i,j,kk,p)) / wmax
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
                  goto 10
 100              continue
               enddo            !p loop     
            endif               ! endif lsmk
c     finalize ny loop
         enddo                  ! j
c     finalize nx loop
      enddo                     ! k
         
c     loop dos pfts encerrado
c     daqui por diante nada util --- vai tudo para env5.f
      return
      end subroutine wbm
      
      
      
      subroutine budget (pft,month,w1,g1,s1,ts,temp,prec,p0,ae,ca,ipar
     $     ,cl1_pft,ca1_pft,cf1_pft,w2,g2,s2,cl2_pft,ca2_pft,cf2_pft
     $     ,smavg,ruavg,evavg,epavg,phavg,aravg,nppavg,laiavg,clavg
     $     ,csavg,hravg,rcavg,rmlavg,rmfavg,rmsavg,rmavg,rglavg,rgfavg
     $     ,rgsavg,rgavg,cleafavg_pft,cawoodavg_pft,cfrootavg_pft)
      
!     Surface water             !soil moisture, snow and ice budget for a single month
!     
!     =============
!     REPARE QUE VC TEM UMA NOVA VARIAVEL PARA A BUDGET ---> pft
!     INPUTS
      
      integer month, pft             !Actual month (1-12)
      real w1                   !Initial (previous month last day) soil moisture storage (mm)
      real g1                   !Initial soil ice storage (mm)
      real s1                   !Initial overland snow storage (mm)
      real cl1_pft              ! initial npp allocation to cleaf compartment
      real cf1_pft              !                           froot
      real ca1_pft              !                           cawood
      real ts                   !Soil temperature (oC)
      real temp                 !Surface air temperature (oC)
      real prec                 !Precipitation (mm/day)
      real p0                   !Surface pressure (mb)
      real ae                   !Available energy (W/m2)
      real ca                   !Atmospheric carbon
      real ipar                 !Incident photosynthetic active radiation
!     
!     Output
!     ------
!     
      real w2                   !Final (last day) soil moisture storage (mm)
      real g2                   !Final soil ice storage (mm)
      real s2                   !Final overland snow storage (mm)
      real cl2_pft              !FINAL npp allocation to veg_pool
      real cf2_pft
      real ca2_pft
      real smavg                !Snowmelt monthly average (mm/day)
      real ruavg                !Runoff monthly average (mm/day)
      real evavg                !Actual evapotranspiration monthly average (mm/day)
      real epavg                !Maximum evapotranspiration monthly average (mm/day)
      real phavg                !Monthly photosynthesis
      real aravg                !Monthly autotrophic respiration
      real nppavg               !Monthly NPP (average between PFTs)
      real laiavg               !Monthly leaf area Index
      real clavg                !Monthly carbon litter
      real csavg                !Monthly carbon soil
      real hravg                !Monthly heterotrophic respiration
      real rcavg                !Monthly canopy resistence

      real rmlavg,rmfavg,rmsavg,rmavg,rglavg,rgfavg
      real rgsavg,rgavg,cleafavg_pft,cawoodavg_pft,cfrootavg_pft
      
!     Internal Variables

      real rh                   !Relative humidity
      real wmax                 !Soil moisture availability (mm)
      real tsnow                !Temperature threshold for snowfall (oC)
      real tice                 !Temperature threshold for soil freezing (oC)
      real psnow                !Snowfall (mm/day)
      real prain                !Rainfall (mm/day)
      real w                    !Daily soil moisture storage (mm)
      real g                    !Daily soil ice storage (mm)
      real s                    !Daily overland snow storage (mm)
      real rimelt               !Runoff due to soil ice melting
      real smelt                !Snowmelt (mm/day)
      real roff                 !Total runoff
      real evap                 !Actual evapotranspiration (mm/day)
      real emax                 !Maximum evapotranspiration
      integer ndmonth(12)       !Number of months
      data ndmonth /31,28,31,30,31,30,31,31,30,31,30,31/ !Number of days for each month
!     
!     Carbon Cycle
!     

      real ph                   !Canopy gross photosynthesis (kgC/m2/yr)
      real ar                   !Autotrophic respiration (kgC/m2/yr)
      real nppa,nppb            !Net primary productivity / auxiliar
      real laia                 !Leaf area index (m2 leaf/m2 area)
      real cl                   !Litter carbon (kgC/m2) ---- anual?
      real cs                   !Soil carbon (kgC/m2) ---- anual?
      real hr                   !Heterotrophic (microbial) respiration (kgC/m2/yr)
      real mgama                !Photo-respiration compensation point (Pa)
      real jc,je,jl,jp          !Factors that limit photosynthesis
      real f4sun,f4shade        !Scaling-up to canopy level(dimensionless)
      real gs                   !Stomatal conductance (mol/m2/s)
      real rc2                  !Canopy resistence (s/m)
      real f1                   !Photosynthesis (mol/m2/s)
      real f1b                  !Photosynthesis (micromol/m2/s)
      real rm,rml,rmf,rms,rg,rgl,rgf,rgs ! maintenance and growth plant respiration
      real cl1, cf1, ca1, cl2,cf2, ca2
!     Initialize Canopy Resistence Parameters
!     ---------------------------------------
!     
      gs  = 0.0
      rc2 = 0.0
      f1  = 0.0
      f1b = 0.0
!     
!     Parameters
!     ----------
!     
      rh    = 0.685
      wmax  = 500.0
      tsnow = -1.0
      tice  = -2.5
!     
!     Precipitation
!     =============
!     
      psnow = 0.0
      prain = 0.0
      if (temp.lt.tsnow) then
         psnow = prec/real(ndmonth(month))
      else
         prain = prec/real(ndmonth(month))
      endif
!     
!     Initialization
!     --------------
!     
      w           = w1
      g           = g1
      s           = s1
      smavg       = 0.0
      ruavg       = 0.0
      evavg       = 0.0
      epavg       = 0.0
      rcavg       = 0.0
      laiavg      = 0.0
      phavg       = 0.0
      aravg       = 0.0
      nppavg      = 0.0
      clavg       = 0.0
      csavg       = 0.0
      hravg       = 0.0
      rmlavg       = 0.0
      rmfavg       = 0.0
      rmsavg       = 0.0
      rmavg       = 0.0

      rglavg       = 0.0
      rgfavg       = 0.0
      rgsavg       = 0.0
      rgavg       = 0.0

      cleafavg_pft = 0.0         !mean monthly leaf biomass for all the PFTs
      cawoodavg_pft = 0.0        !mean monthly aboveground biomass for all the PFTs
      cfrootavg_pft = 0.0        !mean monthly belowground biomass for all the PFTs
      
      
      cl2_pft = 0.0
      ca2_pft = 0.0
      cf2_pft = 0.0
      
      alfa_leaf = 0.0
      alfa_awood = 0.0
      alfa_froot = 0.0

!     Numerical integration
!     ---------------------
!     
      do i=1,ndmonth(month)
!     
         nppa      = 0.0        !Auxiliar_nppa
         ph        = 0.0        !Auxiliar_ph
         ar        = 0.0        !Auxiliar_ar
         laia      = 0.0        !Auxiliar_laia
         f5        = 0.0        !Auxiliar_f5
         f1        = 0.0        !Auxiliar_f1
         vpd       = 0.0
         rc2       = 0.0        !Auxiliar_rc2
         
         rm   = 0.0
         rml  = 0.0
         rmf  = 0.0
         rms  = 0.0
         rg   = 0.0
         rgl  = 0.0
         rgf  = 0.0
         rgs  = 0.0

         cl1 = cl2_pft
         ca1 = ca2_pft
         cf1 = cf2_pft
         
         
         beta_leaf = alfa_leaf
         beta_awood = alfa_awood
         beta_froot = alfa_froot

         
         if ((i.eq.1).and.(month.eq.1)) then
            cl1 = cl1_pft
            ca1 = ca1_pft
            cf1 = cf1_pft
            
            beta_leaf=0.
            beta_awood=0.
            beta_froot=0.

         endif
!     Carbon cycle (photosynthesis, plant respiration and NPP)

         call productivity1 (pft, temp, p0, w, wmax, ca, ipar, tsoil,
     $        cl1, ca1, cf1, beta_leaf, beta_awood, beta_froot, ph,
     $        ar, nppa,laia, f5, f1, vpd, rm, rml, rmf, rms, rg, rgl,
     $        rgf, rgs)
!     
!     carbon allocation (carbon content on each compartment)
         call allocation (pft, nppa, cl1, ca1, cf1, !input
     &        cl2, ca2, cf2)      !output
         
         alfa_leaf  = cl2 - cl1 
         alfa_awood = ca2 - ca1 
         alfa_froot = cf2 - cf1 
         
         cl2_pft = cl2
         ca2_pft = ca2
         cf2_pft = cf2
         
!     Maximum evapotranspiration (emax)
!     =================================
         if (pft.eq.1) then
            call evpot2 (p0,temp,rh,ae,emax)
         endif
!     Snow budget
!     ===========
!     
         smelt = 2.63 + 2.55*temp + 0.0912*temp*prain !Snowmelt (mm/day)
         smelt = amax1(smelt,0.)
         smelt = amin1(smelt,s+psnow)
         ds = psnow - smelt
         s = s + ds
!     
!     Water budget
!     ============
!     
         if (ts.le.tice) then   !Frozen soil
            g = g + w           !Soil moisture freezes
            w = 0.0
            roff = smelt + prain !mm/day
            evap = 0.0
            ph = 0.0
            ar = 0.0
            nppa = 0.0
            laia = 0.0
            cl = 0.0
            cs = 0.0
            hr = 0.0
c     rc2 = 100.0               !Default value, equal to aerodynamic resistance (below)
c     
         else                   !Non-frozen soil
            w = w + g           !Soil ice melts
            g = 0.0
            rimelt = 0.0
            if (w.gt.wmax) then
               rimelt = w - wmax !Runoff due to soil ice melting
               w = wmax
            endif
!     
!     Canopy resistance (Based on Medlyn et al.(2011)
!     (rc2 ; s/m)
!     
            call canopy_resistence (pft,vpd,f1,rc2)     
            call runoff (w,wmax,roff) !Soil moisture runoff (roff, mm/day)
            call penman (p0,temp,rh,ae,rc2,evap) !Actual evapotranspiration (evap, mm/day)
            dw = prain + smelt - evap - roff
            w = w + dw
            if (w.gt.wmax) then
               roff = roff + (w - wmax)
               w = wmax
            endif
            if (w.lt.0.) w = 0.
            roff = roff + rimelt !Total runoff
!     
!     Carbon cycle (Microbial respiration, litter and soil carbon)
!     ============================================================
!     
            call carbon2 (ts,f5,evap,laia, !Inputs
     &           cl,cs,hr)      !Outputs
!     
         endif
!     
!     Updating monthly values
!     =======================
!
         if (pft.eq.1) epavg = epavg + emax   !mm/day
         smavg = smavg + smelt  !mm/day
         ruavg = ruavg + roff   !mm/day
         evavg = evavg + evap   !mm/day
         
         rcavg = rcavg + rc2    !s/m/day
         phavg = phavg + ph/365.0 !kgC/m2/day
         aravg = aravg + ar/365.0 !kgC/m2/day
         nppavg = nppavg + nppa/365.0 !kgC/m2/day
         laiavg = laiavg + laia/365.0 !m2leaf/m2area/day
         clavg = clavg + cl/365.0 !kgC/m2/day
         csavg = csavg + cs/365.0 !kgC/m2/day
         hravg = hravg + hr/365.0 !kgC/m2/day
         rmlavg = rmlavg + rml/365
         rmfavg = rmfavg + rmf/365
         rmsavg = rmsavg + rms/365
         rmavg = rmavg + rm/365
         rglavg = rglavg + rgl/365
         rgfavg = rgfavg + rgf/365
         rgsavg = rgsavg + rgs/365
         rgavg = rgavg + rg/365
         cleafavg_pft = cleafavg_pft + cl2/365
         cawoodavg_pft = cawoodavg_pft + ca2 /365
         cfrootavg_pft = cfrootavg_pft + cf2 /365
         
         
      enddo
!     
!     Final calculations
!     ------------------
!     
      w2 = w
      g2 = g
      s2 = s
      smavg = smavg/real(ndmonth(month))
      ruavg = ruavg/real(ndmonth(month))
      evavg = evavg/real(ndmonth(month))
      if (pft.eq.1) epavg = epavg/real(ndmonth(month))
      rcavg = rcavg/real(ndmonth(month))
      phavg = phavg * 12.0      !kgC/m2/yr
      aravg = aravg * 12.0      !kgC/m2/yr
      nppavg = nppavg * 12.0    !kgC/m2/yr
      laiavg = laiavg * 12.0
      clavg = clavg * 12.0      !kgC/m2
      csavg = csavg * 12.0      !kgC/m2
      hravg = hravg * 12.0      !kgC/m2/yr
      rmlavg = rmlavg * 12.0 
      rmfavg = rmfavg * 12.0
      rmsavg = rmsavg * 12.0
      rmavg = rmavg * 12.0
      rglavg = rglavg * 12.0
      rgfavg = rgfavg * 12.0
      rgsavg = rgsavg * 12.0
      rgavg = rgavg * 12.0
      cleafavg_pft = cleafavg_pft * 12.0
      cawoodavg_pft = cawoodavg_pft * 12.0
      cfrootavg_pft = cfrootavg_pft * 12.0
               
      return
      end subroutine budget
!     
!     =========================================================
!     
      subroutine penman (spre,temp,ur,rn,rc2,evap)
!     
!     Inputs
!     ------
!     
      real spre                 !Surface pressure (mb)
      real temp                 !Temperature (oC)
      real ur                   !Relative humidity (0-1,dimensionless)
      real rn                   !Radiation balance (W/m2)
      real rc2                  !Canopy resistence (s/m)
!     
!     
!     Output
!     ------
!     
      real evap                 !Evapotranspiration (mm/day)
!     
!     Parameters
!     ----------
!     
      ra = 100                  !s/m
      h5 = 0.0275               !mb-1
!     
!     Delta
!     -----
!     
      t1 = temp + 1.
      t2 = temp - 1.
      call tetens(t1,es1)       !Saturation partial pressure of water vapour at temperature T
      call tetens(t2,es2)
      delta = (es1-es2)/(t1-t2) !mb/oC
!     
!     Delta_e
!     -------
!     
      call tetens (temp,es)
      delta_e = es*(1. - ur)    !mb
!     
      if ((delta_e.ge.(1./h5)-0.5).or.(rc2.ge.4500)) evap = 0.
      if ((delta_e.lt.(1./h5)-0.5).or.(rc2.lt.4500)) then
!     
!     Gama and gama2
!     --------------
!     
         gama  = spre*(1004.)/(2.45e6*0.622)
         gama2 = gama*(ra + rc2)/ra
!     
!     Real evapotranspiration
!     -----------------------
!     
         evap = (delta*rn + (1.20*1004./ra)*delta_e)/(delta+gama2) !W/m2
         evap = evap*(86400./2.45e6) !mm/day
         evap = amax1(evap,0.)  !Eliminates condensation
      endif
!     
      return
      end
!     
!     ============================================
!     
      subroutine evpot2 (spre,temp,ur,rn,evap) 
!     
!     Inputs
!     ------
!     
      real spre                 !Surface pressure (mb)
      real temp                 !Temperature (oC)
      real ur                   !Relative humidity (0-1,dimensionless)
      real rn                   !Irradiation balance (W/m2)
!     
!     Output
!     ------
!     
      real evap                 !Potencial evapotranspiration without stress (mm/day)
!     
!     Parameters
!     ----------
!     
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
      end
!     
!     =================================================================
!     ===
!     
      subroutine runoff (w,wmax,roff)
      real w,roff
      roff = 11.5*((w/wmax)**6.6) !From NCEP-NCAR Reanalysis data
      return
      end
!     
!     =================================================================
!     ====
!     
      subroutine tetens (t,es)  !SVP (kpa)!
      real t,es
      if (t.ge.0.) then
         es = 6.1078*exp((7.5*t/(237.3+t))*log(10.))
      else
         es = 6.1078*exp((9.5*t/(265.5+t))*log(10.))
      endif                                                             
!     
      return
      end        
!     =============================================================
