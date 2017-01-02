C23456
      subroutine wbm (prec,temp,lsmk,p0,ca,par, 
     &     cleaf_ini, cawood_ini, cfroot_ini,   
     &     emaxm, tsoil, photo_pft,aresp_pft,npp_pft,lai_pft,
     &     clit_pft,csoil_pft, hresp_pft,rcm_pft,runom_pft,
     &     evapm_pft,wsoil_pft,rml_pft,rmf_pft,rms_pft,rm_pft,rgl_pft
     $     ,rgf_pft,rgs_pft,rg_pft,cleaf_pft,cawood_pft,cfroot_pft)
!     implicit none

c      call wbm (prec,temp,lsmk,p0,ca,par,cleafin, cawoodin,cfrootin,
c     &     emaxm, tsoil, photo_pft,aresp_pft,npp_pft,lai_pft,
c     &     clit_pft,csoil_pft, hresp_pft,rcm_pft,runom_pft,
c     &     evapm_pft,wsoil_pft,rml_pft,rmf_pft,rms_pft,rm_pft,rgl_pft
c     $     ,rgf_pft,rgs_pft,rg_pft,cleaf_pft,cawood_pft,cfroot_pft)   

!     =================================================================
!     Water balance model (WBM).
!     From monthly climatologies of precipitation and surface
!     temperature, 
!     the WBM calculates the environmental variables.
!     
!     05 Jul 2005, MDO: ae, rh & runoff are changed
!     11 Jul 2005, MDO: wsoil2 is written (for testing purpose)
!     31 Ago 2006, DML: carbon cycle is included
!     31 dez 2016, DML, HAP, BR and JP- several changes 
C     MODIFIED BY HAP AND JPDF 15-12-2016 -- ADAPTING TO CAETE TROIA
C       passagem tranquila-JP 
        ! CARBON2 PASSA A SE CHAMAR PRODUCTIVITY1
C     MODIFIED BY BR AND JPDF 27-12-2016 -- ADAPTING TO CAETE BIANCA + TROIA
C       
!     =================================================================
!     ! COMENTARIOS: (27-12-2016)

C     ALGUNS PROCESSOS IMPLEMENTADOS PELA BIANCA (F5 A APARTIR DE RC2 EM CARBON2
C     E LAIA A PARTIR DE SLA) NAO PODEM AINDA SER IMPLEMENTADOS NESTA VERSAO POR PROBLEMAS 
C     PARCIALMENTE DESCONHECIDOS. JP
C     1 - RC2 SO EH CALCULADA DEPÓIS DA CHAMADA DE PRODUCTIVITY1 (ANTIGA CARBON2) 
C     BIA VC ESTAVA USANDO RC2 DENTRO DA PRODUCTIVITY (COMO DENOMINADOR EM UMA DIVISAO)
C     PRECISAMOS REVER A ORDEM DA BUDGET PRA IMPLEMENTAR ISSO.
C     2 - NAO SEI O QUE ESTA HAVENDO COM O CALCULO DE LAIA A PARTIR DA SLA- DA ERRO-NAO SEI   
C     ASSIM ESTES DOIS PROCESSOS (CALCULO DE F5 E LAIA CONTINUAM COMO CARBON2)

      
!     ! COMENTARIOS:


!     Variables
!     =========
!     
      integer,parameter :: nx=720,ny=360,q=3
      real,parameter :: no_data = -9999.0
      
!     
c     --------------------------I N P U T S----------------------------
      real ca                   !CO2 atmospheric concentration (ppmv)
      real p0(nx,ny,12)         !Atmospheric pressure (mb)
      real lsmk(nx,ny)          !Land=1/Ocean=0
      real prec(nx,ny,12)       !Precipitation (mm/month)
      real temp(nx,ny,12)       !Temperature (oC)
      real par(nx,ny,12)        !IPAR (Ein/m2/s)

      real  cleaf_ini(nx,ny,q)  ! Initial carbon content in leaves (kg m-2)
      real  cawood_ini(nx,ny,q) ! Initial carbon content in aboveground wood (kg m-2)
      real  cfroot_ini(nx,ny,q) ! Initial carbon content in fineroots (kg m-2)
C     -----------------------------E N D-------------------------------
      
c     -------------------------O U T P U T S---------------------------
       
      real tsoil (nx,ny,12)     !soil temperature
      
      real photo_pft(nx,ny,12,q) !Monthly photosynthesis   (kgC m-2)
      real aresp_pft(nx,ny,12,q) !Monthly autotrophic res  (kgC m-2)
      real npp_pft  (nx,ny,12,q) !Monthly net primary produ (kgC m-2)
      
      real lai_pft  (nx,ny,12,q) !Monthly leaf area index
      real clit_pft (nx,ny,12,q) !Monthly litter carbon
      real csoil_pft(nx,ny,12,q) !Monthly soil carbon
      real hresp_pft(nx,ny,12,q) !Monthly het resp          (kgC/m2)
      real rcm_pft  (nx,ny,12,q) 

      real emaxm    (nx,ny,12,q) !Maximum evapotranspiration (kg m-2 day-1)
      real runom_pft(nx,ny,12,q) !Runoff 
      real evapm_pft(nx,ny,12,q) !Actual evapotranspiration        
      real wsoil_pft(nx,ny,12,q) !Soil moisture (mm)
      
c     NOVOS OUTPUTS DA BIANCA
      real rml_pft  (nx,ny,12,q)
      real rmf_pft  (nx,ny,12,q)
      real rms_pft  (nx,ny,12,q)
      real rm_pft   (nx,ny,12,q)
      
      real rgl_pft  (nx,ny,12,q)
      real rgf_pft  (nx,ny,12,q)
      real rgs_pft  (nx,ny,12,q)
      real rg_pft   (nx,ny,12,q)
      
      real cleaf_pft (nx,ny,12,q) !monthly  leaf biomass (KgC/m2)
      real cawood_pft(nx,ny,12,q) !monthly aboveground biomass (KgC/m2)
      real cfroot_pft(nx,ny,12,q) !monthly fine root
      
      
c     --------------------------------E N D--------------------------------
c     ------------------------- internal variables-------------------------

      integer i, j, k, kk, mes, nerro, p
      
      real wg0  (nx,ny,12)      !Moisture of the previous year
      real wsoilt(nx,ny,12)     !soil water (used to check wbm equilibrium)
      real gsoilt(nx,ny,12)     !soil ice   (used to check wbm equilibrium)
      real gsoil(nx,ny,12,q)    !Soil ice 
      real ssoil(nx,ny,12,q)    !Soil snow
      real snowm(nx,ny,12,q)    !Snowmelt

      real sini(q), sfim(q)
      real wini(q), wfim(q)
      real gini(q), gfim(q)

      real cleaf1_pft (q), clfim(q)
      real cawood1_pft(q), cafim(q)
      real cfroot1_pft(q), crfim(q)
      
      real rmes(q),phmes(q),smes(q),rcmes(q),hrmes(q)
      real nppmes(q),laimes(q), armes(q),clmes(q),csmes(q)
      real emes(q),epmes(q),rmlmes(q),rmfmes(q),rmsmes(q)
      real rmmes(q),rglmes(q),rgfmes(q),rgsmes(q),rgmes(q)
      real cleafmes(q),cawoodmes(q),cfrootmes(q)

      
      real t0,t1                !Soil temperature aux variables 
      real wsaux1,dwww,wmax     !auxiliar to equilibrium check
      real ae                   !Available energy     
      real pr,spre,ta,td,ipar
      
      real, parameter :: H = 1.0 !Soil layer(m) 
      real, parameter :: diffu = 4.e-7*(30.*86400.0) !Soil thermal diffusivity (m2/month)
      real, parameter :: tau = (H**2)/(2.0*diffu) !E-folding time (months)
      real, parameter :: auxs = -100.0 !Auxiliar for calculation of Snpp
      
            

      
!     Soil temperature
!     ================
      
!     For all grid points
!     -------------------
      do i=1,nx
         do j=1,ny
!     Initialize soil temperature
!     ---------------------------
            do k=1,12
               tsoil(i,j,k) = no_data
            enddo

!     Only for land grid points
!     -------------------------!     
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
!     ============
!     Water budget
!     ============
!     
!     For all grid points
!     -------------------
!     
      do i=1,nx
         do j=1,ny
            do p = 1,q   
               cleaf1_pft (p) =  cleaf_ini(i,j,p)
               cawood1_pft(p) = cawood_ini(i,j,p)
               cfroot1_pft(p) = cfroot_ini(i,j,p)
            enddo
c     Write to track program execution     
            if ((mod(j,ny).eq.0).and.(mod(i,10).eq.0))
     &          write(*,*) 'working:',i           
!     Initialize variables
            do k=1,12
               wg0(i,j,k) = no_data               
               do p=1,q
                  emaxm(i,j,k,p)      = no_data !Maximum evapotranspiration
                  photo_pft(i,j,k,p)  = no_data !Monthly photosynthesis (kgC/m2)
                  aresp_pft(i,j,k,p)  = no_data !Monthly autotrophic respiration (kgC/m2)
                  npp_pft(i,j,k,p)    = no_data !Monthly net primary productivity (average between PFTs) (kgC/m2)
                  lai_pft(i,j,k,p)    = no_data !Monthly leaf area index
                  clit_pft(i,j,k,p)   = no_data !Monthly litter carbon
                  csoil_pft(i,j,k,p)  = no_data !Monthly soil carbon
                  hresp_pft(i,j,k,p)  = no_data !Monthly heterotrophic respiration (kgC/m2)
                  rcm_pft(i,j,k,p)    = no_data
                  gsoil(i,j,k,p)      = no_data !Soil ice
                  ssoil(i,j,k,p)      = no_data !Soil snow
                  runom_pft(i,j,k,p)  = no_data !Runoff
                  evapm_pft(i,j,k,p)  = no_data !Actual evapotranspiration        
                  wsoil_pft(i,j,k,p)  = no_data !Soil moisture (mm)
C     NOVAS VARIAVEIS VEG POOLS
                  rml_pft(i,j,k,p)    = no_data
                  rmf_pft(i,j,k,p)    = no_data
                  rms_pft(i,j,k,p)    = no_data
                  rm_pft(i,j,k,p)     = no_data
                  rgl_pft(i,j,k,p)    = no_data
                  rgf_pft(i,j,k,p)    = no_data
                  rgs_pft(i,j,k,p)    = no_data
                  rg_pft(i,j,k,p)     = no_data
                  cleaf_pft(i,j,k,p)  = no_data !monthly npp alloc to leaf biomass (KgC/m2)
                  cawood_pft(i,j,k,p) = no_data !monthly npp alloc to aboveground biomass (KgC/m2)
                  cfroot_pft(i,j,k,p) = no_data !monthly npp alloc to fine root biomass (KgC/m2)
                  
               enddo              
            enddo
!     
!     Only for land grid points
!     -------------------------
            if (nint(lsmk(i,j)).ne.0) then
               do k=1,12
                  wg0(i,j,k) = -1.0 ! initialize previous year water in soil (w + g)
               enddo
!     Set some variables
!     ------------------
               do p = 1,q
                  wini(p)  = 0.01  !Soil moisture_initial condition (mm)
                  gini(p)  = 0.0   !Soil ice_initial condition (mm)
                  sini(p)  = 0.0   !Overland snow_initial condition (mm)
               enddo
!     START INTEGRATION
               n = 0
 10            continue
               n = n + 1
!     Pre-processing
!     --------------
               k = mod(n,12)
               if (k.eq.0) k = 12
               mes = k
               spre = p0(i,j,k)
               td = tsoil(i,j,k)
               ta = temp(i,j,k)
               pr = prec(i,j,k)
               ipar = par(i,j,k)
               ae = 2.895*ta+52.326 !Available energy (W/m2) - From NCEP-NCAR reanalysis data
!     
!     Monthly water budget
!     ====================
               
               call budget (mes,wini,gini,sini,td,ta,pr,spre,ae
     $              , ca, ipar, cleaf1_pft, cawood1_pft, cfroot1_pft
     $              ,wfim,gfim, sfim, clfim, cafim,crfim,smes,rmes,emes
     $              ,epmes,phmes,armes,nppmes,laimes,clmes,csmes,hrmes
     $              ,rcmes,rmlmes,rmfmes, rmsmes, rmmes, rglmes, rgfmes
     $              ,rgsmes,rgmes, cleafmes, cawoodmes,cfrootmes) 

               do p=1,q
                  emaxm    (i,j,k,p) = epmes(p)
                  gsoil    (i,j,k,p) = gfim(p)
                  ssoil    (i,j,k,p) = sfim(p)
                  wsoil_pft(i,j,k,p) = wfim(p)
                  snowm    (i,j,k,:) = smes(p)
                  runom_pft(i,j,k,p) = rmes(p)
                  evapm_pft(i,j,k,p) = emes(p)
                  rcm_pft  (i,j,k,p) = rcmes(p)
                  lai_pft  (i,j,k,p) = laimes(p)
                  photo_pft(i,j,k,p) = phmes(p)
                  aresp_pft(i,j,k,p) = armes(p)
                  npp_pft  (i,j,k,p) = nppmes(p)
                  clit_pft (i,j,k,p) = clmes(p)
                  csoil_pft(i,j,k,p) = csmes(p)
                  hresp_pft(i,j,k,p) = hrmes(p)
                  rml_pft  (i,j,k,p) = rmlmes(p)
                  rmf_pft  (i,j,k,p) = rmfmes(p)
                  rms_pft  (i,j,k,p) = rmsmes(p)
                  rm_pft   (i,j,k,p) = rmmes(p)
                  rgl_pft   (i,j,k,p) = rglmes(p)
                  rgf_pft   (i,j,k,p) = rgfmes(p)
                  rgs_pft   (i,j,k,p) = rgsmes(p)
                  rg_pft    (i,j,k,p) = rgmes(p)
                  cleaf_pft (i,j,k,p) = cleafmes(p)
                  cawood_pft(i,j,k,p) = cawoodmes(p)
                  cfroot_pft(i,j,k,p) = cfrootmes(p)
                  wini(p)         = wfim(p)
                  gini(p)         = gfim(p)
                  sini(p)         = sfim(p)
                  cleaf1_pft(p)  = clfim(p) 
                  cawood1_pft(p) = cafim(p)
                  cfroot1_pft(p) = crfim(p)
               enddo
!     Check if equilibrium is attained
!     --------------------------------
               if (k.eq.12) then
                  wsoilt(i,j,k) = 0.0
                  gsoilt(i,j,k) = 0.0
                  do p = 1,q
                     wsoilt(i,j,k) = wsoilt(i,j,k) + wsoil_pft(i,j,k,p)
                     gsoilt(i,j,k) = gsoilt(i,j,k) + gsoil(i,j,k,p)
                  enddo
                  wmax = 500.
                  nerro = 0
                  do kk=1,12
                     wsaux1 = wsoilt(i,j,kk) + gsoilt(i,j,kk)
                     
                     dwww = (wsaux1 - wg0(i,j,kk)) / wmax
                     if (abs(dwww).gt.0.001) nerro = nerro + 1
                  enddo
                  if (nerro.ne.0) then
                     do kk=1,12
                        wg0(i,j,kk) = wsoilt(i,j,kk) + gsoilt(i,j,kk)
                     enddo
                  else
                     goto 100
                  endif
               endif
               goto 10               
 100           continue
               
            endif               ! endif lsmk
c     finalize ny loop
         enddo                  ! j
c     finalize nx loop
      enddo                     ! k
      return
      end subroutine wbm
      
      
      
      subroutine budget (month,w1,g1,s1,ts,temp,prec,p0,ae,ca
     $     ,ipar,cl1_pft,ca1_pft,cf1_pft,w2,g2,s2,cl2_pft,ca2_pft
     $     ,cf2_pft,smavg,ruavg,evavg,epavg,phavg,aravg,nppavg,laiavg
     $     ,clavg,csavg,hravg,rcavg,rmlavg,rmfavg,rmsavg,rmavg,rglavg
     $     ,rgfavg,rgsavg,rgavg,cleafavg_pft,cawoodavg_pft
     $     ,cfrootavg_pft)


      integer, parameter :: npft = 3
!     INPUTS
      
      integer month             !Actual month (1-12)
      real w1(npft)             !Initial (previous month last day) soil moisture storage (mm)
      real g1(npft)             !Initial soil ice storage (mm)
      real s1(npft)             !Initial overland snow storage (mm)
      real cl1_pft(npft)        ! initial npp allocation to cleaf compartment
      real cf1_pft(npft)        !                           froot
      real ca1_pft(npft)        !                           cawood
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
      real w2(npft)             !Final (last day) soil moisture storage (mm)
      real g2(npft)             !Final soil ice storage (mm)
      real s2(npft)             !Final overland snow storage (mm)
      real cl2_pft(npft)        !FINAL npp allocation to veg_pool
      real cf2_pft(npft)
      real ca2_pft(npft)
      real smavg(npft)          !Snowmelt monthly average (mm/day)
      real ruavg(npft)          !Runoff monthly average (mm/day)
      real evavg(npft)          !Actual evapotranspiration monthly average (mm/day)
      real epavg(npft)          !Maximum evapotranspiration monthly average (mm/day)
      real phavg(npft)          !Monthly photosynthesis
      real aravg(npft)          !Monthly autotrophic respiration
      real nppavg(npft)         !Monthly NPP (average between PFTs)
      real laiavg(npft)         !Monthly leaf area Index
      real clavg(npft)          !Monthly carbon litter
      real csavg(npft)          !Monthly carbon soil
      real hravg(npft)          !Monthly heterotrophic respiration
      real rcavg(npft)          !Monthly canopy resistence
      real rmlavg(npft),rmfavg(npft) ! maintenance/growth respiration
      real rmsavg(npft),rmavg(npft)  
      real rglavg(npft),rgfavg(npft)
      real rgsavg(npft),rgavg(npft)
      real cleafavg_pft(npft),cawoodavg_pft(npft),cfrootavg_pft(npft) ! vegetal carbon pools
      
!     Internal Variables
      integer p
      
      real alfa_leaf(npft), alfa_awood(npft), alfa_froot(npft)
      real beta_leaf(npft), beta_awood(npft), beta_froot(npft)

      ! RELATED WITH GRIDCELL OCUPATION AND LIGHT LIMITATION
      REAL TOTAL_BIOMASS_PFT(NPFT), OCP_COEFFS(NPFT)
      REAL MAX_VALUE, TOTAL_BIOMASS
      INTEGER LT(NPFT), MAX_INDEX
      
!     WBM COMMUNICATION (water balance)
      real rh                   !Relative humidity
      real wmax                 !Soil moisture availability (mm)
      real tsnow                !Temperature threshold for snowfall (oC)
      real tice                 !Temperature threshold for soil freezing (oC)
      real psnow                !Snowfall (mm/day)
      real prain                !Rainfall (mm/day)
      real rimelt               !Runoff due to soil ice melting
      real smelt                !Snowmelt (mm/day)
      real w   (npft)           !Daily soil moisture storage (mm)
      real g   (npft)           !Daily soil ice storage (mm)
      real s   (npft)           !Daily overland snow storage (mm)
      real ds  (npft)
      real dw  (npft)
      real roff(npft)           !Total runoff
      real evap(npft)           !Actual evapotranspiration (mm/day)
      real emax(npft)           !Maximum evapotranspiration
      
      integer ndmonth(12)       !Number of months
      data ndmonth /31,28,31,30,31,30,31,31,30,31,30,31/ !Number of days for each month
      
c     Carbon Cycle
      
      real ph  (npft)           !Canopy gross photosynthesis (kgC/m2/yr)
      real ar  (npft)           !Autotrophic respiration (kgC/m2/yr)
      real nppa(npft)           !Net primary productivity / auxiliar
      real laia(npft)           !Leaf area index (m2 leaf/m2 area)
      real cl  (npft)           !Litter carbon (kgC/m2) ---- anual?
      real cs  (npft)           !Soil carbon (kgC/m2) ---- anual?
      real hr  (npft)           !Heterotrophic (microbial) respiration (kgC/m2/yr)
      real rc2 (npft)           !Canopy resistence (s/m)
      real f1  (npft)
      real f5  (npft)           !Photosynthesis (mol/m2/s)
      real f1b (npft)           !Photosynthesis (micromol/m2/s)
      real vpd (npft)
      real rm(npft),rml(npft),rmf(npft),rms(npft)
      real rg(npft),rgl(npft),rgf(npft),rgs(npft)
      
      real cl1(npft), cf1(npft), ca1(npft)
      real cl2(npft),cf2(npft), ca2(npft)

      
!     Initialize Canopy Resistence Parameters
!     ---------------------------------------
      do p = 1,npft
         rc2(p) = 0.0
         f1(p)  = 0.0
         f1b(p) = 0.0
         f5(p) = 0.0
      enddo
      
!     Parameters
!     ----------
      rh    = 0.685
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
      
      do p = 1,npft
         LT(P)      = 0
         w(p)       = w1(p)
         g(p)       = g1(p)
         s(p)       = s1(p)
         smavg(p)   = 0.0
         ruavg(p)   = 0.0
         evavg(p)   = 0.0
         epavg(p)   = 0.0
         rcavg(p)   = 0.0
         laiavg(p)  = 0.0
         phavg(p)   = 0.0
         aravg(p)   = 0.0
         nppavg(p)  = 0.0
         clavg(p)   = 0.0
         csavg(p)   = 0.0
         hravg(p)   = 0.0
         rmlavg(p)  = 0.0
         rmfavg(p)  = 0.0
         rmsavg(p)  = 0.0
         rmavg(p)   = 0.0
         rglavg(p)  = 0.0
         rgfavg(p)  = 0.0
         rgsavg(p)  = 0.0
         rgavg(p)   = 0.0
         cleafavg_pft(p)  = 0.0 !mean monthly leaf biomass for all the PFTs
         cawoodavg_pft(p) = 0.0 !mean monthly aboveground biomass for all the PFTs
         cfrootavg_pft(p) = 0.0 !mean monthly belowground biomass for all the PFTs
         cl2_pft(p) = 0.0
         ca2_pft(p) = 0.0
         cf2_pft(p) = 0.0
         alfa_leaf(p)  = 0.0
         alfa_awood(p) = 0.0
         alfa_froot(p) = 0.0
      enddo
      
!     Numerical integration
!     --------------------
      do i=1,ndmonth(month)
         do p=1,npft
            nppa(p)  = 0.0  !Auxiliar_nppa
            ph(p)    = 0.0  !Auxiliar_ph
            ar(p)    = 0.0  !Auxiliar_ar
            laia(p)  = 0.0  !Auxiliar_laia
            f5(p)    = 0.0  !Auxiliar_f5
            f1(p)    = 0.0  !Auxiliar_f1
            vpd(p)   = 0.0
            rc2(p)   = 0.0  !Auxiliar_rc2
            rm(p)    = 0.0
            rml(p)   = 0.0
            rmf(p)   = 0.0
            rms(p)   = 0.0
            rg(p)    = 0.0
            rgl(p)   = 0.0
            rgf(p)   = 0.0
            rgs(p)   = 0.0
            cl1(p) = cl2_pft(p)
            ca1(p) = ca2_pft(p)
            cf1(p) = cf2_pft(p)
            beta_leaf(p) = alfa_leaf(p)
            beta_awood(p) = alfa_awood(p)
            beta_froot(p) = alfa_froot(p)
            
            if ((i.eq.1).and.(month.eq.1)) then
               cl1(p) = cl1_pft(p)
               ca1(p) = ca1_pft(p)
               cf1(p) = cf1_pft(p)
               
               beta_leaf(p)=0.
               beta_awood(p)=0.
               beta_froot(p)=0.
               LT(P) = 1
            endif

            
!     Maximum evapotranspiration (emax)
!     =================================
            call evpot2 (p0,temp,rh,ae,emax(p))
            
!     Productivity (ph, aresp, vpd, rc2)
!     =================================
            call productivity1 (p,LT(P),temp, p0, w(p), wmax, ca, ipar
     $          ,cl1(p), ca1(p), cf1(p), beta_leaf(p), beta_awood(p)
     $          ,beta_froot(p),emax(p), ph(p),ar(p), nppa(p),laia(p),
     $          f5(p),f1(p), vpd(p), rm(p), rml(p),rmf(p), rms(p),
     $          rg(p),rgl(p),rgf(p),rgs(p), rc2(p))

            
            call critical_value(nppa(p))
c     if(nppa(p) .gt. 0.0) print*, nppa(p), 'nppa_prod'
!     Carbon allocation (carbon content on each compartment)
!     =====================================================
            call allocation (p, nppa(p), cl1(p), ca1(p), cf1(p), !input
     &          cl2(p), ca2(p), cf2(p)) !output

            call critical_value(cl2)
            call critical_value(ca2)
            call critical_value(cf2)
c     if(nppa(p) .gt. 0.0) then
c     print*, cl1(p), 'cl1_prod',p,i
c     print*, cl2(p), 'cl2_prod',p,i
c     print*, cf1(p), 'cf1_prod',p,i
c     print*, cf2(p), 'cf2_prod',p,i
c     print*, ca1(p), 'ca1_prod',p,i
c     print*, ca2(p), 'ca2_prod',p,i
c     endif
            
            alfa_leaf(p)  = abs(cl2(p) - cl1(p)) 
            alfa_awood(p) = abs(ca2(p) - ca1(p)) 
            alfa_froot(p) = abs(cf2(p) - cf1(p)) 

!     Snow budget
!     ===========     
            smelt = 2.63 + 2.55*temp + 0.0912*temp*prain !Snowmelt (mm/day)
            smelt = amax1(smelt,0.)
            smelt = amin1(smelt,s(p)+psnow)
            ds(p) = psnow - smelt
            s(p) = s(p) + ds(p)
            
!     Water budget
!     ============
            if (ts.le.tice) then !Frozen soil
               g(p) = g(p) + w(p) !Soil moisture freezes
               w(p) = 0.0
               roff(p) = smelt + prain !mm/day
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
               rimelt = 0.0
               if (w(p).gt.wmax) then
                  rimelt = w(p) - wmax !Runoff due to soil ice melting
                  w(p) = wmax
               endif
     
!     Canopy resistance (Based on Medlyn et al.(2011)
!     (rc2 ; s/m)
!     
!     call canopy_resistence (p,vpd(p),f1(p),rc2(p))
!     canopy_resistence now called inside productivity (L 297)
               
               call runoff (w(p),wmax,roff(p)) !Soil moisture runoff (roff, mm/day)
               call penman (p0,temp,rh,ae,rc2(p),evap(p)) !Actual evapotranspiration (evap, mm/day)
               dw(p) = prain + smelt - evap(p) - roff(p)
               w(p) = w(p) + dw(p)
               if (w(p).gt.wmax) then
                  roff(p) = roff(p) + (w(p) - wmax)
                  w(p) = wmax
               endif
               if (w(p).lt.0.) w(p) = 0.
               roff(p) = roff(p) + rimelt !Total runoff

!     Carbon cycle (Microbial respiration, litter and soil carbon)
!     ============================================================     
               call carbon2 (ts,f5,evap(p),laia(p), !Inputs
     &              cl(p),cs(p),hr(p))   !Outputs
            endif
            
            if(nppa(p) .gt. 0.0)print*, nppa(p)
         enddo                  ! end p loop
         
!     Day i completed ---> calculate gridcell occupation and light limitation factor
         
         TOTAL_BIOMASS = 0.0
         
         do P = 1,NPFT
            
            TOTAL_BIOMASS_PFT(P) = CL2(P) + CF2(P) + CA2(P) ! BIOMASSA TOTAL DO PFT
            TOTAL_BIOMASS = TOTAL_BIOMASS + TOTAL_BIOMASS_PFT(P)
c            if(total_biomass .gt. 0.0) print*, total_biomass,'biomass',
c     &          p

         enddo
 
         DO P = 1,NPFT
            IF(TOTAL_BIOMASS .GT. 0.0) THEN
c               print*, total_biomass_pft(p),'biomass', p
               OCP_COEFFS(P) = TOTAL_BIOMASS_PFT(P) / TOTAL_BIOMASS
               IF(OCP_COEFFS(P) .LT. 0.0) OCP_COEFFS(P) = 0.0
               CALL CRITICAL_VALUE(OCP_COEFFS(P))
            ELSE
               OCP_COEFFS(P) = 0.0
            ENDIF
c            if(ocp_coeffs(p) .gt. 0.0) PRINT*, OCP_COEFFS(P), 'OCP'

         enddo

!     Light  limitation factor
         p = 0
         MAX_VALUE = 0.0
         MAX_INDEX = 0
         
         DO  P = 1,NPFT
c            if(ocp_coeffs(p) .gt. 0.0) PRINT*, OCP_COEFFS(P), 'OCP',
c     &          'light'
            IF(OCP_COEFFS(P) .GT. MAX_VALUE) THEN
               MAX_INDEX = P
               MAX_VALUE = OCP_COEFFS(P)
            ENDIF
            LT(P) = 0
         enddo
         
         
         IF(MAX_INDEX .GT. 0) THEN
            LT(MAX_INDEX) = 1
         ELSE
            DO P = 1,NPFT
               LT(P) = 1
            enddo
         ENDIF
         
         p = 0
!     Accumulate daily budgets weighted by occupation coefficients
         DO P = 1,NPFT

            epavg(p) = epavg(p) + emax(p) * OCP_COEFFS(P) !mm/day
            if (ocp_coeffs(p) .gt. 0.0 )print*, ocp_coeffs(p),
     &          'ocp in weighthing', p, i
            smavg(p) = smavg(p) + smelt
            ruavg(p) = ruavg(p) + roff(p) * OCP_COEFFS(P) !mm/day
            evavg(p) = evavg(p) + evap(p) * OCP_COEFFS(P)!mm/day
            
            rcavg(p) = rcavg(p) + rc2(p) * OCP_COEFFS(P) !s/m/day
            phavg(p) = phavg(p) + (ph(p) * OCP_COEFFS(P))/365.0 !kgC/m2/day
            aravg(p) = aravg(p) + (ar(p) * OCP_COEFFS(P))/365.0 !kgC/m2/day
            nppavg(p) = nppavg(p) + (nppa(p)* OCP_COEFFS(P))/365.0 !kgC/m2/day
            if(nppavg(p) .gt. 0.0) print*, nppavg(p)
     $           
            laiavg(p) = laiavg(p) + laia(p)/365.0 !m2leaf/m2area/day
     $           
            clavg(p) = clavg(p) + (cl(p) * OCP_COEFFS(P))/365.0 !kgC/m2/day
            csavg(p) = csavg(p) + (cs(p) * OCP_COEFFS(P))/365.0 !kgC/m2/day
            hravg(p) = hravg(p) + (hr(p) * OCP_COEFFS(P))/365.0 !kgC/m2/day
            rmlavg(p) = rmlavg(p) + (rml(p) * OCP_COEFFS(P))/365.
            rmfavg(p) = rmfavg(p) + (rmf(p) * OCP_COEFFS(P))/365.
            rmsavg(p) = rmsavg(p) + (rms(p) * OCP_COEFFS(P))/365.
            rmavg(p) = rmavg(p) + (rm(p) * OCP_COEFFS(P))/365.
            rglavg(p) = rglavg(p) + (rgl(p) * OCP_COEFFS(P)) /365.
            rgfavg(p) = rgfavg(p) + (rgf(p) * OCP_COEFFS(P)) /365.
            rgsavg(p) = rgsavg(p) + (rgs(p) * OCP_COEFFS(P))/365.
            rgavg(p) = rgavg(p) + (rg(p) * OCP_COEFFS(P))/365.
            cleafavg_pft(p) = cleafavg_pft(p) +  cl2(p) * OCP_COEFFS(P) 
            cawoodavg_pft(p) = cawoodavg_pft(p) + ca2(p) * OCP_COEFFS(P)
            cfrootavg_pft(p) = cfrootavg_pft(p) +  cf2(p) *
     &          OCP_COEFFS(P)
         enddo
         
      enddo                     ! end ndmonth loop
      

!     Final calculations
!     ------------------
!     monthly values
      do p=1,NPFT
         ocp_coeffs(p) = 0.0
         w2(p) = w(p)
         g2(p) = g(p)
         s2(p) = s(p)
         cl2_pft(p) = cl2(p) 
         ca2_pft(p) = ca2(p)
         cf2_pft(p) = cf2(p)
         smavg(p) = smavg(p)/real(ndmonth(month))
         ruavg(p) = ruavg(p)/real(ndmonth(month))
         evavg(p) = evavg(p)/real(ndmonth(month))
         epavg(p) = epavg(p)/real(ndmonth(month))
         rcavg(p) = rcavg(p)/real(ndmonth(month))
         phavg(p) = phavg(p) * 12.0 !kgC/m2/yr
         aravg(p) = aravg(p)* 12.0 !kgC/m2/yr
         nppavg(p) = nppavg(p) * 12.0 !kgC/m2/yr
         laiavg(p) = laiavg(p) * 12.0
         clavg(p) = clavg(p) * 12.0   !kgC/m2
         csavg(p) = csavg(p) * 12.0   !kgC/m2
         hravg(p) = hravg(p) * 12.0   !kgC/m2/yr
         rmlavg(p) = rmlavg(p) * 12.0 
         rmfavg(p) = rmfavg(p) * 12.0
         rmsavg(p) = rmsavg(p) * 12.0
         rmavg(p) = rmavg(p) * 12.0 
         rglavg(p) = rglavg(p) * 12.0 
         rgfavg(p) = rgfavg(p) * 12.0 
         rgsavg(p) = rgsavg(p) * 12.0 
         rgavg(p) = rgavg(p) * 12.0 
         cleafavg_pft(p) = cleafavg_pft(p)  * 12.0
         cawoodavg_pft(p) = cawoodavg_pft(p) * 12.0
         cfrootavg_pft(p) = cfrootavg_pft(p) * 12.0
c         IF(NPPAVG(P) .GT. 0.0)print*, nppavg(p), 'nppavg', p
      enddo
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
