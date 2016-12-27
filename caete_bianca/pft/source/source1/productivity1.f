c234567
!
!==============================================================================================
! Productivity1: Photosynthesis, Plant Respiration and NPP
! Carbon2: Microbial Respiration
!
! Code written by David Lapola and Helena Alves do Prado
! Last update: 10/2016
!
! A "bug" has been found on Sep/2007: when the water soil is between 0.205 and 0.265, NPP
! unrealistically drops to a level below those when wsoil is lesser than 0.205 (because of f5)
!=============================================================================================
!
      subroutine productivity1 (pft, temp, p0, w, wmax, ca, ipar, tsoil
     $     ,cl1, ca1, cf1, beta_leaf, beta_awood, beta_froot, ph,
     $     ar, nppa, laia, f5, f1, vpd, rm, rml, rmf, rms, rg, rgl, rgf
     $     ,rgs)
!     implicit none
!     
!     Variables
!     =========
!     
!     Input
!     -----
!     
      integer pft
      real temp                 !Mean monthly temperature (oC)
      real p0                   !Mean surface pressure (hPa)
      real wa,w,wmax            !Soil moisture (dimensionless)
      real ca                   !Atmospheric CO2 concentration (Pa)
      real ipar                 !Incident photosynthetic active radiation (w/m2)
!
      real cl1, cf1, ca1
      real beta_leaf, beta_awood, beta_froot
      real rm, rml, rmf, rms, rg, rgl, rgf, rgs
!     Output
!     ------
!     
      real ph                   !Canopy gross photosynthesis (kgC/m2/yr)
      real laia                 !Autotrophic respiration (kgC/m2/yr)
      real ar                   !Leaf area index (m2 leaf/m2 area)
      real nppa                 !Net primary productivity (kgC/m2/yr)
      real f5
!     
!     Internal
!     --------
!     
      real vm                   !Rubisco maximum carboxylaton rate (molCO2/m2/s)
      real mgama                !Photo-respiration compensation point (Pa)
      real es                   !Saturation partial pressure (Pa?)
      real es2                  !Saturation partial pressure (kPa?)
      real rmax                 !Saturated mixing ratio (kg/kg)
      real r                    !Moisture deficit at leaf level (kg/kg)
      real ci                   !Internal leaf CO2 partial pressure (Pa)
      real a,b,c,a2,b2,c2       !Auxiliars
      real delta,delta2         !Auxiliars
      real rl                   !Leaf respiration (kgC/m2/yr)
      real rp                   !Non-leaf parts respiration (kgC/m2/yr)
      real sunlai,shadelai      !Sunlai/Shadelai
      real vpd                  !Vapor pressure deficit (kPa)
!     
!     Rubisco, light and transport limited photosynthesis rate
!     --------------------------------------------------------
!     
      real jc,jl,je,jp,jp1,jp2,j1,j2 !Auxiliars
!     
!     Functions
!     ---------
!     
      real f1                   !Leaf level gross photosynthesis (molCO2/m2/s)
      real f1a                  !auxiliar_f1
      real f1b                  !Leaf level gross photosynthesis (micromolCO2/m2/s)
      real f2                   !Michaelis-Menton CO2 constant (Pa)
      real f3                   !Michaelis-Menton O2 constant (Pa)
      real f4,f4sun,f4shade     !Scaling-up to canopy level (dimensionless)
!

      ! BIANCA ___________________________________________________
      real tleaf(3)             !leaf turnover time (yr)
      data tleaf /1.0,0.5,1.0/  !leaf turnover time for the 3 PFTs
      real sla                  !specific leaf area (m2/kg)
      real cl2                  !leaf compartment's carbon content (kgC/m2)
      real ca2                  !aboveground woody compartment's carbon content (kgC/m2)
      real cf2                  !fineroots compartment's carbon content (kgC/m2)
      real csa                  !sapwood compartment´s carbon content (5% of woody tissues) (kgC/m2)
      real ncl                  !leaf N:C ratio (gN/gC)
      real ncf                  !fine roots N:C ratio (gN/gC)
      real ncs                  !sapwood N:C ratio(gN/gC)
      real tsoil                !soil temperature (ºC)
      real csai 
      real pt                   !taxa potencial de fornecimento para transpiração (mm/dia)
      real csru                 !Specific root water uptake (0.5 mm/gC/dia; based in Pavlick et al (2013))
      real ad                   !atmospheric demand for transpiration (mm/dia;based in Gerten et al. 2004)
      real emax                 !potential evapotranspiration (mm/dia)
      real alfm                 !maximum Priestley-Taylor coefficient (based in Gerten et al. 2004; used to calculate ad)
      real gm                   !scaling conductance (mm/dia)(based in Gerten et al. 2004; used to calculate ad)
      real gc                   !Canopy conductance (mm/dia)(based in Gerten et al. 2004; used to calculate ad)

c     HELENA____________________________________________________
!     Parameters 
!     ----------
!     
      real p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,
     &     p19,p20,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31
      real p21(3)
!     
      a   = 0.83                !Photosynthesis co-limitation coefficient
      a2  = 0.93                !Photosynthesis co-limitation coefficient
      p3  = 21200.0             !Atmospheric oxygen concentration (Pa)
      p4  = 0.08                !Quantum efficiency (mol electrons/Ein)
      p5  = 0.15                !Light scattering rate
      p6  = 2.0                 !Parameter for jl
      p7  = 0.50                !Ratio of light limited photosynthesis to Rubisco carboxylation
      p8  = 5200.0              !Photo-respiration compensation point
      p9  = 0.57                !Photosynthesis co-limitation coefficient
      p10 = 0.10                !Q10 function
      p11 = 25.0                !Q10 function reference temperature (oC)
      p12 = 30.0                !Michaelis-Menten constant for CO2 (Pa)
      p13 = 2.10                !Michaelis-Menten constant for CO2
      p14 = 30000.0             !Michaelis-Menten constant for O2 (Pa)
      p15 = 1.20                !Michaelis-Menten constant for O2
      p19 = 0.90                !Maximum ratio of internal to external CO2
      p20 = 0.10                !Critical humidity deficit (kg/kg)
!     p21 = 0.00004                                                     !Maximum Rubisco carboxylation rate (molCO2/m2/s)
      p22 = 2.0                 !Rubisco carboxylation rate
      p23 = 0.30                !Rubisco carboxylation rate
      p24 = 36.0                !Rubisco carboxylation rate (oC)
      p25 = 0.000008            !Maximum gross photosynthesis rate (molCO2/m2/s)
      p26 = 0.50                !light extinction coefficient for IPAR/sun (0.5/sen90)
      p27 = 1.50                !light extinction coefficient for IPAR/shade (0.5/sen20)
      p28 = 0.500               !Soil moisture at field capacity
      p29 = 0.205               !Soil moisture at wilting point
      p30 = 0.015               !Ratio of respiration to Rubisco carboxylation rates
      p31 = 3.85                !Whole plant to leaf respiration ratio
!     
!     Initialize
!     ----------
!     
      vm       = 0.0
      rmax     = 0.0
      r        = 0.0
      es       = 0.0
      ci       = 0.0
      mgama    = 0.0
      jc       = 0.0
      jl       = 0.0
      je       = 0.0
      jp       = 0.0
      rl       = 0.0
      rp       = 0.0
      f4sun    = 0.0
      f4shade  = 0.0
      f1a      = 0.0
      f1b      = 0.0
      a        = 0.0
      b        = 0.0
      c        = 0.0
      a2       = 0.0
      b2       = 0.0
      c2       = 0.0
      delta    = 0.0
      delta2   = 0.0
      sunlai   = 0.0
      shadelai = 0.0
!     
!     ==============
!     Photosynthesis 
!     ==============
!     
!     Rubisco maximum carboxylaton rate (molCO2/m2/s)
!     -----------------------------------------------
!     
      data p21 /0.00003,0.000058,0.000067/ !Tropical/Temperate/Boreal


      
      vm = (p21(pft)*p22**(p10*(temp-p11)))/ !Free-range parameter --> 0.0358>vm>840 (micromol)
     &     (1.0+exp(p23*(temp-p24)))
!     
!     Photo-respiration compensation point (Pa)
!     -----------------------------------------
!     
      mgama = p3/(p8*(p9**(p10*(temp-p11))))
!     
!     Michaelis-Menton CO2 constant (Pa)
!     ----------------------------------
!     
      f2 = p12*(p13**(p10*(temp-p11)))
!     
!     Michaelis-Menton O2 constant (Pa)
!     ---------------------------------
!     
      f3 = p14*(p15**(p10*(temp-p11)))
!     
!     Saturation partial pressure of water vapour (Pa)                                       
!     ------------------------------------------------
!     
      call tetens (temp,es)     !-71.1<temp<+38.2
!     
      es2 = es*100.0
      vpd = (((100.0-68.5)/100.0)*(es2))/1000.0 !kPa
!     
!     Saturated mixing ratio (kg/kg)
!     ------------------------------
!     
      rmax = 0.622*(es/(p0-es))
!     
!     Moisture deficit at leaf level (kg/kg)
!     --------------------------------------
!     
      r = -0.315*rmax
!     
!     Internal leaf CO2 partial pressure (Pa)
!     ---------------------------------------
!     
      ci = p19*(1-(r/p20))*(ca-mgama)+mgama
!     
!     Rubisco carboxilation limited photosynthesis rate (molCO2/m2/s)
!     ---------------------------------------------------------------
!     
      jc = vm*((ci-mgama)/(ci+(f2*(1+(p3/f3)))))
!     
!     Light limited photosynthesis rate (molCO2/m2/s)
!     -----------------------------------------------
!     
      jl = p4*(1.0-p5)*ipar*((ci-mgama)/(ci+(p6*mgama)))
!     
!     Transport limited photosynthesis rate (molCO2/m2/s)
!     ---------------------------------------------------
!     
      je = p7*vm
!     
!     Jp (minimum between jc and jl)
!     ------------------------------
!     
      a = 0.83
      b = (-1)*(jc+jl)
      c = jc*jl
      delta = (b**2)-4.0*a*c
!     
      jp1=(-b-(sqrt(delta)))/(2.0*a)
      jp2=(-b+(sqrt(delta)))/(2.0*a)
      jp= amin1(jp1,jp2)
!     
!     Leaf level gross photosynthesis (minimum between jc, jl and je)
!     ---------------------------------------------------------------
!     
      a2 = 0.93
      b2 = (-1)*(jp+je)
      c2 = jp*je
      delta2 = (b2**2)-4.0*a2*c2
!     
      j1=(-b2-(sqrt(delta2)))/(2.0*a2)
      j2=(-b2+(sqrt(delta2)))/(2.0*a2)
      f1a = amin1(j1,j2)
c      PRINT*, F1A, 'f1a'
!     Soil water
!     ==========
!     
      wa = w/wmax
!     
!     Water stress response modifier (dimensionless)
!     ----------------------------------------------
!     
      if (wa.gt.0.5) f5 = 1.0   !Not too lower in e.g. Amazonian dry season
      if ((wa.ge.0.205).and.(wa.le.0.5))
     &     f5 = (wa-0.205)/(0.5-0.205)
      if (wa.lt.0.205) f5 = wa  !Below wilting point f5 accompains wa (then Sahara is well represented)
!     
!     ...f1

! temos que mudar aqui no futuro... a rc2 nao eh disponivel para a
!     chamada da productivity, uma vez que ela é calculada pela rotina
!     canopy_resistence... por enquanto fica a moda antiga mesmo 
c      csru = 0.5 
c      pt = csru*(cf1*1000.)*wa  !(based in Pavlick et al. 2013; *1000. converts kgC/m2 to gC/m2)
c      alfm = 1.391
c      gm = 3.26*86400           !(*86400 transform mm/s to mm/dia)
c      gc = (1/rc2)*86400000     !*86400000 transfor m/s to mm/dia)
c      d =(emax*alfm)/(1+gm/gc)  !(based in Gerten et al. 2004)
c      f5 = 1-(exp(-1*(pt/d)))
!     
!     Photosysthesis minimum and maximum temperature
!     ----------------------------------------------
!     
      if ((temp.ge.-10.0).and.(temp.le.50.0)) then
         f1 = f1a*f5            !f5:water stress factor
      else
         f1 = 0.0               !Temperature above/below photosynthesis windown
      endif
c      if (f1 .gt. 0.5) PRINT*, f1, 'f1'
!     Leaf area index (m2 leaf/m2 arEa) bianca
!     ---------------------------------
      sla=(0.030*1000.)*((365/(((tleaf(pft))*365)/12))**(-0.46))
c      PRINT*, sla, 'sla'
      laia  = 0.25*exp(2.5*(f1/p25)) !Adjusted after using observed ipar
c      if(cl1 .gt. 0) print*, cl1, 'cl1'
c      laia = (cl1*365.*sla)
c      if(laia .gt. 0) PRINT*, laia, 'laia'
!     SunLAI
!     ------
!     
      sunlai = (1.0-(exp(-p26*laia)))/p26
!     
!     ShadeLAI
!     --------
!     
      shadelai = laia - sunlai
!     
!     Scaling-up to canopy level (dimensionless)
!     ------------------------------------------
!     
      f4 = (1.0-(exp(-p26*laia)))/p26 !Sun 90 degrees in the whole canopy, to be used for respiration
!     
!     Sun/Shade approach to canopy scaling                                  !Based in de Pury & Farquhar (1997)
!     ------------------------------------
!     
      f4sun = (1.0-(exp(-p26*sunlai)))/p26 !sun 90 degrees
      f4shade = (1.0-(exp(-p27*shadelai)))/p27 !sun ~20 degrees
!     
!     Canopy gross photosynthesis (kgC/m2/yr)
!     =======================================
!     (0.012 converts molCO2 to kgC)
!     (31557600 converts seconds to year [with 365.25 days])
!     
      ph = 0.012*31557600.0*f1*f4sun*f4shade
c      PRINT*, PH, 'ph'
!     Plant respiration
!     =================
!     c Maintenance respiration (kgC/m2/yr) (based in Ryan 1991)
      
      
      csa= 0.05*ca1             !sapwood carbon content (kgC/m2). 5% of woody tissues (Pavlick, 2013)
      ncl = 0.034               !(gN/gC)
      ncf = 0.034               !(gN/gC)
      ncs = 0.003               !(gN/gC)
      rml = (ncl*cl1)*27*(exp(0.03*temp))
      rmf = (ncf*cf1)*27*(exp(0.03*tsoil))
      rms = (ncs*csa)*27*(exp(0.03*temp))
      
      rm = rml + rmf + rms
c      print*, rm, 'rm'
      
c     Growth respiration (KgC/m2/yr)(based in Ryan 1991; Sitch et al. 2003; Levis et al. 2004)         
      
      csai= 0.05*beta_awood
      rgl = (0.25*((beta_leaf)*365))
      rgf = (0.25*((beta_froot)*365))
      rgs = (0.25*(csai)*365)
      
      rg = rgl + rgf + rgs
c      print*, rg, 'rg'
      if (rg.lt.0) then
         rg = 0
      endif
!     c Autotrophic (plant) respiration -ar- (kgC/m2/yr)
!     Respiration minimum and maximum temperature
!     -------------------------------------------
!     
      if ((temp.ge.-10.0).and.(temp.le.50.0)) then
         ar = rm+rg
c         if (ar .gt. 0.)PRINT*, AR ,'ar'
      else
         ar = 0.0               !Temperature above/below respiration windown
      endif
!     
!     ============
!     Productivity
!     ============
!     
!     Net primary productivity(kgC/m2/yr)
!     ===================================
!     
      nppa = ph-ar
        
c       if (nppa .gt. 0.) print*, nppa, 'npp'
c      PRINT*, NPPA
c      if (nppa.lt.0.0) nppa = 0.0 !Maybe this is incorrect, but demands retuning of every biome limits
!     
      return
      end subroutine productivity1
!     
!     ======================
!     Canopy resistence(s/m)
!     ======================
!     
      subroutine canopy_resistence (m,vpd,f1,rc2)
!     
!     Variables
!     =========
!     
!     Inputs
!     ------
!     
      integer m
      real f1                   !Photosynthesis (molCO2/m2/s)
      real vpd                  !kPa
!     
!     Outputs
!     -------
!     
      real rc2                  !Canopy resistence (s/m)
!     
!     Internal
!     --------
!     
      real f1b                  !Photosynthesis (micromolCO2/m2/s)
      real gs2                  !Canopy conductance (m/s)
      real gs                   !Canopy conductance (molCO2/m2/s)
      real g0                   !Residual stomatance conductance
      real g1(3)                !3.00 (boreal) / 2.05 (temperate) / 3.02 (tropical)
!     real g1 
      real D                    !kPA
      real aa
      real rcmax 
!     
      data g1 /6.0,4.0,2.0/ 
      f1b = (f1*10e5)           !maior f1b = 13.38
      aa = (f1b/363)
      g0 = 0.01
!     g1 = 4.9
      D = sqrt(vpd)
!     rcmax = 200.0                                                    !A=2.0
      rcmax = 550.0             !A=0.5                             
!     rcmax = 1800.0                                                   !A=0.1
!     
      if (vpd.lt.0.25) gs = 1.5
      if (vpd.ge.0.25) then 
         gs = g0 + 1.6 * (1 + (g1(m)/D)) * (aa) !Based on Medlyn et al. 2011
!     gs = g0 + 1.6 * (1 + (g1/D)) * (aa)                               !Based on Medlyn et al. 2011
      endif
!     
!     if (gs.le.0.01) gs = 0.0                                        
!     if (gs.le.0.0) then                                             !-0.029<gs<1.29
!     print*,gs
!     endif
!     
      gs2 = gs/41
!     if (gs2.le.0.00025) gs2 = 0.0
!     if (gs2.gt.0.0318) then                                         !-0.00069<gs2<0.0317
!     print*,print
!     endif
!     
      rc2 = 1/gs2
!     if (rc2.ge.151.47) rc2 = rcmax                                    !A=2
      if (rc2.ge.545.43) rc2 = rcmax !A=0.5  
!     if (rc2.ge.1779.97) rc2 = rcmax                                   !A=0.1
!     
!     if (rc2.lt.19.57) then
!     print*,rc2
!     endif
!     
      return
      end subroutine
!     
!     =====================================
!     Microbial (heterotrophic) respiration
!     =====================================
      
      subroutine carbon2 (tsoil,f5,evap,laia, !Inputs
     &     cl,cs,hr)            !Outputs
!     
!     Variables
!     =========
!     
!     Inputs
!     ------
      
      real tsoil                !Mean monthly soil temperature (oC)
      real f5                   !Stress response to soil moisture (dimensionless)
      real evap                 !Actual evapotranspiration (mm/day)
      real laia
!     
!     Outputs 
!     -------
!     
      real cl                   !Litter carbon (kgC/m2)
      real cs                   !Soil carbon (kgC/m2)
      real hr                   !Heterotrophic (microbial) respiration (kgC/m2)
!     
!     Internal
!     --------
!     
      real lf                   !Litterfall (kgC/m2)
      real f6                   !Litter decayment function
      real f7                   !Soil carbon storage function
      real p10                  !Q10 function
      real p11                  !Q10 function reference temperature (oC)
      real p32                  !Q10 parameter of soil respiration sensibility to temperature
      real p33                  !Litterfall (kg/C/yr)
      real p34                  !Average fraction of litter carbon lost to atmosphere
      real p35                  !Carbon soil turnover (1/20yr)
      real p36                  !Specific rate heterotrophic respiration
!     
!     Initialize
!     ----------
!     
      lf  = 0.0
      f6  = 0.0
      f7  = 0.0
      cl  = 0.0
      cs  = 0.0
!     
      p10 = 0.10
      p11 = 25.0
      p32 = 2.00
      p33 = 0.10
      p34 = 0.30
      p35 = 0.05
      p36 = 0.25
!     
!     Litter decayment function                                             !Controlled by annual evapotranspiration
!     -------------------------
!     
      f6 = 1.16*10**(-1.4553+0.0014175*(evap*365.0))
!     
!     Soil carbon storage function                                          !Controlled by temperature
!     ----------------------------
!     
      f7 = p32**(p10*(tsoil-p11))
!     
!     Litterfall (kgC/m2)
!     ------------------
!     
      lf = p33*laia
!     
!     Litter carbon (kgC/m2)
!     ----------------------
!     
      cl = lf/f6
!     
!     Soil carbon(kgC/m2)
!     -------------------
!     
      cs = ((p34*cl)/(p35*f7))*f5
!     
!     Respiration minimum and maximum temperature
!     -------------------------------------------
!     
      if ((tsoil.ge.-10.0).and.(tsoil.le.50.0)) then
         hr = p36*(cl*(f6**2)+(cs*f5*evap*(f7**2))) !Litter and Soil respectively
      else
         hr = 0.0               !Temperature above/below respiration windown
      endif
!     
      return
      end
!     ===================================================
