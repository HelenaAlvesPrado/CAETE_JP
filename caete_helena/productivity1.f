c23456789
!

      subroutine productivity1 (temp,p0,w,wmax,ca,ipar,l,                                      !input
     &                          ph,ar,nppa,laia,f5,vm,mgama,
     &                          f2,f3,rmax,r,ci,f10,jc)                                      !output
!
!==============================================================================================
!
! Productivity1: Photosynthesis, Plant Respiration and NPP
! Carbon2: Microbial Respiration (to be improved)
!
! Code written by David Lapola and associates!
! Last update: 06/2016
!
! A "bug" has been found on Sep/2007: when the water soil is between 0.205 and 0.265, NPP
! unrealistically drops to a level below those when wsoil is lesser than 0.205 (because of f5)
!
!==============================================================================================
!
! Input/Output Variables
! ----------------------
!
! Input
! -----
!
!     temp           : Mean monthly temperature                                                !oC
!     p0             : Mean surface pressure                                                   !hPa
!     wa,w,wmax      : Soil moisture                                                           !Dimensionless
!     ca             : Atmospheric CO2 pressure                                                !Pa
!     ipar           : Incident photosynthetic active radiation                                !w/m^2                   !
!
! Output
! ------
!
!     ph             : Canopy gross photosynthesis                                             !kgC/m2/yr
!     ar             : Autotrophic respiration                                                 !kgC/m2/yr
!     nppa           : Net primary productivity                                                !kgC/m2/yr
!     laia           : Leaf area index                                                         !m2 leaf/m2 area
!
!==============================================================================================
!
! Input/Output Variables
! ----------------------
!
      integer l
      real temp
      real p0
      real wa,w,wmax
      real ca
      real ipar
      
      real ph
      real laia
      real ar
      real nppa
!
! Internal Variables
! ------------------
!
      real vm,mgama,es,rmax,r,ci,a,b,c,a2,b2,c2,
     &     delta,delta2,f1,f1a,f2,f3,f4,f5,k16,k17,rl,rp
! CO2 sensitivity function
      real f10
!
! Rubisco, light and transport limited photosynthesis rate
! --------------------------------------------------------
!
      real jc,jl,je,
     &     jp,jp1,jp2,j1,j2	                                                               !auxiliars
!
! Parameters
! ----------
      real p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p22,
     &     p23,p24,p25,p26,p27,p28,p29,p30,p31
      real p21(3)
!
      p1  = 0.93         !Photosynthesis co-limitation coefficient
      p2  = 0.83         !Photosynthesis co-limitation coefficient
      p3  = 21200        !Atmospheric oxygen concentration                                     Pa
      p4  = 0.08         !Quantum efficiency                                                   mol electrons/Ein
      p5  = 0.15         !Light scattering rate
      p6  = 2.0          !Parameter for jl
      p7  = 0.5          !Ratio of light limited photosynthesis to Rubisco carboxylation
      p8  = 5200         !Photo-respiration compensation point
      p9  = 0.57         !Photosynthesis co-limitation coefficient
      p10 = 0.1          !Q10 function
      p11 = 25           !Q10 function reference temperature                                 !oC
      p12 = 30           !Michaelis-Menten constant for CO2                                  !Pa
      p13 = 2.1          !Michaelis-Menten constant for CO2
      p14 = 30000        !Michaelis-Menten constant for O2                                   !Pa
      p15 = 1.2          !Michaelis-Menten constant for O2
      p16 = 1.23         !CO2 sensitivity function
      p17 = 24.11        !CO2 sensitivity function
      p18 = 8.93         !CO2 sensitivity function
      p19 = 0.9          !Maximum ratio of internal to external CO2
      p20 = 0.1          !Critical humidity deficit                                          !kg/kg
!      p21 = 0.00004      !: Maximum Rubisco carboxylation rate                                 !molCO2/m2/s
      p22 = 2            !Rubisco carboxylation rate
      p23 = 0.3          !Rubisco carboxylation rate
      p24 = 36           !Rubisco carboxylation rate                                         !oC
      p25 = 0.000008     !Maximum gross photosynthesis rate                                  !molCO2/m2/s
      p26 = 0.5          !light extinction coefficient for IPAR/sun
      p27 = 1.5          !light extinction coefficient for IPAR/shade
      p28 = 0.5          !Soil moisture at field capacity
      p29 = 0.205        !Soil moisture at wilting point
      p30 = 0.015        !Ratio of respiration to Rubisco carboxylation rates
      p31 = 3.85         !Whole plant to leaf respiration ratio
!
!==============================================================================================
!
! Photosynthesis
! ==============
!
! Rubisco maximum carboxylaton rate(molCO2/m2/s)
! ----------------------------------------------
!
!      data p21 /0.00003,0.000058,0.000067/
      data p21 /0.00004,0.00005,0.00006/

      vm = (p21(l)*p22**(p10*(temp-p11)))/                                                                              !ok
     &     (1.0+exp(p23*(temp-p24)))

!      print*,vm
!      stop
!
!     Values based on Try,considering 25oC, 345 until 400 ppmCO2 and excluding crazy numbers
!
!     Tropical : 30 micromol/m2/s
!     Temperate : 58 micromol/m2/s
!     Boreal : 67 micromol/m2/s
!
! Photo-respiration compensation point(Pa)
! ----------------------------------------
!
      mgama = 21200.0/(5200.0*(0.57**(0.1*(temp-25.0))))

!      print*,temp,mgama
!      stop
!
! Michaelis-Menton CO2 constant(Pa)
! ---------------------------------
!
      f2 = 30*(2.1**(0.1*(temp-25.0)))
!
!      print*,f2
!      stop

! Michaelis-Menton O2 constant(Pa)
! --------------------------------
!
      f3 = 30000.0*(1.2**(0.1*(temp-25.0)))
      
!      print*,f3
!      stop
!
! Saturation partial pressure of water vapour at temperature 'temp'(Pa)
! ---------------------------------------------------------------------
!
      call tetens (temp,es)
!
! Saturated mixing ratio(kg/kg)
! -----------------------------
!
      rmax = 0.622*(es/(p0-es))                                                                !pO=surface pressure(Pa)

!      print*,rmax
!      stop
!
! Moisture deficit at leaf level(kg/kg)
! -------------------------------------
!
!     relative humidity = 0.685
      r = (-0.315)*rmax                                                                        !rmax=saturated mixing ratio

!      print*,r
!      stop
!
! Internal leaf CO2 partial pressure(Pa)
! --------------------------------------
!      k16 = 0.875	!Generic value for forest types
!      k17 = 0.08	!Generic value for forest types
!      k16 = 0.9	!Value for C3 grass types
!      k17 = 0.1	!Value for C3 grass types
!
      ci = 0.9*(1-(r/0.1))*(ca-mgama)+mgama
!
! Reduction of NPP sensitivity to CO2 based in simulations for pre-industrial period
! ----------------------------------------------------------------------------------
!
      f10 = 1.23/(1+exp(-1*(ca-24.11)/8.93))
!
! Rubisco carboxilation limited photosynthesis rate(molCO2/m2/s)
! ==============================================================
!
      jc = vm*((ci-mgama)/(ci+(f2*f10*(1+(21200.0/f3)))))

!      print*,'temp',temp,'mgama',mgama,'f2',f2,'f3',f3,
!     &       'rmax',rmax,'r',r,'ci',ci,'f10',f10,'jc',jc
!      stop

!
! Light limited photosynthesis rate(molCO2/m2/s)
! ==============================================
!
      jl = 0.08*(1.0-0.15)*ipar*((ci-mgama)/(ci+(2.0*mgama)))
!
! Transport limited photosynthesis rate(molCO2/m2/s)
! ==================================================
!
      je = 0.5*vm                                                                              !p7 = ratio of light limited photosynthesis to Rubisco carboxylation
!
! Jp(minimum between jc and jl)
! -----------------------------
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
! Leaf level gross photosynthesis (minimum between jc, jl and je)                              !f1
! ---------------------------------------------------------------
!
      a2 = 0.93
      b2 = (-1)*(jp+je)
      c2 = jp*je
      delta2 = (b2**2)-4.0*a2*c2
!
      j1=(-b2-(sqrt(delta2)))/(2.0*a2)
      j2=(-b2+(sqrt(delta2)))/(2.0*a2)
         f1a = amin1(j1,j2)
!
! Soil water
! ==========
!
      wa = w/wmax
!
! Water stress response modifier(dimensionless)
! ---------------------------------------------
!
      if (wa.gt.0.5) f5 = 1.0	                                                               !Not too lower in e.g. Amazonian dry season
      if ((wa.ge.0.205).and.(wa.le.0.5))
     &   f5 = (wa-0.205)/(0.5-0.205)
      if (wa.lt.0.205) f5 = wa		                                                       !Below wilting point f5 accompains wa (then Sahara is well represented)
!
! ...f1
!
! Photosysthesis minimum and maximum temperature !f1 = leaf level gross photosynthesis
! ----------------------------------------------
!
      if ((temp.ge.-10.0).and.(temp.le.50.0)) then                                             !!!!!!!!!!!!!!!! Nao interfere !!!!!!!!!!!!!!!!!!
         f1 = f1a*f5
      else
         f1 = 0.0                                                                              !!!!!!!! interfere !!!!!! Mudando para 0.1, por ex, muda btt. Mas fica igual p todos os pfts         !Temperature above/below photosynthesis windown
      endif                                                                                    !!!! Mantendo 0.0, o padrao permanece igual para todos, mas os valores mudam !!!!!!!!
!
! Leaf area index(m2 leaf/m2 area)
! ================================
!
!     laia  = 0.2*exp(2.5*(f1/0.000008))
      laia  = 0.25*exp(2.5*(f1/0.000008))                                                           !adjusted after using observed ipar
!
! SunLAI
! ------
!
      sunlai = (1.0-(exp(-0.5*laia)))/0.5
!
! ShadeLAI
! --------
!
      shadelai = laia - sunlai
!
! Scaling-up to canopy level(dimensionless)
! =========================================
!
      f4 = (1.0-(exp(-0.5*laia)))/0.5	                                               !sun 90° in the whole canopy, to be used for respiration
!
! Sun/Shade approach to canopy scaling                                                         !Based in de Pury & Farquhar 1997
! ------------------------------------
!
      f4sun = (1.0-(exp(-0.5*sunlai)))/0.5	                                               !sun 90°
      f4shade = (1.0-(exp(-1.5*shadelai)))/1.5                                              !sun ~20°
!
! ======================================
! Canopy gross photosynthesis(kgC/m2/yr)
! ======================================
! (0.012 converts molCO2 to kgC)
! [31557600 converts seconds to year (with 365.25 days)]
!
      ph = 0.012*31557600.0*f1*f4sun*f4shade
!
! =================
! Plant respiration
! =================
!
! Maintenance respiration(kgC/m2/yr)
! ----------------------------------
!
! Leaf respiration
! ----------------
!
      rl = 0.012*31557600.0*0.015*vm*f4*f5
!
! Non-leaf parts respiration
! --------------------------
!
      rp = 3.85*rl

! Respiration minimum and maximum temperature
! -------------------------------------------
!
      if ((temp.ge.-10.0).and.(temp.le.50.0)) then                                             !!!!!!! interfere um pouco !!!!!!!!!
         ar = rl+rp
      else
         ar = 0.0	                                                                       !!!!!!!!!!!!!!!!! Nao interefere   !!!!!!!!!!!!!!!!!                   !Temperature above/below respiration windown
      endif                                                                                    ! Essa condi‡ao criada nao interferiu no padrao de npp !
!
! ============
! Productivity
! ============
!
! Net primary productivity(kgC/m2/yr)
! ===================================
!
      nppa = ph-ar
      if (nppa.lt.0.0) nppa = 0.0                                                              !maybe this is incorrect, but demands retuning of every biome limits

      return
      end subroutine productivity1
!
! ===========================
! Canopy resistence(mol/m2/s)
! ===========================
! [NPP*2.64 converts kgC/m2/yr to micromolCO2/m2/s]
!
      subroutine canopy_resistence (nppb,nppa,m,                                               !input
     &                              cres2)                                                     !output
!
      real nppb,nppa,cres2
      integer m
      real g1(3)
      data g1 /3.5,4.4,2.7/
!
!     es = 0.6108*exp((17.27*temp)/(temp+237.3))                                               !temp = 25oC
!     ea = (rh/100)*es                                                                         !rh = 68.5%
!
      vpd = 0.997850189                                                                        !kPa
      g0 = 0.5
!
      nppb = amax1(nppa,0.05)
      gs = g0+(1.6)*(1+(g1(m)/sqrt(vpd)))*((nppb*2.64)/350)                                    !mol/m^2/s
!
      gs2 = gs*1000                                                                            !mmol/m^2/s
      gs3 = gs2/41                                                                             !mm/s
      gs4 = gs3/1000                                                                           !m/s
      cres2 = 1/gs4
!
!     gs        : Optimal stomatal conductance to water vapour                                 !mol/m2/s
!     g0        : Residual conductance to water vapour (minimal stomatal conductance)          !mol/m2/s
!     g1        : Stomatal sensibility to assimilation                                         !kPa^0.5
!     nppb      : NPP                                                                          !micromol/m2/s
!     d         : Vapour Pressure Deficit                                                      !kPa
!     ca        : CO2 concentration                                                            !micromo/mol(ppm)
!
!
!                                                       PVM-2           Medlyn et al.(2001)
!
!                                                        g0                     g1                                      !unidade g0
!
!     Tropical Evergreen Tree (angiosperm)               0.5                   3.5 (n=6)
!     Temperate Deciduous Tree (angiosperm)              0.5                   4,4 (n=16)
!     Boreal Deciduous Tree (angiosperm)                 0.5                   2.7 (n=1)
!
      return
      end
!
! =====================
! Microbial respiration
! =====================
!
      subroutine carbon2 (tsoil,f5,evap,laia,                                                  !input
     &                    cl,cs,hr)                                                            !output
!
! Input/Output Variables
! ----------------------
!
      real tsoil	                                                                       !mean monthly soil temperature (oC)
      real f5		                                                                       !stress response to soil moisture (dimensionless)
      real evap		                                                                       !actual evapotranspiration (mm/day)
      real laia		                                                                       !Leaf area index (m2 leaf/m2 area)
!
      real cl		                                                                       !Litter carbon (kgC/m2)
      real cs		                                                                       !Soil carbon (kgC/m2)
      real hr		                                                                       !Heterotrophic (microbial) respiration (kgC/m2/yr)
!
! Internal Variables
! ------------------
!
      real lf,f6,f7
!
! =========================
! Litter decayment function                                                                    !controlled by annual evapotranspiration
! =========================
!
      f6 = 1.16*10**(-1.4553+0.0014175*(evap*365.0))
!
! Soil carbon storage function                                                                 !controlled by temperature
! ----------------------------
!
      f7 = 2.0**(0.1*(tsoil-25.0))
!
! Litterfall(kgC/m2)
! ------------------
!
      lf = 0.1*laia
!
! Litter carbon(kgC/m2)
! ---------------------
!
      cl = lf/f6
!
! Soil carbon(kgC/m2)
! -------------------
!
      cs = ((0.3*cl)/(0.05*f7))*f5
!
! Respiration minimum and maximum temperature
! -------------------------------------------
!
      if ((tsoil.ge.-10.0).and.(tsoil.le.50.0)) then
         hr = 0.25*(cl*(f6**2)+(cs*f5*evap*(f7**2)))                                           !Litter and Soil respectively
!
      else
         hr = 0.0	                                                                       !Temperature above/below respiration windown
      endif
!
      return
      end
!
!==============================================================================================
