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
      subroutine productivity1 (temp,p0,w,wmax,ca,ipar,l,
     &                          ph,ar,nppa,laia,f5,f1)                                             
!
! Variables
! =========
!
! Input
! -----
!
      integer l
      real temp		                                                
      real p0		                                                
      real wa,w,wmax	                                                
      real ca		                                                
      real ipar		                                                
!
! Output
! ------
!
      real ph		                                               
      real laia		                                          
      real ar		                                             
      real nppa	
      real f5	                                             
!
! Internal
! --------
!
      real vm
      real mgama
      real es
      real rmax
      real r
      real ci
      real a,b,c,a2,b2,c2
      real delta,delta2
      real rl
      real rp
      real sunlai,shadelai
!
! Rubisco, light and transport limited photosynthesis rate
! --------------------------------------------------------
!
      real jc,jl,je,
     &     jp,jp1,jp2,j1,j2	
!
! Functions
! ---------
!
      real f1,f1a
      real f2 
      real f3
      real f4,f4sun,f4shade
!      real f10
!
! Parameters
! ----------
!
      real p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,
     &     p19,p20,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31    
      real p21(3)
      real f10,p16,p17,p18
!
      a   = 0.83         !Photosynthesis co-limitation coefficient
      a2  = 0.93         !Photosynthesis co-limitation coefficient
      p3  = 21200.0      !Atmospheric oxygen concentration (Pa)
      p4  = 0.08         !Quantum efficiency (mol electrons/Ein)
      p5  = 0.15         !Light scattering rate
      p6  = 2.0          !Parameter for jl
      p7  = 0.5          !Ratio of light limited photosynthesis to Rubisco carboxylation
      p8  = 5200.0       !Photo-respiration compensation point
      p9  = 0.57         !Photosynthesis co-limitation coefficient
      p10 = 0.1          !Q10 function
      p11 = 25.0         !Q10 function reference temperature (oC)
      p12 = 30           !Michaelis-Menten constant for CO2 (Pa)
      p13 = 2.1          !Michaelis-Menten constant for CO2 
      p14 = 30000.0      !Michaelis-Menten constant for O2 (Pa)
      p15 = 1.2          !Michaelis-Menten constant for O2
      p16 = 1.23         !CO2 sensitivity function
      p17 = 24.11        !CO2 sensitivity function
      p18 = 8.93         !CO2 sensitivity function
      p19 = 0.9          !Maximum ratio of internal to external CO2
      p20 = 0.1          !Critical humidity deficit (kg/kg)
!      p21 = 0.00004       !Maximum Rubisco carboxylation rate (molCO2/m2/s)
      p22 = 2.0          !Rubisco carboxylation rate
      p23 = 0.3          !Rubisco carboxylation rate
      p24 = 36.0         !Rubisco carboxylation rate (oC)
      p25 = 0.000008     !Maximum gross photosynthesis rate (molCO2/m2/s)
      p26 = 0.5          !light extinction coefficient for IPAR/sun (0.5/sen90)                 
      p27 = 1.50         !light extinction coefficient for IPAR/shade (0.5/sen20)               
      p28 = 0.500          !Soil moisture at field capacity
      p29 = 0.205        !Soil moisture at wilting point
      p30 = 0.015        !Ratio of respiration to Rubisco carboxylation rates
      p31 = 3.85         !Whole plant to leaf respiration ratio
!
! Initialize
! ----------
!
      vm       = 0.
      rmax     = 0.
      r        = 0.
      es       = 0.
      ci       = 0.
      mgama    = 0.
      jc       = 0.
      jl       = 0.
      je       = 0.
      jp       = 0.
      rl       = 0.
      rp       = 0.
      f4sun    = 0.
      f4shade  = 0.
      f1a      = 0.
      a        = 0.
      b        = 0.
      c        = 0.
      a2       = 0.
      b2       = 0.
      c2       = 0.
      delta    = 0.
      delta2   = 0.                                                                      
      sunlai   = 0.
      shadelai = 0.       
!
! ==============
! Photosynthesis 
! ==============
!
! Rubisco maximum carboxylaton rate (molCO2/m2/s)
! -----------------------------------------------
!
      data p21 /0.00003,0.000058,0.000067/      !Tropica/Temperate/Boreal 
      vm = (p21(l)*p22**(p10*(temp-p11)))/                                                     
     &     (1.0+exp(p23*(temp-p24)))
!
! Photo-respiration compensation point (Pa)
! -----------------------------------------
!
      mgama = p3/(p8*(p9**(p10*(temp-p11))))
!
! Michaelis-Menton CO2 constant (Pa)
! ----------------------------------
!
      f2 = p12*(p13**(p10*(temp-p11)))
!
! Michaelis-Menton O2 constant (Pa)
! ---------------------------------
!
      f3 = p14*(p15**(p10*(temp-p11)))
!
! Saturation partial pressure of water vapour at temperature 'temp' (Pa)
! ----------------------------------------------------------------------
!
      call tetens (temp,es)
!
! Saturated mixing ratio (kg/kg)
! ------------------------------
!
      rmax = 0.622*(es/(p0-es)) 
!
! Moisture deficit at leaf level (kg/kg)
! --------------------------------------
!
      r = -0.315*rmax 
!
! Internal leaf CO2 partial pressure (Pa)
! ---------------------------------------
!
      ci = p19*(1-(r/p20))*(ca-mgama)+mgama
!
! Reduction of NPP sensitivity to CO2 based in simulations for pre-industrial period
! ----------------------------------------------------------------------------------
!
c      f10 = p16/(1+exp(-1*(ca-p17)/p18))
!
! Rubisco carboxilation limited photosynthesis rate (molCO2/m2/s)
! ---------------------------------------------------------------
!
c      jc = vm*((ci-mgama)/(ci+(f2*f10*(1+(p3/f3)))))
       jc = vm*((ci-mgama)/(ci+(f2*(1+(p3/f3)))))
!
! Light limited photosynthesis rate (molCO2/m2/s)
! -----------------------------------------------
!
      jl = p4*(1.0-p5)*ipar*((ci-mgama)/(ci+(p6*mgama)))
!
! Transport limited photosynthesis rate (molCO2/m2/s)
! ---------------------------------------------------
!
      je = p7*vm  
!
! Jp (minimum between jc and jl)
! ------------------------------
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
! Leaf level gross photosynthesis (minimum between jc, jl and je)                              
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
! Water stress response modifier (dimensionless)
! ----------------------------------------------
!
      if (wa.gt.0.5) f5 = 1.0 	                                                            !Not too lower in e.g. Amazonian dry season
      if ((wa.ge.0.205).and.(wa.le.0.5))
     &   f5 = (wa-0.205)/(0.5-0.205)
      if (wa.lt.0.205) f5 = wa		                                                      !Below wilting point f5 accompains wa (then Sahara is well represented)
!
! ...f1
!
! Photosysthesis minimum and maximum temperature                                               
! ----------------------------------------------
!
      if ((temp.ge.-10.0).and.(temp.le.50.0)) then
         f1 = f1a*f5
      else
         f1 = 0.0	                                                                              !Temperature above/below photosynthesis windown 
      endif

!      if (f1.lt.0.) then
!      print*,f1,'cu'
!      endif
!
! Leaf area index (m2 leaf/m2 area)
! ---------------------------------
!
      laia  = 0.25*exp(2.5*(f1/p25))  
!
! SunLAI
! ------
!
      sunlai = (1.0-(exp(-p26*laia)))/p26
!
! ShadeLAI
! --------
!
      shadelai = laia - sunlai
!
! Scaling-up to canopy level (dimensionless)
! ------------------------------------------
!
      f4 = (1.0-(exp(-p26*laia)))/p26	                                                      
!
! Sun/Shade approach to canopy scaling                                                       
! ------------------------------------
!
      f4sun = (1.0-(exp(-p26*sunlai)))/p26	                                            
      f4shade = (1.0-(exp(-p27*shadelai)))/p27	                                         
!
! Canopy gross photosynthesis (kgC/m2/yr)
! =======================================
! (0.012 converts molCO2 to kgC)
! (31557600 converts seconds to year [with 365.25 days])
!
      ph = 0.012*31557600.0*f1*f4sun*f4shade
!
! Plant respiration
! =================
!
! Maintenance respiration (kgC/m2/yr)
! ----------------------------------
!
      rl = 0.012*31557600.0*p30*vm*f4*f5  
!
      rp = p31*rl 
!
! Respiration minimum and maximum temperature
! -------------------------------------------
!
      if ((temp.ge.-10.0).and.(temp.le.50.0)) then
         ar = rl+rp
      else
         ar = 0.0	                                                                              !Temperature above/below respiration windown 
      endif
!
! ============
! Productivity
! ============
!
! Net primary productivity(kgC/m2/yr)
! ===================================
!
      nppa = ph-ar
      if (nppa.lt.0.0) nppa = 0.0                                                               !maybe this is incorrect, but demands retuning of every biome limits
!
      return
      end subroutine
!      
! ======================
! Canopy resistence(s/m)
! ======================
!
!      subroutine canopy_resistence (f1,rc2)

! Variables
! =========
!
!      real f1
!      real gs
!      real gs2
!      real rc2 
!
!      rc2 = (ca/(0.9*0.685*(nppb*2.64e-6)*(p0*100)))
!      gs = g0 + 1.6*(1+(g1/sqrt(vpd)))*((nppb*2.64e-6)/ca)
!
!      gs = -0.03 + 1.6 * (1+(4.9/sqrt(2.2))) * ((f1*1000000)/350)  !molCO2/m2/s
!      gs2 = gs/41                                                  !m/s
!      rc2 = 1/gs2                                                  !s/m
!      print*,rc2
!      return
!      end
!      
! =====================================
! Microbial (heterotrophic) respiration
! =====================================

      subroutine carbon2 (tsoil,f5,evap,laia,                                                   !input
     &                    cl,cs,hr)		                                                      !output
!
! Variables
! =========
!
! Inputs
! ------

      real tsoil	                                                                             !Mean monthly soil temperature (oC)
      real f5		                                                                       !Stress response to soil moisture (dimensionless)
      real evap		                                                                       !Actual evapotranspiration (mm/day)
      real laia
!
! Outputs 
! -------
!
      real cl		                                                                       !Litter carbon (kgC/m2)
      real cs		                                                                       !Soil carbon (kgC/m2)
      real hr		                                                                       !Heterotrophic (microbial) respiration (kgC/m2)
!
! Internal
! --------
!
      real lf                                                                                  
      real f6                                                                                  
      real f7                                                                                   
      real p10
      real p11
      real p32
      real p33
      real p34
      real p35
      real p36
!
! Initialize
! ----------
!
      lf  = 0.                                                                                 !Litterfall (kgC/m2)
      f6  = 0.                                                                                 !Litter decayment function 
      f7  = 0.                                                                                 !Soil carbon storage function
      cl  = 0.                                                                                   
      cs  = 0.
!
      p10 = 0.1                                                                                !Q10 function
      p11 = 25.0                                                                               !Q10 function reference temperature (oC)
      p32 = 2.0                                                                                !Q10 parameter of soil respiration sensibility to temperature  
      p33 = 0.1                                                                                !Litterfall (kg/C/yr)  
      p34 = 0.3                                                                                !Average fraction of litter carbon lost to atmosphere  
      p35 = 0.05                                                                               !Carbon soil turnover (1/20yr)
      p36 = 0.25    
!
! Litter decayment function (controlled by annual evapotranspiration)                                                                    !Controlled by annual evapotranspiration
! -------------------------------------------------------------------
!
      f6 = 1.16*10**(-1.4553+0.0014175*(evap*365.0))                                           
!
! Soil carbon storage function                                                                 !Controlled by temperature
! ----------------------------
!
      f7 = p32**(p10*(tsoil-p11))
!
! Litterfall (kgC/m2)
! ------------------
!
      lf = p33*laia
!
! Litter carbon (kgC/m2)
! ---------------------
!
      cl = lf/f6
!
! Soil carbon(kgC/m2)
! -------------------
!
      cs = ((p34*cl)/(p35*f7))*f5
!
! Respiration minimum and maximum temperature
! -------------------------------------------
!
      if ((tsoil.ge.-10.0).and.(tsoil.le.50.0)) then
         hr = p36*(cl*(f6**2)+(cs*f5*evap*(f7**2)))                                            !Litter and Soil respectively
      else
         hr = 0.0	                                                                             !Temperature above/below respiration windown
      endif
!
      return
      end
