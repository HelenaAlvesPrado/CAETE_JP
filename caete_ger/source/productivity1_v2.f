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
      subroutine productivity1 (pft,ocp_pft,LIGTH_LIMIT,temp, p0, w,
     &    wmax, ca,ipar,rh,cl1, ca1, cf1, beta_leaf, beta_awood,        !!!!!!!! rh adicionado
     &    beta_froot,emax, ph,ar, nppa, laia, f5, f1, vpd, rm, rml, 
     &    rmf, rms, rg,rgl,rgf,rgs,rc)
!
!     implicit none
!     
!     Variables
!     =========
!     
!     Input
!     -----
!     
      integer pft
      real ocp_pft,temp         !Mean monthly temperature (oC)
      real p0                   !Mean surface pressure (hPa)
      real wa,w,wmax            !Soil moisture (dimensionless)
      real ca                   !Atmospheric CO2 concentration (Pa)
      real ipar                 !Incident photosynthetic active radiation (w/m2)'
      real rh                   !Relative humidity ---------- adicionado
      real cl1, cf1, ca1
      real beta_leaf, beta_awood, beta_froot
      LOGICAL LIGTH_LIMIT
!     Output
!     ------
!     
      real ph,rc                !Canopy gross photosynthesis (kgC/m2/yr)
      real laia                 !Autotrophic respiration (kgC/m2/yr)
      real ar                   !Leaf area index (m2 leaf/m2 area)
      real nppa, vpd            !Net primary productivity (kgC/m2/yr)
      real f5
      real*16 :: f5_64
      real rm, rml, rmf, rms, rg, rgl, rgf, rgs
      
      
!     Internal
!     --------
!     
      double precision vm      !Rubisco maximum carboxylaton rate (molCO2/m2/s)

      double precision mgama    !Photo-respiration compensation point (Pa)
      REAL es, es2              !Saturation partial pressure (Pa?)
      double precision rmax     !Saturated mixing ratio (kg/kg)
      double precision r        !Moisture deficit at leaf level (kg/kg)
      double precision ci       !Internal leaf CO2 partial pressure (Pa)
      double precision a,b,c,a2,b2,c2 !Auxiliars
      double precision delta,delta2 !Auxiliars
      double precision rl       !Leaf respiration (kgC/m2/yr)
      double precision rp       !Non-leaf parts respiration (kgC/m2/yr)
      double precision sunlai,shadelai !Sunlai/Shadelai
      double precision vpd64,d  !Vapor pressure deficit (kPa)
!     
!     Rubisco, light and transport limited photosynthesis rate
!     --------------------------------------------------------
!     
      double precision jc,jl,je,jp,jp1,jp2,j1,j2 !Auxiliars
!     
!     Functions
!     ---------
!     
      double precision f1       !Leaf level gross photosynthesis (molCO2/m2/s)
      double precision f1a      !auxiliar_f1
      double precision f1b      !Leaf level gross photosynthesis (micromolCO2/m2/s)
      double precision f2       !Michaelis-Menton CO2 constant (Pa)
      double precision f3       !Michaelis-Menton O2 constant (Pa)
      double precision f4,f4sun,f4shade !Scaling-up to canopy level (dimensionless)
      
! BIANCA ___________________________________________________
      real tleaf(7)             !leaf turnover time (yr)
      real p21(7)
      double precision sla      !specific leaf area (m2/kg)      real cl2                  !leaf compartment's carbon content (kgC/m2)
      double precision csa      !sapwood compartment´s carbon content (5% of woody tissues) (kgC/m2)
      double precision ncl      !leaf N:C ratio (gN/gC)
      double precision ncf      !fine roots N:C ratio (gN/gC)
      double precision ncs      !sapwood N:C ratio(gN/gC)
      double precision csai 
      double precision pt                   !taxa potencial de fornecimento para
c     transpiração (mm/dia)
      double precision csru    !Specific root water uptake (0.5 mm/gC/dia; based in Pavlick et al (2013))
      real emax                 !potential evapotranspiration (mm/dia)
      double precision alfm     !maximum Priestley-Taylor coefficient (based in Gerten et al. 2004; used to calculate ad)
      double precision gm       !scaling conductance (mm/dia)(based in Gerten et al. 2004; used to calculate ad)
      double precision gc       !Canopy conductance (mm/dia)(based in Gerten et al. 2004; used to calculate ad)
      double precision laia64, ph64, ar64, rm64, rms64, rml64
      double precision rmf64
      double precision rg64, rgf64, rgs64, rgl64
      
c     HELENA____________________________________________________
!     Parameters 
!     ----------
!     
      DOUBLE PRECISION p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14
     &    ,p15,p19,p20,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31

!     
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
      
!     Rubisco maximum carboxylaton rate (molCO2/m2/s)
!     -----------------------------------------------
!     getting pft parameters
      call pft_par(2, p21)
      call pft_par(6, tleaf)

      vm = (p21(pft)) ! ------------ alterado! vm completamente variante. Nao é mais p21. Coloquei assim momentanemante pra facilitar

!      vm = (p21(pft)*(p22**(p10*(temp-p11))))/ !Free-range parameter --> 0.0358>vm>840 (micromol)
!     &    (1.0+exp(p23*(temp-p24)))
!
c      if (vm .gt. 0.00001) print*, vm
!     
c     call critical_value(vm)
!     Photo-respiration compensation point (Pa)
!     -----------------------------------------
!     
      mgama = p3/(p8*(p9**(p10*(temp-p11))))
!     call critical_value(mgama)
!     Michaelis-Menton CO2 constant (Pa)
!     ----------------------------------
!     
      f2 = p12*(p13**(p10*(temp-p11)))
c     call critical_value(f2)
!     
!     Michaelis-Menton O2 constant (Pa)
!     ---------------------------------
!     
      f3 = p14*(p15**(p10*(temp-p11)))
c     call critical_value(f3)
!     Saturation partial pressure of water vapour (Pa)
!     
!     ------------------------------------------------
!     
      call tetens (temp,es)     !-71.1<temp<+38.2
!     


      es2 = es*100.0
      vpd64 = (rh*es2)/1000. !----------- incorporada no umidade relativa
      vpd = real(vpd64, 4)

C      call critical_value(vpd)
!     Saturated mixing ratio (kg/kg)
!     -----------------------------
!     
      rmax = 0.622*(es/(p0-es))
c     call critical_value(rmax)
!     
!     Moisture deficit at leaf level (kg/kg)
!     --------------------------------------
!     
      r = -0.315*rmax
c     call critical_value(r)
!     
!     Internal leaf CO2 partial pressure (Pa)
!     ---------------------------------------
!     
      ci = p19* (1.-(r/p20)) * (ca-mgama) + mgama
c     call critical_value(ci)
!     Rubisco carboxilation limited photosynthesis rate (molCO2/m2/s)
!     ---------------------------------------------------------------
!     
      jc = vm*((ci-mgama)/(ci+(f2*(1.+(p3/f3)))))
c     call critical_value(jc)
!     
!     Light limited photosynthesis rate (molCO2/m2/s)
!     -----------------------------------------------
      if (LIGTH_LIMIT) then
         aux_ipar= ipar   
      else
         aux_ipar = ipar - (ipar * 0.25)
      endif
      
      jl = p4*(1.0-p5)*aux_ipar*((ci-mgama)/(ci+(p6*mgama)))
      
c     call critical_value(jl)
!     
!     Transport limited photosynthesis rate (molCO2/m2/s)
!     ---------------------------------------------------
!     
      je = p7*vm
c     call critical_value(je)
!     
!     Jp (minimum between jc and jl)
!     ------------------------------
!     
      a = 0.83
      b = (-1.)*(jc+jl)
C     call critical_value(b)
      c = jc*jl
c     call critical_value(c)
      delta = (b**2)-4.0*a*c
c     call critical_value(c)
!     
      jp1=(-b-(sqrt(delta)))/(2.0*a)
      jp2=(-b+(sqrt(delta)))/(2.0*a)
      jp= amin1(jp1,jp2)
      
C     call critical_value2(jp)
      
!     1
!     Leaf level gross photosynthesis (minimum between jc, jl and je)
!     ---------------------------------------------------------------
!     
      a2 = 0.93
      b2 = (-1.)*(jp+je)
      c2 = jp*je
      delta2 = (b2**2)-4.0*a2*c2
!     
      j1=(-b2-(sqrt(delta2)))/(2.0*a2)
      j2=(-b2+(sqrt(delta2)))/(2.0*a2)
      f1a = amin1(j1,j2)
!
!     Soil water
!     ==========      

      wa = (w/wmax) * ocp_pft
      
      call canopy_resistence(pft , vpd64, f1a, ocp_pft, rc)

      
!     Water stress response modifier (dimensionless)
!     ----------------------------------------------     

      csru = 0.5 
      pt = csru * (cf1 * 1000.) * wa  !(based in Pavlick et al. 2013; *1000. converts kgC/m2 to gC/m2)
      alfm = 1.391
      gm = 3.26 * 86400.           !(*86400 transform mm/s to mm/dia)
      
      if(rc .gt. 0.05) then
         gc = (1./rc) * 86400000. !*86400000 transfor m/s to mm/dia)
      else
         gc = (1./0.05) * 86400000. ! BIANCA E HELENA - Mudei este esquema..   
      endif                     ! tentem entender o algoritmo
                                ! e tenham certeza que faz sentido ecologico

      d =(emax*alfm)/(1. + gm/gc) !(based in Gerten et al. 2004)

!     BIanca- Eu alterei a estrutura desta equacao para evitar erros
!     Isso faz a mesma coisa que o calculo que vc implementou - jp
      if(d .gt. 0.0) then
         f5_64 = pt/d
         f5_64 = exp(-1 * f5_64)
         f5_64 = 1.0 - f5_64
      else
         f5_64 = wa  ! eu acrescentei esta parte caso d seja igual a zero
      endif          ! nao sei se faz sentido!

      f5 = real(f5_64,4) !esta funcao transforma o f5 (double precision)
!     em real (32 bits) 
      
!     Photosysthesis minimum and maximum temperature
!     ----------------------------------------------
      
      if ((temp.ge.-10.0).and.(temp.le.50.0)) then
         f1 = f1a*real(f5_64,8) !f5:water stress factor-- Notem que aqui a tranformacao eh de 128 pra 64 bits
      else
         f1 = 0.0               !Temperature above/below photosynthesis windown
      endif

!     Leaf area index (m2/m2)
!     ---------------------------------
      sla=((0.0300*1000.)*((365./(((tleaf(pft))/365.)/12.))**(-0.46)))    
      laia64 = (cl1 * 365.0 * sla) * ocp_pft    !boa!
      laia = real(laia64,4)
      
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

      laia  = real((f4sun + f4shade), 4) ! pra mim faz sentido que a laia final seja
                                         ! a soma da lai em nivel de dossel (sun + shade) - jp
!     Canopy gross photosynthesis (kgC/m2/yr)
!     =======================================
!     (0.012 converts molCO2 to kgC)
!     (31557600 converts seconds to year [with 365.25 days])
      ph64 = 0.012*31557600.0*f1*f4sun*f4shade
      ph = real(ph64, 4)       ! kg m-2 year-1

c     ============================================================
c     Autothrophic respiration
!     ========================
!     Maintenance respiration (kgC/m2/yr) (based in Ryan 1991)
      csa= 0.05 * ca1           !sapwood carbon content (kgC/m2). 5% of woody tissues (Pavlick, 2013)
!     Atencao
!     eu alterei os valores da razao c:n ...os dados foram chutados por mim,
!     com base numa leitura superficial do seguinte paper:
      
!     http://onlinelibrary.wiley.com/wol1/doi/10.1111/gcb.12022/suppinfo
      
!     Notem que eu li apenas o material suplementar.
!     Assumi que a massa de carbono em todos os tecidos
!     eh de 50%! Quem concorda? por favor nao coloquem acentuacao
!     grafica nos comments -- eu sei que eh dificil,
!     eu esqueço as vezes =P. Mas eu desconfio que essa eh a origem do bug que
!     nao deixa a bia editar o codigo no notepad++ 
      
      ncl = 0.021               !(gN/gC) 
      ncf = 0.004               !(gN/gC)
      ncs = 0.003               !(gN/gC)

      
      rml64 = (ncl * cl1) * 27. * exp(0.03*temp)
      rml =  real(rml64,4) * ocp_pft

      rmf64 = (ncf * cf1) * 27. * exp(0.03*temp)
      rmf =  real(rmf64,4) * ocp_pft

      rms64 = (ncs * csa) * 27. * exp(0.03*temp)
      rms = real(rms64,4) * ocp_pft
       
      rm64 = (rml64 + rmf64 + rms64)
      rm = real(rm64, 4)

c     Growth respiration (KgC/m2/yr)(based in Ryan 1991; Sitch et al.
c     2003; Levis et al. 2004)         
       
      csai =  (beta_awood * 0.05)

      rgl64 = 0.25 * beta_leaf * 365.
      rgl = real(rgl64,4) * ocp_pft

      rgf64 =  0.25* beta_froot * 365.
      rgf = real(rgf64,4) * ocp_pft

      rgs64 = (0.25 * csai * 365.)
      rgs = real(rgs64,4) * ocp_pft
     
      rg64 = (rgl64 + rgf64 + rgs64)
      rg = real(rg64,4) * ocp_pft
     
      if (rg.lt.0) then
         rg = 0.0
      endif
      
!     c Autotrophic (plant) respiration -ar- (kgC/m2/yr)
!     Respiration minimum and maximum temperature
!     -------------------------------------------     
      if ((temp.ge.-10.0).and.(temp.le.50.0)) then
         ar64 = rm64 + rg64
         ar = real(ar64,4) 
      else
         ar = 0.0               !Temperature above/below respiration windown
      endif
      
c     -----------------------------------------------------------------
!     NPP
!     ============
!     Productivity
!     ============
!     Net primary productivity(kgC/m2/yr)
!     ====================================
      nppa = ph - ar
      !if(nppa .lt. 0.0) nppa = 0.0 
!     No futuro proximo poderiamos usar uma npp negativa para indicar uma perda de carbono
!     dos tecidos vegetais... e uma possivel extincao do pft/pls da celula de grid.
      return
      end subroutine productivity1



      
!     ==================================================================
      subroutine canopy_resistence(m,vpd_in,f1_in,ocp_pft1,rc2_in)
      
!     Variables
!     =========
      
!     Inputs
!     ------      
      integer m
      double precision f1_in    !Photosynthesis (molCO2/m2/s)
      double precision vpd_in   !kPa
      
!     Outputs
!     -------
      real rc2_in               !Canopy resistence (s/m)
      
!     Internal
!     --------
      real g1(7) 
      double precision f1b      !Photosynthesis (micromolCO2/m2/s)
      double precision gs2      !Canopy conductance (m/s)
      double precision gs       !Canopy conductance (molCO2/m2/s)
      double precision g0       !Residual stomatance conductance
      double precision D1       !kPA
      double precision aa

      call pft_par(1, g1)
      
      f1b = (f1_in*10e5)        ! Helena - Mudei algumas coisas aqui
      aa = (f1b/363.)           ! Entenda o algoritmo e tenha certeza de que  
      g0 = 0.01                 ! condiz com a realidade esperada =)
      
      if(vpd_in .gt. 0.0) then
         goto 10
      else
         gs = 1.5
         goto 100
      endif
 10   continue
      if (vpd_in.lt.0.25) then
         gs = 1.5
      else
         D1 = sqrt(vpd_in)
         gs = g0 + 1.6 * (1. + (g1(m)/D1)) * (aa) !Based on Medlyn et al. 2011
      endif
 100  continue
      gs2 = gs/41.     
      rc2_in = real((1./gs2),4) * ocp_pft1
      
      return
      end subroutine canopy_resistence
      
!     =====================================
!     Microbial (heterotrophic) respiration
!     =====================================
      
      subroutine carbon2 (tsoil,f5,evap,laia, !Inputs
     &    cl,cs,hr)             !Outputs
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

!     xOutputs 
!     -------
!     
      real cl                   !Litter carbon (kgC/m2)
      real cs                   !Soil carbon (kgC/m2)
      real hr                   !Heterotrophic (microbial) respiration (kgC/m2)

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
      f6 = 1.16*10.**(-1.4553+0.0014175*(evap*365.0))
C      call critical_value(f6)
!     
!     Soil carbon storage function                                          !Controlled by temperature
!     ----------------------------
!     
      f7 = p32**(p10*(tsoil-p11))
C      call critical_value(f7)
!     Litterfall (kgC/m2)
!     ------------------
!     
      lf = p33 * laia
C      call critical_value(lf)
!     
!     Litter carbon (kgC/m2)
!     ----------------------
!     
      cl = lf/f6
C      call critical_value(cl)
!     
!     Soil carbon(kgC/m2)
!     -------------------
!     
      cs = ((p34*cl)/(p35*f7))*f5
C      call critical_value(cs)
!     
!     Respiration minimum and maximum temperature
!     -------------------------------------------
!     
      if ((tsoil.ge.-10.0).and.(tsoil.le.50.0)) then
         hr = p36*(cl*(f6**2)+(cs*f5*evap*(f7**2))) !Litter and Soil respectively
C         call critical_value(hr)
      else
         hr = 0.0               !Temperature above/below respiration windown
      endif
!     
      return
      end
!     ===================================================

!     Estas 3 subrotinas critical_value eu utilizei durante o debug
!     Ainda tem algumas espalhadas pelo codigo, pretendo
!     ir retirando aos poucos agora, mas no futuro eu vou utiliza-las
      
      subroutine critical_value(var)
      implicit none
      real var
      
      if(abs(var) .lt. 0.000001) var = 0.0

      return
      end subroutine critical_value


      
      subroutine critical_value2(var1)
      implicit none
      double precision var1
      
      if(abs(var1) .lt. 0.00000001) var1 = 0.0

      return
      end subroutine critical_value2


      
      subroutine critical_value3(var2)
      implicit none
      real*16 var2
      
      if(abs(var2) .lt. 0.000000000001) var2 = 0.0

      return
      end subroutine critical_value3


c23456
      SUBROUTINE PFT_AREA_FRAC(CLEAF, CFROOT, CAWOOD, OCP_COEFFS,
     &    OCP_WOOD)
c     esta subrotina calcula a areA DE CADA pft  
c     a partir do conteudo de carbono nos compartimentos vegetais no 
c     dia anterior 
      INTEGER, PARAMETER :: NPFT = 7
      INTEGER :: P, MAX_INDEX(1), I
      REAL :: CLEAF(NPFT), CFROOT(NPFT), CAWOOD(NPFT)
      REAL :: TOTAL_BIOMASS_PFT(NPFT),OCP_COEFFS(NPFT)
      REAL :: TOTAL_BIOMASS, TOTAL_WOOD, TOTAL_W_PFT(NPFT)
      LOGICAL :: OCP_WOOD(NPFT)
      
      TOTAL_BIOMASS = 0.0
      TOTAL_WOOD = 0.0
      do p = 1,npft
         TOTAL_W_PFT(P) = 0.0
         TOTAL_BIOMASS_PFT(P) = 0.0
         OCP_COEFFS(P) = 0.0
         OCP_WOOD(P) = .FALSE.
      enddo
      

      DO P = 1,NPFT
         TOTAL_BIOMASS_PFT(P) = CLEAF(P) + CFROOT(P) + CAWOOD(P) ! biomassa total no dia i
         TOTAL_BIOMASS = TOTAL_BIOMASS + TOTAL_BIOMASS_PFT(P)
         TOTAL_WOOD = TOTAL_WOOD + CAWOOD(P)
         TOTAL_W_PFT(P) = CAWOOD(P)
      ENDDO                 ! end p loop

!     GRID CELL OCCUPATION COEFFICIENTS
      IF(TOTAL_BIOMASS .GT. 0.0) THEN
         DO P = 1,NPFT   
            OCP_COEFFS(P) = TOTAL_BIOMASS_PFT(P) / TOTAL_BIOMASS
            IF(OCP_COEFFS(P) .LT. 0.0) OCP_COEFFS(P) = 0.0
c            CALL CRITICAL_VALUE(OCP_COEFFS(P))
c            PRINT*, OCP_COEFFS(P), 'ocp'
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
      
C23456
      RETURN
      END SUBROUTINE PFT_AREA_FRAC






!     =========================================================
!     
      subroutine penman (spre,temp,ur,rn,rc2,evap)
!     Desconfio seriamente que temos que revisar esta subrotina
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
      real ra, h5, t1, t2, es, es1, es2, delta_e
      real gama, gama2
      ra = 100.                  !s/m
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
      real ra, rcmin, t1, t2, es, es1, es2, delta_e
      real gama, gama2
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
      subroutine runoff (wa,roff)
      real*8 :: wa
      real :: roff
      real*16 :: roff64
      roff64 = 38.*(wa**11.)
c      roff64 = 11.5*(wa**6.6) * 1000. !From NCEP-NCAR Reanalysis data
      call critical_value3(roff64)
      roff = real(roff64, 4)
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
c=====================================================================
c     
c     subroutine allocation calculates the daily carbon content of each
c     compartment
c     
c     code written by Bianca Rius & David Lapola (27.Ago.2015)
c     
c=====================================================================
      
      subroutine allocation (pft, npp ,scl1,sca1,scf1,
     &    scl2,sca2,scf2)          !output
c     
c     
!     variables
      integer, parameter :: npfts = 7
      integer pft   
      real npp                  !potential npp (KgC/m2/yr)
      real npp_aux              !auxiliary variable to calculate potential npp in KgC/m2/day
      real scl1                  !previous day carbon content on leaf compartment (KgC/m2)
      real scl2                  !final carbon content on leaf compartment (KgC/m2)
      real sca1                  !previous day carbon content on aboveground woody biomass compartment(KgC/m2)
      real sca2                  !final carbon content on aboveground woody biomass compartment (KgC/m2)
      real scf1                  !previous day carbon content on fine roots compartment (KgC/m2)
      real scf2                  !final carbon content on fine roots compartment (KgC/m2)      
      
      real aleaf(npfts)             !npp percentage allocated compartment
      real aawood(npfts)
      real afroot(npfts)
      real tleaf(npfts)             !turnover time (yr)
      real tawood(npfts)
      real tfroot(npfts)            

      call pft_par(3, aleaf)
      call pft_par(4, aawood)
      call pft_par(5, afroot)
      call pft_par(6, tleaf)
      call pft_par(7, tawood)
      call pft_par(8, tfroot)
    
!     Carbon content of each compartment(KgC/m2)
c     
c     
c     initialization
c      if((scl1 .lt. 0.0000001) .or. (scf1 .lt. 0.0000001)) then
c         IF(NPP .lt. 0.0000001) THEN
c            scl2 = 0.0
c            scf2 = 0.0
c            sca2 = 0.0 
c            goto 10
c         ENDIF
c      endif   
      npp_aux = npp/365.0       !transform (KgC/m2/yr) in (KgC/m2/day)
      scl2 = scl1 + (aleaf(pft) * npp_aux) -(scl1 /(tleaf(pft)*365.0))
      scf2 = scf1 +(afroot(pft) * npp_aux)-(scf1 /(tfroot(pft)*365.0))
      if(aawood(pft) .gt. 0.0) then
         sca2 = sca1 +(aawood(pft)*npp_aux)-(sca1/(tawood(pft)*365.0))
      else
         sca2 = 0.0
      endif


      if(scl2 .lt. 0.0) scl2 = 0.0
      if(scf2 .lt. 0.0) scf2 = 0.0
      if(sca2 .lt. 0.0) sca2 = 0.0
      
c     10   continue
      return
      end
c    