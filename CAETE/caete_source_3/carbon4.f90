! this code is based on CPTEC-PVM2 source code. Wrote by Carlos Nobre, 
! Marcos Oyama, David Lapola, Bianca Rius e Helena Alves do Prado


! author :: jpdarela 
! in an attempt to write CAETE calculation routines


subroutine prod (pft, temp, p0, w, wmax, ca, ipar, tsoil,&
     cl1, ca1, cf1, beta_leaf, beta_awood, beta_froot, ph,&
     ar, nppa, laia, f5, f1, vpd, rm, rml, rmf, rms, rg, rgl,&
     rgf,rgs)
  implicit none
  
  
  !c i/o variables
  integer, intent(in) :: pft
  real,intent(in) :: temp  !mean monthly temperature (oC)
  real,intent(in) :: p0  !mean surface pressure (hPa)
  real,intent(in) :: w,wmax !soil moisture (dimensionless)
  real,intent(in) :: ca  !Atmospheric CO2 pressure (Pa)
  real,intent(in) :: ipar  !incident photos. active radiation
  
  real, intent(in) :: beta_leaf
  real, intent(in) :: beta_awood
  real, intent(in) :: beta_froot
  real, intent(in) :: cl1
  real, intent(in) :: ca1
  real, intent(in) :: cf1
  
  
  real, intent(out) :: ph  !Canopy gross photosyntheis(kgC/m2/yr)
  real, intent(out) :: laia !Leaf area index (m2 leaf/ m2 area)
  real, intent(out) :: ar !Autotrophic (plant) respiration (kgC/m2/yr)
  real, intent(out) :: nppa !Net primary productivity (kgC/m2/yr)
  real, intent(out) :: f5
  real, intent(out) :: f1
  real, intent(out) :: vpd
  real, intent(out) :: rm
  real, intent(out) :: rml
  real, intent(out) :: rms
  real, intent(out) :: rmf
  real, intent(out) :: rg
  real, intent(out) :: rgf
  real, intent(out) :: rgl
  real, intent(out) :: rgs
  
  !c internal variables
  
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
  
  
  !  c     HELENA____________________________________________________
  !     Parameters 
  !     ----------
  !     
  real p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p19,p20,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31
  real p21(3)
  data p21 /0.00003,0.000058,0.000067/ !Tropical/Temperate/Boreal   
  real wa
  real vm,mgama,es,es2,rmax,r,ci,a,b,c,a2,b2,c2,delta,delta2,f1a,f2,f3,f4,k16,k17,rl,rp
  real f1b
  real f10 !CO2 sensitivity
  real sunlai, shadelai, f4shade, f4sun
  !c Rubisco, light and transport limited photosynthesis rate respectively
  real jc,jl,je,jp,jp1,jp2,j1,j2 !auxiliars

  
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
  !p21 = 0.00004            !Maximum Rubisco carboxylation rate (molCO2/m2/s)
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
  

  !c
  !c
  !c
  !c=======================================================================
  !c Photosynthesis =======================================================
  !c
  !c Rubisco maximum carboxylaton rate
  !c (vm ; molCO2/m^2s) [Eq. 12]
  vm = (p21(pft)*p22**(p10*(temp-p11)))/(1.0+exp(p23*(temp-p24)))
  !vm = (0.00004*2.0**(0.1*(temp-25.0)))/(1.0+exp(0.3*(temp-36.0)))
  call critical_value(vm)
  !c Photo-respiration compensation point
  !c (mgama ; Pa) [Eq. 8]
  !c
  mgama = p3/(p8*(p9**(p10*(temp-p11))))
  call critical_value(mgama)
  !mgama = 21200.0/(5200.0*(0.57**(0.1*(temp-25.0))))
  
  
  !     Michaelis-Menten CO2 constant (Pa)
  !     ----------------------------------
  !     
  f2 = p12*(p13**(p10*(temp-p11)))
  call critical_value(f2)
  !     
  !     Michaelis-Menten O2 constant (Pa)
  !     ---------------------------------
  !     
  f3 = p14*(p15**(p10*(temp-p11)))
  call critical_value(f3)
  !     Saturation partial pressure of water vapour (Pa)                                       
  !     -----------------------------------------------
  !call tetens (temp,es)
  if (temp.ge.0.) then
     es = 6.1078*exp((7.5*temp/(237.3+temp))*log(10.))
  else
     es = 6.1078*exp((9.5*temp/(265.5+temp))*log(10.))
  endif
  
  call critical_value(es)
  es2 = es*100.0
  vpd = (((100.0-68.5)/100.0)*(es2))/1000.0 !kPa
  call critical_value(vpd)
  !c
  
  !     Saturated mixing ratio (kg/kg)
  !     ------------------------------
  !     
  rmax = 0.622*(es/(p0-es))
  call critical_value(rmax)
  !     
  !     Moisture deficit at leaf level (kg/kg)
  !     --------------------------------------
  !     
  r = -0.315*rmax
  call critical_value(r)
  
  ci = p19*(1-(r/p20))*(ca-mgama)+mgama
  call critical_value(ci)
  !     Rubisco carboxilation limited photosynthesis rate (molCO2/m2/s)
  !     ---------------------------------------------------------------
  !     
  jc = vm*((ci-mgama)/(ci+(f2*(1+(p3/f3)))))
  call critical_value(jc)
  !     
  !     Light limited photosynthesis rate (molCO2/m2/s)
  !     -----------------------------------------------
  !     
  jl = p4*(1.0-p5)*ipar*((ci-mgama)/(ci+(p6*mgama)))
  call critical_value(jl)
  
  
  je = p7*vm
  call critical_value(je)
  !     
  !     Jp (minimum between jc and jl)
  !     ------------------------------
  !     
  a = 0.83
  b = (-1)*(jc+jl)
  call critical_value(b)
  c = jc*jl
  call critical_value(c)
  delta = (b**2)-4.0*a*c
  call critical_value(c)
  !     
  jp1=(-b-(sqrt(delta)))/(2.0*a)
  jp2=(-b+(sqrt(delta)))/(2.0*a)
  jp= amin1(jp1,jp2)
  
  call critical_value(jp)
  
  !     1
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
  call critical_value(f1a)
  !  c      PRINT*, F1A, 'f1a'
  
  !     Soil water
  !     ==========
  !     
  wa = w/wmax  !c
  !c Water stress response modifier (dimensionless)
  !c [f5 ; Eq. 21]
  !c 
  if (wa.gt.0.5) f5 = 1.0  !Not too lower in e.g. Amazonian dry season
  if ((wa.ge.0.205).and.(wa.le.0.5)) f5 = (wa-0.205)/(0.5-0.205)
  if (wa.lt.0.205) f5 = wa  !Below wilting point f5 accompains wa (then Sahara is well represented)
  !c
  !c ...f1
  !c
  !c Photosysthesis minimum and maximum temperature
  if ((temp.ge.-10.0).and.(temp.le.50.0)) then
     f1 = f1a*f5            !f5:water stress factor
     call critical_value(f1)
  else
     f1 = 0.0               !Temperature above/below photosynthesis windown
  endif
  
  
  !     Leaf area index (m2 leaf/m2 arEa) bianca
  !     ---------------------------------
  sla=(0.030*1000.)*((365/(((tleaf(pft))*365)/12))**(-0.46))
  !  c      PRINT*, sla, 'sla'
  laia  = 0.25*exp(2.5*(f1/p25)) !Adjusted after using observed ipar
  !  c      if(cl1 .gt. 0) print*, cl1, 'cl1'
  !  c      laia = (cl1*365.*sla)
  !  c      if(laia .gt. 0) PRINT*, laia, 'laia'
  !     SunLAI
  !     ------
  !     
  sunlai = (1.0-(exp(-p26*laia)))/p26
  call critical_value(sunlai)
  !     ShadeLAI
  !     --------
  !     
  shadelai = laia - sunlai
  call critical_value(shadelai)
  !     
  !     Scaling-up to canopy level (dimensionless)
  !     ------------------------------------------
  !     
  f4 = (1.0-(exp(-p26*laia)))/p26 !Sun 90 degrees in the whole canopy, to be used for respiration
  call critical_value(f4)
  
  !     Sun/Shade approach to canopy scaling                                  !Based in de Pury & Farquhar (1997)
  !     ------------------------------------
  !     
  f4sun = (1.0-(exp(-p26*sunlai)))/p26 !sun 90 degrees
  f4shade = (1.0-(exp(-p27*shadelai)))/p27 !sun ~20 degrees
  !
  call critical_value(f4sun)
  call critical_value(f4shade)
  !     Canopy gross photosynthesis (kgC/m2/yr)
  !     =======================================
  !     (0.012 converts molCO2 to kgC)
  !     (31557600 converts seconds to year [with 365.25 days])
  !     
  ph = 0.012*31557600.0*f1*f4sun*f4shade
  call critical_value(ph)
  !  c      PRINT*, PH, 'ph'
  
  
  !     Plant respiration
  !     =================
  !     c Maintenance respiration (kgC/m2/yr) (based in Ryan 1991)
      
      
  csa= 0.05*ca1             !sapwood carbon content (kgC/m2). 5% of woody tissues (Pavlick, 2013)
  call critical_value(csa)
  ncl = 0.034               !(gN/gC)
  ncf = 0.034               !(gN/gC)
  ncs = 0.003               !(gN/gC)
  rml = (ncl*cl1)*27*(exp(0.03*temp))
  call critical_value(rml)
  rmf = (ncf*cf1)*27*(exp(0.03*tsoil))
  call critical_value(rmf)
  rms = (ncs*csa)*27*(exp(0.03*temp))
  call critical_value(rms)
  
  rm = rml + rmf + rms
  call critical_value(rm)
  !  c      print*, rm, 'rm'
  
  ! c     Growth respiration (KgC/m2/yr)(based in Ryan 1991; Sitch et al. 2003; Levis et al. 2004)         
  
  csai= 0.05*beta_awood
  call critical_value(csai)
  rgl = (0.25*((beta_leaf)*365))
  call critical_value(rgl)
  rgf = (0.25*((beta_froot)*365))
  call critical_value(rgf)
  rgs = (0.25*(csai)*365)
  call critical_value(rgs)
  
  rg = rgl + rgf + rgs
  call critical_value(rg)
  !  c      print*, rg, 'rg'
  if (rg.lt.0) then
     rg = 0.
  endif
  
  
  !     c Autotrophic (plant) respiration -ar- (kgC/m2/yr)
  !     Respiration minimum and maximum temperature
  !     -------------------------------------------
  !     
  if ((temp.ge.-10.0).and.(temp.le.50.0)) then
     ar = rm+rg
     call critical_value(ar)
     !     c         if (ar .gt. 0.)PRINT*, AR ,'ar'
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
  call critical_value(nppa)
  
  !  c       if (nppa .gt. 0.) print*, nppa, 'npp'
  !  c      PRINT*, NPPA
  !  c      if (nppa.lt.0.0) nppa = 0.0 !Maybe this is incorrect, but demands retuning of every biome limits
  !     
  return
  
end subroutine prod






subroutine carbon_hr (tsoil,f5,evap,laia,cl,cs,hr)
  
  implicit none
  !c i/o variables
  real, intent(in) :: tsoil !mean monthly soil temperature (oC)
  real, intent(in) :: f5  !stress response to soil moisture (dimensionless)
  real, intent(in) :: evap  !actual evapotranspiration (mm/day)
  real, intent(in) :: laia  !Leaf area index (m2 leaf/ m2 area)
  !c
  real, intent(out) :: cl  !Litter carbon (kgC/m2)
  real, intent(out) :: cs  !Soil carbon (kgC/m2)
  real, intent(out) :: hr  !Heterotrophic (microbial) respiration (kgC/m2/yr)
  !c
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
  call critical_value(f6)
  !     
  !     Soil carbon storage function                                          !Controlled by temperature
  !     ----------------------------
  !     
  f7 = p32**(p10*(tsoil-p11))
  call critical_value(f7)
  !     Litterfall (kgC/m2)
  !     ------------------
  !     
  lf = p33*laia
  call critical_value(lf)
  !     
  !     Litter carbon (kgC/m2)
  !     ----------------------
  !     
  cl = lf/f6
  call critical_value(cl)
  !     
!     Soil carbon(kgC/m2)
  !     -------------------
  !     
  cs = ((p34*cl)/(p35*f7))*f5
  call critical_value(cs)
  !     
  !     Respiration minimum and maximum temperature
  !     -------------------------------------------
  !     
  if ((tsoil.ge.-10.0).and.(tsoil.le.50.0)) then
     hr = p36*(cl*(f6**2)+(cs*f5*evap*(f7**2))) !Litter and Soil respectively
     call critical_value(hr)
  else
     hr = 0.0               !Temperature above/below respiration windown
  endif
  !     
  return
end subroutine carbon_hr
    
    
    
    

subroutine critical_value(var)
  implicit none
  real,intent(inout) :: var
  
  if(abs(var) .lt. 0.000000001) var = 0.0
  
  return
end subroutine critical_value





subroutine penman (spre,temp,ur,rn,rc2,evap)
  
  implicit none
  
  !c
  !c Entradas
  !!c --------
  !c spre   = pressao aa supeficie (mb)
  !c temp   = temperatura (oC)
  !c w      = grau de saturacao (0-1,adimensional)
  !c ur     = umidade relativa  (0-1,adimensional)
  !c rn     = saldo de radiacao (W m-2)
  !c rc2    = resistencia do dossel (s/m)
  !c
  !c Saida
  !c -----
  !c evap  = evapotranspiracao (mm/dia)

  ! parametros
  real, parameter ::  ra = 100.0   !s/m
  real, parameter ::  h5 = 0.0275 !mb-1

  ! i/o
  real, intent( in) :: spre
  real, intent( in) :: temp
  !real, intent( in) :: w
  !real, intent( in) :: wmax
  real, intent( in) :: ur
  real, intent( in) :: rn
  real, intent( in) :: rc2
  
  real, intent(out) :: evap

  !c internal var
  real :: es, es1, es2, t1, t2
  real :: delta, delta_e, gama, gama2
    
  !c delta
  t1 = temp + 1.
  t2 = temp - 1.

  !call tetens(t1,es1)
  if (t1.ge.0.) then
     es1 = 6.1078*exp((7.5*t1/(237.3+t1))*log(10.))
  else
     es1 = 6.1078*exp((9.5*t1/(265.5+t1))*log(10.))
  endif
  
  !call tetens(t2,es2)
    if (t2.ge.0.) then
     es2 = 6.1078*exp((7.5*t2/(237.3+t2))*log(10.))
  else
     es2 = 6.1078*exp((9.5*t2/(265.5+t2))*log(10.))
  endif
  
  delta = (es1-es2)/(t1-t2) !mb/oC
  
  !c delta_e
  ! call tetens (temp,es)
  if (temp.ge.0.) then
     es = 6.1078*exp((7.5*temp/(237.3+temp))*log(10.))
  else
     es = 6.1078*exp((9.5*temp/(265.5+temp))*log(10.))
  endif
  
  delta_e = es*(1. - ur) !mb
  
  if ((delta_e.ge.(1./h5)-0.5).or.(rc2.ge.4500)) evap = 0.
  
  if ((delta_e.lt.(1./h5)-0.5).or.(rc2.lt.4500)) then
     !c gama e gama2
     gama  = spre*(1004.)/(2.45e6*0.622)
     gama2 = gama*(ra + rc2)/ra
     !c evapotranspiracao real
     evap = (delta*rn + (1.20*1004./ra)*delta_e)/(delta+gama2) ! W/m2
     evap = evap*(86400./2.45e6)!mm/dia
     !evap = amax1(evap,0.) !elimina condensacao
  endif
end subroutine penman






subroutine evpot2 (spre,temp,ur,rn,evap)
  !c
  !c Entradas
  !c --------
  !c spre   = pressao aa supeficie (mb)
  !c temp   = temperatura (oC)
  !c ur     = umidade relativa  (0-1,adimensional)
  !c rn     = saldo de radiacao (W m-2)
  !c
  !c Saida
  !c -----
  !c evap  = evapotranspiracao potencial(mm/dia)
  !c
  
  !c parametros
  real, parameter :: ra =    100.   !s/m
  real, parameter :: rcmin = 100.   !s/m

  ! i/o
  real, intent( in) :: spre
  real, intent( in) :: temp
  real, intent( in) :: ur
  real, intent( in) :: rn
  real, intent(out) :: evap
  
  ! internal vars
  real :: t1,t2,es1,es2,gama, gama2

  
  !c delta
  t1 = temp + 1.
  t2 = temp - 1.
  !call tetens(t1,es1)
  if (t1.ge.0.) then
     es1 = 6.1078*exp((7.5*t1/(237.3+t1))*log(10.))
  else
     es1 = 6.1078*exp((9.5*t1/(265.5+t1))*log(10.))
  endif
  
  !call tetens(t2,es2)
    if (t2.ge.0.) then
     es2 = 6.1078*exp((7.5*t2/(237.3+t2))*log(10.))
  else
     es2 = 6.1078*exp((9.5*t2/(265.5+t2))*log(10.))
  endif
  delta = (es1-es2)/(t1-t2) !mb/oC
  !c
  !c delta_e
  !call tetens (temp,es)
  if (temp.ge.0.) then
     es = 6.1078*exp((7.5*temp/(237.3+temp))*log(10.))
  else
     es = 6.1078*exp((9.5*temp/(265.5+temp))*log(10.))
  endif
  delta_e = es*(1. - ur) !mb
  !c
  !c resistencia estomatica
  rc = rcmin
  !c
  !c gama e gama2
  gama  = spre*(1004.)/(2.45e6*0.622)
  gama2 = gama*(ra + rc)/ra
  !c
  !c evapotranspiracao potencial sem estresse
  evap = (delta*rn + (1.20*1004./ra)*delta_e)/(delta+gama2) ! W/m2
  evap = evap*(86400./2.45e6)                               ! mm/dia
  !      evap = amax1(evap,0.) !elimina condensacao
  !c
end subroutine evpot2






subroutine soil_temp(temp, tsoil)
  ! Calcula a temperatura do solo. Aqui vamos mudar no futuro!
  ! a tsoil deve ter relacao com a et realizada...
  ! a profundidade do solo (H) e o coef de difusao (DIFFU) devem ser variaveis (MAPA DE SOLO?; agua no solo?)  
  implicit none

  !parameters
  real, parameter :: H = 1.0                         ! soil layer thickness (meters)
  real, parameter :: DIFFU = 4.e7 * (30.0 * 86400.0) ! soil thermal diffusivity (m2/mes)
  real, parameter :: TAU = (H ** 2) / (2.0 * DIFFU)  ! e-folding times (months) 
  integer, parameter :: m = 12

  ! i/o
  
  real,dimension(m), intent( in) :: temp ! future __ make temps an allocatable array
  real,dimension(m), intent(out) :: tsoil
 
  ! internal vars
  
  integer :: n, k
  real :: t0 = 0.0
  real :: t1 = 0.0

  tsoil = -9999.0

  do n=1,1200 !run to attain equilibrium
     k = mod(n,m)
     if (k.eq.0) k = 12
     t1 = (t0*exp(-1.0/TAU) + (1.0 - exp(-1.0/TAU)))*temp(k)
     tsoil(k) = (t0 + t1)/2.0
     t0 = t1
  enddo
end subroutine soil_temp






subroutine canopy_resistence(m,vpd,f1,rc2)
  integer, intent(in) :: m
  real, intent( in) :: vpd, f1
  real, intent(out) :: rc2
  
  real f1b                  !Photosynthesis (micromolCO2/m2/s)
  real gs2                  !Canopy conductance (m/s)
  real gs                   !Canopy conductance (molCO2/m2/s)
  real g0                   !Residual stomatance conductance
  real g1(3)                !3.00 (boreal) / 2.05 (temperate) / 3.02 (tropical)
  !     real g1 
  real D                    !kPA
  real aa
  real rcmax 
  
  data g1 /6.0,4.0,2.0/ 
  f1b = (f1*10e5)           !maior f1b = 13.38
  aa = (f1b/363)
  g0 = 0.01
  !     g1 = 4.9
  D = sqrt(vpd)
  !     rcmax = 200.0                                                    !A=2.0
  rcmax = 550.0             !A=0.5                             
  !     rcmax = 1800.0      
  if (vpd.lt.0.25) gs = 1.5
  if (vpd.ge.0.25) then 
     gs = g0 + 1.6 * (1 + (g1(m)/D)) * (aa) !Based on Medlyn et al. 2011
     !     gs = g0 + 1.6 * (1 + (g1/D)) * (aa)                               !Based on Medlyn et al. 2011
     call critical_value(gs)
  endif
  !     
  !     if (gs.le.0.01) gs = 0.0                                        
  !     if (gs.le.0.0) then                                             !-0.029<gs<1.29
  !     print*,gs
  !     endif
  !     
  gs2 = gs/41
  call critical_value(gs2)
  !     if (gs2.le.0.00025) gs2 = 0.0
  !     if (gs2.gt.0.0318) then                                         !-0.00069<gs2<0.0317
  !     print*,print
  !     endif
  !     
  rc2 = 1/gs2
  call critical_value(rc2)
  !     if (rc2.ge.151.47) rc2 = rcmax                                    !A=2
  if (rc2.ge.545.43) rc2 = rcmax !A=0.5  
  !     if (rc2.ge.1779.97) rc2 = rcmax                                   !A=0.1
  !     
  !     if (rc2.lt.19.57) then
  !     print*,rc2
  !     endif
  !     
  return
end subroutine canopy_resistence





subroutine available_energy(temp, ae)
  real, intent( in) :: temp
  real, intent(out) :: ae
  
  ae = 2.895*temp + 52.326 !from NCEP-NCAR Reanalysis data
  
end subroutine available_energy




subroutine runoff_c(w, wmax, roff)
  real, intent( in) :: w, wmax
  real, intent(out) :: roff
  
  roff = 38.*((w/wmax)**11.) ! [Eq. 10]
  roff = 11.5*((w/wmax)**6.6) !from NCEP-NCAR Reanalysis data 
end subroutine runoff_c




!c234567
!c=====================================================================
!c
!c subroutine allocation calculates the daily carbon content of each
!c compartment
!c
!c code written by Bianca Rius & David Lapola (27.Ago.2015)
!c
!c=====================================================================
!      
subroutine allocation (i6,npp,cl1,ca1,cf1, cl2,ca2,cf2)
  implicit none
  !c     
  !!     variables
  
  integer, intent(in) :: i6                !index for calculate the 3 PFTs (i6.eq.1 - TBE; i6.eq.2 - TBD; i6.eq.3 - HERB)
  real, intent( in)   :: npp                  !potential npp (KgC/m2/yr)
  real npp_aux              !auxiliary variable to calculate potential npp in KgC/m2/day
  real, intent( in)   ::  cl1                  !previous day carbon content on leaf compartment (KgC/m2)
  real, intent(out)   ::  cl2                  !final carbon content on leaf compartment (KgC/m2)
  real, intent( in)   ::  ca1                  !previous day carbon content on aboveground woody biomass compartment(KgC/m2)
  real, intent(out)   ::  ca2                  !final carbon content on aboveground woody biomass compartment (KgC/m2)
  real, intent( in)   ::  cf1                  !previous day carbon content on fine roots compartment (KgC/m2)
  real, intent(out)   ::  cf2                  !final carbon content on fine roots compartment (KgC/m2)
  
  
  
  
  
  real aleaf(3)             !npp percentage alocated to leaf compartment
  data aleaf /0.35, 0.25, 0.45/
  real aawood (3)           !npp percentage alocated to aboveground woody biomass compartment
  data aawood /0.40, 0.40, 0.0/
  real afroot(3)            !npp percentage alocated to fine roots compartment
  data afroot /0.25, 0.35, 0.55/ 
  real tleaf(3)             !turnover time of the leaf compartment (yr)
  data tleaf /1.0, 0.5, 1.0/ 
  real tawood (3)           !turnover time of the aboveground woody biomass compartment (yr)
  data tawood /31.0, 33.0, 0.0/
  real tfroot(3)            !turnover time of the fine roots compartment
  data tfroot /3.0, 2.0, 1.0/
  
  
!!$  
!!$!  c     
!!$!  c=====================================================================
!!$!  c     
!!$!  c     
!!$  !     Carbon content of each compartment(KgC/m2)
!!$  c     
!!$  c     
!!$c     initialization
!!$      
  npp_aux = npp/365.     !transform (KgC/m2/yr) in (KgC/m2/day)
  
  
  cl2 = ((aleaf(i6)*npp_aux) - (cl1/((tleaf(i6))*365)))+ cl1
  cf2 = ((afroot(i6)*npp_aux) - (cf1/((tfroot(i6))*365)))+ cf1
  if(aawood(i6) .gt. 0.0) then
     ca2 = ((aawood(i6)*npp_aux) - (ca1/((tawood(i6))*365)))+ ca1
  else
     ca2 = 0.0
  endif
  
  
  
  return
end subroutine allocation





subroutine spinup(nppot,cleafini,cfrootini,cawoodini)
!        cbwoodini,cstoini,cotherini,crepini) 

  IMPLICIT NONE
  
  integer, parameter :: nt=5000
  integer, parameter :: npft=3
  !     inputs
  integer i6, kk, k
  
  real, intent(in) :: nppot
  real :: sensitivity,sensitivity2

      
  !c     outputs
  real, intent(out) :: cleafini(npft)
  real, intent(out) :: cawoodini(npft)
  !real, intent(out) :: cbwoodini(npft)
  real, intent(out) :: cfrootini(npft)
  !real :: cstoini(npft)
  !real :: cotherini(npft)
  !real :: crepini(npft)

!!$c     internal vars
!!$
!!$
  real cleafi_aux(nt)
  real cfrooti_aux(nt)
  real cawoodi_aux(nt)
!!$!c      real cbwoodi_aux(nt)
!!$!c      real cstoi_aux(nt)
!!$!c      real cotheri_aux(nt)
!!$!c      real crepi_aux(nt)
!!$
!!$      
!!$    
  real aleaf(3)             !npp percentage alocated to leaf compartment
  data aleaf /0.35, 0.27, 0.45/
  real aawood (3)           !npp percentage alocated to aboveground woody biomass compartment
  data aawood /0.40, 0.40, 0.0/
  real afroot(3)            !npp percentage alocated to fine roots compartment
  data afroot /0.25, 0.33, 0.55/ 
!!$c      real abwood(3)            !npp percentage alocated to belowground woody biomass compartment
!!$c      data abwood /0.10,0.10,0.001/
!!$c      real asto(3)              !npp percentage alocated to storage compartment
!!$c      data asto /0.10,0.10,0.10/
!!$c      real arep(3)              !npp percentage alocated to reproduction compartment
!!$c      data arep /0.15,0.15,0.10/
!!$c      real aother(3)            !npp percentage alocated to other compartment
!!$c      data aother /0.05,0.05,0.06/ 
  
  real tleaf(3)             !turnover time of the leaf compartment (yr)
  data tleaf /1.0, 0.5, 1.0/ 
  real tawood (3)           !turnover time of the aboveground woody biomass compartment (yr)
  data tawood /31.0, 33.0, 0.0/
  real tfroot(3)            !turnover time of the fine roots compartment
  data tfroot /3.0, 2.0, 1.0/
!!$c      real tbwood (3)           !turnover time of the belowground woody biomass compartment
!!$c      data tbwood /40.0,40.0,40.0/
!!$c      real tsto  (3)            !turnover time of the storage compartmentturn
!!$c      data tsto /5.0,5.0,5.0/ 
!!$c      real trep (3)             !turnover time of the reproduction compartment
!!$c      data trep /0.25,0.25,0.25/ 
!!$c      real tother (3)           !turnover time of the other compartment
!!$c      data tother /0.12,0.12,0.12/

  sensitivity = 1.10
  sensitivity2 = 1.40
  
  do i6=1,npft
     do k=1,nt
        if (k.eq.1) then
           cleafi_aux(k) = aleaf(i6)*(nppot)
           cawoodi_aux(k) = aawood(i6)*(nppot)
           cfrooti_aux(k) = afroot(i6)*(nppot)
           
        else
           cleafi_aux(k) = ((aleaf(i6)*(nppot))-(cleafi_aux(k-1)/(tleaf(i6)))) + cleafi_aux(k-1)
           cawoodi_aux(k) = ((aawood(i6)*(nppot))-(cawoodi_aux(k-1)/(tawood(i6)))) + cawoodi_aux(k-1)
           cfrooti_aux(k) = ((afroot(i6)*(nppot))-(cfrooti_aux(k-1)/(tfroot(i6)))) + cfrooti_aux(k-1)

           
           kk =  int(k*0.66)
           if(aawood(i6) .gt. 0.0) then
              if((cfrooti_aux(k)/cfrooti_aux(kk).lt.sensitivity).and.&
                   (cleafi_aux(k)/cleafi_aux(kk).lt.sensitivity).and.&
                   (cawoodi_aux(k)/cawoodi_aux(kk).lt.sensitivity2)) then              
                 
                 cleafini(i6) = cleafi_aux(k) ! carbon content (kg m-2)
                 
                 cawoodini(i6) = cawoodi_aux(k)
                 cfrootini(i6) = cfrooti_aux(k)
                 exit   
              endif
              
           else
              if((cfrooti_aux(k)/cfrooti_aux(kk).lt.sensitivity).and.&
                   (cleafi_aux(k)/cleafi_aux(kk).lt.sensitivity)) then              
                 
                 cleafini(i6) = cleafi_aux(k) ! carbon content (kg m-2)
                 cawoodini(i6) = 0.0
                 cfrootini(i6) = cfrooti_aux(k)

                 exit
              endif
           endif
           
              
        endif
        
     enddo
     
  enddo
  
  
  
  return
  
end subroutine spinup





subroutine budget (pft,month,w1,g1,s1,tsoil,temp,prec,p0,ae,ca,&
     ipar,cl1_pft,ca1_pft,cf1_pft,w2,g2,s2,cl2_pft,ca2_pft,cf2_pft,smavg,&
     ruavg,evavg,epavg,phavg,aravg,nppavg,laiavg,clavg,csavg,hravg,rcavg,&
     rmlavg,rmfavg,rmsavg,rmavg,rglavg,rgfavg,rgsavg,rgavg,cleafavg_pft,&
     cawoodavg_pft,cfrootavg_pft)
  
  implicit none
  ! Surface water (soil moisture, snow and ice) budget for a single month.
  !c
  !c I/O variables
  !c -------------
  !c input  month : actual month (1-12)
  !c        w1    : initial (previous month last day) soil moisture storage (mm)
  !c        g1    : initial soil ice storage (mm)
  !c        s1    : initial overland snow storage (mm)
  !c        tsoil : soil temperature (oC)
  !c        temp  : surface air temperature (oC)
  !c        prec  : precipitation (mm/day)
  !c        p0    : surface pressure (mb)
  !c        ae    : available energy (W/m2)
  !c output w2    : final (last day) soil moisture storage (mm)
  !c        g2    : final soil ice storage (mm)
  !c        s2    : final overland snow storage (mm)
  !c        smavg : snowmelt monthly average (mm/day)
  !c        ruavg : runoff monthly average (mm/day)
  !c        evavg : actual evapotranspiration monthly average (mm/day)
  !c        epavg : maximum evapotranspiration monthly average (mm/day)
  !c
  !c=======================================================================
  !c
  !c i/o variables
  integer, intent( in) ::  month, pft
  real, intent( in) ::  w1,g1,s1,tsoil,temp,prec,p0,ae,ca,ipar
  real, intent( in) ::  cl1_pft,ca1_pft,cf1_pft
  real, intent(out) ::  w2,g2,s2,cl2_pft,ca2_pft,cf2_pft
  real, intent(out) ::  smavg,ruavg,evavg,epavg,rcavg,phavg,aravg,nppavg,laiavg
  real, intent(out) ::  clavg,csavg,hravg
  real, intent(out) ::  rmlavg,rmfavg,rmsavg,rmavg,rglavg,rgfavg,rgsavg
  real, intent(out) ::  rgavg,cleafavg_pft,cawoodavg_pft,cfrootavg_pft
  !c
  !c internal variables
  real :: rh,wmax,tsnow,tice
  real :: psnow,prain
  real :: w,g,s
  real :: rimelt,smelt,roff,evap,emax, ds, dw, rc2
  integer ::  ndmonth(12) !number of days for each month
  data ndmonth /31,28,31,30,31,30,31,31,30,31,30,31/
  !c carbon cycle
  integer :: i 
  real :: ph,ar,nppa,laia,cl,cs,hr,f5
  real gs                   !Stomatal conductance (mol/m2/s)
  real f1                   !Photosynthesis (mol/m2/s)
  real f1b                  !Photosynthesis (micromol/m2/s)
  real rm,rml,rmf,rms,rg,rgl,rgf,rgs ! maintenance and growth plant respiration
  real cl1, cf1, ca1, cl2,cf2, ca2
  real alfa_leaf, alfa_awood, alfa_froot, beta_leaf, beta_awood, beta_froot
  real vpd
  
  !     Initialize Canopy Resistence Parameters
  !     ---------------------------------------
  !     
  gs  = 0.0
  rc2 = 0.0
  f1  = 0.0
  f1b = 0.0
  !     
  !c
  !c parameters
  !c      rh    = 0.6   !relative humidity (adimensional)
  rh    = 0.685 !from NCEP-NCAR Reanalysis data
  wmax  = 500.0 !soil moisture availability (mm)
  tsnow = -1.0  !temperature threshold for snowfall (oC)
  tice  = -2.5  !temperature threshold for soil freezing (oC)
  !c
  !c precipitation [Eq. 3]
  psnow = 0.0
  prain = 0.0
  if (temp.lt.tsnow) then
     psnow = prec/real(ndmonth(month)) !snowfall (mm/day)
  else
     prain = prec/real(ndmonth(month)) !rainfall (mm/day)
  endif
  !c
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
  !c numerical integration
  do i=1,ndmonth(month)
     
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
     
     
     !c
     !c carbon cycle (photosynthesis, plant respiration and NPP)
     call prod (pft, temp, p0, w, wmax, ca, ipar, tsoil, cl1, ca1,&
          cf1, beta_leaf, beta_awood, beta_froot, ph, ar,&
          nppa,laia, f5, f1, vpd, rm, rml, rmf, rms, rg, rgl,&
          rgf, rgs)
     
     !cl1_pft = cl1
     !cf1_pft = cf1
     !ca1_pft = ca1
     !     
     !     carbon allocation (carbon content on each compartment)
     call allocation (pft, nppa, cl1, ca1, cf1, cl2, ca2, cf2)
     
!!$     C         if(nppa .gt. 0.) PRINT*, NPPA, 'nppa_ after alloc'
!!$     c         if(cl2 .gt. 0.) PRINT*, cl2, 'cl2_ after alloc'
!!$     c         if(ca2 .gt. 0.) PRINT*, ca2, 'ca2_ after alloc'
!!$     c         if(cf2 .gt. 0.) PRINT*, cf2, 'cf2_ after alloc'
     
     alfa_leaf  = cl2 - cl1 
     alfa_awood = ca2 - ca1 
     alfa_froot = cf2 - cf1
     
     
     !     Maximum evapotranspiration (emax)
     !     =================================
     if (pft.eq.1) then
        call evpot2 (p0,temp,rh,ae,emax)
     endif
     
     !c snow budget
     smelt = 2.63 + 2.55*temp + 0.0912*temp*prain !snowmelt (mm/day) [Eq. 4]
     smelt = amax1(smelt,0.)
     smelt = amin1(smelt,s+psnow)
     ds = psnow - smelt ![Eq. 2]
     s = s + ds
     !c
     !c water budget
     if (tsoil.le.tice) then !frozen soil
        g = g + w !soil moisture freezes
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
        rc2 = 100.0 !default value, equal to aerodynamic resistance (below)
     else                    !non-frozen soil
        w = w + g !soil ice melts
        g = 0.0
        rimelt = 0.0
        if (w.gt.wmax) then
           rimelt = w - wmax !runoff due to soil ice melting
           w = wmax
        endif
        
        call canopy_resistence (pft,vpd,f1,rc2)     
        call runoff_c (w,wmax,roff) !Soil moisture runoff (roff, mm/day)
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
        call carbon_hr (tsoil,f5,evap,laia,cl,cs,hr)
        !     
     endif
     
     !     Updating monthly values
     !     =======================
     !
     if (pft.eq.1) epavg = epavg + emax   !mm/day
     smavg = smavg + smelt  !mm/day
     ruavg = ruavg + roff   !mm/day
     evavg = evavg + evap   !mm/day
     
     rcavg = rcavg + rc2    !s/m/day
     phavg = phavg + ph /365.0 !kgC/m2/day
     aravg = aravg + ar /365.0 !kgC/m2/day
     nppavg = nppavg + nppa /365.0 !kgC/m2/day
     laiavg = laiavg + laia /365.0 !m2leaf/m2area/day
     clavg = clavg + cl /365.0 !kgC/m2/day
     csavg = csavg + cs /365.0 !kgC/m2/day
     hravg = hravg + hr /365.0 !kgC/m2/day
     rmlavg = rmlavg + rml /365
     rmfavg = rmfavg + rmf /365
     rmsavg = rmsavg + rms /365
     rmavg = rmavg + rm /365
     rglavg = rglavg + rgl /365
     rgfavg = rgfavg + rgf /365
     rgsavg = rgsavg + rgs /365
     rgavg = rgavg + rg /365
     cleafavg_pft = cleafavg_pft + cl2 /365
     cawoodavg_pft = cawoodavg_pft + ca2 /365
     cfrootavg_pft = cfrootavg_pft + cf2 /365
     
     !c         if(nppavg .gt. 0.) PRINT*, NPPAVG, 
     !c     $    'nppa_ after monthly integration'
     
  enddo
  !     
  !     Final calculations
  !     ------------------
  !     
  w2 = w
  g2 = g
  s2 = s
  cl2_pft = cl2
  ca2_pft = ca2
  cf2_pft = cf2
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
  !  c      print*, phavg, 'phavg' 
  
end subroutine budget

!     subroutine wbm (prec,temp,lsmk,p0,ca,par, 
!    &     cleaf_ini, cawood_ini, cfroot_ini,   
!     &     emaxm, tsoil, photo_pft,aresp_pft,npp_pft,lai_pft,
!     &     clit_pft,csoil_pft, hresp_pft,rcm_pft,runom_pft,
!     &     evapm_pft,wsoil_pft,rml_pft,rmf_pft,rms_pft,rm_pft,rgl_pft
!     $     ,rgf_pft,rgs_pft,rg_pft,cleaf_pft,cawood_pft,cfroot_pft)



subroutine wbm (prec,temp,p0,ca,par, cleaf_ini, cawood_ini, cfroot_ini,&
     emaxm, tsoil, photo, aresp, npp, lai, clit, csoil, hresp,&
     rcm, runom, evapm, wsoil, rml, rmf, rms, rm, rgl, rgf, rgs,&
     rg, cleaf, cawood, cfroot) 
  
  
  implicit none
  !c=======================================================================
  !c
  !c Water balance model (WBM). From monthly climatologies of
  !c precipitation and surface temperature, the WBM calculates the
  !c environmental variables.
  !c
  !c 05Jul2005, MDO: ae, rh & runoff are changed.
  !c 11Jul2005, MDO: wsoil2 is written (for testing purpose).
  !c 31Ago2006, DML: carbon cycle is included
  !c
  !c=======================================================================
  !c
  !c i/o variables
  integer, parameter :: m = 12, q = 3
  integer p
  real, dimension(m), intent( in) :: prec,temp,p0,par
  real, dimension(q), intent(in) :: cleaf_ini, cawood_ini, cfroot_ini 
  real, intent( in) :: ca
  
  real, dimension(m,q), intent(out) :: npp,photo,aresp,rcm , wsoil, runom,evapm,lai, clit, csoil, hresp
  real, dimension(m,q), intent(out) :: rml, rmf, rms, rm, rgl, rgf, rgs, rg, cleaf, cawood, cfroot
  real, dimension(m), intent(out) ::  emaxm,tsoil 
  
  !c internal variables
  real, dimension(m) :: spre
  real, dimension(m,q) :: wg0, gsoil, ssoil, snowm
  
  real ae,wini,gini,sini,ipar,wfim,gfim,sfim,smes,rmes,emes,epmes,phmes
  real rmlmes,rmfmes,rmsmes,rmmes, clfim, cafim, crfim, cawood1_pft, cfroot1_pft, cleaf1_pft 
  real rglmes,rgfmes,rgsmes,rgmes,cleafmes,cawoodmes,cfrootmes
  real armes,nppmes,laimes, clmes,csmes,hrmes,rcmes,ca1,dwww,wmax
  real pr,ps,ta,td
  !c
  integer :: n, k, kk, mes, nerro
  
  !c Soil temperature
  ca1 = ca * 0.101325 ! atm CO2 pressure in Pa
  call soil_temp(temp, tsoil)
  
  !c Water budget
  do k=1,m
     emaxm(k) = 0. !average maximum evapotranspiration (mm/day)
     spre(k) = p0(k) * 0.01
     do p = 1,q
        wsoil(k,p) = 0. !soil moisture (mm)
        gsoil(k,p) = 0. !soil ice (mm)
        ssoil(k,p) = 0. !soil snow (mm)
        snowm(k,p) = 0. !average snowmelt (mm/day)
        runom(k,p) = 0. !average runoff (mm/day)
        evapm(k,p) = 0. !average actual evapotranspiration (mm/day)
        wg0(  k,p) = 0. !soil moisture of the previous year (mm)
        rcm(  k,p) = 0. !average canopy resistance (s/m)
        lai(  k,p) = 0.
        photo(k,p) = 0.
        aresp(k,p) = 0.
        npp(  k,p) = 0.
        clit( k,p) = 0.
        csoil(k,p) = 0.
        hresp(k,p) = 0.
        rml(  k,p) = 0.
        rmf(  k,p) = 0.0
        rms(  k,p) = 0.0
        rm(   k,p) = 0.0  
        rgl(  k,p) = 0.0
        rgf(  k,p) = 0.0
        rgs(  k,p) = 0.0
        rg(   k,p) = 0.0
        cleaf(k,p) = 0.0 !monthly npp alloc to leaf biomass (KgC/m2)
        cawood(k,p) = 0.0 !monthly npp alloc to aboveground biomass (KgC/m2)
        cfroot(k,p) = 0.0
     enddo
  enddo
  
  do p=1,q
     
     wini  = 0.01 !soil moisture initial condition (mm)
     gini  = 0.0  !soil ice initial condition (mm)
     sini  = 0.0  !overland snow initial condition (mm)
     
     cleaf1_pft  =  cleaf_ini(p) !inital leaf biomass for each PFT from spinup (KgC/m2) 
     cawood1_pft = cawood_ini(p)
     cfroot1_pft = cfroot_ini(p)
     !c
     !c initialization
     do k=1,m
        wg0(k,p) = -1.0
     enddo
     !surface pressure (mb)
     !c
     !c start integration
     n = 0
10   continue
     n = n + 1
     !c
     !c pre-processing
     k = mod(n,12)
     if (k.eq.0) k = 12
     mes = k
     ps = spre(k)
     td = tsoil(k)
     ta = temp(k)
     pr = prec(k)
     ipar = par(k) / 2.18e5   !converting to Ein/m2/s
     !c      ae = 2.26457*ta + 67.5876 !available energy (W/m2) [Eq. 8]
     ae = 2.895*ta + 52.326 !from NCEP-NCAR Reanalysis data
     !c
     smes = 0.0
     rmes = 0.0
     emes = 0.0
     epmes = 0.0
     phmes = 0.0
     armes = 0.0
     nppmes = 0.0
     laimes = 0.0
     clmes = 0.0
     csmes = 0.0
     hrmes = 0.0
     rcmes = 0.0
     rmlmes = 0.0
     rmfmes = 0.0
     rmsmes = 0.0
     rmmes = 0.0
     rglmes = 0.0
     rgfmes = 0.0
     rgsmes = 0.0
     rgmes = 0.0
     cleafmes = 0.0 
     cawoodmes  = 0.0
     cfrootmes = 0.0
     !c monthly water budget
     
     !subroutine budget (pft,month,w1,g1,s1,tsoil,temp,prec,p0,ae,ca,&
     !     ipar,cl1_pft,ca1_pft,cf1_pft,w2,g2,s2,cl2_pft,ca2_pft,cf2_pft,smavg,&
     !     ruavg,evavg,epavg,phavg,aravg,nppavg,laiavg,clavg,csavg,hravg,rcavg,&
     !     rmlavg,rmfavg,rmsavg,rmavg,rglavg,rgfavg,rgsavg,rgavg,cleafavg_pft,&
     !     cawoodavg_pft,cfrootavg_pft
     
     call budget (p,mes,wini,gini,sini,td,ta,pr,ps,ae,ca1,ipar,&
          cleaf1_pft, cawood1_pft, cfroot1_pft,wfim,gfim,&
          sfim, clfim, cafim,crfim,smes,rmes,emes,epmes,phmes,&
          armes,nppmes,laimes,clmes,csmes,hrmes,rcmes,rmlmes,&
          rmfmes, rmsmes, rmmes, rglmes, rgfmes,rgsmes,rgmes,&
          cleafmes, cawoodmes,cfrootmes)
     
     !c update variables
     if(p .eq. 1) emaxm(k) = epmes
     wsoil(k,p) = wfim
     gsoil(k,p) = gfim
     ssoil(k,p) = sfim
     snowm(k,p) = smes
     runom(k,p) = rmes
     evapm(k,p) = emes
     
     rcm(k,p) = rcmes
     lai(k,p) = laimes
     photo(k,p) = phmes
     aresp(k,p) = armes
     npp(k,p) = nppmes
     clit(k,p) = clmes
     csoil(k,p) = csmes
     hresp(k,p) = hrmes
     rml(k,p) =rmlmes
     rmf (k,p) =rmfmes
     rms(k,p) = rmsmes
     rm(k,p) =  rmmes
     rgl(k,p) =  rglmes
     rgf(k,p) =  rgfmes
     rgs(k,p) =  rgsmes
     rg(k,p) =  rgmes 
     cleaf(k,p) =  cleafmes  
     cawood(k,p) =  cawoodmes
     cfroot(k,p) =  cfrootmes
     
     wini = wfim
     gini = gfim
     sini = sfim
     cleaf1_pft  =  clfim 
     cawood1_pft = cafim
     cfroot1_pft = crfim
     !c
     !c check if equilibrium is attained (k=12)
     if (k.eq.12) then
        wmax = 500.
        nerro = 0
        do kk=1,12
           dwww = (wsoil(kk,p)+gsoil(kk,p)-wg0(kk,p))/wmax
           if (abs(dwww).gt.0.001) nerro = nerro + 1
        enddo
        if (nerro.ne.0) then
           do kk=1,12
              wg0(kk,p) = wsoil(kk,p) + gsoil(kk,p)
           enddo
        else
           goto 100
        endif
     endif
     
     goto 10
100  continue
  enddo
  return
end subroutine wbm
