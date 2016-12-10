! this code is based on CPTEC-PVM2 source code. Wrote by Carlos Nobre, Marcos Oyama e David Lapola.
! author :: jpdarela 
! in an attempt to write CAETE calculation routines


subroutine prod (temp,p0,w,wmax,ca,ipar,ph,ar,nppa,laia,f5)
  implicit none

  
  !c i/o variables
  real,intent(in) :: temp  !mean monthly temperature (oC)
  real,intent(in) :: p0  !mean surface pressure (hPa)
  real,intent(in) :: w,wmax !soil moisture (dimensionless)
  real,intent(in) :: ca  !Atmospheric CO2 pressure (Pa)
  real,intent(in) :: ipar  !incident photos. active radiation

  
  real, intent(out) :: ph  !Canopy gross photosyntheis(kgC/m2/yr)
  real, intent(out) :: laia !Leaf area index (m2 leaf/ m2 area)
  real, intent(out) :: ar !Autotrophic (plant) respiration (kgC/m2/yr)
  real, intent(out) :: nppa !Net primary productivity (kgC/m2/yr)
  real, intent(out) :: f5
  
  !c internal variables

  real wa
  real vm,mgama,es,rmax,r,ci,a,b,c,a2,b2,c2,delta,delta2,f1,f1a,f2,f3,f4,k16,k17,rl,rp
  real f10 !CO2 sensitivity
  real sunlai, shadelai, f4shade, f4sun
  
  !c Rubisco, light and transport limited photosynthesis rate respectively
  real jc,jl,je,jp,jp1,jp2,j1,j2 !auxiliars
  !c
  !c
  !c
  !c=======================================================================
  !c Photosynthesis =======================================================
  !c
  !c Rubisco maximum carboxylaton rate
  !c (vm ; molCO2/m^2s) [Eq. 12]
  
  vm = (0.00004*2.0**(0.1*(temp-25.0)))/(1.0+exp(0.3*(temp-36.0)))
  
  ! new CO2 uptake leaf level model
  ! old:
  
  !Jc
  
  ! Jc[vm(T,Ca)] <- Vm(T) * (ci(T,Ca) - mgama(T))/
  !                     (ci(T,Ca) +  f2(T) * (1 +(k3/f3(T)))
  
  ! N Limitation goes here
  ! new variables  Vm(Vmax)
  
  ! Vm(T) = is * vs * Na 
  ! 
  ! vm = 
  
  
  
  !c Photo-respiration compensation point
  !c (mgama ; Pa) [Eq. 8]
  !c
  mgama = 21200.0/(5200.0*(0.57**(0.1*(temp-25.0))))
  !c
  !c Michaelis-Menten CO2 constant
  !c (f2 ; Pa) [Eq. 9]
  !c
  f2 = 30*(2.1**(0.1*(temp-25.0)))
  !c
  !c Michaelis-Menten O2 constant
  !c (f3 ; Pa) [Eq. 10]
  !c
  f3 = 30000.0*(1.2**(0.1*(temp-25.0)))
  !c
  !c Saturation partial pressure of water vapour at temperature 'temp' 
  !c (es ; Pa) [from WBM, to be used in Eq. 14]
  !c subroutine tetens (t,es)
  
  !call tetens (temp,es)
  if (temp.ge.0.) then
     es = 6.1078*exp((7.5*temp/(237.3+temp))*log(10.))
  else
     es = 6.1078*exp((9.5*temp/(265.5+temp))*log(10.))
  endif
  !c
  !c Saturated mixing ratio
  !c (rmax ; kg/kg) [Eq. 14]
  rmax = 0.622*(es/(p0-es))
  
  !c Moisture deficit at leaf level
  !c (r ; kg/kg) [Eq. 13]
  !c Considering that relative humidity is constant (0.685)
  !c
  r = (-0.315)*rmax
  !c
  !c Internal leaf CO2 partial pressure
  !c (ci ; Pa) [Eq. 11]
  !c
  !c      k16 = 0.875 !Generic value for forest types
  !c      k17 = 0.08 !Generic value for forest types
  k16 = 0.9 !Value for C3 grass types
  k17 = 0.1 !Value for C3 grass types
  !c
  ci = k16*(1-(r/k17))*(ca-mgama)+mgama
  !c
  !c Empirical function for reduction of NPP sensitivity to CO2
  !c based in simulations for pre-Industrial period (dimensionless)
  f10 = 1.23/(1+exp(-1*(ca-24.11)/8.93))
  !c
  !c Rubisco carboxilation limited photosynthesis rate
  !c (jc ; molCO2/m^2s) [Eq. 4]
  !c
  jc = vm*((ci-mgama)/(ci+(f2*f10*(1+(21200.0/f3)))))
  !c
  !c Light limited photosynthesis rate
  !c (jl ; molCO2/m^2s) [Eq. 5]
  !c
  jl = 0.08*(1.0-0.15)*ipar*((ci-mgama)/(ci+(2.0*mgama)))
  !c
  !c Transport limited photosynthesis rate
  !c (je ; molCO2/m^2s) [Eq. 6]
  !c
  je = 0.5*vm
  !c
  !c jp (minimum between jc and jl)
  !c [Eq. 3]
  !c
  a = 0.83
  b = (-1)*(jc+jl)
  c = jc*jl
  delta = (b**2)-4.0*a*c
  !c
  jp1=(-b-(sqrt(delta)))/(2.0*a)
  jp2=(-b+(sqrt(delta)))/(2.0*a)
  jp= amin1(jp1,jp2)
  !c
  !c Leaf level gross photosynthesis (minimum between jc, jl and je)
  !c [f1, Eq. 2]
  !c
  a2 = 0.93
  b2 = (-1)*(jp+je)
  c2 = jp*je
  delta2 = (b2**2)-4.0*a2*c2
  !c
  j1=(-b2-(sqrt(delta2)))/(2.0*a2)
  j2=(-b2+(sqrt(delta2)))/(2.0*a2)
  f1a = amin1(j1,j2)
  !c
  !c Soil water
  wa = w/wmax
  !c
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
     f1 = f1a*f5
  else
     f1 = 0.0 !Temperature above/below photosynthesis windown 
  endif
  !c
  !c Leaf area index (m2 leaf/m2 area)
  !c [lai ; Eq. 15]
  !c
  !c      laia  = 0.2*exp(2.5*(f1/0.000008))
  laia  = 0.25*exp(2.5*(f1/0.000008)) !adjusted after using observed ipar
  !c
  !c LAI with direct incident sun
  !c [sunlai; Eq. 16]      
  sunlai = (1.0-(exp(-0.5*laia)))/0.5
  !c
  !c LAI below sunlai
  !c [shadelai, Eq. 17]
  shadelai = laia - sunlai
  !c
  !c Scaling-up to canopy level (dimensionless)
  !c [f4 ; Eq. 18]
  !c
  f4 = (1.0-(exp(-0.5*laia)))/0.5 !sun 90 in the whole canopy, to be used for respiration
  !c
  !c Sun/Shade approach to canopy scaling (based in de Pury & Farquhar 1997) 
  !c [f4sun, f4shade ; Eqs. 19, 20]
  !c
  f4sun = (1.0-(exp(-0.5*sunlai)))/0.5 !sun 90
  f4shade = (1.0-(exp(-1.5*shadelai)))/1.5 !sun ~20
  !c
  !c Canopy gross photosynthesis (kgC/m^2/yr)
  !c [ph ; Eq. 1]
  !c (0.012 converts molCO2 to kgC)
  !c [31557600 converts seconds to year (with 365.25 days)]
  !c
  ph = 0.012*31557600.0*f1*f4sun*f4shade
  !!c
  !c=======================================================================
  !c Plant (autotrophic) respiration ======================================
  !c     
  
  !c Leaf respiration (kgC/m2/yr)
  !c [rl ; Eq. 23]
  !c
  rl = 0.012*31557600.0*0.015*vm*f4*f5
  
  !c
  !c Non-leaf parts respiration (kgC/m2/yr)
  !c [rp ; Eq. 24]
  !c
  rp = 3.85*rl !3.85?
  !c
  !c Autotrophic (plant) respiration (kgC/m2/yr)
  !c [ar ; Eq. 22]
  !c
  !c Respiration minimum and maximum temperature
  if ((temp.ge.-10.0).and.(temp.le.50.0)) then
     ar = rl+rp
  else
     ar = 0.0 !Temperature above/below respiration windown 
  endif
  
  !c Productivity =========================================================
  !c
  !c Net primary productivity (kgC/m2/yr)
  !c [npp ; Eq. 25]
  !c
  nppa = ph-ar
  !if (nppa.lt.0.0) nppa = 0.0
  
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
  !c internal variables
  real lf,f6,f7
  !c
  !c========================================================================
  !c Litter decayment function 
  !c (controlled by annual evapotranspiration)
  !c [f6 ; Eq. 29]
  !c
  f6 = 1.16*10**(-1.4553+0.0014175*(evap*365.0))
  !c
  !c Soil carbon storage function
  !c (controlled by temperature)
  !c [f7 ; Eq.31]
  !c
  f7 = 2.0**(0.1*(tsoil-25.0))
  !c Litterfall (kgC/m2)
  !c [lf ; Eq. 27]
  !c
  lf = 0.1*laia
  !c
  !c Litter carbon (kgC/m2)
  !c [cl ; Eq. 28]
  !c
  cl = lf/f6
  !c Soil carbon (kgC/m2)
  !c [cs ; Eq. 30]
  !c
  cs = ((0.3*cl)/(0.05*f7))*f5
  !c
  !c Heterotrophic (Microbial) respiration (kgC/m2/yr)
  !c [hr ; Eq. 26]
  !c
  !c Respiration minimum and maximum temperature
  if ((tsoil.ge.-10.0).and.(tsoil.le.50.0)) then
     hr = 0.25*(cl*(f6**2)+(cs*f5*evap*(f7**2))) !Litter and Soil respectively
  else
     hr = 0.0 !Temperature above/below respiration windown 
  endif
  !c
end subroutine carbon_hr






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






subroutine canopy_resistence(npp, ca, p0, rc2)

  real, intent( in) :: npp, ca, p0
  real, intent(out) :: rc2

  real :: nppb, p1

  p1 = p0 * 100. ! convertendo mbar(hPa) para Pa
  
  nppb = amax1(npp,0.05)
  rc2 = (ca/(0.9*(nppb*2.64e-6)*0.685*p1))

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







subroutine budget (month,w1,g1,s1,tsoil,temp,prec,p0,ae,ca,&
     ipar,w2,g2,s2,smavg,ruavg,evavg,&
     epavg,phavg,aravg,nppavg,laiavg,& 
     clavg,csavg,hravg,rcavg)

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
  integer, intent( in) ::  month
  real, intent( in) ::  w1,g1,s1,tsoil,temp,prec,p0,ae,ca,ipar
  real, intent(out) ::  w2,g2,s2,smavg,ruavg,evavg,epavg,rcavg
  real, intent(out) ::  phavg,aravg,nppavg,laiavg
  real, intent(out) ::  clavg,csavg,hravg
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
  real :: ph,ar,nppa,laia,cl,cs,hr, f5
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
  !c initialization
  w = w1 !w = daily soil moisture storage (mm)
  g = g1 !g = daily soil ice storage (mm)
  s = s1 !s = daily overland snow storage (mm)
  smavg = 0.
  ruavg = 0.
  evavg = 0.
  epavg = 0.
  rcavg = 0.
  laiavg = 0.
  phavg = 0.
  aravg = 0.
  nppavg = 0.
  clavg = 0.
  csavg = 0.
  hravg = 0.
  !c
  !c numerical integration
  do i=1,ndmonth(month)
     !c
     !c carbon cycle (photosynthesis, plant respiration and NPP)
     call prod (temp,p0,w,wmax,ca,ipar, ph,ar,nppa,laia,f5)
     !c maximum evapotranspiration (emax)
     call evpot2 (p0,temp,rh,ae,emax)
     !c
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
        
        
        call canopy_resistence(nppa, ca, p0, rc2)
        call runoff_c (w,wmax,roff) !soil moisture runoff (roff, mm/day) [Eq. 10]
        call penman (p0,temp,rh,ae,rc2,evap) !actual evapotranspiration (evap, mm/day)
        dw = prain + smelt - evap - roff ![Eq. 1]
        w = w + dw
        if (w.gt.wmax) then
           roff = roff + (w - wmax)
           w = wmax
        endif
        if (w.lt.0.) w = 0.
        roff = roff + rimelt !total runoff
        !    c carbon cycle (Microbial respiration, litter and soil carbon)
        call carbon_hr (tsoil,f5,evap,laia, cl,cs,hr)              !output
     endif
     !  c
     !  c updating monthly values
     smavg = smavg + smelt
     ruavg = ruavg + roff
     evavg = evavg + evap
     epavg = epavg + emax
     rcavg = rcavg + rc2
     phavg = phavg + ph/365.0 !kgC/m2
     aravg = aravg + ar/365.0 !kgC/m2
     nppavg = nppavg + nppa/365.0 !kgC/m2
     laiavg = laiavg + laia/365.0
     clavg = clavg + cl/365.0
     csavg = csavg + cs/365.0
     hravg = hravg + hr/365.0 !kgC/m2
     !   c
  enddo
  ! c
  !c final calculations
  w2 = w
  g2 = g
  s2 = s
  smavg = smavg/real(ndmonth(month))
  ruavg = ruavg/real(ndmonth(month))
  evavg = evavg/real(ndmonth(month))
  epavg = epavg/real(ndmonth(month))
  rcavg = rcavg/real(ndmonth(month))
  phavg = phavg*12.0 !kgC/m2/yr
  aravg = aravg*12.0 !kgC/m2/yr
  nppavg = nppavg*12.0 !kgC/m2/yr
  laiavg = laiavg*12.0
  clavg = clavg*12.0 !kgC/m2
  csavg = csavg*12.0 !kgC/m2
  hravg = hravg*12.0 !kgC/m2/yr
end subroutine budget





subroutine wbm (prec,temp,p0,ca,par,&
     npp,photo,aresp,rcm,tsoil, wsoil, runom, evapm, emaxm, lai, clit, csoil, hresp)
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
  integer, parameter :: m = 12
  real, dimension(m), intent( in) :: prec,temp,p0,par
  real, intent( in) :: ca
  
  real, dimension(m), intent(out) :: npp,photo,aresp,rcm,tsoil, wsoil, runom
  real, dimension(m), intent(out) :: evapm, emaxm, lai, clit, csoil, hresp
  
  !c internal variables
  
  real, dimension(m) :: wg0, spre, gsoil, ssoil, snowm
  
  real ae,wini,gini,sini,ipar,wfim,gfim,sfim,smes,rmes,emes,epmes,phmes
  real armes,nppmes,laimes, clmes,csmes,hrmes,rcmes,ca1,dwww,wmax
  real pr,ps,ta,td
  !c
  integer :: n, k, kk, mes, nerro

  !c Soil temperature
  ca1 = ca * 0.101325 ! atm CO2 pressure in Pa
  call soil_temp(temp, tsoil)
 
  !c Water budget
  do k=1,m
     wsoil(k) = 0. !soil moisture (mm)
     gsoil(k) = 0. !soil ice (mm)
     ssoil(k) = 0. !soil snow (mm)
     snowm(k) = 0. !average snowmelt (mm/day)
     runom(k) = 0. !average runoff (mm/day)
     evapm(k) = 0. !average actual evapotranspiration (mm/day)
     emaxm(k) = 0. !average maximum evapotranspiration (mm/day)
     wg0(  k) = 0. !soil moisture of the previous year (mm)
     rcm(  k) = 0. !average canopy resistance (s/m)
     lai(  k) = 0.
     photo(k) = 0.
     aresp(k) = 0.
     npp(  k) = 0.
     clit( k) = 0.
     csoil(k) = 0.
     hresp(k) = 0.
  end do
  
  wini  = 0.01 !soil moisture initial condition (mm)
  gini  = 0.0  !soil ice initial condition (mm)
  sini  = 0.0  !overland snow initial condition (mm)
  !c
  !c initialization
  do k=1,m
     wg0(k) = -1.0
     spre(k) = p0(k) * 0.01
  enddo
  !surface pressure (mb)
  !c
  !c start integration
  n = 0
10 continue
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
  !c monthly water budget
  
  ! budget args:
  !  month,w1,g1,s1,tsoil,temp,prec,p0,ae,ca,&
  !     ipar,w2,g2,s2,smavg,ruavg,evavg,&
  !     epavg,phavg,aravg,nppavg,laiavg,& 
  !     clavg,csavg,hravg,rcavg)
  
  call budget (mes,wini,gini,sini,td,ta,pr,ps,ae,ca1,ipar,&
       wfim,gfim,sfim,smes,rmes,emes,epmes,&
       phmes,armes,nppmes,laimes,&
       clmes,csmes,hrmes,rcmes)
  !c
  !c update variables
  wsoil(k) = wfim
  gsoil(k) = gfim
  ssoil(k) = sfim
  snowm(k) = smes
  runom(k) = rmes
  evapm(k) = emes
  emaxm(k) = epmes
  rcm(k) = rcmes
  lai(k) = laimes
  photo(k) = phmes
  aresp(k) = armes
  npp(k) = nppmes
  clit(k) = clmes
  csoil(k) = csmes
  hresp(k) = hrmes
  wini = wfim
  gini = gfim
  sini = sfim
  !c
  !c check if equilibrium is attained (k=12)
  if (k.eq.12) then
     wmax = 500.
     nerro = 0
     do kk=1,12
        dwww = (wsoil(kk)+gsoil(kk)-wg0(kk))/wmax
        if (abs(dwww).gt.0.001) nerro = nerro + 1
     enddo
     if (nerro.ne.0) then
        do kk=1,12
           wg0(kk) = wsoil(kk) + gsoil(kk)
        enddo
     else
        goto 100
     endif
  endif
  
  goto 10
100 continue
  return
end subroutine wbm
