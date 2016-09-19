! c234567
! c=======================================================================
! c Subroutine carbon1 calculates photosynthesis, plant respiration and
! c net primary productivity. Subroutine carbon2 calculates microbial
! c respiration (yet to be improved).
! c
! c Code written by David Lapola
! c Last update: Aug/2007
! c
! c A "bug" has been found on Sep/2007: when the water soil
! c is between 0.205 and 0.265, NPP unrealistically drops to a level
! c below those when wsoil is lesser than 0.205 (because of f5).
! c
! c=======================================================================
subroutine carbon1 (temp,p0,w,wmax,ca,ipar,&    !input
                    ph,ar,nppa,laia,f5)        !output
    implicit none
    real, intent(in) :: temp, p0, w, wmax, ca, ipar
    ! temp        !mean monthly temperature (oC)
    ! p0        !mean surface pressure (hPa)
    ! wa,w,wmax    !soil moisture (dimensionless)
    ! ca        !Atmospheric CO2 pressure (Pa)
    ! ipar        !incident photos. active radiation

    real, intent(out) :: ph, ar, nppa, laia, f5
    ! ph        !Canopy gross photosyntheis (kgC/m2/yr)tivity
    ! laia        !Leaf area index (m2 leaf/ m2 area)
    ! ar        !Autotrophic (plant) respiration (kgC/m2/yr)
    ! nppa        !Net primary productivity (kgC/m2/yr)
! c internal variables
    real, parameter :: temp_min = -10, temp_max = 50
    real vm, mgama, es, rmax, r, ci, a, b, c, a2, b2, c2,&
         delta, delta2, f1, f1a, f2, f3, f4, k16, k17,&
         rl, rp, f10, wa, sunlai,shadelai, f4sun, f4shade

         !f10 = !CO2 sensitivity

! c Rubisco, light and transport limited photosynthesis rate respectively
    real jc, jl, je, jp, jp1, jp2, j1, j2    !auxiliars
!=========================================================================
! c Photosynthesis =======================================================
! ========================================================================
! c Rubisco maximum carboxylaton rate
! c (vm ; molCO2/m^2s) [Eq. 12]

    vm = (0.00004*2.0**(0.1*(temp-25.0))) / (1.0+exp(0.3*(temp-36.0)))
    !============================================================
!c new CO2 uptake leaf level model

      ! old:

      !Jc

      ! Jc[vm(T,Ca)] <- Vm(T) * (ci(T,Ca) - mgama(T))/
      !                     (ci(T,Ca) +  f2(T) * (1 +(k3/f3(T)))

      ! N Limitation goes here
      ! new variables  Vm(Vmax)

      ! Vm(T) = is * vs * Na
      !
      ! vm =
      !=========================================================

! c Photo-respiration compensation point
! c (mgama ; Pa) [Eq. 8]
! c
     mgama = 21200.0/(5200.0*(0.57**(0.1*(temp-25.0))))
! c
! c Michaelis-Menton CO2 constant
! c (f2 ; Pa) [Eq. 9]
! c
    f2 = 30*(2.1**(0.1*(temp-25.0)))
! c
! c Michaelis-Menton O2 constant
! c (f3 ; Pa) [Eq. 10]
! c
    f3 = 30000.0*(1.2**(0.1*(temp-25.0)))
! c
! c Saturation partial pressure of water vapour at temperature 'temp'
! c (es ; Pa) [from WBM, to be used in Eq. 14]
! c
    call tetens (temp,es)
! c
! c Saturated mixing ratio
! c (rmax ; kg/kg) [Eq. 14]
! c
    rmax = 0.622*(es/(p0-es))
! c
! c Moisture deficit at leaf level
! c (r ; kg/kg) [Eq. 13]
! c Considering that relative humidity is constant (0.685)
! c
    r = (-0.315)*rmax
! c
! c Internal leaf CO2 partial pressure
! c (ci ; Pa) [Eq. 11]
! ! c
! c      k16 = 0.875   !Generic value for forest types
! c      k17 = 0.08    !Generic value for forest types
    k16 = 0.9    !Value for C3 grass types
    k17 = 0.1    !Value for C3 grass types

    ci = k16*(1-(r/k17))*(ca-mgama)+mgama
! c
! c Empirical function for reduction of NPP sensitivity to CO2
! c based in simulations for pre-Industrial period (dimensionless)
    f10 = 1.23/(1+exp(-1*(ca-24.11)/8.93))
! c
! c Rubisco carboxilation limited photosynthesis rate
! c (jc ; molCO2/m^2s) [Eq. 4]
! c
    jc = vm*((ci-mgama)/(ci+(f2*f10*(1+(21200.0/f3)))))
! c
! c Light limited photosynthesis rate
! c (jl ; molCO2/m^2s) [Eq. 5]
! c
    jl = 0.08*(1.0-0.15)*ipar*((ci-mgama)/(ci+(2.0*mgama)))
! c
! c Transport limited photosynthesis rate
! c (je ; molCO2/m^2s) [Eq. 6]
! c
    je = 0.5*vm
! c
! c jp (minimum between jc and jl)
! c [Eq. 3]
! c
    a = 0.83
    b = (-1)*(jc+jl)
    c = jc*jl
    delta = (b**2)-4.0*a*c

    jp1=(-b-(sqrt(delta)))/(2.0*a)
    jp2=(-b+(sqrt(delta)))/(2.0*a)
    jp= amin1(jp1,jp2)
! c
! c Leaf level gross photosynthesis (minimum between jc, jl and je)
! c [f1, Eq. 2]
! c
    a2 = 0.93
    b2 = (-1)*(jp+je)
    c2 = jp*je
    delta2 = (b2**2)-4.0*a2*c2

    j1=(-b2-(sqrt(delta2)))/(2.0*a2)
    j2=(-b2+(sqrt(delta2)))/(2.0*a2)
    f1a = amin1(j1,j2)
! c
! c Soil water
       wa = w/wmax
! c
! c Water stress response modifier (dimensionless)
! c [f5 ; Eq. 21]
! c
    if (wa > 0.5) f5 = 1.0        !Not too lower in e.g. Amazonian dry season
    if ((wa >= 0.205).and.(wa <= 0.5)) f5 = (wa-0.205)/(0.5-0.205)
    if (wa < 0.205) f5 = wa        !Below wilting point f5 accompains wa (then Sahara is well represented)
! c
! c ...f1
! c
! c Photosysthesis minimum and maximum temperature
    if ((temp >= temp_min).and.(temp <= temp_max)) then
        f1 = f1a*f5
    else
        f1 = 0.0    !Temperature above/below photosynthesis windown
    endif
! c
! c Leaf area index (m2 leaf/m2 area)
! c [lai ; Eq. 15]
! c
! c      laia  = 0.2*exp(2.5*(f1/0.000008))
    laia  = 0.25*exp(2.5*(f1/0.000008)) !adjusted after using observed ipar
! c
! c LAI with direct incident sun
! c [sunlai; Eq. 16]
    sunlai = (1.0-(exp(-0.5*laia)))/0.5
! c
! c LAI below sunlai
! c [shadelai, Eq. 17]
    shadelai = laia - sunlai
! c
! c Scaling-up to canopy level (dimensionless)
! c [f4 ; Eq. 18]
! c
    f4 = (1.0-(exp(-0.5*laia)))/0.5    !sun 90� in the whole canopy, to be used for respiration
! c
! c Sun/Shade approach to canopy scaling (based in de Pury & Farquhar 1997)
! c [f4sun, f4shade ; Eqs. 19, 20]
! c
    f4sun = (1.0-(exp(-0.5*sunlai)))/0.5    !sun 90�
    f4shade = (1.0-(exp(-1.5*shadelai)))/1.5    !sun ~20�
! c
! c Canopy gross photosynthesis (kgC/m^2/yr)
! c [ph ; Eq. 1]
! c (0.012 converts molCO2 to kgC)
! c [31557600 converts seconds to year (with 365.25 days)]
! c
    ph = 0.012*31557600.0*f1*f4sun*f4shade
! c=======================================================================
! c Plant (autotrophic) respiration ======================================
! c
! c Leaf respiration (kgC/m2/yr)
! c [rl ; Eq. 23]
! c
    rl = 0.012*31557600.0*0.015*vm*f4*f5
! c
! c Non-leaf parts respiration (kgC/m2/yr)
! c [rp ; Eq. 24]
! c
    rp = 3.85*rl
! c
! c Autotrophic (plant) respiration (kgC/m2/yr)
! c [ar ; Eq. 22]
! c
! c Respiration minimum and maximum temperature
    if ((temp >= temp_min).and.(temp <= temp_max)) then
        ar = rl+rp
    else
         ar = 0.0    !Temperature above/below respiration windown
    endif

! c=======================================================================
! c Productivity =========================================================
! c
! c Net primary productivity (kgC/m2/yr)
! c [npp ; Eq. 25]
! c
    nppa = ph-ar
    if (nppa.lt.0.0) nppa = 0.0 !maybe this is incorrect, but demands retuning of every biome limits
    return
end subroutine carbon1
!
! c=======================================================================
! c=======================================================================
! c Microbial (heterotrophic) respiration ================================
subroutine carbon2 (tsoil,f5,evap,laia,cl,cs,hr)

    implicit none

    real, intent(in) :: tsoil,& !mean monthly soil temperature (oC)
                        f5,&    !stress response to soil moisture (dimensionless)
                        evap,&  !actual evapotranspiration (mm/day)
                        laia    !Leaf area index (m2 leaf/ m2 area)

    real, intent(out) :: cl, cs, hr
      !cl         Litter carbon (kgC/m2)
      !cs         Soil carbon (kgC/m2)
      !hr         Heterotrophic (microbial) respiration (kgC/m2/yr)
! c internal variables
    real lf,f6,f7
! c
! c========================================================================
! c Litter decayment function
! c (controlled by annual evapotranspiration)
! c [f6 ; Eq. 29]
! c
    f6 = 1.16*10**(-1.4553+0.0014175*(evap*365.0))
! c
! c Soil carbon storage function
! c (controlled by temperature)
! c [f7 ; Eq.31]
! c
    f7 = 2.0**(0.1*(tsoil-25.0))
! c
! c Litterfall (kgC/m2)
! c [lf ; Eq. 27]
! c
    lf = 0.1*laia
! c
! c Litter carbon (kgC/m2)
! c [cl ; Eq. 28]
! c
    cl = lf/f6
! c
! c Soil carbon (kgC/m2)
! c [cs ; Eq. 30]
! c
    cs = ((0.3*cl)/(0.05*f7))*f5
! c
! c Heterotrophic (Microbial) respiration (kgC/m2/yr)
! c [hr ; Eq. 26]
! c
! c Respiration minimum and maximum temperature
    if ((tsoil > -10).and.(tsoil<=50)) then
        hr = 0.25*(cl*(f6**2)+(cs*f5*evap*(f7**2))) !Litter and Soil respectively
    else
        hr = 0.0    !Temperature above/below respiration windown
    endif
! c
    return
end subroutine carbon2
