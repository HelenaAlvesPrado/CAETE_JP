module water
  
  ! this module defines functions related to surface water balance
  implicit none
  private

  ! functions defined here:
  
  public ::              &
       soil_temp        ,&
       penman           ,&
       evpot2           ,&
       available_energy ,&
       runoff

  
contains

  !=================================================================
  !=================================================================
  
  function soil_temp(t0,temp) result(tsoil)
    use global_pars, only: i4, r4, H, TAU, DIFFU
    implicit none
    
    real(kind=r4),intent( in) :: temp
    real(kind=r4),intent( in) :: t0! future __ make temps an allocatable array
    real(kind=r4) :: tsoil
 
    real(kind=r4) :: t1 = 0.0
    
    t1 = (t0*exp(-1.0/TAU) + (1.0 - exp(-1.0/TAU)))*temp
    tsoil = (t0 + t1)/2.0
  end function soil_temp
  
  !=================================================================
  !=================================================================
!!$  
!!$  function soil_water(w_state,temp,prec,p0,ipar,rh)
!!$    ! resturns wsoil of a day (mm) ----> kg(H2O) m-2(ground)
!!$    use global_pars
!!$    implicit none
!!$    
!!$    integer(kind=i4),intent(in) :: id                ! >5 = w; 1 = evap; 2 = emax; 3 = runoff; 4 = snow; 5 = ice.
!!$    real(kind=r4),intent(in) :: w0                   !previous day soil water (mm)
!!$    real(kind=r4),intent(in) :: t0                   !previous day temperature
!!$    real(kind=r4),intent(in) :: temp                 !Surface air temperature (oC)
!!$    real(kind=r4),intent(in) :: prec                 !Precipitation (mm/day)
!!$    real(kind=r4),intent(in) :: p0                   !Surface pressure (mb)
!!$    real(kind=r4),intent(in) :: ipar                 !Incident photosynthetic active radiation
!!$    !__future__ real(kind=r4),intent(in) :: rlds                 !Incident longwave radiation (W/m2)
!!$    real(kind=r4),intent(in) :: rh                   !Relative humidity
!!$    
!!$    real(kind=r4) :: w_state                
!!$    !     -----------------------Internal Variables------------------------
!!$    real(kind=r4) :: tsnow                !Temperature threshold for snowfall (oC)
!!$    real(kind=r4) :: tice                 !Temperature threshold for soil freezing (oC)
!!$    real(kind=r4) :: psnow = 0.0                !Snowfall (mm/day)
!!$    real(kind=r4) :: prain = 0.0          !Rainfall (mm/day)
!!$    real(kind=r4) :: rimelt = 0.0         !Runoff due to soil ice melting
!!$    real(kind=r4) :: smelt = 0.0                !Snowmelt (mm/day)
!!$    real(kind=r4) :: g = 0.0                    !Daily soil ice storage (mm)
!!$    real(kind=r4) :: s = 0.0                    !Daily overland snow storage (mm)
!!$    real(kind=r4) :: ds = 0.0
!!$    real(kind=r4) :: dw = 0.0
!!$    real(kind=r4) :: roff = 0.0                 !Total runoff
!!$    real(kind=r4) :: evap = 0.0                 !Actual evapotranspiration (mm/day)
!!$    real(kind=r8) :: wapft = 0.0         
!!$    real(kind=r4) :: rc2 = 0.0
!!$    
!!$    tsnow = -1.0
!!$    tice  = -2.5
!!$    psnow = 0.0
!!$    prain = 0.0
!!$
!!$    if(t0 .lt. tsnow) s = w0
!!$    
!!$    if (temp.lt.tsnow) then
!!$       psnow = prec
!!$    else
!!$       prain = prec
!!$    endif
!!$    
!!$    
!!$    smelt = 2.63 + 2.55* temp + 0.0912 * temp * prain !Snowmelt (mm/day)
!!$    smelt = amax1(smelt,0.)
!!$    smelt = amin1(smelt,(s+psnow))
!!$    ds = psnow - smelt
!!$    s = ds
!!$    
!!$    !     Water budget
!!$    !     ============
!!$    if (ts.le.tice) then !Frozen soil
!!$       g =  w0!Soil moisture freezes
!!$       w = 0.0
!!$       roff = smelt + prain !mm/day
!!$       evap = 0.0
!!$    else                !Non-frozen soil
!!$       w =  g + w0
!!$       g = 0.0
!!$       rimelt = 0.0
!!$       if (w .gt. wmax) then
!!$          rimelt = w - wmax !Runoff due to soil ice melting
!!$          w = wmax
!!$       endif
!!$       wapft = (w/wmax)
!!$       !runoff function
!!$       roff = runoff(wapft)       !Soil moisture runoff (roff, mm/day)
!!$       !# penman function 
!!$       evap = penman(p0,temp,rh,ae,rc2) !Actual evapotranspiration (evap, mm/day)
!!$       dw = prain + smelt - evap - roff
!!$       w = w + dw
!!$       if (w.gt.wmax) then
!!$          roff = roff + (w - wmax)
!!$          w = wmax
!!$       endif
!!$       if (w.lt.0.) w = 0.
!!$       roff = roff + rimelt !Total runoff
!!$    endif
!!$    
!!$    if(id .eq. 0) w = roff
!!$    if(id .eq. 1) w = evap
!!$    if(id .eq. 2) w = emax
!!$    
!!$  end function soil_water

  !=================================================================
  !=================================================================
  
  function penman (spre,temp,ur,rn,rc2) result(evap)
    use global_pars, only: r4, rcmin, rcmax
    use photo, only: tetens
    implicit none
    
    
    real(kind=r4),intent(in) :: spre                 !Surface pressure (mbar)
    real(kind=r4),intent(in) :: temp                 !Temperature (oC)
    real(kind=r4),intent(in) :: ur                   !Relative humidity (0-1,dimensionless)
    real(kind=r4),intent(in) :: rn                   !Radiation balance (W/m2)
    real(kind=r4),intent(in) :: rc2                  !Canopy resistence (s/m)
    
    real(kind=r4) :: evap                             !Evapotranspiration (mm/day)
    !     Parameters
    !     ----------
    real(kind=r4) :: ra, h5, t1, t2, es, es1, es2, delta_e, delta
    real(kind=r4) :: gama, gama2
    
    
    ra = rcmin
    h5 = 0.0275               !mb-1
    
    !     Delta
    !     -----
    t1 = temp + 1.
    t2 = temp - 1.
    es1 = tetens(t1)       !Saturation partial pressure of water vapour at temperature T
    es2 = tetens(t2)
    
    delta = (es1-es2)/(t1-t2) !mbar/oC
    !     
    !     Delta_e
    !     -------
    es = tetens (temp)
    delta_e = es*(1. - ur)    !mbar
    
    if ((delta_e.ge.(1./h5)-0.5).or.(rc2.ge.rcmax)) evap = 0.
    if ((delta_e.lt.(1./h5)-0.5).or.(rc2.lt.rcmax)) then
       !     Gama and gama2
       !     --------------
       gama  = spre*(1004.)/(2.45e6*0.622)
       gama2 = gama*(ra + rc2)/ra
       
       !     Real evapotranspiration
       !     -----------------------     
       evap = (delta* rn + (1.20*1004./ra)*delta_e)/(delta+gama2) !W/m2
       evap = evap*(86400./2.45e6) !mm/day
       evap = amax1(evap,0.)  !Eliminates condensation
    endif
  end function penman
  
  !=================================================================
  !=================================================================

  function available_energy(temp) result(ae)
    use global_pars, only: r4
    implicit none
    
    real(kind=r4),intent(in) :: temp
    real(kind=r4) :: ae
    
    ae = 2.895 * temp + 52.326 !from NCEP-NCAR Reanalysis data  
  end function  available_energy

  !=================================================================
  !=================================================================

  
  function runoff(wa) result(roff)
    use global_pars, only: r4
    implicit none
    
    real(kind=r4),intent(in) :: wa
    real(kind=r4):: roff
    
    !  roff = 38.*((w/wmax)**11.) ! [Eq. 10]
    roff = 11.5*((wa)**6.6) !from NCEP-NCAR Reanalysis data  
  end function  runoff
  
  !=================================================================
  !=================================================================
  
  function evpot2 (spre,temp,ur,rn) result(evap) 
    use global_pars, only: r4, rcmin, rcmax
    use photo, only: tetens
    implicit none
    
    !     Inputs
    
    real(kind=r4),intent(in) :: spre                 !Surface pressure (mb)
    real(kind=r4),intent(in) :: temp                 !Temperature (oC)
    real(kind=r4),intent(in) :: ur                   !Relative humidity (0-1,dimensionless)
    real(kind=r4),intent(in) :: rn                   !Radiation balance (W/m2)
    !     Output
    !     ------
    !     
    real(kind=r4) :: evap                 !Evapotranspiration (mm/day)
    !     Parameters
    !     ----------
    real(kind=r4) :: ra, t1, t2, es, es1, es2, delta_e, delta
    real(kind=r4) :: gama, gama2, rc
    
    ra = rcmin            !s/m
    
    !     Delta
    !     -----     
    t1 = temp + 1.
    t2 = temp - 1.
    es1 = tetens(t1)
    es2 = tetens(t2)
    delta = (es1-es2)/(t1-t2) !mb/oC
    !     
    !     Delta_e
    !     -------
    !     
    es = tetens (temp)
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
    evap =(delta*rn + (1.20*1004./ra)*delta_e)/(delta+gama2) !W/m2
    evap = evap*(86400./2.45e6) !mm/day
    evap = amax1(evap,0.)     !Eliminates condensation
  end function evpot2
  
  !=================================================================
  !=================================================================
  
end module water
