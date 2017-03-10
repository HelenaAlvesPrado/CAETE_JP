module water
  
  ! this module defines functions related to surface water balance
  implicit none
  private

  ! functions defined here:
  
  public ::              &
       soil_temp        ,&
       soil_temp_sub        ,&
       penman           ,&
       evpot2           ,&
       available_energy ,&
       runoff          

  
contains
  
  !=================================================================
  !=================================================================

  subroutine soil_temp_sub(temp, tsoil)
  ! Calcula a temperatura do solo. Aqui vamos mudar no futuro!
  ! a tsoil deve ter relacao com a et realizada...
  ! a profundidade do solo (H) e o coef de difusao (DIFFU) devem ser
  ! variaveis (MAPA DE SOLO?; agua no solo?)
  use global_pars
  implicit none
  integer(kind=i4),parameter :: m = ntimes
  
  real(kind=r4),dimension(m), intent( in) :: temp ! future __ make temps an allocatable array
  real(kind=r4),dimension(m), intent(out) :: tsoil
   
  ! internal vars
  
  integer(kind=i4) :: n, k
  real(kind=r4) :: t0 = 0.0
  real(kind=r4) :: t1 = 0.0

  tsoil = -9999.0

  do n=1,1200 !run to attain equilibrium
     k = mod(n,m)
     if (k.eq.0) k = 12
     t1 = (t0*exp(-1.0/TAU) + (1.0 - exp(-1.0/TAU)))*temp(k)
     tsoil(k) = (t0 + t1)/2.0
     t0 = t1
  enddo
end subroutine soil_temp_sub
  
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
