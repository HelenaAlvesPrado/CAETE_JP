
subroutine photosynthesis_rate(grd)
  !returns instantaneous photosynthesis rate at leaf level (molCO2/m2/s)
  
  use global_pars
  use photo_par
  use grdcell
  
  implicit none

  type(grdcell),intent(inout) :: grd
  integer(kind=i4) :: p
  real(kind=r4),dimension(ntimes) :: temp = grd%temp
  real(kind=r4),dimension(ntimes) :: p0   = grd%p0
  real(kind=r4),dimension(ntimes) :: w    = grd%w
  real(kind=r4),dimension(ntimes) :: rh   = grd%rh
  real(kind=r4),dimension(ntimes) :: ipar = grd%ipar
  
  
  real(kind=r8),dimension(ntimes) :: f2,f3            !Michaelis-Menten CO2/O2 constant (Pa)
  real(kind=r8),dimension(ntimes) :: mgama           !Photo-respiration compensation point (Pa)
  real(kind=r8),dimension(ntimes) :: rmax, r
  real(kind=r8),dimension(ntimes) :: ci
  real(kind=r8),dimension(npls,ntimes) :: jp1, jp2, jp, jc, jl, je
  real(kind=r8),dimension(npls,ntimes) :: b,c,delta, delta2
  real(kind=r8),dimension(npls,ntimes) :: f1a


  !ENVIRONMENTAL STATE (COMMOM TO ALL PFTS)
  !============================================================
  !Photo-respiration compensation point (Pa)
  mgama = p3/(p8*(p9**(p10*(temp-p11))))
  
  !Michaelis-Menten CO2 constant (Pa)     
  f2 = p12*(p13**(p10*(temp-p11)))

  !Michaelis-Menten O2 constant (Pa)
  f3 = p14*(p15**(p10*(temp-p11)))

  !Saturation vapour pressure (hPa)
  call tetens(temp,es)

  !Saturated mixing ratio (kg/kg)     
  rmax = 0.622*(es/(p0-es))
     
  !Moisture deficit at leaf level (kg/kg)    
  r = -0.315*rmax
  
  !Internal leaf CO2 partial pressure (Pa)
  ci = p19* (1.-(r/p20)) * (ca-mgama) + mgama
  !==============================================================

  !Rubisco carboxilation limited photosynthesis rate (molCO2/m2/s)
  do p = 1, npls
     
     jc(p,:) = vm(p)*((ci-mgama)/(ci+(f2*(1.+(p3/f3)))))
  
!!$  !Light limited photosynthesis rate (molCO2/m2/s) 
     jl(p,:) = p4*(1.0-p5)*ipar*((ci-mgama)/(ci+(p6*mgama)))
  
  !     Transport limited photosynthesis rate (molCO2/m2/s)
  !     ---------------------------------------------------
     je(p,:) = p7*vm(p)
  
  !     Jp (minimum between jc and jl)
     !     ------------------------------
  enddo
  b = (-1.)*(jc+jl)
  c = jc*jl
  delta = (b**2)-4.0*a*c
  
  jp1=(-b-(sqrt(delta)))/(2.0*a)
  jp2=(-b+(sqrt(delta)))/(2.0*a)
  jp= amin1(jp1,jp2)
  
  !     Leaf level gross photosynthesis (minimum between jc, jl and je)
  !     ---------------------------------------------------------------
  a2 = 0.93
  b2 = (-1.)*(jp+je)
  c2 = jp*je
  delta2 = (b2**2)-4.0*a2*c2
  
  j1=(-b2-(sqrt(delta2)))/(2.0*a2)
  j2=(-b2+(sqrt(delta2)))/(2.0*a2)
  f1a = amin1(j1,j2)

  grd%photo_rate = f1a
  
  return
end subroutine photosynthesis

subroutine tetens (t,es)  !Saturation Vapor Pressure (hPa)!
  use global_pars
  implicit none
  
  real(KIND=R4),dimension(ntimes),intent( in) :: t
  real(KIND=R4),dimension(ntimes),intent(out) :: es
  
  where (t .ge .0.)
     es = 6.1121 * exp((18.678-(t/234.5))*(t/(257.14+t))) ! mbar == hPa
  elsewhere
     es = 6.1115 * exp((23.036-(t/333.7))*(t/(279.82+t))) ! mbar == hPa
  end where
  
  return
end subroutine tetens
