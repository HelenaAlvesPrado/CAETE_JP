
subroutine budget(wi, gi, si, temp, prec, p0, rh, par, ca, tsoil,&
     wf, gf, sf, runoff, evap, emax, rc2, &
     ph, ar, nppa, laia,&
     cl, cs, hr)
  

  ! i/o
  real, intent( in) :: wi, gi, si ! initial water/ice/snow soil content
  real, intent( in) :: temp, prec, p0, rh, par, ca, tsoil   ! = 0.685 !from NCEP-NCAR Reanalysis data
  real, intent(out) :: wf, gf, sf, runoff, evap, emax, rc2
  real, intent(out) :: ph, ar, nppa, laia
  real, intent(out) :: cl, cs, hr
  ! parameters
  real, parameter ::  wmax  = 500.0  !soil moisture availability (mm)
  real, parameter ::  tsnow = -1.0   !temperature threshold for snowfall (oC)
  real, parameter ::  tice  = -2.5   !temperature threshold for soil freezing (oC)

  ! internal vars
  real :: psnow = 0.0
  real :: prain = 0.0
  real :: w, g, s, smelt, dw, f5, ae
  
  if (temp.lt.tsnow) then
     psnow = prec !snowfall (mm/day)
  else
     prain = prec !rainfall (mm/day)
  endif
  
  !c initialization
  w  = wi   !soil moisture initial condition (mm)
  g  = gi   !soil ice initial condition (mm)
  s  = si   !overland snow initial condition (mm)
  
  !NPP +  potential evapotranspiration
  call prod(temp,p0,w,wmax,ca,par,ph,ar,nppa,laia,f5)
  call evpot2(p0,temp,rh,par,emax)
 
  !c snow budget
  smelt = 2.63 + 2.55*temp + 0.0912*temp*prain !snowmelt (mm/day) [Eq. 4]
  smelt = amax1(smelt,0.0)
  smelt = amin1(smelt,psnow)
  s = psnow - smelt
  
  !c water budget
  if (tsoil.le.tice) then !frozen soil
     g = w !soil moisture freezes
     w = 0.0
     runoff = smelt + prain
     evap = 0.0
     rc2 = 100.0 !default value, equal to aerodynamic resistance (below)
  else
     w = w + g !soil ice melts
     g = 0.0
     rimelt = 0.0
     if (w.gt.wmax) then
        rimelt = w - wmax !runoff due to soil ice melting
        w = wmax
     endif

     call canopy_resistence(nppa, ca, p0, rc2)
     call available_energy(temp, ae)
     call runoff_c(w,wmax,runoff) !soil moisture runoff (roff, mm/day) [Eq. 10]
     call penman (p0,temp,rh,ae,rc2,evap) !actual evapotranspiration (evap, mm/day)
     dw = prain + smelt - evap - runoff ![Eq. 1]
     w = w + dw
     if (w.gt.wmax) then
        runoff = runoff + (w - wmax)
        w = wmax
     endif
     if (w.lt.0.) w = 0.
     runoff = runoff + rimelt !total runoff
  endif
  call carbon_hr(tsoil,f5,evap,laia,cl,cs,hr)

  wf = w
  gf = g
  sf = s
end subroutine budget



subroutine wbm5(prec, temp, p0, ca, par,&
     wsoil, gsoil, ssoil,total_runoff, evapm, rcm, emaxm, &
     photo, aresp, npp, lai,&
     clit, csoil, hresp)

  implicit none
  
  !parameters
  integer, parameter :: m = 12
  
  ! i/o

  real, dimension(m), intent(in) :: prec, temp, p0, par
  real, intent(in) :: ca

  real, dimension(m), intent(out) ::  wsoil, gsoil, ssoil, total_runoff, evapm, rcm
  real, dimension(m), intent(out) ::  emaxm, lai, photo, aresp, npp, clit, csoil, hresp
  !real, dimension(m), intent(out) :: 
  !c set some internal variables
  real, dimension(m) :: wg0, spre, tsoil

  integer :: i, k, d, kk, n, nerro
  real :: td, ta, pr, ipar, dwww, rh, wmax, ca1
  real :: wini  = 0.01  !soil moisture initial condition (mm)
  real :: gini  = 0.0   !soil ice initial condition (mm)
  real :: sini  = 0.0   !overland snow initial condition (mm)
   

  real :: wf, gf, sf, runoff, evap, emax,rc2
  real :: ph, ar, nppa, laia, cl, cs, hr, spr

  ca1  = ca *  101.325
  !c initialization
  do k=1,m
     wg0(k) = -1.0
     spre(k) = p0(k) * 0.01 ! converts Pa to mbar
     tsoil(k) = 0.
     wsoil(k) = 0.
     gsoil(k) = 0.
     ssoil(k) = 0.
     total_runoff(k) = 0.
     evapm(k) = 0.
     rcm(k) = 0.
     emaxm(k) = 0.
     lai(k) = 0.
     photo(k) = 0.
     aresp(k) = 0.
     npp(k) = 0.
     clit(k) = 0.
     csoil(k) = 0.
     hresp(k) = 0.
  enddo
  !c
  !c start integration
  call soil_temp(temp, tsoil)
  n = 0
10 continue
  n = n + 1
  !c
  !c pre-processing
  k = mod(n,12)
  if (k.eq.0) k = 12
  
  select case(k)
  case(1)
     d = 31
  case(2)
     d = 28
  case(3)
     d = 31
  case(4)
     d = 30
  case(5)
     d = 31
  case(6)
     d = 30
  case(7)
     d = 31
  case(8)
     d = 31
  case(9)
     d = 30
  case(10)
     d = 31
  case(11)
     d = 30
  case(12)
     d = 31
  end select
  
  td = tsoil(k)
  ta = temp(k)
  pr = prec(k)/real(d)
  ipar = par(k)
  spr = spre(k)
  
  do i=1,d

     runoff = 0.
     evap = 0.
     emax= 0.
     rc2 = 0.
     wf = 0.
     gf = 0.
     sf = 0.
     ph = 0.
     ar = 0.
     nppa = 0.
     laia = 0.
     cl = 0.
     cs = 0.
     hr = 0.
     rh = 0.700
     
     call budget(wini, gini, sini, ta, pr, spr ,rh, ipar, ca1, td,&
          wf, gf, sf, runoff, evap, emax,rc2, &
          ph, ar, nppa, laia,&
          cl, cs, hr)
     
     wsoil(k) = wsoil(k) + (wf / real(d))
     gsoil(k) = gsoil(k) + (gf / real(d))
     ssoil(k) = ssoil(k) + (sf / real(d))
     total_runoff(k) = total_runoff(k) + (runoff / real(d))
     evapm(k) = evapm(k) + (evap / real(d))
     emaxm(k) = emaxm(k) + (emax / real(d))
     rcm(k) = rcm(k) + (rc2/real(d))
     lai(k) = lai(k) + (laia/real(d))
     photo(k) = (photo(k) + (ph / 365.) * 12)
     aresp(k) = (aresp(k) + (ar / 365.) * 12)
     npp(k) = (npp(k) + (nppa / 365.) * 12)
     clit(k) = (clit(k) + (cl / 365.) * 12)
     csoil(k) = (csoil(k) + (cs / 365.) * 12)
     hresp(k) = (hresp(k) + (hr / 365.))
     wini = wf
     gini = gf
     sini = sf

  end do

  
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
  
end subroutine wbm5
