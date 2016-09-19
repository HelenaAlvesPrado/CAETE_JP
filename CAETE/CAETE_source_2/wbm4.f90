
subroutine wbm (&!in
                prec,&
                temp,&
                lsmk,&
                p0,&
                ca,&
                par,&!out
                 tmin     ,&
                 meanpp   ,&
                 seanpp   ,&
                 mphoto   ,&
                 maresp   ,&
                 meanhr   ,&
                 meancs   ,&
                 wsoil2   ,&
                 evapm    ,&
                 npp      ,&
                 photo    ,&
                 aresp    ,&
                 rcm      ,&
                 ave_wsoil,&
                 ave_evap ,&
                 ave_rc    &
                 )
! c
! c=======================================================================
! c
! c Water balance model (WBM). From monthly climatologies of
! c precipitation and surface temperature, the WBM calculates the
! c environmental variables.
! c
! c 05Jul2005, MDO: ae, rh & runoff are changed.
! c 11Jul2005, MDO: wsoil2 is written (for testing purpose).
! c 31Ago2006, DML: carbon cycle is included
! c
! c=======================================================================
! c
! c i/o variables
    implicit none

    integer, parameter :: nx=720, ny=360, nt=12
    real, parameter :: NO_DATA = -9999.000000
    real, intent(in) :: ca
    real, dimension(nx,ny,nt), intent(in) :: prec,&
                                             temp,&
                                             p0  ,&
                                             par


    real, dimension(nx,ny), intent(in) :: lsmk

    real, dimension(nx,ny), intent(out) :: tmin     ,&
                                           seanpp   ,&
                                           meanpp   ,&
                                           meanhr   ,&
                                           meancs   ,&
                                           mphoto   ,&
                                           maresp   ,&
                                           ave_wsoil,&
                                           ave_evap ,&
                                           ave_rc

    real, dimension(nx,ny,nt), intent(out) :: wsoil2 ,&
                                             photo   ,&
                                             evapm   ,&
                                             aresp   ,&
                                             npp     ,&
                                             rcm
! c
! c internal variables
    real, dimension(nx,ny,nt) :: tsoil, wsoil, gsoil, ssoil,&
                                 snowm, runom, emaxm, lai  ,&
                                 wg0, clit, csoil, hresp

    real, dimension(nx, ny) :: nppmin, nppmax

    real, parameter :: H = 1.0, diffu = 1.0368

    real :: tau, t0, t1, nppmes, laimes, ipar, auxs, wini,&
            gini, sini, ta, td, pr, spre, ae, wfim,&
            gfim, sfim, smes, rmes, emes, epmes, phmes,armes,&
            clmes, csmes, hrmes, rcmes, dwww, wmax

    integer :: i, j, k, nerro, kk, ice, n, mes

! c Soil temperature
! c ----------------
! c
      !H    =  !soil layer (m)
      !diffu = !soil thermal diffusivity (m2/mes)
    tau = (H**2)/(2.0*diffu)    !e-folding time (months)
    auxs = -100.0   !auxiliar for calculation of Snpp

!c for all grid points
    do i=1,nx
        do j=1,ny
        !c initialize soil temperature
            do k=1,nt
                tsoil(i,j,k) = NO_DATA
            enddo
! c
! c jp- spinup ________________________________________
! c only for land grid points

            if (nint(lsmk(i,j)) == 1) then
                t0 = 0.     !initialization
                do n = 1,1200 !100 yr (1200 months) run to attain equilibrium
                    k = mod(n,12)
                    if (k == 0) k = 12
                    t1 = t0 * exp(-1.0/tau) + (1.0 - exp(-1.0/tau)) * temp(i,j,k)
                    tsoil(i,j,k) = (t0 + t1)/2.0
                    t0 = t1
                enddo
            endif

        enddo
    enddo
! c
! c Water budget
! c ------------
! c
! c for all grid points
    do i=1,nx
        do j=1,ny

            if ((mod(j,ny) == 0) .and. (mod(i,10) == 0)) then
                write(*,*) 'water balance:',i
            endif

            tmin(i,j)   =       NO_DATA
            seanpp(i,j) =       NO_DATA
            meanpp(i,j) =       NO_DATA
            mphoto(i,j) =       NO_DATA
            maresp(i,j) =       NO_DATA
            meanhr(i,j) =       NO_DATA
            meancs(i,j) =       NO_DATA
            ave_wsoil(i,j) =    NO_DATA!calculates annual average soil water
            ave_evap(i,j) =     NO_DATA!calculates annual average evapotranspiration
            ave_rc(i,j) =       NO_DATA!calculates annual average canopy resistance
            do k=1,nt
                wsoil(i,j,k) =  NO_DATA!soil moisture (mm)
                gsoil(i,j,k) =  NO_DATA!soil ice (mm)
                ssoil(i,j,k) =  NO_DATA!soil snow (mm)
                snowm(i,j,k) =  NO_DATA!average snowmelt (mm/day)
                runom(i,j,k) =  NO_DATA!average runoff (mm/day)
                evapm(i,j,k) =  NO_DATA!average actual evapotranspiration (mm/day)
                emaxm(i,j,k) =  NO_DATA!average maximum evapotranspiration (mm/day)
                wg0(i,j,k) =    NO_DATA!soil moisture of the previous year (mm)
                wsoil2(i,j,k) = NO_DATA!for testing purpose
                rcm(i,j,k) =    NO_DATA!average canopy resistance (s/m)
                lai(i,j,k) =    NO_DATA
                photo(i,j,k) =  NO_DATA
                aresp(i,j,k) =  NO_DATA
                npp(i,j,k) =    NO_DATA
                clit(i,j,k) =   NO_DATA
                csoil(i,j,k) =  NO_DATA
                hresp(i,j,k) =  NO_DATA
            enddo



            if (nint(lsmk(i,j)) == 1) then
                wini  = 0.01  !soil moisture initial condition (mm)
                gini  = 0.0   !soil ice initial condition (mm)
                sini  = 0.0   !overland snow initial condition (mm)
! c
! c initialization
                do k=1,nt
                    wg0(i,j,k) = -1.0
                enddo
! c
! c start integration
                n = 0
10 continue
                n = n + 1
    ! c
    ! c pre-processing
                k = mod(n,12)
                if (k.eq.0) k = 12
                mes = k
                td = tsoil(i,j,k)
                ta = temp(i,j,k)
                pr = prec(i,j,k)
                ipar = par(i,j,k)
                spre = p0(i,j,k) !surface pressure (mb)
    ! c      ae = 2.26457*ta + 67.5876 !available energy (W/m2) [Eq. 8]
                ae = 2.895*ta + 52.326 !from NCEP-NCAR Reanalysis data
    ! c
    ! c monthly water budget
                call budget (mes,wini,gini,sini,td,ta,pr,&
                             spre,ae,ca,ipar,wfim,gfim,sfim,&
                             smes,rmes,emes,epmes,phmes,armes,&
                             nppmes,laimes,clmes,csmes,hrmes,rcmes)
    ! c
    ! c update variables
                wsoil(i,j,k) = wfim
                gsoil(i,j,k) = gfim
                ssoil(i,j,k) = sfim
                snowm(i,j,k) = smes
                runom(i,j,k) = rmes
                evapm(i,j,k) = emes
                emaxm(i,j,k) = epmes
                rcm(i,j,k) = rcmes
                lai(i,j,k) = laimes
                photo(i,j,k) = phmes
                aresp(i,j,k) = armes
                npp(i,j,k) = nppmes
                clit(i,j,k) = clmes
                csoil(i,j,k) = csmes
                hresp(i,j,k) = hrmes
                wini = wfim
                gini = gfim
                sini = sfim
    ! c
    ! c check if equilibrium is attained (k=12)
                if (k == 12) then
                    wmax = 500.
                    nerro = 0
                    do kk=1,12
                        dwww = (wsoil(i,j,kk) + gsoil(i,j,kk) - wg0(i,j,kk)) / wmax
                        if (abs(dwww) > 0.001) then
                            nerro = nerro + 1
                        endif
                    enddo
                    if (nerro > 0) then
                        do kk=1,12
                            wg0(i,j,kk) = wsoil(i,j,kk) + gsoil(i,j,kk)
                        enddo
                    else
                        goto 100
                    endif
                endif
                goto 10
100 continue
    ! c
    ! c Environmental variables
    ! c -----------------------
    ! c
    ! c initialize
                tmin(i,j) = 100.0
                seanpp(i,j) = 0.0
                meanpp(i,j) = 0.0
                mphoto(i,j) = 0.0
                maresp(i,j) = 0.0
                nppmin(i,j) = 100.0
                nppmax(i,j) = -100.0
                meanhr(i,j) = 0.0
                meancs(i,j) = 0.0
                ave_wsoil(i,j) = 0.0
                ave_evap(i,j) = 0.0
                ave_rc(i,j) = 0.0
    ! c
    ! c calculate tmin, meanpp, seanpp [Eqs 11, ...]
                do k=1,nt

                    if (temp(i,j,k) < tmin(i,j)) then
                        tmin(i,j) = temp(i,j,k)
                    endif

                    meanpp(i,j) = meanpp(i,j) + (npp(i,j,k)/12)
                    mphoto(i,j) = mphoto(i,j) + (photo(i,j,k)/12)
                    maresp(i,j) = maresp(i,j) + (aresp(i,j,k)/12)
                    meanhr(i,j) = meanhr(i,j) + (hresp(i,j,k)/12)
                    meancs(i,j) = meancs(i,j) + (csoil(i,j,k)/12)
                    ave_wsoil(i,j) = ave_wsoil(i,j) + (wsoil(i,j,k)/12)
                    ave_evap(i,j) = ave_evap(i,j) + (evapm(i,j,k)/12)
                    ave_rc(i,j) = ave_rc(i,j) + (rcm(i,j,k)/12)
                    if (npp(i,j,k) < nppmin(i,j)) nppmin(i,j) = npp(i,j,k)
                    if (npp(i,j,k) > nppmax(i,j)) nppmax(i,j) = npp(i,j,k)

                enddo
    ! c
                if (meanpp(i,j) > 0.0) then
                    seanpp(i,j) = (nppmax(i,j)-nppmin(i,j))/(meanpp(i,j))
                else
                    seanpp(i,j) = (nppmax(i,j)-nppmin(i,j))
                endif
                if(seanpp(i,j) > auxs) auxs = seanpp(i,j)    !in order to let Snpp dimensionless
    ! c
    ! c calculate wsoil2 (only for grid points without soil ice)
                ice = 0
                do k=1,12
                    if ((gsoil(i,j,k)/wmax) > 1.e-7) ice = 1
                enddo

                if (ice == 0) then
                    do k=1,12
                        wsoil2(i,j,k) = wsoil(i,j,k)/wmax
                    enddo
                endif
            endif
            ! if ((i.eq.20).and.(j.eq.45)) then
            !     write(*,22) meanpp(i,j), mphoto(i,j), maresp(i,j)
            !     22 format(f6.4)
            !     stop
            ! endif

        enddo
    enddo

!c Final determination of Snpp (loop again cause of auxs...)
    do i=1,nx
        do j=1,ny
            if (nint(lsmk(i,j)) == 1) then
                seanpp(i,j) = seanpp(i,j)/auxs
            endif
        enddo
    enddo

    return

end subroutine wbm

subroutine budget (month,w1,g1,s1,tsoil,temp,prec,p0,ae,ca,& !inputs
                   ipar,w2,g2,s2,smavg,ruavg,evavg,&
                   epavg,phavg,aravg,nppavg,laiavg,&
                   clavg,csavg,hravg,rcavg)             !outputs
! c=======================================================================
! c Surface water (soil moisture, snow and ice) budget for a single month.
! c I/O variables
! c input  month : actual month (1-12)
! c        w1    : initial (previous month last day) soil moisture storage (mm)
! c        g1    : initial soil ice storage (mm)
! c        s1    : initial overland snow storage (mm)
! c        tsoil : soil temperature (oC)
! c        temp  : surface air temperature (oC)
! c        prec  : precipitation (mm/day)
! c        p0    : surface pressure (mb)
! c        ae    : available energy (W/m2)
! c output w2    : final (last day) soil moisture storage (mm)
! c        g2    : final soil ice storage (mm)
! c        s2    : final overland snow storage (mm)
! c        smavg : snowmelt monthly average (mm/day)
! c        ruavg : runoff monthly average (mm/day)
! c        evavg : actual evapotranspiration monthly average (mm/day)
! c        epavg : maximum evapotranspiration monthly average (mm/day)
! c=======================================================================
! c i/o variables
    integer, intent(in) :: month
    real, intent(in) :: w1,g1,s1,tsoil,temp,prec,p0,ae,ca,ipar
    real, intent(out) :: w2,g2,s2,smavg,ruavg,evavg,epavg,rcavg,&
                         phavg,aravg,nppavg,laiavg,&
                         clavg,csavg,hravg
! c
! c internal variables
    real psnow,prain
    real w,g,s
    real rimelt,smelt,roff,evap,emax, dw
    integer ndmonth(12); data ndmonth /31,28,31,30,31,30,31,31,30,31,30,31/
!c carbon cyclE
    real ph,ar,nppa,nppb,laia,cl,cs,hr,f5, rc2
    real, parameter :: rh    = 0.685,& !from NCEP-NCAR Reanalysis data Umidade relativa
                       wmax  = 500.0,& !soil moisture availability (mm)
                       tsnow = -1.0 ,& !temperature threshold for snowfall (oC)
                       tice  = -2.5    !temperature threshold for soil freezing (oC)
! c
! c precipitation [Eq. 3]
    psnow = 0.0
    prain = 0.0
    if (temp.lt.tsnow) then
        psnow = prec/real(ndmonth(month)) !snowfall (mm/day)
    else
        prain = prec/real(ndmonth(month)) !rainfall (mm/day)
    endif
! c
! c initialization
    w = w1     !w = daily soil moisture storage (mm)
    g = g1     !g = daily soil ice storage (mm)
    s = s1     !s = daily overland snow storage (mm)
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
! c
! c numerical integration
    do i=1,ndmonth(month)
! c
! c carbon cycle (photosynthesis, plant respiration and NPP)
        call carbon1 (temp,p0,ca,ipar,&   !inputs
                      ph,ar,nppa,laia,f5)        !output
! c
! c maximum evapotranspiration (emax)
        call evpot2 (p0,temp,rh,ae,emax)
! c
! c snow budget
        smelt = 2.63 + 2.55*temp + 0.0912*temp*prain !snowmelt (mm/day) [Eq. 4]
        smelt = amax1(smelt,0.)
        smelt = amin1(smelt,s+psnow)
        ds = psnow - smelt ![Eq. 2]
        s = s + ds
! c
! c water budget
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
    ! c
! c Canopy resistance (based in Sellers et al. 1996; SiB2)
! c (rc2 ; s/m) [Eq. 32]
! c [NPP*2.64e-6 converts kgC/m2/yr to molCO2/m2/s]
! c [p0*100 convertes hPa (mb) to Pa]
! c
            nppb = amax1(nppa,0.05)
            rc2 = (ca/(0.9*(nppb*2.64e-6)*0.685*(p0*100)))
            call runoff (w,wmax,roff) !soil moisture runoff (roff, mm/day) [Eq. 10]
            call penman (p0,temp,rh,ae,rc2,evap) !actual evapotranspiration (evap, mm/day)
            dw = prain + smelt - evap - roff ![Eq. 1]
            w = w + dw
            if (w.gt.wmax) then
                roff = roff + (w - wmax)
                w = wmax
            endif
            if (w.lt.0.) w = 0.
            roff = roff + rimelt !total runoff
! c carbon cycle (Microbial respiration, litter and soil carbon)
            call carbon2 (tsoil,f5,evap,laia,& !input
                          cl,cs,hr)
        endif
! c
! c updating monthly values
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
    enddo
! c
! c final calculations
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

    return
end subroutine budget

subroutine penman (spre,temp,ur,rn,rc2,evap)
! c
! c Entradas
! c --------
! c spre   = pressao aa supeficie (mb)
! c temp   = temperatura (oC)
! c w      = grau de saturacao (0-1,adimensional)
! c ur     = umidade relativa  (0-1,adimensional)
! c rn     = saldo de radiacao (W m-2)
! c rc2    = resistencia do dossel (s/m)
! c
! c Saida
! c -----
! c evap  = evapotranspiracao (mm/dia)
! c
    real, intent(in) :: spre,temp,ur,rn,rc2
    real, intent(out) :: evap
! c
! c parametros
    real, parameter ::  ra =    100.0 ,&   !s/m
                        h5    = 0.0275 !mb-1

    real t1, t2, es1, es2, delta, delta_e, gama, gama2
! c
! c delta
    t1 = temp + 1.
    t2 = temp - 1.
    call tetens(t1,es1)
    call tetens(t2,es2)
    delta = (es1-es2)/(t1-t2) !mb/oC
! c
! c delta_e
    call tetens (temp,es)
    delta_e = es*(1. - ur) !mb
! c
    if ((delta_e >= (1./h5)-0.5).or.(rc2 >= 4500)) evap = 0.
    if ((delta_e < (1./h5)-0.5).or.(rc2 < 4500)) then
! c gama e gama2
        gama  = spre*(1004.)/(2.45e6*0.622)
        gama2 = gama*(ra + rc2)/ra
! c evapotranspiracao real
        evap = (delta*rn + (1.20*1004./ra)*delta_e)/(delta+gama2) ! W/m2
        evap = evap*(86400./2.45e6)                               ! mm/dia
        evap = amax1(evap,0.) !elimina condensacao
    endif
! c
    return
end subroutine penman


subroutine evpot2 (spre,temp,ur,rn,evap)
! c
! c Entradas
! c --------
! c spre   = pressao aa supeficie (mb)
! c temp   = temperatura (oC)
! c ur     = umidade relativa  (0-1,adimensional)
! c rn     = saldo de radiacao (W m-2)
! c
! c Saida
! c -----
! c evap  = evapotranspiracao potencial sem estresse (mm/dia)
! c
    real, intent(in) :: spre,temp,ur,rn
    real, intent(inout) :: evap
! c
! c parametros
    real, parameter :: ra =    100.0 ,&   !s/m
                       rcmin = 100.0      !s/m

    real t1, t2, es1, es2, delta, delta_e, gama, gama2, rc
! c
! c delta
    t1 = temp + 1.
    t2 = temp - 1.
    call tetens(t1,es1)
    call tetens(t2,es2)
    delta = (es1-es2)/(t1-t2) !mb/oC
! c
! c delta_e
    call tetens (temp,es)
    delta_e = es*(1. - ur) !mb
! c
! c resistencia estomatica
    rc = rcmin
! c
! c gama e gama2
    gama  = spre*(1004.)/(2.45e6*0.622)
    gama2 = gama*(ra + rc)/ra
! c
! c evapotranspiracao potencial sem estresse
    evap = (delta*rn + (1.20*1004./ra)*delta_e)/(delta+gama2) ! W/m2
    evap = evap*(86400./2.45e6)                               ! mm/dia
    evap = amax1(evap,0.) !elimina condensacao
    return
end subroutine evpot2


subroutine runoff (w,wmax,roff)
      real, intent(in) :: w, wmax
      real, intent(inout) :: roff
!c      roff = 38.*((w/wmax)**11.) ! [Eq. 10]
      roff = 11.5*((w/wmax)**6.6) !from NCEP-NCAR Reanalysis data
      return
  end subroutine runoff


subroutine tetens (t,es)
      real, intent(in) :: t
      real, intent(inout) :: es
        if (t > 0.) then
            es = 6.1078*exp((7.5*t/(237.3+t))*log(10.))
        else
            es = 6.1078*exp((9.5*t/(265.5+t))*log(10.))
        endif
      return
end subroutine tetens
