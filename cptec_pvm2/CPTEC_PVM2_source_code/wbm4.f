c234567
      subroutine wbm (prec,temp,lsmk,p0,ca,par,        		!input
     &                tmin,meanpp,seanpp,mphoto,maresp,  	!output
     &                meanhr,meancs,wsoil2,evapm,npp,
     &                photo,aresp,rcm,ave_wsoil,ave_evap,ave_rc)
c
c=======================================================================
c
c Water balance model (WBM). From monthly climatologies of
c precipitation and surface temperature, the WBM calculates the
c environmental variables.
c
c 05Jul2005, MDO: ae, rh & runoff are changed.
c 11Jul2005, MDO: wsoil2 is written (for testing purpose).
c 31Ago2006, DML: carbon cycle is included
c
c=======================================================================
c
c i/o variables
      parameter(nx=192,ny=96)
      real prec(nx,ny,12),temp(nx,ny,12),lsmk(nx,ny),p0(nx,ny)
      real tmin(nx,ny),seanpp(nx,ny),meanpp(nx,ny),meanhr(nx,ny),
     &     meancs(nx,ny),mphoto(nx,ny),maresp(nx,ny),ave_wsoil(nx,ny),
     &     ave_evap(nx,ny),ave_rc(nx,ny)
      real wsoil2(nx,ny,12),par(nx,ny,12)
      real photo(nx,ny,12),aresp(nx,ny,12),npp(nx,ny,12),
     & lai(nx,ny,12),clit(nx,ny,12),csoil(nx,ny,12),hresp(nx,ny,12)
c
c internal variables
      real H,diffu,tau,tsoil(nx,ny,12),t0,t1
      real wsoil(nx,ny,12),gsoil(nx,ny,12),ssoil(nx,ny,12),
     &       snowm(nx,ny,12),runom(nx,ny,12),evapm(nx,ny,12),
     &       emaxm(nx,ny,12),rcm(nx,ny,12)
      real wg0(nx,ny,12)
      real nppmes,laimes,ipar
      real nppmin(nx,ny),nppmax(nx,ny),meant(nx,ny),seat(nx,ny)
c
c Soil temperature
c ----------------
c
      H    = 1.0                  !soil layer (m)
      diffu = 4.e-7*(30.*86400.0) !soil thermal diffusivity (m2/mes)
      tau = (H**2)/(2.0*diffu)    !e-folding time (months)
      auxs=-100.0   !auxiliar for calculation of Snpp
c
c for all grid points
      do i=1,nx
      do j=1,ny
c
c initialize soil temperature
      do k=1,12
        tsoil(i,j,k) = -999.99
      enddo
c
c only for land grid points
      if (int(lsmk(i,j)).ne.0) then
      t0 = 0.     !initialization
      do n=1,1200 !100 yr (1200 months) run to attain equilibrium
        k = mod(n,12)
        if (k.eq.0) k = 12
        t1 = t0*exp(-1.0/tau) + (1.0 - exp(-1.0/tau))*temp(i,j,k)
        tsoil(i,j,k) = (t0 + t1)/2.0
        t0 = t1
      enddo
      endif
c
      enddo
      enddo
c
c Water budget
c ------------
c
c for all grid points
      do i=1,nx
      do j=1,ny
c
c write to track program execution
      if ((mod(j,ny).eq.0).and.(mod(i,10).eq.0))
     &  write(*,*) 'water balance:',i
c
c initialize variables (-999.99 is undef)
      tmin(i,j)   = -999.99
      seanpp(i,j) = -999.99
      meanpp(i,j) = -999.99
      mphoto(i,j) = -999.99
      maresp(i,j) = -999.99
      meanhr(i,j) = -999.99
      meancs(i,j) = -999.99
      ave_wsoil(i,j) = -999.99  !calculates annual average soil water
      ave_evap(i,j) = -999.99   !calculates annual average evapotranspiration
      ave_rc(i,j) = -999.99     !calculates annual average canopy resistance
      do k=1,12
        wsoil(i,j,k) = -999.99 !soil moisture (mm)
        gsoil(i,j,k) = -999.99 !soil ice (mm)
        ssoil(i,j,k) = -999.99 !soil snow (mm)
        snowm(i,j,k) = -999.99 !average snowmelt (mm/day)
        runom(i,j,k) = -999.99 !average runoff (mm/day)
        evapm(i,j,k) = -999.99 !average actual evapotranspiration (mm/day)
        emaxm(i,j,k) = -999.99 !average maximum evapotranspiration (mm/day)
          wg0(i,j,k) = -999.99 !soil moisture of the previous year (mm)
        wsoil2(i,j,k) = -999.99 !for testing purpose
        rcm(i,j,k) = -999.99 !average canopy resistance (s/m)
        lai(i,j,k) = -999.99
	photo(i,j,k) = -999.99
	aresp(i,j,k) = -999.99
        npp(i,j,k) = -999.99
	clit(i,j,k) = -999.99
	csoil(i,j,k) = -999.99
	hresp(i,j,k) = -999.99
      enddo
c
c only for land grid points
      if (int(lsmk(i,j)).ne.0) then
c
c set some variables
      wini  = 0.01  !soil moisture initial condition (mm)
      gini  = 0.0   !soil ice initial condition (mm)
      sini  = 0.0   !overland snow initial condition (mm)
c
c initialization
      do k=1,12
      wg0(i,j,k) = -1.0
      enddo
      spre = p0(i,j) !surface pressure (mb)
c
c start integration
      n = 0
   10 continue
      n = n + 1
c
c pre-processing
      k = mod(n,12)
      if (k.eq.0) k = 12
      mes = k
      td = tsoil(i,j,k)
      ta = temp(i,j,k)
      pr = prec(i,j,k)
      ipar = par(i,j,k)
c      ae = 2.26457*ta + 67.5876 !available energy (W/m2) [Eq. 8]
      ae = 2.895*ta + 52.326 !from NCEP-NCAR Reanalysis data
c
c monthly water budget
      call budget (mes,wini,gini,sini,td,ta,pr,spre,ae,ca,ipar,
     &                 wfim,gfim,sfim,smes,rmes,emes,epmes,
     &                 phmes,armes,nppmes,laimes,
     &                 clmes,csmes,hrmes,rcmes)
c
c update variables
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
c
c check if equilibrium is attained (k=12)
      if (k.eq.12) then
         wmax = 500.
         nerro = 0
         do kk=1,12
            dwww = (wsoil(i,j,kk)+gsoil(i,j,kk)-wg0(i,j,kk))/wmax
            if (abs(dwww).gt.0.001) nerro = nerro + 1
         enddo
         if (nerro.ne.0) then
            do kk=1,12
               wg0(i,j,kk) = wsoil(i,j,kk) + gsoil(i,j,kk)
            enddo
         else
            goto 100
         endif
      endif
      goto 10
  100 continue
c
c Environmental variables
c -----------------------
c
c initialize
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
c
c calculate tmin, meanpp, seanpp [Eqs 11, ...]
      do k=1,12
         if (temp(i,j,k).lt.tmin(i,j)) tmin(i,j) = temp(i,j,k)
	 meanpp(i,j) = meanpp(i,j) + (npp(i,j,k)/12)
	 mphoto(i,j) = mphoto(i,j) + (photo(i,j,k)/12)
	 maresp(i,j) = maresp(i,j) + (aresp(i,j,k)/12)
	 meanhr(i,j) = meanhr(i,j) + (hresp(i,j,k)/12)
	 meancs(i,j) = meancs(i,j) + (csoil(i,j,k)/12)
         ave_wsoil(i,j) = ave_wsoil(i,j) + (wsoil(i,j,k)/12)
         ave_evap(i,j) = ave_evap(i,j) + (evapm(i,j,k)/12)
         ave_rc(i,j) = ave_rc(i,j) + (rcm(i,j,k)/12)
         if (npp(i,j,k).lt.nppmin(i,j)) nppmin(i,j) = npp(i,j,k)
         if (npp(i,j,k).gt.nppmax(i,j)) nppmax(i,j) = npp(i,j,k)
      enddo
c
      if (meanpp(i,j).gt.0.0) then
      seanpp(i,j) = (nppmax(i,j)-nppmin(i,j))/(meanpp(i,j))
      else
      seanpp(i,j) = (nppmax(i,j)-nppmin(i,j))
      endif
      if(seanpp(i,j).gt.auxs) auxs=seanpp(i,j)	!in order to let Snpp dimensionless
c
c calculate wsoil2 (only for grid points without soil ice)
      ice = 0
      do k=1,12
         if (gsoil(i,j,k)/wmax.gt.1.e-7) ice = 1
      enddo
      if (ice.eq.0) then
         do k=1,12
            wsoil2(i,j,k) = wsoil(i,j,k)/wmax
         enddo
      endif
c
c close land "if" and "do" loops
      endif
c
c      if ((i.eq.20).and.(j.eq.45)) then
c      write(*,22) meanpp(i,j), mphoto(i,j), maresp(i,j)
c   22 format(f6.4)
c      stop
c      endif
c
      enddo
      enddo
c
c Final determination of Snpp (loop again cause of auxs...)
      do i=1,nx
      do j=1,ny
        if (int(lsmk(i,j)).ne.0) then
         seanpp(i,j) = seanpp(i,j)/auxs
        endif
      enddo
      enddo

      return
      end
c
c=======================================================================
c234567
      subroutine budget (month,w1,g1,s1,tsoil,temp,prec,p0,ae,ca,
     &                    ipar,w2,g2,s2,smavg,ruavg,evavg,	 	!output
     &                         epavg,phavg,aravg,nppavg,laiavg,	 	!output
     &                         clavg,csavg,hravg,rcavg)		 	!output
c
c=======================================================================
c
c Surface water (soil moisture, snow and ice) budget for a single month.
c
c I/O variables
c -------------
c input  month : actual month (1-12)
c        w1    : initial (previous month last day) soil moisture storage (mm)
c        g1    : initial soil ice storage (mm)
c        s1    : initial overland snow storage (mm)
c        tsoil : soil temperature (oC)
c        temp  : surface air temperature (oC)
c        prec  : precipitation (mm/day)
c        p0    : surface pressure (mb)
c        ae    : available energy (W/m2)
c output w2    : final (last day) soil moisture storage (mm)
c        g2    : final soil ice storage (mm)
c        s2    : final overland snow storage (mm)
c        smavg : snowmelt monthly average (mm/day)
c        ruavg : runoff monthly average (mm/day)
c        evavg : actual evapotranspiration monthly average (mm/day)
c        epavg : maximum evapotranspiration monthly average (mm/day)
c
c=======================================================================
c
c i/o variables
      integer month
      real w1,g1,s1,tsoil,temp,prec,p0,ae,ca,ipar,
     &     w2,g2,s2,smavg,ruavg,evavg,epavg,rcavg,
     &     phavg,aravg,nppavg,laiavg,
     &     clavg,csavg,hravg
c
c internal variables
      real rh,wmax,tsnow,tice
      real psnow,prain
      real w,g,s
      real rimelt,smelt,roff,evap,emax
      integer ndmonth(12) !number of days for each month
      data ndmonth /31,28,31,30,31,30,31,31,30,31,30,31/
c carbon cycle
      real ph,ar,nppa,nppb,laia,cl,cs,hr
c
c parameters
c      rh    = 0.6   !relative humidity (adimensional)
      rh    = 0.685 !from NCEP-NCAR Reanalysis data
      wmax  = 500.0 !soil moisture availability (mm)
      tsnow = -1.0  !temperature threshold for snowfall (oC)
      tice  = -2.5  !temperature threshold for soil freezing (oC)
c
c precipitation [Eq. 3]
      psnow = 0.0
      prain = 0.0
      if (temp.lt.tsnow) then
        psnow = prec/real(ndmonth(month)) !snowfall (mm/day)
      else
        prain = prec/real(ndmonth(month)) !rainfall (mm/day)
      endif
      w1 = 0
c
c initialization
      w = w1 	!w = daily soil moisture storage (mm)
      g = g1 	!g = daily soil ice storage (mm)
      s = s1 	!s = daily overland snow storage (mm)
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
c
c numerical integration
      do i=1,ndmonth(month)
c
c carbon cycle (photosynthesis, plant respiration and NPP)
      call carbon1 (temp,p0,w,wmax,ca,ipar,	!input
     &              ph,ar,nppa,laia,f5)		!output
c
c maximum evapotranspiration (emax)
      call evpot2 (p0,temp,rh,ae,emax)
c
c snow budget
      smelt = 2.63 + 2.55*temp + 0.0912*temp*prain !snowmelt (mm/day) [Eq. 4]
      smelt = amax1(smelt,0.)
      smelt = amin1(smelt,s+psnow)
      ds = psnow - smelt ![Eq. 2]
      s = s + ds
c
c water budget
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
c
c Canopy resistance (based in Sellers et al. 1996; SiB2)
c (rc2 ; s/m) [Eq. 32]
c [NPP*2.64e-6 converts kgC/m2/yr to molCO2/m2/s]
c [p0*100 convertes hPa (mb) to Pa]
c
	nppb = amax1(nppa,0.05)
      	rc2 = (ca/(0.9*(nppb*2.64e-6)*0.685*(p0*100)))
	call runoff (w,wmax,roff) !soil moisture runoff (roff, mm/day) [Eq. 10]
        call penman (p0,temp,w,wmax,rh,ae,rc2,evap) !actual evapotranspiration (evap, mm/day)
	dw = prain + smelt - evap - roff ![Eq. 1]
        w = w + dw
        if (w.gt.wmax) then
          roff = roff + (w - wmax)
          w = wmax
        endif
        if (w.lt.0.) w = 0.
        roff = roff + rimelt !total runoff
c carbon cycle (Microbial respiration, litter and soil carbon)
	call carbon2 (tsoil,f5,evap,laia, !input
     &                cl,cs,hr)              !output
      endif
c
c updating monthly values
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
c
      enddo
c
c final calculations
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
c
      return
      end
c
c======================================================================
c
      subroutine penman (spre,temp,w,wmax,ur,rn,rc2,evap)
c
c Entradas
c --------
c spre   = pressao aa supeficie (mb)
c temp   = temperatura (oC)
c w      = grau de saturacao (0-1,adimensional)
c ur     = umidade relativa  (0-1,adimensional)
c rn     = saldo de radiacao (W m-2)
c rc2    = resistencia do dossel (s/m)
c
c Saida
c -----
c evap  = evapotranspiracao (mm/dia)
c
      real spre,temp,w,wmax,ur,rn,rc2,evap
c
c parametros
      ra =    100.   !s/m
      h5    = 0.0275 !mb-1
c
c delta
      t1 = temp + 1.
      t2 = temp - 1.
      call tetens(t1,es1)
      call tetens(t2,es2)
      delta = (es1-es2)/(t1-t2) !mb/oC
c
c delta_e
      call tetens (temp,es)
      delta_e = es*(1. - ur) !mb
c
      if ((delta_e.ge.(1./h5)-0.5).or.(rc2.ge.4500)) evap = 0.
      if ((delta_e.lt.(1./h5)-0.5).or.(rc2.lt.4500)) then
c gama e gama2
        gama  = spre*(1004.)/(2.45e6*0.622)
        gama2 = gama*(ra + rc2)/ra
c evapotranspiracao real
        evap = (delta*rn + (1.20*1004./ra)*delta_e)/(delta+gama2) ! W/m2
        evap = evap*(86400./2.45e6)                               ! mm/dia
        evap = amax1(evap,0.) !elimina condensacao
      endif
c
      return
      end
c
c======================================================================
c
      subroutine evpot2 (spre,temp,ur,rn,evap)
c
c Entradas
c --------
c spre   = pressao aa supeficie (mb)
c temp   = temperatura (oC)
c ur     = umidade relativa  (0-1,adimensional)
c rn     = saldo de radiacao (W m-2)
c
c Saida
c -----
c evap  = evapotranspiracao potencial sem estresse (mm/dia)
c
      real spre,temp,ur,rn,evap
c
c parametros
      ra =    100.   !s/m
      rcmin = 100.   !s/m
c
c delta
      t1 = temp + 1.
      t2 = temp - 1.
      call tetens(t1,es1)
      call tetens(t2,es2)
      delta = (es1-es2)/(t1-t2) !mb/oC
c
c delta_e
      call tetens (temp,es)
      delta_e = es*(1. - ur) !mb
c
c resistencia estomatica
      rc = rcmin
c
c gama e gama2
      gama  = spre*(1004.)/(2.45e6*0.622)
      gama2 = gama*(ra + rc)/ra
c
c evapotranspiracao potencial sem estresse
      evap = (delta*rn + (1.20*1004./ra)*delta_e)/(delta+gama2) ! W/m2
      evap = evap*(86400./2.45e6)                               ! mm/dia
      evap = amax1(evap,0.) !elimina condensacao
c
      return
      end
c
c=====================================================================
c
      subroutine runoff (w,wmax,roff)
      real w,roff
c      roff = 38.*((w/wmax)**11.) ! [Eq. 10]
      roff = 11.5*((w/wmax)**6.6) !from NCEP-NCAR Reanalysis data 
      return
      end
c
c=====================================================================
c
      subroutine tetens (t,es)
      real t,es
      if (t.ge.0.) then
      es = 6.1078*exp((7.5*t/(237.3+t))*log(10.))
      else
      es = 6.1078*exp((9.5*t/(265.5+t))*log(10.))
      endif
      return
      end

