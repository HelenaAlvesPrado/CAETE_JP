c23456789
!
      subroutine wbm (prec,temp,lsmk,p0,ca,par,        		                               !input
     &                tmin,meanpp,seanpp,mphoto,maresp,           	                       !output
     &                meanhr,meancs,wsoil2,evapm,npp,
     &                photo,aresp,cres,ave_wsoil,ave_evap,ave_cres,
     &                cres1,cres2,cres3,ave_cres1,ave_cres2,ave_cres3
     &                npp1,npp2,npp3,meanpp1,meanpp2,meanpp3,wsoil,
     &                runom,mnpp1_pon,mnpp2_pon,mnpp3_pon)
!
!==============================================================================================
!
! Water balance model (WBM).
! From monthly climatologies of precipitation and surface temperature, the WBM calculates the
! environmental variables.
!
! 05 Jul 2005, MDO: ae, rh & runoff are changed
! 11 Jul 2005, MDO: wsoil2 is written (for testing purpose)
! 31 Ago 2006, DML: carbon cycle is included
!
!==============================================================================================
!
! Input/Output Variables
! ----------------------
!
      integer, parameter :: nx=192, ny=96
      integer counter
      real, parameter :: NO_DATA = -9999.0
      real prec(nx,ny,12),temp(nx,ny,12),lsmk(nx,ny),p0(nx,ny)
      real tmin(nx,ny),seanpp(nx,ny),meanpp(nx,ny),meanhr(nx,ny),
     &     meancs(nx,ny),mphoto(nx,ny),maresp(nx,ny),ave_wsoil(nx,ny),
     &     ave_evap(nx,ny),ave_cres(nx,ny)
      real wsoil2(nx,ny,12),par(nx,ny,12)
      real photo(nx,ny,12),aresp(nx,ny,12),npp(nx,ny,12),
     &     lai(nx,ny,12),clit(nx,ny,12),csoil(nx,ny,12),hresp(nx,ny,12)
      real cres1(nx,ny,12),cres2(nx,ny,12),cres3(nx,ny,12)
      real ave_cres1(nx,ny),ave_cres2(nx,ny),ave_cres3(nx,ny)
      real npp1(nx,ny,12),npp2(nx,ny,12),npp3(nx,ny,12)
      real meanpp1(nx,ny),meanpp2(nx,ny),meanpp3(nx,ny)
      real npptotal(nx,ny)
!
! Internal Variables
! ------------------
!
      real H,diffu,tau,tsoil(nx,ny,12),t0,t1
      real wsoil(nx,ny,12),gsoil(nx,ny,12),ssoil(nx,ny,12),
     &     snowm(nx,ny,12),runom(nx,ny,12),evapm(nx,ny,12),
     &     emaxm(nx,ny,12),cres(nx,ny,12)
      real wg0(nx,ny,12)
      real nppmes,laimes,ipar
      real nppmin(nx,ny),nppmax(nx,ny),meant(nx,ny),seat(nx,ny)
      real nppmes1,nppmes2,nppmes3
      real mnpp1_f(nx,ny),mnpp2_f(nx,ny),mnpp3_f(nx,ny)
      real mnpp1_pon(nx,ny),mnpp2_pon(nx,ny),mnpp3_pon(nx,ny)
      real one
!
! Counting
! --------

      counter = 0
!
! Soil temperature
! ================
!
      H    = 1.0                                                                              !Soil layer(m)
      diffu = 4.e-7*(30.*86400.0)                                                              !Soil thermal diffusivity(m2/month)
      tau = (H**2)/(2.0*diffu)                                                               !E-folding time(months)
      auxs=-100.0                                                                           !Auxiliar for calculation of Snpp
!
! For all grid points
! -------------------
!
      do i=1,nx
      do j=1,ny
!
! Initialize soil temperature
! ---------------------------
!
      do k=1,12
         tsoil(i,j,k) = NO_DATA                                                                   !Undefined
      enddo
!
! Only for land grid points
! -------------------------
!
      if (int(lsmk(i,j)).ne.0) then
      t0 = 0.                                                                                  !Initialization
      do n=1,1200                                                                           !1200 months run to attain equilibrium
        k = mod(n,12)
        if (k.eq.0) k = 12
        t1 = t0*exp(-1.0/tau) + (1.0 - exp(-1.0/tau))*temp(i,j,k)
        tsoil(i,j,k) = (t0 + t1)/2.0
        t0 = t1
      enddo
      endif
!
      enddo
      enddo
!
! ============
! Water budget
! ============
!
! For all grid points
! -------------------
!
      do i=1,nx
      do j=1,ny
!
! Write to track program execution
! --------------------------------
!
      if ((mod(j,ny).eq.0).and.(mod(i,10).eq.0))
     &  write(*,*) 'working:',i
!
! Initialize variables
! --------------------
!
      tmin(i,j)   =  NO_DATA   !Minimal temperature
      seanpp(i,j) =  NO_DATA   
      meanpp(i,j) =  NO_DATA   !NPP
      meanpp1(i,j) = NO_DATA   !NPP_PFT1
      meanpp2(i,j) = NO_DATA   !NPP_PFT2
      meanpp3(i,j) = NO_DATA   !NPP_PFT3
      mphoto(i,j) =  NO_DATA   !Photosynthesis
      maresp(i,j) =  NO_DATA    !Autotrophic Respiration
      meanhr(i,j) =  NO_DATA    !Heterotrophic Respiration
      meancs(i,j) =  NO_DATA    !Carbon Soil
      ave_wsoil(i,j) = NO_DATA !Water Soil
      ave_evap(i,j) =  NO_DATA  !Evapotranspiration
      ave_cres(i,j) =  NO_DATA  !Canopy Resistance
      ave_cres1(i,j) = NO_DATA !Canopy Resistance_PFT1
      ave_cres2(i,j) = NO_DATA !Canopy Resistance_PFT2
      ave_cres3(i,j) = NO_DATA !Canopy Resistance_PFT3
!
      do k=1,12
!
        wsoil(i,j,k)   =       NO_DATA                                                          !Soil Moisture(mm)
        gsoil(i,j,k)   =       NO_DATA                                                          !Soil Ice(mm)
        ssoil(i,j,k)   =       NO_DATA                                                          !Soil Snow(mm)
        snowm(i,j,k)   =       NO_DATA                                                          !Snowmelt(mm/day)
        runom(i,j,k)   =       NO_DATA                                                          !Runoff(mm/day)
        evapm(i,j,k)   =       NO_DATA                                                          !Actual Evapotranspiration(mm/day)
        emaxm(i,j,k)   =       NO_DATA                                                          !Maximum Evapotranspiration(mm/day)
          wg0(i,j,k)   =       NO_DATA                                                          !Moisture of the previous year (mm)
        wsoil2(i,j,k)  =       NO_DATA                                                          !For testing purpose
         cres(i,j,k)   =       NO_DATA                                                          !Canopy Resistance (s/m)
        cres1(i,j,k)   =       NO_DATA                                                          !Canopy Resistance_PFT1 (s/m)
        cres2(i,j,k)   =        NO_DATA                                                         !Canopy Resistance_PFT2 (s/m)
        cres3(i,j,k)   =        NO_DATA                                                          !Canopy Resistance_PFT3 (s/m)
        lai(i,j,k)     =         NO_DATA                                                          !Leaf Area Index
        photo(i,j,k)   =         NO_DATA                                                       !Photosynthesis
        aresp(i,j,k)   =         NO_DATA                                                         !Autotrophic Respiration
        npp(i,j,k)     =         NO_DATA                                                         !NPP
        npp1(i,j,k)    =         NO_DATA                                                         !NPP_PFT1
        npp2(i,j,k)    =         NO_DATA                                                          !NPP_PFT2
        npp3(i,j,k)    =         NO_DATA                                                          !NPP_PFT3
        clit(i,j,k)    =         NO_DATA                                                !Carbon Litter
        csoil(i,j,k)   =         NO_DATA                                                          !Carbon Soil
        hresp(i,j,k)   =         NO_DATA                                                         !Heterotrophic Respiration
!
      enddo
!
! Only for land grid points
! -------------------------
!
      if (int(lsmk(i,j)).ne.0) then
!
! Set some variables
! ------------------
!
      wini  = 0.01                                                                             !Soil Moisture_Initial Condition(mm)
      gini  = 0.0                                                                              !Soil Ice_Initial Condition(mm)
      sini  = 0.0                                                                              !Overland Snow_Initial Condition(mm)
!
! Initialization
! --------------
!
      do k=1,12
      wg0(i,j,k) = -1.0                                                                                              !wb0?
      enddo
      spre = p0(i,j)                                                                           !Surface Pressure(mb)
!
! Start integration
! -----------------
!
      n = 0
   10 continue
      n = n + 1
!
! Pre-processing
! --------------
!
      k = mod(n,12)
      if (k.eq.0) k = 12
      mes = k
      td = tsoil(i,j,k)
      ta = temp(i,j,k)
      pr = prec(i,j,k)
      ipar = par(i,j,k)
!      ae = 2.26457*ta + 67.5876 !available energy (W/m2) [Eq. 8]
      ae = 2.895*ta+52.326                                                                     !Available energy(W/m2)_From NCEP-NCAR Reanalysis Data
!
! Monthly water budget
! ====================
!
      call budget (mes,wini,gini,sini,td,ta,pr,spre,ae,ca,ipar,
     &                 wfim,gfim,sfim,smes,rmes,emes,epmes,
     &                 phmes,armes,nppmes,laimes,
     &                 clmes,csmes,hrmes,cresmes,
     &                 cresmes1,cresmes2,cresmes3,
     &                 nppmes1,nppmes2,nppmes3)
!
!
! Update variables
! ----------------
!
        wsoil(i,j,k) = wfim
        gsoil(i,j,k) = gfim
        ssoil(i,j,k) = sfim
        snowm(i,j,k) = smes
        runom(i,j,k) = rmes
        evapm(i,j,k) = emes
        emaxm(i,j,k) = epmes
        cres(i,j,k)  = cresmes
        cres1(i,j,k) = cresmes1
        cres2(i,j,k) = cresmes2
        cres3(i,j,k) = cresmes3
        lai(i,j,k) = laimes
        photo(i,j,k) = phmes
        aresp(i,j,k) = armes
        npp(i,j,k) = nppmes
        npp1(i,j,k) = nppmes1
        npp2(i,j,k) = nppmes2
        npp3(i,j,k) = nppmes3
        clit(i,j,k) = clmes
        csoil(i,j,k) = csmes
        hresp(i,j,k) = hrmes
        wini = wfim
        gini = gfim
        sini = sfim
!
! Check if equilibrium is attained
! --------------------------------
!
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
!
! Environmental variables
! =======================
!
! Initialize
! ----------
!
      tmin(i,j) = 100.0
      seanpp(i,j) = 0.0
      meanpp(i,j) = 0.0
      meanpp1(i,j) = 0.0
      meanpp2(i,j) = 0.0
      meanpp3(i,j) = 0.0
      mphoto(i,j) = 0.0
      maresp(i,j) = 0.0
      nppmin(i,j) = 100.0
      nppmax(i,j) = -100.0
      meanhr(i,j) = 0.0
      meancs(i,j) = 0.0
      ave_wsoil(i,j) = 0.0
      ave_evap(i,j) = 0.0
      ave_cres(i,j) = 0.0
      ave_cres1(i,j) = 0.0
      ave_cres2(i,j) = 0.0
      ave_cres3(i,j) = 0.0
!
! Calculate tmin, meanpp, seanpp
! ------------------------------
!
      counter = counter + 1
!
      do k=1,12
         if (temp(i,j,k).lt.tmin(i,j)) tmin(i,j) = temp(i,j,k)
         meanpp(i,j) = meanpp(i,j) + (npp(i,j,k)/12)
         meanpp1(i,j) = meanpp1(i,j) + (npp1(i,j,k)/12)
         meanpp2(i,j) = meanpp2(i,j) + (npp2(i,j,k)/12)
         meanpp3(i,j) = meanpp3(i,j) + (npp3(i,j,k)/12)
         mphoto(i,j) = mphoto(i,j) + (photo(i,j,k)/12)
         maresp(i,j) = maresp(i,j) + (aresp(i,j,k)/12)
         meanhr(i,j) = meanhr(i,j) + (hresp(i,j,k)/12)
         meancs(i,j) = meancs(i,j) + (csoil(i,j,k)/12)
         ave_wsoil(i,j) = ave_wsoil(i,j) + (wsoil(i,j,k)/12)
         ave_evap(i,j) = ave_evap(i,j) + (evapm(i,j,k)/12)
         ave_cres(i,j) = ave_cres(i,j) + (cres(i,j,k)/12)
         ave_cres1(i,j) = ave_cres1(i,j) + (cres1(i,j,k)/12)
         ave_cres2(i,j) = ave_cres2(i,j) + (cres2(i,j,k)/12)
         ave_cres3(i,j) = ave_cres3(i,j) + (cres3(i,j,k)/12)
!
         if (npp(i,j,k).lt.nppmin(i,j)) nppmin(i,j) = npp(i,j,k)
         if (npp(i,j,k).gt.nppmax(i,j)) nppmax(i,j) = npp(i,j,k)
!
         if (npp1(i,j,k).lt.nppmin(i,j)) nppmin(i,j) = npp1(i,j,k)
         if (npp1(i,j,k).gt.nppmax(i,j)) nppmax(i,j) = npp1(i,j,k)
!
         if (npp2(i,j,k).lt.nppmin(i,j)) nppmin(i,j) = npp2(i,j,k)
         if (npp2(i,j,k).gt.nppmax(i,j)) nppmax(i,j) = npp2(i,j,k)
!
         if (npp3(i,j,k).lt.nppmin(i,j)) nppmin(i,j) = npp3(i,j,k)
         if (npp3(i,j,k).gt.nppmax(i,j)) nppmax(i,j) = npp3(i,j,k)
!
      enddo
!
      npptotal(i,j) = meanpp1(i,j) + meanpp2(i,j) + meanpp3(i,j)
!
      mnpp1_f(i,j) = ((meanpp1(i,j))/(npptotal(i,j)))
      if (npptotal(i,j).eq.0.and.meanpp1(i,j).eq.0) mnpp1_f(i,j) = 0
!
      mnpp2_f(i,j) = ((meanpp2(i,j)/npptotal(i,j)))
      if (npptotal(i,j).eq.0.and.meanpp2(i,j).eq.0) mnpp2_f(i,j) = 0
!
      mnpp3_f(i,j) = ((meanpp3(i,j)/npptotal(i,j)))
      if (npptotal(i,j).eq.0.and.meanpp3(i,j).eq.0) mnpp3_f(i,j) = 0
!
!      print*,npptotal(i,j),meanpp1(i,j),mnpp1_f(i,j)
      
      one = mnpp1_f(i,j) + mnpp2_f(i,j) + mnpp3_f(i,j)
      if (one.ne.1.) one = 1
!
      mnpp1_pon(i,j) = (mnpp1_f(i,j) * meanpp1(i,j))
      mnpp2_pon(i,j) = (mnpp2_f(i,j) * meanpp2(i,j))
      mnpp3_pon(i,j) = (mnpp3_f(i,j) * meanpp3(i,j))
!
!      if (counter.eq.3197) then                                                                !Manaus = 5245
!      print*,meanpp(i,j)
!      endif
!
!      print*,'npptotal',npptotal(i,j),'meanpp',meanpp(i,j),'meanpp1',
!     & meanpp1(i,j),'meanpp2',meanpp2(i,j),'meanpp3',meanpp3(i,j),
!     & 'mnpp1_f',mnpp1_f(i,j),'mnpp2_f',mnpp2_f(i,j),'mnpp3_f',
!     & mnpp3_f(i,j),'one',one,'mnpp1_pon',mnpp1_pon(i,j),'mnpp2_pon',
!     & mnpp2_pon(i,j),'mnpp3_pon',mnpp3_pon(i,j)
!      endif

      if (meanpp(i,j).gt.0.0) then
      seanpp(i,j) = (nppmax(i,j)-nppmin(i,j))/(meanpp(i,j))
      else
      seanpp(i,j) = (nppmax(i,j)-nppmin(i,j))
      endif
      if(seanpp(i,j).gt.auxs) auxs=seanpp(i,j)	                                       !In order to let Snpp dimensionless
!
      if (meanpp1(i,j).gt.0.0) then
      seanpp(i,j) = (nppmax(i,j) - nppmin(i,j))/(meanpp1(i,j))
      else
      seanpp(i,j) = (nppmax(i,j) - nppmin(i,j))
      endif
!
      if (meanpp2(i,j).gt.0.0) then
      seanpp(i,j) = (nppmax(i,j) - nppmin(i,j))/(meanpp2(i,j))
      else
      seanpp(i,j) = (nppmax(i,j) - nppmin(i,j))
      endif
!
      if (meanpp3(i,j).gt.0.0) then
      seanpp(i,j) = (nppmax(i,j) - nppmin(i,j))/(meanpp3(i,j))
      else
      seanpp(i,j) = (nppmax(i,j) - nppmin(i,j))
      endif
!
! Calculate wsoil2                                                                             !Only for grid points without soil ice
! ----------------
!
      ice = 0
      do k=1,12
         if (gsoil(i,j,k)/wmax.gt.1.e-7) ice = 1
      enddo
      if (ice.eq.0) then
         do k=1,12
            wsoil2(i,j,k) = wsoil(i,j,k)/wmax
         enddo
      endif
!
! Close land "if" and "do" loops
! ------------------------------
!
      endif
!
!     if ((i.eq.20).and.(j.eq.45)) then
!     write(*,22) meanpp(i,j), mphoto(i,j), maresp(i,j)
!  22 format(f6.4)
!     stop
!      endif
!
      enddo
      enddo
!
! Final determination of Snpp                                                                  !Loop again cause of auxs...
! ---------------------------
!
      do i=1,nx
      do j=1,ny
        if (int(lsmk(i,j)).ne.0) then
         seanpp(i,j) = seanpp(i,j)/auxs
        endif
      enddo
      enddo
!
      return
      end
!
!=======================================================================
c23456789
!
      subroutine budget (month,w1,g1,s1,tsoil,temp,prec,p0,ae,ca,                              !input
     &                    ipar,w2,g2,s2,smavg,ruavg,evavg,                                      !output
     &                         epavg,phavg,aravg,nppavg,laiavg,
     &                         clavg,csavg,hravg,cresavg,
     &                   cresavg1,cresavg2,cresavg3,
     &                   nppavg1,nppavg2,nppavg3)
!
! Surface water                                                                                !soil moisture, snow and ice budget for a single month
! =============
!
! Input
! -----
!
!        Month  : actual month (1-12)
!        w1     : initial (previous month last day) soil moisture storage (mm)
!        g1     : initial soil ice storage (mm)
!        s1     : initial overland snow storage (mm)
!        tsoil  : soil temperature (oC)
!        temp   : surface air temperature (oC)
!        prec   : precipitation (mm/day)
!        p0     : surface pressure (mb)
!        ae     : available energy (W/m2)
!        ca     : atmospheric carbon
!
! Output
! ------
!
!        ipar   : incident photosynthetic active radiation
!        w2     : final (last day) soil moisture storage (mm)
!        g2     : final soil ice storage (mm)
!        s2     : final overland snow storage (mm)
!        smavg  : snowmelt monthly average (mm/day)
!        ruavg  : runoff monthly average (mm/day)
!        evavg  : actual evapotranspiration monthly average (mm/day)
!        epavg  : maximum evapotranspiration monthly average (mm/day)
!        phavg  :
!        aravg  :
!        nppavg :
!        laiavg :
!        clavg  :
!        csavg  :
!        hravg  :
!        rcavg  :
!        rcavg1 :
!        rcavg2 :
!
! Input/Output Variables
!
      integer npft
      integer month
      real w1,g1,s1,tsoil,temp,prec,p0,ae,ca,ipar,
     &     w2,g2,s2,smavg,ruavg,evavg,epavg,cresavg,
     &     phavg,aravg,nppavg,laiavg,
     &     clavg,csavg,hravg
      real cresavg1,cresavg2,cresavg3
      real cres_pft1,cres_pft2,cres_pft3
      real nppavg1,nppavg2,nppavg3
      real nppa_pft1,nppa_pft2,nppa_pft3
!
! Internal variables
! ------------------
!
      real rh,wmax,tsnow,tice
      real psnow,prain
      real w,g,s
      real rimelt,smelt,roff,evap,emax
      integer ndmonth(12)
      data ndmonth /31,28,31,30,31,30,31,31,30,31,30,31/                                       !Number of days for each month
!
! Carbon Cycle
! ------------
!
      real ph,ar,nppa,nppb,laia,cl,cs,hr                                                       !Carbon cycle
!
! Parameters
! ----------
!
!     rh    = 0.6                                                                              !Relative humidity (dimensionless)
      rh    = 0.685                                                                            !From NCEP-NCAR Reanalysis data
      wmax  = 500.0                                                                            !Soil moisture availability (mm)
      tsnow = -1.0                                                                             !Temperature threshold for snowfall (oC)
      tice  = -2.5                                                                             !Temperature threshold for soil freezing (oC)
!
! Precipitation
! =============
!
      psnow = 0.0
      prain = 0.0
      if (temp.lt.tsnow) then
        psnow = prec/real(ndmonth(month))                                                        !Snowfall (mm/day)
      else
        prain = prec/real(ndmonth(month))                                                        !Rainfall (mm/day)
      endif
!
! Initialization
! --------------
!
      w = w1 	                                                                       !w = daily soil moisture storage (mm)
      g = g1 	                                                                       !g = daily soil ice storage (mm)
      s = s1 	                                                                       !s = daily overland snow storage (mm)
      smavg = 0.
      ruavg = 0.
      evavg = 0.
      epavg = 0.
      cresavg = 0.
      cresavg1 = 0.
      cresavg2 = 0.
      cresavg3 = 0.
      laiavg = 0.
      phavg = 0.
      aravg = 0.
      nppavg = 0.
      nppavg1 = 0.
      nppavg2 = 0.
      nppavg3 = 0.
      clavg = 0.
      csavg = 0.
      hravg = 0.
!
! Numerical integration
! ---------------------
!
      npft = 3
!
      do i=1,ndmonth(month)
!
      var3 = 0
      var4 = 0
      var5 = 0
      var6 = 0
      var7 = 0
!
      do l=1,3
!
! Carbon cycle (photosynthesis, plant respiration and NPP)
! ========================================================
!
      call productivity1 (temp,p0,w,wmax,ca,ipar,l,                                            !input
     &                    ph,ar,nppa,laia,f5,vm,mgama,
     &                    f2,f3,rmax,r,ci,f10,jc)

      if (l.eq.1) vm_pft1 = vm
      if (l.eq.2) vm_pft2 = vm
      if (l.eq.3) vm_pft3 = vm
!
      if (l.eq.1) nppa_pft1 = nppa
      if (l.eq.2) nppa_pft2 = nppa
      if (l.eq.3) nppa_pft3 = nppa
!
      if (l.eq.1) ph_pft1 = ph
      if (l.eq.2) ph_pft2 = ph
      if (l.eq.3) ph_pft3 = ph
!
      if (l.eq.1) ar_pft1 = ar
      if (l.eq.2) ar_pft2 = ar
      if (l.eq.3) ar_pft3 = ar
!
      if (l.eq.1) laia_pft1 = laia
      if (l.eq.2) laia_pft2 = laia
      if (l.eq.3) laia_pft3 = laia
!
      var3 = var3 + (nppa/npft)
      var4 = var4 + (ph/npft)
      var5 = var5 + (ar/npft)
      var6 = var6 + (laia/npft)
      var7 = var7 + (vm/npft)

!      print*,'temp',temp,'mgama',mgama,'f2',f2,'f3',f3,
!     &       'rmax',rmax,'r',r,'ci',ci,'f10',f10,'jc',jc
!      stop

      enddo
!
      nppa = var3
      ph = var4
      ar = var5
      laia = var6
      vm = var7

!      print*,'temp',temp,'mgama',mgama,'f2',f2,'f3',f3,
!     &       'rmax',rmax,'r',r,'ci',ci,'f10',f10,'jc',jc
!      stop


!      print*,'vm_pft1',vm_pft1,'vm_pft2',vm_pft2,
!     &       'vm_pft3',vm_pft3
!      stop
!
!      print*,temp,mgama
!      stop

!      print*,f2
!      stop
!
!      if (ph.ne.0.) then
!      print*,'media',ph,ph_pft1,ph_pft2,ph_pft3
!      stop
!      endif

!      if (laia.gt.0.25) then
!      print*,'media',laia,laia_pft1,laia_pft2,laia_pft3
!      stop
!      endif

!      if (ar.ne.0.) then
!      print*,'media',ar,ar_pft1,ar_pft2,ar_pft3
!      stop
!      endif

!      if (nppa.ne.0.) then
!      print*,nppa
!      stop
!      endif

!      if (ph.ne.0.) then
!      print*,ph                                        s
!      stop
!      endif

!      if (ar.ne.0.) then
!      print*,ar
!      stop
!      endif
!
!      if(laia.gt.0.25) then
!      print*,laia
!      stop
!      endif

      var2 = 0
!
      do m=1,3
!
! Canopy resistence
! =================

      call canopy_resistence (nppb,nppa,m,                                                     !input
     &                        cres2)                                                           !output
!
      if (m.eq.1) cres_pft1 = cres2
      if (m.eq.2) cres_pft2 = cres2
      if (m.eq.3) cres_pft3 = cres2
!
      var2 = var2 + (cres2/npft)
!
!      print*,'m',m,'npft',npft,'cres2',cres2

!
      enddo
!
      cres2 = var2
!
!      print*,'media',cres2,cres_pft1,cres_pft2,cres_pft3
!
! Maximum evapotranspiration (emax)
! =================================
!
      call evpot2 (p0,temp,rh,ae,emax)
!
! Snow budget
! ===========
!
      smelt = 2.63 + 2.55*temp + 0.0912*temp*prain                                       !Snowmelt (mm/day)
      smelt = amax1(smelt,0.)
      smelt = amin1(smelt,s+psnow)
      ds = psnow - smelt
      s = s + ds
!
! Water budget
! ============
!
      if (tsoil.le.tice) then                                                                  !Frozen soil
        g = g + w                                                                                !Soil moisture freezes
        w = 0.0
        roff = smelt + prain
        evap = 0.0
        ph = 0.0
        ar = 0.0
        nppa = 0.0                                                                                                        !!!!!!!!!!!!
        laia = 0.0
        cl = 0.0
        cs = 0.0
        hr = 0.0
        cres2 = 100.0                                                                            !Default value, equal to aerodynamic resistance (below)
      else                                                                                     !Non-frozen soil
        w = w + g                                                                                !Soil ice melts
        g = 0.0
        rimelt = 0.0
        if (w.gt.wmax) then
          rimelt = w - wmax                                                                        !Runoff due to soil ice melting
          w = wmax
        endif
!
        call runoff (w,wmax,roff)                                                                !soil moisture runoff (roff, mm/day)
        call penman (p0,temp,w,wmax,rh,ae,cres2,evap)                                              !actual evapotranspiration (evap, mm/day)
        dw = prain + smelt - evap - roff
        w = w + dw
        if (w.gt.wmax) then
          roff = roff + (w - wmax)
          w = wmax
        endif
        if (w.lt.0.) w = 0.
        roff = roff + rimelt                                                                     !total runoff
!
! Carbon cycle (Microbial respiration, litter and soil carbon)
! ============================================================
!
        call carbon2 (tsoil,f5,evap,laia,                                                        !input
     &                cl,cs,hr)                                                                  !output
      endif
!
! Updating monthly values
! -----------------------
!
      smavg = smavg + smelt
      ruavg = ruavg + roff
      evavg = evavg + evap
      epavg = epavg + emax                                                                                              !rever esse valor
      cresavg = cresavg + cres2
      cresavg1 = cresavg1 + cres_pft1
      cresavg2 = cresavg2 + cres_pft2
      cresavg3 = cresavg3 + cres_pft3
      phavg = phavg + ph/365.0                                                                 !kgC/m2
      aravg = aravg + ar/365.0                                                                 !kgC/m2
      nppavg = nppavg + nppa/365.0                                                             !kgC/m2
      nppavg1 = nppavg1 + nppa_pft1/365.0                                                             !kgC/m2
      nppavg2 = nppavg2 + nppa_pft2/365.0                                                             !kgC/m2
      nppavg3 = nppavg3 + nppa_pft3/365.0                                                             !kgC/m2
      laiavg = laiavg + laia/365.0
      clavg = clavg + cl/365.0
      csavg = csavg + cs/365.0
      hravg = hravg + hr/365.0                                                                 !kgC/m2

      enddo
!
! Final calculations
! ------------------
!
      w2 = w
      g2 = g
      s2 = s
      smavg = smavg/real(ndmonth(month))
      ruavg = ruavg/real(ndmonth(month))
      evavg = evavg/real(ndmonth(month))
      epavg = epavg/real(ndmonth(month))
      cresavg = cresavg/real(ndmonth(month))
      cresavg1 = cresavg1/real(ndmonth(month))
      cresavg2 = cresavg2/real(ndmonth(month))
      cresavg3 = cresavg3/real(ndmonth(month))
      phavg = phavg*12.0                                                                     !kgC/m2/yr
      aravg = aravg*12.0                                                                     !kgC/m2/yr
      nppavg = nppavg*12.0                                                                   !kgC/m2/yr
      nppavg1 = nppavg1 * 12.0
      nppavg2 = nppavg2 * 12.0
      nppavg3 = nppavg3 * 12.0
      laiavg = laiavg*12.0
      clavg = clavg*12.0                                                                     !kgC/m2
      csavg = csavg*12.0                                                                     !kgC/m2
      hravg = hravg*12.0                                                                     !kgC/m2/yr
!
      return
      end
!
! =========================================================
      subroutine penman (spre,temp,w,wmax,ur,rn,cres2,evap)
! =========================================================
!
! Inputs
! ------
!
!       spre   = surface pressure                                                              !mb
!       temp   = temperature                                                                   !oC
!       w      = saturation                                                                    !0-1,dimensionless
!       ur     = relative humidity                                                             !0-1,dimensionless
!       rn     = radiation balance                                                             !W m-2
!       cres2  = canopy resistence                                                             !s/m
!
! Output
! ------
!
!        evap  = evapotranspiration                                                            !mm/day
!
      real spre,temp,w,wmax,ur,rn,cres2,evap
!
! Parameters
! ----------
!
      ra =    100.                                                      !s/m
      h5    = 0.0275                                                    !mb-1
!
! Delta
! -----
!
      t1 = temp + 1.
      t2 = temp - 1.
      call tetens(t1,es1)
      call tetens(t2,es2)
      delta = (es1-es2)/(t1-t2)                                                            !mb/oC
!
! Delta_e
! -------
!
      call tetens (temp,es)
      delta_e = es*(1. - ur)                                                                  !mb
!
      if ((delta_e.ge.(1./h5)-0.5).or.(cres2.ge.4500)) evap = 0.                                                      !rever
      if ((delta_e.lt.(1./h5)-0.5).or.(cres2.lt.4500)) then                                                           !rever
!
! Gama and gama2
! --------------
!
        gama  = spre*(1004.)/(2.45e6*0.622)
        gama2 = gama*(ra + cres2)/ra
!
! Real evapotranspiration
! -----------------------
!
        evap = (delta*rn + (1.20*1004./ra)*delta_e)/(delta+gama2)                        !W/m2
        evap = evap*(86400./2.45e6)                                                            !mm/day
        evap = amax1(evap,0.)                                                                    !eliminates condensation
      endif
!
      return
      end
!
!=============================================
      subroutine evpot2 (spre,temp,ur,rn,evap)
!=============================================
!
! Inputs
! ------
!
!       spre   = surface pressure                                                              !mb
!       temp   = temperature                                                                   !oC
!       ur     = relative humidity                                                             !0-1,dimensionless
!       rn     = radiation balance                                                             !W m-2
!
! Outputs
! -------
!
!        evap  = potencial evapotranspiration without stress                                   !mm/day
!
      real spre,temp,ur,rn,evap
!
! Parameters
! ----------
!
      ra =    100.                                                      !s/m
      cresmin = 100.                                                      !s/m
!
! Delta
! -----
!
      t1 = temp + 1.
      t2 = temp - 1.
      call tetens(t1,es1)
      call tetens(t2,es2)
      delta = (es1-es2)/(t1-t2)                                                            !mb/oC
!
! Delta_e
! -------
!
      call tetens (temp,es)
      delta_e = es*(1. - ur)                                                                  !mb
!
! Stomatal Conductance
! --------------------
!
      cres = cresmin
!
! Gama and gama2
! --------------
!
      gama  = spre*(1004.)/(2.45e6*0.622)
      gama2 = gama*(ra + cres)/ra
!
! Potencial evapotranspiration(without stress)
! --------------------------------------------
!
      evap = (delta*rn + (1.20*1004./ra)*delta_e)/(delta+gama2)                          !W/m2
      evap = evap*(86400./2.45e6)                                                            !mm/day
      evap = amax1(evap,0.)                                                                    !eliminates condensation
!
      return
      end
!
!=====================================================================
!
      subroutine runoff (w,wmax,roff)
      real w,roff
!     roff = 38.*((w/wmax)**11.)                                        ![Eq. 10]
      roff = 11.5*((w/wmax)**6.6)                                       !from NCEP-NCAR Reanalysis data
      return
      end
!
!=====================================================================
!
      subroutine tetens (t,es)
      real t,es
      if (t.ge.0.) then
      es = 6.1078*exp((7.5*t/(237.3+t))*log(10.))
      else
      es = 6.1078*exp((9.5*t/(265.5+t))*log(10.))
      endif
      return
      end

