     subroutine budget (pft,month,w1,g1,s1,tsoil,temp,prec,p0,ae,ca,ipar,  !input
     &     cl1_pft,ca1_pft,cf1_pft,
     &     w2,g2,s2,cl2_pft,ca2_pft,cf2_pft,smavg,ruavg,
     &     evavg,epavg,phavg,aravg,nppavg,laiavg, !output
     &     rmlavg,rmfavg,rmsavg,rmavg, !output
     &     rglavg,rgfavg,rgsavg,rgavg, !output
     &     clavg,csavg,hravg,rcavg,
     &     cleafavg_pft,cawoodavg_pft,cfrootavg_pft) !output
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
c        cleaf1: initial (previous month last day) carbon content on leaf compartment (kgC/m2)
c output w2    : final (last day) soil moisture storage (mm)
c        g2    : final soil ice storage (mm)
c        s2    : final overland snow storage (mm)
c        smavg : snowmelt monthly average (mm/day)
c        ruavg : runoff monthly average (mm/day)
c        evavg : actual evapotranspiration monthly average (mm/day)
c        epavg : maximum evapotranspiration monthly average (mm/day)
c        cleaf2: month final carbon content on leaf compartment (kgC/m2)
c        cleafavg: average leaf carbon content (kgC/m2)
c
c=======================================================================
c
c i/o variables
      !  MO MAS parameter (npft=3,npls=3)!npls:number of pfts


      !inputs
      integer month, pft
      real w1,g1,s1,tsoil,temp,prec,p0,ae,ca,ipar
      real cleaf1_pft ,cawood1_pft ,cfroot1_pft

      real    w2,g2,s2,smavg,ruavg,evavg,epavg,rcavg,
     &     phavg,aravg,nppavg,laiavg,
     &     clavg,csavg,hravg, cleafavg_pft,cawoodavg_pft,
     &     cfrootavg_pft

      real cl2_pft,ca2_pft,cf2_pft

c internal variables
      real rh,wmax,tsnow,tice
      real psnow,prain
      real w,g,s
      real rimelt,smelt,roff,evap,emax
      integer ndmonth(12) !number of days for each month
      data ndmonth /31,28,31,30,31,30,31,31,30,31,30,31/
c carbon cycle
      real ph,ar,nppa,nppb,laia,cl,cs,hr,
     &     rm,rml,rmf,rms,rg,rgl,rgf,rgs,
     &     rmlavg,rmfavg,rmsavg,rmavg,
     &     rglavg,rgfavg,rgsavg,rgavg
c carbon allocation

      real cl1,cl2,ca1,ca2,cf1,cf2

c
c parameters
c      rh    = 0.6   !relative humidity (adimensional)
      rh    = 0.685 !from NCEP-NCAR Reanalysis data
      wmax  = 500.0 !soil moisture availability (mm)
      tsnow = -1.0  !temperature threshold for snowfall (oC)
      tice  = -2.5  !temperature threshold for soil freezing (oC)
c
!     
!     Precipitation
!     =============
!     
      psnow = 0.0
      prain = 0.0
      if (temp.lt.tsnow) then
         psnow = prec/real(ndmonth(month))
      else
         prain = prec/real(ndmonth(month))
      endif
!     
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


c     
c     
      
c     numerical integration
      do i=1,ndmonth(month)
         
         
         
         cl1 = cl2_pft
         ca1 = ca2_pft
         cf1 = cf2_pft
         
         
         beta_leaf = alfa_leaf
         beta_awood = alfa_awood
         beta_froot = alfa_froot
         
         npp_init = cl1 + ca1 + cf1
         nppb = amax1(npp_init,0.05)
         rc2 = (ca/(9.0*(nppb*2.64e-5)*0.685*(p0*100)))
         
         
         if ((i.eq.1).and.(month.eq.1)) then
            cl1 = cleaf1_pft
            ca1 = cawood1_pft
            cf1 = cfroot1_pft
            
            beta_leaf=0.
            beta_awood=0.
            beta_froot=0.
            
            npp_init = cl1 + ca1 + cf1
            nppb = amax1(npp_init,0.05)
            rc2 = (ca/(9.0*(nppb*2.64e-5)*0.685*(p0*100)))
         endif
         
!     TEM UM PROBLRMA GRAVE COM A VATIAVEL RC2 NA CHAMADE DE CARBON 1/ ELA NAO EXISTE
!     DENTRO DE CARBON1 E POR ISSO ACONTECE UMA POSSIVEL DIVISAO POR 0 DENTRO DA CARBON1
!     vamos resolver este problema calculando rc2 com a npp inicial
         
         
         
c     carbon cycle (photosynthesis, plant respiration and NPP)
         
         
         call carbon1 (pft, temp,p0,w,wmax,ca,ipar,tsoil,emax,rc2, !input
     &        cl1,ca1,cf1,
     &        beta_leaf,beta_awood,beta_froot, !input
     &        ph,ar,nppa,laia,f5,rm,rml,rmf,rms,rg,rgl,rgf,rgs) !output
         
         
         
!     carbon allocation (carbon content on each compartment)
         call allocation (pft, nppa,cl1,ca1,cf1, !input
     &        cl2,ca2,cf2)      !output
         
         alfa_leaf  = cl2 - cl1 
         alfa_awood = ca2 - ca1 
         alfa_froot = cf2 - cf1 
         
         cl2_pft = cl2
         ca2_pft = ca2
         cf2_pft = cf2
         
         
c     
c     maximum evapotranspiration (emax)
         call evpot2 (p0,temp,rh,ae,emax)
c     
c     snow budget
         smelt = 2.63 + 2.55*temp + 0.0912*temp*prain !snowmelt (mm/day) [Eq. 4]
         smelt = amax1(smelt,0.)
         smelt = amin1(smelt,s+psnow)
         ds = psnow - smelt     ![Eq. 2]
         s = s + ds
c     
c     water budget
         if (tsoil.le.tice) then !frozen soil
            g = g + w           !soil moisture freezes
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
            rc2 = 100.0         !default value, equal to aerodynamic resistance (below)
         else                   !non-frozen soil
            w = w + g           !soil ice melts
            g = 0.0
            rimelt = 0.0
            if (w.gt.wmax) then
               rimelt = w - wmax !runoff due to soil ice melting
               w = wmax
            endif
c     
c     Canopy resistance (based in Sellers et al. 1996; SiB2)
c     (rc2 ; s/m) [Eq. 32]
c     [NPP*2.64e-6 converts kgC/m2/yr to molCO2/m2/s]
c     [p0*100 convertes hPa (mb) to Pa]
c     
            nppb = amax1(nppa,0.05)
            rc2 = (ca/(9.0*(nppb*2.64e-5)*0.685*(p0*100)))
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
c     carbon cycle (Microbial respiration, litter and soil carbon)
            call carbon2 (tsoil,f5,evap,laia, !input
     &           cl,cs,hr)      !output
         endif
c     
c     updating monthly values
         smavg = smavg + smelt
         ruavg = ruavg + roff
         evavg = evavg + evap
         epavg = epavg + emax
         rcavg = rcavg + rc2
         phavg = phavg + ph/365.0 !kgC/m2
         aravg = aravg + ar/365.0 !kgC/m2
         rmlavg = rmlavg + rml/365.0
         rmfavg = rmfavg + rmf/365.0
         rmsavg = rmsavg + rms/365.0
         rmavg = rmavg + rm/365.0
         rglavg = rglavg + rgl/365.0
         rgfavg = rgfavg + rgf/365.0
         rgsavg = rgsavg + rgs/365.0
         rgavg = rgavg + rg/365.0
         nppavg = nppavg + nppa/365.0 !kgC/m2
         laiavg = laiavg + laia/365.0
         clavg = clavg + cl/365.0
         csavg = csavg + cs/365.0
         hravg = hravg + hr/365.0 !kgC/m2
         
         cleafavg_pft = cleafavg_pft + cl2_pft
         cawoodavg_pft = cawoodavg_pft + ca2_pft
         cfrootavg_pft = cfrootavg_pft + cf2_pft
c     
         
      enddo
c     
c     final calculations
      
      w2 = w
      g2 = g
      s2 = s
      
      smavg = smavg/real(ndmonth(month))
      ruavg = ruavg/real(ndmonth(month))
      evavg = evavg/real(ndmonth(month))
      epavg = epavg/real(ndmonth(month))
      rcavg = rcavg/real(ndmonth(month))
      phavg = phavg*12.0        !kgC/m2/yr
      aravg = aravg*12.0        !kgC/m2/yr
      nppavg = nppavg*12.0      !kgC/m2/yr
      laiavg = laiavg*12.0
      clavg = clavg*12.0        !kgC/m2
      csavg = csavg*12.0        !kgC/m2
      hravg = hravg*12.0        !kgC/m2/yr

      rmlavg = rmlavg*12.0
      rmfavg = rmfavg*12.0
      rmsavg = rmsavg*12.0
      rmavg = rmavg*12.0
      rglavg = rglavg*12.0
      rgfavg = rgfavg*12.0
      rgsavg = rgsavg*12.0
      rgavg = rgavg*12.0

      cleafavg_pft = cleafavg_pft/real(ndmonth(month))
      cawoodavg_pft = cawoodavg_pft/real(ndmonth(month))
      cfrootavg_pft = cfrootavg_pft/real(ndmonth(month))
      
      return
      end

