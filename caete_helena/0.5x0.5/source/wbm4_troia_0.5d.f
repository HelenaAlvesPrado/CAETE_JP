c234567
!
      subroutine wbm (prec,temp,lsmk,p0,ca,par,        		                             
     &                tmin,meanpp,seanpp,mphoto,maresp,           	                       
     &                meanhr,meancs,wsoil2,evapm,npp,
     &                photo,aresp,rcm,ave_wsoil,ave_evap,ave_rc,
     &                npp1,npp2,npp3,meanpp1,meanpp2,meanpp3)
!
! ===========================================================================================
! Water balance model (WBM).
! From monthly climatologies of precipitation and surface temperature, the WBM calculates the
! environmental variables.
!
! 05 Jul 2005, MDO: ae, rh & runoff are changed
! 11 Jul 2005, MDO: wsoil2 is written (for testing purpose)
! 31 Ago 2006, DML: carbon cycle is included
! ===========================================================================================
!
! Variables
! =========
!
      integer,parameter :: nx=720,ny=360
      real,parameter :: no_data = -9999.0
      integer counter,i,j,k      
!
! Inputs
! ------
!      
      real ca
      real p0(nx,ny,12)                                                                          
      real lsmk(nx,ny)                                                                        
      real prec(nx,ny,12)                                                                      
      real temp(nx,ny,12)                                                                      
      real par(nx,ny,12) 
!
! Outputs
! -------
!            
      real tmin(nx,ny)                                                          
      real meancs(nx,ny)                                                              
      real ave_wsoil(nx,ny)                                                               
      real ave_evap(nx,ny)                                                                     
      real ave_rc(nx,ny)  
      real meanpp(nx,ny)   
      real meanpp1(nx,ny)
      real meanpp2(nx,ny)
      real meanpp3(nx,ny)
      real seanpp(nx,ny) 
      real maresp(nx,ny)
      real meanhr(nx,ny)
      real mphoto(nx,ny)
!
      real aresp(nx,ny,12)                                                         
      real photo(nx,ny,12)                                                             
      real wsoil2(nx,ny,12)                                                                   
      real npp(nx,ny,12)  
      real npp1(nx,ny,12)
      real npp2(nx,ny,12)
      real npp3(nx,ny,12)
      real lai(nx,ny,12)
      real clit(nx,ny,12)
      real csoil(nx,ny,12)                                                               
      real hresp(nx,ny,12)                                                                
      real rcm(nx,ny,12)
!
! Internal Variables
! ------------------
!      
      integer mes
      integer kk
      real H
      real diffu
      real tau
      real t0,t1
      real laimes
      real ipar
      real ae,auxs,dwww,nerro,wmax
      real sini,gfim,sfim,gini,wfim,wini
      real pr,ice,spre,ta,td
      real rmes,phmes,smes,rcmes
      real armes,clmes,csmes,emes,epmes,hrmes
!
      real tsoil(nx,ny,12)
      real wsoil(nx,ny,12)
      real gsoil(nx,ny,12)
      real ssoil(nx,ny,12)
      real snowm(nx,ny,12)
      real runom(nx,ny,12)
      real evapm(nx,ny,12)
      real emaxm(nx,ny,12)
      real wg0(nx,ny,12)
!
! Carbon Cycle      
! ------------
!
      real nppmes,nppmes1,nppmes2,nppmes3
      real nppmin(nx,ny)
      real nppmax(nx,ny)
      real meant(nx,ny)
      real seat(nx,ny)
      real mgama,jc,jl,je,jp,f4sun,f4shade
!
! Counting
! --------
!
      counter = 0
!
! Initialize productivity
! -----------------------
!
      mgama   = 0.
      jc      = 0.
      jl      = 0.
      je      = 0.
      jp      = 0.
      f4sun   = 0.
      f4shade = 0.
!
! Soil temperature
! ================
!
      H     = 1.0 
      diffu = 4.e-7*(30.*86400.0)
      tau   = (H**2)/(2.0*diffu)
      auxs  = -100.0
!
! For all grid points
! -------------------
c234567
!
      do i=1,nx
         do j=1,ny
!
! Initialize soil temperature
! ---------------------------
!
            do k=1,12
               tsoil(i,j,k) = no_data                                                         
            enddo
!
! Only for land grid points
! -------------------------
!
      if (int(lsmk(i,j)).ne.0) then
      t0 = 0.
      do n = 1,1200
         k = mod(n,12)
         if (k.eq.0) k = 12
         t1 = t0*exp(-1.0/tau)+(1.0-exp(-1.0/tau))*temp(i,j,k) 
         tsoil(i,j,k) = (t0 + t1)/2.0
         t0 = t1
      enddo
      endif
!      
      enddo
      enddo
!
c234567
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
     &   write(*,*) 'working:',i
!
! Initialize variables
! --------------------
!
      tmin(i,j)      = no_data
      seanpp(i,j)    = no_data
      meanpp(i,j)    = no_data
      meanpp1(i,j)   = no_data
      meanpp2(i,j)   = no_data
      meanpp3(i,j)   = no_data
      mphoto(i,j)    = no_data                                                          
      maresp(i,j)    = no_data                                                           
      meanhr(i,j)    = no_data                                                           
      meancs(i,j)    = no_data                                                          
      ave_wsoil(i,j) = no_data                                                       
      ave_evap(i,j)  = no_data                                                        
      ave_rc(i,j)    = no_data                                                       
!      
      do k=1,12
!
         wsoil(i,j,k)  = no_data                                                  
         gsoil(i,j,k)  = no_data                                                     
         ssoil(i,j,k)  = no_data                                          
         snowm(i,j,k)  = no_data                                                      
         runom(i,j,k)  = no_data                                                     
         evapm(i,j,k)  = no_data                                                    
         emaxm(i,j,k)  = no_data                                              
         wg0(i,j,k)    = no_data                                                         
         wsoil2(i,j,k) = no_data                                          
         rcm(i,j,k)    = no_data                                                        
         lai(i,j,k)    = no_data                                                         
         photo(i,j,k)  = no_data                                                   
         aresp(i,j,k)  = no_data                                                                 
         npp(i,j,k)    = no_data
         npp1(i,j,k)   = no_data                                                         
         npp2(i,j,k)   = no_data                                                         
         npp3(i,j,k)   = no_data                                                       
         clit(i,j,k)   = no_data                                                        
         csoil(i,j,k)  = no_data                                                       
         hresp(i,j,k)  = no_data                                                         
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
      wini  = 0.01                                              
      gini  = 0.0                                                       
      sini  = 0.0                                                                  
!
! Initialization
! --------------
!
      do k=1,12
      wg0(i,j,k) = -1.0     
      spre = p0(i,j,k) * 0.01  ! getting p0 as (nx,ny,12) and converting from Pa to mbar                                                      
      enddo
                                                                        
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
      ae = 2.895*ta+52.326                                                                                                                                
!
! Monthly water budget
! ====================
!
      call budget (mes,wini,gini,sini,td,ta,pr,spre,ae,ca,ipar,
     &                 wfim,gfim,sfim,smes,rmes,emes,epmes,
     &                 phmes,armes,nppmes,laimes,
     &                 clmes,csmes,hrmes,rcmes,
     &                 nppmes1,nppmes2,nppmes3) 
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
      rcm(i,j,k)   = rcmes
      lai(i,j,k)   = laimes
      photo(i,j,k) = phmes
      aresp(i,j,k) = armes
      npp(i,j,k)   = nppmes
      npp1(i,j,k)  = nppmes1
      npp2(i,j,k)  = nppmes2
      npp3(i,j,k)  = nppmes3
      clit(i,j,k)  = clmes
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
      tmin(i,j)      = 100.0
      seanpp(i,j)    = 0.0
      meanpp(i,j)    = 0.0
      meanpp1(i,j)   = 0.0
      meanpp2(i,j)   = 0.0
      meanpp3(i,j)   = 0.0
      mphoto(i,j)    = 0.0
      maresp(i,j)    = 0.0
      nppmin(i,j)    = 100.0
      nppmax(i,j)    = -100.0
      meanhr(i,j)    = 0.0
      meancs(i,j)    = 0.0
      ave_wsoil(i,j) = 0.0
      ave_evap(i,j)  = 0.0
      ave_rc(i,j)    = 0.0
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
         ave_rc(i,j) = ave_rc(i,j) + (rcm(i,j,k)/12)
!      
         if (npp(i,j,k).lt.nppmin(i,j)) nppmin(i,j) = npp(i,j,k)
         if (npp(i,j,k).gt.nppmax(i,j)) nppmax(i,j) = npp(i,j,k)
  
         if (npp1(i,j,k).lt.nppmin(i,j)) nppmin(i,j) = npp1(i,j,k)
         if (npp1(i,j,k).gt.nppmax(i,j)) nppmax(i,j) = npp1(i,j,k)

         if (npp2(i,j,k).lt.nppmin(i,j)) nppmin(i,j) = npp2(i,j,k)
         if (npp2(i,j,k).gt.nppmax(i,j)) nppmax(i,j) = npp2(i,j,k)
!
         if (npp3(i,j,k).lt.nppmin(i,j)) nppmin(i,j) = npp3(i,j,k)
         if (npp3(i,j,k).gt.nppmax(i,j)) nppmax(i,j) = npp3(i,j,k)
!    
      enddo

!      if (counter.eq.3197.or.counter.eq.5245.or.counter.eq.2000.or.
!     &    counter.eq.5500) then
!      print*,meanpp(i,j),meanpp1(i,j),meanpp2(i,j),
!     &       meanpp3(i,j),ave_rc(i,j),'caete','g0=0.2'
!      endif

!      if (counter.eq.5245) then
!      print*,meanpp(i,j),ave_rc(i,j),'caete'
!      endif
!
      if (meanpp(i,j).gt.0.0) then
      seanpp(i,j) = (nppmax(i,j)-nppmin(i,j))/(meanpp(i,j))
      else
      seanpp(i,j) = (nppmax(i,j)-nppmin(i,j))
      endif
      if (seanpp(i,j).gt.auxs) auxs = seanpp(i,j)	                                    
!
! Calculate wsoil2                                                                         
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
! Final determination of Snpp                                                          
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
c234567
!
      subroutine budget (month,w1,g1,s1,tsoil,temp,prec,p0,ae,ca,ipar,
     &                         w2,g2,s2,smavg,ruavg,evavg,                                    
     &                         epavg,phavg,aravg,nppavg,laiavg,
     &                         clavg,csavg,hravg,rcavg,
     &                         nppavg1,nppavg2,nppavg3)
!
! Surface water                                                                            
! =============
!
      integer month
      integer npft
      real w1
      real g1
      real s1
      real tsoil
      real temp
      real prec
      real p0
      real ae
      real ca
      real ipar
!
! Output
! ------
!      
      real w2
      real g2
      real s2
      real smavg
      real ruavg
      real evavg
      real epavg
      real phavg
      real aravg
      real nppavg
      real nppavg1
      real nppavg2
      real nppavg3
      real laiavg
      real clavg
      real csavg
      real hravg
      real rcavg
!
! Internal Variables
!
      real var2,var3,var4,var5,var6,var7
      real nppa_pft1,nppa_pft2,nppa_pft3
      real ph_pft1,ph_pft2,ph_pft3                                                             !Auxiliars
      real ar_pft1,ar_pft2,ar_pft3                                                             !Auxiliars  
      real laia_pft1,laia_pft2,laia_pft3                                                       !Auxiliars
      real f5_pft1,f5_pft2,f5_pft3                                                             !Auxiliars
      real f1_pft1,f1_pft2,f1_pft3                                                             !Auxiliars
      real rh
      real wmax
      real tsnow
      real tice
      real psnow
      real prain
      real w
      real g
      real s
      real rimelt
      real smelt
      real roff
      real evap
      real emax
      integer ndmonth(12)                                    
      data ndmonth /31,28,31,30,31,30,31,31,30,31,30,31/

! Carbon Cycle
!
      real ph
      real ar
      real nppa,nppb
      real laia
      real cl
      real cs
      real hr
      real mgama
      real jc,je,jl,jp
      real f4sun,f4shade
      real gs
      real rc2
      real aa                                                      
      real f1
      real f1a
!
      gs = 0.0
      rc2 = 0.0
      aa = 0.0
      f1 = 0.0
      f1a = 0.0
      
!
! Parameters
! ----------
!
      rh    = 0.685                                                                            
      wmax  = 500.0                                                                            
      tsnow = -1.0                                                                             
      tice  = -2.5                                                                             
!
! Precipitation
! =============
!
      psnow = 0.0
      prain = 0.0
      if (temp.lt.tsnow) then
        psnow = prec/real(ndmonth(month))                                                      
      else
        prain = prec/real(ndmonth(month))                                                      
      endif
!
! Initialization
! --------------
!
      w = w1 	                                                                             
      g = g1 	                                                                             
      s = s1 	                                                                             
      smavg   = 0.
      ruavg   = 0.
      evavg   = 0.
      epavg   = 0.
      rcavg   = 0.
      laiavg  = 0.
      phavg   = 0.
      aravg   = 0.
      nppavg  = 0.
      nppavg1 = 0.
      nppavg2 = 0.
      nppavg3 = 0.
      clavg   = 0.
      csavg   = 0.
      hravg   = 0.
!
! Numerical integration
! ---------------------
!
      npft = 3
!      
      do i=1,ndmonth(month)
!
      var2 = 0.0        !nppa
      var3 = 0.0        !ph    
      var4 = 0.0        !ar
      var5 = 0.0        !laia
      var6 = 0.0        !f5
      var7 = 0.0        !f1

!
      do l=1,3
!      
! Carbon cycle (photosynthesis, plant respiration and NPP)
! ========================================================
!
      call productivity1 (temp,p0,w,wmax,ca,ipar,l,
     &                    ph,ar,nppa,laia,f5,f1)                                                 
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
      if (l.eq.1) f5_pft1 = f5
      if (l.eq.2) f5_pft1 = f5
      if (l.eq.3) f5_pft1 = f5
!
      if (l.eq.1) f1_pft1 = f1
      if (l.eq.2) f1_pft2 = f1
      if (l.eq.3) f1_pft3 = f1      
!
      var2 = var2 + (nppa/npft)
      var3 = var3 + (ph/npft)
      var4 = var4 + (ar/npft)
      var5 = var5 + (laia/npft)
      var6 = var6 + (f5/npft)
      var7 = var7 + (f1/npft)
!
!      print*,'l',l,'npft',npft,'nppa',nppa      
      enddo
!
      nppa = var2
      ph   = var3
      ar   = var4
      laia = var5
      f5   = var6
      f1   = var7

!
!      print*,'media',nppa,nppa_pft1,nppa_pft2,nppa_pft3
!
! Maximum evapotranspiration (emax)
! =================================
!
      call evpot2 (p0,temp,rh,ae,emax)  

! Snow budget
! ===========
!
      smelt = 2.63 + 2.55*temp + 0.0912*temp*prain                                         
      smelt = amax1(smelt,0.)
      smelt = amin1(smelt,s+psnow)
      ds = psnow - smelt
      s = s + ds
!
! Water budget
! ============
!
      if (tsoil.le.tice) then                                                              
        g = g + w                                                                         
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
!        rc2 = 100.0                                                                           
      else                                                                                     
        w = w + g                                                                             
        g = 0.0
        rimelt = 0.0
        if (w.gt.wmax) then
          rimelt = w - wmax                                                          
          w = wmax
        endif
!
! Canopy resistance (based in Sellers et al. 1996; SiB2)
! (rc2 ; s/m) [Eq. 32]

!      if (f1.lt.10e-7) f1 = 0.0
!      f1a = (f1*10e5)
      f1a = (f1*10e5)
      if (f1a.le.0.0) f1a = 0.0
!
      g0 = -0.03
      g1 = 4.9
      D = sqrt(2.2)
      aa = (f1a/350)
!      bb = f1a/(350*D)
!
      gs = g0 + 1.6 * (1+(g1/D)) * (aa)    !molCO2/m2/s
      if (gs.lt.-0.03) gs = -0.03
!
      gs2 = gs/41                                            !m/s
      if (gs2.lt.-0.00070) gs2 = -0.00070
!
      rc2 = gs2
        call runoff (w,wmax,roff)                                                           
        call penman (p0,temp,w,wmax,rh,ae,rc2,evap)                                       
        dw = prain + smelt - evap - roff                                                        
        w = w + dw
        if (w.gt.wmax) then
          roff = roff + (w - wmax)
          w = wmax
        endif
        if (w.lt.0.) w = 0.
        roff = roff + rimelt                                                            
!
! Carbon cycle (Microbial respiration, litter and soil carbon)
! ============================================================
!
        call carbon2 (tsoil,f5,evap,laia,                                                  
     &                cl,cs,hr)                                                           
      endif
!
! Updating monthly values
! =======================
!
      smavg = smavg + smelt                                                                    
      ruavg = ruavg + roff                                                                     
      evavg = evavg + evap                                                                    
      epavg = epavg + emax                                                                   
      rcavg = rcavg + rc2                                                                   
      phavg = phavg + ph/365.0                                                              
      aravg = aravg + ar/365.0                                                              
      nppavg = nppavg + nppa/365.0 
      nppavg1 = nppavg1 + nppa_pft1/365.0                                                   
      nppavg2 = nppavg2 + nppa_pft2/365.0                                                     
      nppavg3 = nppavg3 + nppa_pft3/365.0                                                         
      laiavg = laiavg + laia/365.0                                                         
      clavg = clavg + cl/365.0                                                                
      csavg = csavg + cs/365.0                                                             
      hravg = hravg + hr/365.0                                                             
!    
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
      rcavg = rcavg/real(ndmonth(month))
      phavg = phavg*12.0                                                                       
      aravg = aravg*12.0                                                                      
      nppavg = nppavg*12.0                                                                   
      nppavg1 = nppavg1 * 12.0
      nppavg2 = nppavg2 * 12.0
      nppavg3 = nppavg3 * 12.0
      laiavg = laiavg*12.0
      clavg = clavg*12.0                                                                
      csavg = csavg*12.0                                                                    
      hravg = hravg*12.0                                                                    
!
      return
      end
!
! =========================================================
!
      subroutine penman (spre,temp,w,wmax,ur,rn,rc2,evap)                                                                                        
!
! Inputs
! ------
!
      real spre                                                                                !Surface pressure (mb)
      real temp                                                                                !Temperature (oC)
      real w                                                                                   !Saturation (0-1,dimensionless)
      real ur                                                                                  !Relative humidity (0-1,dimensionless)
      real rn                                                                                  !Radiation balance (W/m2)
      real rc2                                                                                 !Canopy resistence (s/m)
!
! Output
! ------
!
      real evap
!        
! Parameters
! ----------
!
      ra = 100                                                                       
      h5 = 0.0275                                                                          
!
! Delta
! -----
!
      t1 = temp + 1.
      t2 = temp - 1.
      call tetens(t1,es1)                                                                       
      call tetens(t2,es2)                                                           
      delta = (es1-es2)/(t1-t2)                                                        
!
! Delta_e
! -------
!
      call tetens (temp,es)
      delta_e = es*(1. - ur)                                                             
!
      if ((delta_e.ge.(1./h5)-0.5).or.(rc2.ge.4500)) evap = 0.                               
      if ((delta_e.lt.(1./h5)-0.5).or.(rc2.lt.4500)) then                                    
!
! Gama and gama2
! --------------
!
        gama  = spre*(1004.)/(2.45e6*0.622)
        gama2 = gama*(ra + rc2)/ra
!
! Real evapotranspiration
! -----------------------
!
        evap = (delta*rn + (1.20*1004./ra)*delta_e)/(delta+gama2)                           
        evap = evap*(86400./2.45e6)                                                    
        evap = amax1(evap,0.)                                                                  
      endif
!
      return
      end
!
!=============================================
!
      subroutine evpot2 (spre,temp,ur,rn,evap)                                                                                                                   !Output
!
! Inputs
! ------
!
      real spre                                                                                !Surface pressure (mb)
      real temp                                                                                !Temperature (oC)
      real ur                                                                                  !Relative humidity (0-1,dimensionless)
      real rn                                                                                  !Irradiation balance (W/m2)
!
! Output
! ------
!
      real evap                                                                                !Potencial evapotranspiration without stress (mm/day)
!
! Parameters
! ----------
!
      ra      = 100.                                                                           
      rcmin   = 100.                                                                         
!
! Delta
! -----
!
      t1 = temp + 1.
      t2 = temp - 1.
      call tetens(t1,es1)
      call tetens(t2,es2)
      delta = (es1-es2)/(t1-t2)                                                                
!
! Delta_e
! -------
!
      call tetens (temp,es)
      delta_e = es*(1. - ur)                                                                  
!
! Stomatal Conductance
! --------------------
!
      rc = rcmin
!
! Gama and gama2
! --------------
!
      gama  = spre*(1004.)/(2.45e6*0.622)
      gama2 = gama*(ra + rc)/ra
!
! Potencial evapotranspiration (without stress)
! ---------------------------------------------
!
      evap = (delta*rn + (1.20*1004./ra)*delta_e)/(delta+gama2)                               
      evap = evap*(86400./2.45e6)                                                              
      evap = amax1(evap,0.)                                                                    
!
      return
      end
!
! ====================================================================
!
      subroutine runoff (w,wmax,roff)
      real w,roff
!     roff = 38.*((w/wmax)**11.)                                                            
      roff = 11.5*((w/wmax)**6.6)                                                             
      return
      end
!
! =====================================================================
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
