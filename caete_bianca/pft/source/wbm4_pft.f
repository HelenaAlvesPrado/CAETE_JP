c234567

c       definicao da helena
c      subroutine wbm (prec,temp,lsmk,p0,ca,par,
c     &     emaxm, tsoil, photo_pft,aresp_pft,npp_pft,lai_pft,
c     &     clit_pft,csoil_pft, hresp_pft,rcm_pft,runom_pft,
c     &     evapm_pft,wsoil_pft)



      subroutine wbm (prec,temp,lsmk,p0,ca,par,
     &                cleafini_pft,cawoodini_pft,cfrootini_pft,
     &                tmin,meanpp,seanpp,mphoto,maresp,  	!output
     &                meanrml,meanrmf,meanrms,meanrm,
     &                meanrgl,meanrgf,meanrgs,meanrg,
     &                meanhr,meancs,wsoil2,evapm,npp,
     &                photo,aresp,rcm,ave_wsoil,ave_evap,ave_rc,
     &                monrml,monrmf,monrms,monrm,
     &                monrgl,monrgf,monrgs,monrg,
     &                ctotal_gd,
     &                total_npp,total_ar,total_ph,
     &                ctotal_pft, mean_npp_pft, mean_ar_pft,
     &                mean_photo_pft,
     &                mean_cleaf_pft,mean_cawood_pft,mean_cfroot_pft)

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
      parameter (npft=3,npls=3)
      parameter(nx=192,ny=96)
      real prec(nx,ny,12),temp(nx,ny,12),lsmk(nx,ny),p0(nx,ny)
      real tmin(nx,ny),seanpp(nx,ny),meanpp(nx,ny),meanhr(nx,ny),
     &     meancs(nx,ny),mphoto(nx,ny),maresp(nx,ny),ave_wsoil(nx,ny),
     &     ave_evap(nx,ny),ave_rc(nx,ny),
     &     mean_npp_tbe(nx,ny),mean_npp_tbd(nx,ny),mean_npp_herb(nx,ny),
     &     mean_npp_pft(nx,ny,npls),npp_pft(nx,ny,12,npls),
     &     nppavg_pft(npls),
     &     mean_photo_pft(nx,ny,npls),
     &     photo_pft(nx,ny,12,npls),photoavg_pft(npls),
     &     mean_ar_pft(nx,ny,npls),
     &     ar_pft(nx,ny,12,npls),aravg_pft(npls),
     &     total_npp(nx,ny),total_ph(nx,ny),total_ar(nx,ny)
      real meanrml(nx,ny),meanrmf(nx,ny),meanrms(nx,ny),meanrm(nx,ny),
     &     meanrgl(nx,ny),meanrgf(nx,ny),meanrgs(nx,ny),meanrg(nx,ny)
      real wsoil2(nx,ny,12),par(nx,ny,12)
      real photo(nx,ny,12),aresp(nx,ny,12),npp(nx,ny,12),
     & lai(nx,ny,12),clit(nx,ny,12),csoil(nx,ny,12),hresp(nx,ny,12)
      real monrml(nx,ny,12),monrmf(nx,ny,12),monrms(nx,ny,12),
     &    monrm(nx,ny,12),monrgl(nx,ny,12),monrgf(nx,ny,12),
     &    monrgs(nx,ny,12),monrg(nx,ny,12)
      integer npls2 !number of plant life strategies
	 
c carbon allocation
      real 
     & cleafini_pft(nx,ny,npls),cleaf1_pft(npls),
     & cawoodini_pft(nx,ny,npls),cawood1_pft(npls),
     & cfrootini_pft(nx,ny,npls),cfroot1_pft(npls),
     & mean_cleaf_pft(nx,ny,npls),
     & cleaf_pft(nx,ny,12,npls),
     & cl2_pft(npls),cleafavg_pft(npls),
     & cleafi_pft(npls),cleaff_pft(npls),
     & mean_cawood_pft(nx,ny,npls),
     & cawood_pft(nx,ny,12,npls),
     & ca2_pft(npls),cawoodavg_pft(npls),
     & cawoodi_pft(npls),cawoodf_pft(npls),
     & mean_cfroot_pft(nx,ny,npls),
     & cfroot_pft(nx,ny,12,npls),
     & cf2_pft(npls),cfrootavg_pft(npls),
     & cfrooti_pft(npls),cfrootf_pft(npls),
     & ctotal_gd(nx,ny),ctotal_pft(nx,ny,npls),
     & ocp_pft(nx,ny,npls)
      integer i6 !index for calculate the 3 PFTs (i6.eq.1 - TBE; i6.eq.2 - TBD; i6.eq.3 - HERB)
      integer i1,i2 
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
			
      do k=1,12

      do i6=1,npls
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
		monrml(i,j,k) = -999.99 !mean monthly leaf maintenance respiration (KgC/m2/day)
		monrmf(i,j,k) = -999.99 !mean monthly fine roots maintenance respiration (KgC/m2/day)
		monrms(i,j,k) = -999.99 !mean monthly sapwood maintenance respiration (KgC/m2/day)
		monrm(i,j,k) = -999.99  !mean monthly total maintenance respiration (KgC/m2/day)
		monrgl(i,j,k) = -999.99 !mean monthly leaf growth respiration (KgC/m2/day)
		monrgf(i,j,k) = -999.99 !mean monthly fine roots growth respiration (KgC/m2/day)
		monrgs(i,j,k) = -999.99 !mean monthly sapwood growth respiration (KgC/m2/day)
		monrg(i,j,k) = -999.99  !mean monthly total growth respiration (KgC/m2/day)
        npp(i,j,k) = -999.99
	    clit(i,j,k) = -999.99
	    csoil(i,j,k) = -999.99
	    hresp(i,j,k) = -999.99
		    
        cleaf_pft(i,j,k,i6)= -999.99 !mean monthly leaf biomass for all the PFTs (KgC/m2)
		cawood_pft(i,j,k,i6)= -999.99 !mean monthly aboveground biomass for all the PFTs (KgC/m2)
		cfroot_pft(i,j,k,i6)= -999.99 !mean monthly leaf biomass for all the PFTs (KgC/m2)
		npp_pft(i,j,k,i6)= -999.99!mean monthly npp for all the PFTs (KgC/m2/yr)
		photo_pft(i,j,k,i6)= -999.99!mean monthly photosynthesis for all the PFTs (KgC/m2/yr)
		ar_pft(i,j,k,i6)= -999.99!mean monthly autotrophic respiration for all the PFTs (KgC/m2/yr)
		     enddo
      enddo
	    
c
c only for land grid points
      if (int(lsmk(i,j)).ne.0) then
c
c set some variables
      wini  = 0.01  !soil moisture initial condition (mm)
      gini  = 0.0   !soil ice initial condition (mm)
      sini  = 0.0   !overland snow initial condition (mm)
      do i6=1,npls
		 cleafi_pft(i6) = 0.0 !inital leaf biomass for all the PFTs (just for initialization)
		 cawoodi_pft(i6) = 0.0 !inital aboveground biomass for all the PFTs (just for initialization)
		 cfrooti_pft(i6) = 0.0 !inital leaf biomass for all the PFTs (just for initialization)
		 cleaf1_pft(i6)= cleafini_pft(i,j,i6)!inital leaf biomass for each PFT from spinup (KgC/m2) (just for january 1st)
         cawood1_pft(i6)=cawoodini_pft(i,j,i6)!inital leaf biomass for each PFT from spinup (KgC/m2) (just for january 1st)
		 cfroot1_pft(i6)=cfrootini_pft(i,j,i6)!inital leaf biomass for each PFT from spinup (KgC/m2) (just for january 1st)
	     enddo
	  
               endif
            endif
         enddo
      enddoc initialization
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
!      ae = 2.26457*ta + 67.5876 !available energy (W/m2) [Eq. 8]
      ae = 2.895*ta + 52.326 !from NCEP-NCAR Reanalysis data
c
        
c prepare inputs escalares para budget (a integração que comeca na linha
c 203 precisa ser chamada para cada pft ... )
     
c monthly water budget
      call budget (mes,wini,gini,sini,td,ta,pr,spre,ae,ca,ipar,
     &              cleafi_pft,cawoodi_pft,cfrooti_pft,
     &				cleaf1_pft,cawood1_pft,cfroot1_pft,
     &              wfim,gfim,sfim,smes,rmes,emes,epmes,
     &              phmes,armes,nppmes,laimes,
     &              rmlmes,rmfmes,rmsmes,rmmes,
     &              rglmes,rgfmes,rgsmes,rgmes,
     &              clmes,csmes,hrmes,rcmes,
     &              cleaff_pft,cawoodf_pft,cfrootf_pft, 
     &              cleafavg_pft,cawoodavg_pft,cfrootavg_pft,
     &              nppavg_pft,photoavg_pft,
     &              aravg_pft)  

      
c update variables
       do i6=1,npls
		cleaf_pft(i,j,k,i6) = cleafavg_pft(i6)
		cawood_pft(i,j,k,i6) = cawoodavg_pft(i6)
		cfroot_pft(i,j,k,i6) = cfrootavg_pft(i6)
		npp_pft(i,j,k,i6) = nppavg_pft(i6)
		photo_pft(i,j,k,i6) = photoavg_pft(i6)
		ar_pft(i,j,k,i6) = aravg_pft(i6)
      enddo
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
		monrml(i,j,k) = rmlmes
		monrmf(i,j,k) = rmfmes
		monrms(i,j,k) = rmsmes
		monrm(i,j,k) = rmmes
		monrgl(i,j,k) = rglmes
		monrgf(i,j,k) = rgfmes
		monrgs(i,j,k) = rgsmes
		monrg(i,j,k) = rgmes
        npp(i,j,k) = nppmes
	    clit(i,j,k) = clmes
	    csoil(i,j,k) = csmes
	    hresp(i,j,k) = hrmes
        wini = wfim
        gini = gfim
        sini = sfim
	     do i6=1,npls			
        cleafi_pft(i6) = cleaff_pft(i6)!it turns the carbon content of the last day of the previous month on the carbon content of the 1st day of the next month(for all the PFTs)
	    cawoodi_pft(i6) = cawoodf_pft(i6)!it turns the carbon content of the last day of the previous month on the carbon content of the 1st day of the next month(for all the PFTs)     
        cfrooti_pft(i6) = cfrootf_pft(i6)!it turns the carbon content of the last day of the previous month on the carbon content of the 1st day of the next month(for all the PFTs)    	       
	     enddo

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
c
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
	  meanrml(i,j) = 0.0
	  meanrmf(i,j) = 0.0
	  meanrms(i,j) = 0.0
	  meanrm(i,j) = 0.0
	  meanrgl(i,j) = 0.0
	  meanrgf(i,j) = 0.0
	  meanrgs(i,j) = 0.0
	  meanrg(i,j) = 0.0
      nppmin(i,j) = 100.0
      nppmax(i,j) = -100.0
      meanhr(i,j) = 0.0
      meancs(i,j) = 0.0
      ave_wsoil(i,j) = 0.0
      ave_evap(i,j) = 0.0
      ave_rc(i,j) = 0.0
         do i6=1,npls
         mean_cleaf_pft(i,j,i6)=0.0 !mean annual biomass on leaf compartment for each PFT (kgC/m2)
		 mean_cawood_pft(i,j,i6)=0.0 !mean annual biomass on aboveground woody compartment for each PFT (kgC/m2)
		 mean_cfroot_pft(i,j,i6)=0.0!mean annual biomass on fine roots compartment for each PFT (kgC/m2)
		 mean_npp_pft(i,j,i6)=0.0 !mean annual NPP for each PFT (kgC/m2/year)
		 mean_photo_pft(i,j,i6)=0.0 !mean annual photosynthesis for each PFT (kgC/m2/year)
		 mean_ar_pft(i,j,i6)=0.0 !mean annual autotrophic respiration for each PFT (kgC/m2/year)
	     enddo

c
c calculate tmin, meanpp, seanpp [Eqs 11, ...]
      do k=1,12
         if (temp(i,j,k).lt.tmin(i,j)) tmin(i,j) = temp(i,j,k)
	 meanpp(i,j) = meanpp(i,j) + (npp(i,j,k)/12)
	 mphoto(i,j) = mphoto(i,j) + (photo(i,j,k)/12)
	 maresp(i,j) = maresp(i,j) + (aresp(i,j,k)/12)
	 meanrml(i,j) = meanrml(i,j) + (monrml(i,j,k)/12)
	 meanrmf(i,j) = meanrmf(i,j) + (monrmf(i,j,k)/12)
	 meanrms(i,j) = meanrms(i,j) + (monrms(i,j,k)/12)
	 meanrm(i,j) = meanrm(i,j) + (monrm(i,j,k)/12)
	 meanrgl(i,j) = meanrgl(i,j) + (monrgl(i,j,k)/12)
	 meanrgf(i,j) = meanrgf(i,j) + (monrgf(i,j,k)/12)
	 meanrgs(i,j) = meanrgs(i,j) + (monrgs(i,j,k)/12)
	 meanrg(i,j) = meanrg(i,j) + (monrg(i,j,k)/12)
	 meanhr(i,j) = meanhr(i,j) + (hresp(i,j,k)/12)
	 meancs(i,j) = meancs(i,j) + (csoil(i,j,k)/12)
         ave_wsoil(i,j) = ave_wsoil(i,j) + (wsoil(i,j,k)/12)
         ave_evap(i,j) = ave_evap(i,j) + (evapm(i,j,k)/12)
         ave_rc(i,j) = ave_rc(i,j) + (rcm(i,j,k)/12)
		 
         if (npp(i,j,k).lt.nppmin(i,j)) nppmin(i,j) = npp(i,j,k)
         if (npp(i,j,k).gt.nppmax(i,j)) nppmax(i,j) = npp(i,j,k)
		    
c it calculates mean annual biomass and mean annual carbon fluxes for all the PFTs
			   do i6=1,npls
               mean_cleaf_pft(i,j,i6)=mean_cleaf_pft(i,j,i6)+
     &  		 (cleaf_pft(i,j,k,i6)/12)
	           mean_cawood_pft(i,j,i6)=mean_cawood_pft(i,j,i6)+
     &  		 (cawood_pft(i,j,k,i6)/12)
	           mean_cfroot_pft(i,j,i6)=mean_cfroot_pft(i,j,i6)+
     &  		 (cfroot_pft(i,j,k,i6)/12)
	           mean_npp_pft(i,j,i6)=mean_npp_pft(i,j,i6)+
     &  		 (npp_pft(i,j,k,i6)/12)
	           mean_photo_pft(i,j,i6)=mean_photo_pft(i,j,i6)+
     &  		 (photo_pft(i,j,k,i6)/12)
	           mean_ar_pft(i,j,i6)=mean_ar_pft(i,j,i6)+
     &  		 (ar_pft(i,j,k,i6)/12)
	        
	         enddo
      enddo
	  
c
c it calculates the total biomass for all the PFTs	 
	     do i6=1,npls
	         ctotal_pft(i,j,i6)=mean_cleaf_pft(i,j,i6)+
     &	            mean_cawood_pft(i,j,i6)+mean_cfroot_pft(i,j,i6)
	     enddo
	 
c it calculates the total biomass in a grid cell	    
         do i6=1,npls
		    ctotal_gd(i,j)= ctotal_gd(i,j)+ctotal_pft(i,j,i6)
		    enddo	

c it calculates the proportional occupation of PFT in a grid cell through its total biomass	   
		    do i6=1,npls      
		 ocp_pft(i,j,i6)= ctotal_pft(i,j,i6)/ctotal_gd(i,j)
		     if(ctotal_pft(i,j,i6).le.0.0) then
		        ocp_pft(i,j,i6)=0.0
		     endif
		    enddo
			   
c it recalculates the total biomass and the biomass of each compartment of each PFT from its proportional occupation (Pavlick et al., 2013)		 
		    do i6=1,npls
	         ctotal_pft(i,j,i6)=ctotal_pft(i,j,i6)*ocp_pft(i,j,i6)
			 ctotal_pft(i,j,i6)=ctotal_pft(i,j,i6)
			 mean_cleaf_pft(i,j,i6)=mean_cleaf_pft(i,j,i6)*ocp_pft(i,j,i6)
			 mean_cawood_pft(i,j,i6)=mean_cawood_pft(i,j,i6)*ocp_pft(i,j,i6)
			 mean_cfroot_pft(i,j,i6)=mean_cfroot_pft(i,j,i6)*ocp_pft(i,j,i6)
               enddo

c it recalculates the total biomass in a grid cell 	
             do i6=1,npls
		    ctotal_gd(i,j)= 0
		    enddo	
		 
		    do i6=1,npls
		    ctotal_gd(i,j)= ctotal_gd(i,j)+ctotal_pft(i,j,i6)
		    enddo	
		
c it calculates total npp, autotrophic respiration and photosynthesis in a grid cell 
          do i6=1,npls
		    total_npp(i,j)= total_npp(i,j) + mean_npp_pft(i,j,i6)
			total_ar(i,j)= total_ar(i,j) + mean_ar_pft(i,j,i6)
			total_ph(i,j)= total_ph(i,j) + mean_photo_pft(i,j,i6)
		    enddo	

	 
	 
	 

c
c
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
      subroutine budget (month,w1,g1,s1,tsoil,temp,prec,p0,ae,ca,ipar,  !input
     &                     cleafi_pft,cawoodi_pft,cfrooti_pft,
     &                     cleaf1_pft,cawood1_pft,cfroot1_pft,
     &                     w2,g2,s2,smavg,ruavg,evavg,	 	            !output
     &                     epavg,phavg,aravg,nppavg,laiavg,	            !output
     &                     rmlavg,rmfavg,rmsavg,rmavg,					!output
     &                     rglavg,rgfavg,rgsavg,rgavg,					!output
     &                     clavg,csavg,hravg,rcavg,                     !output
     &                     cleaf2_pft,cawood2_pft,cfroot2_pft,
     &                     cleafavg_pft,cawoodavg_pft,cfrootavg_pft,
     &                     nppavg_pft,photoavg_pft,
     &                     aravg_pft)         			!output
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
      parameter (npft=3,npls=3)!npls:number of pfts
      integer month
      real w1,g1,s1,tsoil,temp,prec,p0,ae,ca,ipar,
     &     w2,g2,s2,smavg,ruavg,evavg,epavg,rcavg,
     &     phavg,aravg,nppavg,laiavg,
     &     clavg,csavg,hravg     
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
     &     rglavg,rgfavg,rgsavg,rgavg,
     &     sla,
     &     nppavg_pft(npls),npp2_pft(npls),
     &     photoavg_pft(npls),photo2_pft(npls),
     &     aravg_pft(npls),ar2_pft(npls)
c carbon allocation
      real cleaf1_pft(npls),cawood1_pft(npls),cfroot1_pft(npls),
     &     cleafi_pft(npls),cawoodi_pft(npls),
     &     cfrooti_pft(npls),
     &     cleaf2_pft(npls),cawood2_pft(npls),cfroot2_pft(npls),
     &     cleafavg_pft(npls),cawoodavg_pft(npls),cfrootavg_pft(npls),	 
     &     cl1,cl2,ca1,ca2,cf1,cf2,cb1,cb2,cs1,cs2,
     &     cr1,cr2,co1,co2,
     &     cl2_pft(npls),ca2_pft(npls),cf2_pft(npls),
     &     alfa_leaf(npls),alfa_awood(npls),alfa_froot(npls),
     &     beta_leaf, beta_awood,beta_froot,
     &     ctotal_pft2(npls),!Daily total biomass for each PFT 
     &     ctotal_gd2,      !Daily total biomass in a grid cell 
     &     ctotali_gd,     !Initial Daily total biomass in a grid cell 
     &     cawood_gd,	 ! Aboveground woody biomass in a grid cell
     &     ocp_wood(npls),
     &     max_value, !to find the greatest value of ocp_wood
     &     max_index, !to find the index of the greatest value of ocp_wood
     &     cawoodi_gd, ! Initial Aboveground woody biomass in a grid cell
     &     ocp_pft2(npls)!PFTs proportional occupation in a grid cell(pft_biomass/gridcell_total_biomass)
      real ocpi_wood_pft(npls)
      real ctotali_pft(npls),ocpi_pft(npls)
      real alc_leaf2,alc_froot2,alc_awood2	 
	 

					 
				 
	  
	  
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
	  rmlavg = 0.
	  rmfavg = 0.
	  rmsavg = 0.
	  rmavg = 0.
	  rglavg = 0.
	  rgfavg = 0.
	  rgsavg = 0.
	  rgavg = 0.
      nppavg = 0.
      clavg = 0.
      csavg = 0.
      hravg = 0.
	  ctotal_gd2=0. !it calculates the total biomass of a grid  cell in daily bases
	  cawood_gd = 0. !it calculates the total aboveground woody biomass in a grid cell (it is used to calculate the proportional uptake of light of each PLS)
	  c=======================================================================
c234567
      subroutine budget (month,w1,g1,s1,tsoil,temp,prec,p0,ae,ca,ipar,  !input
     &                     cleafi_pft,cawoodi_pft,cfrooti_pft,
     &                     cleaf1_pft,cawood1_pft,cfroot1_pft,
     &                     w2,g2,s2,smavg,ruavg,evavg,	 	            !output
     &                     epavg,phavg,aravg,nppavg,laiavg,	            !output
     &                     rmlavg,rmfavg,rmsavg,rmavg,					!output
     &                     rglavg,rgfavg,rgsavg,rgavg,					!output
     &                     clavg,csavg,hravg,rcavg,                     !output
     &                     cleaf2_pft,cawood2_pft,cfroot2_pft,
     &                     cleafavg_pft,cawoodavg_pft,cfrootavg_pft,
     &                     nppavg_pft,photoavg_pft,
     &                     aravg_pft)         			!output
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
      parameter (npft=3,npls=3)!npls:number of pfts
      integer month
      real w1,g1,s1,tsoil,temp,prec,p0,ae,ca,ipar,
     &     w2,g2,s2,smavg,ruavg,evavg,epavg,rcavg,
     &     phavg,aravg,nppavg,laiavg,
     &     clavg,csavg,hravg     
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
     &     rglavg,rgfavg,rgsavg,rgavg,
     &     sla,
     &     nppavg_pft(npls),npp2_pft(npls),
     &     photoavg_pft(npls),photo2_pft(npls),
     &     aravg_pft(npls),ar2_pft(npls)
c carbon allocation
      real cleaf1_pft(npls),cawood1_pft(npls),cfroot1_pft(npls),
     &     cleafi_pft(npls),cawoodi_pft(npls),
     &     cfrooti_pft(npls),
     &     cleaf2_pft(npls),cawood2_pft(npls),cfroot2_pft(npls),
     &     cleafavg_pft(npls),cawoodavg_pft(npls),cfrootavg_pft(npls),	 
     &     cl1,cl2,ca1,ca2,cf1,cf2,cb1,cb2,cs1,cs2,
     &     cr1,cr2,co1,co2,
     &     cl2_pft(npls),ca2_pft(npls),cf2_pft(npls),
     &     alfa_leaf(npls),alfa_awood(npls),alfa_froot(npls),
     &     beta_leaf, beta_awood,beta_froot,
     &     ctotal_pft2(npls),!Daily total biomass for each PFT 
     &     ctotal_gd2,      !Daily total biomass in a grid cell 
     &     ctotali_gd,     !Initial Daily total biomass in a grid cell 
     &     cawood_gd,	 ! Aboveground woody biomass in a grid cell
     &     ocp_wood(npls),
     &     max_value, !to find the greatest value of ocp_wood
     &     max_index, !to find the index of the greatest value of ocp_wood
     &     cawoodi_gd, ! Initial Aboveground woody biomass in a grid cell
     &     ocp_pft2(npls)!PFTs proportional occupation in a grid cell(pft_biomass/gridcell_total_biomass)
      real ocpi_wood_pft(npls)
      real ctotali_pft(npls),ocpi_pft(npls)
      real alc_leaf2,alc_froot2,alc_awood2	 
	 

					 
				 
	  
	  
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
      cawoodi_gd = 0. !it calculates the initial total aboveground woody biomass in a grid cell(it is used to calculate the proportional uptake of light of each PLS; just for january 1st)
	  ctotali_gd=0 !it calculates the initial total biomass in a grid cell(just for january 1st)
	     max_value= -100.
		 
	     do i6 = 1,npls
		    cleafavg_pft(i6) = 0.!mean monthly leaf biomass for all the PFTs
			cawoodavg_pft(i6) = 0.!mean monthly aboveground biomass for all the PFTs
			cfrootavg_pft(i6) = 0.!mean monthly belowground biomass for all the PFTs
			nppavg_pft(i6) = 0. !mean monthly npp for all the PFTs
			photoavg_pft(i6) = 0. !mean monthly photosynthesis for all the PFTs
			aravg_pft(i6) = 0. !mean monthly autotrophic respiration for all the PFTs
			
			cl2_pft(i6)=cleafi_pft(i6)
			ca2_pft(i6)=cawoodi_pft(i6)
			cf2_pft(i6)=cfrooti_pft(i6)

			ctotal_pft2(i6)=0.!total biomass of all PFTs(kgC/m2)
			ocp_wood(i6)=0.
			 ocpi_wood_pft(i6)=0.
			
		     enddo
			
           do i6=1,npls			
          ctotali_pft(i6)=cleaf1_pft(i6)+cawood1_pft(i6)+cfroot1_pft(i6)
		  ctotali_gd = ctotali_gd+ctotali_pft(i6)	  
		     enddo
			 
			   do i6=1,npls
			 ocpi_pft(i6)=ctotali_pft(i6)/ctotali_gd
             	if (ctotali_pft(i6).eq.0.) then
				    ocpi_pft(i6)=0.
					  endif
                enddo			   
			 
              do i6=1,npls
			  cawoodi_gd = cawoodi_gd+cawood1_pft(i6)
			    enddo
				
	          do i6=1,npls
			  ocpi_wood_pft(i6) = cawood1_pft(i6)/cawoodi_gd
			    if (cawood1_pft(i6).eq.0.) then
				    ocpi_wood_pft(i6)=0.
					  endif
			    enddo
					  
			    do i6=1,npls
                 if(ocpi_wood_pft(i6).gt.max_value) then
                     max_value=ocpi_wood_pft(i6)
					 maxi_index=i6
					
                  endif
                  enddo				  
        
		          
c
c	  

c numerical integration
      do i=1,ndmonth(month)
	          ctotal_gd2 = 0
			  cawood_gd = 0
	           i6=0
			    	do i1 = 20,80,40 !e.g. C-leaf
	              do i2 = 20,80,20
	              do i3=0,80,40
				  
				  if(i1+i2+i3.eq.100) then
					alc_leaf2=i1
					alc_froot2=i2
					alc_awood2=i3
					
					i6=i6+1
				   	 	 
			 
			     cl1 = cl2_pft(i6) !transforma o valor do dia anterior no valor atual
			      
				 ca1 = ca2_pft(i6)
				 cf1 = cf2_pft(i6)
				

				 beta_leaf = alfa_leaf(i6)
				 beta_awood = alfa_awood(i6)
				 beta_froot = alfa_froot(i6)

				 
	         if ((i.eq.1).and.(month.eq.1)) then
				     cl1 = cleaf1_pft(i6)
					 ca1 = cawood1_pft(i6)
					 cf1 = cfroot1_pft(i6)

					 beta_leaf=0.
					 beta_awood=0.

					 beta_froot=0.
					 ocp_pft2(i6)=ocpi_pft(i6)
					 max_index=maxi_index
			    endif
				  
! TEM UM PROBLRMA GRAVE COM A VATIAVEL RC2 NA CHAMADE DE CARBON 1/ ELA NAO EXISTE 
! DENTRO DE CARBON1 E POR ISSO ACONTECE UMA POSSIVEL DIVISAO POR 0 DENTRO DA CARBON1              
				
c carbon cycle (photosynthesis, plant respiration and NPP)
      call carbon1 (temp,p0,w,wmax,ca,ipar,i6,tsoil,emax,rc2,  !input
     &               cl1,ca1,cf1,
     &               max_index,
     &               ocp_pft2,
     &               beta_leaf,beta_awood,beta_froot,          !input
     &              ph,ar,nppa,laia,f5,rm,rml,rmf,rms,rg,rgl,rgf,rgs)		        !output

	          
	 
!     carbon allocation (carbon content on each compartment)     
       call allocation (nppa,cl1,ca1,cf1,i6,	   !input
     &                 alc_leaf2,alc_froot2,alc_awood2,  
     &                  cl2,ca2,cf2)  !output 

	            alfa_leaf(i6)  =cl2 - cl1 !it calculates the difference btween the biomass of the last and actual day (for growth respiration)
				alfa_awood(i6) =ca2 - ca1 !it calculates the difference btween the biomass of the last and actual day (for growth respiration)  
				alfa_froot(i6) =cf2 - cf1 !it calculates the difference btween the biomass of the last and actual day (for growth respiration) 
				
				cl2_pft(i6) = cl2
				ca2_pft(i6) = ca2
				cf2_pft(i6) = cf2
	          

				  npp2_pft(i6)=nppa
				  photo2_pft(i6)=ph
				  ar2_pft(i6)=ar
				

c it calculates the total biomass for each PFT in a daily bases 
			ctotal_pft2(i6)=cl2_pft(i6)+ca2_pft(i6)+cf2_pft(i6)
			  
c it calculates the total biomass in a grid cell in a daily bases 
                  ctotal_gd2 = ctotal_gd2 + ctotal_pft2(i6)
				  
c it calculates the total aboveground woody biomass in a grid cell in a daily bases (it is used to calculate the proportional uptake of light of each PLS)				  
				  cawood_gd = cawood_gd + ca2_pft(i6)
				 
				   endif
				  enddo
				  enddo
				  enddo
            			  
		     
c it calculates the proportional occupation of a PFT in a grid cell				  
   				do i6=1, npls
				ocp_pft2(i6) = ctotal_pft2(i6)/ctotal_gd2
				   if(ctotal_pft2(i6).le.0.)then
                       ocp_pft2(i6) = 0.
					  endif		  
                 enddo					  

			 
c it calculates the proportional aboveground woody biomass of a PFT in a grid cell	
					  do i6=1,npls
                       ocp_wood(i6)=ca2_pft(i6)/cawood_gd
					   if ((ca2_pft(i6)).eq.0.) then
					    ocp_wood(i6)=0.
						endif
					   if (ocp_wood(i6).gt.max_value)then
					    max_value=ocp_wood(i6)
						max_index=i6
						
						endif
                      enddo					   
				    	
						
					
			 
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

       do i6 = 1,npls
	       cleafavg_pft(i6) = cleafavg_pft(i6) + cl2_pft(i6)
		   cawoodavg_pft(i6) = cawoodavg_pft(i6) + ca2_pft(i6)
		   cfrootavg_pft(i6) = cfrootavg_pft(i6) + cf2_pft(i6)
		   nppavg_pft(i6) = nppavg_pft(i6) + npp2_pft(i6)/365.0 !kgC/m2
		   photoavg_pft(i6) = photoavg_pft(i6) + photo2_pft(i6)/365.0 !kgC/m2
		   aravg_pft(i6) = aravg_pft(i6) + ar2_pft(i6)/365.0 !kgC/m2
       enddo
    	 
		
c

       enddo
c
c final calculations
        
       do i6=1,npls
	  cleaf2_pft(i6) = cl2_pft(i6)
	  cawood2_pft(i6) = ca2_pft(i6)
	  cfroot2_pft(i6) = cf2_pft(i6)

	      enddo
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
	  rmlavg = rmlavg*12.0
	  rmfavg = rmfavg*12.0
	  rmsavg = rmsavg*12.0
	  rmavg = rmavg*12.0
	  rglavg = rglavg*12.0
	  rgfavg = rgfavg*12.0
	  rgsavg = rgsavg*12.0
	  rgavg = rgavg*12.0
      nppavg = nppavg*12.0 !kgC/m2/yr
      laiavg = laiavg*12.0
      clavg = clavg*12.0 !kgC/m2
      csavg = csavg*12.0 !kgC/m2
      hravg = hravg*12.0 !kgC/m2/yr
	       do i6=1,npls
   	  cleafavg_pft(i6) = cleafavg_pft(i6)/real(ndmonth(month))
	  cawoodavg_pft(i6) = cawoodavg_pft(i6)/real(ndmonth(month))
	  cfrootavg_pft(i6) = cfrootavg_pft(i6)/real(ndmonth(month))
	  nppavg_pft(i6) = nppavg_pft(i6)*12.0 !kgC/m2/yr
	  photoavg_pft(i6) = photoavg_pft(i6)*12.0 !kgC/m2/yr
	  aravg_pft(i6) = aravg_pft(i6)*12.0 !kgC/m2/yr
	        enddo
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

