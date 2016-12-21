c234567
      subroutine wbm (prec,temp,lsmk,p0,ca,par,
     &                cleafini,cfrootini,cawoodini,
     &                cbwoodini,cstoini,cotherini,crepini,
     &                cleafini_tbe,cfrootini_tbe,cawoodini_tbe,
     &				  cbwoodini_tbe,cstoini_tbe,
     &                cotherini_tbe,crepini_tbe,
     &                cleafini_tbd,cfrootini_tbd,cawoodini_tbd,
     &				  cbwoodini_tbd,cstoini_tbd,
     &                cotherini_tbd,crepini_tbd,
     &                cleafini_herb,cfrootini_herb,cawoodini_herb,
     &		          cbwoodini_herb,cstoini_herb,
     &                cotherini_herb,crepini_herb,
     &                cleafini_pft,cawoodini_pft,cfrootini_pft,
     &                tmin,meanpp,seanpp,mphoto,maresp,  	!output
     &                meanrml,meanrmf,meanrms,meanrm,
     &                meanrgl,meanrgf,meanrgs,meanrg,
     &                meanhr,meancs,wsoil2,evapm,npp,
     &                photo,aresp,rcm,ave_wsoil,ave_evap,ave_rc,
     &                monrml,monrmf,monrms,monrm,
     &                monrgl,monrgf,monrgs,monrg,
     &                cleaf,cawood,cfroot,cbwood,csto,crep,cother,
     &                mean_cleaf,mean_cawood,mean_cbwood,mean_cfroot,
     &				  mean_csto,mean_crep,mean_cother,
     &                mean_cleaf_tbe,mean_cawood_tbe,
     &				  mean_cbwood_tbe,mean_cfroot_tbe,
     &				  mean_csto_tbe,mean_crep_tbe,
     &				  mean_cother_tbe,
     &                mean_cleaf_tbd,mean_cawood_tbd,
     &				  mean_cbwood_tbd,mean_cfroot_tbd,
     &                mean_csto_tbd,mean_crep_tbd,
     &				  mean_cother_tbd,
     &                mean_cleaf_herb,mean_cawood_herb,
     &				  mean_cbwood_herb,mean_cfroot_herb,
     &                mean_csto_herb,mean_crep_herb,
     &				  mean_cother_herb,
     &                ctotal_tbe,ctotal_tbd,ctotal_herb,ctotal_gd,
     &                mean_npp_tbe,mean_npp_tbd,mean_npp_herb,
     &                mean_photo_tbe,mean_photo_tbd,mean_photo_herb,
     &	 		      mean_ar_tbe,mean_ar_tbd,mean_ar_herb,
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
c     i/o variables
c     INPUTS________________
      integer, parameter :: npft=3
      integer, parameter :: nx=192 ,ny=96
      real prec(nx,ny,12)
      real temp(nx,ny,12)
      real lsmk(nx,ny)
      real p0(nx,ny)
      real par(nx,ny,12)
      
C     inputs -- initial carbon content in plant tissues      
      real cleaf_pft(nx,ny,12,npft)
      real cawood_pft(nx,ny,12,npft)
      real cbwood_pft(nx,ny,12,npft)
      real cfroot_pft(nx,ny,12,npft)
      real csto_pft(nx,ny,12,npft)
      real crep_pft(nx,ny,12,npft)
      real cother_pft(nx,ny,12,npft)

C     _____________________


C     outputas
      
!     maintenance respiration
      real monrml(nx,ny,12,npft)
      real monrmf(nx,ny,12,npft)
      real monrms(nx,ny,12,npft)
      real monrm(nx,ny,12,npft)
      real monrgl(nx,ny,12,npft)
      real monrgf(nx,ny,12,npft)
      real monrgs(nx,ny,12,npft)
      real monrg(nx,ny,12,npft)

      
C     variaveis originais do wbm
      real tsoil(nx,ny,12)
      real emaxm(nx,ny,12)
      
      real npp(nx,ny,12,npft)
      real photo(nx,ny,12,npft)
      real aresp(nx,ny,12,npft)
      real lai(nx,ny,12,npft)
      real clit(nx,ny,12,npft)
      real csoil(nx,ny,12,npft)
      real hresp(nx,ny,12,npft)
      
      real wsoil(nx,ny,12,npft)
      real runom(nx,ny,12,npft)
      real evapm(nx,ny,12,npft)
      real rcm(nx,ny,12,npft)
      
!     INTERNAL
      
      real gsoil(nx,ny,12)
      real ssoil(nx,ny,12)
      real snowm(nx,ny,12) ,
      
      


      
      
      integer npls2             !number of plant life strategies



      
c carbon allocation
      real cleaf(nx,ny,12),cawood(nx,ny,12),cfroot(nx,ny,12),
     & cbwood(nx,ny,12),csto(nx,ny,12),crep(nx,ny,12),cother(nx,ny,12),
     & cleafini(nx,ny),cfrootini(nx,ny),cawoodini(nx,ny),
     & cbwoodini(nx,ny),cstoini(nx,ny),cotherini(nx,ny),crepini(nx,ny),
     & cleafini_tbe(nx,ny),cfrootini_tbe(nx,ny),cawoodini_tbe(nx,ny),
     & cbwoodini_tbe(nx,ny),cstoini_tbe(nx,ny),
     & cotherini_tbe(nx,ny),crepini_tbe(nx,ny),
     & cleafini_tbd(nx,ny),cfrootini_tbd(nx,ny),cawoodini_tbd(nx,ny),
     & cbwoodini_tbd(nx,ny),cstoini_tbd(nx,ny),
     & cotherini_tbd(nx,ny),crepini_tbd(nx,ny),
     & cleafini_herb(nx,ny),cfrootini_herb(nx,ny),cawoodini_herb(nx,ny),
     & cbwoodini_herb(nx,ny),cstoini_herb(nx,ny),
     & cotherini_herb(nx,ny),crepini_herb(nx,ny),
     & cleafini_pft(nx,ny,npls),cleaf1_pft(npls),
     & cawoodini_pft(nx,ny,npls),cawood1_pft(npls),
     & cfrootini_pft(nx,ny,npls),cfroot1_pft(npls),
     & mean_cleaf(nx,ny),mean_cawood(nx,ny),mean_cbwood(nx,ny),
     & mean_cfroot(nx,ny),mean_csto(nx,ny),mean_crep(nx,ny),
     & mean_cother(nx,ny),
     & mean_cleaf_pft(nx,ny,npls),
     & mean_cleaf_tbe(nx,ny),mean_cleaf_tbd(nx,ny),
     & mean_cleaf_herb(nx,ny),
     & cleaf_pft(nx,ny,12,npls),
     & cl2_pft(npls),cleafavg_pft(npls),
     & cleafi_pft(npls),cleaff_pft(npls),
     & mean_cawood_pft(nx,ny,npls),
     & mean_cawood_tbe(nx,ny),mean_cawood_tbd(nx,ny),
     & mean_cawood_herb(nx,ny),
     & cawood_pft(nx,ny,12,npls),
     & ca2_pft(npls),cawoodavg_pft(npls),
     & cawoodi_pft(npls),cawoodf_pft(npls),
     & mean_cbwood_pft(nx,ny,npls),
     & mean_cbwood_tbe(nx,ny),mean_cbwood_tbd(nx,ny),
     & mean_cbwood_herb(nx,ny),
     & cbwood_pft(nx,ny,12,npls),
     & cb2_pft(npls),cbwoodavg_pft(npls),
     & cbwoodi_pft(npls),cbwoodf_pft(npls),
     & mean_cfroot_pft(nx,ny,npls),
     & mean_cfroot_tbe(nx,ny),mean_cfroot_tbd(nx,ny),
     & mean_cfroot_herb(nx,ny),
     & cfroot_pft(nx,ny,12,npls),
     & cf2_pft(npls),cfrootavg_pft(npls),
     & cfrooti_pft(npls),cfrootf_pft(npls),
     & mean_csto_pft(nx,ny,npls),
     & mean_csto_tbe(nx,ny),mean_csto_tbd(nx,ny),
     & mean_csto_herb(nx,ny),
     & csto_pft(nx,ny,12,npls),
     & cs2_pft(npls),cstoavg_pft(npls),
     & cstoi_pft(npls),cstof_pft(npls),
     & mean_crep_pft(nx,ny,npls),
     & mean_crep_tbe(nx,ny),mean_crep_tbd(nx,ny),
     & mean_crep_herb(nx,ny),
     & crep_pft(nx,ny,12,npls),
     & cr2_pft(npls),crepavg_pft(npls),
     & crepi_pft(npls),crepf_pft(npls),
     & mean_cother_pft(nx,ny,npls),
     & mean_cother_tbe(nx,ny),mean_cother_tbd(nx,ny),
     & mean_cother_herb(nx,ny),	 
     & cother_pft(nx,ny,12,npls),
     & co2_pft(npls),cotheravg_pft(npls),
     & cotheri_pft(npls),cotherf_pft(npls),
     & ctotal_tbe(nx,ny),ctotal_tbd(nx,ny),ctotal_herb(nx,ny),
     & ctotal_gd(nx,ny),ctotal_pft(nx,ny,npls),
     & ocp_tbe(nx,ny),ocp_tbd(nx,ny),ocp_herb(nx,ny),ocp_pft(nx,ny,npls)
      integer i6 !index for calculate the 3 PFTs (i6.eq.1 - TBE; i6.eq.2 - TBD; i6.eq.3 - HERB)
      integer i1,i2 
c
c internal variables
      real H,diffu,tau,tsoil(nx,ny,12),t0,t1

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
c     write to track program execution
            if ((mod(j,ny).eq.0).and.(mod(i,10).eq.0))
     &           write(*,*) 'water balance:',i
c     
c     
c     it calculates the number of plant life strategies (PLS)
            npls2=0    
            do i1 = 20,80,40    !e.g. C-leaf
               do i2 = 20,80,20
                  do i3=0,80,40
                     if(i1+i2+i3.eq.100) then
                        npls2= npls2+1
                        
                     endif
                  enddo
               enddo
            enddo
                   
c     initialize variables (-999.99 is undef)
            tmin(i,j)   = -999.99
            seanpp(i,j) = -999.99
      meanpp(i,j) = -999.99
      mphoto(i,j) = -999.99
      maresp(i,j) = -999.99
	  meanrml(i,j) = -999.99 !mean annual leaf maintenance respiration (KgC/m2/yr)
	  meanrmf(i,j) = -999.99 !mean annual fine roots maintenance respiration (KgC/m2/yr)
	  meanrms(i,j) = -999.99 !mean annual sapwood maintenance respiration (KgC/m2/yr)
	  meanrm(i,j) = -999.99  !mean annual total maintenance respiration (KgC/m2/yr)
	  meanrgl(i,j) = -999.99 !mean annual leaf growth respiration (KgC/m2/yr)
	  meanrgf(i,j) = -999.99 !mean annual fine roots growth respiration (KgC/m2/yr)
	  meanrgs(i,j) = -999.99 !mean annual sapwood maintenance growth spiration (KgC/m2/yr)
	  meanrg(i,j) = -999.99  !mean annual total growth respiration (KgC/m2/yr)
      meanhr(i,j) = -999.99
      meancs(i,j) = -999.99
      ave_wsoil(i,j) = -999.99  !calculates annual average soil water
      ave_evap(i,j) = -999.99   !calculates annual average evapotranspiration
      ave_rc(i,j) = -999.99     !calculates annual average canopy resistance
	  mean_cleaf(i,j)= -999.99
	  mean_cawood(i,j)= -999.99
      mean_cfroot(i,j)= -999.99
!	  mean_cbwood(i,j)= -999.99
	  !mean_csto(i,j)= -999.99
	  !mean_crep(i,j)= -999.99
	  !mean_cother(i,j)= -999.99
	  mean_cleaf_tbe(i,j)=-999.99 !mean annual leaf biomass for the TBE PFT (kgC/m2)
	  mean_cawood_tbe(i,j)=-999.99 !mean annual aboveground biomass for the TBE PFT (kgC/m2)
	  mean_cfroot_tbe(i,j)=-999.99 !mean annual leaf biomass for the TBE PFT (kgC/m2)
!	  mean_cbwood_tbe(i,j)=-999.99 !mean annual belowground biomass for the TBE PFT (kgC/m2)
	  !mean_csto_tbe(i,j)=-999.99 !mean annual leaf biomass for the TBE PFT (kgC/m2)
	  !mean_crep_tbe(i,j)=-999.99 !mean annual aboveground biomass for the TBE PFT (kgC/m2)
	  !mean_cother_tbe(i,j)=-999.99 !mean annual belowground biomass for the TBE PFT (kgC/m2)
	  mean_cleaf_tbd(i,j)=-999.99 !mean annual leaf biomass for the TBD PFT (kgC/m2)
	  mean_cawood_tbd(i,j)=-999.99 !mean annual aboveground biomass for the TBD PFT (kgC/m2)
	  mean_cfroot_tbd(i,j)=-999.99 !mean annual leaf biomass for the TBD PFT (kgC/m2)
!	  mean_cbwood_tbd(i,j)=-999.99 !mean annual belowground biomass for the TBD PFT (kgC/m2)
	  !mean_csto_tbd(i,j)=-999.99 !mean annual aboveground biomass for the TBD PFT (kgC/m2)
	  !mean_crep_tbd(i,j)=-999.99 !mean annual belowground biomass for the TBD PFT (kgC/m2)
	  !mean_cother_tbd(i,j)=-999.99 !mean annual leaf biomass for the TBD PFT (kgC/m2)
	  mean_cleaf_herb(i,j)=-999.99 !mean annual leaf biomass for the HERB PFT (kgC/m2)
	  mean_cawood_herb(i,j)=-999.99 !mean annual aboveground biomass for the HERB PFT (kgC/m2)
	  mean_cfroot_herb(i,j)=-999.99 !mean annual leaf biomass for the HERB PFT (kgC/m2)
!	  mean_cbwood_herb(i,j)=-999.99 !mean annual belowground biomass for the HERB PFT (kgC/m2)
	  !mean_csto_herb(i,j)=-999.99 !mean annual aboveground biomass for the HERB PFT (kgC/m2)
	  !mean_crep_herb(i,j)=-999.99 !mean annual belowground biomass for the HERB PFT (kgC/m2)
	  !mean_cother_herb(i,j)=-999.99 !mean annual leaf biomass for the HERB PFT (kgC/m2)
      ctotal_tbe(i,j) = -999.99 ! it calculates the total biomass of the TBE PFT (kgC/m2)
	  ctotal_tbd(i,j) = -999.99 ! it calculates the total biomass of the TBD PFT (kgC/m2)
	  ctotal_herb(i,j) = -999.99 ! it calculates the total biomass of the HERB PFT (kgC/m2)
	  mean_npp_tbe(i,j) = -999.99 !it calculates the mean annual npp of the TBE PFT (kgC/m2/yr)
	  mean_npp_tbd(i,j) = -999.99 !it calculates the mean annual npp of the TBD PFT (kgC/m2/yr)
	  mean_npp_herb(i,j) = -999.99 !it calculates the mean annual npp of the HERB PFT (kgC/m2/yr)
	  mean_photo_tbe(i,j) = -999.99 !it calculates the mean annual photosynthesis of the TBE PFT (kgC/m2/yr)
	  mean_photo_tbd(i,j) = -999.99 !it calculates the mean annual photosynthesis of the TBD PFT (kgC/m2/yr)
	  mean_photo_herb(i,j) = -999.99 !it calculates the mean annual photosynthesis of the HERB PFT (kgC/m2/yr)
	  mean_ar_tbe(i,j) = -999.99 !it calculates the mean annual autotrophic respiration of the TBE PFT (kgC/m2/yr)
	  mean_ar_tbd(i,j) = -999.99 !it calculates the mean annual autotrophic respiration of the TBD PFT (kgC/m2/yr)
	  mean_ar_herb(i,j) = -999.99 !it calculates the mean annual autotrophic respiration of the HERB PFT (kgC/m2/yr)
   	   do i6=1,npls
	    mean_cleaf_pft(i,j,i6)=-999.99 !mean annual leaf biomass for all the PFTs (kgC/m2)
		mean_cawood_pft(i,j,i6)=-999.99 !mean annual aboveground biomass for all the PFTs (kgC/m2)
		mean_cfroot_pft(i,j,i6)=-999.99 !mean annual leaf biomass for all the PFTs (kgC/m2)
!		mean_cbwood_pft(i,j,i6)=-999.99 !mean annual aboveground biomass for all the PFTs (kgC/m2)
		!mean_csto_pft(i,j,i6)=-999.99 !mean annual leaf biomass for all the PFTs (kgC/m2)
		!mean_crep_pft(i,j,i6)=-999.99 !mean annual aboveground biomass for all the PFTs (kgC/m2)
		!mean_cother_pft(i,j,i6)=-999.99 !mean annual aboveground biomass for all the PFTs (kgC/m2)
		mean_npp_pft(i,j,i6) = -999.99 !it calculates the mean annual npp for all the PFTs (kgC/m2/yr)
		mean_photo_pft(i,j,i6) = -999.99 !it calculates the mean annual npp for all the PFTs (kgC/m2/yr)
		mean_ar_pft(i,j,i6) = -999.99 !it calculates the mean annual autotrophic respiration for all the PFTs (kgC/m2/yr)
		ctotal_pft(i,j,i6)= -999.99
      enddo
			
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
        cleaf(i,j,k) = -999.99  !average carbon content on leaf compartment in a month(KgC/m2)
		cawood(i,j,k) = -999.99 !average carbon content on abovegroung wood compartment in a month(KgC/m2)
		cfroot(i,j,k) = -999.99 !average carbon content on fine roots compartment in a month(KgC/m2)
		!cbwood(i,j,k) = -999.99 !average carbon content on belowground wood compartment in a month(KgC/m2)
		!csto(i,j,k) = -999.99 !average carbon content on storage compartment in a month(KgC/m2)
		!crep(i,j,k) = -999.99 !average carbon content on reproduction compartment in a month(KgC/m2)
		!cother(i,j,k) = -999.99 !average carbon content on other compartment in a month(KgC/m2)
		    do i6=1,npls
        cleaf_pft(i,j,k,i6)= -999.99 !mean monthly leaf biomass for all the PFTs (KgC/m2)
		cawood_pft(i,j,k,i6)= -999.99 !mean monthly aboveground biomass for all the PFTs (KgC/m2)
		cfroot_pft(i,j,k,i6)= -999.99 !mean monthly leaf biomass for all the PFTs (KgC/m2)
		!cbwood_pft(i,j,k,i6)= -999.99 !mean monthly aboveground biomass for all the PFTs (KgC/m2)
		!csto_pft(i,j,k,i6)= -999.99 !mean monthly leaf biomass for all the PFTs (KgC/m2)
		!crep_pft(i,j,k,i6)= -999.99 !mean monthly aboveground biomass for all the PFTs (KgC/m2)
		!cother_pft(i,j,k,i6)= -999.99 !mean monthly aboveground biomass for all the PFTs (KgC/m2)
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
		 !cbwoodi_pft(i6) = 0.0 !inital aboveground biomass for all the PFTs (just for initialization)
		 !cstoi_pft(i6) = 0.0 !inital leaf biomass for all the PFTs (just for initialization)
		 !crepi_pft(i6) = 0.0 !inital aboveground biomass for all the PFTs (just for initialization)
		 !cotheri_pft(i6) = 0.0 !inital aboveground biomass for all the PFTs (just for initialization)
		 cleaf1_pft(i6)= 0.4!cleafini_pft(i,j,i6)!inital leaf biomass for each PFT from spinup (KgC/m2) (just for january 1st)
         cawood1_pft(i6)= 12.!cawoodini_pft(i,j,i6)
		 cfroot1_pft(i6)=0.3 !cfrootini_pft(i,j,i6)
	     enddo
	  cleafi_tbe = cleafini_tbe(i,j)  !inital leaf biomass for TBE PFT (KgC/m2)
      cawoodi_tbe = cawoodini_tbe(i,j)!inital abovegroung woody biomass for TBE PFT (KgC/m2)
	  cfrooti_tbe = cfrootini_tbe(i,j)!inital fine roots biomass for TBE PFT (KgC/m2)
	  !cbwoodi_tbe = cbwoodini_tbe(i,j)!inital belowground woody biomass for TBE PFT (KgC/m2)
	  !cstoi_tbe = cstoini_tbe(i,j)    !inital storage biomass for TBE PFT (KgC/m2)
	  !crepi_tbe = crepini_tbe(i,j)    !inital reproduction biomass for TBE PFT (KgC/m2)
	  !cotheri_tbe = cotherini_tbe(i,j)!inital other biomass for TBE PFT (KgC/m2)
	  cleafi_tbd = cleafini_tbd(i,j)  !inital leaf biomass for TBD PFT (KgC/m2)
      cawoodi_tbd = cawoodini_tbd(i,j)!inital abovegroung woody biomass for TBD PFT (KgC/m2)
	  cfrooti_tbd = cfrootini_tbd(i,j)!inital fine roots biomass for TBD PFT (KgC/m2)
	  !cbwoodi_tbd = cbwoodini_tbd(i,j)!inital belowground woody biomass for TBD PFT (KgC/m2)
	  !cstoi_tbd = cstoini_tbd(i,j)    !inital storage biomass for TBD PFT (KgC/m2)
	  !crepi_tbd = crepini_tbd(i,j)    !inital reproduction biomass for TBD PFT (KgC/m2)
      !cotheri_tbd = cotherini_tbd(i,j)!inital other biomass for TBD PFT (KgC/m2)
	  cleafi_herb = cleafini_herb(i,j)  !inital leaf biomass for HERB PFT (KgC/m2)
      cawoodi_herb = cawoodini_herb(i,j)!inital abovegroung woody biomass for HERB PFT (KgC/m2)
	  cfrooti_herb = cfrootini_herb(i,j)!inital fine roots biomass for HERB PFT (KgC/m2)
	  !cbwoodi_herb = cbwoodini_herb(i,j)!inital belowground woody biomass for HERB PFT (KgC/m2)
	  !cstoi_herb = cstoini_herb(i,j)    !inital storage biomass for HERB PFT (KgC/m2)
	  !crepi_herb = crepini_herb(i,j)    !inital reproduction biomass for HERB PFT (KgC/m2)
      !cotheri_herb = cotherini_herb(i,j)!inital other biomass for HERB PFT (KgC/m2)
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
!      ae = 2.26457*ta + 67.5876 !available energy (W/m2) [Eq. 8]
      ae = 2.895*ta + 52.326 !from NCEP-NCAR Reanalysis data
c
        

     
c monthly water budget
      call budget (mes,wini,gini,sini,td,ta,pr,spre,ae,ca,ipar,
     &              cleafi,cawoodi,cfrooti,
     &              cbwoodi,cstoi,crepi,cotheri,
     &              cleafi_pft,cawoodi_pft,
     & 				cbwoodi_pft,cfrooti_pft,
     &				cstoi_pft,crepi_pft,
     & 				cotheri_pft,
     &              cleafi_tbe,cawoodi_tbe,cfrooti_tbe,
     &              cbwoodi_tbe,cstoi_tbe,crepi_tbe,cotheri_tbe,
     &              cleafi_tbd,cawoodi_tbd,cfrooti_tbd,
     &              cbwoodi_tbd,cstoi_tbd,crepi_tbd,cotheri_tbd,
     &              cleafi_herb,cawoodi_herb,cfrooti_herb,
     &              cbwoodi_herb,cstoi_herb,crepi_herb,cotheri_herb,
     &				cleaf1_pft,cawood1_pft,cfroot1_pft,
     &              wfim,gfim,sfim,smes,rmes,emes,epmes,
     &              phmes,armes,nppmes,laimes,
     &              rmlmes,rmfmes,rmsmes,rmmes,
     &              rglmes,rgfmes,rgsmes,rgmes,
     &              clmes,csmes,hrmes,rcmes,
     &              cleaff,cawoodf,cfrootf,cbwoodf,cstof,crepf,cotherf,
     &              cleafavg,cawoodavg,cfrootavg,
     &              cbwoodavg,cstoavg,crepavg,cotheravg,
     &              cleaff_pft,cawoodf_pft,
     &              cbwoodf_pft,cfrootf_pft, 
     &				cstof_pft,crepf_pft,
     &              cotherf_pft,
     &              cleafavg_pft,cawoodavg_pft,
     &              cbwoodavg_pft,cfrootavg_pft,
     &              cstoavg_pft,crepavg_pft,
     &              cotheravg_pft,nppavg_pft,photoavg_pft,
     &              aravg_pft)  

      
c update variables
       do i6=1,npls
		cleaf_pft(i,j,k,i6) = cleafavg_pft(i6)
		cawood_pft(i,j,k,i6) = cawoodavg_pft(i6)
		cfroot_pft(i,j,k,i6) = cfrootavg_pft(i6)
		!cbwood_pft(i,j,k,i6) = cbwoodavg_pft(i6)
		!csto_pft(i,j,k,i6) = cstoavg_pft(i6)
		!crep_pft(i,j,k,i6) = crepavg_pft(i6)
		!cother_pft(i,j,k,i6) = cotheravg_pft(i6)
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
	    cleaf(i,j,k) = cleafavg
		cawood(i,j,k) = cawoodavg
		cfroot(i,j,k) = cfrootavg
		!cbwood(i,j,k) = cbwoodavg
		!csto(i,j,k) = cstoavg
		!crep(i,j,k) = crepavg
		!cother(i,j,k) = cotheravg
        wini = wfim
        gini = gfim
        sini = sfim
	     do i6=1,npls			
        cleafi_pft(i6) = cleaff_pft(i6)!it turns the carbon content of the last day of the previous month on the carbon content of the 1st day of the next month(for all the PFTs)
	    cawoodi_pft(i6) = cawoodf_pft(i6)!it turns the carbon content of the last day of the previous month on the carbon content of the 1st day of the next month(for all the PFTs)     
        cfrooti_pft(i6) = cfrootf_pft(i6)!it turns the carbon content of the last day of the previous month on the carbon content of the 1st day of the next month(for all the PFTs)    	
		!cbwoodi_pft(i6) = cbwoodf_pft(i6)!it turns the carbon content of the last day of the previous month on the carbon content of the 1st day of the next month(for all the PFTs)
		!cstoi_pft(i6) = cstof_pft(i6)!it turns the carbon content of the last day of the previous month on the carbon content of the 1st day of the next month(for all the PFTs)
	    !crepi_pft(i6) = crepf_pft(i6)!it turns the carbon content of the last day of the previous month on the carbon content of the 1st day of the next month(for all the PFTs) 
        !cotheri_pft(i6) = cotherf_pft(i6)!it turns the carbon content of the last day of the previous month on the carbon content of the 1st day of the next month(for all the PFTs)       
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
	  mean_cleaf(i,j) = 0.0
	  mean_cawood(i,j) = 0.0
	  mean_cfroot(i,j) = 0.0
	  !mean_cbwood(i,j) = 0.0
	  !mean_csto(i,j) = 0.0
	  !mean_crep(i,j) = 0.0
	  !mean_cother(i,j) = 0.0
	  ctotal_tbe(i,j) = 0.0
	  ctotal_tbd(i,j) = 0.0
	  ctotal_herb(i,j) = 0.0
         do i6=1,npls
         mean_cleaf_pft(i,j,i6)=0.0
		 mean_cawood_pft(i,j,i6)=0.0
		 mean_cfroot_pft(i,j,i6)=0.0
		 !mean_cbwood_pft(i,j,i6)=0.0
		! mean_csto_pft(i,j,i6)=0.0
		 !mean_crep_pft(i,j,i6)=0.0
		 !mean_cother_pft(i,j,i6)=0.0
		 mean_npp_pft(i,j,i6)=0.0
		 mean_photo_pft(i,j,i6)=0.0
		 mean_ar_pft(i,j,i6)=0.0
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
		 mean_cleaf(i,j) = mean_cleaf(i,j) + (cleaf(i,j,k)/12)
		 mean_cawood(i,j) = mean_cawood(i,j) + (cawood(i,j,k)/12) 
		 mean_cfroot(i,j) = mean_cfroot(i,j) + (cfroot(i,j,k)/12)
		 !mean_cbwood(i,j) = mean_cbwood(i,j) + (cbwood(i,j,k)/12)
		 !mean_csto(i,j) = mean_csto(i,j) + (csto(i,j,k)/12)
		 !mean_crep(i,j) = mean_crep(i,j) + (crep(i,j,k)/12)
		 !mean_cother(i,j) = mean_cother(i,j) + (cother(i,j,k)/12)
		 
         if (npp(i,j,k).lt.nppmin(i,j)) nppmin(i,j) = npp(i,j,k)
         if (npp(i,j,k).gt.nppmax(i,j)) nppmax(i,j) = npp(i,j,k)
		    
c it calculates mean annual biomass for all the PFTs
			   do i6=1,npls
               mean_cleaf_pft(i,j,i6)=mean_cleaf_pft(i,j,i6)+
     &  		 (cleaf_pft(i,j,k,i6)/12)
	           mean_cawood_pft(i,j,i6)=mean_cawood_pft(i,j,i6)+
     &  		 (cawood_pft(i,j,k,i6)/12)
	           mean_cfroot_pft(i,j,i6)=mean_cfroot_pft(i,j,i6)+
     &  		 (cfroot_pft(i,j,k,i6)/12)
!	           mean_cbwood_pft(i,j,i6)=mean_cbwood_pft(i,j,i6)+
!     &  		 (cbwood_pft(i,j,k,i6)/12)
!	           mean_csto_pft(i,j,i6)=mean_csto_pft(i,j,i6)+
!     &  		 (csto_pft(i,j,k,i6)/12)
!	           mean_crep_pft(i,j,i6)=mean_crep_pft(i,j,i6)+
!     &  		 (crep_pft(i,j,k,i6)/12)
!	           mean_cother_pft(i,j,i6)=mean_cother_pft(i,j,i6)+
!     &  		 (cother_pft(i,j,k,i6)/12)
	           mean_npp_pft(i,j,i6)=mean_npp_pft(i,j,i6)+
     &  		 (npp_pft(i,j,k,i6)/12)
	           mean_photo_pft(i,j,i6)=mean_photo_pft(i,j,i6)+
     &  		 (photo_pft(i,j,k,i6)/12)
	           mean_ar_pft(i,j,i6)=mean_ar_pft(i,j,i6)+
     &  		 (ar_pft(i,j,k,i6)/12)
	        
	         enddo
      enddo
	  
c
c calculate PFTs total biomass
!      ctotal_tbe(i,j)=mean_cleaf_tbe(i,j)+mean_cawood_tbe(i,j)+
!     &	  mean_cbwood_tbe(i,j)+mean_cfroot_tbe(i,j)+mean_csto_tbe(i,j)+
!     &    mean_crep_tbe(i,j)+mean_cother_tbe(i,j)
     
!	  ctotal_tbd(i,j)=mean_cleaf_tbd(i,j)+mean_cawood_tbd(i,j)+
!     &	  mean_cbwood_tbd(i,j)+mean_cfroot_tbd(i,j)+mean_csto_tbd(i,j)+
!     &    mean_crep_tbd(i,j)+mean_cother_tbd(i,j)

!      ctotal_herb(i,j)=mean_cleaf_herb(i,j)+mean_cawood_herb(i,j)+
!     &	  mean_cbwood_herb(i,j)+mean_cfroot_herb(i,j)+
!     &    mean_csto_herb(i,j)+ mean_crep_herb(i,j)+mean_cother_herb(i,j)
	 
	     do i6=1,npls
	         ctotal_pft(i,j,i6)=mean_cleaf_pft(i,j,i6)+
     &	            mean_cawood_pft(i,j,i6)+mean_cfroot_pft(i,j,i6)
	     enddo
	 
	     ctotal_tbe(i,j)=mean_cleaf_tbe(i,j)+mean_cawood_tbe(i,j)+
     &	  mean_cfroot_tbe(i,j)
	       
         ctotal_tbd(i,j)=mean_cleaf_tbd(i,j)+mean_cawood_tbd(i,j)+
     &	  mean_cfroot_tbd(i,j)
	

         ctotal_herb(i,j)=mean_cleaf_herb(i,j)+mean_cawood_herb(i,j)+
     &	  mean_cfroot_herb(i,j)
	    
         do i6=1,npls
		    ctotal_gd(i,j)= ctotal_gd(i,j)+ctotal_pft(i,j,i6)
		    enddo	
		 
!		 ctotal_gd(i,j)=ctotal_tbe(i,j)+ctotal_tbd(i,j)+ctotal_herb(i,j)

		    do i6=1,npls      
		 ocp_pft(i,j,i6)= ctotal_pft(i,j,i6)/ctotal_gd(i,j)
		     if(ctotal_pft(i,j,i6).le.0.0) then
		        ocp_pft(i,j,i6)=0.0
		     endif
		    enddo
		
		 ocp_tbe(i,j)=ctotal_tbe(i,j)/ctotal_gd(i,j)
	      if(ctotal_tbe(i,j).le.0.0) then
		    ocp_tbe(i,j)=0.0
			   endif
			   
			 
	     ocp_tbd(i,j)=ctotal_tbd(i,j)/ctotal_gd(i,j)
	      if(ctotal_tbd(i,j).le.0.0) then
		    ocp_tbd(i,j)=0.0
			   endif
			  
			   
	      ocp_herb(i,j)=ctotal_herb(i,j)/ctotal_gd(i,j)
	      if(ctotal_herb(i,j).le.0.0) then
		    ocp_herb(i,j)=0.0
			   endif
			   
			   
		 ctotal_tbe(i,j)= ctotal_tbe(i,j)*ocp_tbe(i,j)
         ctotal_tbd(i,j)= ctotal_tbd(i,j)*ocp_tbd(i,j)
         ctotal_herb(i,j)= ctotal_herb(i,j)*ocp_herb(i,j)
		 
		    do i6=1,npls
	         ctotal_pft(i,j,i6)=ctotal_pft(i,j,i6)*ocp_pft(i,j,i6)
			 ctotal_pft(i,j,i6)=ctotal_pft(i,j,i6)
			 mean_cleaf_pft(i,j,i6)=mean_cleaf_pft(i,j,i6)*ocp_pft(i,j,i6)
			 mean_cawood_pft(i,j,i6)=mean_cawood_pft(i,j,i6)*ocp_pft(i,j,i6)
			 mean_cfroot_pft(i,j,i6)=mean_cfroot_pft(i,j,i6)*ocp_pft(i,j,i6)
               enddo
		 
!		 ctotal_gd(i,j)=ctotal_tbe(i,j)+ctotal_tbd(i,j)+ctotal_herb(i,j)


             do i6=1,npls
		    ctotal_gd(i,j)= 0
		    enddo	
		 
		    do i6=1,npls
		    ctotal_gd(i,j)= ctotal_gd(i,j)+ctotal_pft(i,j,i6)
		    enddo	
		
!         total_npp(i,j)= mean_npp_tbe(i,j) + mean_npp_tbd(i,j)+
!    &       mean_npp_herb(i,j)
!	     total_ar(i,j)= mean_ar_tbe(i,j) + mean_ar_tbd(i,j)+
 !    &       mean_npp_herb(i,j)
! total_npp(i,j)= mean_photo_tbe(i,j) + mean_photo_tbd(i,j)+
 !    &       mean_photo_herb(i,j)

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
     &                     cleaf1,cawood1,cfroot1,                      !input
     &                     cbwood1,csto1,crep1,cother1,	                !input
     &                     cleafi_pft,cawoodi_pft,
     &    				   cbwood1_pft,cfrooti_pft,
     &	                   csto1_pft,crep1_pft,
     &    				   cother1_pft,
     &                     cleafi_tbe,cawoodi_tbe,cfrooti_tbe,          !input
     &                     cbwoodi_tbe,cstoi_tbe,crepi_tbe,cotheri_tbe, !input
     &                     cleafi_tbd,cawoodi_tbd,cfrooti_tbd,          !input
     &                     cbwoodi_tbd,cstoi_tbd,crepi_tbd,cotheri_tbd, !input
     &                     cleafi_herb,cawoodi_herb,cfrooti_herb,       !input
     &                  cbwoodi_herb,cstoi_herb,crepi_herb,cotheri_herb, !input
     &                     cleaf1_pft,cawood1_pft,cfroot1_pft,
     &                     w2,g2,s2,smavg,ruavg,evavg,	 	            !output
     &                     epavg,phavg,aravg,nppavg,laiavg,	            !output
     &                     rmlavg,rmfavg,rmsavg,rmavg,					!output
     &                     rglavg,rgfavg,rgsavg,rgavg,					!output
     &                     clavg,csavg,hravg,rcavg,                     !output
     &                     cleaf2,cawood2,cfroot2,                      !output
     &                     cbwood2,csto2,crep2,cother2,                 !output
     &                     cleafavg,cawoodavg,cfrootavg,                !output
     &                     cbwoodavg,cstoavg,crepavg,cotheravg,         !output
     &                     cleaf2_pft,cawood2_pft,
     &					   cbwood2_pft,cfroot2_pft,
     &	 				   csto2_pft,crep2_pft,
     &					   cother2_pft,
     &                     cleafavg_pft,cawoodavg_pft,
     &                     cbwoodavg_pft,cfrootavg_pft,
     &                     cstoavg_pft,crepavg_pft,
     &                     cotheravg_pft,nppavg_pft,photoavg_pft,
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
      parameter (npft=3,npls=3,mal=80,maf=80,maw=80,mtl=2)!npls:number of pfts; mal:allocation maximum value for leaf compartment;maf: llocation maximum value for fine roots compartment
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
      real cleaf1,cawood1,cfroot1,cbwood1,csto1,
     &     crep1,cother1,
     &     cleaf1_pft(npls),cawood1_pft(npls),cfroot1_pft(npls),
     &     cleafi_pft(npls),cawoodi_pft(npls),
     &     cbwood1_pft(npls),cfrooti_pft(npls),
     &     csto1_pft(npls),crep1_pft(npls),
     &     cother1_pft(npls),
     &     cleaf2,cawood2,cfroot2,cbwood2,csto2,
     &     crep2,cother2,
     &     cleaf2_pft(npls),cawood2_pft(npls),
     &     cbwood2_pft(npls),cfroot2_pft(npls),
     &	   csto2_pft(npls),crep2_pft(npls),
     &     cother2_pft(npls),
     &     cleafavg,cawoodavg,cfrootavg,
     &     cbwoodavg,cstoavg,crepavg,cotheravg,
     &     cleafavg_pft(npls),cawoodavg_pft(npls),
     &     cbwoodavg_pft(npls),cfrootavg_pft(npls),
     &	   cstoavg_pft(npls),crepavg_pft(npls),
     &     cotheravg_pft(npls), 
     &     cleafi_tbe,cawoodi_tbe,cfrooti_tbe,
     &     cbwoodi_tbe,cstoi_tbe,crepi_tbe,cotheri_tbe,
     &     cleafi_tbd,cawoodi_tbd,cfrooti_tbd,
     &     cbwoodi_tbd,cstoi_tbd,crepi_tbd,cotheri_tbd,
     &     cleafi_herb,cawoodi_herb,cfrooti_herb,
     &     cbwoodi_herb,cstoi_herb,crepi_herb,cotheri_herb,	 
     &     cl1,cl2,ca1,ca2,cf1,cf2,cb1,cb2,cs1,cs2,
     &     cr1,cr2,co1,co2,
     &     cl2_pft(npls),ca2_pft(npls),cb2_pft(npls),
     &     cf2_pft(npls),cs2_pft(npls),cr2_pft(npls),co2_pft(npls),
     &     alfa_leaf(npls),alfa_awood(npls),alfa_bwood(npls),
     &     alfa_froot(npls),
     &     beta_leaf, beta_awood,beta_bwood,beta_froot,
     &     ctotal_pft2(npls),ctotal_tbe2,ctotal_tbd2,ctotal_herb2,!Daily total biomass for each PFT 
     &     ctotal_gd2,      !Daily total biomass in a grid cell 
     &     ocp_tbe2,ocp_tbd2,ocp_herb2, !PFTs proportional occupation in a grid cell(pft_biomass/gridcell_total_biomass) 
     &     ctotali_tbe2,ctotali_tbd2, ctotali_herb2,!Initial daily total biomass for each PFT 
     &     ctotali_gd,     !Initial Daily total biomass in a grid cell 
     &     ocpi_tbe,ocpi_tbd,ocpi_herb, !PFTs proportional occupation in a grid cell(pft_biomass/gridcell_total_biomass)for the 1st day of the yr 
     &     cawood_tbe2,cawood_tbd2,cawood_herb2, !Daily aboveground woody biomass for each PFT
     &     cawood_gd,	 ! Aboveground woody biomass in a grid cell
     &     ocp_wood(npls),
     &     max_value, !to find the greatest value of ocp_wood
     &     max_index, !to find the index of the greatest value of ocp_wood
     &     cawoodi_gd, ! Initial Aboveground woody biomass in a grid cell
     &     ocpi_wood_tbe,ocpi_wood_tbd,ocpi_wood_herb,maxi_index,
     &     ocp_pft2(npls)
      real ocpi_wood_pft(npls)
      real ctotali_pft(npls),ocpi_pft(npls)
      real alc_leaf(mal,maf,maw),alc_froot(mal,maf,maw),
     &     alc_awood(mal,maf,maw)
	 
	 

					 
				 
	  
	  
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
      cleafavg = 0.
	  cawoodavg = 0.
	  cfrootavg = 0.
!	  cbwoodavg = 0.
!	  cstoavg = 0.
!	  crepavg = 0.
!	  cotheravg = 0.
	  ctotal_gd2=0. !it calculates the total biomass of a grid  cell
	  cawood_gd = 0. !it calculates the total aboveground woody biomass in a grid cell
	  cawoodi_gd = 0. !it calculates the initial total aboveground woody biomass in a grid cell
	  ctotali_gd=0
	     max_value= -100.
		 
	     do i6 = 1,npls
		    cleafavg_pft(i6) = 0.!mean monthly leaf biomass for all the PFTs
			cawoodavg_pft(i6) = 0.!mean monthly aboveground biomass for all the PFTs
			cfrootavg_pft(i6) = 0.!mean monthly belowground biomass for all the PFTs
!			cbwoodavg_pft(i6) = 0.!mean monthly belowground biomass for all the PFTs
!			cstoavg_pft(i6) = 0.!mean monthly leaf biomass for all the PFTs
!			crepavg_pft(i6) = 0.!mean monthly aboveground biomass for all the PFTs
!			cotheravg_pft(i6) = 0.!mean monthly belowground biomass for all the PFTs
			nppavg_pft(i6) = 0. !mean monthly npp for all the PFTs
			photoavg_pft(i6) = 0. !mean monthly photosynthesis for all the PFTs
			aravg_pft(i6) = 0. !mean monthly autotrophic respiration for all the PFTs
			
			cl2_pft(i6)=cleafi_pft(i6)
			ca2_pft(i6)=cawoodi_pft(i6)
			cf2_pft(i6)=cfrooti_pft(i6)
!			cb2_pft(i6)=cbwood1_pft(i6)
!			cs2_pft(i6)=csto1_pft(i6)
!			cr2_pft(i6)=crep1_pft(i6)
!			co2_pft(i6)=cother1_pft(i6)
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
				    alc_leaf(i1,i2,i3)=i1
					alc_froot(i1,i2,i3)=i2
				    alc_awood(i1,i2,i3)=i3
					
					i6=i6+1
				   			
			 
			     cl1 = cl2_pft(i6) !transforma o valor do dia anterior no valor atual
			      
				 ca1 = ca2_pft(i6)
				 cf1 = cf2_pft(i6)
				
!				 cb1 = cb2_pft(i6)
!				 cs1 = cs2_pft(i6) !transforma o valor do dia anterior no valor atual
!				 cr1 = cr2_pft(i6)
!				 co1 = co2_pft(i6)
				 beta_leaf = alfa_leaf(i6)
				 beta_awood = alfa_awood(i6)
				 beta_froot = alfa_froot(i6)
!				  beta_bwood = alfa_bwood(i6)
				 
	         if ((i.eq.1).and.(month.eq.1)) then
				     cl1 = cleaf1_pft(i6)
					 ca1 = cawood1_pft(i6)
					 cf1 = cfroot1_pft(i6)
!	                 cb1 = cbwoodi_tbe
!					 cs1 = csto_tbe
!		             co1 = cotheri_tbe
!		             cr1 = crepi_tbe
					 beta_leaf=0.
					 beta_awood=0.
!					 beta_bwood=0.
					 beta_froot=0.
					 ocp_pft2(i6)=ocpi_pft(i6)
					 max_index=maxi_index
			    endif
				  
                
				
c carbon cycle (photosynthesis, plant respiration and NPP)
      call carbon1 (temp,p0,w,wmax,ca,ipar,i6,tsoil,emax,rc2,  !input
     &               ca2,cf2,cb2,cl1,ca1,cf1,cb1,
     &               max_index,
     &               ocp_pft2,
     &               beta_leaf,beta_awood,beta_bwood,beta_froot,          !input
     &              ph,ar,nppa,laia,f5,rm,rml,rmf,rms,rg,rgl,rgf,rgs)		        !output

	          
	 
!     carbon allocation (carbon content on each compartment)     
       call allocation (nppa,cl1,ca1,cf1,cb1,cs1,cr1,co1,i6,	   !input
     &                 alc_leaf,alc_froot,alc_awood,
     &                  cl2,ca2,cf2,cb2,cs2,cr2,co2)  !output 

	            alfa_leaf(i6)  =cl2 - cl1 !it calculates the difference btween the biomass of the last and actual day (for growth respiration)
				alfa_awood(i6) =ca2 - ca1 !it calculates the difference btween the biomass of the last and actual day (for growth respiration)  
			    alfa_bwood(i6) =cb2 - cb1 !it calculates the difference btween the biomass of the last and actual day (for growth respiration) 
				alfa_froot(i6) =cf2 - cf1 !it calculates the difference btween the biomass of the last and actual day (for growth respiration) 
				
				cl2_pft(i6) = cl2
				ca2_pft(i6) = ca2
				cf2_pft(i6) = cf2
	          
!				cb2_pft(i6) = cb2
!				cs2_pft(i6) = cs2
!				cr2_pft(i6) = cr2
!				co2_pft(i6) = co2
				
				  npp2_pft(i6)=nppa
				  photo2_pft(i6)=ph
				  ar2_pft(i6)=ar
				


			ctotal_pft2(i6)=cl2_pft(i6)+ca2_pft(i6)+cf2_pft(i6)
			  
                  ctotal_gd2 = ctotal_gd2 + ctotal_pft2(i6)
				  cawood_gd = cawood_gd + ca2_pft(i6)
				 
				   endif
				  enddo
				  enddo
				  enddo
            			  
		     
				  
   				do i6=1, npls
				ocp_pft2(i6) = ctotal_pft2(i6)/ctotal_gd2
				   if(ctotal_pft2(i6).le.0.)then
                       ocp_pft2(i6) = 0.
					  endif		  
                 enddo					  

			 

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
      cleafavg = cleafavg + cl2 !KgC/m2
	  cawoodavg = cawoodavg + ca2
	  cfrootavg = cfrootavg + cf2
!	  cbwoodavg = cbwoodavg + cb2
!	  cstoavg = cstoavg + cs2
!	  crepavg = crepavg + cr2
!	  cotheravg = cotheravg + co2	  
       do i6 = 1,npls
	       cleafavg_pft(i6) = cleafavg_pft(i6) + cl2_pft(i6)
		   cawoodavg_pft(i6) = cawoodavg_pft(i6) + ca2_pft(i6)
		   cfrootavg_pft(i6) = cfrootavg_pft(i6) + cf2_pft(i6)
!		   cbwoodavg_pft(i6) = cbwoodavg_pft(i6) + cb2_pft(i6)
!		   cstoavg_pft(i6) = cstoavg_pft(i6) + cs2_pft(i6)
!		   crepavg_pft(i6) = crepavg_pft(i6) + cr2_pft(i6)
!		   cotheravg_pft(i6) = cotheravg_pft(i6) + co2_pft(i6)
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
!	  cbwood2_pft(i6) = cb2_pft(i6)
!	   csto2_pft(i6) = cs2_pft(i6)
!	   crep2_pft(i6) = cr2_pft(i6)
!	  cother2_pft(i6) = co2_pft(i6)
	      enddo
      w2 = w
      g2 = g
      s2 = s
      smavg = smavg/real(ndmonth(month))
      ruavg = ruavg/real(ndmonth(month))
      evavg = evavg/real(ndmonth(month))
      epavg = epavg/real(ndmonth(month))
      rcavg = rcavg/real(ndmonth(month))
      cleafavg = cleafavg/real(ndmonth(month)) !monthly carbon content on leaf compart
	  cawoodavg = cawoodavg/real(ndmonth(month))
	  cfrootavg = cfrootavg/real(ndmonth(month))
!	  cbwoodavg = cbwoodavg/real(ndmonth(month))
!	  cstoavg = cstoavg/real(ndmonth(month))
!	  crepavg = crepavg/real(ndmonth(month))
!	  cotheravg = cotheravg/real(ndmonth(month))
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
!	  cbwoodavg_pft(i6) = cbwoodavg_pft(i6)/real(ndmonth(month))
!	  cstoavg_pft(i6) = cstoavg_pft(i6)/real(ndmonth(month))
!	  crepavg_pft(i6) = crepavg_pft(i6)/real(ndmonth(month))
!	  cotheravg_pft(i6) = cotheravg_pft(i6)/real(ndmonth(month))
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

