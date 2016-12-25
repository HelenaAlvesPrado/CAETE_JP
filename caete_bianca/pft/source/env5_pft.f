c234567
      program env
c
c=======================================================================
c CPTEC-PVM2
c Adapted from env4.f. 
c Code written by David Lapola
c Last update: Aug/2007
c
c Compile with: g95 (or gfortran) env5.f wbm4.f carbon2.f -o a.exe (a.exe is the executable)
c Execute with ./a.exe
c Then run g95 (gfortran) pvm5.f -o pvm5.exe
c Execute with ./pvm5.exe
c
c=======================================================================
c
c Parameters and variables
c ------------------------
c
c     implicit none
      integer nx,ny,ni
      integer cell_id(6003),var_id
      parameter(nx=192,ny=96,a=360,b=180,npls=3)
      real lsmk(nx,ny),p0(nx,ny)
      real pr(nx,ny,12),t(nx,ny,12)
      real prec(nx,ny,12),temp(nx,ny,12),par(nx,ny,12),ipar(nx,ny,12),ca
      real anpr(nx,ny,12),ant(nx,ny,12)
      real tmin(nx,ny),seanpp(nx,ny),meanpp(nx,ny),meanhr(nx,ny),
     & meancs(nx,ny),mphoto(nx,ny),maresp(nx,ny),ave_wsoil(nx,ny),
     & ave_evap(nx,ny),ave_rc(nx,ny),nppot(nx,ny),nppot_aux
      
      real bpot(nx,ny)
      real wsoil2(nx,ny,12),evaptr(nx,ny,12),rc(nx,ny,12),waux(nx,ny),
     & waux2(nx,ny),waux3(nx,ny),waux4(nx,ny),waux5(nx,ny),waux6(nx,ny)
      real long(nx,ny),lat(nx,ny),lat_graus,long_graus,aux_lat(6003),
     & aux_long(6003)
      real trans_month(nx,ny),trans_month2(nx,ny)
      real aux_des(nx,ny),aux_des2(nx,ny)
      real meanrml(nx,ny), meanrmf(nx,ny), meanrms(nx,ny),meanrm(nx,ny),
     &     meanrgl(nx,ny), meanrgf(nx,ny), meanrgs(nx,ny), meanrg(nx,ny)
c carbon cycle
      real mp(nx,ny,12),mar(nx,ny,12),mnpp(nx,ny,12),
     &     lai(nx,ny,12),clit(nx,ny,12),csoil(nx,ny,12),
     &     hresp(nx,ny,12)
      real monrml(nx,ny,12),monrmf(nx,ny,12),monrms(nx,ny,12), 
     &     monrm(nx,ny,12), monrgl(nx,ny,12),monrgf(nx,ny,12),
     &     monrgs(nx,ny,12), monrg(nx,ny,12)
      real mean_npp_pft(nx,ny,npls), mean_ar_pft(nx,ny,npls),
     &     mean_photo_pft(nx,ny,npls),
     &     total_npp(nx,ny),total_ar(nx,ny),total_ph(nx,ny)
c carbon allocation
      real cleafini_pft(nx,ny,npls),cawoodini_pft(nx,ny,npls),
     & cfrootini_pft(nx,ny,npls),
     & cleafini_pft2(npls),cawoodini_pft2(npls),cfrootini_pft2(npls)
      real cleaf(nx,ny,12),cawood(nx,ny,12),cfroot(nx,ny,12),
     &      waux7(nx,ny),waux8(nx,ny),waux9(nx,ny),
     &     waux10(nx,ny),waux11(nx,ny),waux12(nx,ny),waux13(nx,ny),
     &     waux14(nx,ny),waux15(nx,ny),waux16(nx,ny),waux17(nx,ny),
     &     waux18(nx,ny),waux19(nx,ny),waux20(nx,ny),waux21(nx,ny),
     &     waux22(nx,ny),waux23(nx,ny),waux24(nx,ny),waux25(nx,ny),
     &     waux26(nx,ny),waux27(nx,ny),waux28(nx,ny),waux29(nx,ny),
     &     waux30(nx,ny), waux31(nx,ny),waux32(nx,ny),
     &     mean_cleaf_pft(nx,ny,npls),mean_cawood_pft(nx,ny,npls),
     &     mean_cfroot_pft(nx,ny,npls),
     & ctotal_gd(nx,ny),ctotal_pft(nx,ny,npls)
      integer i6 ! index for calculate the 3 PFTs (i4.eq.1 - TBE; i4.eq.2 - TBD; i4.eq.3 - HERB)
    
c spinup
      integer supindex,i1,i2,i3,i10 
      real alc_leaf2,alc_froot2,alc_awood2

c Open files
c ----------


      open( 9,file='../inputs/psup.bin',status='old',
     &        form='unformatted',access='direct',recl=4*nx*ny)
      open(10,file='../inputs/lsmk.bin',status='old',
     &        form='unformatted',access='direct',recl=4*nx*ny)
      open(11,file='../inputs/prec.bin',status='old',
     &        form='unformatted',access='direct',recl=4*nx*ny)
      open(12,file='../inputs/temp.bin',status='old',
     &        form='unformatted',access='direct',recl=4*nx*ny)
      open(13,file='../inputs/ipar.bin',status='old',
     &        form='unformatted',access='direct',recl=4*nx*ny)
      open(20,file='../inputs/nppot.bin',status='old',
     &        form='unformatted',access='direct',recl=4*nx*ny)

! Future anomalies, see lines 85 and 86
!      open(14,file='../ipcc_ar4/anomalias/SRESA2/temp/'//
!     &               'HADCM3_A2_2070_2099_anom_t.bin',status='old',
!     &        form='unformatted',access='direct',recl=4*nx*ny)
!      open(15,file='../ipcc_ar4/anomalias/SRESA2/prec/'//
!     &               'HADCM3_A2_2070_2099_anom_pr.bin',status='old',
!     &        form='unformatted',access='direct',recl=4*nx*ny)
!
! Read data
c ---------

      read (9,rec=1) p0     !mb
      read (10,rec=1) lsmk   !land = 1; ocean = 0
      read (20,rec=1) nppot
      call read12 (11,pr)  !mm/month
      call read12 (12,t)  !oC
      call read12 (13,ipar)  !w/m2
!      call read12 (14,ant)  !oC
!      call read12 (15,anpr)  !oC

c Close files
c ---------------

      close (9)
      close (10)
      close (11)
      close (12)
      close (13)
      close (20)
!      close (14)
!      close (15)
c
c
c--------------------------------------------------------------------------------------------
!      do i=1,nx
!         do j=1,ny
!
!      if (i.le.96) then     !inverte hemisf‚rios leste-oeste
!         ni=i+96
!         else
!         ni=(96-i)*(-1)
!      endif
!
!      nj=(j-96)*(-1)    !faz flip na latitude
!
!     pO(i,j)=p0(ni,nj)
!     lsmk(i,j)=lsmk(ni,nj)
!     aux_des2(i,j)=aux_des(ni,nj) !Aux_des com nova configura‡Æo de mapa!
!     trans_month2(i,j)=trans_month(ni,nj) !trans_month com nova configura‡Æo de mapa!
!
!         enddo
!      enddo
c----------------------------------------------------------------------------------------------
c
c
      var_id=0 !define uma vari vel para contabilizar cell_id
      prec_t=0.0 !define a vari vel precipita‡Æo_total
c
c
!      do i=1,nx
!         do j=1,ny
!      aux_des(i,j)=0 !define uma vari vel para contabilizar os meses com precipita‡Æo inferior a 100mm
!      p0(i,j)=1010.0
 !        enddo
 !     enddo
c
c
      do k=1,12
         do i=1,nx
         do j=1,ny
		 

c --- Photosynthetically active radiation reaching canopy --------------
c     (ipar ; Ein/m2/s) [Eq. 7] observed data from ISLSCP2
      par(i,j,k) = ipar(i,j,k)/(2.18e5) !converting to Ein/m2/s
      temp(i,j,k) = t(i,j,k) !+ant(i,j,k) !uncomment to use future anomalies
      prec(i,j,k) = pr(i,j,k)!+anpr(i,j,k) !+pr(i,j,k)*0.2 !uncomment to use future anomalies
      if (prec(i,j,k).lt.0.0) prec (i,j,k) = 0.0

c!!!!!!!!!!!!! CA NÃO PRECISA ESTAR DENTRO DO LOOP

c --- Atmospheric CO2 pressure (Pa) ! ppmv / Pa ------------------------
c     1 Pa CO2 = 9.901 ppmv CO2 (Adams et al. 2004)
!      ca= 18.18 !Pa (=180 ppmv; Last Glacial Maximum)
!      ca= 28.28 !Pa (=280 ppmv; Pre-Industrial Rev.)
!      ca= 35.35 !Pa (=350 ppmv; 1961-1990)
      ca= 400/9.901 !Pa (=350 ppmv; 1961-1990)
!      ca= 54.03 !Pa (=535 ppmv; SRES-B1 2080's)
!      ca= 73.73 !Pa (=730 ppmv; SRES-A2 2080's)
!      ca= ((73.73-35.35)*0.5)+35.35 !Pa half effect!!!

c
c     !exerc¡cio: dura‡Æo da esta‡Æo seca ------------------------------

!      if(prec(i,j,k).lt.100.0) then
!      aux_des(i,j)=aux_des(i,j)+1.0
!      else
!      aux_des(i,j)=aux_des(i,j)+0.0
!      endif
c
c
c     !exerc¡cio: mˆs de transi‡Æo (esta‡Æo seca --> chuvosa) ----------

      if (prec(i,j,k-1).lt.100.and.prec(i,j,k-2).lt.100
     & .and.prec(i,j,k-3).lt.100.and.prec(i,j,k).gt.100)then
      trans_month(i,j)=k-1
      endif
c
      if (trans_month(i,j).eq.0) then
      trans_month(i,j)=12
      endif

c
c     !exercicio: valor de inputs para c‚lulas especificas -------------

      cellx=-60.06 !longitude
      celly=-2.35 !latitude
      long_graus=a/nx !extensÆo longitudinal de cada c‚lula em graus
      lat_graus=b/ny  !extensÆo latitudinal de cada c‚lula em graus

      long(i,j)=long_graus*(i) !coloca em coordenada "0-360"
        if (long(i,j).lt.180) then
        long(i,j)=long(i,j) !transforma para coordenada "-180|180"
        else
        long(i,j)=(360-long(i,j))*(-1)
        endif

      lat(i,j)=(lat_graus*(j)) !coloca em coordenada "0-180'
      lat(i,j)=(lat(i,j)-90.0)


c     if (k.eq.1.and.lat(i,j).ge.celly.and.lat(i,j).lt.(celly+0.5)!para determinado mˆs
c     & .and.long(i,j).ge.cellx.and.long(i,j).lt.(cellx+0.5)) then

c      if (lat(i,j).ge.celly.and.lat(i,j).lt.(celly+0.5) !para o ano inteiro
c     & .and.long(i,j).ge.cellx.and.long(i,j).lt.(cellx+0.5)) then
c      prec_t=prec_t+prec(i,j,k)


c     print*,ipar(i,j,k),prec(i,j,k),t(i,j,k),ave_evap(i,j)
c     endif
         enddo
         enddo
      enddo
c
c
c ------------------------------------------------------
c  Spinup to calculate initial carbon content for the compartments
c -----------------------------------------------------
c


      do i=1,nx
	     do j=1,ny
		    if(lsmk(i,j).eq.1)then
			   nppot_aux = nppot(i,j)
			
!		    print*, nppot(i,j),nppot_aux

	!!! NAO ENTENDO ESTAS VARIAVESI ALC_LEAF2 ETC ... NÃO SAO UTILIZADAS 
    !!! NA SPINUP... 		    

       call spinup (nppot_aux, 
     &          cleafini_pft2,cawoodini_pft2,cfrootini_pft2)
	             
				   do i1=1,npls
				   cleafini_pft(i,j,i1)=cleafini_pft2(i1)
				   cawoodini_pft(i,j,i1)=cawoodini_pft2(i1)
				   cfrootini_pft(i,j,i1)=cfrootini_pft2(i1)
				   enddo
	           endif
			   enddo
			   enddo

         
c -------------------------------------------------------
c Calculate environmental variables (water balance model)
c -------------------------------------------------------
c
      call wbm (prec,temp,lsmk,p0,ca,par,
     &          cleafini_pft,cawoodini_pft,cfrootini_pft,
     &          tmin,meanpp,seanpp,mphoto,maresp,  !carbon
     &          meanrml,meanrmf,meanrms,meanrm,
     &          meanrgl,meanrgf,meanrgs,meanrg,
     &          meanhr,meancs,wsoil2,evaptr,mnpp,
     &          mp,mar,rc,ave_wsoil,ave_evap,ave_rc,
     &          monrml,monrmf,monrms,monrm,
     &          monrgl,monrgf,monrgs,monrg,
     &          ctotal_gd,
     &          total_npp,total_ar,total_ph,
     &          ctotal_pft, mean_npp_pft, mean_ar_pft,
     &          mean_photo_pft,
     &          mean_cleaf_pft,mean_cawood_pft,mean_cfroot_pft)

     c     !exerc¡cio: nomear as celulas de grad ----------------------------

       do i=1,nx
         do j=1,ny

      if (lsmk(i,j).eq.1) then
         var_id = var_id + 1.0
         cell_id(var_id)= var_id
         aux_lat(var_id)=lat(i,j)
         aux_long(var_id)=long(i,j)
		 
	      if (cell_id(var_id) .eq. 5244)  then

      print*,'ctoal_pft', ctotal_pft(i,j,1), ctotal_pft(i,j,2),
     &        ctotal_pft(i,j,3) 
c  100 format (A8,1x,I4,1x,A8,1x,F7.2,1x,A9,f7.2)

      endif
      endif
           enddo
      enddo
	  
c verificar se as variáveis com dimensão i6 estão saindo da subrotina wbm com 3 valores 
       do i=1,nx
	      do j=1,ny
		     if (lsmk(i,j).eq.1) then
			   do i6=1,npls
!			     print*, i6, 'ctotal',ctotal_pft(i,j,i6),
!     &          'npp', mean_npp_pft(i,j,i6), 
!     &	        'ar', mean_ar_pft(i,j,i6),
!     &          'ph', mean_photo_pft(i,j,i6) 
			    enddo
				   endif
			    enddo
			    enddo

c
c write environmental variables
      open(16,file='../outputs/ambientais4_tbe.bin',
!      open(16,file='../ipcc_env_vars2/'//
!     & 'ambientais4_HADCM3_A2_onlyclimate.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      write(16,rec=1) tmin   !oC
      write(16,rec=2) meanpp !kgC/m2/y
      write(16,rec=3) seanpp !dimensionless
      write(16,rec=4) meanhr !kgC/m2/y
      write(16,rec=5) meancs !kgC/m2/y
      write(16,rec=6) mphoto !kgC/m2/y
      write(16,rec=7) maresp !kgC/m2/y
      write(16,rec=8) ave_wsoil !mm
      write(16,rec=9) ave_evap !mm/d/y
      write(16,rec=10) ave_rc !s/m
      write(16,rec=11) meanrml
      write(16,rec=12) meanrmf
      write(16,rec=13) meanrms
   	  write(16,rec=14) meanrm
      write(16,rec=15) meanrgl
      write(16,rec=16) meanrgf
      write(16,rec=17) meanrgs
   	  write(16,rec=18) meanrg
      write(16,rec=19) mean_cleaf
      close(16)
	  
	  
	  
c annual carbon content on compartments
       open(18,file='../outputs/annual_biomass.bin',
!      open(16,file='../ipcc_env_vars2/'//
!     & 'ambientais4_HADCM3_A2_onlyclimate.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      write(18,rec=1) mean_cleaf  !oC
      write(18,rec=2) mean_cawood  !oC
      write(18,rec=3) mean_cbwood  !oC
      write(18,rec=4) mean_cfroot  !oC
      write(18,rec=5) mean_csto  !oC
      write(18,rec=6) mean_crep  !oC
      write(18,rec=7) mean_cother  !oC
      write(18,rec=8) mean_cleaf_tbe  !oC
      write(18,rec=9) mean_cleaf_tbd  !oC
      write(18,rec=10) mean_cawood_tbe  !oC
      write(18,rec=11) mean_cawood_tbd  !oC
      write(18,rec=12) mean_cbwood_tbe  !oC
      write(18,rec=13) mean_cbwood_tbd  !oC
      write(18,rec=14) mean_cfroot_tbe  !oC
      write(18,rec=15) mean_cfroot_tbd  !oC
      write(18,rec=16) mean_csto_tbe  !oC
      write(18,rec=17) mean_csto_tbd  !oC
      write(18,rec=18) mean_crep_tbe  !oC
      write(18,rec=19) mean_crep_tbd  !oC
      write(18,rec=20) mean_cother_tbe  !oC
      write(18,rec=21) mean_cother_tbd  !oC
      write(18,rec=22) mean_cleaf_herb !oC
      write(18,rec=23) mean_cawood_herb  !oC
      write(18,rec=24) mean_cbwood_herb  !oC
      write(18,rec=25) mean_cfroot_herb  !oC
      write(18,rec=26) mean_csto_herb  !o
      write(18,rec=27) mean_crep_herb  !oC
      write(18,rec=28) mean_cother_herb  !oC
      write(18,rec=29) ctotal_tbe  !oC
      write(18,rec=30) ctotal_tbd  !oC
      write(18,rec=31) ctotal_herb  !oC
      write(18,rec=32) ctotal_gd  !oC
      close(18)
	  
c anuual carbon fluxes
       open(18,file='../outputs/carbon_fluxes.bin',
!      open(16,file='../ipcc_env_vars2/'//
!     & 'ambientais4_HADCM3_A2_onlyclimate.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      write(18,rec=1) mean_npp_tbe !oC
      write(18,rec=2) mean_npp_tbd !oC
      write(18,rec=3) mean_npp_herb  !oC
      write(18,rec=4) mean_photo_tbe !oC
      write(18,rec=5) mean_photo_tbd !oC
      write(18,rec=6) mean_photo_herb  !oC
      write(18,rec=7) mean_ar_tbe !oC
      write(18,rec=8) mean_ar_tbd !oC
      write(18,rec=9) mean_ar_herb  !oC
      write(18,rec=10) total_npp !oC
      write(18,rec=11) total_ar !oC
      write(18,rec=12) total_ph  !oC
      close(18)	  


c write carbon content on leaf compartment for each PFT
      open(92,file='../outputs/mean_cleaf_pft.bin',
!      open(17,file='../ipcc_env_vars2/'//
!     & 'soilm_HADCM3_A2_onlyclimate.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do i6=1,npls
         do i=1,nx
         do j=1,ny
            waux30(i,j) = mean_cleaf_pft(i,j,i6)
         enddo
         enddo
         write(92,rec=i6) waux30
      enddo
      close(92)
	  
c write carbon content on awood compartment for each PFT
      open(92,file='../outputs/mean_cawood_pft.bin',
!      open(17,file='../ipcc_env_vars2/'//
!     & 'soilm_HADCM3_A2_onlyclimate.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do i6=1,npls
         do i=1,nx
         do j=1,ny
            waux31(i,j) = mean_cawood_pft(i,j,i6)
         enddo
         enddo
         write(92,rec=i6) waux31
      enddo
      close(92)

c write carbon content on awood compartment for each PFT
      open(92,file='../outputs/mean_cfroot_pft.bin',
!      open(17,file='../ipcc_env_vars2/'//
!     & 'soilm_HADCM3_A2_onlyclimate.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do i6=1,npls
         do i=1,nx
         do j=1,ny
            waux32(i,j) = mean_cfroot_pft(i,j,i6)
         enddo
         enddo
         write(92,rec=i6) waux32
      enddo
      close(92)	 	  
	  
	  
	  
c write total biomass for each PFT
      open(92,file='../outputs/ctotal_pft.bin',
!      open(17,file='../ipcc_env_vars2/'//
!     & 'soilm_HADCM3_A2_onlyclimate.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do i6=1,npls
         do i=1,nx
         do j=1,ny
            waux26(i,j) = ctotal_pft(i,j,i6)
         enddo
         enddo
         write(92,rec=i6) waux26
      enddo
      close(92)
	  
	  
c write mean annual npp for each PFT
      open(92,file='../outputs/mean_npp_pft.bin',
!      open(17,file='../ipcc_env_vars2/'//
!     & 'soilm_HADCM3_A2_onlyclimate.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do i6=1,npls
         do i=1,nx
         do j=1,ny
            waux27(i,j) = mean_npp_pft(i,j,i6)
         enddo
         enddo
         write(92,rec=i6) waux27
      enddo
      close(92)
	  
c write mean annual autotrophic respiration for each PFT
      open(92,file='../outputs/mean_ar_pft.bin',
!      open(17,file='../ipcc_env_vars2/'//
!     & 'soilm_HADCM3_A2_onlyclimate.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do i6=1,npls
         do i=1,nx
         do j=1,ny
            waux28(i,j) = mean_ar_pft(i,j,i6)
         enddo
         enddo
         write(92,rec=i6) waux28
      enddo
      close(92)	  

c write mean annual photosynthesis for each PFT
      open(92,file='../outputs/mean_photo_pft.bin',
!      open(17,file='../ipcc_env_vars2/'//
!     & 'soilm_HADCM3_A2_onlyclimate.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do i6=1,npls
         do i=1,nx
         do j=1,ny
            waux29(i,j) = mean_photo_pft(i,j,i6)
         enddo
         enddo
         write(92,rec=i6) waux29
      enddo
      close(92)	  
	  
c write wsoil2
      open(17,file='../outputs/soilm.bin',
!      open(17,file='../ipcc_env_vars2/'//
!     & 'soilm_HADCM3_A2_onlyclimate.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do k=1,12
         do i=1,nx
         do j=1,ny
            waux(i,j) = wsoil2(i,j,k)
         enddo
         enddo
         write(17,rec=k) waux
      enddo
      close(17)
c
c write evapotranspiration
      open(18,file='../outputs/evaptr.bin',
!      open(18,file='../ipcc_env_vars2/'//
!     & 'evaptr_HADCM3_A2_onlyclimate.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do k=1,12
         do i=1,nx
         do j=1,ny
            waux2(i,j) = evaptr(i,j,k)
         enddo
         enddo
         write(18,rec=k) waux2
      enddo
      close(18)
c
c write monthly NPP
      open(19,file='../outputs/mon_npp.bin',
!      open(19,file='../ipcc_env_vars2/'//
!     & 'mon_npp_HADCM3_A2_onlyclimate.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do k=1,12
         do i=1,nx
         do j=1,ny
            waux3(i,j) = mnpp(i,j,k)
         enddo
         enddo
         write(19,rec=k) waux3
      enddo
      close(19)
c
c write monthly Photosynthesis
      open(20,file='../outputs/mon_photo.bin',
!      open(20,file='../ipcc_env_vars2/'//
!     & 'mon_photo_HADCM3_A2_onlyclimate.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do k=1,12
         do i=1,nx
         do j=1,ny
            waux4(i,j) = mp(i,j,k)
         enddo
         enddo
         write(20,rec=k) waux4
      enddo
      close(20)
      
c
c write monthly plant respiration
      open(21,file='../outputs/mon_ar.bin',
!      open(21,file='../ipcc_env_vars2/'//
!     & 'mon_ar_HADCM3_A2_onlyclimate.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do k=1,12
         do i=1,nx
         do j=1,ny
            waux5(i,j) = mar(i,j,k)
         enddo
         enddo
         write(21,rec=k) waux5
      enddo
      close(21)
      
c
c write canopy resistance
      open(22,file='../outputs/rc.bin',
!      open(22,file='../ipcc_env_vars2/'//
!     & 'rc_HADCM3_A2_onlyclimate.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do k=1,12
         do i=1,nx
         do j=1,ny
            waux6(i,j) = rc(i,j,k)
         enddo
         enddo
         write(22,rec=k) waux6
      enddo
      close(22)
      
c write monthly cleaf 
      open(65,file='../outputs/mon_cleaf.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do k=1,12
         do i=1,nx
         do j=1,ny
            waux7(i,j) = cleaf(i,j,k)
         enddo
         enddo
         write(65,rec=k) waux7
      enddo
      close(65)
	  
c write monthly cawood 
      open(78,file='../outputs/mon_cawood.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do k=1,12
         do i=1,nx
         do j=1,ny
            waux20(i,j) = cawood(i,j,k)
         enddo
         enddo
         write(78,rec=k) waux20
      enddo
      close(78)

c write monthly cfroot 
      open(79,file='../outputs/mon_cfroot.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do k=1,12
         do i=1,nx
         do j=1,ny
            waux21(i,j) = cfroot(i,j,k)
         enddo
         enddo
         write(79,rec=k) waux21
      enddo
      close(79)	
      
c write trans_month2
      open (50,file='../outputs/trans_month.flt',
     &      status='unknown',form='unformatted',
     &      access='direct',recl=4*nx*ny)
      write(50,rec=1) trans_month
      close(50)
      
      
      open (60,file='../outputs/aux_des.flt',
     &      status='unknown',form='unformatted',
     &      access='direct',recl=4*nx*ny)
      write(60,rec=1) aux_des
      close(60)
c


      
c     !criando um arquivo ascii_lat/long/cell_id
c      open (70,file='../outputs/cell_id.txt',
c     &      status='unknown')
c      do i=1,6003
c      write(70,354)'cell_id:',cell_id(i),'lat:',aux_lat(i),'long:',
c     & aux_long(i)
c      enddo
c  354 format (a8,1x,i4,1x,a8,1x,f7.2,1x,a9,f7.2)
c      close (70)

C---------------------------------------------------------------------------------
c     !Criando um arquivo_dura‡Æo esta‡Æo seca
!      open (80,file='../outputs/lsmk_ni_nj.flt',
!     &      status='unknown',form='unformatted',
!     &      access='direct',recl=4*nx*ny)
!      write(80,rec=1) lsmk
!      close(80)
C---------------------------------------------------------------------------------

c
c
c Program end
c -----------
c
      stop
      end
c
c=======================================================================

      subroutine read12(nunit,var)
c auxiliar reading routine
      parameter(nx=192,ny=96)
      integer nunit
      real var(nx,ny,12)
      real aux(nx,ny)
      do k=1,12
        read(nunit,rec=k) aux
        do i=1,nx
          do j=1,ny
!-----------------------------------------------------------------------
!      if (i.le.96) then     !inverte hemisferios leste-oeste
!      ni=i+96
!      else
!      ni=(96-i)*(-1)
!      endif
!      nj=(j-96)*(-1)    !faz flip na latitude
!-----------------------------------------------------------------------
        var(i,j,k) = aux(i,j)

          enddo
        enddo
      enddo
      return
      end
	  
c========================================================================

      subroutine spinup (nppot_aux,
     &          cleafini_pft2,cawoodini_pft2,cfrootini_pft2)	 
      parameter(nx=192,ny=96,nt=5000,npls=3)
      real nppot_aux
      real cleafini_pft2(npls),cawoodini_pft2(npls),
     &	   cfrootini_pft2(npls)
      real cleafi_aux(nt),cawoodi_aux(nt),
     &   cfrooti_aux(nt),cbwoodi_aux(nt),
     &   cstoi_aux(nt),cotheri_aux(nt),crepi_aux(nt)

     
      real sensitivity,sensitivity2
      integer n
	     real aleaf(3) 		 !npp percentage alocated to leaf compartment
	     data aleaf /0.35,0.35,0.45/
         real aawood (3) !npp percentage alocated to aboveground woody biomass compartment
         data aawood /0.40,0.40,0.000000001/
   	     real afroot(3)  !npp percentage alocated to fine roots compartment
	     data afroot /0.25,0.25,0.55/ 
         real tleaf(3)   !turnover time of the leaf compartment (yr)
         data tleaf /1.0,0.5,1.0/ 
         real tawood (3)  !turnover time of the aboveground woody biomass compartment (yr)
	     data tawood /30.0,30.0,30.0/
  	     real tfroot(3)  !turnover time of the fine roots compartment
	     data tfroot /1.0,1.0,1.0/
			sensitivity = 1.10
		    sensitivity2 = 1.40
			
			
			   
			
			   
		    
	            
	          do i10=1,npls
             do k=1,nt
             
			 
             
			      
		     if (k.eq.1) then
			    
              cleafi_aux(k) = aleaf(i10)*(nppot_aux)
			 cawoodi_aux(k) =aawood(i10)*(nppot_aux)
			 cfrooti_aux(k) = afroot(i10)*(nppot_aux)
                 else
				 
             cleafi_aux(k) = ((aleaf(i10)*(nppot_aux))-
     &		(cleafi_aux(k-1)/(tleaf(i10)))) + cleafi_aux(k-1)
	        cawoodi_aux(k) = ((aawood(i10)*(nppot_aux))-
     &		(cawoodi_aux(k-1)/(tawood(i10)))) + cawoodi_aux(k-1)
	         cfrooti_aux(k) = ((afroot(i10)*(nppot_aux))-
     &		(cfrooti_aux(k-1)/(tfroot(i10)))) + cfrooti_aux(k-1)
	         
                 kk =  int(k*0.66)
				  if ((cawoodi_aux(k).ne.0.).and.(i10.ne.1)) then
!				  print*,i10, k,kk,cawoodi_aux(k),cawoodi_aux(kk)		 
	               endif
      if((cfrooti_aux(k)/cfrooti_aux(kk).lt.sensitivity).and.
     &	   (cleafi_aux(k)/cleafi_aux(kk).lt.sensitivity).and.
     &    (cawoodi_aux(k)/cawoodi_aux(kk).lt.sensitivity))then
	                
					cleafini_pft2(i10)=cleafi_aux(k)
					cawoodini_pft2(i10)=cawoodi_aux(k)
					cfrootini_pft2(i10)=cfrooti_aux(k)
					
			     
                 exit
                 endif
               endif
               enddo
               enddo
			   		   
                
			   
             return
               end				   
	 
	 
	 
	 
	 
	 
