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
      parameter(nx=192,ny=96,a=360,b=180)
      real lsmk(nx,ny),p0(nx,ny)
      real pr(nx,ny,12),t(nx,ny,12)
      real prec(nx,ny,12),temp(nx,ny,12),par(nx,ny,12),ipar(nx,ny,12),ca
      real anpr(nx,ny,12),ant(nx,ny,12)
      real tmin(nx,ny),seanpp(nx,ny),meanpp(nx,ny),meanhr(nx,ny),
     & meancs(nx,ny),mphoto(nx,ny),maresp(nx,ny),ave_wsoil(nx,ny),
     & ave_evap(nx,ny),ave_rc(nx,ny),nppot(nx,ny)
      
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
c carbon allocation
      real cleafini(nx,ny),cawoodini(nx,ny),cfrootini(nx,ny),
     & cbwoodini(nx,ny),cstoini(nx,ny),cotherini(nx,ny),crepini(nx,ny),
     & cleafini_tbe(nx,ny),cawoodini_tbe(nx,ny),cfrootini_tbe(nx,ny),
     & cbwoodini_tbe(nx,ny),cstoini_tbe(nx,ny),
     & cotherini_tbe(nx,ny),crepini_tbe(nx,ny),
     & cleafini_tbd(nx,ny),cawoodini_tbd(nx,ny),cfrootini_tbd(nx,ny),
     & cbwoodini_tbd(nx,ny),cstoini_tbd(nx,ny),
     & cotherini_tbd(nx,ny),crepini_tbd(nx,ny),
     & cleafini_herb(nx,ny),cawoodini_herb(nx,ny),cfrootini_herb(nx,ny),
     & cbwoodini_herb(nx,ny),cstoini_herb(nx,ny),
     & cotherini_herb(nx,ny),crepini_herb(nx,ny)
      real cleaf(nx,ny,12),cawood(nx,ny,12),cfroot(nx,ny,12),
     &     cbwood(nx,ny,12),csto(nx,ny,12),crep(nx,ny,12),
     &     cother(nx,ny,12),waux7(nx,ny),waux8(nx,ny),waux9(nx,ny),
     &     waux10(nx,ny),waux11(nx,ny),waux12(nx,ny),waux13(nx,ny),
     &     waux14(nx,ny),waux15(nx,ny),waux16(nx,ny),waux17(nx,ny),
     &     waux18(nx,ny),waux19(nx,ny),waux20(nx,ny),waux21(nx,ny),
     &     waux22(nx,ny),waux23(nx,ny),waux24(nx,ny),waux25(nx,ny),
     &     mean_cleaf(nx,ny),mean_cawood(nx,ny),mean_cbwood(nx,ny),
     &     mean_cfroot(nx,ny),mean_crep(nx,ny),mean_cother(nx,ny),
     &     mean_csto(nx,ny),
     &     mean_cleaf_tbe(nx,ny), mean_cawood_tbe(nx,ny),
     &     mean_cbwood_tbe(nx,ny),mean_cfroot_tbe(nx,ny),
     &     mean_csto_tbe(nx,ny),mean_crep_tbe(nx,ny),
     &     mean_cother_tbe(nx,ny),
     &     mean_cleaf_tbd(nx,ny), mean_cawood_tbd(nx,ny),
     &     mean_cbwood_tbd(nx,ny),mean_cfroot_tbd(nx,ny),
     &     mean_csto_tbd(nx,ny), mean_crep_tbd(nx,ny),
     &     mean_cother_tbd(nx,ny),
     &     mean_cleaf_herb(nx,ny), mean_cawood_herb(nx,ny),
     &     mean_cbwood_herb(nx,ny),mean_cfroot_herb(nx,ny),
     &     mean_csto_herb(nx,ny), mean_crep_herb(nx,ny),
     &     mean_cother_herb(nx,ny),
     &     ctotal_tbe(nx,ny),ctotal_tbd(nx,ny),ctotal_herb(nx,ny),
     &     ctotal_gd(nx,ny)
      integer i3 ! index for calculate the 3 PFTs (i4.eq.1 - TBE; i4.eq.2 - TBD; i4.eq.3 - HERB)
      parameter (npft=3)!it defines the amount of PFTs
c spinup
      integer supindex 


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
      open(20,file='../outputs/nppot.bin',status='old',
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
c
c --- Photosynthetically active radiation reaching canopy --------------
c     (ipar ; Ein/m2/s) [Eq. 7] observed data from ISLSCP2
      par(i,j,k) = ipar(i,j,k)/(2.18e5) !converting to Ein/m2/s
      temp(i,j,k) = t(i,j,k) !+ant(i,j,k) !uncomment to use future anomalies
      prec(i,j,k) = pr(i,j,k)!+anpr(i,j,k) !+pr(i,j,k)*0.2 !uncomment to use future anomalies
      if (prec(i,j,k).lt.0.0) prec (i,j,k) = 0.0

c
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
      
	  supindex = 1 ! it turns on or turns off the spinup to calculate initial carbon content on the compartments
	     if (supindex.eq.1) then
	     		 do i3=1,npft!index to calculate spinup for the three PFTs. i3=1 Tropical broadleaf evergreen; i3=2 Tropical broadleaf deciduous; i3=3 Tropical herbaceous
      
       call spinup (nppot,lsmk,i3,
     &          cleafini,cawoodini,cfrootini,cbwoodini,
     &          cstoini,cotherini,crepini,
     &          cleafini_tbe,cawoodini_tbe,
     &          cfrootini_tbe,cbwoodini_tbe,
     &          cstoini_tbe,cotherini_tbe,
     &          crepini_tbe,
     &          cleafini_tbd,cawoodini_tbd,
     &          cfrootini_tbd,cbwoodini_tbd,
     &          cstoini_tbd,cotherini_tbd,
     &          crepini_tbd,
     &          cleafini_herb,cawoodini_herb,
     &          cfrootini_herb,cbwoodini_herb,
     &          cstoini_herb,cotherini_herb,
     &          crepini_herb)	 


c	 write cleafini_tbe (inital carbon content in the leaf compartment for Tropical broadleaf evergreen PFT)
      open (99,file='../outputs/cleaf_initial_tbe.bin',
     &      status='unknown',form='unformatted',
     &      access='direct',recl=4*nx*ny)
      write(99,rec=1) cleafini_tbe
      close(99)
	  

	  
c write cawoodini_tbe (inital carbon content in the aboveground wood compartment for Tropical broadleaf evergreen PFT)
      open (102,file='../outputs/cawood_initial_tbe.bin',
     &      status='unknown',form='unformatted',
     &      access='direct',recl=4*nx*ny)
      write(102,rec=1) cawoodini_tbe
      close(102)
	  

c write cfrootini_tbe (inital carbon content in the fine roots compartment for Tropical broadleaf evergreen PFT)
      open (105,file='../outputs/cfroot_initial_tbe.bin',
     &      status='unknown',form='unformatted',
     &      access='direct',recl=4*nx*ny)
      write(105,rec=1) cfrootini_tbe
      close(105)
	  

	  
c write cbwoodini_tbe (inital carbon content in the belowground wood compartment for Tropical broadleaf evergreen PFT)
      open (108,file='../outputs/cbwood_initial_tbe.bin',
     &      status='unknown',form='unformatted',
     &      access='direct',recl=4*nx*ny)
      write(108,rec=1) cbwoodini_tbe
      close(108)
	  
	  
c write cstoini_tbe (inital carbon content in the storage compartment for Tropical broadleaf evergreen PFT)
      open (111,file='../outputs/csto_initial_tbe.bin',
     &      status='unknown',form='unformatted',
     &      access='direct',recl=4*nx*ny)
      write(111,rec=1) cstoini_tbe
      close(111)
	  
	  
c write cotherini_tbe (inital carbon content in the other compartment for Tropical broadleaf evergreen PFT)
      open (114,file='../outputs/cother_initial_tbe.bin',
     &      status='unknown',form='unformatted',
     &      access='direct',recl=4*nx*ny)
      write(114,rec=1) cotherini_tbe
      close(114)
	  
	  
c write crepini_tbe (inital carbon content in the reproduction compartment for Tropical broadleaf evergreen PFT)
      open (117,file='../outputs/crep_initial_tbe.bin',
     &      status='unknown',form='unformatted',
     &      access='direct',recl=4*nx*ny)
      write(117,rec=1) crepini_tbe
      close(117)
	  
c	 write cleafini_tbd (inital carbon content in the leaf compartment for Tropical broadleaf deciduous PFT)
      open (99,file='../outputs/cleaf_initial_tbd.bin',
     &      status='unknown',form='unformatted',
     &      access='direct',recl=4*nx*ny)
      write(99,rec=1) cleafini_tbd
      close(99)
	  

	  
c write cawoodini_tbd (inital carbon content in the aboveground wood compartment for Tropical broadleaf deciduous PFT)
      open (102,file='../outputs/cawood_initial_tbd.bin',
     &      status='unknown',form='unformatted',
     &      access='direct',recl=4*nx*ny)
      write(102,rec=1) cawoodini_tbd
      close(102)
	  

c write cfrootini_tbd (inital carbon content in the fine roots compartment for Tropical broadleaf deciduous PFT)
      open (105,file='../outputs/cfroot_initial_tbd.bin',
     &      status='unknown',form='unformatted',
     &      access='direct',recl=4*nx*ny)
      write(105,rec=1) cfrootini_tbd
      close(105)
	  

	  
c write cbwoodini_tbd (inital carbon content in the belowground wood compartment for Tropical broadleaf deciduous PFT)
      open (108,file='../outputs/cbwood_initial_tbd.bin',
     &      status='unknown',form='unformatted',
     &      access='direct',recl=4*nx*ny)
      write(108,rec=1) cbwoodini_tbd
      close(108)
	  
	  
c write cstoini_tbd (inital carbon content in the storage compartment for Tropical broadleaf deciduous PFT)
      open (111,file='../outputs/csto_initial_tbd.bin',
     &      status='unknown',form='unformatted',
     &      access='direct',recl=4*nx*ny)
      write(111,rec=1) cstoini_tbd
      close(111)
	  
	  
c write cotherini_tbd (inital carbon content in the other compartment for Tropical broadleaf deciduous PFT)
      open (114,file='../outputs/cother_initial_tbd.bin',
     &      status='unknown',form='unformatted',
     &      access='direct',recl=4*nx*ny)
      write(114,rec=1) cotherini_tbd
      close(114)
	  
	  
c write crepini_tbd (inital carbon content in the reproduction compartment for Tropical broadleaf deciduous PFT)
      open (117,file='../outputs/crep_initial_tbd.bin',
     &      status='unknown',form='unformatted',
     &      access='direct',recl=4*nx*ny)
      write(117,rec=1) crepini_tbd
      close(117)
	  
	  
c	 write cleafini_herb (inital carbon content in the leaf compartment for Tropical Herbaceous PFT)
      open (99,file='../outputs/cleaf_initial_herb.bin',
     &      status='unknown',form='unformatted',
     &      access='direct',recl=4*nx*ny)
      write(99,rec=1) cleafini_herb
      close(99)
	  

	  
c write cawoodini_herb (inital carbon content in the aboveground wood compartment for Tropical Herbaceous PFT)
      open (102,file='../outputs/cawood_initial_herb.bin',
     &      status='unknown',form='unformatted',
     &      access='direct',recl=4*nx*ny)
      write(102,rec=1) cawoodini_herb
      close(102)
	  

c write cfrootini_herb (inital carbon content in the fine roots compartment for Tropical Herbaceous PFT)
      open (105,file='../outputs/cfroot_initial_herb.bin',
     &      status='unknown',form='unformatted',
     &      access='direct',recl=4*nx*ny)
      write(105,rec=1) cfrootini_herb
      close(105)
	  

	  
c write cbwoodini_herb (inital carbon content in the belowground wood compartment for Tropical Herbaceous PFT)
      open (108,file='../outputs/cbwood_initial_herb.bin',
     &      status='unknown',form='unformatted',
     &      access='direct',recl=4*nx*ny)
      write(108,rec=1) cbwoodini_herb
      close(108)
	  
	  
c write cstoini_herb (inital carbon content in the storage compartment for Tropical Herbaceous PFT)
      open (111,file='../outputs/csto_initial_herb.bin',
     &      status='unknown',form='unformatted',
     &      access='direct',recl=4*nx*ny)
      write(111,rec=1) cstoini_herb
      close(111)
	  
	  
c write cotherini_herb (inital carbon content in the other compartment for Tropical Herbaceous PFT)
      open (114,file='../outputs/cother_initial_herb.bin',
     &      status='unknown',form='unformatted',
     &      access='direct',recl=4*nx*ny)
      write(114,rec=1) cotherini_herb
      close(114)
	  
	  
c write crepini_herb (inital carbon content in the reproduction compartment for Tropical Herbaceous PFT)
      open (117,file='../outputs/crep_initial_herb.bin',
     &      status='unknown',form='unformatted',
     &      access='direct',recl=4*nx*ny)
      write(117,rec=1) crepini_herb
      close(117)	  
	  

	      enddo
	  
	  
	  
	  
	     else if (supindex.eq.2) then 
     
         open (121,file='../outputs/cleaf_initial_tbe.bin',status='old',
     &         form='unformatted',access='direct',recl=4*nx*ny)
	        read (121,rec=1) cleafini_tbe
		    close(121)

         open(123,file='../outputs/cfroot_initial_tbe.bin',status='old',
     &        form='unformatted',access='direct',recl=4*nx*ny)
		    read (123,rec=1) cfrootini_tbe
	      close (123) 
		  
       
        open(125,file='../outputs/cawood_initial_tbe.bin',status='old',
     &        form='unformatted',access='direct',recl=4*nx*ny)
		    read (125,rec=1) cawoodini_tbe
	      close (125)
		  
		  
        open(128,file='../outputs/cbwood_initial_tbe.bin',status='old',
     &        form='unformatted',access='direct',recl=4*nx*ny)
		    read (128,rec=1) cbwoodini_tbe
	      close (128)
	
	     open(131,file='../outputs/csto_initial_tbe.bin',status='old',
     &        form='unformatted',access='direct',recl=4*nx*ny)
		    read (131,rec=1) cstoini_tbe
	      close (131)
		  
		  
        open(134,file='../outputs/cother_initial_tbe.bin',status='old',
     &        form='unformatted',access='direct',recl=4*nx*ny)
		    read (134,rec=1) cotherini_tbe
	      close (134)
		  
        
		    open(137,file='../outputs/crep_initial_tbe.bin',status='old',
     &        form='unformatted',access='direct',recl=4*nx*ny)
		    read (137,rec=1) crepini_tbe
	      close (137)
		  
       open (121,file='../outputs/cleaf_initial_tbd.bin',status='old',
     &         form='unformatted',access='direct',recl=4*nx*ny)
	        read (121,rec=1) cleafini_tbd
		    close(121)

         open(123,file='../outputs/cfroot_initial_tbd.bin',status='old',
     &        form='unformatted',access='direct',recl=4*nx*ny)
		    read (123,rec=1) cfrootini_tbd
	      close (123) 
		  
       
        open(125,file='../outputs/cawood_initial_tbd.bin',status='old',
     &        form='unformatted',access='direct',recl=4*nx*ny)
		    read (125,rec=1) cawoodini_tbd
	      close (125)
		  
		  
        open(128,file='../outputs/cbwood_initial_tbd.bin',status='old',
     &        form='unformatted',access='direct',recl=4*nx*ny)
		    read (128,rec=1) cbwoodini_tbd
	      close (128)
	
	     open(131,file='../outputs/csto_initial_tbd.bin',status='old',
     &        form='unformatted',access='direct',recl=4*nx*ny)
		    read (131,rec=1) cstoini_tbd
	      close (131)
		  
		  
        open(134,file='../outputs/cother_initial_tbd.bin',status='old',
     &        form='unformatted',access='direct',recl=4*nx*ny)
		    read (134,rec=1) cotherini_tbd
	      close (134)
		  
        
		    open(137,file='../outputs/crep_initial_tbd.bin',status='old',
     &        form='unformatted',access='direct',recl=4*nx*ny)
		    read (137,rec=1) crepini_tbd
	      close (137)

      open (121,file='../outputs/cleaf_initial_herb.bin',status='old',
     &         form='unformatted',access='direct',recl=4*nx*ny)
	        read (121,rec=1) cleafini_herb
		    close(121)

      open(123,file='../outputs/cfroot_initial_herb.bin',status='old',
     &        form='unformatted',access='direct',recl=4*nx*ny)
		    read (123,rec=1) cfrootini_herb
	      close (123) 
		  
       
      open(125,file='../outputs/cawood_initial_herb.bin',status='old',
     &        form='unformatted',access='direct',recl=4*nx*ny)
		    read (125,rec=1) cawoodini_herb
	      close (125)
		  
		  
      open(128,file='../outputs/cbwood_initial_herb.bin',status='old',
     &        form='unformatted',access='direct',recl=4*nx*ny)
		    read (128,rec=1) cbwoodini_herb
	      close (128)
	
      open(131,file='../outputs/csto_initial_herb.bin',status='old',
     &        form='unformatted',access='direct',recl=4*nx*ny)
		    read (131,rec=1) cstoini_herb
	      close (131)
		  
		  
      open(134,file='../outputs/cother_initial_herb.bin',status='old',
     &        form='unformatted',access='direct',recl=4*nx*ny)
		    read (134,rec=1) cotherini_herb
	      close (134)
		  
        
      open(137,file='../outputs/crep_initial_herb.bin',status='old',
     &        form='unformatted',access='direct',recl=4*nx*ny)
		    read (137,rec=1) crepini_herb
	      close (137)		  
		 
		  
		    endif
   
   
   
   
c -------------------------------------------------------
c Calculate environmental variables (water balance model)
c -------------------------------------------------------
c
      call wbm (prec,temp,lsmk,p0,ca,par,
     &          cleafini,cfrootini,cawoodini,
     &          cbwoodini,cstoini,cotherini,crepini,
     &          cleafini_tbe,cfrootini_tbe,cawoodini_tbe,
     &			cbwoodini_tbe,cstoini_tbe,
     &          cotherini_tbe,crepini_tbe,
     &          cleafini_tbd,cfrootini_tbd,cawoodini_tbd,
     &		    cbwoodini_tbd,cstoini_tbd,
     &          cotherini_tbd,crepini_tbd,
     &          cleafini_herb,cfrootini_herb,cawoodini_herb,
     &		    cbwoodini_herb,cstoini_herb,
     &          cotherini_herb,crepini_herb,
     &          tmin,meanpp,seanpp,mphoto,maresp,  !carbon
     &          meanrml,meanrmf,meanrms,meanrm,
     &          meanrgl,meanrgf,meanrgs,meanrg,
     &          meanhr,meancs,wsoil2,evaptr,mnpp,
     &          mp,mar,rc,ave_wsoil,ave_evap,ave_rc,
     &          monrml,monrmf,monrms,monrm,
     &          monrgl,monrgf,monrgs,monrg,
     &          cleaf,cawood,cfroot,cbwood,csto,crep,cother,
     &          mean_cleaf,mean_cawood,mean_cbwood,mean_cfroot,
     &          mean_csto,mean_crep,mean_cother,
     &          mean_cleaf_tbe,mean_cawood_tbe,
     &          mean_cbwood_tbe,mean_cfroot_tbe,
     &          mean_csto_tbe,mean_crep_tbe,
     &          mean_cother_tbe,
     &          mean_cleaf_tbd,mean_cawood_tbd,
     &			mean_cbwood_tbd,mean_cfroot_tbd,
     &          mean_csto_tbd,mean_crep_tbd,
     &			mean_cother_tbd,
     &          mean_cleaf_herb,mean_cawood_herb,
     &			mean_cbwood_herb,mean_cfroot_herb,
     &          mean_csto_herb,mean_crep_herb,
     &			mean_cother_herb,
     &          ctotal_tbe,ctotal_tbd,ctotal_herb,ctotal_gd)

     c     !exerc¡cio: nomear as celulas de grad ----------------------------

       do i=1,nx
         do j=1,ny

      if (lsmk(i,j).eq.1) then
         var_id = var_id + 1.0
         cell_id(var_id)= var_id
         aux_lat(var_id)=lat(i,j)
         aux_long(var_id)=long(i,j)
		 
	      if (cell_id(var_id) .eq. 5244)  then

      print*, cell_id(var_id),'cleaf',cleaf(i,j,7),
     &       'cawood',cawood(i,j,7),'cfroot',cfroot(i,j,7),
     &       'cbwood',cbwood(i,j,7),'csto',csto(i,j,7),
     &       'crep',crep(i,j,7),'cother',cother(i,j,7)
	 
c  100 format (A8,1x,I4,1x,A8,1x,F7.2,1x,A9,f7.2)

      endif
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

c write monthly cbwood 
      open(80,file='../outputs/mon_cbwood.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do k=1,12
         do i=1,nx
         do j=1,ny
            waux22(i,j) = cbwood(i,j,k)
         enddo
         enddo
         write(80,rec=k) waux22
      enddo
      close(80)	

c write monthly csto 
      open(81,file='../outputs/mon_csto.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do k=1,12
         do i=1,nx
         do j=1,ny
            waux23(i,j) = csto(i,j,k)
         enddo
         enddo
         write(81,rec=k) waux23
      enddo
      close(81)	

c write monthly crep 
      open(82,file='../outputs/mon_crep.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do k=1,12
         do i=1,nx
         do j=1,ny
            waux24(i,j) = crep(i,j,k)
         enddo
         enddo
         write(82,rec=k) waux24
      enddo
      close(82)

c write monthly cother 
      open(83,file='../outputs/mon_cother.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do k=1,12
         do i=1,nx
         do j=1,ny
            waux25(i,j) = cother(i,j,k)
         enddo
         enddo
         write(83,rec=k) waux25
      enddo
      close(83) 	  
	  
c write january cleaf 
      open(66,file='../outputs/cleaf1.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do k=1,12
         do i=1,nx
         do j=1,ny
            waux8(i,j) = cleaf(i,j,1)
         enddo
         enddo
         write(66,rec=k) waux8
      enddo
      close(66)	  
	  
c write february cleaf 
      open(67,file='../outputs/cleaf2.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do k=1,12
         do i=1,nx
         do j=1,ny
            waux9(i,j) = cleaf(i,j,2)
         enddo
         enddo
         write(67,rec=k) waux9
      enddo
      close(67)	 

c write march cleaf 
      open(68,file='../outputs/cleaf3.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do k=1,12
         do i=1,nx
         do j=1,ny
            waux10(i,j) = cleaf(i,j,3)
         enddo
         enddo
         write(68,rec=k) waux10
      enddo
      close(68)	

c write april cleaf 
      open(69,file='../outputs/cleaf4.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do k=1,12
         do i=1,nx
         do j=1,ny
            waux11(i,j) = cleaf(i,j,4)
         enddo
         enddo
         write(69,rec=k) waux11
      enddo
      close(69)	 

c write may cleaf 
      open(70,file='../outputs/cleaf5.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do k=1,12
         do i=1,nx
         do j=1,ny
            waux12(i,j) = cleaf(i,j,5)
         enddo
         enddo
         write(70,rec=k) waux12
      enddo
      close(70)

c write june cleaf 
      open(71,file='../outputs/cleaf6.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do k=1,12
         do i=1,nx
         do j=1,ny
            waux13(i,j) = cleaf(i,j,6)
         enddo
         enddo
         write(71,rec=k) waux13
      enddo
      close(71)

c write july cleaf 
      open(72,file='../outputs/cleaf7.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do k=1,12
         do i=1,nx
         do j=1,ny
            waux14(i,j) = cleaf(i,j,7)
         enddo
         enddo
         write(72,rec=k) waux14
      enddo
      close(72)

c write august cleaf 
      open(73,file='../outputs/cleaf8.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do k=1,12
         do i=1,nx
         do j=1,ny
            waux15(i,j) = cleaf(i,j,8)
         enddo
         enddo
         write(73,rec=k) waux15
      enddo
      close(73)	

c write september cleaf 
      open(74,file='../outputs/cleaf9.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do k=1,12
         do i=1,nx
         do j=1,ny
            waux16(i,j) = cleaf(i,j,9)
         enddo
         enddo
         write(74,rec=k) waux16
      enddo
      close(74)	

c write october cleaf 
      open(75,file='../outputs/cleaf10.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do k=1,12
         do i=1,nx
         do j=1,ny
            waux17(i,j) = cleaf(i,j,10)
         enddo
         enddo
         write(75,rec=k) waux17
      enddo
      close(75)

c write november cleaf 
      open(76,file='../outputs/cleaf11.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do k=1,12
         do i=1,nx
         do j=1,ny
            waux18(i,j) = cleaf(i,j,11)
         enddo
         enddo
         write(76,rec=k) waux18
      enddo
      close(76)

c write december cleaf 
      open(77,file='../outputs/cleaf12.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do k=1,12
         do i=1,nx
         do j=1,ny
            waux19(i,j) = cleaf(i,j,12)
         enddo
         enddo
         write(77,rec=k) waux19
      enddo
      close(77)	 	  


      
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

c write cleafini (inital carbon content in the leaf compartment)
      open (80,file='../outputs/cleaf_initial.bin',
     &      status='unknown',form='unformatted',
     &      access='direct',recl=4*nx*ny)
      write(80,rec=1) cleafini
      close(80)

c write cawoodini (inital carbon content in the aboveground woody biomass compartment)
      open (88,file='../outputs/cawood_initial.bin',
     &      status='unknown',form='unformatted',
     &      access='direct',recl=4*nx*ny)
      write(88,rec=1) cawoodini
      close(88)
	  
c write cfrootini (inital carbon content in the aboveground woody biomass compartment)
      open (90,file='../outputs/cfroot_initial.bin',
     &      status='unknown',form='unformatted',
     &      access='direct',recl=4*nx*ny)
      write(90,rec=1) cfrootini
      close(90)	 
	  
c write cbwoodini (inital carbon content in the belowground woody biomass compartment)
      open (92,file='../outputs/cbwood_initial.bin',
     &      status='unknown',form='unformatted',
     &      access='direct',recl=4*nx*ny)
      write(92,rec=1) cbwoodini
      close(92)

c write cstoini (inital carbon content in the storage compartment)
      open (94,file='../outputs/csto_initial.bin',
     &      status='unknown',form='unformatted',
     &      access='direct',recl=4*nx*ny)
      write(94,rec=1) cstoini
      close(94)
	  
c write cotherini (inital carbon content in the other compartment)
      open (96,file='../outputs/cother_initial.bin',
     &      status='unknown',form='unformatted',
     &      access='direct',recl=4*nx*ny)
      write(96,rec=1) cotherini
      close(96)	
	  
c write crepini (inital carbon content in the reproduction compartment)
      open (98,file='../outputs/crep_initial.bin',
     &      status='unknown',form='unformatted',
     &      access='direct',recl=4*nx*ny)
      write(98,rec=1) crepini
      close(98)	   	  

      
      
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

      subroutine spinup (nppot,lsmk,i3,
     &          cleafini,cawoodini,cfrootini,cbwoodini,
     &          cstoini,cotherini,crepini,
     &          cleafini_tbe,cawoodini_tbe,
     &          cfrootini_tbe,cbwoodini_tbe,
     &          cstoini_tbe,cotherini_tbe,
     &          crepini_tbe,
     &          cleafini_tbd,cawoodini_tbd,
     &          cfrootini_tbd,cbwoodini_tbd,
     &          cstoini_tbd,cotherini_tbd,
     &          crepini_tbd,
     &          cleafini_herb,cawoodini_herb,
     &          cfrootini_herb,cbwoodini_herb,
     &          cstoini_herb,cotherini_herb,
     &          crepini_herb)	 
      parameter(nx=192,ny=96,nt=2500)
      real nppot(nx,ny),cleafini(nx,ny),lsmk(nx,ny),cawoodini(nx,ny),
     &     cfrootini(nx,ny),cbwoodini(nx,ny),cstoini(nx,ny),
     &     cotherini(nx,ny),crepini(nx,ny)
      real cleafini_tbe(nx,ny),cawoodini_tbe(nx,ny),
     &     cfrootini_tbe(nx,ny),cbwoodini_tbe(nx,ny),
     &     cstoini_tbe(nx,ny),cotherini_tbe(nx,ny),
     &     crepini_tbe(nx,ny),
     &     cleafini_tbd(nx,ny),cawoodini_tbd(nx,ny),
     &     cfrootini_tbd(nx,ny),cbwoodini_tbd(nx,ny),
     &     cstoini_tbd(nx,ny),cotherini_tbd(nx,ny),
     &     crepini_tbd(nx,ny),
     &     cleafini_herb(nx,ny),cawoodini_herb(nx,ny),
     &     cfrootini_herb(nx,ny),cbwoodini_herb(nx,ny),
     &     cstoini_herb(nx,ny),cotherini_herb(nx,ny),
     &     crepini_herb(nx,ny)
      real cleafi_aux(nx,ny,nt),cawoodi_aux(nx,ny,nt),
     &   cfrooti_aux(nx,ny,nt),cbwoodi_aux(nx,ny,nt),
     &   cstoi_aux(nx,ny,nt),cotheri_aux(nx,ny,nt),crepi_aux(nx,ny,nt)
     
      real sensitivity,sensitivity2
      integer n
	     real aleaf(3) 		 !npp percentage alocated to leaf compartment
	     data aleaf /0.23,0.23,0.32/
         real aawood (3) !npp percentage alocated to aboveground woody biomass compartment
         data aawood /0.22,0.22,0.001/
   	     real afroot(3)  !npp percentage alocated to fine roots compartment
	     data afroot /0.15,0.15,0.42/ 
         real abwood(3)  !npp percentage alocated to belowground woody biomass compartment
	     data abwood /0.10,0.10,0.001/
!         real asto(3)    !npp percentage alocated to storage compartment
!         data asto /0.10,0.10,0.10/
!         real arep(3)    !npp percentage alocated to reproduction compartment
!	     data arep /0.15,0.15,0.10/
!         real aother(3)  !npp percentage alocated to other compartment
!	     data aother /0.05,0.05,0.06/ 
         real tleaf(3)   !turnover time of the leaf compartment (yr)
         data tleaf /1.0,0.5,1.0/ 
         real tawood (3)  !turnover time of the aboveground woody biomass compartment (yr)
	     data tawood /20.0,20.0,20.0/
  	     real tfroot(3)  !turnover time of the fine roots compartment
	     data tfroot /1.0,1.0,1.0/
         real tbwood (3)		 !turnover time of the belowground woody biomass compartment
	     data tbwood /30.0,30.0,30.0/
!         real tsto  (3)    !turnover time of the storage compartmentturn
!	     data tsto /5.0,5.0,5.0/ 
!         real trep (3)   !turnover time of the reproduction compartment
!	     data trep /0.25,0.25,0.25/ 
!         real tother (3) !turnover time of the other compartment
!	     data tother /0.12,0.12,0.12/
			sensitivity = 1.10
		    sensitivity2 = 1.40
		    do i=1,nx
             do j=1,ny
             if (lsmk(i,j).eq.1) then
            do k=1,nt
             n = n+1
		     if (k.eq.1) then
              cleafi_aux(i,j,k) = aleaf(i3)*(nppot(i,j))
			 cawoodi_aux(i,j,k) = aawood(i3)*(nppot(i,j))
			 cfrooti_aux(i,j,k) = afroot(i3)*(nppot(i,j))
			 cbwoodi_aux(i,j,k) = abwood(i3)*(nppot(i,j))
!			 cstoi_aux(i,j,k) = asto(i3)*(nppot(i,j))
!			 cotheri_aux(i,j,k) = aother(i3)*(nppot(i,j)/365)
!			 crepi_aux(i,j,k) = arep(i3)*(nppot(i,j)/365)
                 else
             cleafi_aux(i,j,k) = ((aleaf(i3)*(nppot(i,j)))-
     &		(cleafi_aux(i,j,k-1)/(tleaf(i3)))) + cleafi_aux(i,j,k-1)
	        cawoodi_aux(i,j,k) = ((aawood(i3)*(nppot(i,j)))-
     &		(cawoodi_aux(i,j,k-1)/(tawood(i3)))) + cawoodi_aux(i,j,k-1)
	         cfrooti_aux(i,j,k) = ((afroot(i3)*(nppot(i,j)))-
     &		(cfrooti_aux(i,j,k-1)/(tfroot(i3)))) + cfrooti_aux(i,j,k-1)
	         cbwoodi_aux(i,j,k) = ((abwood(i3)*(nppot(i,j)))-
     &		(cbwoodi_aux(i,j,k-1)/(tbwood(i3)))) + cbwoodi_aux(i,j,k-1)
!	        cstoi_aux(i,j,k) = ((asto(i3)*(nppot(i,j)))-
!     &		(cstoi_aux(i,j,k-1)/(tsto(i3)))) + cstoi_aux(i,j,k-1)
!	         cotheri_aux(i,j,k) = ((aother(i3)*(nppot(i,j)))-
!     &		(cotheri_aux(i,j,k-1)/(tother(i3)*365))) + cotheri_aux(i,j,k-1)
!	        crepi_aux(i,j,k) = ((arep(i3)*(nppot(i,j)))-
!     &		(crepi_aux(i,j,k-1)/(trep(i3)*365))) + crepi_aux(i,j,k-1)
                 kk =  int(k*0.66)
!       if((cfrooti_aux(i,j,k)/cfrooti_aux(i,j,kk).lt.sensitivity).and.
!     &	   (cleafi_aux(i,j,k)/cleafi_aux(i,j,kk).lt.sensitivity).and.
!     &     (cawoodi_aux(i,j,k)/cawoodi_aux(i,j,kk).lt.sensitivity2).and.
!     &     (cbwoodi_aux(i,j,k)/cbwoodi_aux(i,j,kk).lt.sensitivity2)then
!     &	   (cstoi_aux(i,j,k)/cstoi_aux(i,j,kk).lt.sensitivity).and.
!     &     (cotheri_aux(i,j,k)/cotheri_aux(i,j,kk).lt.sensitivity).and.
!     &     (crepi_aux(i,j,k)/crepi_aux(i,j,kk).lt.sensitivity))   then
	                
       if((cfrooti_aux(i,j,k)/cfrooti_aux(i,j,kk).lt.sensitivity).and.
     &	 (cleafi_aux(i,j,k)/cleafi_aux(i,j,kk).lt.sensitivity).and.
     &   (cawoodi_aux(i,j,k)/cawoodi_aux(i,j,kk).lt.sensitivity2).and.
     &   (cbwoodi_aux(i,j,k)/cbwoodi_aux(i,j,kk).lt.sensitivity2))then
					cleafini(i,j) = cleafi_aux(i,j,k)
					cawoodini(i,j) = cawoodi_aux(i,j,k)
					cfrootini(i,j) = cfrooti_aux(i,j,k)
					cbwoodini(i,j) = cbwoodi_aux(i,j,k)
!					cstoini(i,j) = cstoi_aux(i,j,k)
!					cotherini(i,j) = cotheri_aux(i,j,k)
!					crepini(i,j) = crepi_aux(i,j,k)
					
					 if (i3.eq.1) then
					  cleafini_tbe(i,j)= cleafini(i,j)
					  cawoodini_tbe(i,j)= cawoodini(i,j) 
					  cfrootini_tbe(i,j)= cfrootini(i,j) 
					  cbwoodini_tbe(i,j)= cbwoodini(i,j) 
!					  cstoini_tbe(i,j)= cstoini(i,j)
!					  cotherini_tbe(i,j)= cotherini(i,j)
!					  crepini_tbe(i,j)= crepini(i,j)
					  else if (i3.eq.2) then
					  cleafini_tbd(i,j)= cleafini(i,j)
					  cawoodini_tbd(i,j)= cawoodini(i,j) 
					  cfrootini_tbd(i,j)= cfrootini(i,j) 
					  cbwoodini_tbd(i,j)= cbwoodini(i,j) 
!					  cstoini_tbd(i,j)= cstoini(i,j)
!					  cotherini_tbd(i,j)= cotherini(i,j)
!					  crepini_tbd(i,j)= crepini(i,j)
					  else if (i3.eq.3) then
					  cleafini_herb(i,j)= cleafini(i,j)
					  cawoodini_herb(i,j)= 0.
					  cfrootini_herb(i,j)= cfrootini(i,j) 
					  cbwoodini_herb(i,j)= 0.
!					  cstoini_herb(i,j)= cstoini(i,j)
!					  cotherini_herb(i,j)= cotherini(i,j)
!					  crepini_herb(i,j)= crepini(i,j)
                     endif					  
					
                 exit
                 endif
               endif
               enddo
               endif
               enddo
               enddo	

             return
               end				   
	 
	 
	 
	 
	 
	 