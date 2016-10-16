c23456789
!
      program env
!
      implicit none

!=======================================================================
! CPTEC-PVM2
! Adapted from env4.f.
! Code written by David Lapola
! Last update: Aug/2007
!
! Compile with: g95 (or gfortran) env5.f wbm4.f productivity1.f -o a.exe
! Execute with ./a.exe
! Then run g95 (gfortran) pvm5.f -o pvm5.exe
! Execute with ./pvm5.exe
!
!=======================================================================
!
! Parameters and variables
! ------------------------
!
      integer nx,ny,a,b,i,j,k
      integer cell_id(6003),var_id,ni
      parameter(nx=192,ny=96,a=360,b=180)
      real, parameter :: NO_DATA = -9999.0
      real cellx,celly
      real lsmk(nx,ny),p0(nx,ny)
      real pr(nx,ny,12),t(nx,ny,12)
      real prec(nx,ny,12),temp(nx,ny,12),par(nx,ny,12),ipar(nx,ny,12),ca
      real anpr(nx,ny,12),ant(nx,ny,12)
      real tmin(nx,ny),seanpp(nx,ny),meanpp(nx,ny),meanhr(nx,ny),
     &     meancs(nx,ny),mphoto(nx,ny),maresp(nx,ny),ave_wsoil(nx,ny),
     &     ave_evap(nx,ny),ave_cres(nx,ny)
      real bpot(nx,ny)
      real wsoil2(nx,ny,12),evaptr(nx,ny,12),cres(nx,ny,12)
      real waux(nx,ny),waux2(nx,ny),waux3(nx,ny),waux4(nx,ny),
     &     waux5(nx,ny),waux6(nx,ny),waux7(nx,ny),waux8(nx,ny),
     &     waux9(nx,ny),waux10(nx,ny),waux11(nx,ny),waux12(nx,ny),
     &     waux13(nx,ny),waux14(nx,ny)
      real long(nx,ny),lat(nx,ny),long_graus,lat_graus,aux_long(6003),
     &     aux_lat(6003),trans_month(nx,ny),trans_month2(nx,ny)
      real aux_des(nx,ny),aux_des2(nx,ny)
      real cres1(nx,ny,12),cres2(nx,ny,12),cres3(nx,ny,12)
      real ave_cres1(nx,ny),ave_cres2(nx,ny),ave_cres3(nx,ny)
      real npp1(nx,ny,12),npp2(nx,ny,12),npp3(nx,ny,12)
      real meanpp1(nx,ny),meanpp2(nx,ny),meanpp3(nx,ny)
      real mnpp1_pon(nx,ny),mnpp2_pon(nx,ny),mnpp3_pon(nx,ny)
      real wsoil(nx,ny,12),runom(nx,ny,12)
!
! Carbon Cycle

      real mp(nx,ny,12),mar(nx,ny,12),mnpp(nx,ny,12),
     &     lai(nx,ny,12),clit(nx,ny,12),csoil(nx,ny,12),
     &     hresp(nx,ny,12)
!
! Open files
! ----------
!
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
!
! =============================================================================================
!
! Future anomalies see lines 85 and 86
! ------------------------------------
!
!     open(14,file='../ipcc_ar4/anomalias/SRESA2/temp/'//
!    &               'HADCM3_A2_2070_2099_anom_t.bin',status='old',
!    &        form='unformatted',access='direct',recl=4*nx*ny)
!     open(15,file='../ipcc_ar4/anomalias/SRESA2/prec/'//
!    &               'HADCM3_A2_2070_2099_anom_pr.bin',status='old',
!    &        form='unformatted',access='direct',recl=4*nx*ny)
!
! =============================================================================================
!
! Read data
! ---------
!
      read ( 9,rec=1) p0                                                                       !mb
      read (10,rec=1) lsmk                                                                     !land = 1; ocean = 0
      call read12 (11,pr)                                                                      !mm/month
      call read12 (12,t)                                                                       !oC
      call read12 (13,ipar)                                                                    !w/m2
!      call read12 (14,ant)                                                                    !oC
!      call read12 (15,anpr)                                                                   !oC
!
! Close files
! -----------
!
      close ( 9)
      close (10)
      close (11)
      close (12)
      close (13)
!      close (14)
!      close (15)
!
!==============================================================================================
!
!     do i=1,nx
!     do j=1,ny
!
!     if (i.le.96) then                                                                        !inverte hemisferios leste-oeste
!     ni=i+96
!     else
!     ni=(96-i)*(-1)
!     endif
!
!     nj=(j-96)*(-1)                                                                           !Flip na latitude
!
!     pO(i,j)=p0(ni,nj)
!     lsmk(i,j)=lsmk(ni,nj)
!     aux_des2(i,j)=aux_des(ni,nj)                                                             !Aux_des com nova configuracao de mapa
!     trans_month2(i,j)=trans_month(ni,nj)                                                     !Trans_month com nova configuracao de mapa
!
!     enddo
!     enddo
!
!==============================================================================================
!
      var_id= 0                                                                                !cell_id
!     prec_t= 0.0                                                                              !total precipitation
!
!     do i=1,nx
!     do j=1,ny
!
!     aux_des(i,j)= 0                                                                          !months with less then 100mm of precipitation
!     p0(i,j)= 1010.0
!
!     enddo
!     enddo
!
!==============================================================================================
!
      do k=1,12
         do i=1,nx
         do j=1,ny

! Photosynthetically active radiation reaching canopy
! ipar; Ein/m2/s (observed data from ISLSCP2)

         par(i,j,k) = (ipar(i,j,k)/(2.18e5))                                                      !converting to Ein/m2/s
         temp(i,j,k) = t(i,j,k)                                                                   !+ant(i,j,k) !uncomment to use future anomalies
         prec(i,j,k) = pr(i,j,k)                                                                  !+anpr(i,j,k) !+pr(i,j,k)*0.2 !uncomment to use future anomalies
         if (prec(i,j,k).lt.0.0) prec (i,j,k) = 0.0
!
! Atmospheric CO2 pressure                                                                     !ppmv/Pa
! 1 Pa CO2 = 9.901 ppmv CO2 (Adams et al. 2004)
!     ca= 18.18                                                                                !Pa (=180 ppmv; Last Glacial Maximum)
!     ca= 28.28                                                                                !Pa (=280 ppmv; Pre-Industrial Revolution)
!     ca= 35.35                                                                                !Pa (=350 ppmv; 1961-1990)
      ca= 350/9.901
!      ca= 363.20/9.901                                                                         !Pa (=363.20 ppmv; average_1980/2010)
!     ca= 54.03                                                                                !Pa (=535 ppmv; SRES-B1 2080's)
!     ca= 73.73                                                                                !Pa (=730 ppmv; SRES-A2 2080's)
!     ca= ((73.73-35.35)*0.5)+35.35                                                            !Pa half effect!!!
!
!==============================================================================================
!
! Duracao:estacao seca
! --------------------
!
!     if(prec(i,j,k).lt.100.0) then
!     aux_des(i,j)=aux_des(i,j)+1.0
!     else
!     aux_des(i,j)=aux_des(i,j)+0.0
!     endif
!
! Transition month (dry_rainy season)
! -----------------------------------
!
!     if (prec(i,j,k-1).lt.100.and.prec(i,j,k-2).lt.100
!    & .and.prec(i,j,k-3).lt.100.and.prec(i,j,k).gt.100) then
!     trans_month(i,j)=k-1
!     endif
!
!     if (trans_month(i,j).eq.0) then
!     trans_month(i,j)=12
!     endif
!
!     prec_t=prec_t+prec(i,j,k)
!
      enddo
      enddo
      enddo
!
!==============================================================================================
!
! Calculate environmental variables (WBM)
! ----  -----------------------------------
!
      call wbm (prec,temp,lsmk,p0,ca,par,
     &          tmin,meanpp,seanpp,mphoto,maresp,
     &          meanhr,meancs,wsoil2,evaptr,mnpp,mp,mar,cres,
     &          ave_wsoil,ave_evap,ave_cres,
     &          cres1,cres2,cres3,ave_cres1,ave_cres2,ave_cres3,
     &          npp1,npp2,npp3,meanpp1,meanpp2,meanpp3,wsoil,runom,
     &          mnpp1_pon,mnpp2_pon,mnpp3_pon)
!
!==============================================================================================
!
! Nomeando celulas de grad
! ------------------------
c!
c      do i=1,nx
c         do j=1,ny
!c
c            if (nint(lsmk(i,j)).eq.1) then
c               var_id = var_id + 1
c               cell_id(var_id) = var_id
c               aux_lat(var_id) = lat(i,j)
c               aux_long(var_id) = long(i,j)
c               do k=1,12
c                  if (cell_id(var_id).eq.5245) then
c!     print*,rc(i,j,k)
c                  endif
c               enddo
c               if (cell_id(var_id).eq.5245) then
!     print*,ave_rc(i,j)
c               endif
!
c               cellx= 134       !longitude
c               celly= -25       !latitude
c               long_graus = a/nx !extensao longitudinal de cada celula em graus
c               lat_graus = b/ny !extensao latitudinal de cada celula em graus
c               long(i,j) = long_graus*(i) !coloca em coordenada "0-360"
c               if (long(i,j).lt.180) then
c                  long(i,j) = long(i,j) !transforma para coordenada "-180|180"
c               else
c                  long(i,j) = (360-long(i,j))*(-1)
c              endif
c               lat(i,j) = (lat_graus*(j)) !coloca em coordenada "0-180'
c               lat(i,j) = (lat(i,j) - 90.0)
!
c      if (lat(i,j).ge.celly.and.lat(i,j).lt.(celly+1.875)            !para o ano inteiro
c     &        .and.long(i,j).ge.cellx.and.long(i,j)
c     &        .lt.(cellx+1.875)) then
!     print*,cell_id(var_id),"oi"
c      endif
!     
c      if (cell_id(var_id).eq.3197) then
!      print*,ave_rc(i,j)
c      endif
!
c      if (cell_id(var_id).eq.3197) then
!      print*,mnpp(i,j,k),meanpp(i,j),lat(i,j),long(i,j)
c      endif
!
c      if (cell_id(var_id).eq.5245) then
!      print*,ave_rc(i,j),lat(i,j),long(i,j)
c      endif
!
c      endif
c      enddo
c      enddo
!
!==============================================================================================
!
! Write environmental variables
! -----------------------------
!
      open(16,file='../outputs/ambientais4.bin',                                               !open(16,file='../ipcc_env_vars2/'//& 'ambientais4_HADCM3_A2_onlyclimate.bin',
!      open(16,file='../ipcc_env_vars2/'//
!     & 'ambientais4_HADCM3_A2_onlyclimate.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      write(16,rec=1) tmin                                                                     !oC
      write(16,rec=2) meanpp                                                                   !kgC/m2/yr
      write(16,rec=3) seanpp                                                                   !dimensionless
      write(16,rec=4) meanhr                                                                   !kgC/m2/y
      write(16,rec=5) meancs                                                                   !kgC/m2/y
      write(16,rec=6) mphoto                                                                   !kgC/m2/y
      write(16,rec=7) maresp                                                                   !kgC/m2/y
      write(16,rec=8) ave_wsoil                                                                !mm
      write(16,rec=9) ave_evap                                                                 !mm/d/y
      write(16,rec=10) ave_cres                                                                !s/m
      close(16)
!
! Write wsoil2
! ------------
!
      open(17,file='../outputs/soilm.bin',                                                     !open(17,file='../ipcc_env_vars2/'//& 'soilm_HADCM3_A2_onlyclimate.bin',
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
!
! Write evapotranspiration
! ------------------------
!
      open(18,file='../outputs/evaptr.bin',                                                    !open(18,file='../ipcc_env_vars2/'//& 'evaptr_HADCM3_A2_onlyclimate.bin',
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
!
! Write monthly NPP
! -----------------
!
      open(19,file='../outputs/mon_npp.bin',                                                   !open(19,file='../ipcc_env_vars2/'//& 'mon_npp_HADCM3_A2_onlyclimate.bin',
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
!
! Write monthly photosynthesis
! ----------------------------
!
      open(20,file='../outputs/mon_photo.bin',                                                 !open(20,file='../ipcc_env_vars2/'//& 'mon_photo_HADCM3_A2_onlyclimate.bin',
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
!
! Write monthly plant respiration
! -------------------------------
!
      open(21,file='../outputs/mon_ar.bin',                                                    !open(21,file='../ipcc_env_vars2/'//& 'mon_ar_HADCM3_A2_onlyclimate.bin',
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
!
! Write monthly canopy resistance
! -------------------------------
!
      open(22,file='../outputs/cres.bin',                                                      !open(22,file='../ipcc_env_vars2/'//& 'rc_HADCM3_A2_onlyclimate.bin',
!      open(22,file='../ipcc_env_vars2/'//
!     & 'rc_HADCM3_A2_onlyclimate.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do k=1,12
         do i=1,nx
         do j=1,ny
            waux6(i,j) = cres(i,j,k)
         enddo
         enddo
         write(22,rec=k) waux6
      enddo
      close(22)
      
! Write annual canopy resistence
! ------------------------------
!
      open(23,file='../outputs/ave_cres.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
         write(23,rec=1) ave_cres
      close(23)
!
! Write monthly canopy resistance_pft1
! ------------------------------------
!
      open(24,file='../outputs/cres1.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do k=1,12
         do i=1,nx
         do j=1,ny
            waux7(i,j) = cres1(i,j,k)
         enddo
         enddo
         write(24,rec=k) waux7
      enddo
      close(24)
!
! Write monthly canopy resistance_pft2
! ------------------------------------
!
      open(25,file='../outputs/cres2.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do k=1,12
         do i=1,nx
         do j=1,ny
            waux8(i,j) = cres2(i,j,k)
         enddo
         enddo
         write(25,rec=k) waux8
      enddo
      close(25)
!
! Write monthly canopy resistance_pft3
! ------------------------------------
!
      open(26,file='../outputs/cres3.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do k=1,12
         do i=1,nx
         do j=1,ny
            waux9(i,j) = cres3(i,j,k)
         enddo
         enddo
         write(26,rec=k) waux9
      enddo
      close(26)
!
! Write annual photosynthesis
! ---------------------------
!
      open(27,file='../outputs/mphoto.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
         write(27,rec=1) mphoto
      close(27)
!
! Write anual NPP
! ---------------
!
      open(28,file='../outputs/meanpp.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
         write(28,rec=1) meanpp
      close(28)
!
! Write anual canopy resistence_pft1
! ----------------------------------
!
      open(29,file='../outputs/ave_cres1.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
         write(29,rec=1) ave_cres1
      close(29)
!
! Write anual canopy resistence_pft2
! ----------------------------------
!
      open(30,file='../outputs/ave_cres2.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
         write(30,rec=1) ave_cres2
      close(30)
!
! Write anual canopy resistence_pft3
! ----------------------------------
!
      open(31,file='../outputs/ave_cres3.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
         write(31,rec=1) ave_cres3
      close(31)
!
      open(32,file='../outputs/npp1.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do k=1,12
         do i=1,nx
         do j=1,ny
            waux10(i,j) = npp1(i,j,k)
         enddo
         enddo
         write(32,rec=k) waux10
      enddo
      close(32)
!
      open(33,file='../outputs/npp2.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do k=1,12
         do i=1,nx
         do j=1,ny
            waux11(i,j) = npp2(i,j,k)
         enddo
         enddo
         write(33,rec=k) waux11
      enddo
      close(33)
!
      open(34,file='../outputs/npp3.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do k=1,12
         do i=1,nx
         do j=1,ny
            waux12(i,j) = npp3(i,j,k)
         enddo
         enddo
         write(34,rec=k) waux12
      enddo
      close(34)
!
      open(35,file='../outputs/meanpp1.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
         write(35,rec=1) meanpp1
      close(35)
!
      open(36,file='../outputs/meanpp2.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
         write(36,rec=1) meanpp2
      close(36)
!
      open(37,file='../outputs/meanpp3.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
         write(37,rec=1) meanpp3
      close(37)
!
      open(38,file='../outputs/maresp_caete.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
         write(38,rec=1) maresp
      close(38)
!
      open(39,file='../outputs/wsoil.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do k=1,12
         do i=1,nx
         do j=1,ny
            waux13(i,j) = wsoil(i,j,k)
         enddo
         enddo
         write(39,rec=k) waux13
      enddo
      close(39)
!
      open(40,file='../outputs/runom.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do k=1,12
         do i=1,nx
         do j=1,ny
            waux14(i,j) = runom(i,j,k)
         enddo
         enddo
         write(40,rec=k) waux14
      enddo
      close(40)
!
      open(41,file='../outputs/mnpp1_pon.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
        write(41,rec=1) mnpp1_pon
      close(41)
!
      open(42,file='../outputs/mnpp2_pon.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
        write(42,rec=1) mnpp2_pon
      close(42)
!
       open(43,file='../outputs/mnpp3_pon.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
         write(43,rec=1) mnpp3_pon
      close(43)
!
! Program end
!
      stop
      end
!
!==============================================================================================
!
      subroutine read12(nunit,var)
! auxiliar reading routine
      parameter(nx=192,ny=96)
      integer nunit
      real var(nx,ny,12)
      real aux(nx,ny)
      do k=1,12
        read(nunit,rec=k) aux
        do i=1,nx
          do j=1,ny
            var(i,j,k) = aux(i,j)
          enddo
        enddo
      enddo
      return
      end
!
!==============================================================================================
!
!     if (i.le.96) then                                                                        !inverte hemisferios leste-oeste
!     ni=i+96
!     else
!     ni=(96-i)*(-1)
!     endif
!     nj=(j-96)*(-1)                                                                           !faz flip na latitude
!
!==============================================================================================
!
!      enddo
!      enddo
!      enddo
!
!      return
!      end
