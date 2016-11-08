c234567
      program env
!
!=======================================================================
! CPTEC-PVM2
! Adapted from env4.f.
! Code written by David Lapola
! Last update: Aug/2007
!
! Compile with: g95 (or gfortran) env5.f wbm4.f carbon2.f -o a.exe (a.exe is the executable)
! Execute with ./a.exe
! Then run g95 (gfortran) pvm5.f -o pvm5.exe
! Execute with ./pvm5.exe
!
!=======================================================================
!
! Parameters and variables
! ------------------------
!
      integer nx, ny
      parameter(nx=720,ny=360)
      real lsmk(nx,ny),p0(nx,ny,12)
      real pr(nx,ny,12),t(nx,ny,12)
      real prec(nx,ny,12),temp(nx,ny,12),par(nx,ny,12),ipar(nx,ny,12),ca
      real anpr(nx,ny,12),ant(nx,ny,12)
      real tmin(nx,ny),seanpp(nx,ny),meanpp(nx,ny),meanhr(nx,ny),
     &     meancs(nx,ny),mphoto(nx,ny),maresp(nx,ny),ave_wsoil(nx,ny),
     &     ave_evap(nx,ny),ave_rc(nx,ny),ave_runom(nx,ny),
     &     ave_gsoil(nx,ny),ave_ssoil(nx,ny),ave_snowm(nx,ny)
      real bpot(nx,ny)
      real wsoil2(nx,ny,12),evaptr(nx,ny,12),rc(nx,ny,12),waux(nx,ny),
     & waux2(nx,ny),waux3(nx,ny),waux4(nx,ny),waux5(nx,ny),waux6(nx,ny),
     & waux7(nx,ny),waux8(nx,ny),waux9(nx,ny),waux10(nx,ny),
     & waux11(nx,ny),waux12(nx,ny),waux13(nx,ny),waux14(nx,ny),
     & waux15(nx,ny)
      real wsoil(nx,ny,12),runom(nx,ny,12),wg0(nx,ny,12)
! carbon cycle
      real mp(nx,ny,12),mar(nx,ny,12),mnpp(nx,ny,12),
     &     lai(nx,ny,12),clit(nx,ny,12),csoil(nx,ny,12),
     &     hresp(nx,ny,12)
!
! Open files
! ----------
!

      open( 9,file='../inputs/ps.bin',status='old',
     &        form='unformatted',access='direct',recl=4*nx*ny)
      open(10,file='../inputs/lsmk.bin',status='old',
     &        form='unformatted',access='direct',recl=4*nx*ny)
      open(11,file='../inputs/pr.bin',status='old',
     &        form='unformatted',access='direct',recl=4*nx*ny)
      open(12,file='../inputs/tas.bin',status='old',
     &        form='unformatted',access='direct',recl=4*nx*ny)
      open(13,file='../inputs/rsds.bin',status='old',
     &        form='unformatted',access='direct',recl=4*nx*ny)
!
! Future anomalies, see lines 85 and 86
!      open(14,file='../ipcc_ar4/anomalias/SRESA2/temp/'//
!     &               'HADCM3_A2_2070_2099_anom_t.bin',status='old',
!     &        form='unformatted',access='direct',recl=4*nx*ny)
!      open(15,file='../ipcc_ar4/anomalias/SRESA2/prec/'//
!     &               'HADCM3_A2_2070_2099_anom_pr.bin',status='old',
!     &        form='unformatted',access='direct',recl=4*nx*ny)
!
! Read data
! ---------
!
         !mb
      read (10,rec=1) lsmk      !land = 1; ocean = 0
      call read12 (9, p0)  !Pa
      call read12 (11,pr)  !mm/month
      call read12 (12,t)  !oC
      call read12 (13,ipar)  !w/m2
!      call read12 (14,ant)  !oC
!      call read12 (15,anpr)  !oC
!
! Close files
c -----------
c
      close ( 9)
      close (10)
      close (11)
      close (12)
      close (13)
!      close (14)
!      close (15)
!
      do k=1,12
         do i=1,nx
         do j=1,ny
! photosynthetically active radiation reaching canopy
! (ipar ; Ein/m2/s) [Eq. 7] observed data from ISLSCP2
      par(i,j,k) = (ipar(i,j,k)/(2.18e5))!converting to Ein/m2/s
!
      temp(i,j,k) = t(i,j,k) !+ant(i,j,k) !uncomment to use future anomalies
      prec(i,j,k) = pr(i,j,k)!+anpr(i,j,k) !+pr(i,j,k)*0.2 !uncomment to use future anomalies
      if (prec(i,j,k).lt.0.0) prec (i,j,k) = 0.0
!
! atmospheric CO2 pressure (Pa) ! ppmv / Pa
! 1 Pa CO2 = 9.901 ppmv CO2 (Adams et al. 2004)
!      ca= 18.18 !Pa (=180 ppmv; Last Glacial Maximum)
!      ca= 28.28 !Pa (=280 ppmv; Pre-Industrial Rev.)
!      ca= 35.35 !Pa (=350 ppmv; 1961-1990)
      ca= 363/9.901 !Pa (=350 ppmv; 1961-1990)
!      ca= 54.03 !Pa (=535 ppmv; SRES-B1 2080's)
!      ca= 73.73 !Pa (=730 ppmv; SRES-A2 2080's)
!      ca= ((73.73-35.35)*0.5)+35.35 !Pa half effect!!!

!
         enddo
         enddo
      enddo

c      do i=1,nx
c      do j=1,ny
c        p0(i,j) = 944.59 !mb             ! media obtida pelo input do ILAMB
c      enddo
c      enddo
      
! Calculate environmental variables (water balance model)
! -------------------------------------------------------
!
      call wbm (prec,temp,lsmk,p0,ca,par,
     &          tmin,meanpp,seanpp,mphoto,maresp,  !carbon
     &          meanhr,meancs,wsoil2,evaptr,mnpp,mp,mar,rc,
     &          ave_wsoil,ave_evap,ave_rc,ave_runom,ave_gsoil,
     &          ave_ssoil,ave_snowm,wsoil,runom,wg0)
!
! write environmental variables
!
!      open(16,file='../outputs/ambientais4.bin',
!      open(16,file='../ipcc_env_vars2/'//
!     & 'ambientais4_HADCM3_A2_onlyclimate.bin',
!     &        status='unknown',form='unformatted',
!     &        access='direct',recl=4*nx*ny)
!      write(16,rec=1) tmin   !oC
!      write(16,rec=2) meanpp !kgC/m2/y
!      write(16,rec=3) seanpp !dimensionless
!      write(16,rec=4) meanhr !kgC/m2/y
!      write(16,rec=5) meancs !kgC/m2/y
!      write(16,rec=6) mphoto !kgC/m2/y
!      write(16,rec=7) maresp !kgC/m2/y
!      write(16,rec=8) ave_wsoil !mm
!      write(16,rec=9) ave_evap !mm/d/y
!      write(16,rec=10) ave_rc !s/m
!      close(16)
c
c write wsoil2
      open(17,file='../outputs/wsoil2.flt',
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
      open(18,file='../outputs/et.flt',
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
!      open(19,file='../outputs/mon_npp.flt',
!      open(19,file='../ipcc_env_vars2/'//
!     & 'mon_npp_HADCM3_A2_onlyclimate.bin',
!     &        status='unknown',form='unformatted',
!     &        access='direct',recl=4*nx*ny)
!      do k=1,12
!         do i=1,nx
!         do j=1,ny
!            waux3(i,j) = mnpp(i,j,k)
!         enddo
!         enddo
!         write(19,rec=k) waux3
!      enddo
!      close(19)
      
c write monthly Photosynthesis
!      open(20,file='../outputs/mon_photo.flt',
!      open(20,file='../ipcc_env_vars2/'//
!     & 'mon_photo_HADCM3_A2_onlyclimate.bin',
!     &        status='unknown',form='unformatted',
!     &        access='direct',recl=4*nx*ny)
!      do k=1,12
!         do i=1,nx
!         do j=1,ny
!            waux4(i,j) = mp(i,j,k)
!         enddo
!         enddo
!         write(20,rec=k) waux4
!      enddo
!      close(20)
c
c write monthly plant respiration
!      open(21,file='../outputs/mon_ar.flt',
!      open(21,file='../ipcc_env_vars2/'//
!     & 'mon_ar_HADCM3_A2_onlyclimate.bin',
!     &        status='unknown',form='unformatted',
!     &        access='direct',recl=4*nx*ny)
!      do k=1,12
!         do i=1,nx
!         do j=1,ny
!            waux5(i,j) = mar(i,j,k)
!         enddo
!         enddo
!         write(21,rec=k) waux5
!      enddo
!      close(21)
c
c write canopy resistance
!      open(22,file='../outputs/rc.flt',
!      open(22,file='../ipcc_env_vars2/'//
!     & 'rc_HADCM3_A2_onlyclimate.bin',
!     &        status='unknown',form='unformatted',
!     &        access='direct',recl=4*nx*ny)
!      do k=1,12
!         do i=1,nx
!         do j=1,ny
!            waux6(i,j) = rc(i,j,k)
!         enddo
!         enddo
!         write(22,rec=k) waux6
!      enddo
!      close(22)
      
!      open(23,file='../outputs/ave_runom.flt',
!     &        status='unknown',form='unformatted',
!     &        access='direct',recl=4*nx*ny)
!         do i=1,nx
!         do j=1,ny
!            waux7(i,j) = ave_runom(i,j)
!         enddo
!         enddo
!         write(23,rec=1) waux7
!      close(23)
      
!      open(24,file='../outputs/ave_gsoil.flt',
!     &        status='unknown',form='unformatted',
!     &        access='direct',recl=4*nx*ny)
!         do i=1,nx
!         do j=1,ny
!            waux8(i,j) = ave_gsoil(i,j)
!         enddo
!         enddo
!         write(24,rec=1) waux8
!      close(24)


!      open(25,file='../outputs/ave_ssoil.flt',
!     &        status='unknown',form='unformatted',
!     &        access='direct',recl=4*nx*ny)
!         do i=1,nx
!         do j=1,ny
!            waux9(i,j) = ave_ssoil(i,j)
!         enddo
!         enddo
!         write(25,rec=1) waux9
!      close(25)

!      open(26,file='../outputs/ave_snowm.flt',
!     &        status='unknown',form='unformatted',
!     &        access='direct',recl=4*nx*ny)
!         do i=1,nx
!         do j=1,ny
!            waux10(i,j) = ave_snowm(i,j)
!         enddo
!         enddo
!         write(26,rec=1) waux10
!      close(26)
      
!      open(27,file='../outputs/ave_wsoil.flt',
!     &        status='unknown',form='unformatted',
!     &        access='direct',recl=4*nx*ny)
!         do i=1,nx
!         do j=1,ny
!            waux11(i,j) = ave_wsoil(i,j)
!         enddo
!         enddo
!         write(27,rec=1) waux11
!      close(27)
      
!      open(28,file='../outputs/ave_evap.flt',
!     &        status='unknown',form='unformatted',
!     &        access='direct',recl=4*nx*ny)
!         do i=1,nx
!         do j=1,ny
!            waux12(i,j) = ave_evap(i,j)
!         enddo
!         enddo
!         write(28,rec=1) waux12
!      close(28)

      
      open(29,file='../outputs/wsoil.flt',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
        do k=1,12 
         do i=1,nx
         do j=1,ny
            waux13(i,j) = wsoil(i,j,k)
         enddo
         enddo
         write(29,rec=k) waux13
         enddo
      close(29)

      open(30,file='../outputs/runoff.flt',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do k=1,12
         do i=1,nx
         do j=1,ny
            waux14(i,j) = runom(i,j,k)
         enddo
         enddo
         write(30,rec=k) waux14
         enddo
      close(30)

!      open(31,file='../outputs/wg0.nc',
!     &        status='unknown',form='unformatted',
!     &        access='direct',recl=4*nx*ny)
!      do k=1,12
!         do i=1,nx
!         do j=1,ny
!            waux15(i,j) = wsoil(i,j,k)
!         enddo
!         enddo
!         write(31,rec=k) waux15
!      enddo
!      close(31)
c
      stop
      end
c=======================================================================

      subroutine read12(nunit,var)
c auxiliar reading routine
      parameter(nx=720,ny=360)
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
      
!      subroutine xrev(
!      do i=1,nx
!      do j=1,ny
!         if (i.le.96) then
!         pr2(i,j) = pr
