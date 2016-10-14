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
      integer nx, ny
      parameter(nx=192,ny=96)
      real lsmk(nx,ny),p0(nx,ny)
      real pr(nx,ny,12),t(nx,ny,12)
      real prec(nx,ny,12),temp(nx,ny,12),par(nx,ny,12),ipar(nx,ny,12),ca
      real anpr(nx,ny,12),ant(nx,ny,12)
      real tmin(nx,ny),seanpp(nx,ny),meanpp(nx,ny),meanhr(nx,ny),
     &     meancs(nx,ny),mphoto(nx,ny),maresp(nx,ny),ave_wsoil(nx,ny),
     &     ave_evap(nx,ny),ave_rc(nx,ny)
      real bpot(nx,ny)
      real wsoil2(nx,ny,12),evaptr(nx,ny,12),rc(nx,ny,12),waux(nx,ny),
     & waux2(nx,ny),waux3(nx,ny),waux4(nx,ny),waux5(nx,ny),waux6(nx,ny)
c carbon cycle
      real mp(nx,ny,12),mar(nx,ny,12),mnpp(nx,ny,12),
     &     lai(nx,ny,12),clit(nx,ny,12),csoil(nx,ny,12),
     &     hresp(nx,ny,12)
c
c Open files
c ----------
c
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
c
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
c
      read ( 9,rec=1) p0     !mb
      read (10,rec=1) lsmk   !land = 1; ocean = 0
      call read12 (11,pr)  !mm/month
      call read12 (12,t)  !oC
      call read12 (13,ipar)  !w/m2
!      call read12 (14,ant)  !oC
!      call read12 (15,anpr)  !oC
c
c Close files
c -----------
c
      close ( 9)
      close (10)
      close (11)
      close (12)
      close (13)
!      close (14)
!      close (15)
c
      do k=1,12
         do i=1,nx
         do j=1,ny
c photosynthetically active radiation reaching canopy
c (ipar ; Ein/m2/s) [Eq. 7] observed data from ISLSCP2
      par(i,j,k) = ipar(i,j,k)/(2.18e5) !converting to Ein/m2/s
c
      temp(i,j,k) = t(i,j,k) !+ant(i,j,k) !uncomment to use future anomalies
      prec(i,j,k) = pr(i,j,k)!+anpr(i,j,k) !+pr(i,j,k)*0.2 !uncomment to use future anomalies
      if (prec(i,j,k).lt.0.0) prec (i,j,k) = 0.0
c
c atmospheric CO2 pressure (Pa) ! ppmv / Pa
c 1 Pa CO2 = 9.901 ppmv CO2 (Adams et al. 2004)
!      ca= 18.18 !Pa (=180 ppmv; Last Glacial Maximum)
!      ca= 28.28 !Pa (=280 ppmv; Pre-Industrial Rev.)
!      ca= 35.35 !Pa (=350 ppmv; 1961-1990)
      ca= 350/9.901 !Pa (=350 ppmv; 1961-1990)
!      ca= 54.03 !Pa (=535 ppmv; SRES-B1 2080's)
!      ca= 73.73 !Pa (=730 ppmv; SRES-A2 2080's)
!      ca= ((73.73-35.35)*0.5)+35.35 !Pa half effect!!!
c
         enddo
         enddo
      enddo

c Calculate environmental variables (water balance model)
c -------------------------------------------------------
c
      call wbm (prec,temp,lsmk,p0,ca,par,
     &          tmin,meanpp,seanpp,mphoto,maresp,  !carbon
     &          meanhr,meancs,wsoil2,evaptr,mnpp,mp,mar,rc,
     &          ave_wsoil,ave_evap,ave_rc)
c
c write environmental variables
      open(16,file='../outputs/ambientais4.bin',
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
      close(16)
c
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
            var(i,j,k) = aux(i,j)
          enddo
        enddo
      enddo
      return
      end
