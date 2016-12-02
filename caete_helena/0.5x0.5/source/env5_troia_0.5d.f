c234567
!
      program env
!
! =======================================================================
! CPTEC-PVM2
! Adapted from env4.f.
! Code written by David Lapola and Helena Alves do Prado
! Last update: Aug/2007
!
! Compile with: g95 (or gfortran) env5.f wbm4.f productivity1.f -o a.exe
! Execute with ./a.exe
! Then run g95 (gfortran) pvm5.f -o pvm5.exe
! Execute with ./pvm5.exe
! =======================================================================
!
! Parameters and variables
! ------------------------
!
      integer nx,ny,i,j,k
      parameter (nx=720,ny=360)
      real,parameter :: no_data = -9999.0
!
! Inputs
! ------
!
      real ca                                                                                 
      real lsmk(nx,ny)                                                                            
      real p0(nx,ny,12)                                                                  
      real prec(nx,ny,12)                                                                     
      real pr(nx,ny,12)                                                                      
      real anpr(nx,ny,12)                                                                            
      real temp(nx,ny,12)                                                                
      real t(nx,ny,12)                                                                   
      real ant(nx,ny,12)                                                                
      real par(nx,ny,12)                                                                
      real ipar(nx,ny,12) 
!
! Outputs
! -------
!
      real tmin(nx,ny)                                                                   
      real meanpp(nx,ny)                                                                             
      real meanpp1(nx,ny)                                                                   
      real meanpp2(nx,ny)                                                                     
      real meanpp3(nx,ny) 
      real seanpp(nx,ny)                                                                                
      real mphoto(nx,ny)                                                                                                                                             
      real maresp(nx,ny)                                                                  
      real meanhr(nx,ny)                                          
      real meancs(nx,ny)                                             
      real ave_wsoil(nx,ny)                                                                      
      real ave_evap(nx,ny)                                                        
      real ave_rc(nx,ny)   
!
      real wsoil2(nx,ny,12)
      real evaptr(nx,ny,12)
      real mnpp(nx,ny,12)
      real npp1(nx,ny,12)                                                               
      real npp2(nx,ny,12)                                                                 
      real npp3(nx,ny,12)
      real mp(nx,ny,12)
      real mar(nx,ny,12)
      real rc(nx,ny,12)
      real wsoil(nx,ny,12)                                                                    
      real runom(nx,ny,12)                                                                    
!
! Auxiliars
! ---------
!
      real waux(nx,ny)                                                                                      
      real waux2(nx,ny)                                                           
      real waux3(nx,ny)                                                              
      real waux4(nx,ny)                                                                                              
      real waux5(nx,ny)                                                                                              
      real waux6(nx,ny)                                                                                        
!
! Internal variables
! ------------------
!
      real lai(nx,ny,12)
      real clit(nx,ny,12)
      real csoil(nx,ny,12)
      real hresp(nx,ny,12)
!
! Open files
! ==========
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
c ===================================================================
c
c Future anomalies see lines 85 and 86
c ------------------------------------
c
c     open(14,file='../ipcc_ar4/anomalias/SRESA2/temp/'//
c    &               'HADCM3_A2_2070_2099_anom_t.bin',status='old',
c    &        form='unformatted',access='direct',recl=4*nx*ny)
c     open(15,file='../ipcc_ar4/anomalias/SRESA2/prec/'//
c    &               'HADCM3_A2_2070_2099_anom_pr.bin',status='old',
c    &        form='unformatted',access='direct',recl=4*nx*ny)
c
c ===================================================================
!
! Read data
! =========
!
      !read ( 9,rec=1) p0                                                                       
      read (10,rec=1) lsmk
      call read12 (9,p0)  ! getting p0 as (nx,ny,12)                                                              
      call read12 (11,pr)                                                                 
      call read12 (12,t)                                                                  
      call read12 (13,ipar)                                                               
c      call read12 (14,ant)                                                                 
c      call read12 (15,anpr)                                                                  
!
! Close files
! ===========
!
      close ( 9)
      close (10)
      close (11)
      close (12)
      close (13)
c      close (14)
c      close (15)
!
! ===========
!
      do k=1,12
         do i=1,nx
         do j=1,ny

! Photosynthetically active radiation reaching canopy (IPAR:Ein/m2/s)
! Observed data from ISLSCP2

         par(i,j,k) = (ipar(i,j,k)/(2.18e5))                                                   
         temp(i,j,k) = t(i,j,k)                                                              
         prec(i,j,k) = pr(i,j,k)                                                               
         if (prec(i,j,k).lt.0.0) prec (i,j,k) = 0.0
!
! atmospheric CO2 pressure (Pa) ! ppmv / Pa
! (1 Pa CO2 = 9.901 ppmv CO2 (Adams et al. 2004))
c      ca= 18.18 !Pa (=180 ppmv; Last Glacial Maximum)
c      ca= 28.28 !Pa (=280 ppmv; Pre-Industrial Rev.)
c      ca= 35.35 !Pa (=350 ppmv; 1961-1990)
      ca= 363/9.901 !Pa (=350 ppmv; 1961-1990)   ! changing pCO2 to mean in 1981-2010 
c      ca= 54.03 !Pa (=535 ppmv; SRES-B1 2080's)
c      ca= 73.73 !Pa (=730 ppmv; SRES-A2 2080's)
c      ca= ((73.73-35.35)*0.5)+35.35 !Pa half effect!!!                                             
!
      enddo
      enddo
      enddo
!
! ===================================================================
!
! Calculate environmental variables (wbm)
! =======================================
!
      call wbm (prec,temp,lsmk,p0,ca,par,                                                    
     &          tmin,meanpp,seanpp,mphoto,maresp,                                          
     &          meanhr,meancs,wsoil2,evaptr,mnpp,
     &          mp,mar,rc,ave_wsoil,ave_evap,ave_rc,
     &          npp1,npp2,npp3,meanpp1,meanpp2,meanpp3)
!     
!
! Write environmental variables
! =============================
!
      open(16,file='../outputs/ambientais4.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      write(16,rec=1) tmin 
      write(16,rec=2) meanpp 
      write(16,rec=3) seanpp 
      write(16,rec=4) meanhr 
      write(16,rec=5) meancs 
      write(16,rec=6) mphoto 
      write(16,rec=7) maresp 
      write(16,rec=8) ave_wsoil 
      write(16,rec=9) ave_evap 
      write(16,rec=10) ave_rc 
      close(16)
!
! Write wsoil2 ---> mm == kg m-2
! ------------
!
c      open(17,file='../outputs/mrso.bin',
c     &        status='unknown',form='unformatted',
c     &        access='direct',recl=4*nx*ny)
c      do k=1,12
c         do i=1,nx
c         do j=1,ny
c            waux(i,j) = wsoil2(i,j,k)
c         enddo
c         enddo
c         write(17,rec=k) waux
c      enddo
c      close(17)
!
! Write evapotranspiration  mm/day == kg m-2 day-1
! ------------------------
!
      open(18,file='../outputs/et.bin',
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
! Write monthly NPP kg C m-2 year-1
! -----------------
!
      open(19,file='../outputs/npp.bin',
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
! Write monthly photosynthesis kg C m-2 year-1
! ----------------------------
!
      open(20,file='../outputs/ph.bin',
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
! Write monthly plant respiration  kg C m-2 y-1
! -------------------------------
!
      open(21,file='../outputs/ar.bin',
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
! Write monthly canopy resistence s m-1
! -------------------------------
!
      open(22,file='../outputs/rc.bin',
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
!
! Write annual canopy resistence 
! ------------------------------
!
!      open(23,file='../outputs/ave_rc_caete_new.flt',
!     &        status='unknown',form='unformatted',
!     &        access='direct',recl=4*nx*ny)
!         write(23,rec=1) ave_rc
!      close(23)
!
! Write annual npp 
! ----------------
!
!      open(24,file='../outputs/ave_npp_caete_brandnew.flt',
!     &        status='unknown',form='unformatted',
!     &        access='direct',recl=4*nx*ny)
!         write(24,rec=1) meanpp
!      close(24)
!
!     Program end
!
      stop
      end
!
!====================================================================
!
      subroutine read12(nunit,var)
! auxiliar reading routine
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
!
c ===================================================================
c
c     if (i.le.96) then                                                                        !inverte hemisferios leste-oeste
c     ni=i+96
c     else
c     ni=(96-i)*(-1)
c     endif
c     nj=(j-96)*(-1)                                                                           !faz flip na latitude
c
c ===================================================================
