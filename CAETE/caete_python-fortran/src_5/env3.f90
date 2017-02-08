program env
  use global_pars
  implicit none
  !     
  !     =======================================================================
  !     CPTEC-PVM2
  !     Adapted from env4.f.
  !     Code written by David Lapola and Helena Alves do Prado e Bianca Rius & jpdarela
  !     Last update: jan/2017
  !     Compile with: gfortran commom.f90 env3.f90 prod3.f90 
  !     Execute with ./a.exe
  !     =======================================================================
  
  !     Parameters and variables
  !     ------------------------     
  integer(kind=i4) :: i,j,k,p
  integer(kind=i4),parameter :: q = npls
  integer(kind=i4),parameter :: nt = ntimes
  real(kind=r4),parameter :: no_data = -9999.0
  !     
  !   Model  Inputs
  !   -------------  
  real(kind=r4) :: ca                   !CO2 concentration (Pa)
  real(kind=r4) :: lsmk(nx,ny)          !Land=1/Ocean=0
  real(kind=r4) :: p0(nt)         !Atmospheric pressure (mb)
  real(kind=r4) :: ps(nx,ny,nt)         ! auxiliar to read atm pressure information
  real(kind=r4) :: prec(nt)       !Precipitation (mm/month)
  real(kind=r4) :: pr(nx,ny,nt)         !Auxiliar_precipitation (mm/month)
  real(kind=r4) :: temp(nt)       !Temperature (oC)
  real(kind=r4) :: t(nx,ny,nt)          !Auxiliar_temperature (oC)
  real(kind=r4) :: par(nt)        !Incident photosynthetic active radiation (Ein/m2/s)
  real(kind=r4) :: ipar(nx,ny,nt)       !Auxiliar_incident photosynthetic active radiation (w/m2)
  real(kind=r4) :: rhs(nt)        !Relative humidity
  real(kind=r4) :: rhaux(nx,ny,nt)      !RHS auxiliar
  
  real(kind=r4),dimension(q) :: cleafin = 0.0
  real(kind=r4),dimension(q) :: cawoodin = 0.0
  real(kind=r4),dimension(q) :: cfrootin = 0.0
  
  !    Model Outputs
  !    -------------
  real(kind=r4),dimension(nt) :: emaxm = 0.0
  real(kind=r4),dimension(nt) :: tsoil = 0.0
  
  real(kind=r4),dimension(nt,q) :: photo_pft = 0.0 !Monthly photosynthesis   (kgC/m2)
  real(kind=r4),dimension(nt,q) :: aresp_pft = 0.0 !Monthly autotrophic res  (kgC/m2)
  real(kind=r4),dimension(nt,q) :: rm_pft    = 0.0
  real(kind=r4),dimension(nt,q) :: rg_pft    = 0.0
  real(kind=r4),dimension(nt,q) :: npp_pft   = 0.0  !Monthly net primary produ (kgC/m2)
      
  real(kind=r4),dimension(nt,q) :: lai_pft = 0.0  !Monthly leaf area index
  real(kind=r4),dimension(nt,q) :: clit_pft = 0.0 !Monthly litter carbon
  real(kind=r4),dimension(nt,q) :: csoil_pft = 0.0 !Monthly soil carbon
  real(kind=r4),dimension(nt,q) :: hresp_pft = 0.0 !Monthly het resp          (kgC/m2)
  real(kind=r4),dimension(nt,q) :: rcm_pft = 0.0 
  
  real(kind=r4),dimension(nt,q) :: runom_pft = 0.0 !Runoff
  real(kind=r4),dimension(nt,q) :: evapm_pft = 0.0 !Actual evapotranspiration        
  real(kind=r4),dimension(nt,q) :: wsoil_pft = 0.0 !Soil moisture (mm)
  
  real(kind=r4),dimension(q) :: cleaf_pft = 0.0
  real(kind=r4),dimension(q) :: cawood_pft = 0.0
  real(kind=r4),dimension(q) :: cfroot_pft = 0.0
  real(kind=r4),dimension(q) :: gridcell_ocp = 0.0 !  final grid cell occupation for each pft (percentage of area)

  real(kind=r4),dimension(nx,ny,q) :: grd_ocp = 0.0
  real(kind=r4),dimension(nx,ny,q) :: clini = 0.0
  real(kind=r4),dimension(nx,ny,q) :: cfini = 0.0
  real(kind=r4),dimension(nx,ny,q) :: cwini = 0.0
  real(kind=r4),dimension(nx,ny,q) :: clfim = 0.0
  real(kind=r4),dimension(nx,ny,q) :: cffim = 0.0
  real(kind=r4),dimension(nx,ny,q) :: cwfim = 0.0


  !     variaveis do spinup
  real(kind=r4) npp_sca 
  real(kind=r4),dimension(nx,ny,nt) :: npp_pot
  real(kind=r4),dimension(q) :: aux1, aux2, aux3
  real(kind=r4),dimension(nx,ny) :: aux_npp 

  
  real(kind=r4),dimension(nx,ny) :: ave_ph = 0.0
  real(kind=r4),dimension(nx,ny) :: ave_ar = 0.0
  real(kind=r4),dimension(nx,ny) :: ave_npp = 0.0
  real(kind=r4),dimension(nx,ny) :: ave_lai = 0.0
  real(kind=r4),dimension(nx,ny) :: ave_clit = 0.0
  real(kind=r4),dimension(nx,ny) :: ave_cs = 0.0
  real(kind=r4),dimension(nx,ny) :: ave_hr = 0.0
  real(kind=r4),dimension(nx,ny) :: ave_rc = 0.0
  real(kind=r4),dimension(nx,ny) :: ave_runom = 0.0
  real(kind=r4),dimension(nx,ny) :: ave_evap = 0.0
  real(kind=r4),dimension(nx,ny) :: ave_wsoil = 0.0    

      
  real(kind=r4), dimension(nx,ny,nt) :: ph = 0.0
  real(kind=r4), dimension(nx,ny,nt) :: ar = 0.0
  real(kind=r4), dimension(nx,ny,nt) :: npp = 0.0
  real(kind=r4), dimension(nx,ny,nt) :: lai = 0.0
  real(kind=r4), dimension(nx,ny,nt) :: clit = 0.0
  real(kind=r4), dimension(nx,ny,nt) :: csoil = 0.0
  real(kind=r4), dimension(nx,ny,nt) :: hr = 0.0
  real(kind=r4), dimension(nx,ny,nt) :: rcm = 0.0
  real(kind=r4), dimension(nx,ny,nt) :: evaptr = 0.0
  real(kind=r4), dimension(nx,ny,nt) :: wsoil = 0.0
  real(kind=r4), dimension(nx,ny,nt) :: runom = 0.0
  
  real(kind=r4), dimension(nx,ny,nt) :: rm = 0.0
  real(kind=r4), dimension(nx,ny,nt) :: rg = 0.0

 ! initialize arrays  
  
  !     Open INPUT files
  !     ==========
  open( 9,file='./inputs/lsmk.bin',status='old',form='unformatted'&
       &,access='direct',recl=4*nx*ny)

  open(10,file='./inputs/ps.bin',status='old',form='unformatted'&
       &,access='direct',recl=4*nx*ny)

  open(11,file='./inputs/pr.bin',status='old',form='unformatted'&
       &,access='direct',recl=4*nx*ny)
      
  open(12,file='./inputs/tas.bin',status='old',form='unformatted'&
       &,access='direct',recl=4*nx*ny)
      
  open(13,file='./inputs/rsds.bin',status='old',form='unformatted'&
       &,access='direct',recl=4*nx*ny)

  open(14,file='./inputs/hurs.bin',status='old',form='unformatted'&
       &,access='direct',recl=4*nx*ny)

  open(26,file='./inputs/npp.bin',status='old',form='unformatted'&
       &,access='direct',recl=4*nx*ny)

  !     Read data
  !     =========
      
  read (9,rec=1) lsmk
  !c     read (26,rec=1) aux_npp 
      
  call readx(10,ps,12)
  call readx(11,pr,12)
  call readx(12,t,12)
  call readx(13,ipar,12)
  call readx(14,rhaux,12)
  call readx(26,npp_pot,12)      
  !     Close files
  !     ===========
  !     
  close( 9)
  close(10)
  close(11)
  close(12)
  close(13)
  close(14)
  close(26)



! apply mask
  do i = 1,nx
     do j = 1,ny
        do k = 1,nt
           if(nint(lsmk(i,j)) .eq. 0) then
              ph(i,j,k) = no_data
              ar(i,j,k) = no_data
              npp(i,j,k) = no_data
              lai(i,j,k) = no_data
              clit(i,j,k) = no_data
              csoil(i,j,k) = no_data
              hr(i,j,k) = no_data
              rcm(i,j,k) = no_data
              runom(i,j,k) = no_data
              evaptr(i,j,k) = no_data
              wsoil(i,j,k) = no_data
              rm(i,j,k)  = no_data
              rg(i,j,k)  = no_data
           endif
        enddo
     enddo
  enddo

  do i = 1,nx
     do j = 1,ny
        do p = 1,q
           if(nint(lsmk(i,j)) .eq. 0) then
              clfim(i,j,p)   = no_data
              cffim(i,j,p)   = no_data
              cwfim(i,j,p)   = no_data
              grd_ocp(i,j,p) = no_data
              clini(i,j,p)   = no_data
              cfini(i,j,p)   = no_data
              cwini(i,j,p)   = no_data
           endif
        enddo
     enddo
  enddo

  
  !c     Calculating annual npp
  do i =1,nx
     do j=1,ny
        if(nint(lsmk(i,j)) .eq. 1) then 
           aux_npp(i,j) = 0.0
           do k = 1,nt
              aux_npp(i,j) = aux_npp(i,j) + (npp_pot(i,j,k)/real(nt)) 
           enddo
        else
           aux_npp(i,j) = no_data
        endif
     enddo
  enddo

  
  !     Atmospheric CO2 pressure (Pa) !Ppmv / Pa
  ca= 363/9.901             !Pa (=363 ppmv; 1981-2010)

!$OMP PARALLEL DO SCHEDULE(STATIC),PRIVATE(I,J,K,P),ORDERED,DEFAULT(SHARED) 
  do i=1,nx
     if(mod(I,72) .eq. 0) print*, nint(real(i)/real(nx) * 100.)
     do j=1,ny
        npp_sca = 0.0
        ! for land grid cells only 
        if(nint(lsmk(i,j)) .eq. 1) then
           
           do k = 1,nt
              ph(i,j,k) = 0.0
              ar(i,j,k) = 0.0
              npp(i,j,k) = 0.0
              lai(i,j,k) = 0.0
              clit(i,j,k) = 0.0
              csoil(i,j,k) = 0.0
              hr(i,j,k) = 0.0
              rcm(i,j,k) = 0.0
              runom(i,j,k) = 0.0
              evaptr(i,j,k) = 0.0
              wsoil(i,j,k) = 0.0
              rm(i,j,k)  = 0.0
              rg(i,j,k)  = 0.0
           enddo
           
           ! run spinup
           npp_sca = aux_npp(i,j)
           
           do p=1,q   
              aux1(p) = 0.0
              aux2(p) = 0.0
              aux3(p) = 0.0
           enddo
           
           call spinup(npp_sca, aux1, aux2, aux3)
           
           do p=1,q
              cleafin(p)  = aux1(p)
              cfrootin(p) = aux2(p)
              cawoodin(p) = aux3(p)
           enddo
           !     ------------------------------------------
           do k=1,nt
              do p=1,q
                 photo_pft(k,p)  = 0.0
                 aresp_pft(k,p)  = 0.0
                 npp_pft(k,p)    = 0.0
                 lai_pft(k,p)    = 0.0
                 clit_pft(k,p)   = 0.0
                 csoil_pft(k,p)  = 0.0
                 hresp_pft(k,p)  = 0.0
                 rcm_pft(k,p)    = 0.0
                 runom_pft(k,p)  = 0.0
                 evapm_pft(k,p)  = 0.0
                 wsoil_pft(k,p)  = 0.0
                 rm_pft(k,p)     = 0.0
                 rg_pft(k,p)     = 0.0
                 cleaf_pft(p)    = 0.0
                 cawood_pft(p)   = 0.0
                 cfroot_pft(p)   = 0.0
                 gridcell_ocp(p) = 0.0
              enddo
              rhs (k) = rhaux(i,j,k) 
              par (k) = ipar  (i,j,k)
              temp(k) = t     (i,j,k) 
              p0  (k) = ps    (i,j,k) 
              prec(k) = pr    (i,j,k) 
           enddo
           print*, 'calling wbm:', i,j
           call wbm (prec,temp,p0,ca,par,rhs,cleafin,cawoodin,cfrootin,&
                &    emaxm, tsoil, photo_pft,aresp_pft,npp_pft,lai_pft,&
                &    clit_pft,csoil_pft, hresp_pft,rcm_pft,runom_pft,&
                &    evapm_pft,wsoil_pft,rm_pft,rg_pft,cleaf_pft,cawood_pft,&
                &    cfroot_pft,gridcell_ocp)
!!$
!!$subroutine wbm (prec,temp,p0,ca,par,rhs,cleaf_ini,cawood_ini&
!!$     &,cfroot_ini,emaxm, tsoil, photo_pft,aresp_pft,npp_pft,lai_pft&
!!$     &,clit_pft,csoil_pft, hresp_pft,rcm_pft,runom_pft,evapm_pft&
!!$     &,wsoil_pft,rm_pft,rg_pft,cleaf_pft,cawood_pft,cfroot_pft,grid_area)
           
           do p = 1,q
              clfim(i,j,p) = cleaf_pft(p)
              cffim(i,j,p) = cfroot_pft(p)
              cwfim(i,j,p) = cawood_pft(p)
              grd_ocp(i,j,p) = gridcell_ocp(p)
              clini(i,j,p) = cleafin(p)
              cfini(i,j,p) = cfrootin(p)
              cwini(i,j,p) = cawoodin(p)
              do k = 1,nt
                 ph(i,j,k) = ph(i,j,k) + photo_pft(k,p)
                 ar(i,j,k) = ar(i,j,k) + aresp_pft(k,p)
                 npp(i,j,k) = npp(i,j,k) + npp_pft(k,p)
                 lai(i,j,k) = lai(i,j,k) + lai_pft(k,p)
                 clit(i,j,k) = clit(i,j,k) + clit_pft(k,p)
                 csoil(i,j,k) = csoil(i,j,k) + csoil_pft(k,p)
                 hr(i,j,k) = hr(i,j,k) + hresp_pft(k,p)
                 rcm(i,j,k) = rcm(i,j,k) + rcm_pft(k,p)
                 runom(i,j,k) = runom(i,j,k) + runom_pft(k,p)
                 evaptr(i,j,k) = evaptr(i,j,k) + evapm_pft(k,p)
                 wsoil(i,j,k) = wsoil(i,j,k) + wsoil_pft(k,p)
                 rm(i,j,k)  = rm(i,j,k) + rm_pft(k,p)
                 rg(i,j,k)  = rg(i,j,k) + rg_pft(k,p)
              enddo
           enddo
        endif
     enddo           
  enddo
!$OMP END PARALLEL DO 

  open(10,file='./spinup/gridcell_ocp.bin',&
       &    status='unknown',form='unformatted',&
       &    access='direct',recl=4*nx*ny)
  call savex(10, grd_ocp, q)
  
  open(10,file='./spinup/cleaf.bin',&
       &    status='unknown',form='unformatted',&
       &    access='direct',recl=4*nx*ny)
  call savex(10, clfim, q)
  
  open(10,file='./spinup/cawood.bin',&
       &    status='unknown',form='unformatted',&
       &    access='direct',recl=4*nx*ny)
  call savex(10,cwfim,q)
  
  open(10,file='./spinup/cfroot.bin',&
       &    status='unknown',form='unformatted',&
       &    access='direct',recl=4*nx*ny)
  call savex(10, cffim,q)
  
  open(10,file='./spinup/clini.bin',&
       &    status='unknown',form='unformatted',&
       &    access='direct',recl=4*nx*ny)
  call savex(10, clini, q)
  
  open(10,file='./spinup/cwini.bin',&
       &    status='unknown',form='unformatted',&
       &    access='direct',recl=4*nx*ny)
  call savex(10,cwini,q)
  
  open(10,file='./spinup/cfini.bin',&
       &    status='unknown',form='unformatted',&
       &    access='direct',recl=4*nx*ny)
  call savex(10, cfini,q)
  
  
  do i = 1,nx
     do j = 1,ny
        if(nint(lsmk(i,j)) .ne. 0) then
           do k = 1,12
              ave_ph(i,j) = ave_ph(i,j) + ph(i,j,k)/12.
              ave_ar(i,j) = ave_ar(i,j) + ar(i,j,k)/12.
              ave_npp(i,j) = ave_npp(i,j) + npp(i,j,k)/12.
              ave_lai(i,j) = ave_lai(i,j) + lai(i,j,k)/12.
              ave_clit(i,j) = ave_clit(i,j) + clit(i,j,k)/12.
              ave_cs(i,j) = ave_cs(i,j) + csoil(i,j,k)/12.
              ave_hr(i,j) = ave_hr(i,j) + hr(i,j,k)/12.
              ave_rc(i,j) = ave_rc(i,j) + rcm(i,j,k)/12.
              ave_runom(i,j) = ave_runom(i,j) + runom(i,j,k)/12.
              ave_evap(i,j) = ave_evap(i,j) + evaptr(i,j,k)/12.
              ave_wsoil(i,j) = ave_wsoil(i,j) + wsoil(i,j,k)/12.
           enddo
        else
           ave_ph(i,j) = no_data
           ave_ar(i,j) = no_data
           ave_npp(i,j) = no_data
           ave_lai(i,j) = no_data
           ave_clit(i,j) = no_data
           ave_cs(i,j) = no_data
           ave_hr(i,j) = no_data
           ave_rc(i,j) = no_data
           ave_runom(i,j) = no_data
           ave_evap(i,j) = no_data
           ave_wsoil(i,j) = no_data
        endif
     enddo
  enddo
  
  open(50,file='./outputs/ambientais.bin',&
       &        status='unknown',form='unformatted',&
       &        access='direct',recl=4*nx*ny)
  write(50,rec=1) ave_npp
  write(50,rec=2) ave_rc
  write(50,rec=3) ave_ar
  write(50,rec=4) ave_lai
  write(50,rec=5) ave_clit
  write(50,rec=6) ave_cs
  write(50,rec=7) ave_hr
  write(50,rec=8) ave_runom
  write(50,rec=9) ave_evap
  write(50,rec=10) ave_wsoil
  write(50,rec=11) ave_ph
  close(50)
  
  
!!$  open(10,file='../outputs/ph.bin',&
!!$       &     status='unknown',form='unformatted',&
!!$       &     access='direct',recl=4*nx*ny)
!!$  call save_file12(10, ph)
!!$  
!!$  open(10,file='../outputs/ar.bin',&
!!$       &     status='unknown',form='unformatted',&
!!$       &     access='direct',recl=4*nx*ny)
!!$  call save_file12(10, ar)
!!$  
!!$  open(10,file='../outputs/npp.bin',&
!!$       &     status='unknown',form='unformatted',&
!!$       &     access='direct',recl=4*nx*ny)
!!$  call save_file12(10, npp)
!!$  
!!$  open(10,file='../outputs/clit.bin',&
!!$       &     status='unknown',form='unformatted',&
!!$       &     access='direct',recl=4*nx*ny)
!!$  call save_file12(10, clit)
!!$  
!!$  open(10,file='../outputs/csoil.bin',&
!!$       &     status='unknown',form='unformatted',&
!!$       &     access='direct',recl=4*nx*ny)
!!$  call save_file12(10, csoil)
!!$  
!!$  open(10,file='../outputs/hr.bin',&
!!$       &     status='unknown',form='unformatted',&
!!$       &     access='direct',recl=4*nx*ny)
!!$  call save_file12(10, hr)
!!$  
!!$  open(10,file='../outputs/rcm.bin',&
!!$       &     status='unknown',form='unformatted',&
!!$       &     access='direct',recl=4*nx*ny)
!!$  call save_file12(10, rcm)
!!$  
!!$  open(10,file='../outputs/runom.bin',&
!!$       &     status='unknown',form='unformatted',&
!!$       &     access='direct',recl=4*nx*ny)
!!$  call save_file12(10, runom)
!!$  
!!$  open(10,file='../outputs/evaptr.bin',&
!!$       &     status='unknown',form='unformatted',&
!!$       &     access='direct',recl=4*nx*ny)
!!$  call save_file12(10, evaptr)
!!$  
!!$  open(10,file='../outputs/wsoil.bin',&
!!$       &     status='unknown',form='unformatted',&
!!$       &     access='direct',recl=4*nx*ny)
!!$  call save_file12(10, wsoil)
!!$  
!!$  open(10,file='../outputs/rm.bin',&
!!$       &     status='unknown',form='unformatted',&
!!$       &     access='direct',recl=4*nx*ny)
!!$  call save_file12(10, rm)
!!$  
!!$  open(10,file='../outputs/rg.bin',&
!!$       &     status='unknown',form='unformatted',&
!!$       &     access='direct',recl=4*nx*ny)
!!$  call save_file12(10, rg)

  stop
contains
  
end program env
    
