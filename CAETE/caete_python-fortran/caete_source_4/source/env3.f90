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
  !     
  real(kind=r4) :: ca                   !CO2 concentration (Pa)
  real(kind=r4) :: lsmk(nx,ny)          !Land=1/Ocean=0
  real(kind=r4) :: p0(nx,ny,nt)         !Atmospheric pressure (mb)
  real(kind=r4) :: ps(nx,ny,nt)         ! auxiliar to read atm pressure information
  real(kind=r4) :: prec(nx,ny,nt)       !Precipitation (mm/month)
  real(kind=r4) :: pr(nx,ny,nt)         !Auxiliar_precipitation (mm/month)
  real(kind=r4) :: temp(nx,ny,nt)       !Temperature (oC)
  real(kind=r4) :: t(nx,ny,nt)          !Auxiliar_temperature (oC)
  real(kind=r4) :: par(nx,ny,nt)        !Incident photosynthetic active radiation (Ein/m2/s)
  real(kind=r4) :: ipar(nx,ny,nt)       !Auxiliar_incident photosynthetic active radiation (w/m2)
  real(kind=r4) :: rhs(nx,ny,nt)        !Relative humidity
  real(kind=r4) :: rhaux(nx,ny,nt)      !RHS auxiliar
  real(kind=r4),dimension(nx,ny,q) :: cleafin = no_data ,cawoodin = no_data ,cfrootin = no_data

  !    Model Outputs
  !    -------------
  real(kind=r4) :: emaxm(nx,ny,nt)
  real(kind=r4) :: tsoil(nx,ny,nt)
  
  real(kind=r4) :: photo_pft(nx,ny,nt,q) !Monthly photosynthesis   (kgC/m2)
  real(kind=r4) :: aresp_pft(nx,ny,nt,q) !Monthly autotrophic res  (kgC/m2)
  real(kind=r4), dimension(nx,ny,nt,q) :: rm_pft,rg_pft
  real(kind=r4) :: npp_pft(nx,ny,nt,q)  !Monthly net primary produ (kgC/m2)
      
  real(kind=r4) :: lai_pft(nx,ny,nt,q)  !Monthly leaf area index
  real(kind=r4) :: clit_pft(nx,ny,nt,q) !Monthly litter carbon
  real(kind=r4) :: csoil_pft(nx,ny,nt,q) !Monthly soil carbon
  real(kind=r4) :: hresp_pft(nx,ny,nt,q) !Monthly het resp          (kgC/m2)
  real(kind=r4) :: rcm_pft(nx,ny,nt,q) 
  
  real(kind=r4) :: runom_pft(nx,ny,nt,q) !Runoff
  real(kind=r4) :: evapm_pft(nx,ny,nt,q) !Actual evapotranspiration        
  real(kind=r4) :: wsoil_pft(nx,ny,nt,q) !Soil moisture (mm)
  
  real(kind=r4),dimension(nx,ny,q) :: cleaf_pft,cawood_pft,cfroot_pft
  real(kind=r4) :: gridcell_ocp(nx,ny,q) !  final grid cell occupation for each pft (percentage of area)


  !     variaveis do spinup
  real(kind=r4) npp_sca 
  real(kind=r4) :: npp_pot(nx,ny,nt)
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

      
  real(kind=r4), dimension(nx,ny,nt) :: ph
  real(kind=r4), dimension(nx,ny,nt) :: ar
  real(kind=r4), dimension(nx,ny,nt) :: npp
  real(kind=r4), dimension(nx,ny,nt) :: lai
  real(kind=r4), dimension(nx,ny,nt) :: clit
  real(kind=r4), dimension(nx,ny,nt) :: csoil
  real(kind=r4), dimension(nx,ny,nt) :: hr
  real(kind=r4), dimension(nx,ny,nt) :: rcm
  real(kind=r4), dimension(nx,ny,nt) :: evaptr
  real(kind=r4), dimension(nx,ny,nt) :: wsoil
  real(kind=r4), dimension(nx,ny,nt) :: runom
  
  real(kind=r4), dimension(nx,ny,nt) :: rm
  real(kind=r4), dimension(nx,ny,nt) :: rg
            
  !     Open INPUT files
  !     ==========
  open( 9,file='../inputs/lsmk.bin',status='old',form='unformatted'&
       &,access='direct',recl=4*nx*ny)

  open(10,file='../inputs/ps.bin',status='old',form='unformatted'&
       &,access='direct',recl=4*nx*ny)

  open(11,file='../inputs/pr.bin',status='old',form='unformatted'&
       &,access='direct',recl=4*nx*ny)
      
  open(12,file='../inputs/tas.bin',status='old',form='unformatted'&
       &,access='direct',recl=4*nx*ny)
      
  open(13,file='../inputs/rsds.bin',status='old',form='unformatted'&
       &,access='direct',recl=4*nx*ny)

  open(14,file='../inputs/hurs.bin',status='old',form='unformatted'&
       &,access='direct',recl=4*nx*ny)

  open(26,file='../inputs/npp.bin',status='old',form='unformatted'&
       &,access='direct',recl=4*nx*ny)

!c      open(27,file='../spinup/clini.bin',status='old',
!c     &    form='unformatted',access='direct',recl=4*nx*ny)
!c
!c      open(28,file='../spinup/cfini.bin',status='old',
!c     &    form='unformatted',access='direct',recl=4*nx*ny)
!c
!c      open(29,file='../spinup/cwini.bin',status='old',
!c     &    form='unformatted',access='direct',recl=4*nx*ny)


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
  !c     call readx(29,cleafin,q)
  !c     call readx(28,cfrootin,q)
  !c     call readx(29,cawoodin,q)
      
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
  !c       close(27)
  !c       close(28)
  !c       close(29)
       

  !c     Calculating annual npp
  do i =1,nx
     do j=1,ny
        if(nint(lsmk(i,j)) .ne. 0) then 
           aux_npp(i,j) = 0.0
           do k = 1,12
              aux_npp(i,j) = aux_npp(i,j) + (npp_pot(i,j,k)/12.) 
           enddo
        else
           aux_npp(i,j) = no_data
        endif
     enddo
  enddo
       
  !     calling spinup
  print*, 'running spinup'
  do i=1,nx
     do j=1,ny
        if (nint(lsmk(i,j)) .ne. 0) then
           npp_sca = aux_npp(i,j)
           
           do p=1,q   
              aux1(p) = 0.0
              aux2(p) = 0.0
              aux3(p) = 0.0
              gridcell_ocp(i,j,p) = 0.0
           enddo
           
           call spinup(npp_sca, aux1, aux2, aux3)
           
           do p=1,q   
              cleafin(i,j,p)  = aux1(p)
              cfrootin(i,j,p) = aux2(p)
              cawoodin(i,j,p) = aux3(p)
           enddo
        else
           do p = 1,q
              cleafin(i,j,p)  = no_data
              cfrootin(i,j,p) = no_data
              cawoodin(i,j,p) = no_data
           enddo
        endif
     enddo
  enddo
      
  open(10,file='../spinup/clini.bin',&
       &    status='unknown',form='unformatted',&
       &    access='direct',recl=4*nx*ny)
  
  call savex(10, cleafin, q)
  
  open(10,file='../spinup/cfini.bin',&
       &    status='unknown',form='unformatted',&
       &    access='direct',recl=4*nx*ny)
  call savex(10, cfrootin, q)
  
  open(10,file='../spinup/cwini.bin',&
       &    status='unknown',form='unformatted',&
       &    access='direct',recl=4*nx*ny)
  call savex(10, cawoodin, q)
  
  do i=1,nx
     do j=1,ny
        
        !     this block set ocean grid cells to no_data 
        if(nint(lsmk(i,j)) .eq. 0) then
           do p = 1,q
              cleaf_pft(i,j,p)  = no_data
              cfroot_pft(i,j,p) = no_data
              cawood_pft(i,j,p) = no_data
              gridcell_ocp(i,j,p) = no_data
           enddo
        endif
        !     ------------------------------------------
        
        do k=1,12
           rhs(i,j,k) =  rhaux(i,j,k) / 100. !(rhaux(60,41,10) / 100.0) !Humidade relativa de manaus em outubro
           par(i,j,k) = ipar(i,j,k)/2.18E5 !Converting to Ein/m2/s
           temp(i,j,k) = t(i,j,k) !+ant(i,j,k) !uncomment to use future anomalies
           p0(i,j,k) = ps(i,j,k) * 0.01 ! transforamando de pascal pra mbar (kPa)
           prec(i,j,k) = pr(i,j,k) !+anpr(i,j,k) !+pr(i,j,k)*0.2 !uncomment to use future anomalies
           !c     if (prec(i,j,k).lt.0.0) prec (i,j,k) = 0.0  
        enddo
     enddo
  enddo
      
  !     Atmospheric CO2 pressure (Pa) !Ppmv / Pa
  ca= 363/9.901             !Pa (=363 ppmv; 1981-2010)
      
      
  !     =======================================
  !     Calculate environmental variables (wbm)
  !     =======================================
      
  call wbm (prec,temp,lsmk,p0,ca,par,rhs,cleafin,cawoodin,cfrootin,&
       &    emaxm, tsoil, photo_pft,aresp_pft,npp_pft,lai_pft,&
       &    clit_pft,csoil_pft, hresp_pft,rcm_pft,runom_pft,&
       &    evapm_pft,wsoil_pft,rml_pft,rmf_pft,rms_pft,rm_pft,rgl_pft&
       &    ,rgf_pft,rgs_pft,rg_pft,cleaf_pft,cawood_pft, cfroot_pft&
       &    ,gridcell_ocp,betal,betaw,betaf)   
      
      
      
      
  !     SAVE RESULTS TO FILES
  open(10,file='../outputs/gridcell_ocp.bin',&
       &    status='unknown',form='unformatted',&
       &    access='direct',recl=4*nx*ny)
  call savex(10, gridcell_ocp, q)
  
  open(10,file='../outputs/cleaf.bin',&
       &    status='unknown',form='unformatted',&
       &    access='direct',recl=4*nx*ny)
  call savex(10, cleaf_pft, q)
  
  open(10,file='../outputs/cawood.bin',&
       &    status='unknown',form='unformatted',&
       &    access='direct',recl=4*nx*ny)
  call savex(10,cawood_pft,q)
      
  open(10,file='../outputs/cfroot.bin',&
       &    status='unknown',form='unformatted',&
       &    access='direct',recl=4*nx*ny)
  call savex(10, cfroot_pft,q)
  
  do i = 1,nx
     do j = 1,ny
        do k = 1,12
           if(nint(lsmk(i,j)) .ne. 0) then
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
              rml(i,j,k)  = 0.0
              rmf(i,j,k)  = 0.0
              rms(i,j,k)  = 0.0
              rm(i,j,k)  = 0.0
              rgl(i,j,k)  = 0.0
              rgf(i,j,k)  = 0.0
              rgs(i,j,k)  = 0.0
              rg(i,j,k)  = 0.0
           else
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
              rml(i,j,k)  = no_data
              rmf(i,j,k)  = no_data
              rms(i,j,k)  = no_data
              rm(i,j,k)  = no_data
              rgl(i,j,k)  = no_data
              rgf(i,j,k)  = no_data
              rgs(i,j,k)  = no_data
              rg(i,j,k)  = no_data
           endif
        enddo
     enddo
  enddo
  do i = 1,nx
     do j = 1,ny
        if(nint(lsmk(i,j)) .ne. 0) then
           do k = 1,nt
              do p = 1,q
                 ph(i,j,k) = ph(i,j,k) + photo_pft(i,j,k,p)
                 ar(i,j,k) = ar(i,j,k) + aresp_pft(i,j,k,p)
                 npp(i,j,k) = npp(i,j,k) + npp_pft(i,j,k,p)
                 lai(i,j,k) = lai(i,j,k) + lai_pft(i,j,k,p)
                 clit(i,j,k) = clit(i,j,k) + clit_pft(i,j,k,p)
                 csoil(i,j,k) = csoil(i,j,k) + csoil_pft(i,j,k,p)
                 hr(i,j,k) = hr(i,j,k) + hresp_pft(i,j,k,p)
                 rcm(i,j,k) = rcm(i,j,k) + rcm_pft(i,j,k,p)
                 runom(i,j,k) = runom(i,j,k) + runom_pft(i,j,k,p)
                 evaptr(i,j,k) = evaptr(i,j,k) + evapm_pft(i,j,k,p)
                 wsoil(i,j,k) = wsoil(i,j,k) + wsoil_pft(i,j,k,p)
                 rml(i,j,k)  = rml(i,j,k) + rml_pft(i,j,k,p)
                 rmf(i,j,k)  = rmf(i,j,k) + rmf_pft(i,j,k,p)
                 rms(i,j,k)  = rms(i,j,k) + rms_pft(i,j,k,p)
                 rm(i,j,k)  = rm(i,j,k) + rm_pft(i,j,k,p)
                 rgl(i,j,k)  = rgl(i,j,k) + rgl_pft(i,j,k,p)
                 rgf(i,j,k)  = rgf(i,j,k) + rgf_pft(i,j,k,p)
                 rgs(i,j,k)  = rgs(i,j,k) + rgs_pft(i,j,k,p)
                 rg(i,j,k)  = rg(i,j,k) + rg_pft(i,j,k,p)
              enddo
           enddo
        endif
     enddo
  enddo

  open(10,file='../outputs/ph.bin',&
       &     status='unknown',form='unformatted',&
       &     access='direct',recl=4*nx*ny)
  call save_file12(10, ph)
  
  open(10,file='../outputs/ar.bin',&
       &     status='unknown',form='unformatted',&
       &     access='direct',recl=4*nx*ny)
  call save_file12(10, ar)
  
  open(10,file='../outputs/npp.bin',&
       &     status='unknown',form='unformatted',&
       &     access='direct',recl=4*nx*ny)
  call save_file12(10, npp)
  
  open(10,file='../outputs/clit.bin',&
       &     status='unknown',form='unformatted',&
       &     access='direct',recl=4*nx*ny)
  call save_file12(10, clit)
  
  open(10,file='../outputs/csoil.bin',&
       &     status='unknown',form='unformatted',&
       &     access='direct',recl=4*nx*ny)
  call save_file12(10, csoil)
  
  open(10,file='../outputs/hr.bin',&
       &     status='unknown',form='unformatted',&
       &     access='direct',recl=4*nx*ny)
  call save_file12(10, hr)
  
  open(10,file='../outputs/rcm.bin',&
       &     status='unknown',form='unformatted',&
       &     access='direct',recl=4*nx*ny)
  call save_file12(10, rcm)
  
  open(10,file='../outputs/runom.bin',&
       &     status='unknown',form='unformatted',&
       &     access='direct',recl=4*nx*ny)
  call save_file12(10, runom)
  
  open(10,file='../outputs/evaptr.bin',&
       &     status='unknown',form='unformatted',&
       &     access='direct',recl=4*nx*ny)
  call save_file12(10, evaptr)
  
  open(10,file='../outputs/wsoil.bin',&
       &     status='unknown',form='unformatted',&
       &     access='direct',recl=4*nx*ny)
  call save_file12(10, wsoil)
  
  open(10,file='../outputs/rml.bin',&
       &     status='unknown',form='unformatted',&
       &     access='direct',recl=4*nx*ny)
  call save_file12(10, rml)
  
  open(10,file='../outputs/rms.bin',&
     &     status='unknown',form='unformatted',&
     &     access='direct',recl=4*nx*ny)
  call save_file12(10, rms)
  
  open(10,file='../outputs/rmf.bin',&
       &     status='unknown',form='unformatted',&
       &     access='direct',recl=4*nx*ny)
  call save_file12(10, rmf)
  
  open(10,file='../outputs/rm.bin',&
       &     status='unknown',form='unformatted',&
       &     access='direct',recl=4*nx*ny)
  call save_file12(10, rm)
  
  open(10,file='../outputs/rgl.bin',&
       &     status='unknown',form='unformatted',&
       &     access='direct',recl=4*nx*ny)
  call save_file12(10, rgl)
  
  open(10,file='../outputs/rgf.bin',&
       &     status='unknown',form='unformatted',&
       &     access='direct',recl=4*nx*ny)
  call save_file12(10, rgf)
  
  open(10,file='../outputs/rgs.bin',&
       &     status='unknown',form='unformatted',&
       &     access='direct',recl=4*nx*ny)
  call save_file12(10, rgs)
  
  open(10,file='../outputs/rg.bin',&
       &     status='unknown',form='unformatted',&
       &     access='direct',recl=4*nx*ny)
  call save_file12(10, rg)
  
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
  
  open(50,file='../outputs/ambientais.bin',&
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
  
  stop
end program env
    
