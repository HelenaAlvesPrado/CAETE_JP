c234567
!     
      program env
!     
      implicit none
!     
!     =======================================================================
!     CPTEC-PVM2
!     Adapted from env4.f.
!     Code written by David Lapola and Helena Alves do Prado e Bianca Rius
c     Reviewed by jpdarela  jan/2017

!     Last update: jan/2017
!     
!     Compile with: g95 (or gfortran) env5TR.f wbm4TR.f productivity1.f 
!     Execute with ./a.exe
!     =======================================================================



!     Parameters and variables
!     ------------------------
!     
      integer i,j,k,p
      integer, parameter :: nx=720,ny=360,q=7
      real,parameter :: no_data = -9999.0
!     
!     Inputs
!     ------
!     
      real ca                   !CO2 concentration (Pa)
      real lsmk(nx,ny)          !Land=1/Ocean=0
      real p0(nx,ny,12)         !Atmospheric pressure (mb)
      real ps(nx,ny,12)         ! auxiliar to read atm pressure information
      real prec(nx,ny,12)       !Precipitation (mm/month)
      real pr(nx,ny,12)         !Auxiliar_precipitation (mm/month)
      real temp(nx,ny,12)       !Temperature (oC)
      real t(nx,ny,12)          !Auxiliar_temperature (oC)
      real par(nx,ny,12)        !Incident photosynthetic active radiation (Ein/m2/s)
      real ipar(nx,ny,12)       !Auxiliar_incident photosynthetic active radiation (w/m2)
      real rhs(nx,ny,12)        !Relative humidity
      real rhaux(nx,ny,12)      !RHS auxiliar
      real, dimension(nx,ny,q) :: cleafin = no_data
     $                          ,cawoodin = no_data
     $                          ,cfrootin = no_data


!     
!     Outputs
!     -------

c     Vou declarar aqui os outputas para a wbm, note que estas varia-
c     eis recebem os mesmos nomes das variaveis declaradas na definicao da wbm
c     porem, elas nao sao as mesmas variveis-- todas as variaveis da wbm sao
c     secretas para o env5.f.
      
      
c     -------------------------O U T P U T S--P A R A--W B M ------------
      
c     primeiro algumas variaveis que nao sao dependentes de pfts

      real emaxm(nx,ny,12)
      real tsoil(nx,ny,12)
      
      
c     agora as variaveis para pfts
      real photo_pft(nx,ny,12,q) !Monthly photosynthesis   (kgC/m2)
      real aresp_pft(nx,ny,12,q) !Monthly autotrophic res  (kgC/m2)
      real npp_pft(nx,ny,12,q)  !Monthly net primary produ (kgC/m2)
      
      real lai_pft(nx,ny,12,q)  !Monthly leaf area index
      real clit_pft(nx,ny,12,q) !Monthly litter carbon
      real csoil_pft(nx,ny,12,q) !Monthly soil carbon
      real hresp_pft(nx,ny,12,q) !Monthly het resp          (kgC/m2)
      real rcm_pft(nx,ny,12,q) 
      
c     VARIAVEIS HIDROLOGICAS IMPORTANTES   
      real runom_pft(nx,ny,12,q) !Runoff
      real evapm_pft(nx,ny,12,q) !Actual evapotranspiration        
      real wsoil_pft(nx,ny,12,q) !Soil moisture (mm)
      
c     variables related to carbon allocation and autothrophic respiration (bianca)
      real, dimension(nx,ny,12,q) :: rml_pft,rmf_pft,rms_pft,rm_pft,
     $rgl_pft,rgf_pft,rgs_pft,rg_pft
      real, dimension(nx,ny,q) :: cleaf_pft,cawood_pft,cfroot_pft


C      WARNING - NEW VARIABLES ---
      
      real gridcell_ocp(nx,ny,q) !  final grid cell occupation for each pft (percentage of area)
      real betal(nx,ny,12,q)
      real, dimension(nx,ny,12):: bl1,bl2,bl3,bl4,bl5,bl6,bl7 ! carbon allocated to growth (monthly sums) 
      real betaw(nx,ny,12,q)
      real, dimension(nx,ny,12):: bw1,bw2,bw3,bw4,bw5,bw6,bw7
      real betaf(nx,ny,12,q)
      real, dimension(nx,ny,12):: bf1,bf2,bf3,bf4,bf5,bf6,bf7

      
c     variaveis do spinup
      real, dimension(q) :: aux1, aux2, aux3
      real npp_pot(nx,ny,12)
      real, dimension(nx,ny) :: aux_npp
      real npp_sca  

c     -----FIM DA DEFINICAO DE VARIAVEIS PARA RODAR O MODELO--
C
c     variaveis anuais -
      real, dimension(nx,ny) :: ave_ph = 0.0
      real, dimension(nx,ny) :: ave_ar = 0.0
      real, dimension(nx,ny) :: ave_npp = 0.0
      real, dimension(nx,ny) :: ave_lai = 0.0
      real, dimension(nx,ny) :: ave_clit = 0.0
      real, dimension(nx,ny) :: ave_cs = 0.0
      real, dimension(nx,ny) :: ave_hr = 0.0
      real, dimension(nx,ny) :: ave_rc = 0.0
      real, dimension(nx,ny) :: ave_runom = 0.0
      real, dimension(nx,ny) :: ave_evap = 0.0
      real, dimension(nx,ny) :: ave_wsoil = 0.0    
      
C     THESE WILL RECEIVE MEANS BETWEEN q PFTs and for each pft (ex. ph to mean; ph1 to pft 1)
      real, dimension(nx,ny,12) :: ph!, ph1, ph2, ph3  
c, ph4, ph5, ph6, ph7
      real, dimension(nx,ny,12) :: ar!, ar1, ar2, ar3, ar4, ar5, ar6, ar7
      real, dimension(nx,ny,12) :: npp!,npp1,npp2,npp3
c,npp4,npp5,npp6,npp7
      real, dimension(nx,ny,12) :: lai!, lai1, lai2, lai3
c lai4, lai5,lai6, lai7
      real, dimension(nx,ny,12) :: clit!, clit1, clit2, clit3
c clit4,   &    clit5, clit6, clit7
      real, dimension(nx,ny,12) :: csoil!, csoil1, csoil2, csoil3
c, csoil4    &    , csoil5, csoil6, csoil7
      real, dimension(nx,ny,12) :: hr!, hr1, hr2, hr3
c, hr4, hr5, hr6, hr7
      real, dimension(nx,ny,12) :: rcm!, rcm1, rcm2, rcm3
c, rcm4, rcm5,    &    rcm6, rcm7
      real, dimension(nx,ny,12) :: evaptr!, et1, et2, et3
c, et4, et5, et6,    &    et7
      real, dimension(nx,ny,12) :: wsoil
      real, dimension(nx,ny,12) :: runom!
      
C     NEW OUTPUTS (AUTOTRF RESPIRATION, ALLOCATION)
      
      real, dimension(nx,ny,12) :: rml!, rml1, rml2, rml3
      real, dimension(nx,ny,12) :: rmf!, rmf1, rmf2, rmf3 
      real, dimension(nx,ny,12) :: rms!, rms1, rms2, rms3
      real, dimension(nx,ny,12) :: rm!,  rm1,  rm2,  rm3
      real, dimension(nx,ny,12) :: rgl!, rgl1, rgl2, rgl3
      real, dimension(nx,ny,12) :: rgf!, rgf1, rgf2, rgf3
      real, dimension(nx,ny,12) :: rgs!, rgs1, rgs2, rgs3
      real, dimension(nx,ny,12) :: rg!,  rg1,  rg2,  rg3 
      
C     -------END DECLARATION----------------------------------------
      
!     Open INPUT files
!     ==========
      open( 9,file='../inputs/lsmk.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)

      open(10,file='../inputs/ps.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)

      open(11,file='../inputs/pr.bin',status='old',
     &    form='unformatted',access='direct',recl=4*nx*ny)
      
      open(12,file='../inputs/tas.bin',status='old',
     &    form='unformatted',access='direct',recl=4*nx*ny)
      
      open(13,file='../inputs/rsds.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)

      open(14,file='../inputs/hurs.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)

      open(26,file='../inputs/npp.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)

c      open(27,file='../spinup/clini.bin',status='old',
c     &    form='unformatted',access='direct',recl=4*nx*ny)
c
c      open(28,file='../spinup/cfini.bin',status='old',
c     &    form='unformatted',access='direct',recl=4*nx*ny)
c
c      open(29,file='../spinup/cwini.bin',status='old',
c     &    form='unformatted',access='direct',recl=4*nx*ny)
!     Read data
!     =========
      
      read (9,rec=1) lsmk
c     read (26,rec=1) aux_npp 
      
      call readx(10,ps,12)
      call readx(11,pr,12)
      call readx(12,t,12)
      call readx(13,ipar,12)
      call readx(14,rhaux,12)
      call readx(26,npp_pot,12)
c     call readx(29,cleafin,q)
c     call readx(28,cfrootin,q)
c     call readx(29,cawoodin,q)
      
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
c       close(27)
c       close(28)
c       close(29)
       

c     Calculating annual npp
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
!     if(mod(nx,10) .eq. 0)print*, (real(i)/real(nx))*100.0, '%'
      enddo
      
      call nan2ndt(cleafin, q)  !!! --------- incorporado essa subroutina
      open(10,file='../spinup/clini.bin',
     &    status='unknown',form='unformatted',
     &    access='direct',recl=4*nx*ny)
      call savex(10, cleafin, q)
      
      call nan2ndt(cfrootin, q)
      open(10,file='../spinup/cfini.bin',
     &    status='unknown',form='unformatted',
     &    access='direct',recl=4*nx*ny)
      call savex(10, cfrootin, q)
      
      call nan2ndt(cawoodin, q)
      open(10,file='../spinup/cwini.bin',
     &    status='unknown',form='unformatted',
     &    access='direct',recl=4*nx*ny)
      call savex(10, cawoodin, q)
!     ===========
      
      
      
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
c     if (prec(i,j,k).lt.0.0) prec (i,j,k) = 0.0  
            enddo
         enddo
      enddo
      
!     Atmospheric CO2 pressure (Pa) !Ppmv / Pa
      ca= 363/9.901             !Pa (=363 ppmv; 1981-2010)
      
      
!     =======================================
!     Calculate environmental variables (wbm)
!     =======================================
      
      call wbm (prec,temp,lsmk,p0,ca,par,rhs,cleafin,cawoodin,cfrootin,
     &    emaxm, tsoil, photo_pft,aresp_pft,npp_pft,lai_pft,
     &    clit_pft,csoil_pft, hresp_pft,rcm_pft,runom_pft,
     &    evapm_pft,wsoil_pft,rml_pft,rmf_pft,rms_pft,rm_pft,rgl_pft
     &    ,rgf_pft,rgs_pft,rg_pft,cleaf_pft,cawood_pft, cfroot_pft
     &    ,gridcell_ocp,betal,betaw,betaf)   
      
      
      
      
!     SAVE RESULTS TO FILES
      call nan2ndt(gridcell_ocp, q)
      open(10,file='../outputs/gridcell_ocp.bin',
     &    status='unknown',form='unformatted',
     &    access='direct',recl=4*nx*ny)
      call savex(10, gridcell_ocp, q)
      
      call nan2ndt(cleaf_pft, q)
      open(10,file='../outputs/cleaf.bin',
     &    status='unknown',form='unformatted',
     &    access='direct',recl=4*nx*ny)
      call savex(10, cleaf_pft, q)
      
      call nan2ndt(cawood_pft, q)
      open(10,file='../outputs/cawood.bin',
     &    status='unknown',form='unformatted',
     &    access='direct',recl=4*nx*ny)
      call savex(10,cawood_pft,q)
      
      call nan2ndt(cfroot_pft, q)
      open(10,file='../outputs/cfroot.bin',
     &    status='unknown',form='unformatted',
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
      
           
C     preparando o terreno pra salvar as variaveis
      
      do i = 1,nx
         do j = 1,ny
            if(nint(lsmk(i,j)) .ne. 0) then
               do k = 1,12
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

      open(10,file='../outputs/ph.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, ph)
      
      open(10,file='../outputs/ar.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, ar)

      open(10,file='../outputs/npp.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, npp)

      open(10,file='../outputs/clit.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, clit)

      open(10,file='../outputs/csoil.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, csoil)

      open(10,file='../outputs/hr.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, hr)

      open(10,file='../outputs/rcm.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, rcm)

      open(10,file='../outputs/runom.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, runom)

      open(10,file='../outputs/evaptr.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, evaptr)

      open(10,file='../outputs/wsoil.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, wsoil)

      open(10,file='../outputs/rml.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, rml)

      open(10,file='../outputs/rms.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, rms)

      open(10,file='../outputs/rmf.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, rmf)

      open(10,file='../outputs/rm.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, rm)

      open(10,file='../outputs/rgl.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, rgl)

      open(10,file='../outputs/rgf.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, rgf)

      open(10,file='../outputs/rgs.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, rgs)

      open(10,file='../outputs/rg.bin',
     &     status='unknown',form='unformatted',
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

      open(50,file='../outputs/ambientais.bin',
     &        status='unknown',form='unformatted',
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

      
      do i=1,nx
         do j=1,ny
            do k=1,12
               if(nint(lsmk(i,j)) .ne. 0) then
c                  npp1(i,j,k) = npp_pft(i,j,k,1)
c                  npp2(i,j,k) = npp_pft(i,j,k,2)
c                  npp3(i,j,k) = npp_pft(i,j,k,3)
cc                  npp4(i,j,k) = npp_pft(i,j,k,4)
cc                  npp5(i,j,k) = npp_pft(i,j,k,5)
cc                  npp6(i,j,k) = npp_pft(i,j,k,6)
cc                  npp7(i,j,k) = npp_pft(i,j,k,7)
c
c                  ph1(i,j,k) = photo_pft(i,j,k,1)
c                  ph2(i,j,k) = photo_pft(i,j,k,2)
c                  ph3(i,j,k) = photo_pft(i,j,k,3)
cc                  ph4(i,j,k) = photo_pft(i,j,k,4)
cc                  ph5(i,j,k) = photo_pft(i,j,k,5)
cc                  ph6(i,j,k) = photo_pft(i,j,k,6)
cc                  ph7(i,j,k) = photo_pft(i,j,k,7)
c
c                  ar1(i,j,k) = aresp_pft(i,j,k,1)
c                  ar2(i,j,k) = aresp_pft(i,j,k,2)
c                  ar3(i,j,k) = aresp_pft(i,j,k,3)
cc                  ar4(i,j,k) = aresp_pft(i,j,k,4)
cc                  ar5(i,j,k) = aresp_pft(i,j,k,5)
cc                  ar6(i,j,k) = aresp_pft(i,j,k,6)
cc                  ar7(i,j,k) = aresp_pft(i,j,k,7)
c
c                  lai1(i,j,k) = lai_pft(i,j,k,1)
c                  lai2(i,j,k) = lai_pft(i,j,k,2)
c                  lai3(i,j,k) = lai_pft(i,j,k,3)
cc                  lai4(i,j,k) = lai_pft(i,j,k,4)
cc                  lai5(i,j,k) = lai_pft(i,j,k,5)
cc                  lai6(i,j,k) = lai_pft(i,j,k,6)
cc                  lai7(i,j,k) = lai_pft(i,j,k,7)
c
c                  hr1(i,j,k) = hresp_pft(i,j,k,1)
c                  hr2(i,j,k) = hresp_pft(i,j,k,2)
c                  hr3(i,j,k) = hresp_pft(i,j,k,3)
cc                  hr4(i,j,k) = hresp_pft(i,j,k,4)
cc                  hr5(i,j,k) = hresp_pft(i,j,k,5)
cc                  hr6(i,j,k) = hresp_pft(i,j,k,6)
cc                  hr7(i,j,k) = hresp_pft(i,j,k,7)
c
c                  clit1(i,j,k) = clit_pft(i,j,k,1)
c                  clit2(i,j,k) = clit_pft(i,j,k,2)
c                  clit3(i,j,k) = clit_pft(i,j,k,3)
cc                  clit4(i,j,k) = clit_pft(i,j,k,4)
cc                  clit5(i,j,k) = clit_pft(i,j,k,5)
cc                  clit6(i,j,k) = clit_pft(i,j,k,6)
cc                  clit7(i,j,k) = clit_pft(i,j,k,7)
c
c                  csoil1(i,j,k) = csoil_pft(i,j,k,1)
c                  csoil2(i,j,k) = csoil_pft(i,j,k,2)
c                  csoil3(i,j,k) = csoil_pft(i,j,k,3)
cc                  csoil4(i,j,k) = csoil_pft(i,j,k,4)
cc                  csoil5(i,j,k) = csoil_pft(i,j,k,5)
cc                  csoil6(i,j,k) = csoil_pft(i,j,k,6)
cc                  csoil7(i,j,k) = csoil_pft(i,j,k,7)
c
c                  et1(i,j,k) = evapm_pft(i,j,k,1)
c                  et2(i,j,k) = evapm_pft(i,j,k,2)
c                  et3(i,j,k) = evapm_pft(i,j,k,3)
cc                  et4(i,j,k) = evapm_pft(i,j,k,4)
cc                  et5(i,j,k) = evapm_pft(i,j,k,5)
cc                  et6(i,j,k) = evapm_pft(i,j,k,6)
cc                  et7(i,j,k) = evapm_pft(i,j,k,7)
c                  
c                  rcm1(i,j,k) = rcm_pft(i,j,k,1)
c                  rcm2(i,j,k) = rcm_pft(i,j,k,2)
c                  rcm3(i,j,k) = rcm_pft(i,j,k,3)
cc                  rcm4(i,j,k) = rcm_pft(i,j,k,4)
cc                  rcm5(i,j,k) = rcm_pft(i,j,k,5)
cc                  rcm6(i,j,k) = rcm_pft(i,j,k,6)
cc                  rcm7(i,j,k) = rcm_pft(i,j,k,7)

                  bl1(i,j,k) = betal(i,j,k,1)
                  bl2(i,j,k) = betal(i,j,k,2)
                  bl3(i,j,k) = betal(i,j,k,3)
                  bl4(i,j,k) = betal(i,j,k,4)
                  bl5(i,j,k) = betal(i,j,k,5)
                  bl6(i,j,k) = betal(i,j,k,6)
                  bl7(i,j,k) = betal(i,j,k,7)

                  bw1(i,j,k) = betaw(i,j,k,1)
                  bw2(i,j,k) = betaw(i,j,k,2)
                  bw3(i,j,k) = betaw(i,j,k,3)
                  bw4(i,j,k) = betaw(i,j,k,4)
                  bw5(i,j,k) = betaw(i,j,k,5)
                  bw6(i,j,k) = betaw(i,j,k,6)
                  bw7(i,j,k) = betaw(i,j,k,7)
                  
                  bf1(i,j,k) = betaf(i,j,k,1)
                  bf2(i,j,k) = betaf(i,j,k,2)
                  bf3(i,j,k) = betaf(i,j,k,3)
                  bf4(i,j,k) = betaf(i,j,k,4)
                  bf5(i,j,k) = betaf(i,j,k,5)
                  bf6(i,j,k) = betaf(i,j,k,6)
                  bf7(i,j,k) = betaf(i,j,k,7)
                  
               else
c                  npp1(i,j,k) = no_data
c                  npp2(i,j,k) = no_data
c                  npp3(i,j,k) = no_data
cc                  npp4(i,j,k) = no_data
cc                  npp5(i,j,k) = no_data
cc                  npp6(i,j,k) = no_data
cc                  npp7(i,j,k) = no_data
c
c                  ph1(i,j,k) = no_data
c                  ph2(i,j,k) = no_data
c                  ph3(i,j,k) = no_data
cc                  ph4(i,j,k) = no_data
cc                  ph5(i,j,k) = no_data
cc                  ph6(i,j,k) = no_data
cc                  ph7(i,j,k) = no_data
c
c                  ar1(i,j,k) = no_data
c                  ar2(i,j,k) = no_data
c                  ar3(i,j,k) = no_data
cc                  ar4(i,j,k) = no_data
cc                  ar5(i,j,k) = no_data
cc                  ar6(i,j,k) = no_data
cc                  ar7(i,j,k) = no_data
c
c                  hr1(i,j,k) = no_data
c                  hr2(i,j,k) = no_data
c                  hr3(i,j,k) = no_data
cc                  hr4(i,j,k) = no_data
cc                  hr5(i,j,k) = no_data
cc                  hr6(i,j,k) = no_data
cc                  hr7(i,j,k) = no_data
c
c                  clit1(i,j,k) = no_data
c                  clit2(i,j,k) = no_data
c                  clit3(i,j,k) = no_data
cc                  clit4(i,j,k) = no_data
cc                  clit5(i,j,k) = no_data
cc                  clit6(i,j,k) = no_data
cc                  clit7(i,j,k) = no_data
c
c                  csoil1(i,j,k) = no_data
c                  csoil2(i,j,k) = no_data
c                  csoil3(i,j,k) = no_data
cc                  csoil4(i,j,k) = no_data
cc                  csoil5(i,j,k) = no_data
cc                  csoil6(i,j,k) = no_data
cc                  csoil7(i,j,k) = no_data
c
c                  et1(i,j,k) = no_data
c                  et2(i,j,k) = no_data
c                  et3(i,j,k) = no_data
cc                  et4(i,j,k) = no_data
cc                  et5(i,j,k) = no_data
cc                  et6(i,j,k) = no_data
cc                  et7(i,j,k) = no_data
c                  
c                  rcm1(i,j,k) = no_data
c                  rcm2(i,j,k) = no_data
c                  rcm3(i,j,k) = no_data
cc                  rcm4(i,j,k) = no_data
cc                  rcm5(i,j,k) = no_data
cc                  rcm6(i,j,k) = no_data
cc                  rcm7(i,j,k) = no_data
                  
                  bl1(i,j,k) = no_data
                  bl2(i,j,k) = no_data
                  bl3(i,j,k) = no_data
                  bl4(i,j,k) = no_data
                  bl5(i,j,k) = no_data
                  bl6(i,j,k) = no_data
                  bl7(i,j,k) = no_data

                  bw1(i,j,k) = no_data
                  bw2(i,j,k) = no_data
                  bw3(i,j,k) = no_data
                  bw4(i,j,k) = no_data
                  bw5(i,j,k) = no_data
                  bw6(i,j,k) = no_data
                  bw7(i,j,k) = no_data
                  
                  bf1(i,j,k) = no_data
                  bf2(i,j,k) = no_data
                  bf3(i,j,k) = no_data
                  bf4(i,j,k) = no_data
                  bf5(i,j,k) = no_data
                  bf6(i,j,k) = no_data
                  bf7(i,j,k) = no_data
               endif
            enddo
         enddo
      enddo

      
c!     NPP
c      open(10,file='../outputs_pft/npp.1.bin',
c     &    status='unknown',form='unformatted',
c     &    access='direct',recl=4*nx*ny)
c      call save_file12(10, npp1)
c      open(10,file='../outputs_pft/npp.2.bin',
c     &    status='unknown',form='unformatted',
c     &    access='direct',recl=4*nx*ny)
c      call save_file12(10, npp2)
c      open(10,file='../outputs_pft/npp.3.bin',
c     &    status='unknown',form='unformatted',
c     &    access='direct',recl=4*nx*ny)
c      call save_file12(10, npp3)
cc      open(10,file='../outputs_pft/npp.4.bin',
cc     &    status='unknown',form='unformatted',
cc     &    access='direct',recl=4*nx*ny)
cc      call save_file12(10, npp4)
cc      open(10,file='../outputs_pft/npp.5.bin',
cc     &     status='unknown',form='unformatted',
cc     &    access='direct',recl=4*nx*ny)
cc      call save_file12(10, npp5)
cc      open(10,file='../outputs_pft/npp.6.bin',
cc     &    status='unknown',form='unformatted',
cc     &    access='direct',recl=4*nx*ny)
cc      call save_file12(10, npp6)
cc      open(10,file='../outputs_pft/npp.7.bin',
cc     &    status='unknown',form='unformatted',
cc     &     access='direct',recl=4*nx*ny)
cc      call save_file12(10, npp7)
c
c!     PHOTO      
c      open(10,file='../outputs_pft/ph.1.bin',
c     &    status='unknown',form='unformatted',
c     &    access='direct',recl=4*nx*ny)
c      call save_file12(10, ph1)
c      open(10,file='../outputs_pft/ph.2.bin',
c     &    status='unknown',form='unformatted',
c     &    access='direct',recl=4*nx*ny)
c      call save_file12(10, ph2)
c      open(10,file='../outputs_pft/ph.3.bin',
c     &    status='unknown',form='unformatted',
c     &    access='direct',recl=4*nx*ny)
c      call save_file12(10, ph3)
cc      open(10,file='../outputs_pft/ph.4.bin',
cc     &    status='unknown',form='unformatted',
cc     &    access='direct',recl=4*nx*ny)
cc      call save_file12(10, ph4)
cc      open(10,file='../outputs_pft/ph.5.bin',
cc     &    status='unknown',form='unformatted',
cc     &    access='direct',recl=4*nx*ny)
cc      call save_file12(10, ph5)
cc      open(10,file='../outputs_pft/ph.6.bin',
cc     &    status='unknown',form='unformatted',
cc     &    access='direct',recl=4*nx*ny)
cc      call save_file12(10, ph6)
cc      open(10,file='../outputs_pft/ph.7.bin',
cc     &    status='unknown',form='unformatted',
cc     &    access='direct',recl=4*nx*ny)
cc      call save_file12(10, ph7)
c
c      
c!     ARESP
c      open(10,file='../outputs_pft/ar.1.bin',
c     &    status='unknown',form='unformatted',
c     &    access='direct',recl=4*nx*ny)
c      call save_file12(10, ar1)
c      open(10,file='../outputs_pft/ar.2.bin',
c     &    status='unknown',form='unformatted',
c     &    access='direct',recl=4*nx*ny)
c      call save_file12(10, ar2)
c      open(10,file='../outputs_pft/ar.3.bin',
c     &    status='unknown',form='unformatted',
c     &    access='direct',recl=4*nx*ny)
c      call save_file12(10, ar3)
cc      open(10,file='../outputs_pft/ar.4.bin',
cc     &    status='unknown',form='unformatted',
cc     &    access='direct',recl=4*nx*ny)
cc      call save_file12(10, ar4)
cc      open(10,file='../outputs_pft/ar.5.bin',
cc     &    status='unknown',form='unformatted',
cc     &    access='direct',recl=4*nx*ny)
cc      call save_file12(10, ar5)
cc      open(10,file='../outputs_pft/ar.6.bin',
cc     &    status='unknown',form='unformatted',
cc     &    access='direct',recl=4*nx*ny)
cc      call save_file12(10, ar6)
cc      open(10,file='../outputs_pft/ar.7.bin',
cc     &    status='unknown',form='unformatted',
cc     &    access='direct',recl=4*nx*ny)
cc      call save_file12(10, ar7)
c
c      
c!     HRESP
c      open(10,file='../outputs_pft/hr.1.bin',
c     &    status='unknown',form='unformatted',
c     &    access='direct',recl=4*nx*ny)
c      call save_file12(10, hr1)
c      open(10,file='../outputs_pft/hr.2.bin',
c     &    status='unknown',form='unformatted',
c     &    access='direct',recl=4*nx*ny)
c      call save_file12(10, hr2)
c      open(10,file='../outputs_pft/hr.3.bin',
c     &    status='unknown',form='unformatted',
c     &    access='direct',recl=4*nx*ny)
c      call save_file12(10, hr3)
cc      open(10,file='../outputs_pft/hr.4.bin',
cc     &    status='unknown',form='unformatted',
cc     &    access='direct',recl=4*nx*ny)
cc      call save_file12(10, hr4)
cc      open(10,file='../outputs_pft/hr.5.bin',
cc     &    status='unknown',form='unformatted',
cc     &    access='direct',recl=4*nx*ny)
cc      call save_file12(10, hr5)
cc      open(10,file='../outputs_pft/hr.6.bin',
cc     &    status='unknown',form='unformatted',
cc     &    access='direct',recl=4*nx*ny)
cc      call save_file12(10, hr6)
cc      open(10,file='../outputs_pft/hr.7.bin',
cc     &    status='unknown',form='unformatted',
cc     &    access='direct',recl=4*nx*ny)
cc      call save_file12(10, hr7)
c
c!     CLIT
c      open(10,file='../outputs_pft/clit.1.bin',
c     &     status='unknown',form='unformatted',
c     &     access='direct',recl=4*nx*ny)
c      call save_file12(10, clit1)
c      open(10,file='../outputs_pft/clit.2.bin',
c     &     status='unknown',form='unformatted',
c     &     access='direct',recl=4*nx*ny)
c      call save_file12(10, clit2)
c      open(10,file='../outputs_pft/clit.3.bin',
c     &     status='unknown',form='unformatted',
c     &     access='direct',recl=4*nx*ny)
c      call save_file12(10, clit3)
cc      open(10,file='../outputs_pft/clit.4.bin',
cc     &     status='unknown',form='unformatted',
cc     &     access='direct',recl=4*nx*ny)
cc      call save_file12(10, clit4)
cc      open(10,file='../outputs_pft/clit.5.bin',
cc     &     status='unknown',form='unformatted',
cc     &     access='direct',recl=4*nx*ny)
cc      call save_file12(10, clit5)
cc      open(10,file='../outputs_pft/clit.6.bin',
cc     &     status='unknown',form='unformatted',
cc     &     access='direct',recl=4*nx*ny)
cc      call save_file12(10, clit6)
cc      open(10,file='../outputs_pft/clit.7.bin',
cc     &     status='unknown',form='unformatted',
cc     &     access='direct',recl=4*nx*ny)
cc      call save_file12(10, clit7)
c
c!     CSOIL
c      open(10,file='../outputs_pft/csoil.1.bin',
c     &     status='unknown',form='unformatted',
c     &     access='direct',recl=4*nx*ny)
c      call save_file12(10, csoil1)
c      open(10,file='../outputs_pft/csoil.2.bin',
c     &     status='unknown',form='unformatted',
c     &     access='direct',recl=4*nx*ny)
c      call save_file12(10, csoil2)
c      open(10,file='../outputs_pft/csoil.3.bin',
c     &     status='unknown',form='unformatted',
c     &     access='direct',recl=4*nx*ny)
c      call save_file12(10, csoil3)
cc      open(10,file='../outputs_pft/csoil.4.bin',
cc     &     status='unknown',form='unformatted',
cc     &     access='direct',recl=4*nx*ny)
cc      call save_file12(10, csoil4)
cc      open(10,file='../outputs_pft/csoil.5.bin',
cc     &     status='unknown',form='unformatted',
cc     &     access='direct',recl=4*nx*ny)
cc      call save_file12(10, csoil5)
cc      open(10,file='../outputs_pft/csoil.6.bin',
cc     &     status='unknown',form='unformatted',
cc     &     access='direct',recl=4*nx*ny)
cc      call save_file12(10, csoil6)
cc      open(10,file='../outputs_pft/csoil.7.bin',
cc     &     status='unknown',form='unformatted',
cc     &     access='direct',recl=4*nx*ny)
cc      call save_file12(10, csoil7)
c
c!     EVAPM
c      open(10,file='../outputs_pft/et.1.bin',
c     &    status='unknown',form='unformatted',
c     &    access='direct',recl=4*nx*ny)
c      call save_file12(10, et1)
c      open(10,file='../outputs_pft/et.2.bin',
c     &    status='unknown',form='unformatted',
c     &    access='direct',recl=4*nx*ny)
c      call save_file12(10, et2)
c      open(10,file='../outputs_pft/et.3.bin',
c     &    status='unknown',form='unformatted',
c     &    access='direct',recl=4*nx*ny)
c      call save_file12(10, et3)
cc      open(10,file='../outputs_pft/et.4.bin',
cc     &    status='unknown',form='unformatted',
cc     &    access='direct',recl=4*nx*ny)
cc      call save_file12(10, et4)
cc      open(10,file='../outputs_pft/et.5.bin',
cc     &     status='unknown',form='unformatted',
cc     &    access='direct',recl=4*nx*ny)
cc      call save_file12(10, et5)
cc      open(10,file='../outputs_pft/et.6.bin',
cc     &    status='unknown',form='unformatted',
cc     &     access='direct',recl=4*nx*ny)
cc      call save_file12(10, et6)
cc      open(10,file='../outputs_pft/et.7.bin',
cc     &    status='unknown',form='unformatted',
cc     &    access='direct',recl=4*nx*ny)
cc      call save_file12(10, et7)
c      
c!     RCM 
c      open(10,file='../outputs_pft/rcm.1.bin',
c     &    status='unknown',form='unformatted',
c     &    access='direct',recl=4*nx*ny)
c      call save_file12(10, rcm1)  
c      open(10,file='../outputs_pft/rcm.2.bin',
c     &    status='unknown',form='unformatted',
c     &    access='direct',recl=4*nx*ny)
c      call save_file12(10, rcm2)
c      open(10,file='../outputs_pft/rcm.3.bin',
c     &    status='unknown',form='unformatted',
c     &    access='direct',recl=4*nx*ny)
c      call save_file12(10, rcm3)
cc      open(10,file='../outputs_pft/rcm.4.bin',
cc     &    status='unknown',form='unformatted',
cc     &    access='direct',recl=4*nx*ny)
cc      call save_file12(10, rcm4)
cc      open(10,file='../outputs_pft/rcm.5.bin',
cc     &    status='unknown',form='unformatted',
cc     &    access='direct',recl=4*nx*ny)
cc      call save_file12(10, rcm5)
cc      open(10,file='../outputs_pft/rcm.6.bin',
cc     &    status='unknown',form='unformatted',
cc     &    access='direct',recl=4*nx*ny)
cc      call save_file12(10, rcm6)
cc      open(10,file='../outputs_pft/rcm.7.bin',
cc     &    status='unknown',form='unformatted',
cc     &    access='direct',recl=4*nx*ny)
cc      call save_file12(10, rcm7)
      
!     BLEAF
      open(10,file='../outputs_pft/bl.1.bin',
     &    status='unknown',form='unformatted',
     &    access='direct',recl=4*nx*ny)
      call save_file12(10, bl1)
      open(10,file='../outputs_pft/bl.2.bin',
     &    status='unknown',form='unformatted',
     &    access='direct',recl=4*nx*ny)
      call save_file12(10, bl2)
      open(10,file='../outputs_pft/bl.3.bin',
     &    status='unknown',form='unformatted',
     &    access='direct',recl=4*nx*ny)
      call save_file12(10, bl3)
      open(10,file='../outputs_pft/bl.4.bin',
     &    status='unknown',form='unformatted',
     &    access='direct',recl=4*nx*ny)
      call save_file12(10, bl4)
      open(10,file='../outputs_pft/bl.5.bin',
     &    status='unknown',form='unformatted',
     &    access='direct',recl=4*nx*ny)
      call save_file12(10, bl5)
      open(10,file='../outputs_pft/bl.6.bin',
     &    status='unknown',form='unformatted',
     &    access='direct',recl=4*nx*ny)
      call save_file12(10, bl6)
      open(10,file='../outputs_pft/bl.7.bin',
     &    status='unknown',form='unformatted',
     &    access='direct',recl=4*nx*ny)
      call save_file12(10, bl7)

!     BAWOOD
      open(10,file='../outputs_pft/bw.1.bin',
     &    status='unknown',form='unformatted',
     &    access='direct',recl=4*nx*ny)
      call save_file12(10, bw1)
      open(10,file='../outputs_pft/bw.2.bin',
     &    status='unknown',form='unformatted',
     &    access='direct',recl=4*nx*ny)
      call save_file12(10, bw2)
      open(10,file='../outputs_pft/bw.3.bin',
     &    status='unknown',form='unformatted',
     &    access='direct',recl=4*nx*ny)
      call save_file12(10, bw3)
      open(10,file='../outputs_pft/bw.4.bin',
     &    status='unknown',form='unformatted',
     &    access='direct',recl=4*nx*ny)
      call save_file12(10, bw4)
      open(10,file='../outputs_pft/bw.5.bin',
     &    status='unknown',form='unformatted',
     &    access='direct',recl=4*nx*ny)
      call save_file12(10, bw5)
      open(10,file='../outputs_pft/bw.6.bin',
     &    status='unknown',form='unformatted',
     &    access='direct',recl=4*nx*ny)
      call save_file12(10, bw6)
      open(10,file='../outputs_pft/bw.7.bin',
     &    status='unknown',form='unformatted',
     &    access='direct',recl=4*nx*ny)
      call save_file12(10, bw7)

!     BFROOT
      open(10,file='../outputs_pft/bf.1.bin',
     &    status='unknown',form='unformatted',
     &    access='direct',recl=4*nx*ny)
      call save_file12(10, bl1)
      open(10,file='../outputs_pft/bf.2.bin',
     &    status='unknown',form='unformatted',
     &    access='direct',recl=4*nx*ny)
      call save_file12(10, bf2)
      open(10,file='../outputs_pft/bf.3.bin',
     &    status='unknown',form='unformatted',
     &    access='direct',recl=4*nx*ny)
      call save_file12(10, bf3)
      open(10,file='../outputs_pft/bf.4.bin',
     &    status='unknown',form='unformatted',
     &    access='direct',recl=4*nx*ny)
      call save_file12(10, bf4)
      open(10,file='../outputs_pft/bf.5.bin',
     &    status='unknown',form='unformatted',
     &    access='direct',recl=4*nx*ny)
      call save_file12(10, bf5)
      open(10,file='../outputs_pft/bf.6.bin',
     &    status='unknown',form='unformatted',
     &    access='direct',recl=4*nx*ny)
      call save_file12(10, bf6)
       open(10,file='../outputs_pft/bf.7.bin',
     &    status='unknown',form='unformatted',
     &    access='direct',recl=4*nx*ny)
      call save_file12(10, bf7)

      
      stop
      end program env

!     ==================================================

      
      subroutine pft_par(par, dt) !!!!!!!!!  mudamos os valores de dt
      implicit none
!     input
      integer, parameter :: vars = 7 
      integer :: par            ! parameter number 
      real, dimension(vars) :: dt,dt1,dt2,dt3,dt4,dt5,dt6
     &    ,dt7,dt8
      
      
!     dt1 = g1
!     dt2 = p21 
!     dt3 = aleaf
!     dt4 = aawood
!     dt5 = afroot
!     dt6 = tleaf
!     dt7 = tawood
!     dt8 = tfroot
      
!     PFTS
      
!     1 = Temperate Broad-leaf Deciduous Tree         33          4.43
!     2 = Tropical Broad-leaf Evergreen Tree          59	      3.72
!     3 = C3 grass	      31          4.50

!     PFT         1       2       3       4       5       6       7      
      data dt1/3.37,   4.8,    2.98,   4.6,    4.5,    4.64,   4.22/    !g1
      data dt2/5.9E-5, 3.4E-5, 3.3E-5, 3.1E-5, 3.1E-5, 3.0E-5, 5.5E-5/  !p21
      data dt3/0.35,   0.35,   0.40,   0.40,   0.45,   0.45,   0.30/    !aleaf
      data dt4/0.30,   0.30,   0.20,   0.15,   0.0,    0.0,    0.30/    !aawood
      data dt5/0.35,   0.35,   0.40,   0.45,   0.55,   0.55,   0.40/    !afroot
      data dt6/6.0,    6.5,    8.3,     3.5,    3.8,    3.8,   8.33/     !tleaf
      data dt7/70.,    68.0,   35.6,    4.0,    0.0,    0.0,    42.5/   !tawood
      data dt8/2.8,    3.0,    3.5,     3.0,    3.5,    3.5,    3.0/    !tfroot



!     PFT         1       2       3      4      5       6        7
c      data dt1/3.37,   4.645,  7.18,   2.98,   2.35,   4.64,   4.22/ !g1
c      data dt2/3.2E-5, 3.1E-5, 3.8E-5, 5.5E-5, 5.4E-5, 4.0E-5, 4.5E-5/ !p21
c      data dt3/0.70,   0.65,   0.70,   0.62,   0.55,   0.82,   0.60/ !aleaf
c      data dt4/0.10,   0.15,   0.15,   0.12,   0.0,    0.0,    0.10/ !aawood
c      data dt5/0.20,   0.20,   0.15,   0.26,   0.45,   0.18,   0.30/ !afroot
c      data dt6/6.8,    6.4,    2.4,     2.2,    2.8,    1.8,   4.9/ !tleaf
c      data dt7/80.,    58.0,   38.6,    8.0,    0.0,    0.0,   42.5/ !tawood
c      data dt8/2.8,    2.0,    2.0,     1.5,    1.4,    1.2,   3.0/ !tfroot
      
      if(par .eq. 1 ) then      ! g1
         dt(:) = dt1(:)
      else if(par .eq. 2) then  ! p21
         dt(:) = dt2(:)
      else if(par .eq. 3) then  ! aleaf
         dt(:) = dt3(:)
      else if(par .eq. 4) then  ! awood
         dt(:) = dt4(:)
      else if(par .eq. 5) then  ! afroot
         dt(:) = dt5(:)
      else if(par .eq. 6) then  ! tleaf
         dt(:) = dt6(:)
      else if(par .eq. 7) then  ! tawood
         dt(:) = dt7(:)
      else if(par .eq. 8) then  ! tfroot
         dt(:) = dt8(:)
      else
         print*, "your search failed"
      endif
      
      return
      end subroutine pft_par
      
c     ==================================================
      subroutine spinup(nppot,
     &     cleafini,cfrootini,cawoodini)
c     &     cbwoodini,cstoini,cotherini,crepini) 
      IMPLICIT NONE

      integer, parameter :: nt=10000
      integer, parameter :: npfts=7
      
c     inputs
      integer i6, kk, k
      
      real :: nppot
      real :: sensitivity

c     outputs
      real :: cleafini(npfts)
      real :: cawoodini(npfts)
      real :: cfrootini(npfts)
      real :: ocp(npfts)
      logical :: inutil(npfts)
      real*8 cleafi_aux(nt)
      real*8 cfrooti_aux(nt)
      real*8 cawoodi_aux(nt)

    
      real aleaf(npfts)             !npp percentage alocated to leaf compartment
      real aawood (npfts)           !npp percentage alocated to aboveground woody biomass compartment
      real afroot(npfts)            !npp percentage alocated to fine roots compartmentc 
      real tleaf(npfts)             !turnover time of the leaf compartment (yr)
      real tawood (npfts)           !turnover time of the aboveground woody biomass compartment (yr)
      real tfroot(npfts)            !turnover time of the fine roots compartment


      call pft_par(3, aleaf)
      call pft_par(4, aawood)
      call pft_par(5, afroot)
      call pft_par(6, tleaf)
      call pft_par(7, tawood)
      call pft_par(8, tfroot)

 
      sensitivity = 1.10
      if(nppot .le. 0.0) goto 200
      do i6=1,npfts
         do k=1,nt
            if (k.eq.1) then
               cleafi_aux (k) =  aleaf(i6)*(nppot)
               cawoodi_aux(k) = aawood(i6)*(nppot)
               cfrooti_aux(k) = afroot(i6)*(nppot)

            else
               if(aawood(i6) .gt. 0.0) then
                  cleafi_aux(k) = ((aleaf(i6)*(nppot))-
     &                (cleafi_aux(k-1)/(tleaf(i6)))) + cleafi_aux(k-1)
                  cawoodi_aux(k) = ((aawood(i6)*(nppot))-
     &                (cawoodi_aux(k-1)/(tawood(i6)))) + cawoodi_aux(k
     &                -1)
                  cfrooti_aux(k) = ((afroot(i6)*(nppot))-
     &                (cfrooti_aux(k-1)/(tfroot(i6)))) + cfrooti_aux(k
     &                -1)
               else
                  cleafi_aux(k) = ((aleaf(i6)*(nppot))-
     &                (cleafi_aux(k-1)/(tleaf(i6)))) + cleafi_aux(k-1)
                  cawoodi_aux(k) = 0.0
                  cfrooti_aux(k) = ((afroot(i6)*(nppot))-
     &                (cfrooti_aux(k-1)/(tfroot(i6)))) + cfrooti_aux(k
     &                -1)
               endif
                  

               
               kk =  nint(k*0.66)
               if(cawoodi_aux(kk) .gt. 0.0) then
                  if((cfrooti_aux(k)/cfrooti_aux(kk).lt.sensitivity)
     $                .and.(cleafi_aux(k)/cleafi_aux(kk).lt.sensitivity)
     $                .and.(cawoodi_aux(k)/cawoodi_aux(kk).lt.
     $                sensitivity)) then
                     
                     cleafini(i6) = real(cleafi_aux(k),4) ! carbon content (kg m-2)
                     cfrootini(i6) = real(cfrooti_aux(k),4)
                     cawoodini(i6) = real(cawoodi_aux(k),4)
                     exit
                  ENDIF
               else
                  if((cfrooti_aux(k)/cfrooti_aux(kk).lt.sensitivity)
     $                .and.(cleafi_aux(k)
     &                /cleafi_aux(kk).lt.sensitivity)) then
                     
                     cleafini(i6) = real(cleafi_aux(k),4) ! carbon content (kg m-2)
                     cfrootini(i6) = real(cfrooti_aux(k),4)
                     cawoodini(i6) = 0.0
                     exit
                  endif   
               endif
            endif
         enddo                  !nt
      enddo                     ! npfts 
      call pft_area_frac(cleafini, cfrootini, cawoodini, ocp, inutil)
      do i6 = 1,npfts
         cleafini(i6) = cleafini(i6) * ocp(i6)
         cleafini(i6) = cleafini(i6) * ocp(i6)
         cleafini(i6) = cleafini(i6) * ocp(i6)
      enddo
 200  continue
      return
      end subroutine spinup
!     ================================

      subroutine readx(nunit,var,x)
!     auxiliar reading routine
      integer,parameter :: nx=720,ny=360
      integer nunit,x
      real var(nx,ny,x)
      real aux(nx,ny)
      do k=1,x
         read(nunit,rec=k) aux
         do i=1,nx
            do j=1,ny
               var(i,j,k) = aux(i,j)
            enddo
         enddo
      enddo
      return
      end
!     ================================

      subroutine save_file12(nunit, var)
      parameter(nx = 720, ny = 360, nt = 12)
      integer i, j, k
      real var(nx,ny,nt)
      real waux(nx,ny)
      do k=1,nt
         do i=1,nx
            do j=1,ny
               waux(i,j) = var(i,j,k)
            enddo
         enddo
         write(nunit,rec=k) waux
      enddo
      close(nunit)
      return
      end subroutine save_file12
!     ================================

      subroutine savex(nunit, var, x)
      integer, parameter :: nx=720,ny=360
      integer i, j, k, x
      real var(nx,ny,x)
      real waux(nx,ny)
      do k=1,x
         do i=1,nx
            do j=1,ny
               waux(i,j) = var(i,j,k)
            enddo
         enddo
         write(nunit,rec=k) waux
      enddo
      close(nunit)
      return
      end subroutine savex

      
!     ===============================
      subroutine nan2ndt(var, nl)

      integer nl
      integer, parameter :: nx=720, ny=360
      integer i,j
      real, parameter :: no_data = -9999.0
      real var(nx,ny,nl)
      
      do i = 1,nx
         do j = 1,ny
            do k = 1,nl
               
               if(var(i,j,k) .gt. 1E12)then
                  var(i,j,k) = 0.0
c                  print*, 'var .gt. 1e32',i,j,k
               endif
               
               if(var(i,j,k) .lt. 1E-12) then
                  var(i,j,k) = 0.0
c                  print*, 'var .lt. 1e-32',i,j,k
               endif
               
               if(isnan(var(i,j,k))) then
                  var(i,j,k) = no_data
c                  print*, 'NaN found', '---place',i,j,k
               endif
               
            enddo
         enddo
      enddo  
      
      return
      end subroutine nan2ndt
