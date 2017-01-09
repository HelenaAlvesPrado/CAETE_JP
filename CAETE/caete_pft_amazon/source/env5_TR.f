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
      integer, parameter :: nx=120,ny=160,q=7
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

c variables related to carbon allocation and autothrophic respiration (bianca)
     
      real, dimension(nx,ny,12,q) :: rml_pft,rmf_pft,rms_pft,rm_pft,
     $rgl_pft,rgf_pft,rgs_pft,rg_pft
      real, dimension(nx,ny,q) :: cleaf_pft,cawood_pft,cfroot_pft


c     variaveis do spinup
      real, dimension(q) :: aux1, aux2, aux3
      real npp_pot(nx,ny,12)
      real, dimension(nx,ny) :: aux_npp
      real npp_sca
      
c     -----FIM DA DEFINICAO DE VARIAVEIS PARA RODAR O MODELO--
C     
      
C     THESE WILL RECEIVE MEANS BETWEEN q PFTs and for each pft (ex. ph to mean; ph1 to pft 1)
      real, dimension(nx,ny,12) :: ph!, ph1, ph2, ph3
      real, dimension(nx,ny,12) :: ar!, ar1, ar2, ar3
      real, dimension(nx,ny,12) :: npp!, npp1, npp2, npp3
      real, dimension(nx,ny,12) :: lai!, lai1, lai2, lai3
      real, dimension(nx,ny,12) :: clit!, clit1, clit2, clit3
      real, dimension(nx,ny,12) :: csoil!, csoil1, csoil2, csoil3
      real, dimension(nx,ny,12) :: hr!, hr1, hr2, hr3
      real, dimension(nx,ny,12) :: rcm!, rcm1, rcm2, rcm3
      real, dimension(nx,ny,12) :: runom!, runom1, runom2, runom3
      real, dimension(nx,ny,12) :: evaptr!, evaptr1, evaptr2, evaptr3
      real, dimension(nx,ny,12) :: wsoil!, wsoil1, wsoil2, wsoil3
      
C     NEW OUTPUTS (AUTOTRF RESPIRATION, ALLOCATION)
      
      real, dimension(nx,ny,12) :: rml!, rml1, rml2, rml3
      real, dimension(nx,ny,12) :: rmf!, rmf1, rmf2, rmf3 
      real, dimension(nx,ny,12) :: rms!, rms1, rms2, rms3
      real, dimension(nx,ny,12) :: rm!,  rm1,  rm2,  rm3
      real, dimension(nx,ny,12) :: rgl!, rgl1, rgl2, rgl3
      real, dimension(nx,ny,12) :: rgf!, rgf1, rgf2, rgf3
      real, dimension(nx,ny,12) :: rgs!, rgs1, rgs2, rgs3
      real, dimension(nx,ny,12) :: rg!,  rg1,  rg2,  rg3 
      real, dimension(nx,ny,12) :: cleaf!,  cleaf1,  cleaf2,  cleaf3
      real, dimension(nx,ny,12) :: cawood!, cawood1, cawood2, cawood3
      real, dimension(nx,ny,12) :: cfroot!, cfroot1, cfroot2, cfroot3
      
C     -------END DECLARATION----------------------------------------
      
!     Open INPUT files
!     ==========
      open( 9,file='../inputs/lsmk_sa.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)
      open(10,file='../inputs/ps_sa.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)
      open(11,file='../inputs/pr_sa.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)
      open(12,file='../inputs/tas_sa.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)
      open(13,file='../inputs/rsds_sa.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)
      open(26,file='../inputs/npp_sa.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)
      
!     Read data
!     =========
      
       read (9,rec=1) lsmk

       call readx(10,ps,12)
       call readx(11,pr,12)
       call readx(12,t,12)
       call readx(13,ipar,12)
       call readx(26,npp_pot,12)
     
!     Close files
!     ===========
!     
       close ( 9)
       close (10)
       close (11)
       close (12)
       close (13)
       close (26)

c      Calculating annual npp
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
      do i=1,nx
         do j=1,ny
            if (nint(lsmk(i,j)) .ne. 0) then
               npp_sca = aux_npp(i,j)
               
               do p=1,q   
                  aux1(p) = 0.0
                  aux2(p) = 0.0
                  aux3(p) = 0.0
               enddo
               
               call spinup(npp_sca, aux1, aux2, aux3)

               do p=1,q   
                  cleafin(i,j,p)  = aux1(p)
                  cfrootin(i,j,p) = aux2(p)
                  cawoodin(i,j,p) = aux3(p)
               enddo
            endif
         enddo
         print*, (real(i)/real(nx))*100.0, '%'
      enddo
!     ===========



      do i=1,nx
         do j=1,ny
            do k=1,12
!     Photosynthetically active radiation (IPAR:Ein/m2/s)
!     Observed data from ISLSCP2
               
               par(i,j,k) = ipar(i,j,k)/2.18e5 !Converting to Ein/m2/s
               temp(i,j,k) = t(i,j,k) !+ant(i,j,k) !uncomment to use future anomalies
               p0(i,j,k) = ps(i,j,k) * 0.01 ! transforamando de pascal pra mbar (kPa)
               prec(i,j,k) = pr(i,j,k) !+anpr(i,j,k) !+pr(i,j,k)*0.2 !uncomment to use future anomalies
               if (prec(i,j,k).lt.0.0) prec (i,j,k) = 0.0
               
            enddo
         enddo
      enddo
      
!     Atmospheric CO2 pressure (Pa) !Ppmv / Pa
      ca= 363/9.901             !Pa (=363 ppmv; 1981-2010)

      
!     =======================================
!     Calculate environmental variables (wbm)
!     =======================================
    
      call wbm (prec,temp,lsmk,p0,ca,par,cleafin,cawoodin,cfrootin,
     &     emaxm, tsoil, photo_pft,aresp_pft,npp_pft,lai_pft,
     &     clit_pft,csoil_pft, hresp_pft,rcm_pft,runom_pft,
     &     evapm_pft,wsoil_pft,rml_pft,rmf_pft,rms_pft,rm_pft,rgl_pft
     &     ,rgf_pft,rgs_pft,rg_pft,cleaf_pft,cawood_pft, cfroot_pft)   

!     SAVE RESULTS TO FILES
      open(10,file='../outputs/cleaf.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call savex(10, cleaf_pft, q)

      open(10,file='../outputs/cawood.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call savex(10,cawood_pft,q)

      open(10,file='../outputs/cfroot.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
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
                  cleaf(i,j,k)  = 0.0
                  cawood(i,j,k)  = 0.0
                  cfroot(i,j,k)  = 0.0
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
c$$$                 print*, npp(i,j,k)
               enddo
            endif

         enddo
      enddo
c$$$  

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

c
c
c
c      do i=1,nx
c         do j=1,ny
c            do k=1,12
c               if(nint(lsmk(i,j)) .ne. 0) then
c                  npp1(i,j,k) = npp_pft(i,j,k,1)
c                  npp2(i,j,k) = npp_pft(i,j,k,2)
c                  npp3(i,j,k) = npp_pft(i,j,k,3)
c                  
c                  rcm1(i,j,k) = rcm_pft(i,j,k,1)
c                  rcm2(i,j,k) = rcm_pft(i,j,k,2)
c                  rcm3(i,j,k) = rcm_pft(i,j,k,3)
c
c                  rm1(i,j,k) =rm_pft(i,j,k,1)
c                  rm2(i,j,k) = rm_pft(i,j,k,2)
c                  rm3(i,j,k) = rm_pft(i,j,k,3)
c
c                  rg1(i,j,k) = rg_pft(i,j,k,1)
c                  rg2(i,j,k) = rg_pft(i,j,k,2)
c                  rg3(i,j,k) = rg_pft(i,j,k,3)
c
c
c               else
c                  npp1(i,j,k) = no_data
c                  npp2(i,j,k) = no_data
c                  npp3(i,j,k) = no_data
c                  
c                  rcm1(i,j,k) = no_data
c                  rcm2(i,j,k) = no_data
c                  rcm3(i,j,k) = no_data
c                  
c                  rm1(i,j,k) = no_data
c                  rm2(i,j,k) = no_data
c                  rm3(i,j,k) = no_data
c
c                  rg1(i,j,k) = no_data
c                  rg2(i,j,k) = no_data
c                  rg3(i,j,k) = no_data
c
c               endif
c            enddo
c         enddo
c      enddo
c      
c      open(10,file='../outputs_pft/npp.1.bin',
c     &     status='unknown',form='unformatted',
c     &     access='direct',recl=4*nx*ny)
c      call save_file12(10, npp1)
c
c      open(10,file='../outputs_pft/npp.2.bin',
c     &     status='unknown',form='unformatted',
c     &     access='direct',recl=4*nx*ny)
c      call save_file12(10, npp2)
c
c      open(10,file='../outputs_pft/npp.3.bin',
c     &     status='unknown',form='unformatted',
c     &     access='direct',recl=4*nx*ny)
c      call save_file12(10, npp3)
c
c      open(10,file='../outputs_pft/rcm.1.bin',
c     &     status='unknown',form='unformatted',
c     &     access='direct',recl=4*nx*ny)
c      call save_file12(10, rcm1)
c
c      open(10,file='../outputs_pft/rcm.2.bin',
c     &     status='unknown',form='unformatted',
c     &     access='direct',recl=4*nx*ny)
c      call save_file12(10, rcm2)
c
c      open(10,file='../outputs_pft/rcm.3.bin',
c     &     status='unknown',form='unformatted',
c     &     access='direct',recl=4*nx*ny)
c      call save_file12(10, rcm3)
c
c      open(10,file='../outputs_pft/rm.1.bin',
c     &     status='unknown',form='unformatted',
c     &     access='direct',recl=4*nx*ny)
c      call save_file12(10, rm1)
c
c      open(10,file='../outputs_pft/rm.2.bin',
c     &     status='unknown',form='unformatted',
c     &     access='direct',recl=4*nx*ny)
c      call save_file12(10, rm2)
c
c      open(10,file='../outputs_pft/rm.3.bin',
c     &     status='unknown',form='unformatted',
c     &     access='direct',recl=4*nx*ny)
c      call save_file12(10, rm3)
c
c      open(10,file='../outputs_pft/rg.1.bin',
c     &     status='unknown',form='unformatted',
c     &     access='direct',recl=4*nx*ny)
c      call save_file12(10, rg1)
c
c      open(10,file='../outputs_pft/rg.2.bin',
c     &     status='unknown',form='unformatted',
c     &     access='direct',recl=4*nx*ny)
c      call save_file12(10, rg2)
c
c      open(10,file='../outputs_pft/rg.3.bin',
c     &     status='unknown',form='unformatted',
c     &     access='direct',recl=4*nx*ny)
c      call save_file12(10, rg3)

      stop
      end program env

!     ==================================================

      
      subroutine pft_par(par, dt)
      
!     input
      integer, parameter :: vars = 7 
      integer :: par            ! parameter number 
      real, dimension(vars) :: dt,dt1,dt2,dt3,dt4,dt5,dt6,dt7,dt8
      
      
!     dt1 = g1
!     dt2 = p21 
!     dt3 = aleaf
!     dt4 = aawood
!     dt5 = afroot
!     dt6 = tleaf
!     dt7 = tawood
!     dt8 = tfroot
      
!     PFTS
      
!     1 = tropical evergreen tree
!     2 = tropical deciduous-forest-tree
!     3 = tropical deciduous-savana-tree
!     4 = tropical herb
!     5 = tropical grass
!     6 = temperate tree
!     7 = temperate herb

!     PFT         1       2       3       4       5       6       7 
      data dt1 /3.04,   2.67,   2.0,    3.0,    1.4,    2.05,   1.95/
      data dt2 /3.2e-5, 3.1e-5, 3.0e-5, 3.3e-5, 2.8e-5, 5.0e-5, 4.0e-5/      
      data dt3 /0.30,   0.38,   0.36,   0.45,   0.55,   0.35,   0.55/
      data dt4 /0.28,   0.30,   0.27,   0.10,   0.0,    0.40,   0.12/
      data dt5 /0.42,   0.32,   0.37,   0.55,   0.45,   0.25,   0.23/
      data dt6 /7.0,    2.0,    1.0,    2.0,    1.0,    1.0,    2.0 /
      data dt7 /35.0,   30.0,   25.0,   2.0,    0.0,    32.0,   2.5/
      data dt8 /4.0,    3.5,    3.8,    1.5,    1.0,    3.8,    2.2/ 
     
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
      real :: sensitivity,sensitivity2

c     outputs
      real :: cleafini(npfts)
      real :: cawoodini(npfts)
      real :: cfrootini(npfts)

      real cleafi_aux(nt)
      real cfrooti_aux(nt)
      real cawoodi_aux(nt)

    
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
      sensitivity2 = 1.40

      do i6=1,npfts
         do k=1,nt
            !print*, k, i6
            if (k.eq.1) then
               cleafi_aux(k) = aleaf(i6)*(nppot)
               cawoodi_aux(k) = aawood(i6)*(nppot)
               cfrooti_aux(k) = afroot(i6)*(nppot)
c               cbwoodi_aux(k) = abwood(i6)*(nppot)
c               cstoi_aux(k) = asto(i6)*(nppot)
c               cotheri_aux(k) = aother(i6)*(nppot/365)
c               crepi_aux(k) = arep(i6)*(nppot/365)
            else
               cleafi_aux(k) = ((aleaf(i6)*(nppot))-
     &              (cleafi_aux(k-1)/(tleaf(i6)))) + cleafi_aux(k-1)
               cawoodi_aux(k) = ((aawood(i6)*(nppot))-
     &              (cawoodi_aux(k-1)/(tawood(i6)))) + cawoodi_aux(k-1)
               cfrooti_aux(k) = ((afroot(i6)*(nppot))-
     &              (cfrooti_aux(k-1)/(tfroot(i6)))) + cfrooti_aux(k-1)
c               cbwoodi_aux(k) = ((abwood(i6)*(nppot))-
c     &              (cbwoodi_aux(k-1)/(tbwood(i6)))) + cbwoodi_aux(k-1)
c               cstoi_aux(k) = ((asto(i6)*(nppot))-
c     &              (cstoi_aux(k-1)/(tsto(i6)))) + cstoi_aux(k-1)
c               cotheri_aux(k) = ((aother(i6)*(nppot))-
c     &              (cotheri_aux(k-1)/(tother(i6)*365))) + cotheri_aux(k
c     $              -1)
c               crepi_aux(k) = ((arep(i6)*(nppot))-
c     &              (crepi_aux(k-1)/(trep(i6)*365))) + crepi_aux(k-1)
               
               kk =  int(k*0.66)

               if(aawood(i6) .gt. 0.0) then
               
                  if((cfrooti_aux(k)
     $               /cfrooti_aux(kk).lt.sensitivity).and.(cleafi_aux(k)
     $               /cleafi_aux(kk).lt.sensitivity).and.(cawoodi_aux(k)
     $               /cawoodi_aux(kk).lt.sensitivity2)) then
                  
                  
                  
                     cleafini(i6) = cleafi_aux(k) ! carbon content (kg m-2) 
                     cawoodini(i6) = cawoodi_aux(k)
                     cfrootini(i6) = cfrooti_aux(k)
                     !print*, 'saindo', 'pft:', i6, k
                     exit
                     
                  endif
               else
                  if((cfrooti_aux(k)/cfrooti_aux(kk).lt.
     $                 sensitivity).and.(cleafi_aux(k)
     $                 /cleafi_aux(kk).lt.sensitivity)) then
                     
                  
                     
                     cleafini(i6) = cleafi_aux(k) ! carbon content (kg m-2) 
                     cawoodini(i6) = 0.0
                     cfrootini(i6) = cfrooti_aux(k)
                     !print*, 'saindo', 'pft:', i6, k
                     exit
                  endif
               ENDIF
            endif
            
         enddo
      enddo
      
      
      return
      end subroutine spinup
!     ================================

      subroutine readx(nunit,var,x)
!     auxiliar reading routine
      integer,parameter :: nx=120,ny=160
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
      parameter(nx = 120, ny = 160, nt = 12)
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
      integer, parameter :: nx=120,ny=160
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
