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
c     Reviewed by jpdarela  dec/2016

!     Last update: Dec/2016
!     
!     Compile with: g95 (or gfortran) env5TR.f wbm4TR.f productivity1.f allocation1.f -o a.exe
!     Execute with ./a.exe
!     =======================================================================
!     COMENTARIOS:
C     NOTEM QUE EU SALVEI NA MINHA PASTA INPUTS OS VALORES INICIAIS DE ALOCAÇÃO DE NPP (L 218-231)
C     PARA NAO TER QUE RODAR A SPINUP TODA VEZ QUE FOR RODAR O MODELO. (L 143-148)
C     A PARTE DO CALCULO DO SPINUP ESTA COMENTADA NESTE CODIGO (L 180-216) - JP
!     Parameters and variables
!     ------------------------
!     
      integer i,j,k,p
      integer, parameter :: nx=720,ny=360,q=3
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
c     secretas para o env5.f, exceto aquelas que são outputs na chamada.
      
      
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
     $rgl_pft,rgf_pft,rgs_pft,rg_pft,cleaf_pft,cawood_pft,cfroot_pft


c     variaveis do spinup
      real, dimension(q) :: aux1, aux2, aux3
      real npp_pot(nx,ny,12)
      real, dimension(nx,ny) :: aux_npp
      real npp_sca
      
c     -----FIM DA DEFINICAO DE VARIAVEIS PARA RODAR O MODELO--
C     
      
C     THESE WILL RECEIVE MEANS BETWEEN q PFTs and for each pft (ex. ph to mean; ph1 to pft 1)
      real, dimension(nx,ny,12) :: ph, ph1, ph2, ph3
      real, dimension(nx,ny,12) :: ar, ar1, ar2, ar3
      real, dimension(nx,ny,12) :: npp, npp1, npp2, npp3
      real, dimension(nx,ny,12) :: lai, lai1, lai2, lai3
      real, dimension(nx,ny,12) :: clit, clit1, clit2, clit3
      real, dimension(nx,ny,12) :: csoil, csoil1, csoil2, csoil3
      real, dimension(nx,ny,12) :: hr, hr1, hr2, hr3
      real, dimension(nx,ny,12) :: rcm, rcm1, rcm2, rcm3
      real, dimension(nx,ny,12) :: runom, runom1, runom2, runom3
      real, dimension(nx,ny,12) :: evaptr, evaptr1, evaptr2, evaptr3
      real, dimension(nx,ny,12) :: wsoil, wsoil1, wsoil2, wsoil3
      
C     NEW OUTPUTS (AUTOTRF RESPIRATION, ALLOCATION)
      
      real, dimension(nx,ny,12) :: rml, rml1, rml2, rml3
      real, dimension(nx,ny,12) :: rmf, rmf1, rmf2, rmf3 
      real, dimension(nx,ny,12) :: rms, rms1, rms2, rms3
      real, dimension(nx,ny,12) :: rm,  rm1,  rm2,  rm3
      real, dimension(nx,ny,12) :: rgl, rgl1, rgl2, rgl3
      real, dimension(nx,ny,12) :: rgf, rgf1, rgf2, rgf3
      real, dimension(nx,ny,12) :: rgs, rgs1, rgs2, rgs3
      real, dimension(nx,ny,12) :: rg,  rg1,  rg2,  rg3 
      real, dimension(nx,ny,12) :: cleaf,  cleaf1,  cleaf2,  cleaf3
      real, dimension(nx,ny,12) :: cawood, cawood1, cawood2, cawood3
      real, dimension(nx,ny,12) :: cfroot, cfroot1, cfroot2, cfroot3

      
C     -------END DECLARATION----------------------------------------

      
!     Open INPUT files
!     ==========
      open( 9,file='../inputs/lsmk.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)
      open(10,file='../inputs/ps.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)
      open(11,file='../inputs/pr.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)
      open(12,file='../inputs/tas.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)
      open(13,file='../inputs/rsds.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)
      open(26,file='../inputs/npp.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)
      
      open(21,file='../inputs/cleaf_ini.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)
      open(22,file='../inputs/cfroot_ini.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)
      open(23,file='../inputs/cawood_ini.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)
      

      
!     Read data
!     =========
      
       read (9,rec=1) lsmk
       call read12 (10,ps)
       call read12 (11,pr)
       call read12 (12,t)
       call read12 (13,ipar)
      call read12(26,npp_pot)

      call read3(21, cleafin)
      call read3(22,cfrootin)
      call read3(23,cawoodin)
     
!     Close files
!     ===========
!     
       close ( 9)
       close (10)
       close (11)
       close (12)
       close (13)
       close (26)
       close (21)
       close (22)
       close (23)
!     

c     fazendo medias da npp
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
c$$$       
c$$$       open(10,file='../inputs/npp2.bin',
c$$$     &      status='unknown',form='unformatted',
c$$$     &      access='direct',recl=4*nx*ny)
c$$$       write(10,rec=1) aux_npp
c$$$      close(10) 
c$$$      print*, 'npp_saved'
c$$$      
c$$$      
c$$$!     calling spinup
c$$$      do i=1,nx
c$$$         do j=1,ny
c$$$            if (nint(lsmk(i,j)) .ne. 0) then
c$$$               npp_sca = aux_npp(i,j)
c$$$               
c$$$               do p=1,q   
c$$$                  aux1(p) = 0.0
c$$$                  aux2(p) = 0.0
c$$$                  aux3(p) = 0.0
c$$$               enddo
c$$$               
c$$$               call spinup(npp_sca, aux1, aux2, aux3)
c$$$
c$$$               do p=1,q   
c$$$                  cleafin(i,j,p)  = aux1(p)
c$$$                  cfrootin(i,j,p) = aux2(p)
c$$$                  cawoodin(i,j,p) = aux3(p)
c$$$               enddo
c$$$            endif
c$$$         enddo
c$$$         print*, (real(i)/real(nx))*100.0, '%'
c$$$      enddo
c$$$
c$$$      open(10,file='../inputs/cleaf_ini.bin',
c$$$     &     status='unknown',form='unformatted',
c$$$     &     access='direct',recl=4*nx*ny)
c$$$      call save_file3(10, cleafin)
c$$$
c$$$      open(10,file='../inputs/cfroot_ini.bin',
c$$$     &     status='unknown',form='unformatted',
c$$$     &     access='direct',recl=4*nx*ny)
c$$$      call save_file3(10, cfrootin)
c$$$
c$$$      open(10,file='../inputs/cawood_ini.bin',
c$$$     &     status='unknown',form='unformatted',
c$$$     &     access='direct',recl=4*nx*ny)
c$$$      call save_file3(10, cawoodin)

!     ===========
!     
      do i=1,nx
         do j=1,ny
            do k=1,12
!     Photosynthetically active radiation (IPAR:Ein/m2/s)
!     Observed data from ISLSCP2
               
               par(i,j,k) = (ipar(i,j,k)/(2.18e5)) !Converting to Ein/m2/s
               temp(i,j,k) = t(i,j,k) !+ant(i,j,k) !uncomment to use future anomalies
               p0(i,j,k) = ps(i,j,k) * 0.01 ! transforamando de pascal pra mbar (kPa)
               prec(i,j,k) = pr(i,j,k) !+anpr(i,j,k) !+pr(i,j,k)*0.2 !uncomment to use future anomalies
               if (prec(i,j,k).lt.0.0) prec (i,j,k) = 0.0
               
            enddo
         enddo
      enddo
      
!     Atmospheric CO2 pressure (Pa) !Ppmv / Pa
      ca= 363/9.901             !Pa (=350 ppmv; 1961-1990)

      
!     =======================================
!     Calculate environmental variables (wbm)
!     =======================================

!     wbm definition(soh pra lembrar a ordem dos argumentos)
c     esta eh a ordem que esta na definicao em wbm4.f:
      
c      subroutine wbm (prec,temp,lsmk,p0,ca,par, 
c     &     cleaf_ini, cawood_ini, cfroot_ini   
c     &     emaxm, tsoil, photo_pft,aresp_pft,npp_pft,lai_pft,
c     &     clit_pft,csoil_pft, hresp_pft,rcm_pft,runom_pft,
c     &     evapm_pft,wsoil_pft,rml_pft,rmf_pft,rms_pft,rm_pft,rgl_pft
c     $     ,rgf_pft,rgs_pft,rg_pft,cleaf_pft,cawood_pft,cfroot_pft)
      
      call wbm (prec,temp,lsmk,p0,ca,par,cleafin,cawoodin,cfrootin,
     &     emaxm, tsoil, photo_pft,aresp_pft,npp_pft,lai_pft,
     &     clit_pft,csoil_pft, hresp_pft,rcm_pft,runom_pft,
     &     evapm_pft,wsoil_pft,rml_pft,rmf_pft,rms_pft,rm_pft,rgl_pft
     $     ,rgf_pft,rgs_pft,rg_pft,cleaf_pft,cawood_pft,cfroot_pft)   



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
                  cleaf(i,j,k)  = no_data
                  cawood(i,j,k)  = no_data
                  cfroot(i,j,k)  = no_data
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
                     cleaf(i,j,k)  = cleaf(i,j,k) + cleaf_pft(i,j,k,p)
                     cawood(i,j,k)  = cawood(i,j,k) + cawood_pft(i,j,k
     $                    ,p)
                     cfroot(i,j,k)  = cfroot(i,j,k) + cfroot_pft(i,j,k
     $                    ,p)
                     
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

      open(10,file='../outputs/cleaf.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, cleaf)

      open(10,file='../outputs/cawood.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, cawood)

      open(10,file='../outputs/cfroot.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, cfroot)



      do i=1,nx
         do j=1,ny
            do k=1,12
               if(nint(lsmk(i,j)) .ne. 0) then
                  npp1(i,j,k) = npp_pft(i,j,k,1)
                  npp2(i,j,k) = npp_pft(i,j,k,2)
                  npp3(i,j,k) = npp_pft(i,j,k,3)
                  
                  rcm1(i,j,k) = rcm_pft(i,j,k,1)
                  rcm2(i,j,k) = rcm_pft(i,j,k,2)
                  rcm3(i,j,k) = rcm_pft(i,j,k,3)

                  rm1(i,j,k) =rm_pft(i,j,k,1)
                  rm2(i,j,k) = rm_pft(i,j,k,2)
                  rm3(i,j,k) = rm_pft(i,j,k,3)

                  rg1(i,j,k) = rg_pft(i,j,k,1)
                  rg2(i,j,k) = rg_pft(i,j,k,2)
                  rg3(i,j,k) = rg_pft(i,j,k,3)

                  cleaf1(i,j,k) = cleaf_pft(i,j,k,1)
                  cleaf2(i,j,k) = cleaf_pft(i,j,k,2)
                  cleaf3(i,j,k) = cleaf_pft(i,j,k,3)

                  cfroot1(i,j,k) = cfroot_pft(i,j,k,1)
                  cfroot2(i,j,k) = cfroot_pft(i,j,k,2)
                  cfroot3(i,j,k) = cfroot_pft(i,j,k,3)
                  
                  cawood1(i,j,k) = cawood_pft(i,j,k,1)
                  cawood2(i,j,k) = cawood_pft(i,j,k,2)
                  cawood3(i,j,k) = cawood_pft(i,j,k,3)
               else
                  npp1(i,j,k) = no_data
                  npp2(i,j,k) = no_data
                  npp3(i,j,k) = no_data
                  
                  rcm1(i,j,k) = no_data
                  rcm2(i,j,k) = no_data
                  rcm3(i,j,k) = no_data
                  
                  rm1(i,j,k) = no_data
                  rm2(i,j,k) = no_data
                  rm3(i,j,k) = no_data

                  rg1(i,j,k) = no_data
                  rg2(i,j,k) = no_data
                  rg3(i,j,k) = no_data

                  cleaf1(i,j,k) = no_data 
                  cleaf2(i,j,k) = no_data
                  cleaf3(i,j,k) = no_data

                  cfroot1(i,j,k) = no_data
                  cfroot2(i,j,k) = no_data
                  cfroot3(i,j,k) = no_data
                  
                  cawood1(i,j,k) = no_data
                  cawood2(i,j,k) = no_data
                  cawood3(i,j,k) = no_data
               endif
            enddo
         enddo
      enddo
      
      open(10,file='../outputs_pft/npp.1.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, npp1)

      open(10,file='../outputs_pft/npp.2.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, npp2)

      open(10,file='../outputs_pft/npp.3.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, npp3)

      open(10,file='../outputs_pft/rcm.1.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, rcm1)

      open(10,file='../outputs_pft/rcm.2.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, rcm2)

      open(10,file='../outputs_pft/rcm.3.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, rcm3)

      open(10,file='../outputs_pft/rm.1.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, rm1)

      open(10,file='../outputs_pft/rm.2.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, rm2)

      open(10,file='../outputs_pft/rm.3.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, rm3)

      open(10,file='../outputs_pft/rg.1.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, rg1)

      open(10,file='../outputs_pft/rg.2.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, rg2)

      open(10,file='../outputs_pft/rg.3.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, rg3)

      open(10,file='../outputs_pft/cleaf.1.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, cleaf1)

      open(10,file='../outputs_pft/cleaf.2.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, cleaf2)

      open(10,file='../outputs_pft/cleaf.3.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, cleaf3)

      open(10,file='../outputs_pft/cawood.1.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, cawood1)

      open(10,file='../outputs_pft/cawood.2.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, cawood2)

      open(10,file='../outputs_pft/cawood.3.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, cawood3)

      open(10,file='../outputs_pft/cfroot.1.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, cfroot1)

      open(10,file='../outputs_pft/cfroot.2.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, cfroot2)

      open(10,file='../outputs_pft/cfroot.3.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, cfroot3)

      stop
      end program env
      

      subroutine spinup(nppot,
     &     cleafini,cfrootini,cawoodini)
c     &     cbwoodini,cstoini,cotherini,crepini) 

      IMPLICIT NONE

      integer, parameter :: nt=5000
      integer, parameter :: npft=3
c     inputs
      integer i6, kk, k
      
      real :: nppot
      real :: sensitivity,sensitivity2

      
c     outputs
      real :: cleafini(npft)
      real :: cawoodini(npft)
c      real :: cbwoodini(npft)
      real :: cfrootini(npft)
c      real :: cstoini(npft)
c      real :: cotherini(npft)
c      real :: crepini(npft)

c     internal vars

      real cleafi_aux(nt)
      real cfrooti_aux(nt)
      real cawoodi_aux(nt)
c      real cbwoodi_aux(nt)
c      real cstoi_aux(nt)
c      real cotheri_aux(nt)
c      real crepi_aux(nt)

c$$$      real aleaf(3)             !npp percentage alocated to leaf compartment
c$$$      data aleaf /0.40, 0.25, 0.45/
c$$$      real aawood (3)           !npp percentage alocated to aboveground woody biomass compartment
c$$$      data aawood /0.35, 0.40, 0.0/
c$$$      real afroot(3)            !npp percentage alocated to fine roots compartment
c$$$      data afroot /0.30, 0.25, 0.55/ 
c$$$      real tleaf(3)             !turnover time of the leaf compartment (yr)
c$$$      data tleaf /2.0, 0.7, 1.0/ 
c$$$      real tawood (3)           !turnover time of the aboveground woody biomass compartment (yr)
c$$$      data tawood /30.0, 3.0, 0.0/
c$$$      real tfroot(3)            !turnover time of the fine roots compartment
c$$$      data tfroot /3.0, 2.0, 1.0/

    
      real aleaf(3)             !npp percentage alocated to leaf compartment
      data aleaf /0.40,0.25,0.45/
      real aawood (3)           !npp percentage alocated to aboveground woody biomass compartment
      data aawood /0.35,0.40,0.0/
      real afroot(3)            !npp percentage alocated to fine roots compartment
      data afroot /0.30,0.25,0.55/ 
c      real abwood(3)            !npp percentage alocated to belowground woody biomass compartment
c      data abwood /0.10,0.10,0.001/
c      real asto(3)              !npp percentage alocated to storage compartment
c      data asto /0.10,0.10,0.10/
c      real arep(3)              !npp percentage alocated to reproduction compartment
c      data arep /0.15,0.15,0.10/
c      real aother(3)            !npp percentage alocated to other compartment
c      data aother /0.05,0.05,0.06/ 
c 
      real tleaf(3)             !turnover time of the leaf compartment (yr)
      data tleaf /2.0, 0.7, 1.0/ 
      real tawood (3)           !turnover time of the aboveground woody biomass compartment (yr)
      data tawood /30.0, 3.0, 0.0/
      real tfroot(3)            !turnover time of the fine roots compartment
      data tfroot /3.0,2.0,1.0/
c      real tbwood (3)           !turnover time of the belowground woody biomass compartment
c      data tbwood /40.0,40.0,40.0/
c      real tsto  (3)            !turnover time of the storage compartmentturn
c      data tsto /5.0,5.0,5.0/ 
c      real trep (3)             !turnover time of the reproduction compartment
c      data trep /0.25,0.25,0.25/ 
c      real tother (3)           !turnover time of the other compartment
c      data tother /0.12,0.12,0.12/

      sensitivity = 1.10
      sensitivity2 = 1.40

      do i6=1,npft
         do k=1,nt
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
                     exit
                  endif
               else
                  if((cfrooti_aux(k)/cfrooti_aux(kk).lt.
     $                 sensitivity).and.(cleafi_aux(k)
     $                 /cleafi_aux(kk).lt.sensitivity)) then
                     
                  
                     
                     cleafini(i6) = cleafi_aux(k) ! carbon content (kg m-2) 
                     cawoodini(i6) = 0.0
                     cfrootini(i6) = cfrooti_aux(k)
                     exit
                  endif
               ENDIF
            endif
         enddo
      enddo
      
      
      return
      end subroutine spinup				   
      
!     
      subroutine read12(nunit,var)
!     auxiliar reading routine
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
!     ================================
      subroutine read3(nunit,var)
!     auxiliar reading routine
      parameter(nx=720,ny=360)
      integer nunit
      real var(nx,ny,3)
      real aux(nx,ny)
      do k=1,3
         read(nunit,rec=k) aux
         do i=1,nx
            do j=1,ny
               var(i,j,k) = aux(i,j)
            enddo
         enddo
      enddo
      return
      end

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



      subroutine save_file3(nunit, var)
      parameter(nx = 720, ny = 360, nt = 3)
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
      end subroutine save_file3
