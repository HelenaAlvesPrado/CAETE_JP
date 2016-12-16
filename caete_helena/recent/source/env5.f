c234567
!     
      program env
!     
      implicit none
!     
!     =======================================================================
!     CPTEC-PVM2
!     Adapted from env4.f.
!     Code written by David Lapola and Helena Alves do Prado
c     Reviewed by jpdarela  dec/2016
!     Last update: Dec/2016
!     
!     Compile with: g95 (or gfortran) env5.f wbm4.f productivity1.f -o a.exe
!     Execute with ./a.exe
!     Then run g95 (gfortran) pvm5.f -o pvm5.exe
!     Execute with ./pvm5.exe
!     =======================================================================
!     
!     Parameters and variables
!     ------------------------
!     
      integer nx,ny,i,j,k,n,p,q
      parameter (nx=720,ny=360,q=3)
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


C     Aqui vc pode declarar as variaveis para armazenar as medias anuais
!     ... ex

c     real annual_csoil_pft1(nx,ny)
c     real csoil_pft1(nx,ny,12)
c     e assim por diante...
      
c     depois que vc chamar wbm vc faz estes calculos e salva os resultados
      ! tudo aqui no env5.f!



      
!     Open files
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


      
!     Read the data
!     =========
      
      read (9,rec=1) lsmk
      call read12 (10,ps)
      call read12 (11,pr)
      call read12 (12,t)
      call read12 (13,ipar)
c     call read12 (14,ant)
c     call read12 (15,anpr)
!     
!     Close files
!     ===========
!     
      close ( 9)
      close (10)
      close (11)
      close (12)
      close (13)
c     close (14)
c     close (15)
!     
!     ===========
!     
      do i=1,nx
            do j=1,ny
          
               do k=1,12
!     Photosynthetically active radiation reaching canopy (IPAR:Ein/m2/s)
!     Observed data from ISLSCP2
                  
                  par(i,j,k) = (ipar(i,j,k)/(2.18e5)) !Converting to Ein/m2/s
                  temp(i,j,k) = t(i,j,k) !+ant(i,j,k) !uncomment to use future anomalies
                  p0(i,j,k) = ps(i,j,k) * 0.01 ! transforamando de pascal pra mbar (kPa)
                  prec(i,j,k) = pr(i,j,k) !+anpr(i,j,k) !+pr(i,j,k)*0.2 !uncomment to use future anomalies
                  if (prec(i,j,k).lt.0.0) prec (i,j,k) = 0.0
                  
               enddo
            enddo
         enddo
      
!     
!     Atmospheric CO2 pressure (Pa)                                         !Ppmv / Pa
!     (1 Pa CO2 = 9.901 ppmv CO2 (Adams et al. 2004))
c     ca= 18.18                                                        !Pa (=180 ppmv; Last Glacial Maximum)
c     ca= 28.28                                                        !Pa (=280 ppmv; Pre-Industrial Rev.)
c     ca= 35.35                                                        !Pa (=350 ppmv; 1961-1990)
      ca= 363/9.901   !Pa (=350 ppmv; 1961-1990)
c     ca= 54.03                                                        !Pa (=535 ppmv; SRES-B1 2080's)
c     ca= 73.73                                                        !Pa (=730 ppmv; SRES-A2 2080's)
c     ca= ((73.73-35.35)*0.5)+35.35                                    !Pa half effect!
!     
         
!     
!     =======================================
!     Calculate environmental variables (wbm)
!     =======================================
!     
         
         
!     wbm definition(soh pra lembrar a ordem dos argumentos)
c     esta eh a ordem que esta na definicao em wbm4.f:
      
c            subroutine wbm (prec,temp,lsmk,p0,ca,par,
c     &     emaxm, tsoil, photo_pft,aresp_pft,npp_pft,lai_pft,
c     &     clit_pft,csoil_pft, hresp_pft,rcm_pft,runom_pft,
c     &     evapm_pft,wsoil_pft)   
         
         call wbm (prec,temp,lsmk,p0,ca,par,
     &     emaxm, tsoil, photo_pft,aresp_pft,npp_pft,lai_pft,
     &     clit_pft,csoil_pft, hresp_pft,rcm_pft,runom_pft,
     &     evapm_pft,wsoil_pft)   




C     AGORA JA TEMOS TODOS OS RESULTADOS DO MODELO
C     PODEMOS REALIZAR TODOS OS CALCULOS DE MEDIAS ANUAIS
C     E SALVAR OS RESULTADOS EM ARQUIVOS DE RASTER =P
!     
!     Write environmental variables
!     =============================
!     
!     open(16,file='../outputs/ambientais4_teste_pft.flt',
!     &        status='unknown',form='unformatted',
!     &        access='direct',recl=4*nx*ny)
!     write(16,rec=1) tmin
!     write(16,rec=2) meanpp
!     write(16,rec=3) seanpp
!     write(16,rec=4) meanhr
!     write(16,rec=5) meancs
!     write(16,rec=6) mphoto
!     write(16,rec=7) maresp
!     write(16,rec=8) ave_wsoil
!     write(16,rec=9) ave_evap
!     write(16,rec=10) ave_rc
c     close(16)
!
! Write monthly NPP 
! -----------------
!
!      open(17,file='../outputs/mnpp_caete_teste.flt',
!     &        status='unknown',form='unformatted',
!     &        access='direct',recl=4*nx*ny)
!      do k=1,12
!         do i=1,nx
!         do j=1,ny
!            waux(i,j) = mnpp(i,j,k)
!         enddo
!         enddo
!         write(17,rec=k) waux
!      enddo
!      close(17)
!      
! Write monthly photosynthesis 
! ----------------------------
!
!      open(18,file='../outputs/mon_photo_caete_teste.flt',
!     &        status='unknown',form='unformatted',
!     &        access='direct',recl=4*nx*ny)
!      do k=1,12
!         do i=1,nx
!         do j=1,ny
!            waux2(i,j) = mp(i,j,k)
!         enddo
!         enddo
!        write(18,rec=k) waux2
!      enddo
!      close(18)
!
! Write monthly plant respiration 
! -------------------------------
!
!      open(19,file='../outputs/mon_ar_caete_teste.flt',
!     &        status='unknown',form='unformatted',
!     &        access='direct',recl=4*nx*ny)
!      do k=1,12
!         do i=1,nx
!         do j=1,ny
!            waux3(i,j) = mar(i,j,k)
!         enddo
!         enddo
!         write(19,rec=k) waux3
!      enddo
!      close(19)
!
! Write monthly canopy resistence 
! -------------------------------
!
!      open(20,file='../outputs/mon_rc_caete_teste.flt',
!     &        status='unknown',form='unformatted',
!     &        access='direct',recl=4*nx*ny)
!      do k=1,12
!         do i=1,nx
!         do j=1,ny
!            waux4(i,j) = rc(i,j,k)
!         enddo
!         enddo
!         write(20,rec=k) waux4
!      enddo
!      close(20)    
!
! Write monthly canopy resistence_pft1 
! ------------------------------------
!
!      open(21,file='../outputs/mon_rc1_teste_pft.flt',
!     &        status='unknown',form='unformatted',
!     &        access='direct',recl=4*nx*ny)
!      do k=1,12
!         do i=1,nx
!         do j=1,ny
!            waux5(i,j) = rc_pft1(i,j,k)
!         enddo
!         enddo
!         write(21,rec=k) waux5
!      enddo
!      close(21)    
!
! Write monthly canopy resistence_pft2 
! ------------------------------------
!
!      open(22,file='../outputs/mon_rc2_teste_pft.flt',
!     &        status='unknown',form='unformatted',
!     &        access='direct',recl=4*nx*ny)
!      do k=1,12
!         do i=1,nx
!         do j=1,ny
!            waux6(i,j) = rc_pft2(i,j,k)
!         enddo
!         enddo
!         write(22,rec=k) waux6
!      enddo
!      close(22)    
!
! Write monthly canopy resistence_pft3 
! ------------------------------------
!
!      open(23,file='../outputs/mon_rc3_teste_pft.flt',
!     &        status='unknown',form='unformatted',
!     &        access='direct',recl=4*nx*ny)
!      do k=1,12
!         do i=1,nx
!         do j=1,ny
!            waux7(i,j) = rc_pft3(i,j,k)
!         enddo
!         enddo
!         write(23,rec=k) waux7
!      enddo
!      close(23)    
!
! Write annual canopy resistence 
! ------------------------------
!
!      open(24,file='../outputs/ave_rc_teste_550.flt',
!     &        status='unknown',form='unformatted',
!     &        access='direct',recl=4*nx*ny)
!         write(24,rec=1) ave_rc
!      close(24)
!
! Write annual npp 
! ----------------
!
!      open(25,file='../outputs/ave_npp_teste_pft.flt',
!     &        status='unknown',form='unformatted',
!     &        access='direct',recl=4*nx*ny)
!         write(25,rec=1) meanpp
!      close(25)
!
!      open(26,file='../outputs/meanpp1_teste_pft.flt',
!     &        status='unknown',form='unformatted',
!     &        access='direct',recl=4*nx*ny)
!         write(26,rec=1) meanpp1
!      close(26)
!
!      open(27,file='../outputs/meanpp2_teste_pft.flt',
!     &        status='unknown',form='unformatted',
!     &        access='direct',recl=4*nx*ny)
!         write(27,rec=1) meanpp2
!      close(27)
!
!      open(28,file='../outputs/meanpp3_teste_pft.flt',
!     &        status='unknown',form='unformatted',
!     &        access='direct',recl=4*nx*ny)
!         write(28,rec=1) meanpp3
!      close(28)            
!
!      open(29,file='../outputs/ave_rc1_teste_pft.flt',
!     &        status='unknown',form='unformatted',
!     &        access='direct',recl=4*nx*ny)
!         write(29,rec=1) ave_rc1
!      close(29) 
!
!      open(30,file='../outputs/ave_rc2_teste_pft.flt',
!     &        status='unknown',form='unformatted',
!     &        access='direct',recl=4*nx*ny)
!         write(30,rec=1) ave_rc2
!      close(30) 
!
!      open(31,file='../outputs/ave_rc3_teste_pft.flt',
!     &        status='unknown',form='unformatted',
!     &        access='direct',recl=4*nx*ny)
!         write(31,rec=1) ave_rc3
!      close(31) 
!
!      open(32,file='../outputs/meanpp1_p_teste_pft.flt',
!     &        status='unknown',form='unformatted',
!     &        access='direct',recl=4*nx*ny)
!         write(32,rec=1) meanpp1_p
!      close(32)
!
!      open(33,file='../outputs/meanpp2_p_teste_pft.flt',
!     &        status='unknown',form='unformatted',
!     &        access='direct',recl=4*nx*ny)
!         write(33,rec=1) meanpp2_p
!      close(33)
!
!      open(34,file='../outputs/meanpp3_p_teste_pft.flt',
!     &        status='unknown',form='unformatted',
!     &        access='direct',recl=4*nx*ny)
!         write(34,rec=1) meanpp3_p
!      close(34)            
!
!      open(35,file='../outputs/meanpp_final_teste_pft.flt',
!     &        status='unknown',form='unformatted',
!     &        access='direct',recl=4*nx*ny)
!         write(35,rec=1) meanpp_final
!      close(35)
!
!      open(36,file='../outputs/meanpp_final_novoteste.flt',
!     &        status='unknown',form='unformatted',
!     &        access='direct',recl=4*nx*ny)
!         write(36,rec=1) meanpp_final
!      close(36)
!      
!      open(37,file='../outputs/meanpp_pft_novoteste.flt',
!     &        status='unknown',form='unformatted',
!     &        access='direct',recl=4*nx*ny)
!      do q=1,3
!         waux8(i,j) = meanpp_pft(i,j,q)
!      enddo
!         write(37,rec=q) waux8
!      close(21) 
!
!     Program end
!
      stop
      end program env
      
!
! ================================
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
! ================================
