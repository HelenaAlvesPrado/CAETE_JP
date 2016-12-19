c     12345
      program spin_test
c      implicit none


      integer, parameter :: nx = 720, ny = 360, npft = 3
      integer i, j, k
      real npp_pot(nx,ny,12)
      real lsmk(nx,ny), aux_npp(nx,ny)

      real, dimension(npft) :: leaf_allocation_coeffs,
     $     leaf_turnover_coeffs

      data leaf_allocation_coeffs  /0.35,0.35,0.45/
      data leaf_turnover_coeffs /1.0,0.5,1.0/ 
      
c
c     carbon allocation in plant tissues
      real, dimension(nx,ny,npft) :: cleaf_ini = 0.0
c      real, dimension(nx,ny,npft) :: cfrootini1 = 0.0
c      real, dimension(nx,ny,npft) :: cawoodini1 = 0.0
c      real, dimension(nx,ny,npft) :: cbwoodini1 = 0.0
c      real, dimension(nx,ny,npft) :: cstoini1 = 0.0
c      real, dimension(nx,ny,npft) :: cotherini1 = 0.0
c      real, dimension(nx,ny,npft) :: crepini1 = 0.0


c     ler arquivo npp
c     call spinup

      
      open(10,file='../inputs/lsmk_30min.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)
      
      open(20,file='../inputs/npp.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)

      read(10, rec=1) lsmk
      call read12(20,npp_pot)
      close(10)
      close(20)

!     fazendo medias da npp
      do i =1,nx
         do j=1,ny
            aux_npp(i,j) = 0.0
            do k = 1,12
               aux_npp(i,j) = aux_npp(i,j) + (npp_pot(i,j,k)/12.) 
            enddo
         enddo
      enddo
      

       !cleafini_pft,cawoodini_pft,cfrootini_pft
c      call spinup2 (npp_pot,lsmk,
c     &     cleafini1,cawoodini1,cfrootini1)
c     &     cbwoodini1,cstoini1,cotherini1,crepini1)

      call spinup3(aux_npp, lsmk,leaf_allocation_coeffs,
     $     leaf_turnover_coeffs, .false.,cleaf_ini)

      
      open(23,file='../outputs/cleaf_ini_spinup30m_ok.bin',status='new'
     $     ,form='unformatted',access='direct',recl=4*nx*ny)
      
c      open(24,file='./outputs/cawood_ini_pft.bin',status='new',
c     &     form='unformatted',access='direct',recl=4*nx*ny)
c
c      open(25,file='./outputs/cbwood_ini_pft.bin',status='new',
c     &     form='unformatted',access='direct',recl=4*nx*ny)
c
c      open(26,file='./outputs/cfroot_ini_pft.bin',status='new',
c     &     form='unformatted',access='direct',recl=4*nx*ny)
c      
c      open(27,file='./outputs/crep_ini_pft.bin',status='new',
c     &     form='unformatted',access='direct',recl=4*nx*ny)
c      
c      open(28,file='./outputs/cother_ini_pft.bin',status='new',
c     &     form='unformatted',access='direct',recl=4*nx*ny)
c      
c      open(29,file='./outputs/csto_ini_pft.bin',status='new',
c     &     form='unformatted',access='direct',recl=4*nx*ny)


      write(23, rec=1) cleaf_ini(:,:,1)
      write(23, rec=2) cleaf_ini(:,:,2)
      write(23, rec=3) cleaf_ini(:,:,3)
c     write(24, rec=1)
c      write(24, rec=2)
c      write(24, rec=3)

c      write(25, rec=1)
c      write(25, rec=2)
c      write(25, rec=3)

c      write(26, rec=1)
c      write(26, rec=2)
c      write(26, rec=3)

c      write(27, rec=1)
c      write(27, rec=2)
c      write(27, rec=3)

 
c      write(28, rec=1)
c      write(28, rec=2)
c      write(28, rec=3)

c      write(29, rec=1)
c      write(29, rec=2)
c      write(29, rec=3)



      
      
      close(23)
c     close(24)
c      close(25)
c      close(26)
c      close(27)
c      close(28)
c      close(29)
      
      end program spin_test

      subroutine read12(nunit,var)
c     auxiliar reading routine
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
	  
