c     12345
      program spin_test
c      implicit none


      integer, parameter :: nx = 192, ny = 96, npft = 3
      
      real npp_pot(nx,ny)
      real lsmk(nx,ny)
      
c     carbon allocation in plant tissues
      real, dimension(nx,ny,npft) :: cleafini1 = 0.0
      real, dimension(nx,ny,npft) :: cfrootini1 = 0.0
      real, dimension(nx,ny,npft) :: cawoodini1 = 0.0
      real, dimension(nx,ny,npft) :: cbwoodini1 = 0.0
      real, dimension(nx,ny,npft) :: cstoini1 = 0.0
      real, dimension(nx,ny,npft) :: cotherini1 = 0.0
      real, dimension(nx,ny,npft) :: crepini1 = 0.0


c     ler arquivo npp
c     call spinup

      
      open(10,file='../inputs/lsmk.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)
      
      open(20,file='../inputs/nppot.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)

      read(10, rec=1) lsmk
      read(20, rec=1) npp_pot
      close(10)
      close(20)


       !cleafini_pft,cawoodini_pft,cfrootini_pft
      call spinup2 (npp_pot,lsmk,
     &     cleafini1,cawoodini1,cfrootini1)
c     &     cbwoodini1,cstoini1,cotherini1,crepini1)


      open(23,file='../outputs/cleaf_ini_pft.bin',status='new',
     &     form='unformatted',access='direct',recl=4*nx*ny)
      
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


      write(23, rec=1) cleafini1(:,:,1)
      write(23, rec=2) cleafini1(:,:,2)
      write(23, rec=3) cleafini1(:,:,3)

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
