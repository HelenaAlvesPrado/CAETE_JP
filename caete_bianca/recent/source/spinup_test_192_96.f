c     12345
      program spin_test
c      implicit none


      integer, parameter :: nx = 192, ny = 96, npft = 3
      real, parameter :: no_data = -9999.0
      integer i, j, k
      real npp_pot(nx,ny)
      real lsmk(nx,ny), aux_npp(nx,ny)

      
      real, dimension(npft) :: aux1, aux2, aux3, aux4
      real, dimension(npft) :: aux5, aux6, aux7

c
c     carbon allocation in plant tissues
      
      real, dimension(nx,ny,npft) :: cleaf_ini = no_data
      real, dimension(nx,ny,npft) :: cfroot_ini = no_data
      real, dimension(nx,ny,npft) :: cawood_ini = no_data
      real, dimension(nx,ny,npft) :: cbwood_ini = no_data
      real, dimension(nx,ny,npft) :: csto_ini = no_data
      real, dimension(nx,ny,npft) :: cother_ini = no_data
      real, dimension(nx,ny,npft) :: crep_ini = no_data
      
      
      
      open(10,file='../inputs/lsmk.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)
      
      open(20,file='../inputs/nppot.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)
      
      read(10, rec=1) lsmk
      read(20, rec=1) npp_pot
      close(10)
      close(20)
      

      do i=1,nx
         do j=1,ny
            
            if (nint(lsmk(i,j)) .ne. 0) then
               
               npp_sca = aux_npp(i,j)
               
               aux1 = 0.0
               aux2 = 0.0
               aux3 = 0.0
               aux4 = 0.0
               aux5 = 0.0
               aux6 = 0.0
               aux7 = 0.0
               
               call spinup(npp_sca, aux1, aux2, aux3, aux4, aux5, aux6,
     $              aux7)
               
               cleaf_ini(nx,ny,:) = aux1
               cfroot_ini(nx,ny,:) = aux2
               cawood_ini(nx,ny,:) = aux3
               cbwood_ini(nx,ny,:) = aux4
               csto_ini(nx,ny,:) = aux5
               cother_ini(nx,ny,:) = aux6
               crep_ini(nx,ny,:) = aux7
            endif
         enddo
         print*, (real(i)/real(nx))*100.0, '%'
      enddo
      


      
      open(23,file='../outputs192/cleaf_ini_pft.bin',status='new'
     $     ,form='unformatted',access='direct',recl=4*nx*ny)
      
      open(24,file='../outputs192/cawood_ini_pft.bin',status='new',
     &     form='unformatted',access='direct',recl=4*nx*ny)

      open(25,file='../outputs192/cbwood_ini_pft.bin',status='new',
     &     form='unformatted',access='direct',recl=4*nx*ny)

      open(26,file='../outputs192/cfroot_ini_pft.bin',status='new',
     &     form='unformatted',access='direct',recl=4*nx*ny)
      
      open(27,file='../outputs192/crep_ini_pft.bin',status='new',
     &     form='unformatted',access='direct',recl=4*nx*ny)
      
      open(28,file='../outputs192/cother_ini_pft.bin',status='new',
     &     form='unformatted',access='direct',recl=4*nx*ny)
      
      open(29,file='../outputs192/csto_ini_pft.bin',status='new',
     &     form='unformatted',access='direct',recl=4*nx*ny)


            
         

      write(23, rec=1) cleaf_ini(:,:,1)
      write(23, rec=2) cleaf_ini(:,:,2)
      write(23, rec=3) cleaf_ini(:,:,3)

      write(26, rec=1) cfroot_ini(:,:,1)
      write(26, rec=2) cfroot_ini(:,:,2)
      write(26, rec=3) cfroot_ini(:,:,3)

      write(24, rec=1) cawood_ini(:,:,1)
      write(24, rec=2) cawood_ini(:,:,2)
      write(24, rec=3) cawood_ini(:,:,3)

      write(25, rec=1) cbwood_ini(:,:,1)
      write(25, rec=2) cbwood_ini(:,:,2)
      write(25, rec=3) cbwood_ini(:,:,3)

      write(29, rec=1) csto_ini(:,:,1)
      write(29, rec=2) csto_ini(:,:,2)
      write(29, rec=3) csto_ini(:,:,3)
 
      write(28, rec=1) cother_ini(:,:,1)
      write(28, rec=2) cother_ini(:,:,2)
      write(28, rec=3) cother_ini(:,:,3)

      write(27, rec=1) crep_ini(:,:,1)
      write(27, rec=2) crep_ini(:,:,2)
      write(27, rec=3) crep_ini(:,:,3)


      
      
      close(23)
      close(24)
      close(25)
      close(26)
      close(27)
      close(28)
      close(29)


   
      end program spin_test

      subroutine read12(nunit,var)
c     auxiliar reading routine
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
	  
