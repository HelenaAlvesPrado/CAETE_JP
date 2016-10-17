program flip_image

implicit none

!VARIABLES TO GET CMDLINE ARGS
character *100 BUFFER  !HOLDS ALL CMD ARGS AS A STRING
character *20  mode
character *100 file_in, file_out
character *10 layers

integer nlay, aux_nlay
integer,parameter :: nx = 196
integer,parameter :: ny = 92
integer,parameter :: strd = int(4*nx*ny)
real*4 :: input_data(nx,ny), out_data(nx,ny)


!GET THE PARAMETERS FROM THE COMMAND LINE ARGUMENT
    call GETARG(1,BUFFER)
    read(BUFFER,*) mode
    
    call GETARG(2,BUFFER)
    read(BUFFER,*) layers

    read(layers,*,iostat=aux_nlay) nlay
    
    call GETARG(3,BUFFER)
    read(BUFFER,*) file_in

    if(mode.eq.'3d')then
       call GETARG(4,BUFFER)
       read(BUFFER,*) file_out
    endif
    
    if(mode.eq.'2d') then

       open(14,file=file_in,status='old', access='direct',recl=strd)
       read(14, rec=1) input_data
       close(14)

       call flip_image1(input_data, out_data, nx, ny)
       
       open(15,file=file_in,status='old',access='direct',recl=strd)
       write(15,rec=1) out_data
       close(15)
       
    else if(mode.eq.'3d') then

       open(14,file=file_in,status='old', access='direct',recl=strd)
       read(14, rec=1) input_data
       close(14)
       
       open(15,file=file_out,status='new',access='direct',recl=strd)
       write(15,rec=1) out_data
       close(15)

       call read_and_flip(14,15,nlay)

    else
       continue
       print*, "NADA FEITO"
       
    endif
    


 contains

      subroutine flip_image1(input, output, a, b)
        integer :: i, j

        integer, intent(in) :: a,b
        real*4, intent(in), dimension (a,b) :: input
        real*4, intent(out), dimension (a,b) :: output
        do i = 1, a
          do j = b, 1, -1
             output(i, j) = input(i, (b + 1 - j))
          end do
        end do
      end subroutine flip_image1

      subroutine read_and_flip(nunit_in, nunit_out, layers)
        ! auxiliar reading routine
        integer nx,ny,k
        parameter(nx=192,ny=96)
        integer, intent(in) :: nunit_in, nunit_out
        integer, intent(in) :: layers

        real aux(nx,ny),aux1(nx,ny)
        
        do k=1,layers
           read(nunit_in,rec=k) aux
           call flip_image1(aux, aux1, nx, ny)
           write(nunit_out,rec=k) aux1
        enddo
        return
      end subroutine read_and_flip
      
end program flip_image

