program flip_image

implicit none

!VARIABLES TO GET CMDLINE ARGS
character *100 BUFFER  !HOLDS ALL CMD ARGS AS A STRING
character *100 file_in

integer,parameter :: nx = 720
integer,parameter :: ny = 360
integer,parameter :: strd = int(4*nx*ny)
real*4 :: input_data(nx,ny), out_data(nx,ny)


!GET THE PARAMETERS FROM THE COMMAND LINE ARGUMENT
    call GETARG(1,BUFFER)
    read(BUFFER,*) file_in

open(14,file=file_in,status='old', access='direct',recl=strd)
    read(14, rec=1) input_data
close(14)

call flip_image1(input_data, out_data, nx, ny)

open(15,file=file_in,status='old',access='direct',recl=strd)
    write(15,rec=1) out_data
close(15)


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
end program flip_image
