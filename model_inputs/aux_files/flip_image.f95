subroutine flip_image(input, output, a, b)

  implicit none

  integer, intent (in) :: a, b
  real,    intent (in), dimension (a,b) :: input
  real,    intent(out), dimension (a,b) :: output
  
  integer :: i, j
  do i = 1, a

     do j = b, 1, -1

        output(i, j) = input(i, (b + 1 - j))

     end do

  end do

end subroutine flip_image

