subroutine flip_image_1(input, a, b, output)

  implicit none

  integer, intent (in) :: a, b
  real(kind=4), intent (in), dimension (a,b) :: input
  real(kind=4), intent(out), dimension (a,b) :: output
  
  integer :: i, j
  do i = 1, a

     do j = b, 1, -1

        output(i, j) = input(i, (b + 1 - j))

     end do

  end do

end subroutine flip_image_1

subroutine flip_image_2(input, a, b, c, output)

  implicit none

  integer, intent (in) :: a, b, c
  real(kind=4),    intent (in), dimension (a, b, c) :: input
  real(kind=4),    intent(out), dimension (a, b, c) :: output
  
  integer :: i, j, k
  do k = 1, c
     do i = 1, a
        do j = b, 1, -1

           output(i, j, k) = input(i, (b + 1 - j), k)

        end do
     end do
  enddo
  

end subroutine flip_image_2



