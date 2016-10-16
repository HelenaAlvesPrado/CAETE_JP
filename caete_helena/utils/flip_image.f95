subroutine flip_image_1(input, output, a, b)

  implicit none

  integer, intent (in) :: a, b
  real,    intent (in), dimension (0:a-1,0:b-1) :: input
  real,    intent(out), dimension (0:a-1,0:b-1) :: output
  
  integer :: i, j
  do i = 0, a-1

     do j = b, 0, -1

        output(i, j) = input(i, (b + 1 - j))

     end do

  end do

end subroutine flip_image_1

subroutine flip_image_2(input, output, a, b, c)

  implicit none

  integer, intent (in) :: a, b, c
  real,    intent (in), dimension (0:a-1,0:b-1,0:c-1) :: input
  real,    intent(out), dimension (0:a-1,0:b-1,0:c-1) :: output
  
  integer :: i, j, k
  do k = 0, c-1
     do i = 0, a-1
        do j = b, 0, -1

           output(i, j, k) = input(i, (b + 1 - j), k)

        end do
     end do
  enddo
  

end subroutine flip_image_2



