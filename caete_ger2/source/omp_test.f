      program omp_test

      integer i,j,k
      real start, finish
      real arr(12,345,345)
      call cpu_time(start)

      do i = 1,12
         do j=1,345
            do k=1,345
               arr(i,j,k) = real(k)/322 - 234 * real(i)
               !print*, arr(i,j,k)
            enddo
         enddo
      enddo

      call cpu_time(finish)
      print*, finish - start
      stop
      end
      
