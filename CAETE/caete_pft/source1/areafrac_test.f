      program ocptest
      
      integer, parameter :: q = 3
      
      real c1(q), c2(q), c3(q)
      
      data c1 /1.23, 4.54, 3.54/
      data c2 /2.43, 2.43, 2.55/      
      data c3 /1.34, 5.56, 2.44/



      real ocp_coeffs(q)


      call pft_area_frac(c1,c2,c3,ocp_coeffs)

      do k = 1,q
         print*, ocp_coeffs(k), 'ocp', k

      enddo
      stop
      
      end program ocptest
