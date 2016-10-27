c     
      program maior_de_3
      
      integer  len
      parameter (len=3)         ! numero de suas proporcoes de biomassa
      integer i
      real  x(len)              ! aqui o seu vetor com as proporcoes
      real  max_value           ! variavel que vai guardar o maior valor
      
      max_value = 0.            ! criando com um valor mais baixo 
      
      data x /0.25, 0.25, 0.50/ ! criando valores teoricos pras suas prop.
      
      do i = 1,len
         if(x(i) .gt. max_value) max_value = x(i)
      enddo
      
      print*, max_value
      end program maior_de_3
      
      
      
      
