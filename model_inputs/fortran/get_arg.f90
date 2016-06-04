program get_args
  
  
  CHARACTER *10 file_in, file_out 

!DEFINE BUFFER HOLDS THE COMMAND LINE ARGUMENT
  CHARACTER *20 BUFFER  

!GET THE PARAMETERS FROM THE COMMAND LINE ARGUMENT
    CALL GETARG(1,BUFFER)
    READ(BUFFER,*) file_in
    CALL GETARG(2,BUFFER)
    READ(BUFFER,*) file_out
    
    WRITE(*,*) file_in
    WRITE(*,*) file_out
end program get_args