# MAKEFILE FOR CAETE
    # Define the symbols we might want to change:
    FC	        := gfortran
    FCFLAGS     := -g -Wall -Warray-bounds -trap-ffpe=overflow,underflow
    SOURCES     := env3.f wbm3.f prod3.f

    caete_env : $(SOURCES)
	$(FC) $(FCLAGS) $(SOURCES) -o caete_env.out




