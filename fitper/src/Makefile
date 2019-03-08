#IFORT
#F77 = ifort
#FFLAGS = -132 -r8 -save -free #-zero
#FFLAGS = -132 -r8 -save -w -noarg_temp_created -zero -fpconstant -traceback -g

#GFORTRAN
F77 = gfortran
FFLAGS = -ffree-line-length-132 -fdefault-real-8 -fdefault-double-8 -fno-automatic -w -finit-local-zero -fbacktrace -g -fbounds-check

fitperiods: modules.f90 subroutines.f90 fitseq_subroutine.f90 match_periods_subroutine.f90 fit.subroutine.f90 main.f90
	$(F77) $(FFLAGS) -o ../../$@ $^ 

clean:		
		rm *.o *.mod *~

