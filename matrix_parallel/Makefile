ALLSOURCE=$(wildcard *.f90)
LIBDEVICE_OBJ=$(patsubst %.f90,%.o,$(ALLSOURCE))
# Files included by the sources here: every object is rebuilt when one changes.
DEPENDALL=size_parallel.h ../Fourier.inc ../unit_number.inc
MODFILES=mtmod.mod
include ../Makefile.inc

libmparallel.a: ${MODFILES} ${LIBDEVICE_OBJ}
	${AR} ${ARFLAGS} $@ $^
	${RANLIB} $@

mtmod.mod: mt19937.f90
	${FC} -c ${FCFLAGS} $<

# The unrolled .f kernels are included only by these two sources.
Multiply_Dirac.o: multiplication_direct_1.f multiplication_direct_2.f \
	multiplication_direct_3.f multiplication_direct_4.f
Multiply_Dirac_dagger.o: multiplication_dagger_direct_1.f multiplication_dagger_direct_2.f \
	multiplication_dagger_direct_3.f multiplication_dagger_direct_4.f

clean:
	echo cleaning up in hostsource.
	$(RM) -f *.o *.a *.mod
