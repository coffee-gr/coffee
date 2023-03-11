SRC     = rrfdata.F90 rrflist.f90 factorials.f90 start.f90 \
          arithmetic.f90 conversions.f90 wigner.F90 input.F90

OBJ     = ${subst .f90,.o,${subst .F90,.o,${SRC}}}

LIB     = $(patsubst %.o,librrf.a(%.o),${OBJ})

all: rrfcalc test3j test6j test9j
	${MAKE} library

rrfcalc: modules
	${FC} ${FFLAGS} -o $@ $@.f90 ${OBJ}

test3j test6j test9j: modules
	${FC} ${FFLAGS} -o $@ $@.f90 ${OBJ}
	mv -f $@ ../bin

library: ${LIB}

modules: ${SRC}
	./compile.sh ${FC} ${FFLAGS}
	touch modules

librrf.a(%.o): %.o
	ar -r librrf.a $<

clean:
	rm -f *.o *.mod

force:
