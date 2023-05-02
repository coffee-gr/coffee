#  Select the default compiler.
ifndef COMPILER
  COMPILER := gfortran
endif
# COMPILER = pgf
# COMPILER = gfortran
# COMPILER = ifc
# COMPILER = ifort
# COMPILER = g95
# COMPILER = alpha
# COMPILER = sgi
# COMPILER = ifc

#  Uncomment this line if the compiler doesn't support 16-byte reals
# DEFS       := NOR16

FC       := ${COMPILER}
FFLAGS   :=

all:
	cd src; ${MAKE} FC=${FC} FFLAGS="${FFLAGS}"

test:
	cd tests; ../bin/rrfcalc < test.data

test3j:
	cd tests; ../bin/test3j < test3j.data > test3j.out

test6j:
	cd tests; ../bin/test6j < test6j.data > test6j.out

test9j:
	cd tests; ../bin/test9j < test9j.data > test9j.out

clean:
	cd src; ${MAKE} clean
	cd tests; rm *.out

