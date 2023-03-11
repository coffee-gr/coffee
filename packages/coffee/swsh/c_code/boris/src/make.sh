#!/usr/bin/env bash



#Generic Paths:
codeDir="$PWD"
compDir="$PWD"
libDir="$PWD/../lib"

#includes
incl="-I$codeDir"
#flags:
#Code should be position indep. anyway
fla="-fPIC -std=c99"

#Create compiled object code
echo "Compiling genArrManip.c"
gcc $fla -g -c -Wall genArrManip.c \
    -o "$compDir/genArrManip.o"
    
echo "Compiling genArrSorting.c"
gcc $fla -g -c -Wall genArrSorting.c \
    -ltime \
    -o "$compDir/genArrSorting.o"
    
echo "Compiling hdf5Interface.c"
gcc $fla -g -c -Wall hdf5Interface.c \
    -lhdf5 \
    -o "$compDir/hdf5Interface.o"
    
echo "Compiling hdf5rwDataStructs.c"
gcc $fla -g -c -Wall hdf5rwDataStructs.c \
    -o "$compDir/hdf5rwDataStructs.o" \
    -lmath

echo "Compiling mathWigner3jRecursion.c"
gcc $fla -g -c -Wall mathWigner3jRecursion.c \
    -o "$compDir/mathWigner3jRecursion.o"
    
echo "Compiling mathWigner3j.c"
gcc $fla -g -c -Wall mathWigner3j.c \
    -o "$compDir/mathWigner3j.o"

echo "Compiling mathClebschGordan.c"
gcc $fla -g -c -Wall mathClebschGordan.c \
    -o "$compDir/mathClebschGordan.o"

echo "Making shared library"   
gcc -std=c99 $fla -shared -Wl,-install_name,libboris.so.1 \
    -o "libboris.so.1.0.1" *.o -lhdf5 -lc
   
mv "libboris.so.1.0.1" "../../lib/"

#echo "Compiling main.c"
#gcc $fla -g -Wall main.c -L. -lw3jboris\
#    -Wl,-rpath,.\
#    -o "$compDir/main"
   

#-lmath -ltime

