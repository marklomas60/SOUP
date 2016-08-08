To run the cohort version of SDGVM from this directory in linux using gfortran.

$ cd data
$ gunzip land_mask.bin
$ gunzip land_mask.bin
$ cd ..

$ cd source/f03
$ make clean
$ make
$ ./bin/sdgvm.exe ../../output/f03.dat

This will hopefuly run 9 sites and put the output in the
directory output/f03

NOTE
Although the source files have extensions of 'f90', the code uses some fortran 2003 so a fortran compiler supporting 2003 is needed.
Version 070607 is included and is in 'source/f90'.

