To run the cohort version of SDGVM from this directory in linux using gfortran.

$ cd data
$ gunzip land_mask.dat.gz
$ gunzip land_mask.bin.gz
$ cd ..

$ cd source/f03
$ make clean
$ make
$ ./bin/sdgvm.exe ../../output/f03.dat

This will hopefuly run 9 sites and put the output in the
directory output/f03.

A link to the documentation is given in. 
'source/f03/docs/sdgvm.html'

NOTE
Although the source files have extensions of 'f90', the code uses some fortran 2003 features, so a fortran compiler supporting 2003 is needed.
Version 070607 is included and is in 'source/f90'.

