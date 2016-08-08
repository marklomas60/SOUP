FFLAGS= -O3 -fcheck=all -g
#FFLAGS= -g -fbounds-check -Wuninitialized -O -ftrapv -fimplicit-none -fno-automatic
FF = gfortran

OBJ = sdgvm0.o sdgvm1.o data.o growth.o parameter_adjustment.o hydrology.o phenology.o func.o doly.o soil.o nppcalc.o light.o sunshade.o weathergenerator.o metdos.o

sdgvm0:	$(OBJ)
	$(FF) $(FFLAGS)  -o sdgvm0 $(OBJ)

%.o: %.f
	$(FF)  $(FFLAGS)  -c $<



clean:
	rm -f *.o *.il *.OBJ *~ *.exe sdgvm0 a.out *.MOD

