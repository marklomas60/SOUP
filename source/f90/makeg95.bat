g95 sdgvm0.f %1
g95 data.f  %1
g95 growth.f %1
g95 hydrology.f %1
g95 phenology.f %1
g95 func.f %1
g95 doly.f %1
g95 soil.f %1
g95 nppcalc.f %1
g95 light.f %1
g95 weathergenerator.f %1
g95 metdos.f %1
g95 sunshade.f %1
g95 sdgvm1.f %1
g95 parameter_adjustment.f %1

g95 -o sdgvm0 sdgvm0.o sdgvm1.o data.o growth.o parameter_adjustment.o hydrology.o phenology.o func.o doly.o soil.o nppcalc.o light.o weathergenerator.o metdos.o sunshade.o

