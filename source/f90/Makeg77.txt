ftn95 sdgvm0.f %1
ftn95 data.f  %1
ftn95 growth.f %1
ftn95 hydrology.f %1
ftn95 phenology.f %1
ftn95 func.f %1
ftn95 doly.f %1
ftn95 soil.f %1
ftn95 nppcalc.f %1
ftn95 light.f %1
ftn95 weathergenerator.f %1
ftn95 metdos.f %1
ftn95 sunshade.f %1
ftn95 sdgvm1.f %1
ftn95 parameter_adjustment.f %1

slink sdgvm0.obj sdgvm1.obj data.obj growth.obj parameter_adjustment.obj hydrology.obj phenology.obj func.obj doly.obj soil.obj nppcalc.obj light.obj weathergenerator.obj metdos.obj sunshade.obj

