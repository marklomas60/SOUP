del *.obj *.exe *.mod

set st1=/dreal /debug

ftn95 misc_values.f90 %st1%
ftn95 dims.f90 %st1%
ftn95 system_state.f90 %st1%
ftn95 pft_parameters.f90 %st1%
ftn95 site_parameters.f90 %st1%
ftn95 bfire.f90 %st1%
ftn95 sdgvm0.f90 %st1%
ftn95 sdgvm1.f90 %st1%
ftn95 read_input.f90 %st1%
ftn95 open_files.f90 %st1%
ftn95 data.f90 %st1%
ftn95 veg_dynamics.f90 %st1%
ftn95 hydrology.f90 %st1%
ftn95 phenology.f90 %st1%
ftn95 func.f90 %st1%
ftn95 doly.f90 %st1%
ftn95 soil.f90 %st1%
ftn95 nppcalc.f90 %st1%
ftn95 light.f90 %st1%
ftn95 weather_generator.f90 %st1%
ftn95 metdos.f90 %st1%
ftn95 sunshade.f90 %st1%
ftn95 state_methods.f90 %st1%

slink sdgvm0.obj misc_values.obj bfire.obj state_methods.obj dims.obj system_state.obj site_parameters.obj pft_parameters.obj sdgvm1.obj read_input.obj open_files.obj data.obj veg_dynamics.obj hydrology.obj phenology.obj func.obj doly.obj soil.obj nppcalc.obj light.obj weather_generator.obj metdos.obj sunshade.obj

