# ============================================================================
# Name        : Makefile
# Author      : 
# Version     :
# Copyright   : Your copyright notice
# Description : Makefile for Hello World in Fortran
# ============================================================================

.PHONY: all clean

# Change this line if you are using a different Fortran compiler
#FORTRAN_COMPILER = gfortran
FC = gfortran -O3 -fcheck=all -fdefault-real-8
#FC = gfortran -O3 -fcheck=all

s = src
m = mod
o = obj
b = bin
e = f90

all: $(o)/real_precision.o $(o)/dims.o $(o)/pft_parameters.o $(o)/system_state.o \
 $(o)/site_parameters.o $(o)/tuning_parameters.o $(o)/misc_parameters.o \
 $(o)/misc_values.o $(o)/metdos.o $(o)/func.o $(o)/state_methods.o \
 $(o)/open_files.o $(o)/data.o $(o)/weather_generator.o $(o)/veg_dynamics.o \
 $(o)/sdgvm1.o $(o)/read_input.o $(o)/hydrology_methods.o \
 $(o)/phenology_methods.o $(o)/sunshade.o $(o)/productivity_methods.o \
 $(o)/doly.o $(o)/soil_methods.o $(o)/light_methods.o $(o)/sdgvm.o
	$(FC) -o $b/sdgvm.exe \
 $(o)/real_precision.o $(o)/dims.o $(o)/pft_parameters.o $(o)/system_state.o \
 $(o)/site_parameters.o $(o)/tuning_parameters.o $(o)/misc_parameters.o \
 $(o)/misc_values.o $(o)/metdos.o $(o)/func.o $(o)/state_methods.o \
 $(o)/open_files.o $(o)/data.o $(o)/weather_generator.o $(o)/veg_dynamics.o \
 $(o)/sdgvm1.o $(o)/read_input.o $(o)/hydrology_methods.o \
 $(o)/phenology_methods.o $(o)/sunshade.o $(o)/productivity_methods.o \
 $(o)/doly.o $(o)/soil_methods.o $(o)/light_methods.o $(o)/sdgvm.o -Jmod

$(o)/real_precision.o: $(s)/real_precision.$(e)
	$(FC) -o $(o)/real_precision.o -c $(s)/real_precision.$(e) -J$(m)

$(o)/dims.o: $(s)/dims.$(e)
	$(FC) -o $(o)/dims.o -c $(s)/dims.$(e) -J$(m)

$(o)/pft_parameters.o: $(s)/pft_parameters.$(e) $(s)/real_precision.$(e) $(s)/dims.$(e)
	$(FC) -o $(o)/pft_parameters.o -c $(s)/pft_parameters.$(e) -J$(m)

$(o)/system_state.o: $(s)/system_state.$(e) $(s)/real_precision.$(e) $(s)/dims.$(e)
	$(FC) -o $(o)/system_state.o -c $(s)/system_state.$(e) -J$(m)

$(o)/site_parameters.o: $(s)/site_parameters.$(e) $(s)/real_precision.$(e) $(s)/dims.$(e)
	$(FC) -o $(o)/site_parameters.o -c $(s)/site_parameters.$(e) -J$(m)

$(o)/tuning_parameters.o: $(s)/tuning_parameters.$(e) $(s)/real_precision.$(e)
	$(FC) -o $(o)/tuning_parameters.o -c $(s)/tuning_parameters.$(e) -J$(m)

$(o)/misc_parameters.o: $(s)/misc_parameters.$(e) $(s)/real_precision.$(e)
	$(FC) -o $(o)/misc_parameters.o -c $(s)/misc_parameters.$(e) -J$(m)

$(o)/misc_values.o: $(s)/misc_values.$(e) $(s)/real_precision.$(e)
	$(FC) -o $(o)/misc_values.o -c $(s)/misc_values.$(e) -J$(m)

$(o)/metdos.o: $(s)/metdos.$(e) $(s)/real_precision.$(e)
	$(FC) -o $(o)/metdos.o -c $(s)/metdos.$(e) -J$(m)

$(o)/func.o: $(s)/func.$(e) $(s)/real_precision.$(e) $(s)/pft_parameters.$(e)
	$(FC) -o $(o)/func.o -c $(s)/func.$(e) -J$(m)

$(o)/state_methods.o: $(s)/state_methods.$(e) $(s)/real_precision.$(e) $(s)/dims.$(e) $(s)/system_state.$(e) $(s)/site_parameters.$(e) $(s)/pft_parameters.$(e) $(s)/misc_values.$(e) $(s)/tuning_parameters.$(e)
	$(FC) -o $(o)/state_methods.o -c $(s)/state_methods.$(e) -J$(m)

$(o)/open_files.o: $(s)/open_files.$(e) $(s)/real_precision.$(e) $(s)/func.$(e) $(s)/pft_parameters.$(e)
	$(FC) -o $(o)/open_files.o -c $(s)/open_files.$(e) -J$(m)

$(o)/data.o: $(s)/data.$(e) $(s)/real_precision.$(e) $(s)/dims.$(e) $(s)/func.$(e)
	$(FC) -o $(o)/data.o -c $(s)/data.$(e) -J$(m)

$(o)/weather_generator.o: $(s)/weather_generator.$(e) $(s)/real_precision.$(e) $(s)/data.$(e) $(s)/metdos.$(e)
	$(FC) -o $(o)/weather_generator.o -c $(s)/weather_generator.$(e) -J$(m)

$(o)/veg_dynamics.o: $(s)/veg_dynamics.$(e) $(s)/real_precision.$(e) $(s)/dims.$(e)  $(s)/system_state.$(e) $(s)/site_parameters.$(e) $(s)/state_methods.$(e) $(s)/func.$(e) $(s)/tuning_parameters.$(e)
	$(FC) -o $(o)/veg_dynamics.o -c $(s)/veg_dynamics.$(e) -J$(m)

$(o)/sdgvm1.o: $(s)/sdgvm1.$(e) $(s)/real_precision.$(e) $(s)/dims.$(e) $(s)/data.$(e) $(s)/weather_generator.$(e) $(s)/veg_dynamics.$(e) $(s)/func.$(e)
	$(FC) -o $(o)/sdgvm1.o -c $(s)/sdgvm1.$(e) -J$(m)

$(o)/read_input.o: $(s)/read_input.$(e) $(s)/real_precision.$(e) $(s)/pft_parameters.$(e) $(s)/site_parameters.$(e) $(s)/sdgvm1.$(e) $(s)/func.$(e) $(s)/misc_parameters.$(e) $(s)/tuning_parameters.$(e) $(s)/open_files.$(e)
	$(FC) -o $(o)/read_input.o -c $(s)/read_input.$(e) -J$(m)

$(o)/hydrology_methods.o: $(s)/hydrology_methods.$(e) $(s)/real_precision.$(e) $(s)/system_state.$(e) $(s)/site_parameters.$(e) $(s)/tuning_parameters.$(e) $(s)/pft_parameters.$(e)
	$(FC) -o $(o)/hydrology_methods.o -c $(s)/hydrology_methods.$(e) -J$(m)

$(o)/phenology_methods.o: $(s)/phenology_methods.$(e) $(s)/real_precision.$(e) $(s)/dims.$(e) $(s)/system_state.$(e) $(s)/pft_parameters.$(e) $(s)/site_parameters.$(e) $(s)/misc_values.$(e)
	$(FC) -o $(o)/phenology_methods.o -c $(s)/phenology_methods.$(e) -J$(m)

$(o)/sunshade.o: $(s)/sunshade.$(e) $(s)/real_precision.$(e)
	$(FC) -o $(o)/sunshade.o -c $(s)/sunshade.$(e) -J$(m)

$(o)/productivity_methods.o: $(s)/productivity_methods.$(e) $(s)/real_precision.$(e) $(s)/dims.$(e) $(s)/sunshade.$(e) $(s)/system_state.$(e) $(s)/pft_parameters.$(e) $(s)/site_parameters.$(e) $(s)/tuning_parameters.$(e)
	$(FC) -o $(o)/productivity_methods.o -c $(s)/productivity_methods.$(e) -J$(m)

$(o)/doly.o: $(s)/doly.$(e) $(s)/real_precision.$(e) $(s)/phenology_methods.$(e) $(s)/productivity_methods.$(e) $(s)/misc_values.$(e) $(s)/site_parameters.$(e) $(s)/tuning_parameters.$(e)
	$(FC) -o $(o)/doly.o -c $(s)/doly.$(e) -J$(m)

$(o)/soil_methods.o: $(s)/soil_methods.$(e) $(s)/real_precision.$(e) $(s)/soil_methods.$(e) $(s)/misc_values.$(e) $(s)/tuning_parameters.$(e)
	$(FC) -o $(o)/soil_methods.o -c $(s)/soil_methods.$(e) -J$(m)

$(o)/light_methods.o: $(s)/light_methods.$(e) $(s)/real_precision.$(e)
	$(FC) -o $(o)/light_methods.o -c $(s)/light_methods.$(e) -J$(m)

$(o)/sdgvm.o: $(s)/sdgvm.$(e) $(s)/real_precision.$(e) $(s)/dims.$(e) $(s)/system_state.$(e) $(s)/pft_parameters.$(e) $(s)/state_methods.$(e) $(s)/sdgvm1.$(e) $(s)/read_input.$(e) $(s)/phenology_methods.$(e) $(s)/doly.$(e) $(s)/hydrology_methods.$(e) $(s)/light_methods.$(e) $(s)/soil_methods.$(e) $(s)/misc_values.$(e)
	$(FC) -o $(o)/sdgvm.o -c $(s)/sdgvm.$(e) -J$(m)

clean:
	rm -f $b/sdgvm.exe $b/sdgvm $m/*.mod $o/*.o $s/*.obj
