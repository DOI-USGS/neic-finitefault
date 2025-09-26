MY_BIN=.

modules = constants.o modelling_inputs.o model_parameters.o wavelet_param.o wavelets.o get_stations_data.o retrieve_surf_gf.o rad_pattern.o geodesics.o rise_time.o retrieve_gf.o save_forward.o
modules2 = $(modules) static_data.o imagery_data.o misfit_eval.o regularization.o random_gen.o annealing.o annealing_static.o

#Linux
#F95= ifort -shared-intel -mcmodel=large -O2 -assume byterecl -p#-p -132 -shared-intel
#F77= ifort -O3 -132  -assume byterecl
F95= gfortran -O3 -fno-range-check -mcmodel=medium -g #-mavx2 -Wall -ftree-loop-vectorize -ftree-loop-distribution#-fno-omit-frame-pointer -ffpe-trap=overflow,invalid#-ffpe-trap=overflow -funroll-loops -fcheck=all
F77= gfortran -O1 -fno-range-check -mcmodel=large -g -fbacktrace


all: run_forward run_modelling green_tele clean


run_forward: $(modules) static_data.o imagery_data.o run_forward.o
	$(F95) -o run_forward $(modules) static_data.o imagery_data.o run_forward.o

run_modelling: $(modules2) ffm_methods.o run_modelling.o
	$(F95) -o run_modelling -fopenmp $(modules2) ffm_methods.o run_modelling.o #-fopenmp

green_tele: bpfilter.o ddis.o green_tele.o
	$(F77) -o green_tele bpfilter.o ddis.o green_tele.o

ffm_methods.mod: ffm_methods.o ffm_methods.f95
	$(F95) -fopenmp -c ffm_methods.f95 #-fopenmp

annealing.mod: annealing.o annealing.f95
	$(F95) -fopenmp -c annealing.f95 #-fopenmp

%.mod: %.o %.f95
	$(F95) -c $<


bpfilter.o : bpfilter.f
	$(F77) -c bpfilter.f

green_tele.o : green_tele.f
	$(F77) -c green_tele.f

ddis.o : ddis.f
	$(F77) -c ddis.f

ffm_methods.o: ffm_methods.f95
	$(F95) -fopenmp -c ffm_methods.f95

annealing.o: annealing.f95
	$(F95) -fopenmp -c annealing.f95

%.o: %.f95
	$(F95) -c $<

 
clean:
	\rm -f *.o
	\rm -f *.mod
