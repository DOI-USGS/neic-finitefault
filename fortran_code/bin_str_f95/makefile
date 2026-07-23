BIN=./
OBJS3= fft.o mult_single.o interporation.o cmodel.o ddis.o  Read_green.o bpfilter.o
OBJS5= fft.o green_cgps.o interporation.o cmodel.o ddis.o  Read_green.o bpfilter.o
OBJS4= fft.o get_green_fun.o interporation.o cmodel.o ddis.o  Read_green.o bpfilter.o
modules = constants.o vel_model_data.o model_parameters.o
modules2 = $(modules) geodesics.o retrieve_gf.o wave_travel.o bpfilter.o fast_fourier_trans.o rad_pattern.o store_gf.o
# -mcmodel is x86-64-specific; not supported (and not needed) on macOS arm64
ifeq ($(shell uname -s),Darwin)
MCMODEL_LARGE =
else
MCMODEL_LARGE = -mcmodel=large
endif
F95= gfortran -Wall -O2 -ffixed-line-length-132 -fno-range-check $(MCMODEL_LARGE) -fbacktrace -g -fno-rerun-loop-opt

 
all: get_strong_motion clean


get_strong_motion: $(modules2) get_strong_motion.o
	$(F95) -o get_strong_motion $(modules2) get_strong_motion.o


%.mod: $.o $.f95
	$(f95) -c $<

#constants.mod: constants.o constants.f95
#	$(F95) -c constants.f95
#
#retrieve_gf.mod: retrieve_gf.o retrieve_gf.f95
#	$(F95) -c retrieve_gf.f95
#
#vel_model_data.mod: vel_model_data.o vel_model_data.f95
#	$(F95) -c vel_model_data.f95
#
#wave_travel.mod: wave_travel.o wave_travel.f95
#	$(F95) -c wave_travel.f95
#
#geodesics.mod: geodesics.o geodesics.f95
#	$(F95) -c geodesics.f95
#
#bpfilter.mod: bpfilter.o bpfilter.f95
#	$(F95) -c bpfilter.f95
#
#fast_fourier_trans.mod: fast_fourier_trans.o fast_fourier_trans.f95
#	$(F95) -c fast_fourier_trans.f95
#
#rad_pattern.mod: rad_pattern.o rad_pattern.f95
#	$(F95) -c rad_pattern.f95


%.o: %.f95
	$(F95) -c $<

#get_strong_motion.o: get_strong_motion.f95
#	$(F95) -c get_strong_motion.f95
#
#retrieve_gf.o: retrieve_gf.f95
#	$(F95) -c retrieve_gf.f95
#
#constants.o: constants.f95
#	$(F95) -c constants.f95
#
#vel_model_data.o: vel_model_data.f95
#	$(F95) -c vel_model_data.f95
#
#wave_travel.o: wave_travel.f95
#	$(F95) -c wave_travel.f95
#
#geodesics.o: geodesics.f95
#	$(F95) -c geodesics.f95
#
#bpfilter.o: bpfilter.f95
#	$(F95) -c bpfilter.f95
#
#fast_fourier_trans.o: fast_fourier_trans.f95
#	$(F95) -c fast_fourier_trans.f95
#
#rad_pattern.o: rad_pattern.f95
#	$(F95) -c rad_pattern.f95


clean:
	rm -f *.mod 
	rm -f *.o 
