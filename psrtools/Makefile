#PROGS = autotoa normalize_rms running_mean_sub
PROGS = polifitgsl running_mean_sub
all: $(PROGS)
install: $(PROGS)
	cp $(PROGS) $(PSR)/bin

#autotoa: autotoa.C
#	$(CXX) -g `psrchive --cflags` -o autotoa autotoa.C `psrchive --libs` -I/usr/local/include/ 

#parallactic: parallactic.C
#	$(CXX) -g `psrchive --cflags` -o parallactic parallactic.C `psrchive --libs` -I/usr/local/include/

#check_model_phase: check_model_phase.C
#	$(CXX) -g `psrchive --cflags` -o check_model_phase check_model_phase.C `psrchive --libs` -I/usr/local/include/

#normalize_rms: normalize_rms.C
#	$(CXX) -g `psrchive --cflags` -o normalize_rms normalize_rms.C `psrchive --libs` -I/usr/local/include/

polifitgsl: polifitgsl.C
	$(CXX)  -c polifitgsl.C -I/usr/local/include/ `pkg-config --cflags gsl`

running_mean_sub: running_mean_sub.C
	$(CXX) -g `psrchive --cflags` -o running_mean_sub running_mean_sub.C polifitgsl.o  `psrchive --libs` `pkg-config --libs gsl` -lgsl -lgslcblas  -lm

#dm_stat_search: dm_stat_search.C
#	$(CXX) -O2 `dspsr_cflags` `psrchive --cflags` -o $@ $< `dspsr_ldflags` 

clean:
	rm -f $(PROGS) *.o
