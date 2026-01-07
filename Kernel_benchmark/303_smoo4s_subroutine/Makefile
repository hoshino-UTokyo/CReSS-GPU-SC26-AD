# Include common compiler configuration
include ../Makefile.inc

# Makefile for smoo4s kernel benchmark


# For Fujitsu compiler (comment out gfortran lines and uncomment these)
# FC = frtpx
# FFLAGS = -Kfast,openmp

# Targets
BENCHMARK = kernel_benchmark

all: $(BENCHMARK)

$(BENCHMARK): kernel_benchmark.f90
	$(FC) $(FFLAGS) -o $@ $<

run: $(BENCHMARK)
	OMP_NUM_THREADS=4 ./$(BENCHMARK)

clean:
	rm -f $(BENCHMARK) *.o *.mod

.PHONY: all run clean
