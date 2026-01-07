# Include common compiler configuration
include ../Makefile.inc

all: kernel_benchmark

kernel_benchmark: kernel_benchmark.f90
	$(FC) $(FFLAGS) -o $@ $<

generate_data: generate_data.f90
	$(FC) -O2 -o $@ $<

data: generate_data
	./generate_data

run: kernel_benchmark
	./kernel_benchmark

clean:
	rm -f kernel_benchmark generate_data *.o *.mod

.PHONY: all data run clean
