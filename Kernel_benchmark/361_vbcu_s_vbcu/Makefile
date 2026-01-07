# Include common compiler configuration
include ../Makefile.inc

#======================================================================
# Makefile for vbcu Kernel Benchmark
#======================================================================


BENCHMARK = kernel_benchmark
GENERATOR = generate_data

all: $(BENCHMARK)

$(BENCHMARK): kernel_benchmark.f90
	$(FC) $(FFLAGS) -o $@ $<

$(GENERATOR): generate_data.f90
	$(FC) -O2 -o $@ $<

data: $(GENERATOR)
	./$(GENERATOR)

run: $(BENCHMARK)
	./$(BENCHMARK)

clean:
	rm -f $(BENCHMARK) $(GENERATOR) *.o *.mod

cleanall: clean
	rm -rf data

.PHONY: all data run clean cleanall
