# Include common compiler configuration
include ../Makefile.inc

# Kernel Benchmark Makefile: trilat
#
# Source: Src/trilat.f90
# Description: Calculate Coriolis parameters from latitude

# Targets
BENCHMARK = kernel_benchmark

all: $(BENCHMARK)

$(BENCHMARK): kernel_benchmark.f90
	$(FC) $(FFLAGS) -o $@ $<

run: $(BENCHMARK)
	./$(BENCHMARK)

clean:
	rm -f $(BENCHMARK) *.o *.mod

.PHONY: all run clean
