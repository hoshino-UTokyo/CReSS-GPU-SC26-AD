# Include common compiler configuration
include ../Makefile.inc

# Kernel Benchmark Makefile: coriuv
#
# Source: Src/coriuv.f90
# Description: Calculate Coriolis force for u and v velocity equations
# NOTE: Uses dump data from test_real/kernel_dump/coriuv (symlinked to ./data)

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
