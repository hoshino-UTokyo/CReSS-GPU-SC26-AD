# Include common compiler configuration
include ../Makefile.inc

# Kernel Benchmark Makefile: var8uvw
#
# Source: Src/var8uvw.f90
# Description: Average variable to u, v, and w staggered grid points

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
