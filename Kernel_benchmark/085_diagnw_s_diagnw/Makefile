# Include common compiler configuration
include ../Makefile.inc

# Kernel Benchmark Makefile: diagnw
#
# Source: Src/diagnw.f90
# Description: Compute diagnostic concentrations of water hydrometeors

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
