# Include common compiler configuration
include ../Makefile.inc

# Kernel Benchmark Makefile: advbspt
#
# Source: Src/advbspt.f90
# Description: Advection of base-state potential temperature
# NOTE: Uses dump data from test_real/kernel_dump/advbspt (symlinked to ./data)

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
