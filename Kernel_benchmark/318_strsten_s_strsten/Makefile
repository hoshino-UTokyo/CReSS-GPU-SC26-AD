# Include common compiler configuration
include ../Makefile.inc

# Kernel Benchmark Makefile: strsten
#
# Source: Src/strsten.f90
# Description: Calculate strain tensor components for stress calculation
# NOTE: Uses dump data from test_real/kernel_dump/strsten (symlinked to ./data)

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
