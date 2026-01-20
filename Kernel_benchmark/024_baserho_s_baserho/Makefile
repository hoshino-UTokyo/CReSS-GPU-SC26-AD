# Include common compiler configuration
include ../Makefile.inc

# Kernel Benchmark Makefile: baserho
#
# Source: Src/baserho.f90
# Description: Calculate base-state air density
# NOTE: Uses dump data from test_real/kernel_dump/baserho (symlinked to ./data)

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
