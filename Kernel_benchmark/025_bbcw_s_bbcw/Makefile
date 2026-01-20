# Include common compiler configuration
include ../Makefile.inc

# Kernel Benchmark Makefile: bbcw
#
# Source: Src/bbcw.f90
# Description: Bottom boundary condition for vertical velocity
# NOTE: Uses dump data from test_real/kernel_dump/bbcw (symlinked to ./data)

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
