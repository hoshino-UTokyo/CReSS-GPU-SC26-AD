# Include common compiler configuration
include ../Makefile.inc

# Kernel Benchmark Makefile: adjstq
#
# Source: Src/adjstq.f90
# Description: Adjust hydrometeor mixing ratios
# NOTE: Uses dump data from test_real/kernel_dump/adjstq (symlinked to ./data)

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
