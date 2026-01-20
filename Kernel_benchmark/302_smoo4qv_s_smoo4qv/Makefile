# Include common compiler configuration
include ../Makefile.inc

# Kernel Benchmark Makefile: smoo4qv
#
# Source: Src/smoo4qv.f90
# Description: Fourth-order smoothing for moisture variables

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
