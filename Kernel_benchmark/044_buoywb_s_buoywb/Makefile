# Include common compiler configuration
include ../Makefile.inc

# Kernel Benchmark Makefile: buoywb
#
# Source: Src/buoywb.f90
# Description: Calculate buoyancy forcing for vertical velocity

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
