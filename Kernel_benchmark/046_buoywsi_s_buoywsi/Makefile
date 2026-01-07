# Include common compiler configuration
include ../Makefile.inc

# Kernel Benchmark Makefile: buoywsi
#
# Source: Src/buoywsi.f90
# Description: Calculate buoyancy for small time steps


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
