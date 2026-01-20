# Include common compiler configuration
include ../Makefile.inc

# Kernel Benchmark Makefile: sheartke
#
# Source: Src/sheartke.f90
# Description: Calculate shear production in TKE equation

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
