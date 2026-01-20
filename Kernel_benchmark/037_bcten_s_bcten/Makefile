# Include common compiler configuration
include ../Makefile.inc

# Kernel Benchmark Makefile: bcten
#
# Source: Src/bcten.f90
# Description: Set bottom/top boundary conditions for tensor

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
