# Include common compiler configuration
include ../Makefile.inc

# Kernel Benchmark Makefile: bc8u
#
# Source: Src/bc8u.f90
# Description: Set west/east boundary conditions for variable at u points

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
