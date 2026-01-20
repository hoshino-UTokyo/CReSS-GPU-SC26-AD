# Include common compiler configuration
include ../Makefile.inc

# Kernel Benchmark Makefile: bc8w
#
# Source: Src/bc8w.f90
# Description: Set bottom/top boundary conditions for variable at w points

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
