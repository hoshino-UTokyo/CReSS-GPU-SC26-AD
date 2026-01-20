# Include common compiler configuration
include ../Makefile.inc

# Kernel Benchmark Makefile: bc8v
#
# Source: Src/bc8v.f90
# Description: Set south/north boundary conditions for variable at v points

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
