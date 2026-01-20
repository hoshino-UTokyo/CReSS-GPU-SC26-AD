# Include common compiler configuration
include ../Makefile.inc

# Kernel Benchmark Makefile: vbcp
#
# Source: Src/vbcp.f90
# Description: Set vertical boundary conditions for pressure

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
