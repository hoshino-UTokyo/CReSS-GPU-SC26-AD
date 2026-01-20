# Include common compiler configuration
include ../Makefile.inc

# Kernel Benchmark Makefile: bcycley
#
# Source: Src/bcycley.f90
# Description: Set periodic boundary conditions in y direction

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
