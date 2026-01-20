# Include common compiler configuration
include ../Makefile.inc

# Kernel Benchmark Makefile: advs
#
# Source: Src/advs.f90
# Description: Calculate scalar advection


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
