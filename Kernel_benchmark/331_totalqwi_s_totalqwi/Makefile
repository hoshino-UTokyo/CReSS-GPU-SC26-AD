# Include common compiler configuration
include ../Makefile.inc

# Kernel Benchmark Makefile: totalqwi
#
# Source: Src/totalqwi.f90
# Description: Calculate total water and ice mixing ratio

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
