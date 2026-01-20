# Include common compiler configuration
include ../Makefile.inc

# Kernel Benchmark Makefile: chksat
#
# Source: Src/chksat.f90
# Description: Check and avoid super saturation mixing ratio

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
