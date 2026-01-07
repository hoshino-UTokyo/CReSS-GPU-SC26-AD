# Include common compiler configuration
include ../Makefile.inc

# Kernel Benchmark Makefile: advbspi
#
# Source: Src/advbspi.f90
# Description: Calculate base state pressure advection


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
