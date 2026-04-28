# Include common compiler configuration
include ../Makefile.inc

# Kernel Benchmark Makefile: tkeflx
#
# Source: Src/tkeflx.f90

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
