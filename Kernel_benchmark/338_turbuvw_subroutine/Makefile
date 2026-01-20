# Include common compiler configuration
include ../Makefile.inc

# Kernel Benchmark Makefile: turbuvw
#
# Source: Src/turbuvw.f90
# Description: Calculate velocity turbulent mixing
# NOTE: Uses dump data from test_real/kernel_dump/turbuvw (symlinked to ./data)

BENCHMARK = kernel_benchmark

all: $(BENCHMARK)

$(BENCHMARK): kernel_benchmark.f90
	$(FC) $(FFLAGS) -o $@ $<

run: $(BENCHMARK)
	./$(BENCHMARK)

clean:
	rm -f $(BENCHMARK) *.o *.mod

.PHONY: all run clean
