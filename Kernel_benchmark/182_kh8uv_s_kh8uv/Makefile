# Include common compiler configuration
include ../Makefile.inc

# Kernel Benchmark Makefile: kh8uv
#
# Source: Src/kh8uv.f90
# Description: Set horizontal eddy diffusivity at u and v points

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
