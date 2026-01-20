# Include common compiler configuration
include ../Makefile.inc

#======================================================================
# Makefile for Kernel Benchmark
#======================================================================

# Targets
BENCHMARK = kernel_benchmark

# Default target
all: $(BENCHMARK)

# Benchmark executable
$(BENCHMARK): kernel_benchmark.f90
	$(FC) $(FFLAGS) -o $@ $<

# Run benchmark
run: $(BENCHMARK)
	./$(BENCHMARK)

# Clean build files
clean:
	rm -f $(BENCHMARK) *.o *.mod

.PHONY: all run clean
