# Include common compiler configuration
include ../Makefile.inc

#======================================================================
# Makefile for setsfc Kernel Benchmark
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

# Run with different thread counts
run-scaling: $(BENCHMARK)
	@echo "=== Thread Scaling Test ==="
	@for t in 1 2 4 8 16 32; do \
		echo "--- Threads: $$t ---"; \
		OMP_NUM_THREADS=$$t ./$(BENCHMARK); \
	done

# Clean build files
clean:
	rm -f $(BENCHMARK) *.o *.mod

.PHONY: all run run-scaling clean
