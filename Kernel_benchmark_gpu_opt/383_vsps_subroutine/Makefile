#======================================================================
# Makefile for Kernel Benchmark (OpenACC GPU version)
#======================================================================

# Compiler
FC = nvfortran

# Flags for OpenACC GPU with unified memory
FFLAGS = -fast -acc -gpu=managed -Minfo=accel -Mbyteswapio -mcmodel=medium

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
