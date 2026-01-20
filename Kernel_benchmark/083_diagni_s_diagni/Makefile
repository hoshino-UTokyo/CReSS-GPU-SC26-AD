# Include common compiler configuration
include ../Makefile.inc

#======================================================================
# Makefile for diagni Kernel Benchmark
#======================================================================

# Compiler settings (adjust for your environment)

# For Intel compiler:
# FC = ifort
# FFLAGS = -O3 -qopenmp -xHost

# For NVIDIA HPC SDK (OpenACC):
# FC = nvfortran
# FFLAGS = -O3 -mp -Minfo=mp

# Targets
BENCHMARK = kernel_benchmark
GENERATOR = generate_data

# Default target
all: $(BENCHMARK)

# Benchmark executable
$(BENCHMARK): kernel_benchmark.f90
	$(FC) $(FFLAGS) -o $@ $<

# Data generator executable (if exists)
$(GENERATOR): generate_data.f90
	$(FC) -O2 -o $@ $<

# Generate test data
data: $(GENERATOR)
	./$(GENERATOR)

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
	rm -f $(BENCHMARK) $(GENERATOR) *.o *.mod

# Clean all including data
cleanall: clean
	rm -rf data benchmark.conf

# Help
help:
	@echo "Kernel Benchmark: diagni"
	@echo ""
	@echo "Targets:"
	@echo "  all          - Build benchmark"
	@echo "  data         - Generate test data (requires generate_data.f90)"
	@echo "  run          - Run benchmark"
	@echo "  run-scaling  - Run thread scaling test"
	@echo "  clean        - Remove executables"
	@echo "  cleanall     - Remove executables and data"
	@echo ""
	@echo "Environment variables:"
	@echo "  OMP_NUM_THREADS - Number of OpenMP threads"
	@echo "  OMP_SCHEDULE    - OpenMP schedule (e.g., 'static', 'dynamic')"

.PHONY: all data run run-scaling clean cleanall help
