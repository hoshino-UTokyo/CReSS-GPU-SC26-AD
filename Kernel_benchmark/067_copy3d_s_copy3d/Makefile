# Include common compiler configuration
include ../Makefile.inc

#======================================================================
# Makefile for copy3d Kernel Benchmark
#======================================================================


BENCHMARK = kernel_benchmark
GENERATOR = generate_data

all: $(BENCHMARK)

$(BENCHMARK): kernel_benchmark.f90
	$(FC) $(FFLAGS) -o $@ $<

$(GENERATOR): generate_data.f90
	$(FC) -O2 -o $@ $<

data: $(GENERATOR)
	./$(GENERATOR)

run: $(BENCHMARK)
	./$(BENCHMARK)

run-scaling: $(BENCHMARK)
	@echo "=== Thread Scaling Test ==="
	@for t in 1 2 4 8 16 32; do \
		echo "--- Threads: $$t ---"; \
		OMP_NUM_THREADS=$$t ./$(BENCHMARK); \
	done

clean:
	rm -f $(BENCHMARK) $(GENERATOR) *.o *.mod

cleanall: clean
	rm -rf data

.PHONY: all data run run-scaling clean cleanall
