FC = gfortran
FFLAGS = -O3 -fopenmp -march=native
BENCHMARK = kernel_benchmark
all: $(BENCHMARK)
$(BENCHMARK): kernel_benchmark.f90
	$(FC) $(FFLAGS) -o $@ $<
run: $(BENCHMARK)
	./$(BENCHMARK)
clean:
	rm -f $(BENCHMARK) *.o *.mod
.PHONY: all run clean
