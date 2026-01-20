# Makefile for sfcflx kernel benchmark

# Compiler settings
FC = gfortran
FFLAGS = -O3 -fopenmp -march=native -fconvert=big-endian

# For Fujitsu compiler (comment out gfortran and use these)
# FC = frtpx
# FFLAGS = -O3 -Kopenmp

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
