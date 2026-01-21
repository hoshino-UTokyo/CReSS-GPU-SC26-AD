# Makefile for 4D dump segfault reproduction tests
# Use same compiler as solver (NVIDIA Fortran)

FC = nvfortran
FFLAGS = -O2

.PHONY: all clean run_basic run_cumulative

all: test_4d_dump test_cumulative_dump

test_4d_dump: test_4d_dump.f90
	$(FC) $(FFLAGS) -o $@ $<

test_cumulative_dump: test_cumulative_dump.f90
	$(FC) $(FFLAGS) -o $@ $<

# Basic test - single 4D dump (should pass)
run_basic: test_4d_dump
	./test_4d_dump

# Cumulative test - simulates solver with many dumps
# Usage: make run_cumulative N=700
N ?= 100
run_cumulative: test_cumulative_dump
	./test_cumulative_dump $(N)

clean:
	rm -f test_4d_dump test_cumulative_dump *.mod *.o
	rm -rf test_output
