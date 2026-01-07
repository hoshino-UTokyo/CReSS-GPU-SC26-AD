# Include common compiler configuration
include ../Makefile.inc

# Targets
BENCHMARK = kernel_benchmark
DATAGEN = generate_data

all: $(BENCHMARK) $(DATAGEN)

$(BENCHMARK): kernel_benchmark.f90
	$(FC) $(FFLAGS) -o $@ $<

$(DATAGEN): generate_data.f90
	$(FC) $(FFLAGS) -o $@ $<

data: $(DATAGEN)
	./$(DATAGEN)

run: $(BENCHMARK) data
	./$(BENCHMARK)

clean:
	rm -f $(BENCHMARK) $(DATAGEN) *.o *.mod

cleandata:
	rm -rf data

.PHONY: all data run clean cleandata
