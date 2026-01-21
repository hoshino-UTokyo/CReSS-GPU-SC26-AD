#-----------------------------------------------------------------------
# Makefile for Kernel I/O Library
#-----------------------------------------------------------------------

# Compiler settings (can be overridden)
FC ?= gfortran
FFLAGS ?= -O2 -Wall

# Library name
LIB = libkernel_io.a

# Module objects
OBJS = m_kernel_io.o m_dump_kernel.o

# Default target
all: $(LIB)

# Build static library
$(LIB): $(OBJS)
	ar rcs $@ $^
	@echo "Library $(LIB) created successfully"

# Compile modules
m_kernel_io.o: m_kernel_io.f90
	$(FC) $(FFLAGS) -c $<

m_dump_kernel.o: m_dump_kernel.f90 m_kernel_io.o
	$(FC) $(FFLAGS) -c $<

# Test program
test: test_kernel_io
	./test_kernel_io

test_kernel_io: test_kernel_io.o $(LIB)
	$(FC) $(FFLAGS) -o $@ $< -L. -lkernel_io

test_kernel_io.o: test_kernel_io.f90 m_kernel_io.o
	$(FC) $(FFLAGS) -c $<

# Clean
clean:
	rm -f *.o *.mod $(LIB) test_kernel_io
	rm -rf test_dump_output

# Install (copy to specified location)
install: $(LIB)
	@if [ -z "$(INSTALL_DIR)" ]; then \
		echo "Usage: make install INSTALL_DIR=/path/to/install"; \
	else \
		mkdir -p $(INSTALL_DIR)/lib $(INSTALL_DIR)/include; \
		cp $(LIB) $(INSTALL_DIR)/lib/; \
		cp *.mod $(INSTALL_DIR)/include/; \
		echo "Installed to $(INSTALL_DIR)"; \
	fi

.PHONY: all clean test install
